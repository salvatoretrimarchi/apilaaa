use crate::say;
use anyhow::{anyhow, Context, Result};
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::process::Command;
use tiff::encoder::{colortype, Rational, SRational, TiffEncoder};
use tiff::tags::Tag;

/// Stratified samples (~1M) of the requested channel, positive and finite.
fn sample_channel(rgb: &[f32], channel: usize) -> Vec<f32> {
    let n_pix = rgb.len() / 3;
    let stride = (n_pix / 1_000_000).max(1);
    (0..n_pix)
        .step_by(stride)
        .map(|i| rgb[i * 3 + channel])
        .filter(|v| v.is_finite() && *v > 0.0)
        .collect()
}

/// Percentile of an RGB channel over ~1M stratified samples.
fn percentile_channel(rgb: &[f32], channel: usize, percentile: f32) -> f32 {
    let mut samples = sample_channel(rgb, channel);
    if samples.is_empty() {
        return 0.0;
    }
    let p = percentile.clamp(0.0, 100.0) / 100.0;
    let k = ((samples.len() as f32 * p) as usize).min(samples.len() - 1);
    let (_, nth, _) = samples.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
    *nth
}

/// Median and MAD (median absolute deviation) of an RGB channel.
fn median_mad_channel(rgb: &[f32], channel: usize) -> (f32, f32) {
    let mut samples = sample_channel(rgb, channel);
    if samples.is_empty() {
        return (0.0, 0.0);
    }
    let k = samples.len() / 2;
    let (_, med, _) = samples.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
    let med = *med;
    let mut devs: Vec<f32> = samples.iter().map(|v| (v - med).abs()).collect();
    let km = devs.len() / 2;
    let (_, mad, _) = devs.select_nth_unstable_by(km, |a, b| a.partial_cmp(b).unwrap());
    (med, *mad)
}

/// Parameters of the linear histogram stretch, derived automatically from
/// the analysis of the raws. Both ends are computed **per channel** so that
/// every colour clips the same fraction of its own distribution and no tints
/// appear from imbalance (if `white` were global, the dominant channel —
/// green in astro — would see its range compressed against a foreign one and
/// would lose midtones that R/B would keep).
#[derive(Clone, Copy, Debug)]
pub struct StretchParams {
    pub white: [f32; 3],
    pub blacks: [f32; 3],
}

/// `WHITE` at 99.9 saturates only the brightest 0.1% of each channel.
///
/// `BLACKS` is anchored at the **left foot of the dominant peak** (the sky in
/// astro), not at its mode. For a unimodal ETTR astro-like distribution the
/// median falls practically on the mode and the MAD describes its width, so
/// `median − k·MAD` lands on the left skirt without biting into the body of
/// the peak. k≈3 in MAD units is equivalent to ~2σ for a Gaussian — it cuts
/// bias and pure noise but preserves sky, milky way and nebula.
const WHITE_PERCENTILE: f32 = 99.9;
const BLACKS_K_MAD: f32 = 3.0;

/// Analyses an interleaved RGB buffer and returns the stretch parameters.
/// It is applied to the final stack (already balanced via `frame.wb` in
/// `stack::add`) so that the per-channel peaks fall in comparable ranges.
pub fn analyze_stretch(rgb: &[f32]) -> StretchParams {
    let (m0, mad0) = median_mad_channel(rgb, 0);
    let (m1, mad1) = median_mad_channel(rgb, 1);
    let (m2, mad2) = median_mad_channel(rgb, 2);
    StretchParams {
        white: [
            percentile_channel(rgb, 0, WHITE_PERCENTILE),
            percentile_channel(rgb, 1, WHITE_PERCENTILE),
            percentile_channel(rgb, 2, WHITE_PERCENTILE),
        ],
        blacks: [
            (m0 - BLACKS_K_MAD * mad0).max(0.0),
            (m1 - BLACKS_K_MAD * mad1).max(0.0),
            (m2 - BLACKS_K_MAD * mad2).max(0.0),
        ],
    }
}

/// XYZ→camera matrix (D65) ×10000, row-major 3x3. Official LibRaw/dcraw
/// values. The generic fallback is the identity, which produces flat but
/// valid colours.
fn color_matrix1(camera_model: &str) -> [(i32, i32); 9] {
    let d = 10000;
    let raw: [i32; 9] = match camera_model {
        "SONY ILCE-7M3" => [7374, -2389, -551, -5435, 13162, 2519, -1006, 1795, 6552],
        "SONY ILCE-7M4" => [7460, -2365, -588, -5687, 13442, 2474, -624, 1156, 6584],
        "SONY ILCE-7RM3" | "SONY ILCE-7RM3A" => [6640, -1847, -503, -5238, 13010, 2474, -993, 1859, 6861],
        "SONY ILCE-7RM4" | "SONY ILCE-7RM4A" => [7662, -2686, -660, -5361, 13391, 2221, -1150, 1826, 7494],
        "SONY ILCE-6400" => [7657, -2847, -607, -4083, 11966, 2389, -684, 1418, 5844],
        _ => [10000, 0, 0, 0, 10000, 0, 0, 0, 10000],
    };
    [
        (raw[0], d), (raw[1], d), (raw[2], d),
        (raw[3], d), (raw[4], d), (raw[5], d),
        (raw[6], d), (raw[7], d), (raw[8], d),
    ]
}

/// Writes a linear RGB 16-bit DNG in [0, 65535] (black=0, white=65535)
/// preserving the original linear mapping of the sensor. Adds the tags
/// needed for darktable/Siril/Lightroom to process the file through their
/// raw pipeline (colour matrix, demosaic skipped).
///
/// The `rgb` buffer is assumed to be **already white-balanced** (the WB is
/// baked into the stack via `frame.wb`), so `AsShotNeutral` is written as
/// `[1,1,1]` and the developer does not apply WB again on top.
///
/// `stretch`: if Some, applies the linear stretch derived from the automatic
/// analysis of the balanced stack. The `blacks[c]` and `white[c]` points are
/// **baked into the data**: every pixel is remapped as
/// `(v - blacks[c]) / (white[c] - blacks[c]) · 65535`, clamped to [0, 65535].
/// BlackLevel is written as 0 because the normalization is already applied;
/// this way the clipping is visible even if the developer ignores
/// BlackLevel[3] in LinearRaw DNGs. None keeps the sensor-native scale
/// without clipping.
pub fn write_dng(
    path: &Path,
    rgb: &[f32],
    w: usize,
    h: usize,
    camera_model: &str,
    stretch: Option<StretchParams>,
) -> Result<()> {
    write_dng_impl(path, rgb, w, h, camera_model, stretch, true)
}

/// Same as `write_dng` but without printing the stretch parameters (for
/// batch exports where the stretch is shared and has already been reported).
pub fn write_dng_quiet(
    path: &Path,
    rgb: &[f32],
    w: usize,
    h: usize,
    camera_model: &str,
    stretch: Option<StretchParams>,
) -> Result<()> {
    write_dng_impl(path, rgb, w, h, camera_model, stretch, false)
}

fn write_dng_impl(
    path: &Path,
    rgb: &[f32],
    w: usize,
    h: usize,
    camera_model: &str,
    stretch: Option<StretchParams>,
    verbose: bool,
) -> Result<()> {
    let (blacks, white) = match stretch {
        None => ([0.0f32; 3], [1.0f32; 3]),
        Some(s) => (s.blacks, [
            s.white[0].max(1e-6),
            s.white[1].max(1e-6),
            s.white[2].max(1e-6),
        ]),
    };
    let inv_range = [
        65535.0 / (white[0] - blacks[0]).max(1e-6),
        65535.0 / (white[1] - blacks[1]).max(1e-6),
        65535.0 / (white[2] - blacks[2]).max(1e-6),
    ];
    if verbose && stretch.is_some() {
        let pct = |c: usize| 100.0 * blacks[c] / white[c];
        say!(
            "stretch: white=R {:.6} G {:.6} B {:.6}  blacks=R {:.6} G {:.6} B {:.6} ({:.1}% / {:.1}% / {:.1}% of the range)",
            white[0], white[1], white[2],
            blacks[0], blacks[1], blacks[2],
            pct(0), pct(1), pct(2)
        );
    }

    let mut u16_buf = Vec::with_capacity(rgb.len());
    for chunk in rgb.chunks_exact(3) {
        for c in 0..3 {
            let s = ((chunk[c] - blacks[c]) * inv_range[c]).clamp(0.0, 65535.0);
            u16_buf.push((s + 0.5) as u16);
        }
    }

    let f = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let bw = BufWriter::new(f);
    let mut enc = TiffEncoder::new(bw)?;
    let mut img = enc.new_image::<colortype::RGB16>(w as u32, h as u32)?;
    let de = img.encoder();

    // --- TIFF/DNG structure ---
    // NewSubfileType = 0 (main image, full-res)
    de.write_tag(Tag::Unknown(254), 0u32)?;
    // PhotometricInterpretation = 34892 (LinearRaw) — enables the raw
    // pipeline. Overrides the RGB=2 value set by new_image().
    de.write_tag(Tag::PhotometricInterpretation, 34892u16)?;

    // --- Mandatory DNG tags ---
    de.write_tag(Tag::Unknown(50706), &[1u8, 4, 0, 0][..])?;   // DNGVersion 1.4.0.0
    de.write_tag(Tag::Unknown(50707), &[1u8, 4, 0, 0][..])?;   // DNGBackwardVersion 1.4.0.0
    de.write_tag(Tag::Unknown(50708), camera_model)?;          // UniqueCameraModel

    // --- Sensor levels (after processing) ---
    // BlackLevelRepeatDim LONG[2] = [1, 1] (one value per plane)
    de.write_tag(Tag::Unknown(50713), &[1u32, 1][..])?;
    // BlackLevel LONG[3] = [0, 0, 0]: the black clipping (if any) is baked
    // directly into the u16 data, not delegated to the developer.
    de.write_tag(Tag::Unknown(50714), &[0u32, 0, 0][..])?;
    de.write_tag(Tag::Unknown(50717), &[65535u32][..])?;       // WhiteLevel

    // --- Active area and crop ---
    // ActiveArea LONG[4] (top, left, bottom, right)
    de.write_tag(Tag::Unknown(50829), &[0u32, 0, h as u32, w as u32][..])?;
    de.write_tag(
        Tag::Unknown(50719),
        &[Rational { n: 0, d: 1 }, Rational { n: 0, d: 1 }][..],
    )?;
    de.write_tag(
        Tag::Unknown(50720),
        &[Rational { n: w as u32, d: 1 }, Rational { n: h as u32, d: 1 }][..],
    )?;

    // --- Color ---
    // CalibrationIlluminant1 = 21 (D65)
    de.write_tag(Tag::Unknown(50778), 21u16)?;
    // ColorMatrix1 SRATIONAL[9] — XYZ (D65) → native camera
    let m = color_matrix1(camera_model);
    let color_matrix: Vec<SRational> = m.iter().map(|&(n, d)| SRational { n, d }).collect();
    de.write_tag(Tag::Unknown(50721), &color_matrix[..])?;
    // AnalogBalance RATIONAL[3] = [1, 1, 1]
    de.write_tag(
        Tag::Unknown(50727),
        &[
            Rational { n: 1, d: 1 },
            Rational { n: 1, d: 1 },
            Rational { n: 1, d: 1 },
        ][..],
    )?;
    // AsShotNeutral RATIONAL[3] = [1, 1, 1]: the WB is already baked into
    // the data (see stack::add), so we tell the developer not to balance
    // again on top.
    de.write_tag(
        Tag::Unknown(50728),
        &[
            Rational { n: 1, d: 1 },
            Rational { n: 1, d: 1 },
            Rational { n: 1, d: 1 },
        ][..],
    )?;
    // BaselineExposure = 0 (no correction by default)
    de.write_tag(Tag::Unknown(50730), SRational { n: 0, d: 1 })?;

    img.write_data(&u16_buf)?;
    Ok(())
}

/// Copies EXIF + MakerNotes + XMP + ICC from the source ARW to the
/// destination DNG using exiftool, excluding the structural and DNG tags we
/// have just written.
pub fn copy_metadata(source_arw: &Path, dest_dng: &Path) -> Result<()> {
    let out = Command::new("exiftool")
        .args([
            "-TagsFromFile",
        ])
        .arg(source_arw)
        .args([
            "-EXIF:all",
            "-MakerNotes:all",
            "-XMP:all",
            "-ICC_Profile:all",
            // Exclude tags we control ourselves or that are structural to the DNG.
            "--PhotometricInterpretation",
            "--NewSubfileType",
            "--SubfileType",
            "--ImageWidth",
            "--ImageHeight",
            "--ImageLength",
            "--BitsPerSample",
            "--Compression",
            "--SamplesPerPixel",
            "--RowsPerStrip",
            "--StripOffsets",
            "--StripByteCounts",
            "--PlanarConfiguration",
            "--SampleFormat",
            "--CFAPattern",
            "--CFAPattern2",
            "--CFARepeatPatternDim",
            "--CFAPlaneColor",
            "--CFALayout",
            "--BlackLevel",
            "--WhiteLevel",
            "--BlackLevelRepeatDim",
            "--DNGVersion",
            "--DNGBackwardVersion",
            "--UniqueCameraModel",
            "--LocalizedCameraModel",
            "--AsShotNeutral",
            "--ColorMatrix1",
            "--ColorMatrix2",
            "--CalibrationIlluminant1",
            "--CalibrationIlluminant2",
            "--AnalogBalance",
            "--BaselineExposure",
            "--ActiveArea",
            "--DefaultCropOrigin",
            "--DefaultCropSize",
            "-overwrite_original",
        ])
        .arg(dest_dng)
        .output()
        .context("running exiftool (is it installed?)")?;

    if !out.status.success() {
        return Err(anyhow!(
            "exiftool failed ({}): {}",
            out.status,
            String::from_utf8_lossy(&out.stderr)
        ));
    }
    Ok(())
}
