use crate::raw::CameraProfile;
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
    /// The background the ends were placed around, per channel, in the input
    /// scale. Kept so the run can report where it lands once stretched,
    /// which is the number the neutrality of the file rests on.
    pub median: [f32; 3],
}

impl StretchParams {
    /// Where each channel's background lands in the output range, as a
    /// fraction: `(median − black) / (white − black)`. The three are equal by
    /// construction (see `analyze_stretch`), so this is the level the sky
    /// comes out at whatever colour the night was.
    pub fn background_level(&self) -> [f32; 3] {
        let mut out = [0.0f32; 3];
        for c in 0..3 {
            let span = self.white[c] - self.blacks[c];
            out[c] = if span > 1e-9 { (self.median[c] - self.blacks[c]) / span } else { 0.0 };
        }
        out
    }
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
/// Ceiling on the margin below the background, in units of the span between
/// the background and the white point. Reached only by an image with almost
/// no highlights above its own sky — there is no stretch worth writing for
/// that, and this stops the black point from running away below zero.
const BLACKS_MAX_SPAN: f32 = 1.0;

/// Analyses an interleaved RGB buffer and returns the stretch parameters.
///
/// The buffer is already white-balanced (`frame.wb`, baked in `stack::add`),
/// but a camera's white balance is a statement about daylight, not about the
/// light the sky was actually under. What the pipeline preserves on purpose —
/// the pedestal, the sky's own level — carries that light's colour, and a
/// per-channel stretch that treats each channel on its own terms leaves it
/// there: on the Canon test session the background came out at 7.6 % of the
/// range in red, 4.9 % in green and 10.2 % in blue, which is not a stretch
/// artefact but a cast, and one no global curve in a raw developer can undo
/// because it lives in the black points.
///
/// So the two ends are chosen to say the same thing in all three channels.
/// Each channel's white point still saturates the brightest 0.1 % of itself,
/// and each black point is placed so that the **background lands at the same
/// height in the output range for R, G and B**:
///
/// ```text
/// black_c = median_c − s · (white_c − median_c)
/// ```
///
/// with a single `s` shared by the three. That `s` is the largest any channel
/// needs to keep its `BLACKS_K_MAD` margin below its own background, so no
/// channel is clipped more tightly than it would have been on its own, and
/// the two that needed less simply keep more of their skirt. What comes out
/// is neutral by construction: the sky sits at one level, the stars keep the
/// colour they had relative to it, and the developer opens a file whose
/// channels are already aligned.
pub fn analyze_stretch(rgb: &[f32]) -> StretchParams {
    let mut median = [0.0f32; 3];
    let mut mad = [0.0f32; 3];
    let mut white = [0.0f32; 3];
    for c in 0..3 {
        let (m, d) = median_mad_channel(rgb, c);
        median[c] = m;
        mad[c] = d;
        white[c] = percentile_channel(rgb, c, WHITE_PERCENTILE);
    }
    let mut span = 0.0f32;
    for c in 0..3 {
        let above = white[c] - median[c];
        if above > 1e-9 {
            span = span.max(BLACKS_K_MAD * mad[c] / above);
        }
    }
    let span = span.min(BLACKS_MAX_SPAN);
    let mut blacks = [0.0f32; 3];
    for c in 0..3 {
        blacks[c] = (median[c] - span * (white[c] - median[c])).max(0.0);
    }
    StretchParams { white, blacks, median }
}

/// The camera's XYZ→camera matrix (D65) as the nine SRATIONALs the DNG
/// `ColorMatrix1` tag takes, in dcraw's ×10000 fixed point.
///
/// The values come from the body the frames were read from — `rawloader`
/// carries the official LibRaw/dcraw matrix of every camera it decodes — so
/// this is a unit conversion and nothing else. `raw::load` substitutes the
/// identity when a file offers no matrix at all, which produces flat but
/// valid colours.
fn color_matrix1(xyz_to_cam: [[f32; 3]; 3]) -> Vec<SRational> {
    let d = 10000;
    xyz_to_cam
        .iter()
        .flatten()
        .map(|v| SRational { n: (v * d as f32).round() as i32, d })
        .collect()
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
    camera: &CameraProfile,
    stretch: Option<StretchParams>,
) -> Result<()> {
    write_dng_impl(path, rgb, w, h, camera, stretch, true)
}

/// Same as `write_dng` but without printing the stretch parameters (for
/// batch exports where the stretch is shared and has already been reported).
pub fn write_dng_quiet(
    path: &Path,
    rgb: &[f32],
    w: usize,
    h: usize,
    camera: &CameraProfile,
    stretch: Option<StretchParams>,
) -> Result<()> {
    write_dng_impl(path, rgb, w, h, camera, stretch, false)
}

fn write_dng_impl(
    path: &Path,
    rgb: &[f32],
    w: usize,
    h: usize,
    camera: &CameraProfile,
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
    if verbose && let Some(st) = stretch {
        let pct = |c: usize| 100.0 * blacks[c] / white[c];
        let bg = st.background_level();
        say!(
            "stretch: white=R {:.6} G {:.6} B {:.6}  blacks=R {:.6} G {:.6} B {:.6} ({:.1}% / {:.1}% / {:.1}% of the range); background at {:.1}% / {:.1}% / {:.1}% — the three aligned",
            white[0], white[1], white[2],
            blacks[0], blacks[1], blacks[2],
            pct(0), pct(1), pct(2),
            100.0 * bg[0], 100.0 * bg[1], 100.0 * bg[2]
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
    de.write_tag(Tag::Unknown(50708), camera.model.as_str())?; // UniqueCameraModel

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
    let color_matrix = color_matrix1(camera.xyz_to_cam);
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

/// Copies EXIF + MakerNotes + XMP + ICC from the source RAW file to the
/// destination DNG using exiftool, excluding the structural and DNG tags we
/// have just written.
pub fn copy_metadata(source_raw: &Path, dest_dng: &Path) -> Result<()> {
    let out = Command::new("exiftool")
        .args([
            "-TagsFromFile",
        ])
        .arg(source_raw)
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
