//! Clean sequence export for timelapse (`--export-clean`).
//!
//! For every frame, in chronological order:
//! 1. Load + clean in sensor coordinates (`flatten::clean_frame` with the
//!    global model plus the frame's own temporal anomaly,
//!    `flatten::fit_frame_corr_ex`: horizon glow/twilight and halo variation
//!    with respect to the session median).
//! 2. **Stabilization**: resampled to the reference system with the same
//!    similarity transform (rotation + translation) used by the stacking,
//!    and cropped to the stack's common rectangle. The sky stays fixed;
//!    the guiding errors disappear.
//! 3. **Deflickering**: per channel, the frame's median (sky, foreground
//!    excluded) and high percentile (star cores, physically constant) are
//!    matched to those of the reference stack via a linear gain + offset.
//!    This way tone, brightness and contrast do not drift over the session,
//!    and the neighbouring frames end up at the same level before being
//!    combined.
//! 4. **Temporal noise reduction**: trimmed mean of a sliding window of
//!    `window` aligned frames centred on the frame (the maximum and the
//!    minimum per pixel are discarded → rejects satellite/plane trails and
//!    hot pixels). Noise ≈ /√(window−2). The **foreground** pixels (trees,
//!    horizon; per-cell mask from pass 1) do not enter the combination:
//!    where the current frame is occluded it is shown as is, and where the
//!    neighbours are, they are not used. The **transients** (meteors,
//!    satellites, planes) that the trimmed mean would erase are detected as
//!    an excess of the frame over the combination and are kept with the
//!    value from the original frame.
//! 5. A linear DNG is written with the same stretch as the stack.

use crate::say;
use crate::align::Similarity;
use crate::exif;
use crate::fixed;
use crate::flatten::{self, CellMask};
use crate::output::{self, StretchParams};
use crate::raw::{self, CameraProfile};
use crate::FrameInfo;
use anyhow::{anyhow, Context, Result};
use rayon::prelude::*;
use std::collections::{HashMap, VecDeque};
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Robust per-channel statistics for the deflickering.
#[derive(Clone, Copy, Debug)]
pub struct Stats {
    pub med: [f32; 3],
    pub hi: [f32; 3],
}

const HI_PERCENTILE: f32 = 99.9;

/// Per-channel sky level and star-core amplitude, over ~1M stratified
/// samples. `mask` (optional, one byte per pixel, ≠0 = foreground) excludes
/// those pixels: the tree/horizon is not sky.
///
/// Neither statistic may depend on how much of the frame a cloud covers,
/// because the deflickering divides one by the other and a light-polluted
/// cloud is the brightest thing in the frame:
///
/// * `med` is the `SKY_PERCENTILE`-th percentile of the channel after a
///   light blur, not its median. A cloud only adds light, so the level of
///   the sky is read from below it.
/// * `hi` is that level plus the high percentile of the **high-pass**
///   (channel − its own blur): the amplitude of a star core over the sky it
///   sits on. A cloud is smooth and contributes nothing to it, so what gets
///   measured is always the brightest star cores the frame still shows.
pub fn stats(rgb: &[f32], w: usize, h: usize, mask: Option<&[u8]>) -> Stats {
    let n_pix = w * h;
    let stride = (n_pix / 1_000_000).max(1);
    let mut s = Stats { med: [0.0; 3], hi: [0.0; 3] };
    let pick = |v: &mut Vec<f32>, pct: f32| -> f32 {
        let k = ((v.len() as f32 * pct / 100.0) as usize).min(v.len() - 1);
        let (_, x, _) = v.select_nth_unstable_by(k, |a: &f32, b: &f32| a.partial_cmp(b).unwrap());
        *x
    };
    for c in 0..3 {
        let ch: Vec<f32> = (0..n_pix).into_par_iter().map(|i| rgb[i * 3 + c]).collect();
        let b = box_blur(&ch, w, h, OCC_STAR_RADIUS);
        let keep = |i: &usize| mask.map_or(true, |m| m[*i] == 0);
        let mut lv: Vec<f32> = (0..n_pix).step_by(stride).filter(keep).map(|i| b[i]).filter(|x| x.is_finite()).collect();
        if lv.is_empty() {
            continue;
        }
        // Below the sky percentile there must be sky, not landscape: the
        // stack's reference statistics are taken without a mask, and a
        // horizon covering a fifth of the frame is exactly what the 20th
        // percentile would land on. Anything clearly darker than the median
        // is not sky (`SKY_DARK_MIN`, as in the transient detection) and is
        // dropped first.
        let med = pick(&mut lv, 50.0);
        let dark = SKY_DARK_MIN * med;
        let mut sky: Vec<f32> = lv.iter().copied().filter(|&v| v >= dark).collect();
        if sky.is_empty() {
            sky = lv;
        }
        s.med[c] = pick(&mut sky, SKY_PERCENTILE);
        let mut hv: Vec<f32> = (0..n_pix)
            .step_by(stride)
            .filter(keep)
            .filter(|&i| b[i] >= dark)
            .map(|i| ch[i] - b[i])
            .filter(|x| x.is_finite())
            .collect();
        if hv.is_empty() {
            continue;
        }
        s.hi[c] = s.med[c] + pick(&mut hv, HI_PERCENTILE).max(0.0);
    }
    s
}

/// Per-channel gain + offset to bring `cur` to `reference`.
/// `v' = (v − med_cur)·g + med_ref`, `g = (hi_ref − med_ref)/(hi_cur − med_cur)`:
/// the offset matches the level and the gain the star-core contrast.
///
/// With `guard` on the gain may compress the contrast but never expand it.
/// A frame whose star cores come out below the reference's has almost
/// always lost them to the atmosphere — cloud, haze, dew — and stretching
/// them back up to the reference is exactly what makes the stars of a
/// clouded frame burn brighter than the stars of a clear one, and what
/// blows the cloud they sit on past the white of the stretch along the way.
pub fn normalize(rgb: &mut [f32], cur: &Stats, reference: &Stats, guard: bool) -> [f32; 3] {
    let mut gain = [1.0f32; 3];
    for c in 0..3 {
        let num = reference.hi[c] - reference.med[c];
        let den = (cur.hi[c] - cur.med[c]).max(1e-6);
        let hi = if guard { 1.0 } else { 4.0 };
        gain[c] = (num / den).clamp(0.25, hi);
    }
    rgb.par_chunks_mut(3).for_each(|px| {
        for c in 0..3 {
            px[c] = ((px[c] - cur.med[c]) * gain[c] + reference.med[c]).max(0.0);
        }
    });
    gain
}

/// Resamples `rgb` (fw×fh, frame system) to the reference system cropped to
/// `(x0, y0, ow, oh)`: for every pixel p of the reference the frame is
/// sampled at `q = m(p)` (bilinear). Outside the frame → 0.
pub fn warp_to_ref(
    rgb: &[f32],
    fw: usize,
    fh: usize,
    m: &Similarity,
    x0: usize,
    y0: usize,
    ow: usize,
    oh: usize,
) -> Vec<f32> {
    let mut out = vec![0.0f32; ow * oh * 3];
    let fwf = fw as f32;
    let fhf = fh as f32;
    out.par_chunks_mut(ow * 3).enumerate().for_each(|(ry, row)| {
        let y = (y0 + ry) as f32;
        for rx in 0..ow {
            let x = (x0 + rx) as f32;
            let (qx, qy) = m.apply(x, y);
            if qx < 0.0 || qy < 0.0 || qx >= fwf - 1.0 || qy >= fhf - 1.0 {
                continue;
            }
            let ix = qx.floor() as usize;
            let iy = qy.floor() as usize;
            let dx = qx - ix as f32;
            let dy = qy - iy as f32;
            let i00 = (iy * fw + ix) * 3;
            let i10 = i00 + 3;
            let i01 = ((iy + 1) * fw + ix) * 3;
            let i11 = i01 + 3;
            let w00 = (1.0 - dx) * (1.0 - dy);
            let w10 = dx * (1.0 - dy);
            let w01 = (1.0 - dx) * dy;
            let w11 = dx * dy;
            let d = rx * 3;
            for c in 0..3 {
                row[d + c] = rgb[i00 + c] * w00 + rgb[i10 + c] * w10 + rgb[i01 + c] * w01 + rgb[i11 + c] * w11;
            }
        }
    });
    out
}

/// Foreground mask of the frame brought to the cropped reference system
/// (one byte per pixel, 1 = foreground; outside the frame → 1, no data).
/// Nearest neighbour over the cells.
pub fn warp_mask_to_ref(
    mask: &CellMask,
    fw: usize,
    fh: usize,
    m: &Similarity,
    x0: usize,
    y0: usize,
    ow: usize,
    oh: usize,
) -> Vec<u8> {
    let mut out = vec![0u8; ow * oh];
    let fwf = fw as f32;
    let fhf = fh as f32;
    out.par_chunks_mut(ow).enumerate().for_each(|(ry, row)| {
        let y = (y0 + ry) as f32;
        for rx in 0..ow {
            let x = (x0 + rx) as f32;
            let (qx, qy) = m.apply(x, y);
            row[rx] = if qx < 0.0 || qy < 0.0 || qx >= fwf - 1.0 || qy >= fhf - 1.0 || mask.at_px(qx, qy) { 1 } else { 0 };
        }
    });
    out
}

/// Per-pixel trimmed mean of several aligned frames of the same size: with
/// ≥ 4 samples the maximum and the minimum are discarded (trails, hot
/// pixels); with fewer, a plain mean. `masks[f]` (one byte per pixel,
/// ≠0 = foreground) excludes those samples; if the current frame (`cur`) is
/// occluded at the pixel, its own value is returned: the tree/horizon is
/// seen as is and is not smeared with the neighbouring frames, and the sky
/// that is occluded in other frames is not darkened.
pub fn combine_window(frames: &[&[f32]], masks: &[&[u8]], cur: usize) -> Vec<f32> {
    let n = frames.len();
    assert!(n >= 1 && masks.len() == n && cur < n);
    let len = frames[0].len();
    let mut out = vec![0.0f32; len];
    out.par_chunks_mut(3 * 1024).enumerate().for_each(|(chunk_i, chunk)| {
        let base = chunk_i * 1024;
        let mut buf = [0.0f32; 64];
        for (k, px) in chunk.chunks_mut(3).enumerate() {
            let i = base + k;
            if masks[cur][i] != 0 {
                px.copy_from_slice(&frames[cur][i * 3..i * 3 + 3]);
                continue;
            }
            for c in 0..3 {
                let mut m = 0usize;
                for f in 0..n.min(64) {
                    if masks[f][i] != 0 { continue; }
                    buf[m] = frames[f][i * 3 + c];
                    m += 1;
                }
                if m == 0 {
                    px[c] = frames[cur][i * 3 + c];
                    continue;
                }
                let vals = &buf[..m];
                if m >= 4 {
                    let mut mn = f32::MAX;
                    let mut mx = f32::MIN;
                    let mut sum = 0.0f32;
                    for &v in vals.iter() {
                        mn = mn.min(v);
                        mx = mx.max(v);
                        sum += v;
                    }
                    px[c] = (sum - mn - mx) / (m - 2) as f32;
                } else {
                    px[c] = vals.iter().sum::<f32>() / m as f32;
                }
            }
        }
    });
    out
}

/// Same trimmed mean as `combine_window`, except that the neighbours are
/// **not** in the same system as the current frame: neighbour `f` is
/// sampled at `ts[f](p)` (bilinear) instead of at `p`.
///
/// This is what makes temporal noise reduction possible on an untracked
/// sequence. There the camera is fixed and the sky rotates, so combining
/// the window as it stands would drag every star into a short trail and
/// dim it. `ts[f]` takes the current frame's sky to the neighbour's, so
/// the stars line up while the landscape deliberately does not — and the
/// landscape is precisely what the mask makes the current frame supply on
/// its own. A neighbour's own foreground, once displaced by `ts[f]`, lands
/// on sky it does not belong to; the same mask lookup (in the neighbour's
/// coordinates, where its landscape still is) throws those samples away,
/// and so does a sample falling outside the frame.
pub fn combine_window_warped(
    frames: &[&[f32]],
    masks: &[&[u8]],
    ts: &[Similarity],
    cur: usize,
    w: usize,
    h: usize,
) -> (Vec<f32>, Vec<u8>) {
    let n = frames.len();
    assert!(n >= 1 && masks.len() == n && ts.len() == n && cur < n);
    let (wf, hf) = (w as f32, h as f32);
    let mut out = vec![0.0f32; w * h * 3];
    let mut used = vec![0u8; w * h];
    out.par_chunks_mut(w * 3).zip(used.par_chunks_mut(w)).enumerate().for_each(|(y, (row, urow))| {
        let yf = y as f32;
        let mut buf = [[0.0f32; 64]; 3];
        for x in 0..w {
            let i = y * w + x;
            if masks[cur][i] != 0 {
                row[x * 3..x * 3 + 3].copy_from_slice(&frames[cur][i * 3..i * 3 + 3]);
                urow[x] = 1;
                continue;
            }
            let xf = x as f32;
            let mut m = 0usize;
            for f in 0..n.min(64) {
                if f == cur {
                    for c in 0..3 {
                        buf[c][m] = frames[cur][i * 3 + c];
                    }
                    m += 1;
                    continue;
                }
                let (qx, qy) = ts[f].apply(xf, yf);
                if qx < 0.0 || qy < 0.0 || qx >= wf - 1.0 || qy >= hf - 1.0 {
                    continue;
                }
                let ix = qx as usize;
                let iy = qy as usize;
                if masks[f][iy * w + ix] != 0 {
                    continue;
                }
                let (dx, dy) = (qx - ix as f32, qy - iy as f32);
                let i00 = (iy * w + ix) * 3;
                let i10 = i00 + 3;
                let i01 = ((iy + 1) * w + ix) * 3;
                let i11 = i01 + 3;
                let w00 = (1.0 - dx) * (1.0 - dy);
                let w10 = dx * (1.0 - dy);
                let w01 = (1.0 - dx) * dy;
                let w11 = dx * dy;
                let src = frames[f];
                for c in 0..3 {
                    buf[c][m] = src[i00 + c] * w00 + src[i10 + c] * w10 + src[i01 + c] * w01 + src[i11 + c] * w11;
                }
                m += 1;
            }
            for c in 0..3 {
                let vals = &buf[c][..m];
                row[x * 3 + c] = if m >= 4 {
                    let mut mn = f32::MAX;
                    let mut mx = f32::MIN;
                    let mut sum = 0.0f32;
                    for &v in vals {
                        mn = mn.min(v);
                        mx = mx.max(v);
                        sum += v;
                    }
                    (sum - mn - mx) / (m - 2) as f32
                } else if m > 0 {
                    vals.iter().sum::<f32>() / m as f32
                } else {
                    frames[cur][i * 3 + c]
                };
            }
            urow[x] = m.min(255) as u8;
        }
    });
    (out, used)
}

/// Ramp (frame sky level / session median) between the night-sky cleaning
/// and the natural dawn/twilight version.
const DAWN_R0: f32 = 1.5;
const DAWN_R1: f32 = 2.3;

/// Threshold (in σ of the noise of the smoothed difference) to consider a
/// pixel a transient, and dilation radius of the mask.
const TRANSIENT_K: f32 = 4.0;
const TRANSIENT_DILATE: usize = 2;
/// Minimum extent (px, longer side of the box) of a connected component
/// for it to be considered a real trail and not a star residue or noise.
const TRANSIENT_MIN_EXTENT: usize = 40;
/// Minimum ratio between a component's long and short axes, measured from
/// its second moments. What the extent alone cannot tell apart: a meteor or a
/// satellite is long and thin, while a patch of cloud or a blob left by the
/// high-pass filter is about as wide as it is long, and taking those from the
/// frame untouched is what put square patches into the export.
const TRANSIENT_MIN_ELONGATION: f32 = 3.0;
/// Fewest samples a pixel's temporal combination must have had for a
/// transient to be worth rescuing there. Below the trimmed mean's own
/// threshold of four the combination is a plain mean of two or three frames,
/// noisy enough to be taken for a transient wherever the window loses
/// coverage — which is the whole border of an untracked export.
const TRANSIENT_MIN_SAMPLES: usize = 4;
/// Radius (px) of the high-pass applied to the difference before the
/// threshold: a trail is narrow; a smooth and extensive change of the sky
/// between frames (twilight, horizon glow, a faint cloud passing) is not a
/// transient and must not cancel the noise reduction over half the image.
const TRANSIENT_HP_RADIUS: usize = 32;
/// Fraction of the sky level below which a pixel is not considered sky for
/// transient purposes (foreground and its penumbra).
const TRANSIENT_SKY_MIN: f32 = 0.6;

/// Radius (px) of the high-pass that separates the stars from the sky they
/// sit on, and radius over which its modulus is averaged to measure how
/// much star signal a region carries.
const OCC_STAR_RADIUS: usize = 8;
const OCC_STAR_AVG: usize = 48;
/// Sigmas of the high-pass noise subtracted before a region's star signal
/// is measured.
const OCC_NOISE_K: f32 = 3.0;
/// Radius (px) the occlusion weight is smoothed over.
const OCC_SMOOTH: usize = 128;
/// Factor the occlusion measurement is downscaled by. Every radius above is
/// given in full-size pixels and divided by it.
const OCC_SCALE: usize = 4;
/// Ratio between the star signal of the frame and that of its temporal
/// combination below which the frame counts as occluded there — a cloud
/// passing in front of that patch of sky. The ramp runs from `OCC_R1` (the
/// combination is left alone) to `OCC_R0` (the frame is shown as it is), so
/// that the edge of a cloud does not become an edge in the image.
const OCC_R0: f32 = 0.40;
const OCC_R1: f32 = 0.75;
/// Percentile of the (lightly blurred) channel taken as the frame's sky
/// level. Not the median: a cloud only ever adds light, so the level of the
/// sky is read below it, low enough that a frame the cloud covers for the
/// most part still measures its clear sky and not its cloud.
const SKY_PERCENTILE: f32 = 20.0;
/// Fraction of the median below which a pixel is not sky (foreground and
/// its penumbra), as in `TRANSIENT_SKY_MIN`.
const SKY_DARK_MIN: f32 = 0.6;

/// Separable box moving average (radius `r`, clipped border) — O(n).
pub(crate) fn box_blur(src: &[f32], w: usize, h: usize, r: usize) -> Vec<f32> {
    let mut tmp = vec![0.0f32; w * h];
    tmp.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        let s = &src[y * w..(y + 1) * w];
        let mut acc = 0.0f32;
        let mut cnt = 0usize;
        for x in 0..r.min(w) { acc += s[x]; cnt += 1; }
        for x in 0..w {
            if x + r < w { acc += s[x + r]; cnt += 1; }
            if x > r { acc -= s[x - r - 1]; cnt -= 1; }
            row[x] = acc / cnt as f32;
        }
    });
    let mut out = vec![0.0f32; w * h];
    // Columns: parallelize over blocks of columns.
    let cols: Vec<Vec<f32>> = (0..w)
        .into_par_iter()
        .map(|x| {
            let mut col = vec![0.0f32; h];
            let mut acc = 0.0f32;
            let mut cnt = 0usize;
            for y in 0..r.min(h) { acc += tmp[y * w + x]; cnt += 1; }
            for y in 0..h {
                if y + r < h { acc += tmp[(y + r) * w + x]; cnt += 1; }
                if y > r { acc -= tmp[(y - r - 1) * w + x]; cnt -= 1; }
                col[y] = acc / cnt as f32;
            }
            col
        })
        .collect();
    out.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        for x in 0..w { row[x] = cols[x][y]; }
    });
    out
}

/// Three box passes ≈ a Gaussian. The shape of the kernel reaches the
/// image: a single box is square, and a square response around every bright
/// star tiles the export with square patches.
pub(crate) fn blur_gauss(src: &[f32], w: usize, h: usize, r: usize) -> Vec<f32> {
    let mut b = box_blur(src, w, h, r);
    b = box_blur(&b, w, h, (r / 2).max(1));
    box_blur(&b, w, h, (r / 2).max(1))
}

/// Undoes the temporal combination wherever the cloud covers the frame's
/// sky, and returns how much of the frame that was.
///
/// The window is combined **on the sky**: the neighbours are warped so that
/// their stars fall on the current frame's stars. A cloud does not follow
/// the sky — it crosses it — so a star this frame has behind a cloud is
/// clear in the neighbours, and the trimmed mean hands it back at its full
/// brightness over a cloud that has meanwhile been averaged smooth. That is
/// what puts a field of hard white stars on top of a cloud that in the
/// frame hides them, and it is not what the sequence looked like.
///
/// Occlusion is measured as the loss of star signal, not as brightness: the
/// modulus of the high-pass with its noise floor subtracted, averaged over
/// `OCC_STAR_AVG`, in the frame and in the combination. Where the frame keeps as much as the
/// combination the two agree; where a cloud has taken the stars away the
/// frame's is only its noise and the ratio collapses. The frame is faded
/// back in over the ramp `OCC_R1` → `OCC_R0` so the edge of a cloud does not
/// come out as an edge in the image. Brightness is deliberately not part of
/// the test: the Milky Way is broad and bright too, and it is full of stars.
pub fn preserve_occlusions(cur: &[f32], combined: &mut [f32], fg: &[u8], w: usize, h: usize) -> f32 {
    // The whole measurement runs on a quarter-scale luminance. What comes
    // out of it is a weight smoothed over a hundred-odd pixels, so the
    // detail thrown away is detail it would have blurred away anyway — and
    // at full size the dozen box passes below cost more than the temporal
    // combination they correct.
    let (sw, sh) = ((w / OCC_SCALE).max(1), (h / OCC_SCALE).max(1));
    let sn = sw * sh;
    let small_lum = |img: &[f32]| -> Vec<f32> {
        let mut out = vec![0.0f32; sn];
        out.par_chunks_mut(sw).enumerate().for_each(|(y, row)| {
            for x in 0..sw {
                let mut acc = 0.0f32;
                let mut cnt = 0.0f32;
                for yy in y * OCC_SCALE..((y + 1) * OCC_SCALE).min(h) {
                    for xx in x * OCC_SCALE..((x + 1) * OCC_SCALE).min(w) {
                        let i = (yy * w + xx) * 3;
                        acc += img[i] + img[i + 1] + img[i + 2];
                        cnt += 3.0;
                    }
                }
                row[x] = if cnt > 0.0 { acc / cnt } else { 0.0 };
            }
        });
        out
    };
    // Star signal: the high-pass with its noise floor taken off, averaged
    // over a neighbourhood. The floor matters — the modulus of the
    // high-pass of a 25 s frame is mostly noise, and a measure that noise
    // dominates says the same thing under a cloud as under a clear sky.
    // What is left after subtracting OCC_NOISE_K sigmas is the stars.
    // The floor is the **frame's**, and it is the one used on both sides.
    // Each image measured against its own would compare unequal things: the
    // combination has already had its noise divided by the window, so a
    // floor of its own would leave it more signal than the frame keeps
    // wherever the sky is empty, and that difference alone — not a cloud —
    // would undo the noise reduction over a good part of a clear frame.
    let hp = |lum: &[f32]| -> Vec<f32> {
        let b = blur_gauss(lum, sw, sh, (OCC_STAR_RADIUS / OCC_SCALE).max(1));
        lum.par_iter().zip(b.par_iter()).map(|(l, s)| (l - s).abs()).collect()
    };
    let hc = hp(&small_lum(cur));
    let hk = hp(&small_lum(combined));
    let sigma = {
        let stride = (sn / 1_000_000).max(1);
        let mut sample: Vec<f32> = hc.iter().step_by(stride).copied().collect();
        let k = sample.len() / 2;
        let (_, med, _) = sample.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
        // median(|x|) = 0.6745 sigma for a centred normal.
        *med * 1.4826
    };
    let star = |h: &[f32]| -> Vec<f32> {
        let core: Vec<f32> = h.par_iter().map(|v| (v - OCC_NOISE_K * sigma).max(0.0)).collect();
        blur_gauss(&core, sw, sh, (OCC_STAR_AVG / OCC_SCALE).max(1))
    };
    let sc = star(&hc);
    let sk = star(&hk);
    // Scale below which there is no star signal to speak of on either side:
    // the ratio there must come out as 1 (nothing to bring back), not as
    // the quotient of two noises.
    let eps = 0.5 * (sk.par_iter().map(|&v| v as f64).sum::<f64>() / sn as f64) as f32;
    let raw: Vec<f32> = (0..sn)
        .into_par_iter()
        .map(|i| {
            let r = (sc[i] + eps) / (sk[i] + eps).max(1e-9);
            let t = ((OCC_R1 - r) / (OCC_R1 - OCC_R0)).clamp(0.0, 1.0);
            t * t * (3.0 - 2.0 * t)
        })
        .collect();
    // The weight is a map of how much cloud lies over the frame, not a
    // per-star verdict: smoothed over `OCC_SMOOTH` it follows the body of a
    // cloud and dilutes the disc a single dimmed star would otherwise stamp
    // on the export.
    //
    // The foreground is left out of that average, and the average is
    // normalized by what did take part. A tree hides the sky exactly as a
    // cloud does and reads the same here, but it is the mask that already
    // gives those pixels to the frame; letting them into the blur would
    // spread the tree's verdict over the sky around it and undo the noise
    // reduction along the whole horizon.
    let valid: Vec<f32> = (0..sn)
        .into_par_iter()
        .map(|i| {
            let (sx, sy) = (i % sw, i / sw);
            let mut any_fg = false;
            for yy in sy * OCC_SCALE..((sy + 1) * OCC_SCALE).min(h) {
                for xx in sx * OCC_SCALE..((sx + 1) * OCC_SCALE).min(w) {
                    if fg[yy * w + xx] != 0 {
                        any_fg = true;
                    }
                }
            }
            if any_fg { 0.0 } else { 1.0 }
        })
        .collect();
    let ws = {
        let r = (OCC_SMOOTH / OCC_SCALE).max(1);
        let num = blur_gauss(&raw.par_iter().zip(valid.par_iter()).map(|(a, b)| a * b).collect::<Vec<f32>>(), sw, sh, r);
        let den = blur_gauss(&valid, sw, sh, r);
        num.par_iter()
            .zip(den.par_iter())
            .zip(valid.par_iter())
            .map(|((n, d), v)| if *v == 0.0 || *d <= 1e-6 { 0.0 } else { n / d })
            .collect::<Vec<f32>>()
    };
    let area = ws.par_iter().map(|&v| v as f64).sum::<f64>() / sn as f64;
    if let Some(dir) = std::env::var_os("APILAAA_DEBUG_DIR") {
        // Occlusion weight (8-bit PGM at the measurement's own scale).
        static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let k = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let mut buf = format!("P5\n{} {}\n255\n", sw, sh).into_bytes();
        buf.extend(ws.iter().map(|&v| (v.clamp(0.0, 1.0) * 255.0) as u8));
        let _ = std::fs::write(std::path::Path::new(&dir).join(format!("occluded_{k:04}.pgm")), buf);
    }
    // Back to full size, bilinear on the small grid (node = block centre).
    let sample_w = |x: usize, y: usize| -> f32 {
        let fx = (x as f32 + 0.5) / OCC_SCALE as f32 - 0.5;
        let fy = (y as f32 + 0.5) / OCC_SCALE as f32 - 0.5;
        let x0 = fx.floor().clamp(0.0, (sw - 1) as f32) as usize;
        let y0 = fy.floor().clamp(0.0, (sh - 1) as f32) as usize;
        let x1 = (x0 + 1).min(sw - 1);
        let y1 = (y0 + 1).min(sh - 1);
        let dx = (fx - x0 as f32).clamp(0.0, 1.0);
        let dy = (fy - y0 as f32).clamp(0.0, 1.0);
        ws[y0 * sw + x0] * (1.0 - dx) * (1.0 - dy)
            + ws[y0 * sw + x1] * dx * (1.0 - dy)
            + ws[y1 * sw + x0] * (1.0 - dx) * dy
            + ws[y1 * sw + x1] * dx * dy
    };
    combined.par_chunks_mut(w * 3).enumerate().for_each(|(y, row)| {
        for x in 0..w {
            let wgt = sample_w(x, y);
            if wgt <= 0.0 {
                continue;
            }
            let i = (y * w + x) * 3;
            for c in 0..3 {
                row[x * 3 + c] = row[x * 3 + c] * (1.0 - wgt) + cur[i + c] * wgt;
            }
        }
    });
    area as f32
}

/// Preserves the transient objects of the current frame (meteors,
/// satellites, planes, anything that is only in this frame) over the
/// combined image of the window: the trimmed mean erases them, so they are
/// detected as a positive excess of the frame over the combination
/// (luminance difference smoothed 3×3 and high-passed with radius
/// `TRANSIENT_HP_RADIUS` > `TRANSIENT_K`·σ, σ robust via MAD; mask dilated
/// `TRANSIENT_DILATE` px so the edges of the trail are not cut off) and at
/// those pixels the value from the original frame is used as is. Stars and
/// hot pixels are static → difference ≈ 0 → they are left alone. Returns
/// the number of preserved pixels.
pub fn preserve_transients(
    cur: &[f32],
    combined: &mut [f32],
    fg: &[u8],
    used: Option<&[u8]>,
    w: usize,
    h: usize,
) -> usize {
    let n = w * h;
    // Luminance difference and luminance of the combination (sky).
    let mut d = vec![0.0f32; n];
    let mut lum = vec![0.0f32; n];
    d.par_chunks_mut(w).zip(lum.par_chunks_mut(w)).enumerate().for_each(|(y, (row, lrow))| {
        for x in 0..w {
            let i = (y * w + x) * 3;
            let lc = (combined[i] + combined[i + 1] + combined[i + 2]) / 3.0;
            lrow[x] = lc;
            row[x] = (cur[i] + cur[i + 1] + cur[i + 2]) / 3.0 - lc;
        }
    });
    // Sky level (median of the combination outside the foreground): pixels
    // clearly darker than the sky (body and out-of-focus penumbra of
    // trees/horizon) are not sky and are not evaluated — there the
    // deflickering (gain around the median) moves the values from one frame
    // to another and there are no trails to preserve.
    let sky_med = {
        let stride = (n / 1_000_000).max(1);
        let mut v: Vec<f32> = (0..n).step_by(stride).filter(|&i| fg[i] == 0).map(|i| lum[i]).collect();
        if v.is_empty() { 0.0 } else {
            let k = v.len() / 2;
            let (_, m, _) = v.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
            *m
        }
    };
    let dark_thr = TRANSIENT_SKY_MIN * sky_med;
    // Separable 3×3 (mean) smoothing.
    let mut tmp = vec![0.0f32; n];
    tmp.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        let src = &d[y * w..(y + 1) * w];
        for x in 0..w {
            let l = if x > 0 { src[x - 1] } else { src[x] };
            let r = if x + 1 < w { src[x + 1] } else { src[x] };
            row[x] = (l + src[x] + r) / 3.0;
        }
    });
    let mut ds = vec![0.0f32; n];
    ds.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        let y0 = y.saturating_sub(1);
        let y1 = (y + 1).min(h - 1);
        for x in 0..w {
            row[x] = (tmp[y0 * w + x] + tmp[y * w + x] + tmp[y1 * w + x]) / 3.0;
        }
    });
    drop(tmp);
    // High-pass: drop the smooth and extensive component of the difference.
    //
    // Three box passes rather than one. A single box is square, so around a
    // strong local excess its mean rises over the whole 2r+1 square and the
    // difference goes negative across all of it — a square plateau that the
    // deficit tail below then takes for a dark moving occluder, and that
    // comes out of the export as a square patch of untouched frame. Three
    // passes approximate a Gaussian, whose response is round and decays.
    let lp = {
        let mut b = box_blur(&ds, w, h, TRANSIENT_HP_RADIUS);
        b = box_blur(&b, w, h, TRANSIENT_HP_RADIUS / 2);
        box_blur(&b, w, h, TRANSIENT_HP_RADIUS / 2)
    };
    ds.par_iter_mut().zip(lp.par_iter()).for_each(|(d, l)| *d -= *l);
    drop(lp);
    // Robust σ via MAD over a stratified sample.
    let stride = (n / 1_000_000).max(1);
    let mut sample: Vec<f32> = ds.iter().step_by(stride).copied().collect();
    let k = sample.len() / 2;
    let (_, med, _) = sample.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
    let med = *med;
    let mut abs: Vec<f32> = sample.iter().map(|v| (v - med).abs()).collect();
    let (_, mad, _) = abs.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
    let sigma = (*mad * 1.4826).max(1e-9);
    let thr_hi = med + TRANSIENT_K * sigma;
    let thr_lo = med - TRANSIENT_K * sigma;
    // Initial mask and connected-component filter: only the components
    // spanning ≥ TRANSIENT_MIN_EXTENT px survive (a trail in a long
    // exposure is always elongated; a star residue or noise is not).
    // Two tails: excess (bright trails) and deficit (dark moving occluders:
    // foreground, clouds) — in both cases the frame wins.
    // Only where the trimmed mean actually ran. Towards the edges of an
    // untracked export the neighbouring frames warp off the sensor, so a
    // pixel there is left with two or three samples and its combination is
    // both noisy and far from the frame — which reads as a transient, and
    // turned whole strips along the border into raw frame pasted over the
    // export. There is nothing to preserve where nothing was averaged away.
    let mut mask: Vec<u8> = (0..n)
        .map(|i| {
            let enough = used.map_or(true, |u| u[i] as usize >= TRANSIENT_MIN_SAMPLES);
            let v = ds[i];
            (enough && (v > thr_hi || v < thr_lo) && lum[i] >= dark_thr) as u8
        })
        .collect();
    drop(lum);
    {
        let mut visited = vec![false; n];
        let mut stack: Vec<usize> = Vec::new();
        let mut comp: Vec<usize> = Vec::new();
        for start in 0..n {
            if mask[start] == 0 || visited[start] { continue; }
            comp.clear();
            stack.push(start);
            visited[start] = true;
            let (mut minx, mut maxx, mut miny, mut maxy) = (w, 0usize, h, 0usize);
            while let Some(i) = stack.pop() {
                comp.push(i);
                let x = i % w;
                let y = i / w;
                minx = minx.min(x); maxx = maxx.max(x); miny = miny.min(y); maxy = maxy.max(y);
                let nb = [
                    if x > 0 { Some(i - 1) } else { None },
                    if x + 1 < w { Some(i + 1) } else { None },
                    if y > 0 { Some(i - w) } else { None },
                    if y + 1 < h { Some(i + w) } else { None },
                ];
                for j in nb.into_iter().flatten() {
                    if mask[j] != 0 && !visited[j] { visited[j] = true; stack.push(j); }
                }
            }
            // Elongation from the second moments, so it does not depend on
            // which way the trail runs: the bounding box alone calls a
            // diagonal trail square and a square plateau elongated. A trail
            // in a long exposure is long and thin whichever way it points; a
            // blob left by the filter, or a patch of cloud, is not.
            let extent = (maxx - minx + 1).max(maxy - miny + 1);
            let elong = {
                let cnt = comp.len() as f64;
                let (mut sx, mut sy) = (0.0f64, 0.0f64);
                for &i in &comp { sx += (i % w) as f64; sy += (i / w) as f64; }
                let (mx, my) = (sx / cnt, sy / cnt);
                let (mut cxx, mut cyy, mut cxy) = (0.0f64, 0.0f64, 0.0f64);
                for &i in &comp {
                    let dx = (i % w) as f64 - mx;
                    let dy = (i / w) as f64 - my;
                    cxx += dx * dx; cyy += dy * dy; cxy += dx * dy;
                }
                cxx /= cnt; cyy /= cnt; cxy /= cnt;
                let tr = cxx + cyy;
                let det = cxx * cyy - cxy * cxy;
                let disc = (0.25 * tr * tr - det).max(0.0).sqrt();
                let l1 = 0.5 * tr + disc;
                let l2 = (0.5 * tr - disc).max(1e-9);
                (l1 / l2).sqrt() as f32
            };
            if extent < TRANSIENT_MIN_EXTENT || elong < TRANSIENT_MIN_ELONGATION {
                for &i in &comp { mask[i] = 0; }
            }
        }
    }
    let r = TRANSIENT_DILATE;
    if r > 0 {
        let mut t2 = vec![0u8; n];
        t2.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
            let src = &mask[y * w..(y + 1) * w];
            for x in 0..w {
                let lo = x.saturating_sub(r);
                let hi = (x + r).min(w - 1);
                row[x] = src[lo..=hi].iter().copied().max().unwrap_or(0);
            }
        });
        let mut m2 = vec![0u8; n];
        m2.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
            let lo = y.saturating_sub(r);
            let hi = (y + r).min(h - 1);
            for x in 0..w {
                let mut m = 0u8;
                for yy in lo..=hi {
                    m = m.max(t2[yy * w + x]);
                }
                row[x] = m;
            }
        });
        mask = m2;
    }
    let kept: usize = mask.iter().map(|&m| m as usize).sum();
    if let Some(dir) = std::env::var_os("APILAAA_DEBUG_DIR") {
        // Transient mask (binary PGM, 1/8 resolution) for inspection.
        static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let k = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let (sw, sh) = (w / 8, h / 8);
        let mut buf = format!("P5\n{} {}\n255\n", sw, sh).into_bytes();
        for y in 0..sh { for x in 0..sw {
            let mut any = 0u8;
            for yy in y * 8..y * 8 + 8 { for xx in x * 8..x * 8 + 8 { any = any.max(mask[yy * w + xx]); } }
            buf.push(if any != 0 { 255 } else { 0 });
        } }
        let _ = std::fs::write(std::path::Path::new(&dir).join(format!("transients_{k:04}.pgm")), buf);
    }
    combined.par_chunks_mut(w * 3).enumerate().for_each(|(y, row)| {
        for x in 0..w {
            if mask[y * w + x] != 0 {
                let i = (y * w + x) * 3;
                row[x * 3] = cur[i];
                row[x * 3 + 1] = cur[i + 1];
                row[x * 3 + 2] = cur[i + 2];
            }
        }
    });
    kept
}

pub struct ExportOpts<'a> {
    pub dir: &'a Path,
    pub window: usize,
    pub stabilize: bool,
    pub deflicker: bool,
    pub keep_transients: bool,
    pub camera: &'a CameraProfile,
    pub stretch: Option<StretchParams>,
    pub n_workers: usize,
    /// Contrast compensation β (see `flatten::clean_frame`).
    pub scatter_comp: f32,
    /// Untracked sequence (fixed tripod): the links of the sky chain
    /// (`links[k]` = frame k → frame k+1), so the temporal window is
    /// combined on the sky rather than on the sensor. None = tracked
    /// sequence, where the frames are already aligned on the sky and the
    /// window combines them as they stand.
    pub sky_links: Option<&'a [Option<Similarity>]>,
    /// Freedom given to the per-frame anomaly (see `flatten::AnomalyMode`).
    pub anomaly: flatten::AnomalyMode,
    /// Treat the sky the cloud covers as the frame's own (see
    /// `cloud_mask`): keep it out of the temporal window and out of the
    /// frame's level statistics.
    pub cloud_guard: bool,
}

/// Exports the sequence. `infos` = aligned frames (similarity ref→cur,
/// background map and foreground mask from pass 1; all of them are
/// exported, including the ones excluded from the stack); the unaligned
/// ones are skipped with a warning if stabilization is on. `crop` =
/// (x0, y0, w, h) of the stack. `bg_med` = median map the `model` was
/// fitted on (reference for the per-frame anomaly). `reference` = stats of
/// the final stack for the deflickering.
pub fn export_sequence(
    paths: &[PathBuf],
    infos: &[FrameInfo],
    model: &flatten::GlareModel,
    bg_med: &flatten::BgMap,
    crop: (usize, usize, usize, usize),
    reference: Stats,
    opts: &ExportOpts,
) -> Result<()> {
    std::fs::create_dir_all(opts.dir).with_context(|| format!("creating {}", opts.dir.display()))?;
    let (x0, y0, ow, oh) = if opts.stabilize { crop } else { (0, 0, model.width, model.height) };
    let window = opts.window.max(1) | 1; // odd
    let half = window / 2;
    let by_idx: HashMap<usize, &FrameInfo> = infos.iter().map(|f| (f.idx, f)).collect();
    // ref → frame k: the geometry the export warped with (a star-based
    // similarity when tracked, the tripod drift when not). Identity when
    // exporting in sensor coordinates, where nothing was warped.
    let lref = |idx: usize| -> Similarity {
        if opts.stabilize {
            by_idx.get(&idx).map(|f| f.m).unwrap_or_else(Similarity::identity)
        } else {
            Similarity::identity()
        }
    };
    let crop_in = Similarity::translation(x0 as f32, y0 as f32);
    let crop_out = Similarity::translation(-(x0 as f32), -(y0 as f32));

    // List of exportable frames in chronological order (the aligned ones;
    // without stabilization all of them are exported and the unaligned ones
    // go without a mask). End of the sequence: once the frame goes in its
    // natural version (full dawn, ratio ≥ DAWN_R1) and its sky is above the
    // stretch white, the DNG would come out entirely white.
    let white_g = opts.stretch.map_or(f32::INFINITY, |s| s.white[1]);
    let med_level = {
        let mut v: Vec<f32> = infos.iter().filter(|f| f.aligned && f.level > 0.0).map(|f| f.level).collect();
        if v.is_empty() { 0.0 } else {
            let k = v.len() / 2;
            let (_, m, _) = v.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
            *m
        }
    };
    let order: Vec<usize> = (0..paths.len())
        .filter(|i| match by_idx.get(i) {
            Some(f) => f.export && !(f.level >= white_g && med_level > 0.0 && f.level / med_level >= DAWN_R1),
            None => !opts.stabilize,
        })
        .collect();
    let skipped = paths.len() - order.len();
    say!(
        "exporting {} frames to {} (window {} frames, stabilize={}, deflicker={}, {}×{}){}",
        order.len(),
        opts.dir.display(),
        window,
        opts.stabilize,
        opts.deflicker,
        ow,
        oh,
        if skipped > 0 { format!("; {skipped} skipped (unaligned or sky saturated/above the white)") } else { String::new() }
    );
    let t0 = Instant::now();
    crate::ui::task_begin("export", "export", order.len() as u64);

    // Producer: load + clean + warp in parallel batches; consumer: sliding
    // window in order.
    // (position in `order`, image, foreground mask, deflicker text, natural
    // dawn version with its weight)
    type Entry = (usize, Vec<f32>, Vec<u8>, String, Option<(Vec<f32>, f32)>);
    let mut buffer: VecDeque<Entry> = VecDeque::new();
    let mut next_load = 0usize;
    let mut next_emit = 0usize;
    let mut meta_failed = false;
    let mut exif_failed = false;
    let batch = opts.n_workers.max(1);

    while next_emit < order.len() {
        if crate::ui::aborted() {
            return Err(anyhow!("stopped by the user"));
        }
        // Load until the window of the next frame to emit is covered.
        while next_load < order.len() && next_load <= next_emit + half + batch {
            let end = (next_load + batch).min(order.len());
            let loaded: Vec<Result<Entry>> = (next_load..end)
                .into_par_iter()
                .map(|pos| -> Result<Entry> {
                    let idx = order[pos];
                    let p = &paths[idx];
                    let frame = raw::load(p).with_context(|| format!("loading {}", p.display()))?;
                    if frame.width != model.width || frame.height != model.height {
                        return Err(anyhow!("{}: different dimensions", p.display()));
                    }
                    let info = by_idx.get(&idx).copied();
                    // Map and mask: the ones from pass 1 if the frame is
                    // aligned; otherwise (only without stabilization), they
                    // are computed here.
                    let (bg, mask) = match info {
                        Some(f) => (None, f.mask.clone()),
                        None => {
                            let bg = flatten::block_background(&frame);
                            let mask = flatten::foreground_mask(&bg);
                            (Some(bg), mask)
                        }
                    };
                    let fc = match (&bg, info) {
                        (Some(bg), _) => flatten::fit_frame_corr_ex(bg, bg_med, model, Some(&mask), opts.anomaly),
                        (None, Some(f)) => flatten::fit_frame_corr_ex(&f.bg, bg_med, model, Some(&mask), opts.anomaly),
                        (None, None) => unreachable!(),
                    };
                    // Full dawn/twilight: when the frame's sky exceeds
                    // DAWN_R0 times the median, the "night sky" cleaning
                    // (anomaly + deflicker) stops making sense (the real sky
                    // is a different one and it saturates at the horizon).
                    // It transitions continuously (ramp DAWN_R0→DAWN_R1) to
                    // the **natural** version: the defect model only, with no
                    // anomaly and no deflicker, and the sky is allowed to
                    // brighten up to white.
                    let w_dawn = {
                        let t = ((fc.level_ratio - DAWN_R0) / (DAWN_R1 - DAWN_R0)).clamp(0.0, 1.0);
                        // The natural version is ~2× brighter: so that the
                        // transition looks gradual, the weight starts off
                        // slowly (w²).
                        let w = (t * t * (3.0 - 2.0 * t)).powi(2);
                        if w < 0.02 { 0.0 } else { w }
                    };
                    let to_ref = |img: Vec<f32>| -> Vec<f32> {
                        if opts.stabilize {
                            let m = &info.unwrap().m;
                            warp_to_ref(&img, frame.width, frame.height, m, x0, y0, ow, oh)
                        } else {
                            img
                        }
                    };
                    let pmask = if opts.stabilize {
                        let m = &info.unwrap().m;
                        warp_mask_to_ref(&mask, frame.width, frame.height, m, x0, y0, ow, oh)
                    } else {
                        let mut pm = vec![0u8; ow * oh];
                        pm.par_chunks_mut(ow).enumerate().for_each(|(y, row)| {
                            for x in 0..ow {
                                row[x] = mask.at_px(x as f32 + 0.5, y as f32 + 0.5) as u8;
                            }
                        });
                        pm
                    };
                    let mut gain_txt = format!("  defects ×{:.2}/{:.2}/{:.2}  dome ×{:.2}/{:.2}/{:.2} (lv {:.2})  anomaly {:.1}%", fc.gain[0], fc.gain[1], fc.gain[2], fc.dome[0], fc.dome[1], fc.dome[2], fc.level_ratio, 100.0 * fc.range[1] / reference.med[1].max(1e-9));
                    // The normalized clean version is always produced (it is
                    // the one that enters the temporal window of every
                    // frame); the natural one is blended after combining.
                    let mut img = to_ref(flatten::clean_frame(model, &frame, Some(&fc), true, opts.scatter_comp));
                    if opts.deflicker {
                        let cur = stats(&img, ow, oh, Some(&pmask));
                        let g = normalize(&mut img, &cur, &reference, opts.cloud_guard);
                        gain_txt += &format!("  gain R/G/B {:.3}/{:.3}/{:.3}  sky {:+.2}%", g[0], g[1], g[2],
                            100.0 * (cur.med[1] - reference.med[1]) / reference.med[1].max(1e-9));
                    }
                    let mut pmask = pmask;
                    let natural = if w_dawn > 0.0 {
                        gain_txt += &format!("  dawn {:.0}%", 100.0 * w_dawn);
                        // A dawn frame is no use as a neighbour in the
                        // temporal window (its normalized version is no
                        // longer the same sky and would drag the frames next
                        // to it): it is marked entirely as non-sky; for
                        // itself the window returns its own value.
                        pmask.iter_mut().for_each(|m| *m = 1);
                        Some((to_ref(flatten::clean_frame(model, &frame, None, true, opts.scatter_comp)), w_dawn))
                    } else {
                        None
                    };
                    Ok((pos, img, pmask, gain_txt, natural))
                })
                .collect();
            for r in loaded {
                buffer.push_back(r?);
            }
            buffer.make_contiguous().sort_by_key(|e| e.0);
            next_load = end;
        }
        // Emit frame `next_emit` with the available window.
        let lo = next_emit.saturating_sub(half);
        let hi = (next_emit + half).min(order.len() - 1);
        let mut win_entries: Vec<&Entry> = buffer
            .iter()
            .filter(|e| e.0 >= lo && e.0 <= hi)
            .collect();
        if win_entries.is_empty() {
            return Err(anyhow!("empty window at frame {}", next_emit));
        }
        // Untracked sequence: each neighbour gets the transform that puts
        // the current frame's sky onto its own, in the cropped reference
        // system the buffered images live in. A neighbour the star chain
        // cannot reach — a frame in between with too few stars, or one that
        // failed to load — drops out of the window instead of being
        // combined out of register.
        let mut ts: Option<Vec<Similarity>> = None;
        if let Some(links) = opts.sky_links {
            let c_idx = order[next_emit];
            let l_c = lref(c_idx);
            let mut keep: Vec<(&Entry, Similarity)> = Vec::with_capacity(win_entries.len());
            for e in win_entries.drain(..) {
                let f_idx = order[e.0];
                if f_idx == c_idx {
                    keep.push((e, Similarity::identity()));
                } else if let Some(s) = fixed::sky_between(links, c_idx, f_idx) {
                    let t = crop_out.compose(
                        &lref(f_idx).inverse().compose(&s.compose(&l_c.compose(&crop_in))),
                    );
                    keep.push((e, t));
                }
            }
            win_entries = keep.iter().map(|(e, _)| *e).collect();
            ts = Some(keep.iter().map(|(_, t)| *t).collect());
        }
        let win: Vec<&[f32]> = win_entries.iter().map(|e| e.1.as_slice()).collect();
        let wmasks: Vec<&[u8]> = win_entries.iter().map(|e| e.2.as_slice()).collect();
        let cur_k = win_entries
            .iter()
            .position(|e| e.0 == next_emit)
            .ok_or_else(|| anyhow!("frame {} is not in the buffer", next_emit))?;
        let cur: &[f32] = win[cur_k];
        let gain_txt = win_entries[cur_k].3.clone();
        let (mut img, used) = match (&ts, win.len()) {
            (_, 1) => (cur.to_vec(), None),
            (Some(ts), _) => {
                let (v, u) = combine_window_warped(&win, &wmasks, ts, cur_k, ow, oh);
                (v, Some(u))
            }
            (None, _) => (combine_window(&win, &wmasks, cur_k), None),
        };
        let mut kept_txt = String::new();
        if win.len() > 1 && opts.cloud_guard {
            // Before the transients: where the cloud covers this frame the
            // combination has already been undone, so there is no excess
            // left there for the transient detection to find.
            let occ = preserve_occlusions(cur, &mut img, wmasks[cur_k], ow, oh);
            if occ > 0.001 {
                kept_txt = format!("  occluded {:.1}%", 100.0 * occ);
            }
        }
        if win.len() > 1 && opts.keep_transients {
            let kept = preserve_transients(cur, &mut img, wmasks[cur_k], used.as_deref(), ow, oh);
            kept_txt += &format!("  transients {} px", kept);
        }
        // Full dawn/twilight: blend with the natural version of the frame.
        if let Some((natural, w)) = &win_entries[cur_k].4 {
            let w = *w;
            img.par_iter_mut().zip(natural.par_iter()).for_each(|(a, b)| *a = *a * (1.0 - w) + *b * w);
        }
        let idx = order[next_emit];
        let stem = paths[idx].file_stem().unwrap().to_string_lossy();
        let out = opts.dir.join(format!("{stem}_clean.dng"));
        output::write_dng_quiet(&out, &img, ow, oh, opts.camera, opts.stretch)?;
        // Each exported frame carries its **own** source frame's metadata,
        // not the reference's: the capture time is what turns the sequence
        // back into a timelapse in the developer.
        match exif::read_source(&paths[idx]) {
            Ok(e) => exif::embed(&out, &e)
                .with_context(|| format!("writing EXIF into {}", out.display()))?,
            Err(_) => exif_failed = true,
        }
        if output::copy_metadata(&paths[idx], &out).is_err() {
            meta_failed = true;
        }
        say!(
            "  [{}/{}] {}  window {} frames{}{}",
            next_emit + 1,
            order.len(),
            out.file_name().unwrap().to_string_lossy(),
            win.len(),
            kept_txt,
            gain_txt
        );
        next_emit += 1;
        crate::ui::task_add("export", 1);
        // Drop frames that no longer fall into any future window.
        let min_keep = next_emit.saturating_sub(half);
        while let Some(e) = buffer.front() {
            if e.0 < min_keep {
                buffer.pop_front();
            } else {
                break;
            }
        }
    }
    if exif_failed {
        say!("  WARNING: could not read the source EXIF on some frames");
    }
    if meta_failed {
        say!("  WARNING: could not copy MakerNotes/XMP with exiftool on some frames (is it installed?) — the DNGs still carry the EXIF written natively");
    }
    say!("exported {} frames in {:.1}s", order.len(), t0.elapsed().as_secs_f32());
    crate::ui::task_end(
        "export",
        format!("{} frames, {:.1}s", order.len(), t0.elapsed().as_secs_f32()),
    );
    Ok(())
}
