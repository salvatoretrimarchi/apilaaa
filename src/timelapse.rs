//! Clean sequence export for timelapse (`--export-clean`).
//!
//! For every frame, in chronological order:
//! 1. Load + clean in sensor coordinates (`flatten::clean_frame` with the
//!    global model plus the frame's own temporal anomaly,
//!    `flatten::fit_frame_corr`: horizon glow/twilight and halo variation
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

use crate::align::Similarity;
use crate::flatten::{self, CellMask};
use crate::output::{self, StretchParams};
use crate::raw;
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

/// Per-channel median and high percentile over ~1M stratified samples.
/// `mask` (optional, one byte per pixel, ≠0 = foreground) excludes those
/// pixels: the tree/horizon must not move the sky median.
pub fn stats(rgb: &[f32], mask: Option<&[u8]>) -> Stats {
    let n_pix = rgb.len() / 3;
    let stride = (n_pix / 1_000_000).max(1);
    let mut s = Stats { med: [0.0; 3], hi: [0.0; 3] };
    for c in 0..3 {
        let mut v: Vec<f32> = (0..n_pix)
            .step_by(stride)
            .filter(|&i| mask.map_or(true, |m| m[i] == 0))
            .map(|i| rgb[i * 3 + c])
            .filter(|x| x.is_finite())
            .collect();
        if v.is_empty() {
            continue;
        }
        let k = v.len() / 2;
        let (_, m, _) = v.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
        s.med[c] = *m;
        let kh = ((v.len() as f32 * HI_PERCENTILE / 100.0) as usize).min(v.len() - 1);
        let (_, h, _) = v.select_nth_unstable_by(kh, |a, b| a.partial_cmp(b).unwrap());
        s.hi[c] = *h;
    }
    s
}

/// Per-channel gain + offset to bring `cur` to `reference`.
/// `v' = (v − med_cur)·g + med_ref`, `g = (hi_ref − med_ref)/(hi_cur − med_cur)`.
pub fn normalize(rgb: &mut [f32], cur: &Stats, reference: &Stats) -> [f32; 3] {
    let mut gain = [1.0f32; 3];
    for c in 0..3 {
        let den = (cur.hi[c] - cur.med[c]).max(1e-6);
        gain[c] = ((reference.hi[c] - reference.med[c]) / den).clamp(0.25, 4.0);
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
/// Radius (px) of the high-pass applied to the difference before the
/// threshold: a trail is narrow; a smooth and extensive change of the sky
/// between frames (twilight, horizon glow, a faint cloud passing) is not a
/// transient and must not cancel the noise reduction over half the image.
const TRANSIENT_HP_RADIUS: usize = 32;
/// Fraction of the sky level below which a pixel is not considered sky for
/// transient purposes (foreground and its penumbra).
const TRANSIENT_SKY_MIN: f32 = 0.6;

/// Separable box moving average (radius `r`, clipped border) — O(n).
fn box_blur(src: &[f32], w: usize, h: usize, r: usize) -> Vec<f32> {
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
pub fn preserve_transients(cur: &[f32], combined: &mut [f32], fg: &[u8], w: usize, h: usize) -> usize {
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
    let lp = box_blur(&ds, w, h, TRANSIENT_HP_RADIUS);
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
    let mut mask: Vec<u8> = ds
        .iter()
        .zip(lum.iter())
        .map(|(&v, &l)| ((v > thr_hi || v < thr_lo) && l >= dark_thr) as u8)
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
            let extent = (maxx - minx + 1).max(maxy - miny + 1);
            if extent < TRANSIENT_MIN_EXTENT {
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
    pub camera_model: &'a str,
    pub stretch: Option<StretchParams>,
    pub n_workers: usize,
    /// Contrast compensation β (see `flatten::clean_frame`).
    pub scatter_comp: f32,
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
    println!(
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

    // Producer: load + clean + warp in parallel batches; consumer: sliding
    // window in order.
    // (position in `order`, image, foreground mask, deflicker text, natural
    // dawn version with its weight)
    type Entry = (usize, Vec<f32>, Vec<u8>, String, Option<(Vec<f32>, f32)>);
    let mut buffer: VecDeque<Entry> = VecDeque::new();
    let mut next_load = 0usize;
    let mut next_emit = 0usize;
    let mut meta_failed = false;
    let batch = opts.n_workers.max(1);

    while next_emit < order.len() {
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
                        (Some(bg), _) => flatten::fit_frame_corr(bg, bg_med, model, Some(&mask)),
                        (None, Some(f)) => flatten::fit_frame_corr(&f.bg, bg_med, model, Some(&mask)),
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
                    let mut gain_txt = format!("  defects ×{:.2}/{:.2}/{:.2}  anomaly {:.1}%", fc.gain[0], fc.gain[1], fc.gain[2], 100.0 * fc.range[1] / reference.med[1].max(1e-9));
                    // The normalized clean version is always produced (it is
                    // the one that enters the temporal window of every
                    // frame); the natural one is blended after combining.
                    let mut img = to_ref(flatten::clean_frame(model, &frame, Some(&fc), true, opts.scatter_comp));
                    if opts.deflicker {
                        let cur = stats(&img, Some(&pmask));
                        let g = normalize(&mut img, &cur, &reference);
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
        let win_entries: Vec<&Entry> = buffer
            .iter()
            .filter(|e| e.0 >= lo && e.0 <= hi)
            .collect();
        if win_entries.is_empty() {
            return Err(anyhow!("empty window at frame {}", next_emit));
        }
        let win: Vec<&[f32]> = win_entries.iter().map(|e| e.1.as_slice()).collect();
        let wmasks: Vec<&[u8]> = win_entries.iter().map(|e| e.2.as_slice()).collect();
        let cur_k = win_entries
            .iter()
            .position(|e| e.0 == next_emit)
            .ok_or_else(|| anyhow!("frame {} is not in the buffer", next_emit))?;
        let cur: &[f32] = win[cur_k];
        let gain_txt = win_entries[cur_k].3.clone();
        let mut img = if win.len() == 1 { cur.to_vec() } else { combine_window(&win, &wmasks, cur_k) };
        let mut kept_txt = String::new();
        if win.len() > 1 && opts.keep_transients {
            let kept = preserve_transients(cur, &mut img, wmasks[cur_k], ow, oh);
            kept_txt = format!("  transients {} px", kept);
        }
        // Full dawn/twilight: blend with the natural version of the frame.
        if let Some((natural, w)) = &win_entries[cur_k].4 {
            let w = *w;
            img.par_iter_mut().zip(natural.par_iter()).for_each(|(a, b)| *a = *a * (1.0 - w) + *b * w);
        }
        let idx = order[next_emit];
        let stem = paths[idx].file_stem().unwrap().to_string_lossy();
        let out = opts.dir.join(format!("{stem}_clean.dng"));
        output::write_dng_quiet(&out, &img, ow, oh, opts.camera_model, opts.stretch)?;
        if output::copy_metadata(&paths[idx], &out).is_err() {
            meta_failed = true;
        }
        println!(
            "  [{}/{}] {}  window {} frames{}{}",
            next_emit + 1,
            order.len(),
            out.file_name().unwrap().to_string_lossy(),
            win.len(),
            kept_txt,
            gain_txt
        );
        next_emit += 1;
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
    if meta_failed {
        println!("  WARNING: could not copy EXIF with exiftool on some frames (is it installed?)");
    }
    println!("exported {} frames in {:.1}s", order.len(), t0.elapsed().as_secs_f32());
    Ok(())
}
