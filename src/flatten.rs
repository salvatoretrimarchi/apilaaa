//! Automatic detection and removal of the lens' central halo/glare.
//!
//! The defect (internal reflection between sensor and rear lens element,
//! "hot spot", out-of-focus pupil ring) is **additive** and **fixed in
//! sensor coordinates**, whereas the sky drifts between frames. That is why
//! it is modelled in the sensor system and subtracted before averaging.
//!
//! Pipeline:
//! 1. Every frame produces a low-resolution background map (`BgMap`):
//!    clipped median over `BLOCK`×`BLOCK` px blocks, in balanced scale
//!    (×wb). Stars take up a minuscule fraction of the block and do not move
//!    the median. The **foreground mask** (`foreground_mask`) is derived
//!    from that same map: trees, mountains, horizon — cells far below the
//!    sky predicted by a smooth fit of the frame itself. Those cells enter
//!    neither any fit nor the stack; in the timelapse sequence they are
//!    shown as is.
//! 2. The maps of the frames selected for the stack are combined by
//!    temporal median ignoring the foreground (`temporal_median_masked`;
//!    rejects planes, satellites, passing clouds).
//! 3. Over that map `poly3(x, y) + Σ f_p(r)·cos(h_p θ)` is fitted: a cubic
//!    polynomial (light-pollution gradient, non-linear) plus free radial
//!    profiles centred on the optical axis (which is searched for, not
//!    assumed at the geometric centre): the pure radial term captures the
//!    halo/ring and the vignetting over the sky; the cos2θ/cos4θ harmonics
//!    capture the non-circular vignetting. All in a single robust least
//!    squares fit (IRLS). Structure that does not fit (Milky Way, nebulae)
//!    is rejected as an outlier and stays out of the model, that is, it is
//!    preserved.
//! 3b. Second stage: over the residual `map − model` a smooth bilinear
//!    surface is fitted (nodes every 64 px, conjugate gradient) of "lower
//!    envelope" type (asymmetric IRLS: it discards whatever stays clearly
//!    above — nebulae, Milky Way — and follows the dark part: bands, flare
//!    wedges, lens hood shadow). Only its **fine** component is subtracted
//!    (surface − Gaussian smoothing σ=`SURF_COARSE_SIGMA_PX`) plus the part
//!    of the coarse component that is **shared with its left/right mirror**
//!    about the vertical axis through the optical centre (system defects are
//!    symmetric; the sky is not). This way symmetric wide wedges and shadows
//!    are removed and extensive sky patches (dark nebulae) are preserved.
//!    Only towards the edges (ramp r/r_full 0.5→0.9).
//!    `--no-residual-surface` disables it.
//! 3c. Third stage: two 1D patterns of fixed geometry (not of the sky) are
//!    extracted from the remaining residual: bands per sensor row (median
//!    per row, high-pass in y) and flare spokes per angle around the optical
//!    axis (median per θ bin, high-pass in θ). They are subtracted over the
//!    whole image.
//! 3d. Fourth stage: horizon glow — robust 1D profile of the residual along
//!    the direction of maximum gradient (searched over 0–360°; the horizon
//!    is not assumed to be at the bottom, the photo may be rotated).
//! 4. Since stacking is linear, the correction is applied at the end:
//!    `out(p) = Σ_i frame_i(m_i p)/N − Σ_i model(m_i p)/N`. The second sum is
//!    evaluated on a coarse grid and interpolated bilinearly.
//!
//! The pedestal (median sky level) is preserved: only the spatial variation
//! is removed, not the absolute level.

use crate::align::Similarity;
use crate::raw::Frame;
use rayon::prelude::*;

/// Block side for the per-frame background map (sensor px).
pub const BLOCK: usize = 32;
/// Step of the coarse grid where the average correction is evaluated (ref px).
pub const CORR_STEP: usize = 16;

/// Low-resolution background map in sensor coordinates, interleaved RGB,
/// balanced scale.
#[derive(Clone)]
pub struct BgMap {
    pub gw: usize,
    pub gh: usize,
    pub block: usize,
    pub data: Vec<f32>,
    /// Cells with real data (false = filled in by extrapolation because it
    /// was always occluded; it enters no fit). Empty = all of them.
    pub valid: Vec<bool>,
}

impl BgMap {
    /// Does cell `i` (row-major order) hold real data?
    #[inline]
    pub fn is_valid(&self, i: usize) -> bool {
        self.valid.is_empty() || self.valid[i]
    }

    /// Exclusion vector (`true` = no data) for the fits.
    fn invalid_vec(&self) -> Vec<bool> {
        (0..self.gw * self.gh).map(|i| !self.is_valid(i)).collect()
    }

    #[inline]
    fn at(&self, gx: usize, gy: usize, c: usize) -> f32 {
        self.data[(gy * self.gw + gx) * 3 + c]
    }

    /// Coordinates (sensor px) of the cell centre.
    #[inline]
    fn cell_center(&self, gx: usize, gy: usize) -> (f32, f32) {
        (
            (gx * self.block) as f32 + self.block as f32 * 0.5,
            (gy * self.block) as f32 + self.block as f32 * 0.5,
        )
    }
}

fn median_inplace(v: &mut [f32]) -> f32 {
    if v.is_empty() {
        return 0.0;
    }
    let k = v.len() / 2;
    let (_, m, _) = v.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
    *m
}

/// Robust median with one clipping pass: discards values above
/// `median + 3·MAD` (stars, hot pixels) and recomputes. It does not clip
/// from below: the background is the foot of the histogram, we do not want
/// to bias it upwards.
fn clipped_median(buf: &mut Vec<f32>, tmp: &mut Vec<f32>) -> f32 {
    let med = median_inplace(buf);
    tmp.clear();
    tmp.extend(buf.iter().map(|v| (v - med).abs()));
    let mad = median_inplace(tmp);
    if mad <= 0.0 {
        return med;
    }
    let hi = med + 3.0 * mad;
    tmp.clear();
    tmp.extend(buf.iter().copied().filter(|&v| v <= hi));
    if tmp.len() < buf.len() / 2 {
        return med;
    }
    median_inplace(tmp)
}

/// Background map of the frame: clipped median per block and channel, ×wb.
pub fn block_background(frame: &Frame) -> BgMap {
    let b = BLOCK;
    let gw = frame.width / b;
    let gh = frame.height / b;
    let w = frame.width;
    let wb = frame.wb;
    let rows: Vec<Vec<f32>> = (0..gh)
        .into_par_iter()
        .map(|gy| {
            let mut out = vec![0.0f32; gw * 3];
            let mut buf: Vec<f32> = Vec::with_capacity(b * b);
            let mut tmp: Vec<f32> = Vec::with_capacity(b * b);
            for gx in 0..gw {
                for c in 0..3 {
                    buf.clear();
                    for y in gy * b..(gy + 1) * b {
                        let row = y * w;
                        for x in gx * b..(gx + 1) * b {
                            buf.push(frame.rgb[(row + x) * 3 + c]);
                        }
                    }
                    out[gx * 3 + c] = clipped_median(&mut buf, &mut tmp) * wb[c];
                }
            }
            out
        })
        .collect();
    let mut data = Vec::with_capacity(gw * gh * 3);
    for r in rows {
        data.extend_from_slice(&r);
    }
    BgMap { gw, gh, block: b, data, valid: Vec::new() }
}

/// Cell-by-cell temporal median of several maps (all the same size).
pub fn temporal_median(maps: &[BgMap]) -> BgMap {
    assert!(!maps.is_empty());
    let gw = maps[0].gw;
    let gh = maps[0].gh;
    let n = gw * gh * 3;
    let data: Vec<f32> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut v: Vec<f32> = maps.iter().map(|m| m.data[i]).collect();
            median_inplace(&mut v)
        })
        .collect();
    BgMap { gw, gh, block: maps[0].block, data, valid: Vec::new() }
}

/// **Foreground** mask per cell of the background map (sensor
/// coordinates): `true` = the cell is taken up by something that is not sky
/// (trees, mountains, horizon, buildings). See `foreground_mask`.
#[derive(Clone)]
pub struct CellMask {
    pub gw: usize,
    pub gh: usize,
    pub block: usize,
    pub data: Vec<bool>,
}

impl CellMask {
    #[inline]
    pub fn at_cell(&self, gx: usize, gy: usize) -> bool {
        self.data[gy * self.gw + gx]
    }

    /// Query by sensor pixel coordinate (outside the map: the nearest
    /// border cell).
    #[inline]
    pub fn at_px(&self, x: f32, y: f32) -> bool {
        let gx = ((x.max(0.0) / self.block as f32) as usize).min(self.gw - 1);
        let gy = ((y.max(0.0) / self.block as f32) as usize).min(self.gh - 1);
        self.data[gy * self.gw + gx]
    }

    /// **Sky** weight at sensor pixel (x, y): 1 on free cells, 0 on
    /// foreground cells, with a bilinear ramp between cell centres (so the
    /// edge of the mask leaves no steps when stacking).
    #[inline]
    pub fn sky_weight(&self, x: f32, y: f32) -> f32 {
        let b = self.block as f32;
        let fx = (x / b - 0.5).clamp(0.0, (self.gw - 1) as f32 - 1e-3);
        let fy = (y / b - 0.5).clamp(0.0, (self.gh - 1) as f32 - 1e-3);
        let i0 = fx.floor() as usize;
        let j0 = fy.floor() as usize;
        let tx = fx - i0 as f32;
        let ty = fy - j0 as f32;
        let i1 = (i0 + 1).min(self.gw - 1);
        let j1 = (j0 + 1).min(self.gh - 1);
        let m = |i: usize, j: usize| if self.data[j * self.gw + i] { 0.0f32 } else { 1.0f32 };
        m(i0, j0) * (1.0 - tx) * (1.0 - ty) + m(i1, j0) * tx * (1.0 - ty) + m(i0, j1) * (1.0 - tx) * ty + m(i1, j1) * tx * ty
    }

    /// Fraction of cells taken up by foreground.
    pub fn fraction(&self) -> f32 {
        self.data.iter().filter(|&&m| m).count() as f32 / self.data.len().max(1) as f32
    }

    pub fn any(&self) -> bool {
        self.data.iter().any(|&m| m)
    }
}

/// Foreground threshold: a cell is foreground if its background (G) falls
/// below `FG_RATIO` × the sky predicted by the smooth fit of the frame
/// itself. An opaque object against the sky is nearly black (≪ 0.5×); a
/// dark nebula or the real vignetting never drop that far relative to a fit
/// that already includes the radial profile.
const FG_RATIO: f32 = 0.55;
/// Minimum size (cells, 4-connected) of a foreground component.
const FG_MIN_CELLS: usize = 2;
/// Dilation (cells) of the mask: out-of-focus edges of leaves/branches and
/// the half-cells the threshold leaves out.
const FG_DILATE: usize = 1;

/// Detects the foreground of a frame from its background map, without
/// needing the global model: the same smooth `poly3 + radial profiles`
/// model (geometric centre, IRLS) is fitted to the map itself in G, which
/// follows halo, vignetting and gradient but cannot follow a silhouette
/// cutout; cells far below that prediction (`< FG_RATIO`×) are foreground.
/// Two rounds (the second one already without those cells in the fit), a
/// small-component filter and a dilation of `FG_DILATE` cells.
pub fn foreground_mask(map: &BgMap) -> CellMask {
    let n = map.gw * map.gh;
    let w = (map.gw * map.block) as f32;
    let h = (map.gh * map.block) as f32;
    let r_step = (map.block * 2) as f32;
    let mut mask = vec![false; n];
    for _round in 0..2 {
        let fit = fit_channel_ex(map, 1, w * 0.5, h * 0.5, r_step, Some(&mask));
        for i in 0..n {
            let pred = fit.eval_cell(i, r_step);
            let v = map.data[i * 3 + 1];
            mask[i] = pred > 0.0 && v < FG_RATIO * pred;
        }
    }
    // Drop small connected components.
    {
        let (gw, gh) = (map.gw, map.gh);
        let mut visited = vec![false; n];
        let mut stack: Vec<usize> = Vec::new();
        let mut comp: Vec<usize> = Vec::new();
        for start in 0..n {
            if !mask[start] || visited[start] { continue; }
            comp.clear();
            stack.push(start);
            visited[start] = true;
            while let Some(i) = stack.pop() {
                comp.push(i);
                let x = i % gw;
                let y = i / gw;
                let nb = [
                    if x > 0 { Some(i - 1) } else { None },
                    if x + 1 < gw { Some(i + 1) } else { None },
                    if y > 0 { Some(i - gw) } else { None },
                    if y + 1 < gh { Some(i + gw) } else { None },
                ];
                for j in nb.into_iter().flatten() {
                    if mask[j] && !visited[j] { visited[j] = true; stack.push(j); }
                }
            }
            if comp.len() < FG_MIN_CELLS {
                for &i in &comp { mask[i] = false; }
            }
        }
    }
    // Dilation.
    if FG_DILATE > 0 && mask.iter().any(|&m| m) {
        let (gw, gh) = (map.gw, map.gh);
        let r = FG_DILATE;
        let mut out = vec![false; n];
        for gy in 0..gh {
            for gx in 0..gw {
                if !mask[gy * gw + gx] { continue; }
                for yy in gy.saturating_sub(r)..=(gy + r).min(gh - 1) {
                    for xx in gx.saturating_sub(r)..=(gx + r).min(gw - 1) {
                        out[yy * gw + xx] = true;
                    }
                }
            }
        }
        mask = out;
    }
    CellMask { gw: map.gw, gh: map.gh, block: map.block, data: mask }
}

/// Sky level of the frame: median of G over the map cells that are not
/// foreground (balanced scale).
pub fn sky_level(map: &BgMap, mask: &CellMask) -> f32 {
    let mut v: Vec<f32> = (0..map.gw * map.gh)
        .filter(|&i| !mask.data[i])
        .map(|i| map.data[i * 3 + 1])
        .collect();
    if v.is_empty() {
        return 0.0;
    }
    median_inplace(&mut v)
}

/// Minimum number of valid samples per cell in the masked temporal median;
/// below it the cell is filled in from its neighbours.
const TM_MIN_SAMPLES: usize = 5;

/// Cell-by-cell temporal median ignoring, in every frame, the foreground
/// cells. Cells with fewer than `TM_MIN_SAMPLES` samples (always occluded)
/// are filled in iteratively with the mean of their valid neighbours. Also
/// returns the number of filled cells.
pub fn temporal_median_masked(maps: &[BgMap], masks: &[CellMask]) -> (BgMap, usize) {
    assert!(!maps.is_empty() && maps.len() == masks.len());
    let gw = maps[0].gw;
    let gh = maps[0].gh;
    let ncell = gw * gh;
    let cells: Vec<([f32; 3], bool)> = (0..ncell)
        .into_par_iter()
        .map(|i| {
            let mut out = [0.0f32; 3];
            let mut cnt = 0usize;
            for c in 0..3 {
                let mut v: Vec<f32> = maps
                    .iter()
                    .zip(masks)
                    .filter(|(_, mk)| !mk.data[i])
                    .map(|(m, _)| m.data[i * 3 + c])
                    .collect();
                cnt = v.len();
                out[c] = if v.is_empty() { f32::NAN } else { median_inplace(&mut v) };
            }
            (out, cnt >= TM_MIN_SAMPLES)
        })
        .collect();
    let mut data = vec![0.0f32; ncell * 3];
    let mut valid = vec![false; ncell];
    for (i, (v, ok)) in cells.into_iter().enumerate() {
        data[i * 3..i * 3 + 3].copy_from_slice(&v);
        valid[i] = ok && v.iter().all(|x| x.is_finite());
    }
    let filled = ncell - valid.iter().filter(|&&v| v).count();
    let real_valid = valid.clone();
    // Iterative fill-in from valid neighbours (8-connected).
    let mut pending = filled;
    while pending > 0 {
        let mut new_valid = valid.clone();
        let mut new_data = data.clone();
        let mut progress = false;
        for gy in 0..gh {
            for gx in 0..gw {
                let i = gy * gw + gx;
                if valid[i] { continue; }
                let mut acc = [0.0f32; 3];
                let mut k = 0usize;
                for yy in gy.saturating_sub(1)..=(gy + 1).min(gh - 1) {
                    for xx in gx.saturating_sub(1)..=(gx + 1).min(gw - 1) {
                        let j = yy * gw + xx;
                        if !valid[j] { continue; }
                        for c in 0..3 { acc[c] += data[j * 3 + c]; }
                        k += 1;
                    }
                }
                if k > 0 {
                    for c in 0..3 { new_data[i * 3 + c] = acc[c] / k as f32; }
                    new_valid[i] = true;
                    progress = true;
                    pending -= 1;
                }
            }
        }
        valid = new_valid;
        data = new_data;
        if !progress {
            // No valid cell in the whole map: leave zeros.
            for (i, v) in valid.iter().enumerate() {
                if !v { for c in 0..3 { data[i * 3 + c] = 0.0; } }
            }
            break;
        }
    }
    let valid = if filled > 0 { real_valid } else { Vec::new() };
    (BgMap { gw, gh, block: maps[0].block, data, valid }, filled)
}

/// Number of terms of the background polynomial (full degree 3:
/// `1, u, v, u², uv, v², u³, u²v, uv², v³`) with `u, v ∈ [−1, 1]`
/// normalized sensor coordinates.
const N_POLY: usize = 10;

/// Per-channel background model in sensor coordinates:
/// `poly3(u, v) + f(r)`, with `f` a free radial profile (piecewise linear)
/// around the optical centre out to the farthest corner. The polynomial
/// absorbs the light-pollution gradient (non-linear); `f(r)` the lens'
/// halo/ring and the bowl of the vignetting over the sky.
pub struct GlareModel {
    pub cx: f32,
    pub cy: f32,
    pub poly: [[f32; N_POLY]; 3],
    /// Per-channel profiles (index 0 pure radial, then `HARMONICS`),
    /// sampled every `r_step` px starting at r=0.
    pub profiles: [Vec<Vec<f32>>; 3],
    pub r_step: f32,
    pub r_full: f32,
    /// Level that is preserved (median of the model over the frame) per channel.
    pub pedestal: [f32; 3],
    /// Robust residual (MAD of the fit) per channel — diagnostic.
    pub residual_mad: [f32; 3],
    /// "Lower envelope" residual surface (see `fit_residual_surface`);
    /// None if disabled.
    pub surface: Option<Surface>,
    /// Coarse surface (same lower envelope at a ~600 px scale). It is
    /// subtracted from the fine one so that only **fine** structure (lines,
    /// wedges, bands) is removed and not extensive sky patches (dark
    /// nebulae).
    pub surface_coarse: Option<Surface>,
    /// Fixed 1D patterns estimated from the final residual (see
    /// `fit_fixed_lines`): horizontal sensor bands (one value per cell row,
    /// step `block` px, high-passed) and flare spokes (one value per
    /// angular bin of `SPOKE_BIN_DEG`, high-passed in θ). None if disabled.
    pub lines: Option<FixedLines>,
    /// Light-pollution glow from the horizon: smooth 1D profile along the
    /// direction of maximum residual gradient (see `fit_horizon`). None if
    /// disabled.
    pub horizon: Option<Horizon>,
    /// Sensor width/height.
    pub width: usize,
    pub height: usize,
}

/// Smooth bilinear surface over nodes every `step` px (sensor
/// coordinates), RGB per node. Models the structured residual the
/// parametric model does not capture (bands, flare wedges, mechanical
/// vignetting).
pub struct Surface {
    pub nx: usize,
    pub ny: usize,
    pub step: f32,
    pub data: Vec<f32>, // nx*ny*3
    /// Peak−trough amplitude per channel (diagnostic).
    pub range: [f32; 3],
}

impl Surface {
    #[inline]
    pub fn eval(&self, x: f32, y: f32) -> [f32; 3] {
        let fx = (x / self.step).clamp(0.0, (self.nx - 1) as f32 - 1e-3);
        let fy = (y / self.step).clamp(0.0, (self.ny - 1) as f32 - 1e-3);
        let i0 = fx.floor() as usize;
        let j0 = fy.floor() as usize;
        let tx = fx - i0 as f32;
        let ty = fy - j0 as f32;
        let idx = |i: usize, j: usize| (j * self.nx + i) * 3;
        let (a, b, c, d) = (idx(i0, j0), idx(i0 + 1, j0), idx(i0, j0 + 1), idx(i0 + 1, j0 + 1));
        let mut out = [0.0f32; 3];
        for k in 0..3 {
            out[k] = self.data[a + k] * (1.0 - tx) * (1.0 - ty)
                + self.data[b + k] * tx * (1.0 - ty)
                + self.data[c + k] * (1.0 - tx) * ty
                + self.data[d + k] * tx * ty;
        }
        out
    }
}

/// Fills NaN by linear interpolation between valid neighbours and constant
/// extension at the ends.
fn fill_nan(p: &mut [f32]) {
    let n = p.len();
    let Some(first) = p.iter().position(|v| !v.is_nan()) else {
        for v in p.iter_mut() { *v = 0.0; }
        return;
    };
    let last = p.iter().rposition(|v| !v.is_nan()).unwrap();
    for i in 0..first { p[i] = p[first]; }
    for i in last + 1..n { p[i] = p[last]; }
    let mut i = first;
    while i < last {
        if p[i + 1].is_nan() {
            let j = (i + 1..=last).find(|&j| !p[j].is_nan()).unwrap();
            for k in i + 1..j {
                let t = (k - i) as f32 / (j - i) as f32;
                p[k] = p[i] * (1.0 - t) + p[j] * t;
            }
            i = j;
        } else {
            i += 1;
        }
    }
}

/// Step (px) of the horizon profile nodes and smoothing window.
const HORIZON_STEP_PX: f32 = 96.0;
const HORIZON_SMOOTH: usize = 7;

/// Horizon glow: 1D profile `p(d)` with `d` the projection onto the `phi`
/// direction (the way the glow grows), piecewise linear.
pub struct Horizon {
    pub cos_phi: f32,
    pub sin_phi: f32,
    pub d_min: f32,
    pub step: f32,
    pub profile: Vec<[f32; 3]>,
    pub phi_deg: f32,
    pub range: [f32; 3],
}

impl Horizon {
    #[inline]
    pub fn eval(&self, dx: f32, dy: f32) -> [f32; 3] {
        let d = dx * self.cos_phi + dy * self.sin_phi;
        let t = ((d - self.d_min) / self.step).clamp(0.0, (self.profile.len() - 1) as f32 - 1e-3);
        let i = t.floor() as usize;
        let f = t - i as f32;
        let mut out = [0.0f32; 3];
        for c in 0..3 {
            out[c] = self.profile[i][c] * (1.0 - f) + self.profile[i + 1][c] * f;
        }
        out
    }
}

/// Estimates the horizon glow from the residual `map − model`: for a
/// candidate direction, a robust median of the residual per projection bin
/// `d` (non-outlier cells: |res| ≤ 2·MAD, so the Milky Way and bright
/// nebulae do not lift it; a localized dark patch does not lower it because
/// the median is taken along the whole perpendicular line), smoothed by a
/// moving median. The `phi` (step 2°) minimizing the residual MAD over the
/// G channel is chosen, and then each channel's profile is estimated. The
/// profile minus its median is subtracted (the pedestal is preserved
/// separately).
fn fit_horizon(map: &BgMap, model: &GlareModel) -> Horizon {
    let n = map.gw * map.gh;
    let mut resid = vec![[0.0f32; 3]; n];
    for gy in 0..map.gh {
        for gx in 0..map.gw {
            let (x, y) = map.cell_center(gx, gy);
            let m = model.eval_raw_rgb(x, y);
            let i = gy * map.gw + gx;
            for c in 0..3 { resid[i][c] = map.at(gx, gy, c) - m[c]; }
        }
    }
    let ok = |i: usize| map.is_valid(i);
    fit_horizon_resid(map, model, &resid, &ok, None)
}

/// Core of `fit_horizon` over an already computed per-cell residual. `ok`
/// excludes cells a priori (no data, foreground). `phi_hint` = (central φ,
/// half-width) in degrees restricts the direction search.
fn fit_horizon_resid(
    map: &BgMap,
    model: &GlareModel,
    resid: &[[f32; 3]],
    ok_pre: &dyn Fn(usize) -> bool,
    phi_hint: Option<(f32, f32)>,
) -> Horizon {
    let n = map.gw * map.gh;
    let mut xs = vec![0.0f32; n];
    let mut ys = vec![0.0f32; n];
    for gy in 0..map.gh {
        for gx in 0..map.gw {
            let (x, y) = map.cell_center(gx, gy);
            let i = gy * map.gw + gx;
            xs[i] = x - model.cx;
            ys[i] = y - model.cy;
        }
    }
    let step = HORIZON_STEP_PX;
    let r_corner = corner_radius(model.cx, model.cy, model.width, model.height);
    let d_min = -r_corner;
    let n_bins = (2.0 * r_corner / step).ceil() as usize + 2;

    // Valid cells per channel (bright and very dark ones are rejected).
    let mut ok = vec![[true; 3]; n];
    for c in 0..3 {
        let mut abs: Vec<f32> = (0..n).filter(|&i| ok_pre(i)).map(|i| resid[i][c].abs()).collect();
        let mad = median_inplace(&mut abs).max(1e-9);
        for i in 0..n { ok[i][c] = ok_pre(i) && resid[i][c] <= 2.0 * mad && resid[i][c] >= -5.0 * mad; }
    }

    let profile_for = |cs: f32, sn: f32, c: usize| -> Vec<f32> {
        let mut bins: Vec<Vec<f32>> = vec![Vec::new(); n_bins];
        for i in 0..n {
            if !ok[i][c] { continue; }
            let d = xs[i] * cs + ys[i] * sn;
            let b = (((d - d_min) / step) as usize).min(n_bins - 1);
            bins[b].push(resid[i][c]);
        }
        let mut raw = vec![f32::NAN; n_bins];
        for (b, v) in bins.iter_mut().enumerate() {
            if v.len() >= 6 { raw[b] = median_inplace(v); }
        }
        fill_nan(&mut raw);
        // Smoothing by moving median.
        let half = HORIZON_SMOOTH / 2;
        let mut out = vec![0.0f32; n_bins];
        let mut buf = Vec::with_capacity(HORIZON_SMOOTH);
        for b in 0..n_bins {
            buf.clear();
            for k in b.saturating_sub(half)..=(b + half).min(n_bins - 1) { buf.push(raw[k]); }
            out[b] = median_inplace(&mut buf);
        }
        // Centre on its median: the level is carried by the pedestal.
        let mut tmp = out.clone();
        let med = median_inplace(&mut tmp);
        for v in out.iter_mut() { *v -= med; }
        out
    };
    let score = |cs: f32, sn: f32| -> f32 {
        let prof = profile_for(cs, sn, 1);
        let mut abs: Vec<f32> = (0..n)
            .filter(|&i| ok[i][1])
            .map(|i| {
                let d = xs[i] * cs + ys[i] * sn;
                let t = ((d - d_min) / step).clamp(0.0, (n_bins - 1) as f32 - 1e-3);
                let k = t.floor() as usize;
                let f = t - k as f32;
                (resid[i][1] - (prof[k] * (1.0 - f) + prof[k + 1] * f)).abs()
            })
            .collect();
        median_inplace(&mut abs)
    };
    let cands: Vec<f32> = match phi_hint {
        None => (0..180).map(|k| k as f32 * 2.0).collect(),
        Some((c0, half)) => {
            let n_c = (half / 2.0).ceil() as i32;
            (-n_c..=n_c).map(|k| (c0 + k as f32 * 2.0).rem_euclid(360.0)).collect()
        }
    };
    let scored: Vec<(f32, f32)> = cands
        .par_iter()
        .map(|&deg| {
            let (sn, cs) = deg.to_radians().sin_cos();
            (deg, score(cs, sn))
        })
        .collect();
    let (phi_deg, _) = scored.into_iter().fold((0.0f32, f32::INFINITY), |b, s| if s.1 < b.1 { s } else { b });
    let (sin_phi, cos_phi) = phi_deg.to_radians().sin_cos();
    let mut profile = vec![[0.0f32; 3]; n_bins];
    let mut range = [0.0f32; 3];
    for c in 0..3 {
        let p = profile_for(cos_phi, sin_phi, c);
        for b in 0..n_bins { profile[b][c] = p[b]; }
        range[c] = p.iter().cloned().fold(f32::MIN, f32::max) - p.iter().cloned().fold(f32::MAX, f32::min);
    }
    Horizon { cos_phi, sin_phi, d_min, step, profile, phi_deg, range }
}

/// Angular bin of the flare spokes (degrees) and minimum radius (in
/// r_full) from which they are applied (there is no statistics near the
/// centre).
const SPOKE_BIN_DEG: f32 = 0.5;
const SPOKE_R_MIN: f32 = 0.15;
/// High-pass window (in samples) for rows and for θ.
const LINES_HP_ROWS: usize = 9;
const LINES_HP_THETA: usize = 41;

/// Fixed 1D patterns of the residual: bands per sensor row and flare
/// spokes per angle around the optical centre. Both are fixed geometry of
/// the system (sensor rows, optical axis), not of the sky, so they are
/// subtracted over the whole image without a ramp.
pub struct FixedLines {
    pub row_step: f32,
    pub rows: Vec<[f32; 3]>,   // per cell row, high-passed
    pub spokes: Vec<[f32; 3]>, // per angular bin, high-passed
    /// Peak−trough amplitude (diagnostic) per channel: rows, spokes.
    pub rows_range: [f32; 3],
    pub spokes_range: [f32; 3],
}

impl FixedLines {
    #[inline]
    pub fn eval(&self, y: f32, r: f32, theta: f32, r_full: f32) -> [f32; 3] {
        // Rows: linear interpolation between row centres.
        let fy = (y / self.row_step - 0.5).clamp(0.0, (self.rows.len() - 1) as f32 - 1e-3);
        let j = fy.floor() as usize;
        let t = fy - j as f32;
        let mut out = [0.0f32; 3];
        for c in 0..3 {
            out[c] = self.rows[j][c] * (1.0 - t) + self.rows[j + 1][c] * t;
        }
        // Spokes: only beyond r > SPOKE_R_MIN·r_full with a short ramp.
        let w = ((r / r_full - SPOKE_R_MIN) / SPOKE_R_MIN).clamp(0.0, 1.0);
        if w > 0.0 {
            let nb = self.spokes.len();
            let mut deg = theta.to_degrees();
            if deg < 0.0 { deg += 360.0; }
            let fb = deg / SPOKE_BIN_DEG;
            let b0 = (fb.floor() as usize) % nb;
            let b1 = (b0 + 1) % nb;
            let tb = fb - fb.floor();
            for c in 0..3 {
                out[c] += (self.spokes[b0][c] * (1.0 - tb) + self.spokes[b1][c] * tb) * w;
            }
        }
        out
    }
}

/// 1D high-pass median: `v − moving median(v, window)`. `wrap` = circular.
fn highpass_median(v: &[f32], window: usize, wrap: bool) -> Vec<f32> {
    let n = v.len();
    let half = window / 2;
    let mut out = vec![0.0f32; n];
    let mut buf = Vec::with_capacity(window);
    for i in 0..n {
        buf.clear();
        for d in 0..window {
            let k = i as isize + d as isize - half as isize;
            let idx = if wrap {
                Some(k.rem_euclid(n as isize) as usize)
            } else if k >= 0 && (k as usize) < n {
                Some(k as usize)
            } else {
                None
            };
            if let Some(idx) = idx {
                if v[idx].is_finite() { buf.push(v[idx]); }
            }
        }
        out[i] = if v[i].is_finite() && !buf.is_empty() { v[i] - median_inplace(&mut buf) } else { 0.0 };
    }
    out
}

/// Estimates `FixedLines` from the residual `map − model(so far)`. Per
/// row: median over x of the non-outlier cells (|res| ≤ 3·global MAD); per
/// angular bin: median over r ≥ SPOKE_R_MIN·r_full. Both high-passed so as
/// not to touch smooth structure (already modelled, or real).
fn fit_fixed_lines(map: &BgMap, model: &GlareModel) -> FixedLines {
    let n = map.gw * map.gh;
    let mut resid = vec![[0.0f32; 3]; n];
    let mut rs = vec![0.0f32; n];
    let mut ths = vec![0.0f32; n];
    for gy in 0..map.gh {
        for gx in 0..map.gw {
            let (x, y) = map.cell_center(gx, gy);
            let m = model.eval_raw_rgb(x, y);
            let i = gy * map.gw + gx;
            for c in 0..3 { resid[i][c] = map.at(gx, gy, c) - m[c]; }
            let dx = x - model.cx;
            let dy = y - model.cy;
            rs[i] = (dx * dx + dy * dy).sqrt();
            ths[i] = dy.atan2(dx);
        }
    }
    let nb = (360.0 / SPOKE_BIN_DEG) as usize;
    let mut rows = vec![[0.0f32; 3]; map.gh];
    let mut spokes = vec![[0.0f32; 3]; nb];
    let mut rows_range = [0.0f32; 3];
    let mut spokes_range = [0.0f32; 3];
    for c in 0..3 {
        let mut abs: Vec<f32> = (0..n).filter(|&i| map.is_valid(i)).map(|i| resid[i][c].abs()).collect();
        let mad = median_inplace(&mut abs).max(1e-9);
        let ok = |i: usize| map.is_valid(i) && resid[i][c].abs() <= 3.0 * mad;
        // Rows
        let mut raw = vec![f32::NAN; map.gh];
        let mut buf = Vec::with_capacity(map.gw);
        for gy in 0..map.gh {
            buf.clear();
            for gx in 0..map.gw {
                let i = gy * map.gw + gx;
                if ok(i) { buf.push(resid[i][c]); }
            }
            if buf.len() >= map.gw / 4 { raw[gy] = median_inplace(&mut buf); }
        }
        let hp = highpass_median(&raw, LINES_HP_ROWS, false);
        for gy in 0..map.gh { rows[gy][c] = hp[gy]; }
        rows_range[c] = hp.iter().cloned().fold(f32::MIN, f32::max) - hp.iter().cloned().fold(f32::MAX, f32::min);
        // Spokes (over the residual with the rows already removed)
        let mut bins: Vec<Vec<f32>> = vec![Vec::new(); nb];
        for i in 0..n {
            if !ok(i) || rs[i] < SPOKE_R_MIN * model.r_full { continue; }
            let gy = i / map.gw;
            let mut deg = ths[i].to_degrees();
            if deg < 0.0 { deg += 360.0; }
            let b = ((deg / SPOKE_BIN_DEG) as usize).min(nb - 1);
            bins[b].push(resid[i][c] - rows[gy][c]);
        }
        let mut raw = vec![f32::NAN; nb];
        for (b, v) in bins.iter_mut().enumerate() {
            if v.len() >= 6 { raw[b] = median_inplace(v); }
        }
        let hp = highpass_median(&raw, LINES_HP_THETA, true);
        for b in 0..nb { spokes[b][c] = hp[b]; }
        spokes_range[c] = hp.iter().cloned().fold(f32::MIN, f32::max) - hp.iter().cloned().fold(f32::MAX, f32::min);
    }
    FixedLines { row_step: map.block as f32, rows, spokes, rows_range, spokes_range }
}

/// Node spacing of the residual surface, in map cells (2 cells = 64 px:
/// resolves thin flare lines and wedges).
const SURF_STEP_CELLS: usize = 2;
/// Smoothness of the surface (second difference between nodes, in
/// equivalent cells) and maximum conjugate gradient iterations.
const SURF_LAMBDA: f64 = 2.0;
/// Sigma (px) of the Gaussian smoothing that defines the "coarse"
/// component of the surface. Only `fine − coarse` is subtracted: structure
/// wider than ~2·sigma (extensive patches, dark nebulae) is preserved.
const SURF_COARSE_SIGMA_PX: f32 = 150.0;

/// Smoothed copy (separable Gaussian over the node grid, replicated
/// border) of a surface.
fn smooth_surface(sf: &Surface, sigma_px: f32) -> Surface {
    let sigma = sigma_px / sf.step;
    let rad = (3.0 * sigma).ceil() as isize;
    let kern: Vec<f32> = (-rad..=rad).map(|k| (-(k * k) as f32 / (2.0 * sigma * sigma)).exp()).collect();
    let ksum: f32 = kern.iter().sum();
    let (nx, ny) = (sf.nx, sf.ny);
    let mut tmp = vec![0.0f32; nx * ny * 3];
    for j in 0..ny {
        for i in 0..nx {
            for c in 0..3 {
                let mut acc = 0.0;
                for (kk, &kw) in kern.iter().enumerate() {
                    let ii = (i as isize + kk as isize - rad).clamp(0, nx as isize - 1) as usize;
                    acc += kw * sf.data[(j * nx + ii) * 3 + c];
                }
                tmp[(j * nx + i) * 3 + c] = acc / ksum;
            }
        }
    }
    let mut data = vec![0.0f32; nx * ny * 3];
    for j in 0..ny {
        for i in 0..nx {
            for c in 0..3 {
                let mut acc = 0.0;
                for (kk, &kw) in kern.iter().enumerate() {
                    let jj = (j as isize + kk as isize - rad).clamp(0, ny as isize - 1) as usize;
                    acc += kw * tmp[(jj * nx + i) * 3 + c];
                }
                data[(j * nx + i) * 3 + c] = acc / ksum;
            }
        }
    }
    Surface { nx, ny, step: sf.step, data, range: sf.range }
}
const SURF_CG_ITERS: usize = 400;
/// Radial ramp over which the surface is applied, in units of `r_full`.
const SURF_R_START: f32 = 0.5;
const SURF_R_END: f32 = 0.9;

/// Fits to the residual `map − parametric model` a smooth bilinear surface
/// that follows the **lower envelope** of the sky: IRLS excludes the cells
/// clearly above (nebulae, Milky Way, fat stars: residual > +2·MAD) and the
/// ones far below (foreground objects: < −5·MAD). What is left — dark
/// bands, flare wedges, lens hood shadow — is either a defect fixed to the
/// sensor or featureless sky, and is subtracted. Bright structure stays
/// above the surface and is preserved. Second-difference regularization
/// between nodes.
fn fit_residual_surface(map: &BgMap, model: &GlareModel, step_cells: usize, lambda: f64) -> Surface {
    let n = map.gw * map.gh;
    let mut resid = vec![[0.0f32; 3]; n];
    for gy in 0..map.gh {
        for gx in 0..map.gw {
            let (x, y) = map.cell_center(gx, gy);
            let m = model.eval_raw_rgb(x, y);
            let i = gy * map.gw + gx;
            for c in 0..3 {
                resid[i][c] = map.at(gx, gy, c) - m[c];
            }
        }
    }
    let ok = |i: usize| map.is_valid(i);
    fit_surface_core(map.gw, map.gh, map.block, &resid, &ok, step_cells, lambda, false, None, true).0
}

/// Core of the smooth bilinear surface fit to a per-cell residual. `ok`
/// excludes cells a priori. `symmetric`: rejection |res| > 2.5·MAD (to
/// follow extensive structure above and below alike); otherwise, lower
/// envelope (rejects > +2·MAD and < −5·MAD). `extra`: an additional
/// per-cell basis (RGB) fitted jointly with one coefficient per channel
/// (returned; 0 if there is no `extra`).
fn fit_surface_core(
    gw: usize,
    gh: usize,
    block: usize,
    resid: &[[f32; 3]],
    ok: &(dyn Fn(usize) -> bool + Sync),
    step_cells: usize,
    lambda: f64,
    symmetric: bool,
    extra: Option<&[[f32; 3]]>,
    robust: bool,
) -> (Surface, [f32; 3]) {
    let step = (step_cells * block) as f32;
    let nx = gw / step_cells + 2;
    let ny = gh / step_cells + 2;
    let nn = nx * ny;
    let n = gw * gh;

    // Bilinear basis (4 nonzeros per cell).
    let mut basis: Vec<([usize; 4], [f64; 4])> = Vec::with_capacity(n);
    for gy in 0..gh {
        for gx in 0..gw {
            let x = (gx * block) as f32 + block as f32 * 0.5;
            let y = (gy * block) as f32 + block as f32 * 0.5;
            let fx = x / step;
            let fy = y / step;
            let i0 = (fx.floor() as usize).min(nx - 2);
            let j0 = (fy.floor() as usize).min(ny - 2);
            let tx = (fx - i0 as f32) as f64;
            let ty = (fy - j0 as f32) as f64;
            basis.push((
                [j0 * nx + i0, j0 * nx + i0 + 1, (j0 + 1) * nx + i0, (j0 + 1) * nx + i0 + 1],
                [(1.0 - tx) * (1.0 - ty), tx * (1.0 - ty), (1.0 - tx) * ty, tx * ty],
            ));
        }
    }

    // Sparse normal system solved by matrix-free conjugate gradient:
    // (Aᵀ W A + λ Lᵀ L + ε I) s = Aᵀ W r. Every cell touches 4 nodes; L is
    // the second difference in x and y between nodes.
    let lambda_smooth = lambda;
    // Unknowns: nn nodes (+1 coefficient for `extra`, index nn).
    let nu = nn + 1;
    let matvec = |v: &[f64], out: &mut [f64], weight: &[bool], ex: &[f64]| {
        for o in out.iter_mut() { *o = 0.0; }
        for i in 0..n {
            if !weight[i] { continue; }
            let (ids, vs) = &basis[i];
            let mut sv = v[nn] * ex[i];
            for a in 0..4 { sv += v[ids[a]] * vs[a]; }
            for a in 0..4 { out[ids[a]] += vs[a] * sv; }
            out[nn] += ex[i] * sv;
        }
        out[nn] += 1e-9 * v[nn];
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                if i >= 1 && i + 1 < nx {
                    let d = v[k - 1] - 2.0 * v[k] + v[k + 1];
                    out[k - 1] += lambda_smooth * d;
                    out[k] -= 2.0 * lambda_smooth * d;
                    out[k + 1] += lambda_smooth * d;
                }
                if j >= 1 && j + 1 < ny {
                    let d = v[k - nx] - 2.0 * v[k] + v[k + nx];
                    out[k - nx] += lambda_smooth * d;
                    out[k] -= 2.0 * lambda_smooth * d;
                    out[k + nx] += lambda_smooth * d;
                }
                out[k] += 1e-6 * v[k];
            }
        }
    };
    let mut data = vec![0.0f32; nn * 3];
    let mut range = [0.0f32; 3];
    let mut coef = [0.0f32; 3];
    let per_channel: Vec<(Vec<f32>, f32, f32)> = (0..3)
        .into_par_iter()
        .map(|c| {
            let ex: Vec<f64> = (0..n).map(|i| extra.map_or(0.0, |e| e[i][c] as f64)).collect();
            let mut weight: Vec<bool> = (0..n).map(|i| ok(i)).collect();
            let mut sol = vec![0.0f64; nu];
            let mut rhs = vec![0.0f64; nu];
            let (mut rv, mut pv, mut ap) = (vec![0.0f64; nu], vec![0.0f64; nu], vec![0.0f64; nu]);
            let n_iter = if robust { 6 } else { 1 };
            for iter in 0..n_iter {
                for v in rhs.iter_mut() { *v = 0.0; }
                for i in 0..n {
                    if !weight[i] { continue; }
                    let (ids, vs) = &basis[i];
                    let y = resid[i][c] as f64;
                    for a in 0..4 { rhs[ids[a]] += vs[a] * y; }
                    rhs[nn] += ex[i] * y;
                }
                // CG starting from the previous solution.
                matvec(&sol, &mut ap, &weight, &ex);
                for k in 0..nu { rv[k] = rhs[k] - ap[k]; pv[k] = rv[k]; }
                let mut rr: f64 = rv.iter().map(|x| x * x).sum();
                let rr0 = rr.max(1e-30);
                for _ in 0..SURF_CG_ITERS {
                    if rr / rr0 < 1e-12 { break; }
                    matvec(&pv, &mut ap, &weight, &ex);
                    let pap: f64 = pv.iter().zip(&ap).map(|(a, b)| a * b).sum();
                    if pap <= 0.0 { break; }
                    let alpha = rr / pap;
                    for k in 0..nu { sol[k] += alpha * pv[k]; rv[k] -= alpha * ap[k]; }
                    let rr_new: f64 = rv.iter().map(|x| x * x).sum();
                    let beta = rr_new / rr;
                    rr = rr_new;
                    for k in 0..nu { pv[k] = rv[k] + beta * pv[k]; }
                }
                // Residual with respect to the surface; clipping.
                let mut res = vec![0.0f32; n];
                for i in 0..n {
                    let (ids, vs) = &basis[i];
                    let mut sv = sol[nn] * ex[i];
                    for a in 0..4 { sv += sol[ids[a]] * vs[a]; }
                    res[i] = resid[i][c] - sv as f32;
                }
                let mut abs: Vec<f32> = (0..n).filter(|&i| ok(i)).map(|i| res[i].abs()).collect();
                let mad = median_inplace(&mut abs).max(1e-9);
                if iter < 5 {
                    for i in 0..n {
                        weight[i] = ok(i)
                            && if symmetric {
                                res[i].abs() <= 2.5 * mad
                            } else {
                                res[i] <= 2.0 * mad && res[i] >= -5.0 * mad
                            };
                    }
                }
            }
            let sol_f: Vec<f32> = sol[..nn].iter().map(|&v| v as f32).collect();
            let mn = sol_f.iter().cloned().fold(f32::MAX, f32::min);
            let mx = sol_f.iter().cloned().fold(f32::MIN, f32::max);
            (sol_f, mx - mn, sol[nn] as f32)
        })
        .collect();
    for (c, (sol, rg, cf)) in per_channel.into_iter().enumerate() {
        for k in 0..nn { data[k * 3 + c] = sol[k]; }
        range[c] = rg;
        coef[c] = cf;
    }
    (Surface { nx, ny, step, data, range }, coef)
}

struct Fit {
    poly: [f32; N_POLY],
    profiles: Vec<Vec<f32>>,
    resid_mad: f32,
    terms: Vec<[f32; N_POLY]>,
    rs: Vec<f32>,
    ths: Vec<f32>,
}

impl Fit {
    /// Value of the fit at cell `i` (same order as the map).
    fn eval_cell(&self, i: usize, r_step: f32) -> f32 {
        let mut m = eval_poly(&self.poly, self.terms[i][1], self.terms[i][2]);
        for p in 0..N_PROF {
            m += eval_radial(&self.profiles[p], self.rs[i], r_step) * harmonic_weight(p, self.ths[i]);
        }
        m
    }
}

#[inline]
fn poly_terms(u: f32, v: f32) -> [f32; N_POLY] {
    [1.0, u, v, u * u, u * v, v * v, u * u * u, u * u * v, u * v * v, v * v * v]
}

#[inline]
fn eval_poly(p: &[f32; N_POLY], u: f32, v: f32) -> f32 {
    let t = poly_terms(u, v);
    let mut s = 0.0;
    for k in 0..N_POLY {
        s += p[k] * t[k];
    }
    s
}

/// Maximum radius with a full ring inside the sensor for a given centre.
fn full_radius(cx: f32, cy: f32, w: usize, h: usize) -> f32 {
    cx.min(cy).min(w as f32 - cx).min(h as f32 - cy).max(1.0)
}

/// Distance from the centre to the farthest corner.
fn corner_radius(cx: f32, cy: f32, w: usize, h: usize) -> f32 {
    let dx = cx.max(w as f32 - cx);
    let dy = cy.max(h as f32 - cy);
    (dx * dx + dy * dy).sqrt()
}

/// Angular harmonics in addition to the pure radial term: `cos(2θ)` and
/// `cos(4θ)` (θ relative to the sensor axes) capture non-circular
/// vignetting (lens hood, filter holder, rectangular sensor) which is
/// indeed fixed to the sensor. `sin` terms omitted: we assume mirror
/// symmetry about the axes.
const HARMONICS: [u32; 2] = [2, 4];
/// Radial profiles per channel: index 0 = pure radial, 1.. = harmonics.
const N_PROF: usize = 1 + HARMONICS.len();
/// Nonzeros per data row: polynomial + 2 nodes per profile.
const ROW_NZ: usize = N_POLY + 2 * N_PROF;

#[inline]
fn harmonic_weight(p: usize, theta: f32) -> f32 {
    if p == 0 {
        1.0
    } else {
        (HARMONICS[p - 1] as f32 * theta).cos()
    }
}

/// Fits `poly + Σ_p f_p(r)·cos(h_p θ)` to one channel of the map for a
/// given centre as a single least squares problem: unknowns = polynomial
/// coefficients + values of each profile at nodes every `r_step`
/// (piecewise linear basis). Regularization: second-difference smoothness
/// on the profiles (rings with few cells — centre and corners — do not
/// inject noise) and a ridge that resolves the polynomial/radial degeneracy
/// and anchors the harmonics to zero where there is no evidence. Robust via
/// IRLS: cells with |residual| > 3·MAD (nebulae, fat stars, Milky Way
/// bands) are discarded in the following iterations.
fn fit_channel(map: &BgMap, c: usize, cx: f32, cy: f32, r_step: f32) -> Fit {
    if map.valid.is_empty() {
        fit_channel_ex(map, c, cx, cy, r_step, None)
    } else {
        let ex = map.invalid_vec();
        fit_channel_ex(map, c, cx, cy, r_step, Some(&ex))
    }
}

/// Like `fit_channel`, with cells excluded a priori (`exclude[i] = true`:
/// foreground, cells without data) that never enter the fit.
fn fit_channel_ex(map: &BgMap, c: usize, cx: f32, cy: f32, r_step: f32, exclude: Option<&[bool]>) -> Fit {
    let n = map.gw * map.gh;
    let w = (map.gw * map.block) as f32;
    let h = (map.gh * map.block) as f32;
    let r_max = corner_radius(cx, cy, map.gw * map.block, map.gh * map.block);
    let n_bins = (r_max / r_step).ceil() as usize + 2;
    let nu = N_POLY + N_PROF * n_bins;
    let prof_off = |p: usize| N_POLY + p * n_bins;

    let mut terms = Vec::with_capacity(n);
    let mut rs = Vec::with_capacity(n);
    let mut ths = Vec::with_capacity(n);
    let mut vals = Vec::with_capacity(n);
    for gy in 0..map.gh {
        for gx in 0..map.gw {
            let (x, y) = map.cell_center(gx, gy);
            terms.push(poly_terms(x / w * 2.0 - 1.0, y / h * 2.0 - 1.0));
            rs.push(((x - cx).powi(2) + (y - cy).powi(2)).sqrt());
            ths.push((y - cy).atan2(x - cx));
            vals.push(map.at(gx, gy, c));
        }
    }

    let mut tmp = vals.clone();
    let level = median_inplace(&mut tmp).abs().max(1e-6);
    let lambda_smooth = 200.0f64; // in "equivalent cells"
    let lambda_ridge = [1e-4f64, 2.0, 2.0]; // pure radial nearly free; harmonics anchored

    let mut poly = [0.0f32; N_POLY];
    let mut profiles: Vec<Vec<f32>> = vec![vec![0.0f32; n_bins]; N_PROF];
    let mut weight: Vec<bool> = match exclude {
        Some(e) => e.iter().map(|&x| !x).collect(),
        None => vec![true; n],
    };
    let mut resid_mad = 0.0f32;

    let mut ata = vec![0.0f64; nu * nu];
    let mut atb = vec![0.0f64; nu];
    for iter in 0..5 {
        for v in ata.iter_mut() { *v = 0.0; }
        for v in atb.iter_mut() { *v = 0.0; }
        let mut idx = [0usize; ROW_NZ];
        let mut val = [0.0f64; ROW_NZ];
        for i in 0..n {
            if !weight[i] {
                continue;
            }
            let t = rs[i] / r_step;
            let k = (t.floor() as usize).min(n_bins - 2);
            let f = (t - k as f32).clamp(0.0, 1.0) as f64;
            for a in 0..N_POLY {
                idx[a] = a;
                val[a] = terms[i][a] as f64;
            }
            for p in 0..N_PROF {
                let hw = harmonic_weight(p, ths[i]) as f64;
                idx[N_POLY + 2 * p] = prof_off(p) + k;
                val[N_POLY + 2 * p] = (1.0 - f) * hw;
                idx[N_POLY + 2 * p + 1] = prof_off(p) + k + 1;
                val[N_POLY + 2 * p + 1] = f * hw;
            }
            let y = vals[i] as f64;
            for a in 0..ROW_NZ {
                let ia = idx[a];
                atb[ia] += val[a] * y;
                let row = &mut ata[ia * nu..(ia + 1) * nu];
                for b in 0..ROW_NZ {
                    row[idx[b]] += val[a] * val[b];
                }
            }
        }
        for p in 0..N_PROF {
            let off = prof_off(p);
            for k in 1..n_bins - 1 {
                let ids = [off + k - 1, off + k, off + k + 1];
                let cf = [1.0, -2.0, 1.0];
                for a in 0..3 {
                    for b in 0..3 {
                        ata[ids[a] * nu + ids[b]] += lambda_smooth * cf[a] * cf[b];
                    }
                }
            }
            for k in 0..n_bins {
                ata[(off + k) * nu + off + k] += lambda_ridge[p];
            }
        }
        let Some(sol) = solve_dense(&mut ata, &mut atb, nu) else {
            break;
        };
        for a in 0..N_POLY {
            poly[a] = sol[a] as f32;
        }
        for p in 0..N_PROF {
            for k in 0..n_bins {
                profiles[p][k] = sol[prof_off(p) + k] as f32;
            }
        }
        // --- robust cell rejection: |res| > 3·MAD ---
        let res: Vec<f32> = (0..n)
            .map(|i| {
                let mut m = eval_poly(&poly, terms[i][1], terms[i][2]);
                for p in 0..N_PROF {
                    m += eval_radial(&profiles[p], rs[i], r_step) * harmonic_weight(p, ths[i]);
                }
                vals[i] - m
            })
            .collect();
        let mut abs: Vec<f32> = (0..n)
            .filter(|&i| exclude.map_or(true, |e| !e[i]))
            .map(|i| res[i].abs())
            .collect();
        resid_mad = median_inplace(&mut abs);
        if iter < 4 && resid_mad > 0.0 {
            let thr = 3.0 * resid_mad.max(1e-4 * level);
            for i in 0..n {
                weight[i] = res[i].abs() <= thr && exclude.map_or(true, |e| !e[i]);
            }
        }
    }
    Fit { poly, profiles, resid_mad, terms, rs, ths }
}

/// Gauss with partial pivoting over a dense `n×n` matrix (row-major).
/// Destroys `a` and `b`.
fn solve_dense(a: &mut [f64], b: &mut [f64], n: usize) -> Option<Vec<f64>> {
    for col in 0..n {
        let mut piv = col;
        for r in col + 1..n {
            if a[r * n + col].abs() > a[piv * n + col].abs() {
                piv = r;
            }
        }
        if a[piv * n + col].abs() < 1e-18 {
            return None;
        }
        if piv != col {
            for k in 0..n {
                a.swap(col * n + k, piv * n + k);
            }
            b.swap(col, piv);
        }
        let d = a[col * n + col];
        for r in col + 1..n {
            let f = a[r * n + col] / d;
            if f == 0.0 {
                continue;
            }
            for k in col..n {
                a[r * n + k] -= f * a[col * n + k];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0f64; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for k in r + 1..n {
            s -= a[r * n + k] * x[k];
        }
        x[r] = s / a[r * n + r];
    }
    Some(x)
}

/// Interpolates the profile; beyond the last bin it stays constant.
#[inline]
fn eval_radial(prof: &[f32], r: f32, r_step: f32) -> f32 {
    let t = r / r_step;
    let i = t.floor() as usize;
    if i + 1 >= prof.len() {
        return *prof.last().unwrap_or(&0.0);
    }
    let f = t - i as f32;
    prof[i] * (1.0 - f) + prof[i + 1] * f
}

impl GlareModel {
    /// Fits the model to the map: searches for the optical centre (minimum
    /// robust residual in G, coarse→fine search within the central ±25%)
    /// and fits polynomial + profile per channel.
    pub fn fit(map: &BgMap, width: usize, height: usize, with_surface: bool) -> Self {
        let r_step = (map.block * 2) as f32;
        let w = width as f32;
        let h = height as f32;

        let mut best = (w * 0.5, h * 0.5, f32::INFINITY);
        let mut half_x = w * 0.25;
        let mut half_y = h * 0.25;
        let mut steps = 5i32;
        for _level in 0..4 {
            let (bx, by, _) = best;
            let cands: Vec<(f32, f32)> = (-steps..=steps)
                .flat_map(|iy| {
                    (-steps..=steps).map(move |ix| {
                        (
                            bx + ix as f32 * half_x / steps as f32,
                            by + iy as f32 * half_y / steps as f32,
                        )
                    })
                })
                .collect();
            let scored: Vec<(f32, f32, f32)> = cands
                .par_iter()
                .map(|&(cx, cy)| (cx, cy, fit_channel(map, 1, cx, cy, r_step).resid_mad))
                .collect();
            for s in scored {
                if s.2 < best.2 {
                    best = s;
                }
            }
            half_x /= steps as f32;
            half_y /= steps as f32;
            steps = 3;
        }
        let (cx, cy, _) = best;

        let fits: Vec<Fit> = (0..3)
            .into_par_iter()
            .map(|c| fit_channel(map, c, cx, cy, r_step))
            .collect();

        let mut model = GlareModel {
            cx,
            cy,
            poly: [[0.0; N_POLY]; 3],
            profiles: [Vec::new(), Vec::new(), Vec::new()],
            r_step,
            r_full: full_radius(cx, cy, map.gw * map.block, map.gh * map.block),
            pedestal: [0.0; 3],
            residual_mad: [0.0; 3],
            surface: None,
            surface_coarse: None,
            lines: None,
            horizon: None,
            width,
            height,
        };
        for (c, f) in fits.into_iter().enumerate() {
            model.poly[c] = f.poly;
            model.profiles[c] = f.profiles;
            model.residual_mad[c] = f.resid_mad;
        }
        if with_surface {
            model.surface = Some(fit_residual_surface(map, &model, SURF_STEP_CELLS, SURF_LAMBDA));
            model.surface_coarse = model.surface.as_ref().map(|s| smooth_surface(s, SURF_COARSE_SIGMA_PX));
            model.lines = Some(fit_fixed_lines(map, &model));
            model.horizon = Some(fit_horizon(map, &model));
        }
        for c in 0..3 {
            let mut v: Vec<f32> = (0..map.gh)
                .flat_map(|gy| (0..map.gw).map(move |gx| (gx, gy)))
                .filter(|&(gx, gy)| map.is_valid(gy * map.gw + gx))
                .map(|(gx, gy)| {
                    let (x, y) = map.cell_center(gx, gy);
                    model.eval_raw(x, y, c)
                })
                .collect();
            model.pedestal[c] = median_inplace(&mut v);
        }
        model
    }

    /// Value of the model (polynomial + radial) without subtracting the pedestal.
    #[inline]
    fn eval_raw(&self, x: f32, y: f32, c: usize) -> f32 {
        self.eval_raw_rgb(x, y)[c]
    }

    /// Model over the three channels (geometry computed only once).
    #[inline]
    fn eval_raw_rgb(&self, x: f32, y: f32) -> [f32; 3] {
        let dx = x - self.cx;
        let dy = y - self.cy;
        let r = (dx * dx + dy * dy).sqrt();
        let theta = dy.atan2(dx);
        let u = x / self.width as f32 * 2.0 - 1.0;
        let v = y / self.height as f32 * 2.0 - 1.0;
        let mut hw = [0.0f32; N_PROF];
        for p in 0..N_PROF {
            hw[p] = harmonic_weight(p, theta);
        }
        let mut out = [0.0f32; 3];
        for c in 0..3 {
            let mut m = eval_poly(&self.poly[c], u, v);
            for p in 0..N_PROF {
                m += eval_radial(&self.profiles[c][p], r, self.r_step) * hw[p];
            }
            out[c] = m;
        }
        if let Some(sf) = &self.surface {
            // The surface only acts towards the edges (r > 0.5·r_full,
            // smooth ramp up to 0.9·r_full): that is where flare, wedges and
            // the lens hood shadow live; at the centre the parametric model
            // already covers the halo, and this way dark Milky Way nebulae
            // are not flattened.
            let t = ((r / self.r_full - SURF_R_START) / (SURF_R_END - SURF_R_START)).clamp(0.0, 1.0);
            let wgt = t * t * (3.0 - 2.0 * t);
            if wgt > 0.0 {
                let sv = sf.eval(x, y);
                // Coarse component: only the part **shared with the
                // mirror** about the vertical axis through the optical
                // centre. System defects (lens hood, flare wedges,
                // mechanical vignetting) are left/right symmetric; the sky
                // is not. From the pair (value, mirrored value) the one with
                // the smaller magnitude is taken if they share a sign;
                // otherwise zero. A dark nebula on one side only is left
                // untouched; a symmetric wedge is removed entirely.
                let cv = match &self.surface_coarse {
                    Some(cs) => {
                        let a = cs.eval(x, y);
                        let b = cs.eval(2.0 * self.cx - x, y);
                        let mut k = [0.0f32; 3];
                        for c in 0..3 {
                            k[c] = if a[c] * b[c] > 0.0 {
                                a[c].signum() * a[c].abs().min(b[c].abs())
                            } else {
                                0.0
                            };
                        }
                        k
                    }
                    None => [0.0; 3],
                };
                let cfull = self.surface_coarse.as_ref().map(|c| c.eval(x, y)).unwrap_or([0.0; 3]);
                for c in 0..3 {
                    // fine − coarse (fine structure) + mirror-shared coarse
                    out[c] += (sv[c] - cfull[c] + cv[c]) * wgt;
                }
            }
        }
        if let Some(hz) = &self.horizon {
            let hv = hz.eval(x - self.cx, y - self.cy);
            for c in 0..3 {
                out[c] += hv[c];
            }
        }
        if let Some(fl) = &self.lines {
            let lv = fl.eval(y, r, theta, self.r_full);
            for c in 0..3 {
                out[c] += lv[c];
            }
        }
        out
    }

    /// RGB correction to subtract at sensor (x, y).
    #[inline]
    pub fn correction_rgb(&self, x: f32, y: f32) -> [f32; 3] {
        let m = self.eval_raw_rgb(x, y);
        [m[0] - self.pedestal[0], m[1] - self.pedestal[1], m[2] - self.pedestal[2]]
    }

    /// Radial profile of channel `c` in % of the pedestal, relative to
    /// r=0, one value per `r_step` px out to `r_full`.
    pub fn radial_profile_pct(&self, c: usize) -> Vec<f32> {
        let ped = self.pedestal[c].max(1e-9);
        let prof = &self.profiles[c][0];
        prof.iter().map(|v| 100.0 * (v - prof[0]) / ped).collect()
    }

    /// Human-readable diagnostic. Amplitudes relative to the pedestal (sky level).
    pub fn report(&self) -> String {
        let c = 1;
        let ped = self.pedestal[c].max(1e-9);
        let prof = &self.profiles[c][0];
        let f_max = prof.iter().cloned().fold(f32::MIN, f32::max);
        let f_min = prof.iter().cloned().fold(f32::MAX, f32::min);
        let ang_max = self.profiles[c][1..]
            .iter()
            .flat_map(|p| p.iter())
            .fold(0.0f32, |m, v| m.max(v.abs()));
        // Range of the polynomial over the sensor (corners and centre).
        let mut p_min = f32::MAX;
        let mut p_max = f32::MIN;
        for (u, v) in [(-1.0, -1.0), (1.0, -1.0), (-1.0, 1.0), (1.0, 1.0), (0.0, 0.0), (0.0, -1.0), (0.0, 1.0), (-1.0, 0.0), (1.0, 0.0)] {
            let p = eval_poly(&self.poly[c], u, v);
            p_min = p_min.min(p);
            p_max = p_max.max(p);
        }
        format!(
            "optical centre ({:.0},{:.0}) [Δ={:+.0},{:+.0} px vs geom. centre]  \
             halo(G): peak−trough {:.2}% (r≤{:.0}px)  angular max {:.2}%  gradient(G): {:.2}%  \
             residual surface(G): {:.2}%  row bands(G): {:.2}%  spokes(G): {:.2}%  horizon(G): {:.2}% towards {:.0}°  residual MAD R/G/B {:.3}/{:.3}/{:.3}%",
            self.cx,
            self.cy,
            self.cx - self.width as f32 * 0.5,
            self.cy - self.height as f32 * 0.5,
            100.0 * (f_max - f_min) / ped,
            self.r_full,
            100.0 * ang_max / ped,
            100.0 * (p_max - p_min) / ped,
            100.0 * self.surface.as_ref().map(|s| s.range[c]).unwrap_or(0.0) / ped,
            100.0 * self.lines.as_ref().map(|l| l.rows_range[c]).unwrap_or(0.0) / ped,
            100.0 * self.lines.as_ref().map(|l| l.spokes_range[c]).unwrap_or(0.0) / ped,
            100.0 * self.horizon.as_ref().map(|h| h.range[c]).unwrap_or(0.0) / ped,
            self.horizon.as_ref().map(|h| h.phi_deg).unwrap_or(0.0),
            100.0 * self.residual_mad[0] / self.pedestal[0].max(1e-9),
            100.0 * self.residual_mad[1] / ped,
            100.0 * self.residual_mad[2] / self.pedestal[2].max(1e-9),
        )
    }
}

/// Image of the removed layer in cropped reference coordinates:
/// `correction(x, y) + pedestal`, that is, exactly what was subtracted from
/// the stack plus the preserved sky level. Useful to check that only
/// optical/sensor defect is removed and not sky structure.
pub fn correction_layer(
    grid: &CorrGrid,
    pedestal: [f32; 3],
    x0: usize,
    y0: usize,
    w: usize,
    h: usize,
) -> Vec<f32> {
    let mut out = vec![0.0f32; w * h * 3];
    out.par_chunks_mut(w * 3).enumerate().for_each(|(ry, row)| {
        for rx in 0..w {
            let k = grid.at(x0 + rx, y0 + ry);
            for c in 0..3 {
                row[rx * 3 + c] = (k[c] + pedestal[c]).max(0.0);
            }
        }
    });
    out
}

/// Single clean frame: applies the frame's WB (balanced scale, like the
/// stack) and subtracts the model correction (plus the frame's own anomaly,
/// `extra`) in sensor coordinates.
///
/// `scatter_comp` (β): contrast compensation. The stray light that forms
/// the halo/veil is **scattered** light: where the defect is strong, a
/// fraction f of the image light has ended up in the veil and the structure
/// (stars, Milky Way) arrives attenuated by (1−f); subtracting the veil
/// leaves that area smoother/with less contrast than the centre.
/// f ≈ β·max(k, 0)/pedestal (veil relative to the sky) is estimated and the
/// deviation from the sky is rescaled by 1/(1−f) (capped at ×2). β = 0
/// disables it. The result is ready for `output::write_dng`. With `clamp`,
/// negative values → 0 (for unbiased stacking it is left unclipped and
/// clipped at the end).
pub fn clean_frame(model: &GlareModel, frame: &Frame, extra: Option<&FrameCorr>, clamp: bool, scatter_comp: f32) -> Vec<f32> {
    let w = frame.width;
    let wb = frame.wb;
    let mut out = vec![0.0f32; frame.rgb.len()];
    let model_block_half = extra.map_or(0.0, |e| e.k_lp.step * 0.5);
    out.par_chunks_mut(w * 3)
        .enumerate()
        .for_each(|(y, row)| {
            let src = &frame.rgb[y * w * 3..(y + 1) * w * 3];
            let yf = y as f32 + 0.5;
            for x in 0..w {
                let xf = x as f32 + 0.5;
                let mut k = model.correction_rgb(xf, yf);
                let mut cg = [1.0f32; 3];
                if scatter_comp > 0.0 {
                    for c in 0..3 {
                        let f = (scatter_comp * k[c].max(0.0) / model.pedestal[c].max(1e-9)).min(0.5);
                        cg[c] = 1.0 / (1.0 - f);
                    }
                }
                if let Some(e) = extra {
                    let sv = e.surface.eval(xf, yf);
                    // k_lp lives on the cell grid (node = cell centre).
                    let hb = model_block_half;
                    let lp = e.k_lp.eval(xf - hb, yf - hb);
                    let rn = ((xf - model.cx).powi(2) + (yf - model.cy).powi(2)).sqrt() / model.r_full;
                    let t = ((rn - FRAME_GAIN_R0) / (FRAME_GAIN_R1 - FRAME_GAIN_R0)).clamp(0.0, 1.0);
                    let wg = t * t * (3.0 - 2.0 * t);
                    for c in 0..3 {
                        k[c] = k[c] + (e.gain[c] - 1.0) * (k[c] - lp[c]) * wg + sv[c] - e.pedestal[c];
                    }
                }
                for c in 0..3 {
                    let mut v = src[x * 3 + c] * wb[c] - k[c];
                    if cg[c] != 1.0 {
                        v = model.pedestal[c] + (v - model.pedestal[c]) * cg[c];
                    }
                    row[x * 3 + c] = if clamp { v.max(0.0) } else { v };
                }
            }
        });
    out
}

/// A frame's own residual correction: **temporal anomaly**. Everything
/// that is constant over the session (sensor defects, median halo, Milky
/// Way and nebulae — with a tracking mount the sky barely moves on the
/// sensor) cancels out in `frame map − session median map`; what is left is
/// what changes over time: horizon glow and twilight (which decay and
/// rotate), variation of the halo amplitude with sky brightness, faint
/// clouds. A smooth bilinear surface is fitted to that anomaly (nodes every
/// `FRAME_SURF_STEP_CELLS` cells, symmetric 2.5·MAD IRLS, no foreground and
/// no cells without data) and subtracted together with the global model.
/// This way every frame's sky ends up like the session median's, and the
/// global model (fitted over that median) removes the defects. `pedestal` =
/// median of the surface over the sky (the level that is preserved).
pub struct FrameCorr {
    pub surface: Surface,
    pub pedestal: [f32; 3],
    pub range: [f32; 3],
    /// Per-channel gain over the model correction (halo, vignetting,
    /// bands, spokes, horizon): defects are stray light and scale with the
    /// light coming in (twilight/dawn: more brightness, more halo). It is
    /// fitted jointly with the anomaly surface (smooth surface +
    /// gain·correction explain the anomaly). It acts only on the fine
    /// structure of the correction (`k − k_lp`) and outside the core (ramp
    /// `FRAME_GAIN_R0..R1`).
    pub gain: [f32; 3],
    /// Smooth component of the model correction (moving average over
    /// `FRAME_SURF_STEP_CELLS` cells), on the cell grid.
    pub k_lp: Surface,
    /// Frame sky level / session median level.
    pub level_ratio: f32,
}

/// Radii (in r_full) of the gain ramp (0 in the core → 1 outside).
const FRAME_GAIN_R0: f32 = 0.20;
const FRAME_GAIN_R1: f32 = 0.35;

/// Node spacing (cells) and smoothness of the per-frame anomaly surface:
/// coarse (256 px), extensive structure only (it does not follow the drift
/// of the Milky Way).
const FRAME_SURF_STEP_CELLS: usize = 8;
const FRAME_SURF_LAMBDA: f64 = 100.0;
/// Reduced smoothness and sky level threshold (frame/median) for the
/// strong twilight/dawn regime.
const FRAME_SURF_LAMBDA_BRIGHT: f64 = 10.0;
const FRAME_BRIGHT_RATIO: f32 = 1.6;

/// Fits the frame's anomaly relative to the session median map (`ref_map`,
/// the very one the model to be applied was fitted on). `mask` = the
/// frame's foreground.
pub fn fit_frame_corr(map: &BgMap, ref_map: &BgMap, model: &GlareModel, mask: Option<&CellMask>) -> FrameCorr {
    let n = map.gw * map.gh;
    let sky = |i: usize| map.is_valid(i) && ref_map.is_valid(i) && mask.map_or(true, |m| !m.data[i]);
    let (gw, gh) = (map.gw, map.gh);
    // Raw anomaly and model correction per cell.
    let mut resid0 = vec![[0.0f32; 3]; n];
    let mut kk = vec![[0.0f32; 3]; n];
    for gy in 0..gh {
        for gx in 0..gw {
            let i = gy * gw + gx;
            let (x, y) = map.cell_center(gx, gy);
            kk[i] = model.correction_rgb(x, y);
            for c in 0..3 { resid0[i][c] = map.data[i * 3 + c] - ref_map.data[i * 3 + c]; }
        }
    }
    // Defect gain fitted jointly with the surface. Only over the **fine
    // structure** of the correction (k minus its moving average over
    // FRAME_SURF_STEP_CELLS cells: halo rings, bands, spokes) and outside
    // the core (ramp FRAME_GAIN_R0→R1 in r_full): the smooth component is
    // followed by the surface anyway, and the central point of the halo
    // does not grow with the horizon light (leaving it at ×1 avoids a dark
    // hole at the centre). Anomaly ≈ surface + a·k_fine.
    let mut k_lp_data = vec![0.0f32; n * 3];
    let mut k_fine = vec![[0.0f32; 3]; n];
    {
        let r = FRAME_SURF_STEP_CELLS;
        let mut integ = vec![[0.0f64; 4]; (gw + 1) * (gh + 1)];
        for gy in 0..gh {
            for gx in 0..gw {
                let i = gy * gw + gx;
                let a = integ[gy * (gw + 1) + gx + 1];
                let b = integ[(gy + 1) * (gw + 1) + gx];
                let d = integ[gy * (gw + 1) + gx];
                let mut cur = [a[0] + b[0] - d[0], a[1] + b[1] - d[1], a[2] + b[2] - d[2], a[3] + b[3] - d[3]];
                if sky(i) {
                    for c in 0..3 { cur[c] += kk[i][c] as f64; }
                    cur[3] += 1.0;
                }
                integ[(gy + 1) * (gw + 1) + gx + 1] = cur;
            }
        }
        for i in 0..n {
            let (gx, gy) = (i % gw, i / gw);
            let x0 = gx.saturating_sub(r);
            let y0 = gy.saturating_sub(r);
            let x1 = (gx + r + 1).min(gw);
            let y1 = (gy + r + 1).min(gh);
            let si = |x: usize, y: usize| integ[y * (gw + 1) + x];
            let (a, b, c2, d) = (si(x1, y1), si(x0, y1), si(x1, y0), si(x0, y0));
            let cnt = a[3] - b[3] - c2[3] + d[3];
            let (x, y) = map.cell_center(gx, gy);
            let rn = ((x - model.cx).powi(2) + (y - model.cy).powi(2)).sqrt() / model.r_full;
            let t = ((rn - FRAME_GAIN_R0) / (FRAME_GAIN_R1 - FRAME_GAIN_R0)).clamp(0.0, 1.0);
            let w = t * t * (3.0 - 2.0 * t);
            for c in 0..3 {
                let sum = a[c] - b[c] - c2[c] + d[c];
                let lp = if cnt > 0.0 { (sum / cnt) as f32 } else { kk[i][c] };
                k_lp_data[i * 3 + c] = lp;
                k_fine[i][c] = (kk[i][c] - lp) * w;
            }
        }
    }
    let k_lp = Surface { nx: gw, ny: gh, step: map.block as f32, data: k_lp_data, range: [0.0; 3] };
    // Frame sky level relative to the median: in strong twilight/dawn
    // (> FRAME_BRIGHT_RATIO) the anomaly is a huge gradient (sky lit from
    // the horizon, ×2–3): the robust fit would reject the edges and the
    // defect gain is not reliable. There: least squares surface without
    // rejection, more flexible, and gain 1 (the deep-night model only).
    let level_ratio = {
        let lv = |m: &BgMap| -> f32 {
            let mut v: Vec<f32> = (0..n).filter(|&i| sky(i)).map(|i| m.data[i * 3 + 1]).collect();
            if v.is_empty() { 0.0 } else { median_inplace(&mut v) }
        };
        let (a, b) = (lv(map), lv(ref_map));
        if b > 0.0 { a / b } else { 1.0 }
    };
    let bright = level_ratio > FRAME_BRIGHT_RATIO;
    let (surface, coef) = if bright {
        fit_surface_core(
            gw, gh, map.block, &resid0, &sky, FRAME_SURF_STEP_CELLS, FRAME_SURF_LAMBDA_BRIGHT, true, None, false,
        )
    } else {
        fit_surface_core(
            gw, gh, map.block, &resid0, &sky, FRAME_SURF_STEP_CELLS, FRAME_SURF_LAMBDA, true, Some(&k_fine), true,
        )
    };
    let mut gain = [1.0f32; 3];
    // Never < 1: the stray light cannot be lower than in the deep-night
    // median (a fit < 1 —twilight sky with the tree covering half the
    // image— would let the already removed veil reappear).
    if !bright {
        for c in 0..3 { gain[c] = (1.0 + coef[c]).clamp(1.0, 3.0); }
    }
    let mut pedestal = [0.0f32; 3];
    let mut range = [0.0f32; 3];
    let mut vals: Vec<[f32; 3]> = Vec::with_capacity(n);
    for gy in 0..map.gh {
        for gx in 0..map.gw {
            if !sky(gy * map.gw + gx) { continue; }
            let (x, y) = map.cell_center(gx, gy);
            vals.push(surface.eval(x, y));
        }
    }
    for c in 0..3 {
        let mut v: Vec<f32> = vals.iter().map(|p| p[c]).collect();
        if v.is_empty() { v.push(0.0); }
        let mx = v.iter().cloned().fold(f32::MIN, f32::max);
        let mn = v.iter().cloned().fold(f32::MAX, f32::min);
        range[c] = mx - mn;
        pedestal[c] = median_inplace(&mut v);
    }
    if std::env::var_os("APILAAA_DEBUG_CORR").is_some() {
        let lv = {
            let mut v: Vec<f32> = (0..n).filter(|&i| sky(i)).map(|i| map.data[i * 3 + 1]).collect();
            if v.is_empty() { 0.0 } else { median_inplace(&mut v) }
        };
        // Variant without gain (as before) for comparison.
        let (surf_old, _) = fit_surface_core(gw, gh, map.block, &resid0, &sky, FRAME_SURF_STEP_CELLS, FRAME_SURF_LAMBDA, true, None, true);
        let ped_old = {
            let mut v: Vec<f32> = (0..n).filter(|&i| sky(i)).map(|i| { let (x, y) = map.cell_center(i % gw, i / gw); surf_old.eval(x, y)[1] }).collect();
            if v.is_empty() { 0.0 } else { median_inplace(&mut v) }
        };
        let mut out = format!("CORRDBG lv={lv:.4} gain={:.3} ped={:.4} ped_old={ped_old:.4}\n", gain[1], pedestal[1]);
        let (ny, nx) = (12usize, 18usize);
        for (name, f) in [
            ("anomaly", Box::new(|i: usize| resid0[i][1]) as Box<dyn Fn(usize) -> f32>),
            ("k_model", Box::new(|i: usize| kk[i][1])),
            ("surface", Box::new(|i: usize| { let (x, y) = map.cell_center(i % gw, i / gw); surface.eval(x, y)[1] - pedestal[1] })),
            ("total", Box::new(|i: usize| { let (x, y) = map.cell_center(i % gw, i / gw); kk[i][1] + (gain[1] - 1.0) * k_fine[i][1] + surface.eval(x, y)[1] - pedestal[1] })),
            ("total_old", Box::new(|i: usize| { let (x, y) = map.cell_center(i % gw, i / gw); kk[i][1] + surf_old.eval(x, y)[1] - ped_old })),
            ("resid_new", Box::new(|i: usize| { let (x, y) = map.cell_center(i % gw, i / gw); resid0[i][1] - (k_fine[i][1] * (gain[1] - 1.0) + surface.eval(x, y)[1]) })),
            ("resid_old", Box::new(|i: usize| { let (x, y) = map.cell_center(i % gw, i / gw); resid0[i][1] - surf_old.eval(x, y)[1] })),
        ] {
            out += &format!("  {name} (% of the sky):\n");
            for j in 0..ny {
                let mut line = String::from("   ");
                for i in 0..nx {
                    let cells: Vec<f32> = (j * gh / ny..(j + 1) * gh / ny)
                        .flat_map(|gy| (i * gw / nx..(i + 1) * gw / nx).map(move |gx| gy * gw + gx))
                        .filter(|&c| sky(c))
                        .map(|c| f(c))
                        .collect();
                    let m = if cells.is_empty() { f32::NAN } else { let mut v = cells; median_inplace(&mut v) };
                    line += &format!("{:+5.0}", 100.0 * m / lv.max(1e-9));
                }
                out += &line;
                out += "\n";
            }
        }
        print!("{out}");
    }
    FrameCorr { surface, pedestal, range, gain, k_lp, level_ratio }
}


/// Diagnostic: temporal medians of the first and second chronological half
/// (`bg_first.csv`, `bg_last.csv`) and the mean transforms of each half
/// (`halves.txt`), to separate a pattern fixed to the sensor from structure
/// fixed to the sky using the drift.
pub fn debug_dump_halves(
    maps: &[BgMap],
    transforms: &[Similarity],
    frame_idx: &[usize],
    dir: &std::path::Path,
) -> anyhow::Result<()> {
    use std::io::Write;
    std::fs::create_dir_all(dir)?;
    let mut order: Vec<usize> = (0..maps.len()).collect();
    order.sort_by_key(|&i| frame_idx[i]);
    let half = order.len() / 2;
    let mut ht = std::fs::File::create(dir.join("halves.txt"))?;
    for (name, sl) in [("bg_first.csv", &order[..half]), ("bg_last.csv", &order[half..])] {
        let sub: Vec<BgMap> = sl.iter().map(|&i| maps[i].clone()).collect();
        let med = temporal_median(&sub);
        let mut f = std::io::BufWriter::new(std::fs::File::create(dir.join(name))?);
        for gy in 0..med.gh {
            for gx in 0..med.gw {
                writeln!(f, "{gx},{gy},{},{},{}", med.at(gx, gy, 0), med.at(gx, gy, 1), med.at(gx, gy, 2))?;
            }
        }
        let n = sl.len() as f32;
        let tx: f32 = sl.iter().map(|&i| transforms[i].tx).sum::<f32>() / n;
        let ty: f32 = sl.iter().map(|&i| transforms[i].ty).sum::<f32>() / n;
        let ang: f32 = sl.iter().map(|&i| transforms[i].angle_deg()).sum::<f32>() / n;
        writeln!(ht, "{name}: n={} tx_mean={tx:.1} ty_mean={ty:.1} ang_mean={ang:.3}", sl.len())?;
    }
    Ok(())
}

/// Diagnostic dump: `bg.csv` (measured map) and `model.csv` (model
/// evaluated at the cells), both `gw` columns × `gh` rows × 3 channels (one
/// CSV row per cell: gx,gy,R,G,B).
pub fn debug_dump(map: &BgMap, model: &GlareModel, dir: &std::path::Path) -> anyhow::Result<()> {
    use std::io::Write;
    std::fs::create_dir_all(dir)?;
    let mut fb = std::io::BufWriter::new(std::fs::File::create(dir.join("bg.csv"))?);
    let mut fm = std::io::BufWriter::new(std::fs::File::create(dir.join("model.csv"))?);
    for gy in 0..map.gh {
        for gx in 0..map.gw {
            let (x, y) = map.cell_center(gx, gy);
            writeln!(fb, "{gx},{gy},{},{},{}", map.at(gx, gy, 0), map.at(gx, gy, 1), map.at(gx, gy, 2))?;
            writeln!(
                fm,
                "{gx},{gy},{},{},{}",
                model.eval_raw(x, y, 0),
                model.eval_raw(x, y, 1),
                model.eval_raw(x, y, 2)
            )?;
        }
    }
    Ok(())
}

/// Average correction in the reference system, on a coarse grid of step
/// `CORR_STEP`. Node (i,j) ↔ ref pixel (i·STEP, j·STEP).
pub struct CorrGrid {
    pub nx: usize,
    pub ny: usize,
    pub step: usize,
    pub data: Vec<f32>, // nx*ny*3
}

impl CorrGrid {
    /// `Σ_i model(m_i(p)) / N` for every node of the grid.
    pub fn build(model: &GlareModel, transforms: &[Similarity], w: usize, h: usize) -> Self {
        let step = CORR_STEP;
        let nx = w / step + 2;
        let ny = h / step + 2;
        let inv_n = 1.0 / transforms.len().max(1) as f32;
        let data: Vec<f32> = (0..ny)
            .into_par_iter()
            .flat_map_iter(|j| {
                let y = (j * step) as f32;
                (0..nx).flat_map(move |i| {
                    let x = (i * step) as f32;
                    let mut acc = [0.0f32; 3];
                    for m in transforms {
                        let (qx, qy) = m.apply(x, y);
                        let k = model.correction_rgb(qx, qy);
                        for c in 0..3 {
                            acc[c] += k[c];
                        }
                    }
                    [acc[0] * inv_n, acc[1] * inv_n, acc[2] * inv_n]
                })
            })
            .collect();
        CorrGrid { nx, ny, step, data }
    }

    /// Correction interpolated bilinearly at reference pixel (x, y).
    #[inline]
    pub fn at(&self, x: usize, y: usize) -> [f32; 3] {
        let fx = x as f32 / self.step as f32;
        let fy = y as f32 / self.step as f32;
        let i0 = (fx.floor() as usize).min(self.nx - 2);
        let j0 = (fy.floor() as usize).min(self.ny - 2);
        let tx = fx - i0 as f32;
        let ty = fy - j0 as f32;
        let idx = |i: usize, j: usize| (j * self.nx + i) * 3;
        let a = idx(i0, j0);
        let b = idx(i0 + 1, j0);
        let c = idx(i0, j0 + 1);
        let d = idx(i0 + 1, j0 + 1);
        let mut out = [0.0f32; 3];
        for k in 0..3 {
            out[k] = self.data[a + k] * (1.0 - tx) * (1.0 - ty)
                + self.data[b + k] * tx * (1.0 - ty)
                + self.data[c + k] * (1.0 - tx) * ty
                + self.data[d + k] * tx * ty;
        }
        out
    }
}
