//! Support for **untracked** sequences: the camera stayed fixed on a
//! tripod for the whole timelapse, so the landscape is stationary in
//! sensor coordinates and the sky rotates across it.
//!
//! That inverts the assumption the tracked pipeline is built on. There,
//! the sky is what stays put on the sensor and the defects are separated
//! from it by their being fixed; here both the defects *and* the landscape
//! are fixed, and it is the sky that moves. Two pieces replace the star
//! alignment:
//!
//! - `Template` / `shift` — the residual **tripod drift** (wind, the legs
//!   settling, someone brushing past) measured by normalized
//!   cross-correlation of the high-passed luminance, restricted to the
//!   foreground cells so that the stars, which do move, cannot pull the
//!   estimate towards the sky. Two scales: a coarse search over the whole
//!   drift range, then a fine one with parabolic sub-pixel refinement.
//! - `sky_between` — the transform taking one frame's sky to another's,
//!   composed from the star fit between **consecutive** frames. It is only
//!   ever asked for across the few frames of the export's temporal window,
//!   so the similarity stays a good approximation of the sky's motion and
//!   the composition never accumulates a session-long error.

use crate::align::Similarity;
use crate::flatten::CellMask;
use crate::timelapse::box_blur;
use rayon::prelude::*;

/// Downsampling factors of the two correlation stages.
const COARSE: usize = 8;
const FINE: usize = 2;
/// Radius (in each stage's own px) of the high-pass applied before
/// correlating: only the landscape's edges should drive the match, not the
/// overall brightness, which changes over the session.
const HP_RADIUS: usize = 4;
/// Minimum correlation at the peak for a drift estimate to be accepted.
const MIN_NCC: f32 = 0.25;
/// Minimum fraction of the frame taken up by landscape for the estimate to
/// be attempted at all: a framing that is all sky has nothing fixed to
/// correlate against.
const MIN_FG_FRACTION: f32 = 0.02;
/// Cap on the samples used per correlation stage.
const MAX_SAMPLES: usize = 150_000;

/// A frame's luminance reduced to what the correlation needs, at both
/// scales: box-downsampled and then high-passed.
pub struct Pyramid {
    pub coarse: Vec<f32>,
    pub cw: usize,
    pub ch: usize,
    pub fine: Vec<f32>,
    pub fw: usize,
    pub fh: usize,
}

fn reduce(lum: &[f32], w: usize, h: usize, d: usize) -> (Vec<f32>, usize, usize) {
    let (dw, dh) = (w / d, h / d);
    let mut out = vec![0.0f32; dw * dh];
    out.par_chunks_mut(dw).enumerate().for_each(|(gy, row)| {
        for gx in 0..dw {
            let mut s = 0.0f32;
            for y in gy * d..(gy + 1) * d {
                let off = y * w;
                for x in gx * d..(gx + 1) * d {
                    s += lum[off + x];
                }
            }
            row[gx] = s / (d * d) as f32;
        }
    });
    (out, dw, dh)
}

fn highpass(mut v: Vec<f32>, w: usize, h: usize) -> Vec<f32> {
    let lp = box_blur(&v, w, h, HP_RADIUS);
    v.par_iter_mut().zip(lp.par_iter()).for_each(|(a, b)| *a -= *b);
    v
}

pub fn pyramid(lum: &[f32], w: usize, h: usize) -> Pyramid {
    let (c, cw, ch) = reduce(lum, w, h, COARSE);
    let (f, fw, fh) = reduce(lum, w, h, FINE);
    Pyramid {
        coarse: highpass(c, cw, ch),
        cw,
        ch,
        fine: highpass(f, fw, fh),
        fw,
        fh,
    }
}

/// The reference frame's pyramid plus the sample positions — the landscape
/// cells — every shift is scored over.
pub struct Template {
    pyr: Pyramid,
    coarse_s: Vec<usize>,
    fine_s: Vec<usize>,
    /// Whether there is enough landscape in frame to correlate against.
    pub usable: bool,
    /// Fraction of the frame the landscape takes up (diagnostic).
    pub fg_fraction: f32,
}

pub fn template(lum: &[f32], w: usize, h: usize, mask: &CellMask) -> Template {
    let pyr = pyramid(lum, w, h);
    let pick = |dw: usize, dh: usize, d: usize| -> Vec<usize> {
        let all: Vec<usize> = (0..dw * dh)
            .filter(|i| {
                let x = ((i % dw) * d + d / 2) as f32;
                let y = ((i / dw) * d + d / 2) as f32;
                mask.at_px(x, y)
            })
            .collect();
        let step = (all.len() / MAX_SAMPLES).max(1);
        all.into_iter().step_by(step).collect()
    };
    let coarse_s = pick(pyr.cw, pyr.ch, COARSE);
    let fine_s = pick(pyr.fw, pyr.fh, FINE);
    let fg_fraction = mask.fraction();
    let usable = fg_fraction >= MIN_FG_FRACTION && coarse_s.len() >= 64 && fine_s.len() >= 64;
    Template { pyr, coarse_s, fine_s, usable, fg_fraction }
}

/// Normalized cross-correlation of the template against `b` over every
/// integer shift within `radius` of `(cx, cy)`, all in this stage's own px.
/// Returns the peak refined to sub-pixel by a parabola through the three
/// scores around it on each axis, plus the peak correlation.
fn peak(
    a: &[f32],
    b: &[f32],
    w: usize,
    h: usize,
    samples: &[usize],
    cx: i32,
    cy: i32,
    radius: i32,
) -> Option<(f32, f32, f32)> {
    let n = (2 * radius + 1) as usize;
    let scores: Vec<f32> = (0..n * n)
        .into_par_iter()
        .map(|k| {
            let sx = cx + (k % n) as i32 - radius;
            let sy = cy + (k / n) as i32 - radius;
            let (mut num, mut sa, mut sb) = (0.0f64, 0.0f64, 0.0f64);
            for &i in samples {
                let x = (i % w) as i32 + sx;
                let y = (i / w) as i32 + sy;
                if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 {
                    continue;
                }
                let va = a[i] as f64;
                let vb = b[y as usize * w + x as usize] as f64;
                num += va * vb;
                sa += va * va;
                sb += vb * vb;
            }
            let den = (sa * sb).sqrt();
            if den > 0.0 { (num / den) as f32 } else { 0.0 }
        })
        .collect();
    let (kb, &best) = scores
        .iter()
        .enumerate()
        .max_by(|p, q| p.1.partial_cmp(q.1).unwrap())?;
    if !best.is_finite() {
        return None;
    }
    let (bi, bj) = ((kb % n) as i32, (kb / n) as i32);
    let at = |i: i32, j: i32| scores[j as usize * n + i as usize];
    // Vertex of the parabola through (−1, l), (0, m), (1, r).
    let vertex = |m: f32, l: f32, r: f32| -> f32 {
        let d = 2.0 * (2.0 * m - l - r);
        if d.abs() < 1e-12 { 0.0 } else { ((r - l) / d).clamp(-1.0, 1.0) }
    };
    let ox = if bi > 0 && bi + 1 < n as i32 { vertex(best, at(bi - 1, bj), at(bi + 1, bj)) } else { 0.0 };
    let oy = if bj > 0 && bj + 1 < n as i32 { vertex(best, at(bi, bj - 1), at(bi, bj + 1)) } else { 0.0 };
    Some((
        (cx + bi - radius) as f32 + ox,
        (cy + bj - radius) as f32 + oy,
        best,
    ))
}

/// Tripod drift of `cur` with respect to the template, in sensor px: the
/// landscape the reference shows at `p` is found at `p + (dx, dy)` in this
/// frame, which is exactly the convention the rest of the pipeline uses for
/// the ref → frame similarity. Also returns the peak correlation. None when
/// there is not enough landscape or the peak is too weak to trust — the
/// caller then falls back to identity, which for a tripod is the honest
/// answer.
pub fn shift(tpl: &Template, cur: &Pyramid, search: usize) -> Option<(f32, f32, f32)> {
    if !tpl.usable {
        return None;
    }
    let rc = (search.div_ceil(COARSE)).max(1) as i32;
    let (cx, cy, _) = peak(
        &tpl.pyr.coarse, &cur.coarse, tpl.pyr.cw, tpl.pyr.ch, &tpl.coarse_s, 0, 0, rc,
    )?;
    // Fine stage: the same peak re-searched at FINE resolution, within one
    // coarse cell of where the coarse stage put it.
    let scale = (COARSE / FINE) as f32;
    let (fx, fy) = ((cx * scale).round() as i32, (cy * scale).round() as i32);
    let rf = (COARSE / FINE) as i32 + 1;
    let (dx, dy, ncc) = peak(
        &tpl.pyr.fine, &cur.fine, tpl.pyr.fw, tpl.pyr.fh, &tpl.fine_s, fx, fy, rf,
    )?;
    if ncc < MIN_NCC {
        return None;
    }
    Some((dx * FINE as f32, dy * FINE as f32, ncc))
}

/// Similarity taking frame `c`'s sensor coordinates to frame `f`'s **for
/// the sky**, composed from the consecutive links (`links[k]` = frame k →
/// frame k+1). None if any link in between is missing, in which case that
/// neighbour simply drops out of the temporal window.
pub fn sky_between(links: &[Option<Similarity>], c: usize, f: usize) -> Option<Similarity> {
    if c == f {
        return Some(Similarity::identity());
    }
    let (a, b) = (c.min(f), c.max(f));
    let mut m = Similarity::identity();
    for k in a..b {
        m = links.get(k)?.as_ref()?.compose(&m);
    }
    Some(if f > c { m } else { m.inverse() })
}
