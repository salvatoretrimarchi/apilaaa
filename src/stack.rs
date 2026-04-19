use crate::align::Similarity;
use crate::flatten::{CellMask, CorrGrid};
use rayon::prelude::*;

pub struct Accumulator {
    pub width: usize,
    pub height: usize,
    pub sum: Vec<f32>,
    /// Accumulated weight per pixel (sky samples; the foreground weighs 0
    /// and its edge is ramped, see `CellMask::sky_weight`).
    pub count: Vec<f32>,
    /// Geometric coverage per pixel: frames whose sensor covers the pixel,
    /// covered or not. Defines the valid crop (`valid_bounds`).
    pub cover: Vec<u32>,
}

impl Accumulator {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            sum: vec![0.0; width * height * 3],
            count: vec![0.0; width * height],
            cover: vec![0; width * height],
        }
    }

    /// Accumulates the image `rgb` (`fw`×`fh`, frame system) into the
    /// reference system by applying the inverse of `m`. `m` is the
    /// transform ref -> cur: q = R·p + t. So for every pixel p of the
    /// reference we sample the current frame at q = R·p + t.
    ///
    /// Applies `wb` while sampling (**balanced** scale: R·G·B aligned on
    /// grey; pass `[1,1,1]` if the image is already balanced). This serves
    /// two purposes: (1) the later stretch analysis sees the three
    /// distributions in the same range and the per-channel cuts are
    /// comparable, (2) the final DNG is written with a neutral
    /// AsShotNeutral, avoiding a double WB application in the developer.
    ///
    /// `mask` (optional, sensor coordinates): pixels whose cell is
    /// foreground count as geometric coverage but do **not** contribute a
    /// sample (so the tree/horizon of some frames does not stain the sky
    /// that other frames do see at that pixel).
    pub fn add(
        &mut self,
        rgb: &[f32],
        fw: usize,
        fh: usize,
        wb: [f32; 3],
        m: &Similarity,
        mask: Option<&CellMask>,
    ) {
        let w = self.width;
        let h = self.height;
        let fwf = fw as f32;
        let fhf = fh as f32;

        let rows: Vec<(usize, Vec<f32>, Vec<f32>)> = (0..h)
            .into_par_iter()
            .map(|y| {
                let mut row_sum = vec![0.0f32; w * 3];
                // weight: <0 = outside, 0 = covered but occluded, >0 = sample
                let mut row_wt = vec![-1.0f32; w];
                for x in 0..w {
                    let (qx, qy) = m.apply(x as f32, y as f32);
                    if qx < 0.0 || qy < 0.0 || qx >= fwf - 1.0 || qy >= fhf - 1.0 {
                        continue;
                    }
                    let wt = match mask {
                        Some(mk) => mk.sky_weight(qx, qy),
                        None => 1.0,
                    };
                    if wt <= 0.0 {
                        row_wt[x] = 0.0;
                        continue;
                    }
                    let x0 = qx.floor() as usize;
                    let y0 = qy.floor() as usize;
                    let dx = qx - x0 as f32;
                    let dy = qy - y0 as f32;
                    let i00 = (y0 * fw + x0) * 3;
                    let i10 = i00 + 3;
                    let i01 = ((y0 + 1) * fw + x0) * 3;
                    let i11 = i01 + 3;
                    let w00 = (1.0 - dx) * (1.0 - dy);
                    let w10 = dx * (1.0 - dy);
                    let w01 = (1.0 - dx) * dy;
                    let w11 = dx * dy;
                    let dst = x * 3;
                    for c in 0..3 {
                        let v = rgb[i00 + c] * w00
                            + rgb[i10 + c] * w10
                            + rgb[i01 + c] * w01
                            + rgb[i11 + c] * w11;
                        row_sum[dst + c] = v * wb[c] * wt;
                    }
                    row_wt[x] = wt;
                }
                (y, row_sum, row_wt)
            })
            .collect();

        for (y, row_sum, row_wt) in rows {
            let row_off = y * w;
            for x in 0..w {
                let wt = row_wt[x];
                if wt < 0.0 { continue; }
                self.cover[row_off + x] += 1;
                if wt == 0.0 { continue; }
                let dst = (row_off + x) * 3;
                let src = x * 3;
                self.sum[dst]     += row_sum[src];
                self.sum[dst + 1] += row_sum[src + 1];
                self.sum[dst + 2] += row_sum[src + 2];
                self.count[row_off + x] += wt;
            }
        }
    }

    /// Largest axis-aligned rectangle within which **every** pixel was
    /// geometrically covered by the maximum number of aligned frames
    /// (occluded by foreground or not). Returns (x0, y0, x1, y1) with
    /// `x1`/`y1` exclusive. If there is no coverage, returns the whole
    /// frame.
    pub fn valid_bounds(&self) -> (usize, usize, usize, usize) {
        let max_c = self.cover.iter().copied().max().unwrap_or(0);
        if max_c == 0 {
            return (0, 0, self.width, self.height);
        }
        let w = self.width;
        let h = self.height;
        let mut heights = vec![0u32; w];
        let mut best_area = 0u64;
        let mut best = (0usize, 0usize, w, h);
        for y in 0..h {
            for x in 0..w {
                if self.cover[y * w + x] == max_c {
                    heights[x] += 1;
                } else {
                    heights[x] = 0;
                }
            }
            let (area, x0, x1, hh) = largest_rect_in_histogram(&heights);
            if area > best_area {
                best_area = area;
                let y1 = y + 1;
                let y0 = y1 - hh as usize;
                best = (x0, y0, x1, y1);
            }
        }
        best
    }

    /// Averages and crops the result to `valid_bounds()` to remove the
    /// edges where the drift left some corner without full coverage.
    /// If `corr` is Some, subtracts the average glare/gradient correction
    /// (see `flatten`) from every pixel; negative values are clamped to 0.
    /// Returns `(rgb, width, height)` of the crop.
    pub fn finalize_cropped(&self, corr: Option<&CorrGrid>) -> (Vec<f32>, usize, usize) {
        let (x0, y0, x1, y1) = self.valid_bounds();
        let new_w = x1 - x0;
        let new_h = y1 - y0;
        let mut out = vec![0.0f32; new_w * new_h * 3];
        out.par_chunks_mut(new_w * 3)
            .enumerate()
            .for_each(|(ry, row)| {
                let y = y0 + ry;
                for x in x0..x1 {
                    let src_idx = y * self.width + x;
                    let c = self.count[src_idx];
                    if c <= 0.0 { continue; }
                    let inv = 1.0 / c;
                    let src = src_idx * 3;
                    let dst = (x - x0) * 3;
                    let k = match corr {
                        Some(g) => g.at(x, y),
                        None => [0.0; 3],
                    };
                    row[dst]     = (self.sum[src]     * inv - k[0]).max(0.0);
                    row[dst + 1] = (self.sum[src + 1] * inv - k[1]).max(0.0);
                    row[dst + 2] = (self.sum[src + 2] * inv - k[2]).max(0.0);
                }
            });
        (out, new_w, new_h)
    }
}

/// Largest rectangle under a histogram (stack algorithm, O(n)). Returns
/// `(area, x0, x1, height)` with `x1` exclusive.
fn largest_rect_in_histogram(heights: &[u32]) -> (u64, usize, usize, u32) {
    let mut stack: Vec<(usize, u32)> = Vec::new();
    let mut best = (0u64, 0usize, 0usize, 0u32);
    let n = heights.len();
    for i in 0..=n {
        let cur = if i == n { 0 } else { heights[i] };
        let mut start = i;
        while let Some(&(prev_i, prev_h)) = stack.last() {
            if prev_h <= cur { break; }
            stack.pop();
            let area = prev_h as u64 * (i - prev_i) as u64;
            if area > best.0 {
                best = (area, prev_i, i, prev_h);
            }
            start = prev_i;
        }
        if stack.last().map_or(true, |&(_, h)| h < cur) {
            stack.push((start, cur));
        }
    }
    best
}
