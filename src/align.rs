use crate::stars::Star;

#[derive(Clone, Copy, Debug)]
pub struct Similarity {
    pub cos_t: f32,
    pub sin_t: f32,
    pub tx: f32,
    pub ty: f32,
}

impl Similarity {
    pub fn identity() -> Self {
        Self { cos_t: 1.0, sin_t: 0.0, tx: 0.0, ty: 0.0 }
    }

    /// Pure translation (no rotation): the residual drift of a fixed
    /// tripod, and the crop offsets composed with it.
    pub fn translation(tx: f32, ty: f32) -> Self {
        Self { cos_t: 1.0, sin_t: 0.0, tx, ty }
    }

    pub fn apply(&self, x: f32, y: f32) -> (f32, f32) {
        (
            self.cos_t * x - self.sin_t * y + self.tx,
            self.sin_t * x + self.cos_t * y + self.ty,
        )
    }

    /// `self ∘ first`: applies `first` and then `self`.
    pub fn compose(&self, first: &Similarity) -> Similarity {
        Similarity {
            cos_t: self.cos_t * first.cos_t - self.sin_t * first.sin_t,
            sin_t: self.sin_t * first.cos_t + self.cos_t * first.sin_t,
            tx: self.cos_t * first.tx - self.sin_t * first.ty + self.tx,
            ty: self.sin_t * first.tx + self.cos_t * first.ty + self.ty,
        }
    }

    /// Inverse transform: `q = R·p + t` → `p = Rᵀ·(q − t)`.
    pub fn inverse(&self) -> Similarity {
        let (c, s) = (self.cos_t, self.sin_t);
        Similarity {
            cos_t: c,
            sin_t: -s,
            tx: -(c * self.tx + s * self.ty),
            ty: -(c * self.ty - s * self.tx),
        }
    }

    pub fn angle_deg(&self) -> f32 {
        self.sin_t.atan2(self.cos_t).to_degrees()
    }
}

struct Triangle {
    i: usize,
    j: usize,
    k: usize,
    inv: (f32, f32),
}

fn build_triangles(stars: &[Star]) -> Vec<Triangle> {
    let n = stars.len();
    let mut tris = Vec::with_capacity(n * (n - 1) * (n - 2) / 6);
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                let mut sides = [
                    (dist(&stars[i], &stars[j]), i, j),
                    (dist(&stars[j], &stars[k]), j, k),
                    (dist(&stars[i], &stars[k]), i, k),
                ];
                sides.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                let c = sides[2].0;
                if c < 1e-3 { continue; }
                let inv = (sides[0].0 / c, sides[1].0 / c);
                let (vi, vj, vk) = order_vertices(&sides);
                tris.push(Triangle { i: vi, j: vj, k: vk, inv });
            }
        }
    }
    tris
}

fn dist(a: &Star, b: &Star) -> f32 {
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    (dx * dx + dy * dy).sqrt()
}

fn order_vertices(sides: &[(f32, usize, usize); 3]) -> (usize, usize, usize) {
    let ids = [sides[0].1, sides[0].2, sides[1].1, sides[1].2, sides[2].1, sides[2].2];
    let mut uniq: Vec<usize> = Vec::new();
    for id in ids {
        if !uniq.contains(&id) { uniq.push(id); }
    }
    let va = find_opposite(&uniq, sides[0].1, sides[0].2);
    let vb = find_opposite(&uniq, sides[1].1, sides[1].2);
    let vc = find_opposite(&uniq, sides[2].1, sides[2].2);
    (va, vb, vc)
}

fn find_opposite(uniq: &[usize], e1: usize, e2: usize) -> usize {
    *uniq.iter().find(|&&u| u != e1 && u != e2).unwrap()
}

pub fn fit(ref_stars: &[Star], cur_stars: &[Star]) -> Option<(Similarity, usize)> {
    fit_ex(ref_stars, cur_stars).map(|(m, n, _)| (m, n))
}

/// Same as `fit`, also returning the RMS residual (px) of the refined
/// similarity over its own inliers. On an untracked sequence that residual
/// is the diagnostic that matters: the sky's apparent motion is a rotation
/// on the celestial sphere, and a similarity only approximates it over a
/// small angle, so a residual that grows past a pixel means the temporal
/// window of the export is reaching too far in time.
pub fn fit_ex(ref_stars: &[Star], cur_stars: &[Star]) -> Option<(Similarity, usize, f32)> {
    if ref_stars.len() < 3 || cur_stars.len() < 3 {
        return None;
    }
    let ref_tris = build_triangles(ref_stars);
    let cur_tris = build_triangles(cur_stars);
    if ref_tris.is_empty() || cur_tris.is_empty() {
        return None;
    }

    let mut votes = vec![vec![0u32; cur_stars.len()]; ref_stars.len()];
    let tol = 0.004f32; // tolerance on the invariants

    for rt in &ref_tris {
        for ct in &cur_tris {
            let d0 = rt.inv.0 - ct.inv.0;
            let d1 = rt.inv.1 - ct.inv.1;
            if d0.abs() < tol && d1.abs() < tol {
                votes[rt.i][ct.i] += 1;
                votes[rt.j][ct.j] += 1;
                votes[rt.k][ct.k] += 1;
            }
        }
    }

    let mut pairs: Vec<(usize, usize, u32)> = Vec::new();
    for i in 0..ref_stars.len() {
        let (best_j, best_v) = votes[i].iter().enumerate().max_by_key(|(_, v)| **v).unwrap();
        if *best_v >= 3 {
            pairs.push((i, best_j, *best_v));
        }
    }

    if pairs.len() < 3 {
        return None;
    }

    // RANSAC similarity from 2 pairs
    let mut best: Option<(Similarity, usize)> = None;
    let k_iter = 200usize;
    let inlier_tol = 3.0f32;
    let seed_pairs: Vec<(Star, Star)> = pairs
        .iter()
        .map(|(i, j, _)| (ref_stars[*i], cur_stars[*j]))
        .collect();

    let mut rng_state: u32 = 0x12345678;
    let n = seed_pairs.len();
    for _ in 0..k_iter {
        let a = next_rand(&mut rng_state) as usize % n;
        let mut b = next_rand(&mut rng_state) as usize % n;
        if b == a { b = (b + 1) % n; }
        if let Some(m) = similarity_from_two(&seed_pairs[a], &seed_pairs[b]) {
            let inliers = count_inliers(&m, &seed_pairs, inlier_tol);
            if best.as_ref().map_or(true, |(_, c)| inliers > *c) {
                best = Some((m, inliers));
            }
        }
    }

    let (m0, inliers0) = best?;
    if inliers0 < 3 {
        return None;
    }

    // Refine with least squares over the inliers
    let inlier_set: Vec<(Star, Star)> = seed_pairs
        .iter()
        .filter(|pq| {
            let (p, q) = **pq;
            let (qx, qy) = m0.apply(p.x, p.y);
            let dx = qx - q.x;
            let dy = qy - q.y;
            (dx * dx + dy * dy).sqrt() < inlier_tol
        })
        .copied()
        .collect();

    let refined = fit_least_squares(&inlier_set)?;
    let rms = {
        let s: f32 = inlier_set
            .iter()
            .map(|(p, q)| {
                let (qx, qy) = refined.apply(p.x, p.y);
                (qx - q.x).powi(2) + (qy - q.y).powi(2)
            })
            .sum();
        (s / inlier_set.len() as f32).sqrt()
    };
    Some((refined, inlier_set.len(), rms))
}

fn next_rand(state: &mut u32) -> u32 {
    *state = state.wrapping_mul(1664525).wrapping_add(1013904223);
    *state
}

fn similarity_from_two(a: &(Star, Star), b: &(Star, Star)) -> Option<Similarity> {
    let (p1, q1) = a;
    let (p2, q2) = b;
    let dpx = p2.x - p1.x;
    let dpy = p2.y - p1.y;
    let dqx = q2.x - q1.x;
    let dqy = q2.y - q1.y;
    let len_p = (dpx * dpx + dpy * dpy).sqrt();
    let len_q = (dqx * dqx + dqy * dqy).sqrt();
    if len_p < 1e-3 || len_q < 1e-3 { return None; }
    let ratio = len_q / len_p;
    if (ratio - 1.0).abs() > 0.05 { return None; }
    let ang_p = dpy.atan2(dpx);
    let ang_q = dqy.atan2(dqx);
    let theta = ang_q - ang_p;
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let tx = q1.x - (cos_t * p1.x - sin_t * p1.y);
    let ty = q1.y - (sin_t * p1.x + cos_t * p1.y);
    Some(Similarity { cos_t, sin_t, tx, ty })
}

fn count_inliers(m: &Similarity, pairs: &[(Star, Star)], tol: f32) -> usize {
    pairs
        .iter()
        .filter(|(p, q)| {
            let (qx, qy) = m.apply(p.x, p.y);
            let dx = qx - q.x;
            let dy = qy - q.y;
            (dx * dx + dy * dy).sqrt() < tol
        })
        .count()
}

fn fit_least_squares(pairs: &[(Star, Star)]) -> Option<Similarity> {
    if pairs.len() < 2 { return None; }
    let n = pairs.len() as f32;
    let mut px = 0.0; let mut py = 0.0;
    let mut qx = 0.0; let mut qy = 0.0;
    for (p, q) in pairs {
        px += p.x; py += p.y;
        qx += q.x; qy += q.y;
    }
    px /= n; py /= n; qx /= n; qy /= n;
    let mut num = 0.0f32;
    let mut den = 0.0f32;
    for (p, q) in pairs {
        let ppx = p.x - px; let ppy = p.y - py;
        let qqx = q.x - qx; let qqy = q.y - qy;
        num += ppx * qqy - ppy * qqx;
        den += ppx * qqx + ppy * qqy;
    }
    let theta = num.atan2(den);
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let tx = qx - (cos_t * px - sin_t * py);
    let ty = qy - (sin_t * px + cos_t * py);
    Some(Similarity { cos_t, sin_t, tx, ty })
}
