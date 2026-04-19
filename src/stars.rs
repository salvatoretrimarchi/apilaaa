#[derive(Clone, Copy, Debug)]
pub struct Star {
    pub x: f32,
    pub y: f32,
    pub flux: f32,
}

pub struct Background {
    pub median: f32,
    pub noise: f32,
}

pub fn estimate_background(img: &[f32]) -> Background {
    let mut sample: Vec<f32> = img.iter().step_by(991).copied().collect();
    sample.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if sample.is_empty() {
        return Background { median: 0.0, noise: 1e-6 };
    }
    let median = sample[sample.len() / 2];
    let mut dev: Vec<f32> = sample.iter().map(|v| (v - median).abs()).collect();
    dev.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mad = dev[dev.len() / 2].max(1e-6);
    Background { median, noise: 1.4826 * mad }
}

pub fn detect(img: &[f32], w: usize, h: usize, max_stars: usize) -> Vec<Star> {
    let bg = estimate_background(img);
    let threshold = bg.median + 5.0 * bg.noise;
    let radius: i32 = 2;
    let mut candidates: Vec<Star> = Vec::new();

    for y in (radius as usize)..(h - radius as usize) {
        for x in (radius as usize)..(w - radius as usize) {
            let c = img[y * w + x];
            if c < threshold {
                continue;
            }
            let mut is_max = true;
            'outer: for dy in -radius..=radius {
                for dx in -radius..=radius {
                    if dx == 0 && dy == 0 { continue; }
                    let nx = (x as i32 + dx) as usize;
                    let ny = (y as i32 + dy) as usize;
                    let v = img[ny * w + nx];
                    if v > c || (v == c && (dy < 0 || (dy == 0 && dx < 0))) {
                        is_max = false;
                        break 'outer;
                    }
                }
            }
            if !is_max { continue; }

            let l = img[y * w + (x - 1)];
            let r = img[y * w + (x + 1)];
            let u = img[(y - 1) * w + x];
            let d = img[(y + 1) * w + x];
            let dx_sub = parabola_offset(l, c, r);
            let dy_sub = parabola_offset(u, c, d);

            let mut flux = 0.0;
            for dy in -radius..=radius {
                for dx in -radius..=radius {
                    let nx = (x as i32 + dx) as usize;
                    let ny = (y as i32 + dy) as usize;
                    flux += (img[ny * w + nx] - bg.median).max(0.0);
                }
            }

            candidates.push(Star {
                x: x as f32 + dx_sub,
                y: y as f32 + dy_sub,
                flux,
            });
        }
    }

    candidates.sort_by(|a, b| b.flux.partial_cmp(&a.flux).unwrap());
    candidates.truncate(max_stars);
    candidates
}

fn parabola_offset(l: f32, c: f32, r: f32) -> f32 {
    let denom = l - 2.0 * c + r;
    if denom.abs() < 1e-12 {
        return 0.0;
    }
    (0.5 * (l - r) / denom).clamp(-1.0, 1.0)
}
