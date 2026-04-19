use anyhow::{anyhow, Context, Result};
use rawloader::{CFA, RawImageData};
use std::path::Path;

pub struct Frame {
    pub width: usize,
    pub height: usize,
    pub rgb: Vec<f32>,
    pub wb: [f32; 3],
    pub camera: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Cfa {
    Rggb,
    Bggr,
    Grbg,
    Gbrg,
}

impl Cfa {
    fn from(cfa: &CFA) -> Result<Self> {
        let s: String = (0..2)
            .flat_map(|y| (0..2).map(move |x| cfa.color_at(y, x)))
            .map(|c| match c {
                0 => 'R',
                1 => 'G',
                2 => 'B',
                _ => '?',
            })
            .collect();
        match s.as_str() {
            "RGGB" => Ok(Cfa::Rggb),
            "BGGR" => Ok(Cfa::Bggr),
            "GRBG" => Ok(Cfa::Grbg),
            "GBRG" => Ok(Cfa::Gbrg),
            other => Err(anyhow!("unsupported CFA pattern: {other}")),
        }
    }

    fn is_red(self, x: usize, y: usize) -> bool {
        let (rx, ry) = match self {
            Cfa::Rggb => (0, 0),
            Cfa::Bggr => (1, 1),
            Cfa::Grbg => (1, 0),
            Cfa::Gbrg => (0, 1),
        };
        x % 2 == rx && y % 2 == ry
    }

    fn is_blue(self, x: usize, y: usize) -> bool {
        let (bx, by) = match self {
            Cfa::Rggb => (1, 1),
            Cfa::Bggr => (0, 0),
            Cfa::Grbg => (0, 1),
            Cfa::Gbrg => (1, 0),
        };
        x % 2 == bx && y % 2 == by
    }

}

pub fn load(path: &Path) -> Result<Frame> {
    let raw = rawloader::decode_file(path)
        .with_context(|| format!("decode {}", path.display()))?;
    let cfa = Cfa::from(&raw.cfa)?;
    let w = raw.width;
    let h = raw.height;

    let data: Vec<u16> = match raw.data {
        RawImageData::Integer(v) => v,
        RawImageData::Float(v) => v.into_iter().map(|f| f.clamp(0.0, 65535.0) as u16).collect(),
    };
    if data.len() != w * h {
        return Err(anyhow!("unexpected size: {} != {}*{}", data.len(), w, h));
    }

    let black = raw.blacklevels.iter().map(|&b| b as f32).sum::<f32>() / 4.0;
    let white = raw.whitelevels.iter().map(|&b| b as f32).sum::<f32>() / 4.0;
    let scale = 1.0 / (white - black).max(1.0);

    let mut bayer = vec![0.0f32; w * h];
    for (i, &v) in data.iter().enumerate() {
        bayer[i] = ((v as f32 - black) * scale).clamp(0.0, 1.0);
    }

    let rgb = demosaic_bilinear(&bayer, w, h, cfa);

    let wb = [
        raw.wb_coeffs[0] / raw.wb_coeffs[1].max(1.0),
        1.0,
        raw.wb_coeffs[2] / raw.wb_coeffs[1].max(1.0),
    ];

    let camera = format!("{} {}", raw.make.trim(), raw.model.trim());

    Ok(Frame { width: w, height: h, rgb, wb, camera })
}

fn demosaic_bilinear(bayer: &[f32], w: usize, h: usize, cfa: Cfa) -> Vec<f32> {
    let mut rgb = vec![0.0f32; w * h * 3];
    for y in 0..h {
        for x in 0..w {
            let v = bayer[y * w + x];
            let (mut r, mut g, mut b) = (0.0, 0.0, 0.0);
            let (mut nr, mut ng, mut nb) = (0, 0, 0);
            if cfa.is_red(x, y) {
                r = v; nr = 1;
            } else if cfa.is_blue(x, y) {
                b = v; nb = 1;
            } else {
                g = v; ng = 1;
            }
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    if dx == 0 && dy == 0 { continue; }
                    let nx = x as i32 + dx;
                    let ny = y as i32 + dy;
                    if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 { continue; }
                    let nu = nx as usize;
                    let nv = ny as usize;
                    let s = bayer[nv * w + nu];
                    if cfa.is_red(nu, nv) { r += s; nr += 1; }
                    else if cfa.is_blue(nu, nv) { b += s; nb += 1; }
                    else { g += s; ng += 1; }
                }
            }
            let i = (y * w + x) * 3;
            rgb[i]     = if nr > 0 { r / nr as f32 } else { 0.0 };
            rgb[i + 1] = if ng > 0 { g / ng as f32 } else { 0.0 };
            rgb[i + 2] = if nb > 0 { b / nb as f32 } else { 0.0 };
        }
    }
    rgb
}

pub fn luminance(frame: &Frame) -> Vec<f32> {
    let n = frame.width * frame.height;
    let mut y = vec![0.0f32; n];
    for i in 0..n {
        let j = i * 3;
        y[i] = 0.299 * frame.rgb[j] + 0.587 * frame.rgb[j + 1] + 0.114 * frame.rgb[j + 2];
    }
    y
}
