//! Native EXIF: carries the photographic metadata of the source RAW into
//! the DNG we write, without depending on any external tool.
//!
//! `output::copy_metadata` already does a far more complete job through
//! exiftool — EXIF, MakerNotes, XMP and ICC — but only when exiftool is
//! installed. What this module writes is the subset a raw developer needs
//! to identify the shot on its own: camera make and model, lens make and
//! model, focal length, aperture, exposure, ISO and the capture time. That
//! is exactly what lensfun matches a lens profile on, so distortion,
//! vignetting and chromatic aberration corrections light up in darktable
//! on a bare DNG. exiftool, when present, then adds everything else on top.
//!
//! The DNG is not rewritten: its IFD0 is *relocated*. TIFF puts no
//! constraint on where an IFD lives — the header's first-IFD offset simply
//! points at it — so a new IFD0 carrying the extra tags, the new Exif IFD
//! and their out-of-line values are appended to the end of the file and the
//! header is repointed. The image strips, and every value the original IFD0
//! referenced, stay exactly where they were.

use anyhow::{anyhow, Context, Result};
use std::fs::OpenOptions;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;

/// A tag's value, decoded from the source's byte order so it can be
/// re-encoded in the destination's.
#[derive(Clone, Debug)]
enum Value {
    /// ASCII (2) — the trailing NUL is part of the value.
    Ascii(Vec<u8>),
    /// BYTE (1) or UNDEFINED (7): opaque, byte order does not apply.
    Bytes(u16, Vec<u8>),
    Short(Vec<u16>),
    Long(Vec<u32>),
    Rational(Vec<(u32, u32)>),
    SRational(Vec<(i32, i32)>),
}

impl Value {
    fn type_id(&self) -> u16 {
        match self {
            Value::Ascii(_) => 2,
            Value::Bytes(t, _) => *t,
            Value::Short(_) => 3,
            Value::Long(_) => 4,
            Value::Rational(_) => 5,
            Value::SRational(_) => 10,
        }
    }

    fn count(&self) -> u32 {
        (match self {
            Value::Ascii(v) => v.len(),
            Value::Bytes(_, v) => v.len(),
            Value::Short(v) => v.len(),
            Value::Long(v) => v.len(),
            Value::Rational(v) => v.len(),
            Value::SRational(v) => v.len(),
        }) as u32
    }

    /// The value serialized in `little` byte order, as it goes either
    /// inline in the entry or out of line elsewhere in the file.
    fn bytes(&self, little: bool) -> Vec<u8> {
        let u16b = |v: u16| if little { v.to_le_bytes() } else { v.to_be_bytes() };
        let u32b = |v: u32| if little { v.to_le_bytes() } else { v.to_be_bytes() };
        let i32b = |v: i32| if little { v.to_le_bytes() } else { v.to_be_bytes() };
        let mut out = Vec::new();
        match self {
            Value::Ascii(v) => out.extend_from_slice(v),
            Value::Bytes(_, v) => out.extend_from_slice(v),
            Value::Short(v) => v.iter().for_each(|&x| out.extend_from_slice(&u16b(x))),
            Value::Long(v) => v.iter().for_each(|&x| out.extend_from_slice(&u32b(x))),
            Value::Rational(v) => v.iter().for_each(|&(n, d)| {
                out.extend_from_slice(&u32b(n));
                out.extend_from_slice(&u32b(d));
            }),
            Value::SRational(v) => v.iter().for_each(|&(n, d)| {
                out.extend_from_slice(&i32b(n));
                out.extend_from_slice(&i32b(d));
            }),
        }
        out
    }
}

/// The tags read from the source RAW, split by the IFD they belong in.
/// Photographic tags live in the Exif IFD, not in IFD0: that is where
/// exiv2 — and therefore darktable and lensfun — look for them.
pub struct SourceExif {
    ifd0: Vec<(u16, Value)>,
    exif: Vec<(u16, Value)>,
}

impl SourceExif {
    pub fn is_empty(&self) -> bool {
        self.ifd0.is_empty() && self.exif.is_empty()
    }

    /// Camera and lens as read, for reporting.
    pub fn describe(&self) -> String {
        let ascii = |list: &[(u16, Value)], tag: u16| -> Option<String> {
            list.iter().find(|(t, _)| *t == tag).and_then(|(_, v)| match v {
                Value::Ascii(b) => Some(
                    String::from_utf8_lossy(b).trim_end_matches('\0').trim().to_string(),
                ),
                _ => None,
            })
        };
        let mut parts: Vec<String> = Vec::new();
        if let Some(m) = ascii(&self.ifd0, TAG_MODEL) {
            parts.push(m);
        }
        if let Some(l) = ascii(&self.exif, TAG_LENS_MODEL) {
            parts.push(l);
        }
        if parts.is_empty() { String::from("no identifying tags") } else { parts.join(" + ") }
    }
}

/// When the frame was taken, in seconds, from `DateTimeOriginal` (or
/// `DateTime` if that is all there is).
///
/// The zero point is arbitrary — only differences are ever asked for — so
/// the civil date is turned into days with the usual proleptic-Gregorian
/// count and no timezone is applied. A session that crosses midnight or a
/// DST change is still monotonic in the only sense that matters here.
pub fn capture_seconds(path: &Path) -> Option<f64> {
    let src = read_source(path).ok()?;
    let text = |list: &[(u16, Value)], tag: u16| -> Option<String> {
        list.iter().find(|(t, _)| *t == tag).and_then(|(_, v)| match v {
            Value::Ascii(b) => Some(String::from_utf8_lossy(b).trim_end_matches('\0').trim().to_string()),
            _ => None,
        })
    };
    let s = text(&src.exif, 36867)
        .or_else(|| text(&src.exif, 36868))
        .or_else(|| text(&src.ifd0, 306))?;
    // "YYYY:MM:DD HH:MM:SS"
    let b = s.as_bytes();
    if b.len() < 19 { return None; }
    let num = |a: usize, z: usize| -> Option<i64> { s.get(a..z)?.trim().parse::<i64>().ok() };
    let (y, mo, d) = (num(0, 4)?, num(5, 7)?, num(8, 10)?);
    let (h, mi, se) = (num(11, 13)?, num(14, 16)?, num(17, 19)?);
    // Days from the civil date (the usual days_from_civil).
    let y2 = y - if mo <= 2 { 1 } else { 0 };
    let era = if y2 >= 0 { y2 } else { y2 - 399 } / 400;
    let yoe = y2 - era * 400;
    let doy = (153 * (mo + if mo > 2 { -3 } else { 9 }) + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    let days = era * 146_097 + doe - 719_468;
    Some((days * 86_400 + h * 3600 + mi * 60 + se) as f64)
}

const TAG_MODEL: u16 = 272;
const TAG_LENS_MODEL: u16 = 42036;
const TAG_EXIF_IFD: u16 = 34665;

/// Tags copied into the DNG's IFD0.
const IFD0_TAGS: &[u16] = &[
    271,   // Make
    272,   // Model
    305,   // Software
    306,   // DateTime
    315,   // Artist
    33432, // Copyright
];

/// Tags copied into the DNG's Exif IFD. Everything a lens profile is
/// matched on (lens model, focal length, aperture) plus the exposure
/// triangle and the capture time.
const EXIF_TAGS: &[u16] = &[
    33434, // ExposureTime
    33437, // FNumber
    34850, // ExposureProgram
    34855, // ISOSpeedRatings
    36867, // DateTimeOriginal
    36868, // DateTimeDigitized
    37377, // ShutterSpeedValue
    37378, // ApertureValue
    37380, // ExposureBiasValue
    37381, // MaxApertureValue
    37383, // MeteringMode
    37384, // LightSource
    37385, // Flash
    37386, // FocalLength
    41986, // ExposureMode
    41987, // WhiteBalance
    41989, // FocalLengthIn35mmFilm
    42032, // CameraOwnerName
    42033, // BodySerialNumber
    42034, // LensSpecification
    42035, // LensMake
    42036, // LensModel
    42037, // LensSerialNumber
];

struct Reader<'a> {
    d: &'a [u8],
    little: bool,
}

impl<'a> Reader<'a> {
    fn u16(&self, off: usize) -> Result<u16> {
        let b = self.d.get(off..off + 2).ok_or_else(|| anyhow!("truncated at {off}"))?;
        Ok(if self.little {
            u16::from_le_bytes([b[0], b[1]])
        } else {
            u16::from_be_bytes([b[0], b[1]])
        })
    }

    fn u32(&self, off: usize) -> Result<u32> {
        let b = self.d.get(off..off + 4).ok_or_else(|| anyhow!("truncated at {off}"))?;
        Ok(if self.little {
            u32::from_le_bytes([b[0], b[1], b[2], b[3]])
        } else {
            u32::from_be_bytes([b[0], b[1], b[2], b[3]])
        })
    }

    fn i32(&self, off: usize) -> Result<i32> {
        Ok(self.u32(off)? as i32)
    }

    /// Decodes one IFD entry's value. Returns None for a type we do not
    /// carry across rather than failing the whole read.
    fn value(&self, typ: u16, count: u32, field_off: usize) -> Result<Option<Value>> {
        let elem = match typ {
            1 | 2 | 6 | 7 => 1usize,
            3 | 8 => 2,
            4 | 9 | 11 => 4,
            5 | 10 | 12 => 8,
            _ => return Ok(None),
        };
        let total = elem * count as usize;
        // Up to four bytes live inline in the entry's value field.
        let base = if total <= 4 { field_off } else { self.u32(field_off)? as usize };
        if self.d.len() < base + total {
            return Ok(None);
        }
        let v = match typ {
            2 => Value::Ascii(self.d[base..base + total].to_vec()),
            1 | 7 => Value::Bytes(typ, self.d[base..base + total].to_vec()),
            3 => Value::Short((0..count as usize).map(|i| self.u16(base + 2 * i)).collect::<Result<_>>()?),
            4 => Value::Long((0..count as usize).map(|i| self.u32(base + 4 * i)).collect::<Result<_>>()?),
            5 => Value::Rational(
                (0..count as usize)
                    .map(|i| Ok((self.u32(base + 8 * i)?, self.u32(base + 8 * i + 4)?)))
                    .collect::<Result<_>>()?,
            ),
            10 => Value::SRational(
                (0..count as usize)
                    .map(|i| Ok((self.i32(base + 8 * i)?, self.i32(base + 8 * i + 4)?)))
                    .collect::<Result<_>>()?,
            ),
            _ => return Ok(None),
        };
        Ok(Some(v))
    }

    /// Collects the requested tags of the IFD at `off`. `exif_ptr`, when
    /// given, receives the offset of the Exif IFD if this IFD points to one.
    fn collect(&self, off: usize, want: &[u16], out: &mut Vec<(u16, Value)>, exif_ptr: &mut Option<usize>) -> Result<()> {
        let n = self.u16(off)? as usize;
        for i in 0..n {
            let e = off + 2 + 12 * i;
            let tag = self.u16(e)?;
            let typ = self.u16(e + 2)?;
            let count = self.u32(e + 4)?;
            if tag == TAG_EXIF_IFD {
                *exif_ptr = Some(self.u32(e + 8)? as usize);
                continue;
            }
            if !want.contains(&tag) {
                continue;
            }
            if let Some(v) = self.value(typ, count, e + 8)? {
                out.push((tag, v));
            }
        }
        Ok(())
    }
}

/// Reads the identifying and photographic tags from a source RAW.
pub fn read_source(path: &Path) -> Result<SourceExif> {
    // The tags we want all sit in IFD0 and the Exif IFD, both near the head
    // of the file; the values they point at are close behind. Reading the
    // first few MiB avoids pulling a 50 MiB raw into memory for a handful
    // of strings.
    const HEAD: usize = 4 << 20;
    let mut f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut d = vec![0u8; HEAD];
    let n = read_upto(&mut f, &mut d)?;
    d.truncate(n);
    if d.len() < 8 {
        return Err(anyhow!("{}: too short to be a TIFF", path.display()));
    }
    let little = match &d[0..2] {
        b"II" => true,
        b"MM" => false,
        _ => return Err(anyhow!("{}: not a TIFF/RAW container", path.display())),
    };
    let r = Reader { d: &d, little };
    let ifd0_off = r.u32(4)? as usize;
    let mut ifd0 = Vec::new();
    let mut exif = Vec::new();
    let mut exif_ptr = None;
    r.collect(ifd0_off, IFD0_TAGS, &mut ifd0, &mut exif_ptr)?;
    if let Some(p) = exif_ptr {
        let mut ignored = None;
        r.collect(p, EXIF_TAGS, &mut exif, &mut ignored)?;
    }
    Ok(SourceExif { ifd0, exif })
}

fn read_upto(f: &mut std::fs::File, buf: &mut [u8]) -> Result<usize> {
    let mut n = 0usize;
    while n < buf.len() {
        match f.read(&mut buf[n..])? {
            0 => break,
            k => n += k,
        }
    }
    Ok(n)
}

/// Writes `src`'s tags into the DNG at `path` by relocating its IFD0 to the
/// end of the file, with the Exif IFD appended alongside. Tags the DNG
/// already carries are left alone: what the writer chose about the file's
/// own structure and colour always wins over what the source RAW said.
pub fn embed(path: &Path, src: &SourceExif) -> Result<()> {
    if src.is_empty() {
        return Ok(());
    }
    let mut f = OpenOptions::new()
        .read(true)
        .write(true)
        .open(path)
        .with_context(|| format!("opening {}", path.display()))?;
    let mut header = [0u8; 8];
    f.read_exact(&mut header)?;
    let little = match &header[0..2] {
        b"II" => true,
        b"MM" => false,
        _ => return Err(anyhow!("{}: not a TIFF", path.display())),
    };
    let rd_u16 = |b: &[u8]| if little { u16::from_le_bytes([b[0], b[1]]) } else { u16::from_be_bytes([b[0], b[1]]) };
    let rd_u32 = |b: &[u8]| {
        if little {
            u32::from_le_bytes([b[0], b[1], b[2], b[3]])
        } else {
            u32::from_be_bytes([b[0], b[1], b[2], b[3]])
        }
    };
    let wr_u16 = |v: u16| if little { v.to_le_bytes() } else { v.to_be_bytes() };
    let wr_u32 = |v: u32| if little { v.to_le_bytes() } else { v.to_be_bytes() };

    let ifd0_off = rd_u32(&header[4..8]) as u64;
    f.seek(SeekFrom::Start(ifd0_off))?;
    let mut cnt = [0u8; 2];
    f.read_exact(&mut cnt)?;
    let n_existing = rd_u16(&cnt) as usize;
    let mut existing = vec![0u8; 12 * n_existing];
    f.read_exact(&mut existing)?;
    let mut next_ifd = [0u8; 4];
    f.read_exact(&mut next_ifd)?;

    let present: Vec<u16> = (0..n_existing).map(|i| rd_u16(&existing[12 * i..])).collect();
    let new_ifd0: Vec<&(u16, Value)> = src.ifd0.iter().filter(|(t, _)| !present.contains(t)).collect();
    let have_exif = !src.exif.is_empty() && !present.contains(&TAG_EXIF_IFD);
    if new_ifd0.is_empty() && !have_exif {
        return Ok(());
    }

    // An IFD entry holds its value inline when it fits in four bytes;
    // otherwise the entry holds an offset and the value lives elsewhere.
    // Lay everything out before writing, since the offsets have to be known
    // in advance.
    let inline = |v: &Value| v.bytes(little).len() <= 4;
    let n0 = n_existing + new_ifd0.len() + have_exif as usize;
    let n1 = src.exif.len();
    let base = {
        let end = f.seek(SeekFrom::End(0))?;
        end + (end & 1) // IFDs and their values are word-aligned
    };
    let ifd0_size = 2 + 12 * n0 as u64 + 4;
    let ifd0_vals_off = base + ifd0_size;
    let ifd0_vals_len: u64 = new_ifd0
        .iter()
        .filter(|(_, v)| !inline(v))
        .map(|(_, v)| v.bytes(little).len() as u64)
        .sum();
    let exif_off = ifd0_vals_off + ifd0_vals_len;
    let exif_size = 2 + 12 * n1 as u64 + 4;
    let mut exif_vals_off = exif_off + exif_size;

    // --- new IFD0 ---
    let mut entries: Vec<(u16, Vec<u8>)> = Vec::with_capacity(n0);
    for i in 0..n_existing {
        entries.push((present[i], existing[12 * i..12 * i + 12].to_vec()));
    }
    let mut vals0: Vec<u8> = Vec::new();
    for (tag, v) in &new_ifd0 {
        let raw = v.bytes(little);
        let mut e = Vec::with_capacity(12);
        e.extend_from_slice(&wr_u16(*tag));
        e.extend_from_slice(&wr_u16(v.type_id()));
        e.extend_from_slice(&wr_u32(v.count()));
        if raw.len() <= 4 {
            let mut field = [0u8; 4];
            field[..raw.len()].copy_from_slice(&raw);
            e.extend_from_slice(&field);
        } else {
            e.extend_from_slice(&wr_u32((ifd0_vals_off + vals0.len() as u64) as u32));
            vals0.extend_from_slice(&raw);
        }
        entries.push((*tag, e));
    }
    if have_exif {
        let mut e = Vec::with_capacity(12);
        e.extend_from_slice(&wr_u16(TAG_EXIF_IFD));
        e.extend_from_slice(&wr_u16(4)); // LONG
        e.extend_from_slice(&wr_u32(1));
        e.extend_from_slice(&wr_u32(exif_off as u32));
        entries.push((TAG_EXIF_IFD, e));
    }
    // TIFF requires the entries of an IFD in ascending tag order.
    entries.sort_by_key(|(t, _)| *t);

    // --- Exif IFD ---
    let mut exif_entries: Vec<(u16, Vec<u8>)> = Vec::with_capacity(n1);
    let mut vals1: Vec<u8> = Vec::new();
    if have_exif {
        for (tag, v) in &src.exif {
            let raw = v.bytes(little);
            let mut e = Vec::with_capacity(12);
            e.extend_from_slice(&wr_u16(*tag));
            e.extend_from_slice(&wr_u16(v.type_id()));
            e.extend_from_slice(&wr_u32(v.count()));
            if raw.len() <= 4 {
                let mut field = [0u8; 4];
                field[..raw.len()].copy_from_slice(&raw);
                e.extend_from_slice(&field);
            } else {
                e.extend_from_slice(&wr_u32((exif_vals_off + vals1.len() as u64) as u32));
                vals1.extend_from_slice(&raw);
            }
            exif_entries.push((*tag, e));
        }
        exif_entries.sort_by_key(|(t, _)| *t);
    } else {
        exif_vals_off = exif_off;
    }

    let mut out: Vec<u8> = Vec::new();
    out.extend_from_slice(&wr_u16(entries.len() as u16));
    entries.iter().for_each(|(_, e)| out.extend_from_slice(e));
    out.extend_from_slice(&next_ifd);
    out.extend_from_slice(&vals0);
    if have_exif {
        out.extend_from_slice(&wr_u16(exif_entries.len() as u16));
        exif_entries.iter().for_each(|(_, e)| out.extend_from_slice(e));
        out.extend_from_slice(&wr_u32(0)); // no IFD after the Exif one
        out.extend_from_slice(&vals1);
    }
    debug_assert_eq!(base + out.len() as u64, exif_vals_off + vals1.len() as u64);

    f.seek(SeekFrom::Start(base))?;
    f.write_all(&out)?;
    // Repoint the header at the relocated IFD0. Everything the old one
    // referenced is still where it was, including the image strips.
    f.seek(SeekFrom::Start(4))?;
    f.write_all(&wr_u32(base as u32))?;
    f.flush()?;
    Ok(())
}
