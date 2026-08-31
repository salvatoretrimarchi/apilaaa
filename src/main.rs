mod align;
mod exif;
mod fixed;
mod flatten;
mod output;
mod raw;
mod stack;
mod stars;
mod timelapse;
mod ui;

use crate::align::Similarity;
use crate::flatten::AnomalyMode;
use anyhow::{anyhow, Context, Result};
use clap::{Parser, ValueEnum};
use rayon::prelude::*;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::sync_channel;
use std::thread;
use std::time::Instant;

/// What `--version` reports. `CARGO_PKG_VERSION` is the floor the manifest
/// declares; a release build overrides it with `APILAAA_VERSION`, which is
/// the tag CI is publishing. Nothing is written back into `Cargo.toml` or
/// `Cargo.lock` for it, so `--locked` always compares the two files exactly
/// as they were committed.
const VERSION: &str = match option_env!("APILAAA_VERSION") {
    Some(v) => v,
    None => env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(version = VERSION, about = "RAW stacker with translation and rotation correction")]
pub struct Args {
    /// Directory holding the input RAW files (any Bayer format rawloader
    /// reads: .ARW, .CR2, .NEF, .ORF, .RW2, .DNG …)
    #[arg(short, long, default_value = "res")]
    pub input: PathBuf,

    /// Output DNG file
    #[arg(short, long, default_value = "stacked.dng")]
    pub output: PathBuf,

    /// Maximum number of stars used per frame for alignment
    #[arg(long, default_value_t = 40)]
    pub max_stars: usize,

    /// Process at most N frames (useful for testing)
    #[arg(long)]
    pub limit: Option<usize>,

    /// Disable the automatic histogram stretch; keeps the native sensor
    /// scale (useful if you want to do the whole development in darktable
    /// without any prior gain).
    #[arg(long)]
    pub no_stretch: bool,

    /// Disable the automatic removal of the lens halo/glare and of the
    /// background gradient (plane model + radial profile in sensor
    /// coordinates, see `flatten`). Enabled by default.
    #[arg(long)]
    pub no_flatten: bool,

    /// Disable the second correction stage (lower-envelope residual surface
    /// that removes non-radial bands/wedges); keeps only the parametric
    /// halo + gradient model.
    #[arg(long)]
    pub no_residual_surface: bool,

    /// Also write a DNG holding the **removed layer** (halo + gradient +
    /// surface + bands/radii, averaged like the stack, plus the sky
    /// pedestal) with the same crop and stretch as the output, so you can
    /// check in the developer that only defect was removed, not sky.
    #[arg(long, value_name = "DNG")]
    pub dump_correction: Option<PathBuf>,

    /// Skip the stack: do not load, clean and average the selected frames,
    /// and write no output DNG. Only the crop is still worked out, from the
    /// frames' geometry alone, because `--export-clean` needs it. Requires
    /// `--export-clean` (there would otherwise be nothing to produce) and
    /// `--fixed-tripod`: on a tracked sequence the stack is what every
    /// exported frame is levelled against, and there is no substitute for it.
    #[arg(long)]
    pub no_stack: bool,

    /// Directory to export the clean timelapse sequence to: every aligned
    /// frame (including the ones excluded from the stack) with the defect
    /// model + its own gradient and horizon glow subtracted, stabilized to
    /// the reference frame system (same alignment as the stack, same crop),
    /// with deflickering (per-channel median and high percentile matched to
    /// the stack's, foreground excluded) and temporal noise reduction
    /// (sliding window, trimmed mean; the foreground —trees, horizon— is
    /// kept as is, without smearing it with the neighbouring frames). Same
    /// stretch as the stack — except under --fixed-tripod, where the stack
    /// is a star-trail image and both the levels reference and the stretch
    /// are taken from a representative frame instead. Requires the
    /// correction to be enabled.
    #[arg(long, value_name = "DIR")]
    pub export_clean: Option<PathBuf>,

    /// Frames in the temporal window used to reduce noise in --export-clean
    /// (odd; 1 = no reduction). Noise ≈ /√(N−2).
    #[arg(long, default_value_t = 7, value_name = "N")]
    pub export_window: usize,

    /// Do not stabilize in --export-clean (exports in sensor coordinates,
    /// uncropped; temporal noise reduction is disabled too because it
    /// requires aligned frames).
    #[arg(long)]
    pub export_no_stabilize: bool,

    /// Do not apply deflickering in --export-clean.
    #[arg(long)]
    pub export_no_deflicker: bool,

    /// In --export-clean, do NOT preserve transients (meteors, satellites,
    /// planes) over the temporal combination. By default they are kept with
    /// the value from the original frame.
    #[arg(long)]
    pub export_no_transients: bool,

    /// In --export-clean, do NOT treat the sky covered by cloud as the
    /// frame's own. By default a cloud is marked like the foreground is: it
    /// does not enter the temporal window (the cloud is the one thing in
    /// the frame that does not follow the sky, so combining it would show
    /// the stars of the neighbouring frames through it) and it does not
    /// enter the frame's level statistics (a light-polluted cloud is
    /// brighter than any star and takes over the high percentile the
    /// deflickering reads the star cores from).
    #[arg(long)]
    pub export_no_cloud_guard: bool,

    /// Contrast compensation in the halo/veil areas (β): the veil light is
    /// scattered light missing from the structure of that area; after
    /// subtracting the veil, the deviation from the sky is rescaled by
    /// 1/(1 − β·veil/sky) to match the contrast at the centre. 0 = off.
    #[arg(long, default_value_t = 1.0, value_name = "BETA")]
    pub scatter_comp: f32,

    /// Frame selection for the stack: frames whose sky level (median of the
    /// background without foreground) falls outside [median/F, median·F] of
    /// the session (twilight, clouds, moon) are excluded. Excluded frames
    /// are still exported by --export-clean (unless the sky is saturated).
    #[arg(long, default_value_t = 1.6, value_name = "F")]
    pub stack_sky_tolerance: f32,

    /// Frame selection for the stack: frames with more than this fraction of
    /// the frame taken up by foreground (trees, horizon) are excluded.
    /// Defaults to 0: only frames with a clear sky are stacked (frames with
    /// foreground usually also carry the horizon glow, which varies over
    /// time and leaves bands when averaged). With a value > 0, frames below
    /// it are included and only those pixels are masked out. If fewer than
    /// 20 % of the aligned frames are left without foreground, frames are
    /// automatically admitted with a mask up to 60 %. In --fixed-tripod the
    /// default is 1: the landscape is in every frame and is part of the
    /// picture, so it is never a reason to drop a frame.
    #[arg(long, value_name = "FRAC")]
    pub stack_max_foreground: Option<f32>,

    /// **Untracked sequence**: the camera stayed fixed on a tripod for the
    /// whole timelapse, so the landscape is stationary on the sensor and
    /// the stars move across it. The star alignment is replaced by the
    /// measurement of the residual tripod drift, the per-frame anomaly is
    /// restricted so that the drifting Milky Way is not mistaken for a
    /// defect, and the temporal window of --export-clean is combined on
    /// the sky instead of on the sensor. The stack (--output) then comes
    /// out as a star-trail image, which is what averaging an untracked
    /// sequence means.
    #[arg(long)]
    pub fixed_tripod: bool,

    /// In --fixed-tripod, do not measure the tripod drift: assume the
    /// camera was perfectly still and export in sensor coordinates,
    /// uncropped.
    #[arg(long)]
    pub fixed_no_stabilize: bool,

    /// In --fixed-tripod, maximum tripod drift searched for, in sensor px.
    #[arg(long, default_value_t = 64, value_name = "PX")]
    pub fixed_search: usize,

    /// In --fixed-tripod, how much freedom the per-frame temporal anomaly
    /// is given. `coarse` (default) lets it follow only structure far
    /// broader than the Milky Way — horizon glow, twilight, the moon
    /// rising; `none` drops it entirely and applies the defect model
    /// alone; `full` uses the same surface as a tracked sequence, which on
    /// an untracked one will subtract part of the drifting Milky Way.
    #[arg(long, default_value = "coarse", value_name = "MODE")]
    pub fixed_anomaly: AnomalyArg,
}

/// CLI spelling of `flatten::AnomalyMode`.
#[derive(Copy, Clone, PartialEq, Eq, Debug, ValueEnum)]
pub enum AnomalyArg {
    Coarse,
    None,
    Full,
}

impl From<AnomalyArg> for AnomalyMode {
    fn from(a: AnomalyArg) -> Self {
        match a {
            AnomalyArg::Coarse => AnomalyMode::Coarse,
            AnomalyArg::None => AnomalyMode::None,
            AnomalyArg::Full => AnomalyMode::Full,
        }
    }
}

fn main() -> Result<()> {
    // No arguments at all and a terminal on both ends: the setup screen,
    // then the dashboard. Anything else — one flag, a pipe, a CI job — is
    // the command line as it always was, reporting line by line.
    let (args, interactive) = if ui::should_prompt() {
        match ui::wizard::run()? {
            Some(a) => (a, true),
            None => return Ok(()),
        }
    } else {
        (Args::parse(), false)
    };
    if interactive {
        ui::start_dashboard(summary(&args));
    }
    let result = run(args);
    ui::shutdown(result.as_ref().err().map(|e| format!("{e:#}")));
    result
}

/// What the dashboard states about the run at the top of the screen: the
/// handful of settings that decide what comes out of it.
fn summary(args: &Args) -> Vec<(String, String)> {
    let mut v = vec![
        (String::from("input"), args.input.display().to_string()),
        (
            String::from("mode"),
            String::from(if args.fixed_tripod { "untracked (fixed tripod)" } else { "tracked" }),
        ),
    ];
    v.push((
        String::from("output"),
        if args.no_stack { String::from("none (--no-stack)") } else { args.output.display().to_string() },
    ));
    v.push((
        String::from("export"),
        match &args.export_clean {
            Some(d) => format!("{} (window {})", d.display(), args.export_window),
            None => String::from("off"),
        },
    ));
    v.push((String::from("frames"), String::from("counting...")));
    v.push((
        String::from("clean"),
        String::from(if args.no_flatten { "off (--no-flatten)" } else { "halo + gradient" }),
    ));
    v
}

/// Stops between frames when the dashboard's user asks for it. Never in the
/// middle of writing one: a half-written DNG is worse than no DNG.
fn check_abort() -> Result<()> {
    if ui::aborted() {
        return Err(anyhow!("stopped by the user"));
    }
    Ok(())
}

fn run(args: Args) -> Result<()> {
    if args.no_stack {
        if args.export_clean.is_none() {
            return Err(anyhow!("--no-stack requires --export-clean: without either one the run would produce nothing"));
        }
        if !args.fixed_tripod {
            return Err(anyhow!("--no-stack requires --fixed-tripod: on a tracked sequence the stack is the reference every exported frame is levelled against"));
        }
    }
    let mut paths = list_raw(&args.input)?;
    if let Some(n) = args.limit {
        paths.truncate(n);
    }
    if paths.is_empty() {
        return Err(anyhow!(
            "no RAW file found in {} (accepted extensions: {})",
            args.input.display(),
            raw::RAW_EXTENSIONS.join(", ")
        ));
    }
    say!("found {} frames", paths.len());
    ui::head("frames", paths.len().to_string());

    let t0 = Instant::now();
    let ref_frame = raw::load(&paths[0])
        .with_context(|| format!("loading reference {}", paths[0].display()))?;
    say!(
        "reference: {} ({} x {}) in {:.2}s",
        paths[0].file_name().unwrap().to_string_lossy(),
        ref_frame.width,
        ref_frame.height,
        t0.elapsed().as_secs_f32()
    );

    let untracked = args.fixed_tripod;
    let anomaly: AnomalyMode = if untracked { args.fixed_anomaly.into() } else { AnomalyMode::Full };
    // The temporal window of the export needs, on an untracked sequence,
    // the star chain between consecutive frames; without it there is
    // nothing to combine the neighbours on.
    let needs_sky_chain =
        untracked && args.export_clean.is_some() && (args.export_window.max(1) | 1) > 1;

    let lum_ref = raw::luminance(&ref_frame);
    let stars_ref = stars::detect(&lum_ref, ref_frame.width, ref_frame.height, args.max_stars);
    say!("stars in reference: {}", stars_ref.len());
    if stars_ref.len() < 6 {
        if !untracked {
            return Err(anyhow!("too few stars in the reference ({}) — check exposure or max_stars", stars_ref.len()));
        }
        say!("  WARNING: too few stars for the sky-aligned temporal window; --export-window will fall back to 1");
    }

    let camera = ref_frame.camera.clone();
    let w_ref = ref_frame.width;
    let h_ref = ref_frame.height;
    let flatten_on = !args.no_flatten;

    // ---------------------------------------------------------------
    // Pass 1: alignment + background map + foreground mask for every
    // frame. Nothing is stacked yet: first we need to know which frames
    // are usable (sky level, foreground) and fit the defect model.
    // ---------------------------------------------------------------
    let mut infos: Vec<FrameInfo> = Vec::new();
    // Untracked sequence: the reference's landscape, high-passed at two
    // scales, is what every frame's tripod drift is measured against.
    let mut drift_template: Option<fixed::Template> = None;
    {
        let bg = flatten::block_background(&ref_frame);
        let mask = flatten::foreground_mask(&bg);
        let level = flatten::sky_level(&bg, &mask);
        if untracked && !args.fixed_no_stabilize {
            let tpl = fixed::template(&lum_ref, w_ref, h_ref, &mask);
            if tpl.usable {
                say!(
                    "tripod drift: measured against the reference's landscape ({:.1}% of the frame), up to {} px",
                    100.0 * tpl.fg_fraction, args.fixed_search
                );
            } else {
                say!(
                    "tripod drift: NOT measurable (landscape {:.1}% of the frame) — every frame is left at identity",
                    100.0 * tpl.fg_fraction
                );
            }
            drift_template = Some(tpl);
        }
        infos.push(FrameInfo { idx: 0, m: Similarity::identity(), bg, fit_mask: mask.clone(), mask, cloud: 0.0, level, in_stack: true, in_model: true, aligned: true, inliers: usize::MAX, export: true });
    }
    drop(lum_ref);
    drop(ref_frame);
    say!("  [1/{}] reference: foreground {:.1}%", paths.len(), 100.0 * infos[0].mask.fraction());
    // Untracked sequence: each frame's stars, kept so the chain of
    // consecutive fits can be built after this pass.
    let mut star_lists: Vec<Vec<stars::Star>> = vec![Vec::new(); paths.len()];
    if needs_sky_chain {
        star_lists[0] = stars_ref.clone();
    }

    let n_cpus = thread::available_parallelism().map(|n| n.get()).unwrap_or(2);
    let (n_workers, chan_cap, ram_info) = compute_pipeline_params(w_ref, h_ref, n_cpus, paths.len() - 1);
    say!(
        "pipeline: {} workers, up to {} frames in flight (total RAM {:.1} GiB, budget 80% = {:.1} GiB, {:.0} MiB/frame)",
        n_workers,
        n_workers + chan_cap + 1,
        ram_info.total_bytes as f64 / (1u64 << 30) as f64,
        ram_info.budget_bytes as f64 / (1u64 << 30) as f64,
        ram_info.frame_bytes as f64 / (1u64 << 20) as f64,
    );

    let total = paths.len();
    let max_stars = args.max_stars;
    let next_idx = AtomicUsize::new(1);
    let (tx, rx) = sync_channel::<Msg>(chan_cap);
    ui::phase("pass 1: alignment and background maps");
    ui::task_begin("align", "align", total as u64);
    ui::task_set("align", 1);

    let mut aligned = 1usize;
    let mut skipped = 0usize;
    // Frames that loaded but have no reliable alignment (too few stars, a
    // satellite/plane trail confusing the detector): recovered by
    // interpolation.
    let mut pending: Vec<FrameInfo> = Vec::new();

    thread::scope(|s| {
        for _ in 0..n_workers {
            let tx = tx.clone();
            let next_idx = &next_idx;
            let paths = &paths;
            let stars_ref = &stars_ref;
            let drift_template = &drift_template;
            let fixed_search = args.fixed_search;
            s.spawn(move || {
                loop {
                    let idx = next_idx.fetch_add(1, Ordering::Relaxed);
                    if idx >= total || ui::aborted() {
                        break;
                    }
                    let p = &paths[idx];
                    let ti = Instant::now();
                    let msg = match raw::load(p) {
                        Err(e) => Msg::LoadErr {
                            idx,
                            error: format!("{e}"),
                        },
                        Ok(frame) => {
                            if frame.width != w_ref || frame.height != h_ref {
                                Msg::DimMismatch { idx }
                            } else {
                                let lum = raw::luminance(&frame);
                                let stars_cur =
                                    stars::detect(&lum, frame.width, frame.height, max_stars);
                                let n_stars = stars_cur.len();
                                // Tracked: the similarity taking the
                                // reference's sky onto this frame's.
                                // Untracked: the sky is not what to align
                                // on — the landscape is, and only to undo
                                // the tripod's own drift.
                                let (fit, shift) = if untracked {
                                    let s = drift_template.as_ref().and_then(|tpl| {
                                        let pyr = fixed::pyramid(&lum, frame.width, frame.height);
                                        fixed::shift(tpl, &pyr, fixed_search)
                                    });
                                    (None, s)
                                } else {
                                    (align::fit(stars_ref, &stars_cur), None)
                                };
                                drop(lum);
                                let stars = if needs_sky_chain { stars_cur } else { Vec::new() };
                                // Map, mask and level are computed even without
                                // an alignment: the frame may still be recovered
                                // by interpolating its neighbours' alignment.
                                let bg = flatten::block_background(&frame);
                                let mask = flatten::foreground_mask(&bg);
                                let level = flatten::sky_level(&bg, &mask);
                                let bg = Some((bg, mask, level));
                                Msg::Processed {
                                    idx,
                                    n_stars,
                                    fit,
                                    shift,
                                    stars,
                                    bg,
                                    elapsed: ti.elapsed().as_secs_f32(),
                                }
                            }
                        }
                    };
                    if tx.send(msg).is_err() {
                        break;
                    }
                }
            });
        }
        drop(tx);

        while let Ok(msg) = rx.recv() {
            ui::task_add("align", 1);
            match msg {
                Msg::LoadErr { idx, error } => {
                    say!(
                        "  [{}/{}] {}: LOAD FAILED {error}",
                        idx + 1,
                        total,
                        paths[idx].file_name().unwrap().to_string_lossy()
                    );
                    skipped += 1;
                }
                Msg::DimMismatch { idx } => {
                    say!(
                        "  [{}/{}] {}: different dimensions, skipped",
                        idx + 1,
                        total,
                        paths[idx].file_name().unwrap().to_string_lossy()
                    );
                    skipped += 1;
                }
                Msg::Processed {
                    idx,
                    n_stars,
                    shift,
                    stars,
                    bg,
                    elapsed,
                    ..
                } if untracked => {
                    // Untracked: there is no such thing as a frame that
                    // "fails to align". The camera did not move, so the
                    // worst case is identity — the very assumption the
                    // tripod is meant to guarantee.
                    let (bg, mask, level) = bg.expect("map computed for every loaded frame");
                    if needs_sky_chain {
                        star_lists[idx] = stars;
                    }
                    let (m, txt) = match shift {
                        Some((dx, dy, ncc)) => (
                            Similarity::translation(dx, dy),
                            format!("drift ({dx:+.2},{dy:+.2}) px, corr {ncc:.2}"),
                        ),
                        None => (
                            Similarity::identity(),
                            String::from(if args.fixed_no_stabilize { "identity" } else { "drift not measurable, identity" }),
                        ),
                    };
                    infos.push(FrameInfo { idx, m, bg, fit_mask: mask.clone(), mask, cloud: 0.0, level, in_stack: true, in_model: true, aligned: true, inliers: n_stars, export: true });
                    aligned += 1;
                    say!(
                        "  [{}/{}] {}: {} stars, {}, sky {:.4} in {:.2}s",
                        idx + 1,
                        total,
                        paths[idx].file_name().unwrap().to_string_lossy(),
                        n_stars,
                        txt,
                        level,
                        elapsed
                    );
                }
                Msg::Processed {
                    idx,
                    n_stars,
                    fit,
                    bg,
                    elapsed,
                    ..
                } => {
                    match (fit, bg) {
                        (Some((m, inliers)), Some((bg, mask, level))) if inliers >= MIN_INLIERS => {
                            let fg = mask.fraction();
                            infos.push(FrameInfo { idx, m, bg, fit_mask: mask.clone(), mask, cloud: 0.0, level, in_stack: true, in_model: true, aligned: true, inliers, export: true });
                            aligned += 1;
                            say!(
                                "  [{}/{}] {}: {} stars, {} inliers, θ={:.3}°, t=({:.2},{:.2}), sky {:.4}{} in {:.2}s",
                                idx + 1,
                                total,
                                paths[idx].file_name().unwrap().to_string_lossy(),
                                n_stars,
                                inliers,
                                m.angle_deg(),
                                m.tx,
                                m.ty,
                                level,
                                if fg > 0.0 { format!(", foreground {:.1}%", 100.0 * fg) } else { String::new() },
                                elapsed
                            );
                        }
                        (fit, Some((bg, mask, level))) => {
                            let inl = fit.map_or(0, |(_, i)| i);
                            say!(
                                "  [{}/{}] {}: NOT ALIGNED ({} stars detected, {} inliers); will try to interpolate the neighbours' alignment",
                                idx + 1,
                                total,
                                paths[idx].file_name().unwrap().to_string_lossy(),
                                n_stars,
                                inl
                            );
                            pending.push(FrameInfo {
                                idx, m: Similarity::identity(), bg, fit_mask: mask.clone(), mask, cloud: 0.0,
                                level, in_stack: false, in_model: false, aligned: false, inliers: inl, export: true,
                            });
                        }
                        (_, None) => unreachable!(),
                    }
                }
            }
        }
    });
    infos.sort_by_key(|f| f.idx);

    // ---------------------------------------------------------------
    // Alignment plausibility: the tracking drift is smooth, so every
    // frame's similarity must resemble the interpolation of its
    // neighbours. A fit that departs > ALIGN_MAX_DEV_PX at any corner
    // (too few stars, hot pixels or a trail taken for stars) is a
    // failure: it is replaced by the interpolation and the frame does not
    // enter the stack. Frames without an alignment are recovered the same
    // way.
    // ---------------------------------------------------------------
    if !untracked {
        let mut n_fixed = 0usize;
        let mut bad: Vec<usize> = Vec::new();
        for k in 0..infos.len() {
            let Some(pred) = predict_transform(&infos, infos[k].idx, Some(k)) else { continue; };
            let dev = corner_deviation(&infos[k].m, &pred, w_ref, h_ref);
            if dev > ALIGN_MAX_DEV_PX {
                bad.push(k);
                say!(
                    "  {}: implausible alignment (departs {:.1} px from the neighbours' drift, {} inliers) — replaced by interpolation, out of the stack",
                    paths[infos[k].idx].file_name().unwrap().to_string_lossy(), dev, infos[k].inliers
                );
            }
        }
        for &k in &bad {
            if let Some(pred) = predict_transform(&infos, infos[k].idx, Some(k)) {
                infos[k].m = pred;
                infos[k].aligned = false;
                infos[k].in_stack = false;
                n_fixed += 1;
            }
        }
        // Frames without an alignment: interpolate from the reliable neighbours.
        for mut f in pending.drain(..) {
            match predict_transform(&infos, f.idx, None) {
                Some(m) => {
                    say!(
                        "  {}: alignment interpolated from the neighbours (θ={:.3}°, t=({:.2},{:.2})); exported, out of the stack",
                        paths[f.idx].file_name().unwrap().to_string_lossy(), m.angle_deg(), m.tx, m.ty
                    );
                    f.m = m;
                    infos.push(f);
                    n_fixed += 1;
                }
                None => {
                    say!(
                        "  {}: no aligned neighbours nearby, skipped",
                        paths[f.idx].file_name().unwrap().to_string_lossy()
                    );
                    skipped += 1;
                }
            }
        }
        infos.sort_by_key(|f| f.idx);
        if n_fixed > 0 {
            say!("  {} frames with interpolated alignment", n_fixed);
        }
    }

    say!(
        "aligned: {}  interpolated: {}  skipped: {}  time: {:.1}s",
        aligned,
        infos.iter().filter(|f| !f.aligned).count(),
        skipped,
        t0.elapsed().as_secs_f32()
    );
    ui::task_end(
        "align",
        format!("{aligned} aligned, {skipped} skipped, {:.1}s", t0.elapsed().as_secs_f32()),
    );
    check_abort()?;

    // ---------------------------------------------------------------
    // Untracked sequence: one landscape for the whole session, and the
    // chain of consecutive sky fits the export's temporal window runs on.
    // ---------------------------------------------------------------
    let mut sky_links: Vec<Option<Similarity>> = Vec::new();
    if untracked {
        {
            let masks: Vec<flatten::CellMask> = infos.iter().map(|f| f.mask.clone()).collect();
            let consensus = flatten::consensus_mask(&masks);
            drop(masks);
            say!(
                "landscape: single consensus mask over {} frames, {:.1}% of the frame",
                infos.len(),
                100.0 * consensus.fraction()
            );
            for f in infos.iter_mut() {
                f.mask = consensus.clone();
            }
            // What the landscape mask cannot answer: which cells this
            // particular frame had cloud over. The landscape is the same in
            // every frame and the model has to keep it out; a cloud is in
            // one frame and not the next, and the model has to keep it out
            // too — for the opposite reason.
            let tc = Instant::now();
            let maps: Vec<flatten::BgMap> = infos.iter().map(|f| f.bg.clone()).collect();
            let clouds = flatten::cloud_masks(&maps, &consensus);
            if let Some(dir) = std::env::var_os("APILAAA_DEBUG_DIR") {
                flatten::debug_dump_masks(&maps, &consensus, &clouds, Path::new(&dir))?;
            }
            drop(maps);
            let mut n_clouded = 0usize;
            let mut cover: Vec<f32> = Vec::with_capacity(infos.len());
            for (f, c) in infos.iter_mut().zip(clouds) {
                f.cloud = c.fraction();
                if f.cloud > CLOUD_REPORT_FRAC { n_clouded += 1; }
                cover.push(f.cloud);
                f.fit_mask = flatten::union_mask(&f.mask, &c);
                f.level = flatten::sky_level(&f.bg, &f.fit_mask);
            }
            cover.sort_by(|a, b| a.partial_cmp(b).unwrap());
            say!(
                "cloud: {} of {} frames with more than {:.0}% covered (median cover {:.1}%, worst {:.1}%) in {:.1}s — those cells enter no fit",
                n_clouded, cover.len(), 100.0 * CLOUD_REPORT_FRAC,
                100.0 * cover[cover.len() / 2], 100.0 * cover[cover.len() - 1],
                tc.elapsed().as_secs_f32()
            );
        }
        if needs_sky_chain && stars_ref.len() >= 6 {
            let tc = Instant::now();
            let fits: Vec<Option<(Similarity, f32)>> = (0..paths.len().saturating_sub(1))
                .into_par_iter()
                .map(|k| {
                    let (a, b) = (&star_lists[k], &star_lists[k + 1]);
                    if a.len() < 3 || b.len() < 3 {
                        return None;
                    }
                    align::fit_ex(a, b)
                        .and_then(|(m, inl, rms)| (inl >= MIN_INLIERS).then_some((m, rms)))
                })
                .collect();
            let links: Vec<Option<Similarity>> = fits.iter().map(|f| f.map(|(m, _)| m)).collect();
            // The residual of a one-frame link is the honest measure of how
            // far a similarity can stand in for the sky's real motion. Past
            // a pixel it is already bending, and the window should not
            // reach as far in time as it is being asked to.
            let mut rms: Vec<f32> = fits.iter().filter_map(|f| f.map(|(_, r)| r)).collect();
            let med_rms = if rms.is_empty() {
                f32::NAN
            } else {
                let k = rms.len() / 2;
                let (_, m, _) = rms.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
                *m
            };
            let ok = links.iter().filter(|l| l.is_some()).count();
            say!(
                "sky chain: {} of {} consecutive links fitted in {:.1}s (median residual {:.2} px)",
                ok, links.len(), tc.elapsed().as_secs_f32(), med_rms
            );
            if med_rms > 1.0 {
                say!(
                    "  WARNING: a similarity already leaves {:.2} px between consecutive frames — the sky's motion is not one at this field of view or interval. Lower --export-window (or set it to 1) if the stars come out soft.",
                    med_rms
                );
            }
            sky_links = links;
        }
        drop(star_lists);
    }

    // ---------------------------------------------------------------
    // Frame selection for stacking: drop the ones with too much
    // foreground and the ones whose sky is anomalously bright or dark
    // compared to the session (twilight, clouds, moon rising).
    // ---------------------------------------------------------------
    let med_level;
    {
        let mut lv: Vec<f32> = infos.iter().map(|f| f.level).filter(|v| *v > 0.0).collect();
        med_level = if lv.is_empty() {
            0.0
        } else {
            let k = lv.len() / 2;
            let (_, m, _) = lv.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
            *m
        };
        let tol = args.stack_sky_tolerance.max(1.0);
        // Untracked: the landscape is in every frame and is part of the
        // picture, never a reason to drop one.
        let mut max_fg = args
            .stack_max_foreground
            .unwrap_or(if untracked { 1.0 } else { 0.0 })
            .max(0.0);
        let n_clear = infos.iter().filter(|f| f.aligned && f.mask.fraction() <= max_fg).count();
        if n_clear < (infos.len() / 5).max(1) && max_fg < 0.6 {
            say!(
                "  only {} of {} frames have foreground ≤ {:.0}%: frames with foreground (masked) are admitted up to 60%",
                n_clear, infos.len(), 100.0 * max_fg
            );
            max_fg = 0.6;
        }
        let mut n_fg = 0usize;
        let mut n_bright = 0usize;
        let mut n_dark = 0usize;
        let mut n_cloud = 0usize;
        let mut n_with_fg = 0usize;
        let mut n_sat = 0usize;
        for f in infos.iter_mut() {
            if f.mask.any() { n_with_fg += 1; }
            // Everything is exported except an already saturated sky (full
            // dawn): there is no image left to clean there.
            if f.level >= EXPORT_MAX_LEVEL {
                f.export = false;
                n_sat += 1;
            }
            if !f.aligned {
                f.in_stack = false;
            } else if f.mask.fraction() > max_fg {
                f.in_stack = false;
                n_fg += 1;
            } else if f.cloud > CLOUD_MAX_FRAC {
                f.in_stack = false;
                n_cloud += 1;
            } else if med_level > 0.0 && f.level > med_level * tol {
                f.in_stack = false;
                n_bright += 1;
            } else if med_level > 0.0 && f.level < med_level / tol {
                f.in_stack = false;
                n_dark += 1;
            }
        }
        let n_stack = infos.iter().filter(|f| f.in_stack).count();
        say!(
            "stack selection: {} of {} aligned frames (median sky {:.4}, tolerance ×{:.2}); {} with foreground (masked); excluded: {} for foreground > {:.0}%, {} for cloud > {:.0}%, {} for bright sky, {} for dark sky; {} with saturated sky (≥ {:.0}%) are not exported",
            n_stack, infos.len(), med_level, tol, n_with_fg, n_fg, 100.0 * max_fg,
            n_cloud, 100.0 * CLOUD_MAX_FRAC, n_bright, n_dark, n_sat, 100.0 * EXPORT_MAX_LEVEL
        );
        if n_stack == 0 {
            return Err(anyhow!("no frame passes the stack selection (tune --stack-sky-tolerance / --stack-max-foreground)"));
        }
        if std::env::var_os("APILAAA_DEBUG").is_some() {
            for f in infos.iter().filter(|f| !f.in_stack) {
                say!("  excluded from the stack: {} (sky {:.4}, foreground {:.1}%)",
                    paths[f.idx].file_name().unwrap().to_string_lossy(), f.level, 100.0 * f.mask.fraction());
            }
        }
    }

    // ---------------------------------------------------------------
    // Which frames are allowed to say what the background is. Every frame
    // that is stacked, minus the ones with cloud in them: a cloud is the one
    // thing in the sky that looks like a defect from a single frame, and the
    // only way to be sure it never enters the model is to not let a doubtful
    // frame near it.
    // ---------------------------------------------------------------
    {
        for f in infos.iter_mut() {
            f.in_model = f.in_stack && f.cloud <= CLOUD_FIT_FRAC;
        }
        let mut n_model = infos.iter().filter(|f| f.in_model).count();
        let n_candidates = infos.iter().filter(|f| f.in_stack).count();
        if n_model < MODEL_MIN_FRAMES.min(n_candidates) {
            // Every frame of this session has some cloud in it. Take the
            // clearest ones rather than none: a background measured on the
            // least clouded frames of a clouded night is still the best
            // answer available, and it is said out loud.
            let mut order: Vec<usize> = (0..infos.len()).filter(|&i| infos[i].in_stack).collect();
            order.sort_by(|&a, &b| infos[a].cloud.partial_cmp(&infos[b].cloud).unwrap());
            let keep = MODEL_MIN_FRAMES.min(order.len());
            for &i in &order[..keep] { infos[i].in_model = true; }
            n_model = infos.iter().filter(|f| f.in_model).count();
            say!(
                "  WARNING: only {} frames are clear of cloud; the background is determined from the {} least clouded instead (up to {:.0}% cover)",
                infos.iter().filter(|f| f.in_stack && f.cloud <= CLOUD_FIT_FRAC).count(),
                n_model,
                100.0 * order[..keep].iter().map(|&i| infos[i].cloud).fold(0.0, f32::max)
            );
        }
        if untracked {
            say!(
                "background from {} of {} stacked frames (cloud ≤ {:.0}%)",
                n_model, n_candidates, 100.0 * CLOUD_FIT_FRAC
            );
        }
    }

    // ---------------------------------------------------------------
    // Defect model: (masked) temporal median of the maps of the selected
    // frames.
    // ---------------------------------------------------------------
    // The camera's multiplicative field, measured before anything additive
    // is fitted: what scales with the sky level is the lens, and dividing it
    // out is right at every level, which subtracting a fixed amount is not.
    // Every map is divided by it so the model that follows only has to
    // account for what is genuinely added to the frame.
    // Only where the sky sweeps across the sensor. Separating a field fixed
    // to the sensor from the structure of the sky needs the two to move with
    // respect to each other, and a tracked mount exists precisely to stop
    // that: there the regression absorbs the sky instead, reproducibly enough
    // to pass every self-consistency check and still be wrong. Measured on
    // the tracked test session, dividing by the field took the stack's radial
    // residual from 0.56% to 4.65%.
    let flat = if !untracked {
        flatten::FlatField::flat(1, 1, flatten::BLOCK, 1.0)
    } else {
        let sel: Vec<&FrameInfo> = infos.iter().filter(|f| f.in_model).collect();
        let maps: Vec<flatten::BgMap> = sel.iter().map(|f| f.bg.clone()).collect();
        let masks: Vec<flatten::CellMask> = sel.iter().map(|f| f.fit_mask.clone()).collect();
        let levels: Vec<f32> = sel.iter().map(|f| f.level).collect();
        flatten::fit_flat_field(&maps, &masks, &levels)
    };
    if untracked {
        say!("flat field: {}", flat.report());
    }
    if flat.usable {
        for f in infos.iter_mut() {
            let mask = f.fit_mask.clone();
            flat.apply_map(&mut f.bg, &mask);
            f.level = flatten::sky_level(&f.bg, &f.fit_mask);
        }
    }

    // How far the sky turned between the two halves of an untracked
    // session, at half the frame's radius: the lever the residual surface's
    // verification works with. The sky turns once a sidereal day, and a
    // point at radius r from the pole moves r·Δθ; taking r from the image
    // centre understates that whenever the pole is off frame, which is the
    // safe direction to be wrong in.
    let sky_motion_px = if untracked {
        let used: Vec<usize> = infos.iter().filter(|f| f.in_model).map(|f| f.idx).collect();
        match (used.first(), used.last()) {
            (Some(&a), Some(&b)) if b > a => {
                match (exif::capture_seconds(&paths[a]), exif::capture_seconds(&paths[b])) {
                    (Some(ta), Some(tb)) => {
                        // The two halves sit, on average, half a session apart.
                        let half = (tb - ta).abs() * 0.5;
                        let omega = std::f64::consts::TAU / 86_164.0; // sidereal day
                        let r_half = 0.25 * ((w_ref * w_ref + h_ref * h_ref) as f32).sqrt();
                        let px = (half * omega) as f32 * r_half;
                        say!(
                            "sky motion between the halves of the session: {:.0} px at r={:.0} px ({:.0} min apart)",
                            px, r_half, half / 60.0
                        );
                        px
                    }
                    _ => {
                        say!("  no capture time in the frames: the residual surface cannot be verified against the session");
                        0.0
                    }
                }
            }
            _ => 0.0,
        }
    } else {
        0.0
    };

    let stack_infos: Vec<&FrameInfo> = infos.iter().filter(|f| f.in_stack).collect();
    let stack_transforms: Vec<Similarity> = stack_infos.iter().map(|f| f.m).collect();
    // Fits model + median map over a set of frames.
    let fit_model = |sel: &[&FrameInfo], label: &str| -> Result<(flatten::BgMap, flatten::GlareModel)> {
        let _ = label;
        let tf = Instant::now();
        let maps: Vec<flatten::BgMap> = sel.iter().map(|f| f.bg.clone()).collect();
        let masks: Vec<flatten::CellMask> = sel.iter().map(|f| f.fit_mask.clone()).collect();
        let (bg, filled) = flatten::temporal_median_masked(&maps, &masks);
        if filled > 0 {
            say!("  [{label}] {} map cells always covered, filled in from neighbours", filled);
        }
        if let Some(dir) = std::env::var_os("APILAAA_DEBUG_DIR") {
            let idxs: Vec<usize> = sel.iter().map(|f| f.idx).collect();
            let tfs: Vec<Similarity> = sel.iter().map(|f| f.m).collect();
            let sub = Path::new(&dir).join(label);
            fs::create_dir_all(&sub)?;
            flatten::debug_dump_halves(&maps, &tfs, &idxs, &sub)?;
        }
        // The same maps, split down the middle of the session: what both
        // halves see in the same sensor cells is background, and what only
        // one of them sees is sky that moved on.
        let halves = if untracked && maps.len() >= 2 * MIN_HALF_FRAMES && sky_motion_px > 0.0 {
            let h = maps.len() / 2;
            let (a, _) = flatten::temporal_median_masked(&maps[..h], &masks[..h]);
            let (b, _) = flatten::temporal_median_masked(&maps[h..], &masks[h..]);
            Some((a, b))
        } else {
            None
        };
        drop(maps);
        drop(masks);
        let mut model = flatten::GlareModel::fit(
            &bg,
            w_ref,
            h_ref,
            !args.no_residual_surface,
            halves.as_ref().map(|(a, b)| (a, b, sky_motion_px)),
        );
        model.flat = if flat.usable { Some(flat.clone()) } else { None };
        say!("glare/gradient [{label}]: {}", model.report());
        if let Some(dir) = std::env::var_os("APILAAA_DEBUG_DIR") {
            flatten::debug_dump(&bg, &model, &Path::new(&dir).join(label))?;
        }
        if std::env::var_os("APILAAA_DEBUG").is_some() {
            for (c, name) in ["R", "G", "B"].iter().enumerate() {
                let p: Vec<String> = model
                    .radial_profile_pct(c)
                    .iter()
                    .map(|v| format!("{v:+.1}"))
                    .collect();
                say!("  radial profile {name} (every {:.0}px, % vs r=0): {}", model.r_step, p.join(" "));
            }
        }
        say!(
            "  model [{label}] fitted over {} frames in {:.2}s",
            sel.len(),
            tf.elapsed().as_secs_f32()
        );
        Ok((bg, model))
    };
    // A single model, fitted over the median of the selected frames (clear
    // sky, no twilight): it also serves as the reference for the temporal
    // anomaly of every frame in the sequence (`fit_frame_corr_ex`), which
    // absorbs whatever that frame has more or less of compared to that
    // median (horizon glow, twilight, halo amplitude).
    ui::phase("fitting the halo + gradient model");
    let model_infos: Vec<&FrameInfo> = infos.iter().filter(|f| f.in_model).collect();
    let stack_model = if flatten_on { Some(fit_model(&model_infos, "stack")?) } else { None };
    check_abort()?;

    // ---------------------------------------------------------------
    // Pass 2: stacking. Every selected frame is loaded again, cleaned
    // (model + its own gradient) in sensor coordinates and accumulated
    // aligned, skipping the foreground pixels.
    // ---------------------------------------------------------------
    let mut acc = stack::Accumulator::new(w_ref, h_ref);
    if args.no_stack {
        // Only the geometry, which is all `--export-clean` needs from this
        // pass: no frame is loaded, demosaiced or cleaned here.
        let ts = Instant::now();
        for info in &stack_infos {
            acc.add_coverage(w_ref, h_ref, &info.m);
        }
        say!(
            "no stack: coverage of {} frames only, in {:.1}s",
            stack_infos.len(),
            ts.elapsed().as_secs_f32()
        );
        ui::task_begin("stack", "stack", stack_infos.len() as u64);
        ui::task_end("stack", String::from("not written (--no-stack)"));
    } else {
        let ts = Instant::now();
        let n_stack = stack_infos.len();
        say!(
            "stacking {} frames{}...",
            n_stack,
            if untracked { " (untracked: the mean of a fixed sequence is a star-trail image)" } else { "" }
        );
        ui::phase("pass 2: cleaning and stacking");
        ui::task_begin("stack", "stack", n_stack as u64);
        let next_k = AtomicUsize::new(0);
        let scatter_comp = args.scatter_comp;
        let (tx, rx) = sync_channel::<Result<(usize, Vec<f32>)>>(chan_cap);
        let mut done = 0usize;
        thread::scope(|s| {
            for _ in 0..n_workers {
                let tx = tx.clone();
                let next_k = &next_k;
                let paths = &paths;
                let stack_infos = &stack_infos;
                let model = stack_model.as_ref();
                s.spawn(move || {
                    loop {
                        let k = next_k.fetch_add(1, Ordering::Relaxed);
                        if k >= n_stack || ui::aborted() {
                            break;
                        }
                        let info = stack_infos[k];
                        let p = &paths[info.idx];
                        let r = raw::load(p)
                            .with_context(|| format!("loading {}", p.display()))
                            .and_then(|frame| {
                                if frame.width != w_ref || frame.height != h_ref {
                                    return Err(anyhow!("{}: different dimensions", p.display()));
                                }
                                let img = match model {
                                    Some((bg_med, m)) => {
                                        let fc = flatten::fit_frame_corr_ex(&info.bg, bg_med, m, Some(&info.fit_mask), anomaly);
                                        flatten::clean_frame(m, &frame, Some(&fc), false, scatter_comp)
                                    }
                                    None => {
                                        let wb = frame.wb;
                                        let mut v = frame.rgb.clone();
                                        v.chunks_mut(3).for_each(|px| {
                                            for c in 0..3 { px[c] *= wb[c]; }
                                        });
                                        v
                                    }
                                };
                                Ok((k, img))
                            });
                        if tx.send(r).is_err() {
                            break;
                        }
                    }
                });
            }
            drop(tx);
            while let Ok(r) = rx.recv() {
                let (k, img) = r?;
                let info = stack_infos[k];
                // Untracked: the landscape is stationary on the sensor, so
                // masking it would leave it with no sample from any frame at
                // all — a black cutout instead of the ground the trails are
                // meant to stand over. It is stacked like everything else.
                let mask = if untracked { None } else { Some(&info.mask) };
                acc.add(&img, w_ref, h_ref, [1.0; 3], &info.m, mask);
                done += 1;
                ui::task_add("stack", 1);
                if done % 25 == 0 || done == n_stack {
                    say!("  [{}/{}] stacked ({:.1}s)", done, n_stack, ts.elapsed().as_secs_f32());
                }
            }
            Ok::<(), anyhow::Error>(())
        })?;
        say!("stacked: {}  total time: {:.1}s", done, t0.elapsed().as_secs_f32());
        ui::task_end("stack", format!("{done} frames, {:.1}s", ts.elapsed().as_secs_f32()));
        check_abort()?;
    }

    let corr = stack_model
        .as_ref()
        .map(|(_, m)| flatten::CorrGrid::build(m, &stack_transforms, acc.width, acc.height));
    let (out, out_w, out_h) = if args.no_stack {
        let (x0, y0, x1, y1) = acc.valid_bounds();
        (Vec::new(), x1 - x0, y1 - y0)
    } else {
        acc.finalize_cropped(None)
    };
    if out_w != acc.width || out_h != acc.height {
        say!(
            "crop from drift: {}×{} → {}×{} (−{} px height, −{} px width)",
            acc.width, acc.height, out_w, out_h,
            acc.height - out_h, acc.width - out_w
        );
    }
    if !args.no_stack {
        let (x0, y0, _, _) = acc.valid_bounds();
        let mut holes = 0usize;
        let mut min_count = f32::MAX;
        for y in y0..y0 + out_h {
            for x in x0..x0 + out_w {
                let c = acc.count[y * acc.width + x];
                if c <= 0.0 { holes += 1; }
                min_count = min_count.min(c);
            }
        }
        if holes > 0 || min_count < stack_infos.len() as f32 - 0.5 {
            say!(
                "  stack coverage: minimum {:.0} samples/pixel out of {} frames{}",
                min_count,
                stack_infos.len(),
                if holes > 0 { format!("; {} px with no sample at all (always covered)", holes) } else { String::new() }
            );
        }
    }
    let stretch = if args.no_stretch || args.no_stack {
        None
    } else if let (true, Some((bg_med, model))) = (untracked, &stack_model) {
        // The stack of an untracked session is a star-trail image, and its
        // own percentiles cannot set its levels: a star lands on any given
        // pixel in only a few frames, so the mean divides its flux by the
        // length of the session and the high percentile taken as white falls
        // back onto the sky. The window then collapses to about the width of
        // the noise — measured at 11% of the sky level on the test session —
        // and throws every residual of a per cent across a fifth of the
        // display range, which is what made a correctly flattened stack look
        // like it was full of defects. The levels come instead from one
        // representative frame, the same one and for the same reason the
        // export already uses.
        let rep = infos
            .iter()
            .filter(|f| f.in_stack)
            .min_by(|a, b| {
                let d = |f: &FrameInfo| (f.level - med_level).abs();
                d(a).partial_cmp(&d(b)).unwrap()
            })
            .ok_or_else(|| anyhow!("no frame available as the stack's levels reference"))?;
        say!(
            "stretch from {} (sky {:.4}, closest to the session median {:.4}) — a star-trail stack cannot set its own levels",
            paths[rep.idx].file_name().unwrap().to_string_lossy(), rep.level, med_level
        );
        let frame = raw::load(&paths[rep.idx])
            .with_context(|| format!("loading the stack's levels reference {}", paths[rep.idx].display()))?;
        let fc = flatten::fit_frame_corr_ex(&rep.bg, bg_med, model, Some(&rep.fit_mask), anomaly);
        let img = flatten::clean_frame(model, &frame, Some(&fc), true, args.scatter_comp);
        Some(output::analyze_stretch(&img))
    } else {
        say!("stretch analysis over the balanced stack ({}×{})", out_w, out_h);
        Some(output::analyze_stretch(&out))
    };
    ui::phase("writing the DNG");
    if args.no_stack {
        say!("no stack written (--no-stack); crop {}×{}", out_w, out_h);
    } else {
        output::write_dng(&args.output, &out, out_w, out_h, &camera, stretch)?;
        say!("DNG written: {}", args.output.display());
    }

    // Photographic metadata of the source RAW, written natively into every
    // DNG we produce: camera, lens, focal length, aperture, exposure, ISO
    // and capture time. That is what a developer matches a lens profile on,
    // and unlike `copy_metadata` below it needs no external tool. exiftool,
    // when installed, then adds MakerNotes, XMP and ICC on top.
    let src_exif = match exif::read_source(&paths[0]) {
        Ok(e) => {
            say!("EXIF from the source: {}", e.describe());
            Some(e)
        }
        Err(e) => {
            say!("WARNING: could not read the source EXIF: {e:#}");
            None
        }
    };
    if let (Some(e), false) = (&src_exif, args.no_stack) {
        exif::embed(&args.output, e)
            .with_context(|| format!("writing EXIF into {}", args.output.display()))?;
    }

    if let Some(path) = &args.dump_correction {
        match (&corr, &stack_model) {
            (Some(grid), Some((_, m))) => {
                let (x0, y0, _, _) = acc.valid_bounds();
                let layer = flatten::correction_layer(grid, m.pedestal, x0, y0, out_w, out_h);
                output::write_dng_quiet(path, &layer, out_w, out_h, &camera, stretch)?;
                if let Some(e) = &src_exif {
                    exif::embed(path, e)?;
                }
                say!("removed layer written: {}", path.display());
            }
            _ => say!("WARNING: --dump-correction ignored (correction disabled)"),
        }
    }

    if !args.no_stack {
        let src = paths[0].file_name().unwrap().to_string_lossy().to_string();
        match output::copy_metadata(&paths[0], &args.output) {
            Ok(()) => say!("MakerNotes + XMP + ICC copied from {src}"),
            Err(e) => say!("WARNING: could not copy MakerNotes/XMP/ICC from {src}: {e:#} — the DNG keeps the EXIF written natively"),
        }
    }

    if let Some(dir) = &args.export_clean {
        ui::phase("exporting the clean sequence");
        check_abort()?;
        let Some((bg_med, model)) = &stack_model else {
            return Err(anyhow!("--export-clean requires the correction to be enabled (drop --no-flatten)"));
        };
        let (cx0, cy0, _, _) = acc.valid_bounds();
        // Untracked: the stabilization is the tripod-drift correction, and
        // the temporal window works on the sky chain rather than on frames
        // that are already aligned, so it survives exporting in sensor
        // coordinates.
        let stabilize = !args.export_no_stabilize && !(untracked && args.fixed_no_stabilize);
        let window = if stabilize || untracked { args.export_window } else { 1 };
        // Levels reference. On a tracked sequence the stack is the right
        // one: it is the same sky as every frame, only deeper. On an
        // untracked one it is a star-trail image — a star lands on a given
        // pixel in only a few frames, so the mean divides its flux by the
        // length of the session and the stack's high percentile no longer
        // stands for a star core. Matching the exported frames to that
        // would crush their contrast frame after frame, so both the
        // deflicker reference and the stretch come instead from one
        // representative frame, cleaned exactly as the export cleans it.
        let (reference, export_stretch) = if untracked {
            let rep = infos
                .iter()
                .filter(|f| f.in_stack)
                .min_by(|a, b| {
                    let d = |f: &FrameInfo| (f.level - med_level).abs();
                    d(a).partial_cmp(&d(b)).unwrap()
                })
                .ok_or_else(|| anyhow!("no frame available as the export's levels reference"))?;
            say!(
                "levels reference: {} (sky {:.4}, closest to the session median {:.4}) — the star-trail stack is not one",
                paths[rep.idx].file_name().unwrap().to_string_lossy(), rep.level, med_level
            );
            let frame = raw::load(&paths[rep.idx])
                .with_context(|| format!("loading the levels reference {}", paths[rep.idx].display()))?;
            let fc = flatten::fit_frame_corr_ex(&rep.bg, bg_med, model, Some(&rep.fit_mask), anomaly);
            let img = flatten::clean_frame(model, &frame, Some(&fc), true, args.scatter_comp);
            drop(frame);
            let (img, pmask) = if stabilize {
                (
                    timelapse::warp_to_ref(&img, w_ref, h_ref, &rep.m, cx0, cy0, out_w, out_h),
                    timelapse::warp_mask_to_ref(&rep.mask, w_ref, h_ref, &rep.m, cx0, cy0, out_w, out_h),
                )
            } else {
                let mut pm = vec![0u8; w_ref * h_ref];
                for y in 0..h_ref {
                    for x in 0..w_ref {
                        pm[y * w_ref + x] = rep.mask.at_px(x as f32 + 0.5, y as f32 + 0.5) as u8;
                    }
                }
                (img, pm)
            };
            let (sw, sh) = if stabilize { (out_w, out_h) } else { (w_ref, h_ref) };
            let stats = timelapse::stats(&img, sw, sh, Some(&pmask));
            let st = if args.no_stretch { None } else { Some(output::analyze_stretch(&img)) };
            (stats, st)
        } else {
            (timelapse::stats(&out, out_w, out_h, None), stretch)
        };
        let opts = timelapse::ExportOpts {
            dir,
            window,
            stabilize,
            deflicker: !args.export_no_deflicker,
            keep_transients: !args.export_no_transients,
            camera: &camera,
            stretch: export_stretch,
            n_workers,
            scatter_comp: args.scatter_comp,
            sky_links: if sky_links.is_empty() { None } else { Some(&sky_links) },
            anomaly,
            cloud_guard: !args.export_no_cloud_guard,
        };
        drop(out);
        drop(acc);
        timelapse::export_sequence(&paths, &infos, model, bg_med, (cx0, cy0, out_w, out_h), reference, &opts)?;
    }
    Ok(())
}

/// Everything pass 1 knows about an aligned frame.
pub struct FrameInfo {
    /// Index into the path list (chronological order).
    pub idx: usize,
    /// Similarity ref → frame.
    pub m: Similarity,
    /// Blockwise background map (sensor).
    pub bg: flatten::BgMap,
    /// Per-cell foreground mask (sensor).
    pub mask: flatten::CellMask,
    /// The cells that take no part in any fit: the foreground, plus — on an
    /// untracked sequence — whatever cloud covered this frame. The
    /// foreground alone is what the export and the stack still use, because
    /// the landscape is part of the picture and a cloud is part of the
    /// frame; it is only the *model* that neither may enter.
    pub fit_mask: flatten::CellMask,
    /// Fraction of the frame this frame's cloud covers.
    pub cloud: f32,
    /// Sky level (median G of the map without foreground).
    pub level: f32,
    /// Whether it enters the stack.
    pub in_stack: bool,
    /// Whether it helps determine the background: the flat field, the
    /// session median map and the model fitted over it. A frame with cloud
    /// in it is still stacked and still exported — it is only barred from
    /// saying what the sky underneath looks like.
    pub in_model: bool,
    /// true = similarity fitted from stars; false = interpolated from the
    /// neighbours (exported only).
    pub aligned: bool,
    pub inliers: usize,
    /// Whether it is exported by --export-clean (false: anomalously bright
    /// sky — saturated dawn/twilight — where cleaning makes no sense).
    pub export: bool,
}

/// Cloud cover above which a frame is counted as clouded in the report.
const CLOUD_REPORT_FRAC: f32 = 0.05;
/// Cloud cover above which a frame no longer helps determine the
/// background.
///
/// Deliberately close to zero. A cell-by-cell mask is the right instrument
/// for the export, where every frame has to be cleaned whatever the weather
/// put in front of it, but not for the model: there the question is not "how
/// much of this frame is cloud" but "is this frame *certainly* clear", and
/// a detector answering the first question wrongly on 2 % of the cells is
/// answering the second one with a no. Frames enough are left — a night of
/// timelapse is hundreds — so doubt is cheap to act on and expensive to
/// ignore: what a wrongly kept frame does is put its cloud in the session
/// median, from where it is subtracted out of every frame of the night.
const CLOUD_FIT_FRAC: f32 = 0.02;
/// Fewest frames each half of the session needs for its own median map to
/// mean anything (see the residual surface's verification).
const MIN_HALF_FRAMES: usize = 8;
/// Fewest frames the background may be determined from before the threshold
/// above is given up on and the cleanest frames are taken instead. Two
/// dozen is what the flat field's per-cell regression needs to have a
/// median at all (`FLAT_MIN_SAMPLES` in either half).
const MODEL_MIN_FRAMES: usize = 24;
/// Cloud cover above which a frame no longer takes part in the model or the
/// stack. What is left of its sky is too little to measure a level on, and
/// the level is what every per-frame correction is scaled by. It is still
/// exported: cleaning it with the session's model is exactly what it needs.
const CLOUD_MAX_FRAC: f32 = 0.5;

/// Sky level (fraction of the sensor range, balanced scale) above which a
/// frame is considered saturated and is not exported.
const EXPORT_MAX_LEVEL: f32 = 0.85;

/// Minimum inliers to accept a star-based alignment.
const MIN_INLIERS: usize = 6;
/// Maximum deviation (px, at the corners) between the fitted similarity and
/// the one interpolated from the neighbours for it to be accepted.
const ALIGN_MAX_DEV_PX: f32 = 6.0;
/// Maximum distance (frames) to an aligned neighbour for interpolation.
const ALIGN_INTERP_RANGE: usize = 6;

/// Interpolates/extrapolates the similarity of frame `idx` from the nearest
/// star-aligned frames (before and after, up to `ALIGN_INTERP_RANGE`),
/// linearly in angle and translation. `skip` = position in `infos` to ignore
/// (the frame itself). None if there are no neighbours.
fn predict_transform(infos: &[FrameInfo], idx: usize, skip: Option<usize>) -> Option<Similarity> {
    let ok = |k: usize| infos[k].aligned && Some(k) != skip && infos[k].idx.abs_diff(idx) <= ALIGN_INTERP_RANGE;
    let prev = (0..infos.len()).filter(|&k| ok(k) && infos[k].idx < idx).max_by_key(|&k| infos[k].idx);
    let next = (0..infos.len()).filter(|&k| ok(k) && infos[k].idx > idx).min_by_key(|&k| infos[k].idx);
    let lerp = |a: &Similarity, b: &Similarity, t: f32| -> Similarity {
        let ta = a.sin_t.atan2(a.cos_t);
        let mut tb = b.sin_t.atan2(b.cos_t);
        // shortest path
        while tb - ta > std::f32::consts::PI { tb -= 2.0 * std::f32::consts::PI; }
        while ta - tb > std::f32::consts::PI { tb += 2.0 * std::f32::consts::PI; }
        let th = ta + (tb - ta) * t;
        Similarity {
            cos_t: th.cos(),
            sin_t: th.sin(),
            tx: a.tx + (b.tx - a.tx) * t,
            ty: a.ty + (b.ty - a.ty) * t,
        }
    };
    match (prev, next) {
        (Some(p), Some(n)) => {
            let t = (idx - infos[p].idx) as f32 / (infos[n].idx - infos[p].idx) as f32;
            Some(lerp(&infos[p].m, &infos[n].m, t))
        }
        (Some(p), None) => {
            // Extrapolate with the one before the previous, if it exists.
            let pp = (0..infos.len())
                .filter(|&k| ok(k) && infos[k].idx < infos[p].idx)
                .max_by_key(|&k| infos[k].idx);
            match pp {
                Some(pp) => {
                    let t = (idx - infos[pp].idx) as f32 / (infos[p].idx - infos[pp].idx) as f32;
                    Some(lerp(&infos[pp].m, &infos[p].m, t))
                }
                None => Some(infos[p].m),
            }
        }
        (None, Some(n)) => {
            let nn = (0..infos.len())
                .filter(|&k| ok(k) && infos[k].idx > infos[n].idx)
                .min_by_key(|&k| infos[k].idx);
            match nn {
                Some(nn) => {
                    let t = (idx as f32 - infos[n].idx as f32) / (infos[nn].idx - infos[n].idx) as f32;
                    Some(lerp(&infos[n].m, &infos[nn].m, t))
                }
                None => Some(infos[n].m),
            }
        }
        (None, None) => None,
    }
}

/// Maximum distance (px) between the images of the four sensor corners under
/// two similarities.
fn corner_deviation(a: &Similarity, b: &Similarity, w: usize, h: usize) -> f32 {
    let corners = [(0.0, 0.0), (w as f32, 0.0), (0.0, h as f32), (w as f32, h as f32)];
    corners
        .iter()
        .map(|&(x, y)| {
            let (ax, ay) = a.apply(x, y);
            let (bx, by) = b.apply(x, y);
            ((ax - bx).powi(2) + (ay - by).powi(2)).sqrt()
        })
        .fold(0.0, f32::max)
}

enum Msg {
    LoadErr { idx: usize, error: String },
    DimMismatch { idx: usize },
    Processed {
        idx: usize,
        n_stars: usize,
        fit: Option<(Similarity, usize)>,
        /// Untracked sequence: (dx, dy, correlation) of the tripod drift
        /// against the reference's landscape. None when it could not be
        /// measured, or on a tracked sequence.
        shift: Option<(f32, f32, f32)>,
        /// Untracked sequence: the frame's stars, kept so the chain of
        /// consecutive sky fits can be built after the pass. Empty when
        /// the chain is not needed.
        stars: Vec<stars::Star>,
        bg: Option<(flatten::BgMap, flatten::CellMask, f32)>,
        elapsed: f32,
    },
}

struct RamInfo {
    total_bytes: u64,
    budget_bytes: u64,
    frame_bytes: u64,
}

/// Reads /proc/meminfo (Linux) and returns MemTotal in bytes. On failure,
/// assumes a conservative value of 4 GiB.
fn detect_total_ram_bytes() -> u64 {
    const FALLBACK: u64 = 4 * (1 << 30);
    let Ok(meminfo) = fs::read_to_string("/proc/meminfo") else {
        return FALLBACK;
    };
    for line in meminfo.lines() {
        let Some(rest) = line.strip_prefix("MemTotal:") else {
            continue;
        };
        let kb: u64 = rest
            .split_whitespace()
            .next()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        if kb > 0 {
            return kb * 1024;
        }
    }
    FALLBACK
}

/// Computes the number of workers and the channel capacity, keeping RAM
/// usage below 80% of the total. Assumes:
/// - resident accumulator: w*h*20 B (sum f32×3 + count u32 + cover u32)
/// - frame in flight: w*h*24 B (the frame's rgb f32×3 plus its clean copy in
///   the stacking pass, or rgb + luminance in the alignment pass).
fn compute_pipeline_params(
    width: usize,
    height: usize,
    n_cpus: usize,
    remaining_frames: usize,
) -> (usize, usize, RamInfo) {
    let total = detect_total_ram_bytes();
    let budget = (total as f64 * 0.8) as u64;
    let px = (width as u64) * (height as u64);
    let acc_bytes = px * 20;
    let frame_bytes = px * 24;

    let available = budget.saturating_sub(acc_bytes);
    let max_in_flight = ((available / frame_bytes).max(1) as usize).min(remaining_frames.max(1));

    let n_workers = n_cpus.min(max_in_flight).max(1);
    let chan_cap = max_in_flight.saturating_sub(n_workers + 1);

    (
        n_workers,
        chan_cap,
        RamInfo {
            total_bytes: total,
            budget_bytes: budget,
            frame_bytes,
        },
    )
}

fn list_raw(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut v: Vec<PathBuf> = fs::read_dir(dir)
        .with_context(|| format!("listing {}", dir.display()))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| raw::is_raw_file(p))
        .collect();
    v.sort();
    Ok(v)
}
