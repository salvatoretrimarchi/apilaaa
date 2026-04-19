mod align;
mod flatten;
mod output;
mod raw;
mod stack;
mod stars;
mod timelapse;

use crate::align::Similarity;
use anyhow::{anyhow, Context, Result};
use clap::Parser;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::sync_channel;
use std::thread;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(about = "RAW stacker with translation and rotation correction")]
struct Args {
    /// Directory holding the input RAW files (.ARW)
    #[arg(short, long, default_value = "res")]
    input: PathBuf,

    /// Output DNG file
    #[arg(short, long, default_value = "stacked.dng")]
    output: PathBuf,

    /// Maximum number of stars used per frame for alignment
    #[arg(long, default_value_t = 40)]
    max_stars: usize,

    /// Process at most N frames (useful for testing)
    #[arg(long)]
    limit: Option<usize>,

    /// Disable the automatic histogram stretch; keeps the native sensor
    /// scale (useful if you want to do the whole development in darktable
    /// without any prior gain).
    #[arg(long)]
    no_stretch: bool,

    /// Disable the automatic removal of the lens halo/glare and of the
    /// background gradient (plane model + radial profile in sensor
    /// coordinates, see `flatten`). Enabled by default.
    #[arg(long)]
    no_flatten: bool,

    /// Disable the second correction stage (lower-envelope residual surface
    /// that removes non-radial bands/wedges); keeps only the parametric
    /// halo + gradient model.
    #[arg(long)]
    no_residual_surface: bool,

    /// Also write a DNG holding the **removed layer** (halo + gradient +
    /// surface + bands/radii, averaged like the stack, plus the sky
    /// pedestal) with the same crop and stretch as the output, so you can
    /// check in the developer that only defect was removed, not sky.
    #[arg(long, value_name = "DNG")]
    dump_correction: Option<PathBuf>,

    /// Directory to export the clean timelapse sequence to: every aligned
    /// frame (including the ones excluded from the stack) with the defect
    /// model + its own gradient and horizon glow subtracted, stabilized to
    /// the reference frame system (same alignment as the stack, same crop),
    /// with deflickering (per-channel median and high percentile matched to
    /// the stack's, foreground excluded) and temporal noise reduction
    /// (sliding window, trimmed mean; the foreground —trees, horizon— is
    /// kept as is, without smearing it with the neighbouring frames). Same
    /// stretch as the stack. Requires the correction to be enabled.
    #[arg(long, value_name = "DIR")]
    export_clean: Option<PathBuf>,

    /// Frames in the temporal window used to reduce noise in --export-clean
    /// (odd; 1 = no reduction). Noise ≈ /√(N−2).
    #[arg(long, default_value_t = 7, value_name = "N")]
    export_window: usize,

    /// Do not stabilize in --export-clean (exports in sensor coordinates,
    /// uncropped; temporal noise reduction is disabled too because it
    /// requires aligned frames).
    #[arg(long)]
    export_no_stabilize: bool,

    /// Do not apply deflickering in --export-clean.
    #[arg(long)]
    export_no_deflicker: bool,

    /// In --export-clean, do NOT preserve transients (meteors, satellites,
    /// planes) over the temporal combination. By default they are kept with
    /// the value from the original frame.
    #[arg(long)]
    export_no_transients: bool,

    /// Contrast compensation in the halo/veil areas (β): the veil light is
    /// scattered light missing from the structure of that area; after
    /// subtracting the veil, the deviation from the sky is rescaled by
    /// 1/(1 − β·veil/sky) to match the contrast at the centre. 0 = off.
    #[arg(long, default_value_t = 1.0, value_name = "BETA")]
    scatter_comp: f32,

    /// Frame selection for the stack: frames whose sky level (median of the
    /// background without foreground) falls outside [median/F, median·F] of
    /// the session (twilight, clouds, moon) are excluded. Excluded frames
    /// are still exported by --export-clean (unless the sky is saturated).
    #[arg(long, default_value_t = 1.6, value_name = "F")]
    stack_sky_tolerance: f32,

    /// Frame selection for the stack: frames with more than this fraction of
    /// the frame taken up by foreground (trees, horizon) are excluded.
    /// Defaults to 0: only frames with a clear sky are stacked (frames with
    /// foreground usually also carry the horizon glow, which varies over
    /// time and leaves bands when averaged). With a value > 0, frames below
    /// it are included and only those pixels are masked out. If fewer than
    /// 20 % of the aligned frames are left without foreground, frames are
    /// automatically admitted with a mask up to 60 %.
    #[arg(long, default_value_t = 0.0, value_name = "FRAC")]
    stack_max_foreground: f32,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let mut paths = list_arw(&args.input)?;
    if let Some(n) = args.limit {
        paths.truncate(n);
    }
    if paths.is_empty() {
        return Err(anyhow!("no ARW found in {}", args.input.display()));
    }
    println!("found {} frames", paths.len());

    let t0 = Instant::now();
    let ref_frame = raw::load(&paths[0])
        .with_context(|| format!("loading reference {}", paths[0].display()))?;
    println!(
        "reference: {} ({} x {}) in {:.2}s",
        paths[0].file_name().unwrap().to_string_lossy(),
        ref_frame.width,
        ref_frame.height,
        t0.elapsed().as_secs_f32()
    );

    let lum_ref = raw::luminance(&ref_frame);
    let stars_ref = stars::detect(&lum_ref, ref_frame.width, ref_frame.height, args.max_stars);
    drop(lum_ref);
    println!("stars in reference: {}", stars_ref.len());
    if stars_ref.len() < 6 {
        return Err(anyhow!("too few stars in the reference ({}) — check exposure or max_stars", stars_ref.len()));
    }

    let camera_model = ref_frame.camera.clone();
    let w_ref = ref_frame.width;
    let h_ref = ref_frame.height;
    let flatten_on = !args.no_flatten;

    // ---------------------------------------------------------------
    // Pass 1: alignment + background map + foreground mask for every
    // frame. Nothing is stacked yet: first we need to know which frames
    // are usable (sky level, foreground) and fit the defect model.
    // ---------------------------------------------------------------
    let mut infos: Vec<FrameInfo> = Vec::new();
    {
        let bg = flatten::block_background(&ref_frame);
        let mask = flatten::foreground_mask(&bg);
        let level = flatten::sky_level(&bg, &mask);
        infos.push(FrameInfo { idx: 0, m: Similarity::identity(), bg, mask, level, in_stack: true, aligned: true, inliers: usize::MAX, export: true });
    }
    drop(ref_frame);
    println!("  [1/{}] reference: foreground {:.1}%", paths.len(), 100.0 * infos[0].mask.fraction());

    let n_cpus = thread::available_parallelism().map(|n| n.get()).unwrap_or(2);
    let (n_workers, chan_cap, ram_info) = compute_pipeline_params(w_ref, h_ref, n_cpus, paths.len() - 1);
    println!(
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
            s.spawn(move || {
                loop {
                    let idx = next_idx.fetch_add(1, Ordering::Relaxed);
                    if idx >= total {
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
                                drop(lum);
                                let n_stars = stars_cur.len();
                                let fit = align::fit(stars_ref, &stars_cur);
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
            match msg {
                Msg::LoadErr { idx, error } => {
                    println!(
                        "  [{}/{}] {}: LOAD FAILED {error}",
                        idx + 1,
                        total,
                        paths[idx].file_name().unwrap().to_string_lossy()
                    );
                    skipped += 1;
                }
                Msg::DimMismatch { idx } => {
                    println!(
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
                    fit,
                    bg,
                    elapsed,
                } => {
                    match (fit, bg) {
                        (Some((m, inliers)), Some((bg, mask, level))) if inliers >= MIN_INLIERS => {
                            let fg = mask.fraction();
                            infos.push(FrameInfo { idx, m, bg, mask, level, in_stack: true, aligned: true, inliers, export: true });
                            aligned += 1;
                            println!(
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
                            println!(
                                "  [{}/{}] {}: NOT ALIGNED ({} stars detected, {} inliers); will try to interpolate the neighbours' alignment",
                                idx + 1,
                                total,
                                paths[idx].file_name().unwrap().to_string_lossy(),
                                n_stars,
                                inl
                            );
                            pending.push(FrameInfo {
                                idx, m: Similarity::identity(), bg, mask, level, in_stack: false, aligned: false, inliers: inl, export: true,
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
    {
        let mut n_fixed = 0usize;
        let mut bad: Vec<usize> = Vec::new();
        for k in 0..infos.len() {
            let Some(pred) = predict_transform(&infos, infos[k].idx, Some(k)) else { continue; };
            let dev = corner_deviation(&infos[k].m, &pred, w_ref, h_ref);
            if dev > ALIGN_MAX_DEV_PX {
                bad.push(k);
                println!(
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
                    println!(
                        "  {}: alignment interpolated from the neighbours (θ={:.3}°, t=({:.2},{:.2})); exported, out of the stack",
                        paths[f.idx].file_name().unwrap().to_string_lossy(), m.angle_deg(), m.tx, m.ty
                    );
                    f.m = m;
                    infos.push(f);
                    n_fixed += 1;
                }
                None => {
                    println!(
                        "  {}: no aligned neighbours nearby, skipped",
                        paths[f.idx].file_name().unwrap().to_string_lossy()
                    );
                    skipped += 1;
                }
            }
        }
        infos.sort_by_key(|f| f.idx);
        if n_fixed > 0 {
            println!("  {} frames with interpolated alignment", n_fixed);
        }
    }

    println!(
        "aligned: {}  interpolated: {}  skipped: {}  time: {:.1}s",
        aligned,
        infos.iter().filter(|f| !f.aligned).count(),
        skipped,
        t0.elapsed().as_secs_f32()
    );

    // ---------------------------------------------------------------
    // Frame selection for stacking: drop the ones with too much
    // foreground and the ones whose sky is anomalously bright or dark
    // compared to the session (twilight, clouds, moon rising).
    // ---------------------------------------------------------------
    {
        let mut lv: Vec<f32> = infos.iter().map(|f| f.level).filter(|v| *v > 0.0).collect();
        let med_level = if lv.is_empty() {
            0.0
        } else {
            let k = lv.len() / 2;
            let (_, m, _) = lv.select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap());
            *m
        };
        let tol = args.stack_sky_tolerance.max(1.0);
        let mut max_fg = args.stack_max_foreground.max(0.0);
        let n_clear = infos.iter().filter(|f| f.aligned && f.mask.fraction() <= max_fg).count();
        if n_clear < (infos.len() / 5).max(1) && max_fg < 0.6 {
            println!(
                "  only {} of {} frames have foreground ≤ {:.0}%: frames with foreground (masked) are admitted up to 60%",
                n_clear, infos.len(), 100.0 * max_fg
            );
            max_fg = 0.6;
        }
        let mut n_fg = 0usize;
        let mut n_bright = 0usize;
        let mut n_dark = 0usize;
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
            } else if med_level > 0.0 && f.level > med_level * tol {
                f.in_stack = false;
                n_bright += 1;
            } else if med_level > 0.0 && f.level < med_level / tol {
                f.in_stack = false;
                n_dark += 1;
            }
        }
        let n_stack = infos.iter().filter(|f| f.in_stack).count();
        println!(
            "stack selection: {} of {} aligned frames (median sky {:.4}, tolerance ×{:.2}); {} with foreground (masked); excluded: {} for foreground > {:.0}%, {} for bright sky, {} for dark sky; {} with saturated sky (≥ {:.0}%) are not exported",
            n_stack, infos.len(), med_level, tol, n_with_fg, n_fg, 100.0 * max_fg, n_bright, n_dark, n_sat, 100.0 * EXPORT_MAX_LEVEL
        );
        if n_stack == 0 {
            return Err(anyhow!("no frame passes the stack selection (tune --stack-sky-tolerance / --stack-max-foreground)"));
        }
        if std::env::var_os("APILAAA_DEBUG").is_some() {
            for f in infos.iter().filter(|f| !f.in_stack) {
                println!("  excluded from the stack: {} (sky {:.4}, foreground {:.1}%)",
                    paths[f.idx].file_name().unwrap().to_string_lossy(), f.level, 100.0 * f.mask.fraction());
            }
        }
    }

    // ---------------------------------------------------------------
    // Defect model: (masked) temporal median of the maps of the selected
    // frames.
    // ---------------------------------------------------------------
    let stack_infos: Vec<&FrameInfo> = infos.iter().filter(|f| f.in_stack).collect();
    let stack_transforms: Vec<Similarity> = stack_infos.iter().map(|f| f.m).collect();
    // Fits model + median map over a set of frames.
    let fit_model = |sel: &[&FrameInfo], label: &str| -> Result<(flatten::BgMap, flatten::GlareModel)> {
        let _ = label;
        let tf = Instant::now();
        let maps: Vec<flatten::BgMap> = sel.iter().map(|f| f.bg.clone()).collect();
        let masks: Vec<flatten::CellMask> = sel.iter().map(|f| f.mask.clone()).collect();
        let (bg, filled) = flatten::temporal_median_masked(&maps, &masks);
        if filled > 0 {
            println!("  [{label}] {} map cells always covered, filled in from neighbours", filled);
        }
        if let Some(dir) = std::env::var_os("APILAAA_DEBUG_DIR") {
            let idxs: Vec<usize> = sel.iter().map(|f| f.idx).collect();
            let tfs: Vec<Similarity> = sel.iter().map(|f| f.m).collect();
            let sub = Path::new(&dir).join(label);
            fs::create_dir_all(&sub)?;
            flatten::debug_dump_halves(&maps, &tfs, &idxs, &sub)?;
        }
        drop(maps);
        drop(masks);
        let model = flatten::GlareModel::fit(&bg, w_ref, h_ref, !args.no_residual_surface);
        println!("glare/gradient [{label}]: {}", model.report());
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
                println!("  radial profile {name} (every {:.0}px, % vs r=0): {}", model.r_step, p.join(" "));
            }
        }
        println!(
            "  model [{label}] fitted over {} frames in {:.2}s",
            sel.len(),
            tf.elapsed().as_secs_f32()
        );
        Ok((bg, model))
    };
    // A single model, fitted over the median of the selected frames (clear
    // sky, no twilight): it also serves as the reference for the temporal
    // anomaly of every frame in the sequence (`fit_frame_corr`), which
    // absorbs whatever that frame has more or less of compared to that
    // median (horizon glow, twilight, halo amplitude).
    let stack_model = if flatten_on { Some(fit_model(&stack_infos, "stack")?) } else { None };

    // ---------------------------------------------------------------
    // Pass 2: stacking. Every selected frame is loaded again, cleaned
    // (model + its own gradient) in sensor coordinates and accumulated
    // aligned, skipping the foreground pixels.
    // ---------------------------------------------------------------
    let mut acc = stack::Accumulator::new(w_ref, h_ref);
    {
        let ts = Instant::now();
        let n_stack = stack_infos.len();
        println!("stacking {} frames...", n_stack);
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
                        if k >= n_stack {
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
                                        let fc = flatten::fit_frame_corr(&info.bg, bg_med, m, Some(&info.mask));
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
                acc.add(&img, w_ref, h_ref, [1.0; 3], &info.m, Some(&info.mask));
                done += 1;
                if done % 25 == 0 || done == n_stack {
                    println!("  [{}/{}] stacked ({:.1}s)", done, n_stack, ts.elapsed().as_secs_f32());
                }
            }
            Ok::<(), anyhow::Error>(())
        })?;
        println!("stacked: {}  total time: {:.1}s", done, t0.elapsed().as_secs_f32());
    }

    let corr = stack_model
        .as_ref()
        .map(|(_, m)| flatten::CorrGrid::build(m, &stack_transforms, acc.width, acc.height));
    let (out, out_w, out_h) = acc.finalize_cropped(None);
    if out_w != acc.width || out_h != acc.height {
        println!(
            "crop from drift: {}×{} → {}×{} (−{} px height, −{} px width)",
            acc.width, acc.height, out_w, out_h,
            acc.height - out_h, acc.width - out_w
        );
    }
    {
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
            println!(
                "  stack coverage: minimum {:.0} samples/pixel out of {} frames{}",
                min_count,
                stack_infos.len(),
                if holes > 0 { format!("; {} px with no sample at all (always covered)", holes) } else { String::new() }
            );
        }
    }
    let stretch = if args.no_stretch {
        None
    } else {
        println!("stretch analysis over the balanced stack ({}×{})", out_w, out_h);
        Some(output::analyze_stretch(&out))
    };
    output::write_dng(&args.output, &out, out_w, out_h, &camera_model, stretch)?;
    println!("DNG written: {}", args.output.display());

    if let Some(path) = &args.dump_correction {
        match (&corr, &stack_model) {
            (Some(grid), Some((_, m))) => {
                let (x0, y0, _, _) = acc.valid_bounds();
                let layer = flatten::correction_layer(grid, m.pedestal, x0, y0, out_w, out_h);
                output::write_dng_quiet(path, &layer, out_w, out_h, &camera_model, stretch)?;
                println!("removed layer written: {}", path.display());
            }
            _ => println!("WARNING: --dump-correction ignored (correction disabled)"),
        }
    }

    print!("copying EXIF + MakerNotes from {}... ", paths[0].file_name().unwrap().to_string_lossy());
    use std::io::Write;
    std::io::stdout().flush().ok();
    match output::copy_metadata(&paths[0], &args.output) {
        Ok(()) => println!("ok"),
        Err(e) => println!("WARNING: {e:#} — valid DNG but without EXIF"),
    }

    if let Some(dir) = &args.export_clean {
        let Some((bg_med, model)) = &stack_model else {
            return Err(anyhow!("--export-clean requires the correction to be enabled (drop --no-flatten)"));
        };
        let (cx0, cy0, _, _) = acc.valid_bounds();
        let stabilize = !args.export_no_stabilize;
        let opts = timelapse::ExportOpts {
            dir,
            window: if stabilize { args.export_window } else { 1 },
            stabilize,
            deflicker: !args.export_no_deflicker,
            keep_transients: !args.export_no_transients,
            camera_model: &camera_model,
            stretch,
            n_workers,
            scatter_comp: args.scatter_comp,
        };
        let reference = timelapse::stats(&out, None);
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
    /// Sky level (median G of the map without foreground).
    pub level: f32,
    /// Whether it enters the stack and the model fit.
    pub in_stack: bool,
    /// true = similarity fitted from stars; false = interpolated from the
    /// neighbours (exported only).
    pub aligned: bool,
    pub inliers: usize,
    /// Whether it is exported by --export-clean (false: anomalously bright
    /// sky — saturated dawn/twilight — where cleaning makes no sense).
    pub export: bool,
}

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

fn list_arw(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut v: Vec<PathBuf> = fs::read_dir(dir)
        .with_context(|| format!("listing {}", dir.display()))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("ARW"))
                .unwrap_or(false)
        })
        .collect();
    v.sort();
    Ok(v)
}
