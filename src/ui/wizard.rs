//! The setup screen shown when `apilaaa` is started with no arguments.
//!
//! It is a form over exactly the options the command line takes — nothing
//! is reachable here that is not reachable there, and the equivalent
//! command line is written out at the bottom as it is edited, so the run
//! can be repeated or scripted without going through this screen again.
//! Options that do not apply are not drawn: the tripod section only exists
//! for an untracked sequence, the timelapse settings only once the export
//! is on.

use crate::{AnomalyArg, Args};
use anyhow::Result;
use clap::Parser;
use ratatui::crossterm::event::{self, Event, KeyCode, KeyEventKind, KeyModifiers};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Paragraph};
use ratatui::Frame;
use std::path::{Path, PathBuf};

/// Runs the screen. `Ok(None)` means the user left without starting a run.
pub fn run() -> Result<Option<Args>> {
    let mut form = Form::new();
    let mut term = ratatui::init();
    let outcome = loop {
        if let Err(e) = term.draw(|f| form.draw(f)) {
            ratatui::restore();
            return Err(e.into());
        }
        match event::read() {
            Err(e) => {
                ratatui::restore();
                return Err(e.into());
            }
            Ok(Event::Key(k)) if k.kind != KeyEventKind::Release => match form.key(k.code, k.modifiers) {
                Step::Continue => {}
                Step::Quit => break None,
                Step::Start => break Some(form.cfg.to_args()),
            },
            Ok(_) => {}
        }
    };
    ratatui::restore();
    Ok(outcome)
}

enum Step {
    Continue,
    Quit,
    Start,
}

// ---------------------------------------------------------------------
// The settings themselves
// ---------------------------------------------------------------------

/// The form's own shape of the options: booleans stated positively (the
/// command line spells most of them as `--no-…`, which is right for a flag
/// and wrong for a checkbox), and the two "off or a path" options split
/// into a switch and a path.
struct Cfg {
    input: String,
    untracked: bool,
    limit: usize,
    stack: bool,
    output: String,
    stretch: bool,
    flatten: bool,
    residual_surface: bool,
    scatter_comp: f32,
    dump_correction: String,
    export: bool,
    export_dir: String,
    export_window: usize,
    export_stabilize: bool,
    export_deflicker: bool,
    export_transients: bool,
    export_cloud_guard: bool,
    fixed_stabilize: bool,
    fixed_search: usize,
    fixed_anomaly: AnomalyArg,
    sky_tolerance: f32,
    /// Negative means "leave it to the default for the sequence type",
    /// which is what omitting `--stack-max-foreground` does.
    max_foreground: f32,
    max_stars: usize,
}

impl Cfg {
    /// Starts from clap's own defaults rather than a second copy of them,
    /// so the screen and `--help` can never drift apart.
    fn new() -> Self {
        let d = Args::parse_from(["apilaaa"]);
        Cfg {
            input: d.input.display().to_string(),
            untracked: d.fixed_tripod,
            limit: d.limit.unwrap_or(0),
            stack: !d.no_stack,
            output: d.output.display().to_string(),
            stretch: !d.no_stretch,
            flatten: !d.no_flatten,
            residual_surface: !d.no_residual_surface,
            scatter_comp: d.scatter_comp,
            dump_correction: String::new(),
            export: false,
            export_dir: String::new(),
            export_window: d.export_window,
            export_stabilize: !d.export_no_stabilize,
            export_deflicker: !d.export_no_deflicker,
            export_transients: !d.export_no_transients,
            export_cloud_guard: !d.export_no_cloud_guard,
            fixed_stabilize: !d.fixed_no_stabilize,
            fixed_search: d.fixed_search,
            fixed_anomaly: d.fixed_anomaly,
            sky_tolerance: d.stack_sky_tolerance,
            max_foreground: -1.0,
            max_stars: d.max_stars,
        }
    }

    fn to_args(&self) -> Args {
        let mut a = Args::parse_from(["apilaaa"]);
        a.input = PathBuf::from(self.input.trim());
        a.fixed_tripod = self.untracked;
        a.limit = (self.limit > 0).then_some(self.limit);
        a.no_stack = !self.stack;
        a.output = PathBuf::from(self.output.trim());
        a.no_stretch = !self.stretch;
        a.no_flatten = !self.flatten;
        a.no_residual_surface = !self.residual_surface;
        a.scatter_comp = self.scatter_comp;
        a.dump_correction = (!self.dump_correction.trim().is_empty())
            .then(|| PathBuf::from(self.dump_correction.trim()));
        a.export_clean = self.export.then(|| PathBuf::from(self.export_dir.trim()));
        a.export_window = self.export_window;
        a.export_no_stabilize = !self.export_stabilize;
        a.export_no_deflicker = !self.export_deflicker;
        a.export_no_transients = !self.export_transients;
        a.export_no_cloud_guard = !self.export_cloud_guard;
        a.fixed_no_stabilize = !self.fixed_stabilize;
        a.fixed_search = self.fixed_search;
        a.fixed_anomaly = self.fixed_anomaly;
        a.stack_sky_tolerance = self.sky_tolerance;
        a.stack_max_foreground = (self.max_foreground >= 0.0).then_some(self.max_foreground);
        a.max_stars = self.max_stars;
        a
    }

    /// The command line that would do the same thing. Only what departs
    /// from the defaults is written out, which is also what makes it worth
    /// copying.
    fn command_line(&self) -> String {
        let d = Args::parse_from(["apilaaa"]);
        let mut s = String::from("apilaaa");
        let q = |v: &str| {
            if v.contains(' ') { format!("'{v}'") } else { v.to_string() }
        };
        if self.input.trim() != d.input.display().to_string() {
            s += &format!(" -i {}", q(self.input.trim()));
        }
        if self.stack && self.output.trim() != d.output.display().to_string() {
            s += &format!(" -o {}", q(self.output.trim()));
        }
        if self.untracked {
            s += " --fixed-tripod";
        }
        if self.limit > 0 {
            s += &format!(" --limit {}", self.limit);
        }
        if !self.stack {
            s += " --no-stack";
        }
        if !self.stretch {
            s += " --no-stretch";
        }
        if !self.flatten {
            s += " --no-flatten";
        }
        if self.flatten && !self.residual_surface {
            s += " --no-residual-surface";
        }
        if self.flatten && (self.scatter_comp - d.scatter_comp).abs() > 1e-6 {
            s += &format!(" --scatter-comp {:.2}", self.scatter_comp);
        }
        if !self.dump_correction.trim().is_empty() {
            s += &format!(" --dump-correction {}", q(self.dump_correction.trim()));
        }
        if self.export {
            s += &format!(" --export-clean {}", q(self.export_dir.trim()));
            if self.export_window != d.export_window {
                s += &format!(" --export-window {}", self.export_window);
            }
            if !self.export_stabilize {
                s += " --export-no-stabilize";
            }
            if !self.export_deflicker {
                s += " --export-no-deflicker";
            }
            if !self.export_transients {
                s += " --export-no-transients";
            }
            if !self.export_cloud_guard {
                s += " --export-no-cloud-guard";
            }
        }
        if self.untracked {
            if !self.fixed_stabilize {
                s += " --fixed-no-stabilize";
            }
            if self.fixed_search != d.fixed_search {
                s += &format!(" --fixed-search {}", self.fixed_search);
            }
            if self.fixed_anomaly != d.fixed_anomaly {
                s += &format!(" --fixed-anomaly {}", anomaly_name(self.fixed_anomaly));
            }
        }
        if (self.sky_tolerance - d.stack_sky_tolerance).abs() > 1e-6 {
            s += &format!(" --stack-sky-tolerance {:.2}", self.sky_tolerance);
        }
        if self.max_foreground >= 0.0 {
            s += &format!(" --stack-max-foreground {:.2}", self.max_foreground);
        }
        if self.max_stars != d.max_stars {
            s += &format!(" --max-stars {}", self.max_stars);
        }
        s
    }
}

fn anomaly_name(a: AnomalyArg) -> &'static str {
    match a {
        AnomalyArg::Coarse => "coarse",
        AnomalyArg::None => "none",
        AnomalyArg::Full => "full",
    }
}

// ---------------------------------------------------------------------
// Fields
// ---------------------------------------------------------------------

#[derive(Copy, Clone, PartialEq, Eq)]
enum Id {
    Input,
    Sequence,
    Limit,
    Stack,
    Output,
    Stretch,
    FixedStabilize,
    FixedSearch,
    FixedAnomaly,
    Flatten,
    ResidualSurface,
    ScatterComp,
    DumpCorrection,
    Export,
    ExportDir,
    ExportWindow,
    ExportStabilize,
    ExportDeflicker,
    ExportTransients,
    ExportCloudGuard,
    SkyTolerance,
    MaxForeground,
    MaxStars,
}

enum Row {
    Section(&'static str),
    Field(Id),
}

/// How a field is edited. `Text` and `Number` take typing; the rest only
/// answer to ← and →.
enum Kind {
    Text,
    Number,
    Switch,
    Choice,
}

fn kind(id: Id) -> Kind {
    use Id::*;
    match id {
        Input | Output | DumpCorrection | ExportDir => Kind::Text,
        Limit | FixedSearch | ScatterComp | ExportWindow | SkyTolerance | MaxForeground
        | MaxStars => Kind::Number,
        FixedAnomaly => Kind::Choice,
        _ => Kind::Switch,
    }
}

fn label(id: Id) -> &'static str {
    use Id::*;
    match id {
        Input => "frames directory",
        Sequence => "sequence",
        Limit => "frame limit",
        Stack => "write the stack",
        Output => "output DNG",
        Stretch => "auto stretch",
        FixedStabilize => "measure tripod drift",
        FixedSearch => "drift search",
        FixedAnomaly => "per-frame anomaly",
        Flatten => "remove halo + gradient",
        ResidualSurface => "residual surface",
        ScatterComp => "scatter compensation",
        DumpCorrection => "dump removed layer",
        Export => "export timelapse",
        ExportDir => "export directory",
        ExportWindow => "temporal window",
        ExportStabilize => "stabilize",
        ExportDeflicker => "deflicker",
        ExportTransients => "keep transients",
        ExportCloudGuard => "cloud guard",
        SkyTolerance => "sky tolerance",
        MaxForeground => "max foreground",
        MaxStars => "stars per frame",
    }
}

fn help(id: Id) -> &'static str {
    use Id::*;
    match id {
        Input => "Directory holding the .ARW files. They are taken in name order, which for a camera is capture order.",
        Sequence => "Tracked: an equatorial mount held the stars still, so the frames are aligned on the stars. Untracked: the camera sat on a fixed tripod, the landscape is what stays still and the stack comes out as star trails.",
        Limit => "Process only the first N frames. 0 uses the whole directory. A few frames are enough to see whether the settings are right before committing to the session.",
        Stack => "Off skips loading, cleaning and averaging entirely: only the crop is worked out, for the timelapse export. Needs the export on, and an untracked sequence.",
        Output => "The stacked linear DNG. Camera, lens, exposure and capture time are written into it from the first source frame.",
        Stretch => "Off keeps the native sensor scale, for doing the whole development in darktable with no prior gain.",
        FixedStabilize => "Measures the residual tripod drift against the reference's landscape and undoes it. Off assumes the camera never moved and exports in sensor coordinates, uncropped.",
        FixedSearch => "Largest tripod drift searched for, in sensor pixels.",
        FixedAnomaly => "coarse: the per-frame anomaly may only follow structure far broader than the Milky Way — horizon glow, twilight, a rising moon. none: the defect model alone. full: the tracked-sequence surface, which on an untracked one takes part of the drifting Milky Way with it.",
        Flatten => "Fits and subtracts the lens halo/glare and the background gradient. Off leaves both in, and the timelapse export cannot run.",
        ResidualSurface => "The second stage: a lower-envelope surface that removes the bands and wedges the parametric halo + gradient model leaves behind.",
        ScatterComp => "Contrast compensation in the veiled areas: after the veil is subtracted, the deviation from the sky is rescaled to match the contrast at the centre. 0 turns it off.",
        DumpCorrection => "Also write a DNG of the layer that was removed, cropped and stretched like the output, so a developer can confirm it took defect and not sky. Empty = off.",
        Export => "Export every frame cleaned, stabilized, deflickered and temporally denoised — the timelapse. Excluded frames are exported too.",
        ExportDir => "Where the cleaned sequence is written. One DNG per frame, so count on around 145 MB each.",
        ExportWindow => "Frames combined per exported frame, odd. Noise falls roughly as the square root of N-2. 1 turns the temporal reduction off.",
        ExportStabilize => "Off exports in sensor coordinates, uncropped, and turns off the temporal reduction with it: that one needs aligned frames.",
        ExportDeflicker => "Matches each frame's per-channel median and high percentile to the reference's, ignoring the foreground, so the sequence stops pulsing.",
        ExportTransients => "Keeps meteors, satellites and planes at their original frame's value instead of letting the temporal combination average them away.",
        ExportCloudGuard => "Treats cloud like foreground: out of the temporal window, since cloud is the one thing that does not follow the sky, and out of the frame's level statistics, since a light-polluted cloud outshines every star.",
        SkyTolerance => "Frames whose sky level falls outside [median/F, median x F] of the session — twilight, cloud, moonrise — stay out of the stack. They are still exported.",
        MaxForeground => "Frames with more of the frame than this taken up by trees or horizon stay out of the stack. 'auto' is 0 for a tracked sequence and 1 for an untracked one, where the landscape is part of the picture.",
        MaxStars => "Stars used per frame for the alignment fit.",
    }
}

/// The visible rows. What does not apply is not drawn at all, which is the
/// difference between a form and a wall of switches.
fn rows(cfg: &Cfg) -> Vec<Row> {
    use Id::*;
    let mut v = vec![Row::Section("frames"), Row::Field(Input), Row::Field(Sequence), Row::Field(Limit)];
    v.push(Row::Section("stack"));
    v.push(Row::Field(Stack));
    if cfg.stack {
        v.push(Row::Field(Output));
        v.push(Row::Field(Stretch));
    }
    if cfg.untracked {
        v.push(Row::Section("tripod"));
        v.push(Row::Field(FixedStabilize));
        if cfg.fixed_stabilize {
            v.push(Row::Field(FixedSearch));
        }
        v.push(Row::Field(FixedAnomaly));
    }
    v.push(Row::Section("correction"));
    v.push(Row::Field(Flatten));
    if cfg.flatten {
        v.push(Row::Field(ResidualSurface));
        v.push(Row::Field(ScatterComp));
        if cfg.stack {
            v.push(Row::Field(DumpCorrection));
        }
    }
    v.push(Row::Section("timelapse"));
    v.push(Row::Field(Export));
    if cfg.export {
        v.push(Row::Field(ExportDir));
        v.push(Row::Field(ExportStabilize));
        if cfg.export_stabilize || cfg.untracked {
            v.push(Row::Field(ExportWindow));
        }
        v.push(Row::Field(ExportDeflicker));
        v.push(Row::Field(ExportTransients));
        v.push(Row::Field(ExportCloudGuard));
    }
    v.push(Row::Section("frame selection"));
    v.push(Row::Field(SkyTolerance));
    v.push(Row::Field(MaxForeground));
    if !cfg.untracked {
        v.push(Row::Field(MaxStars));
    }
    v
}

// ---------------------------------------------------------------------
// The screen
// ---------------------------------------------------------------------

struct Form {
    cfg: Cfg,
    /// Index into `rows(&cfg)`; always lands on a field, never a heading.
    cursor: usize,
    /// First visible row, moved only as much as it takes to keep the
    /// cursor on screen.
    top: usize,
    /// The focused field's text while it is being typed into. Committed on
    /// Enter, on leaving the field, and before starting the run.
    edit: Option<String>,
    /// Cached count of .ARW files in `cfg.input`, and what it was counted
    /// for.
    scan: (String, std::result::Result<usize, String>),
    /// Why the run cannot start yet, shown under the command line.
    error: Option<String>,
}

impl Form {
    fn new() -> Self {
        let cfg = Cfg::new();
        let scan = (cfg.input.clone(), count_arw(Path::new(&cfg.input)));
        let mut f = Form { cfg, cursor: 0, top: 0, edit: None, scan, error: None };
        f.cursor = f.next_field(0, 1).unwrap_or(0);
        f
    }

    /// First row at or after `from` in direction `dir` that is a field.
    fn next_field(&self, from: usize, dir: isize) -> Option<usize> {
        let rows = rows(&self.cfg);
        let mut i = from as isize;
        while i >= 0 && (i as usize) < rows.len() {
            if matches!(rows[i as usize], Row::Field(_)) {
                return Some(i as usize);
            }
            i += dir;
        }
        None
    }

    fn current(&self) -> Option<Id> {
        match rows(&self.cfg).get(self.cursor) {
            Some(Row::Field(id)) => Some(*id),
            _ => None,
        }
    }

    fn key(&mut self, code: KeyCode, mods: KeyModifiers) -> Step {
        if mods.contains(KeyModifiers::CONTROL) {
            return match code {
                KeyCode::Char('c') => Step::Quit,
                KeyCode::Char('r') => self.start(),
                _ => Step::Continue,
            };
        }
        match code {
            KeyCode::Esc => {
                if self.edit.is_some() {
                    self.edit = None;
                    Step::Continue
                } else {
                    Step::Quit
                }
            }
            KeyCode::F(5) => self.start(),
            KeyCode::Enter => {
                if self.edit.is_some() {
                    self.commit();
                    Step::Continue
                } else {
                    self.start()
                }
            }
            KeyCode::Up => {
                self.commit();
                if self.cursor > 0 {
                    if let Some(i) = self.next_field(self.cursor - 1, -1) {
                        self.cursor = i;
                    }
                }
                Step::Continue
            }
            KeyCode::Down | KeyCode::Tab => {
                self.commit();
                if let Some(i) = self.next_field(self.cursor + 1, 1) {
                    self.cursor = i;
                }
                Step::Continue
            }
            KeyCode::Left => {
                self.commit();
                self.bump(-1);
                Step::Continue
            }
            KeyCode::Right => {
                self.commit();
                self.bump(1);
                Step::Continue
            }
            KeyCode::Char(' ') if !matches!(self.current().map(kind), Some(Kind::Text)) => {
                self.commit();
                self.bump(1);
                Step::Continue
            }
            KeyCode::Backspace => {
                if let Some(id) = self.current() {
                    if matches!(kind(id), Kind::Text | Kind::Number) {
                        let buf = self.edit.get_or_insert_with(|| value(&self.cfg, id, &self.scan).0);
                        buf.pop();
                    }
                }
                Step::Continue
            }
            KeyCode::Char(c) => {
                if let Some(id) = self.current() {
                    match kind(id) {
                        Kind::Text => {
                            let buf = self.edit.get_or_insert_with(String::new);
                            buf.push(c);
                        }
                        Kind::Number if c.is_ascii_digit() || c == '.' => {
                            // Typing over a number replaces it: nobody
                            // means to append a digit to a default.
                            let buf = self.edit.get_or_insert_with(String::new);
                            buf.push(c);
                        }
                        _ => {}
                    }
                }
                Step::Continue
            }
            _ => Step::Continue,
        }
    }

    /// Writes the edit buffer back into the settings, if there is one.
    /// A value that does not parse is dropped and the old one stands.
    fn commit(&mut self) {
        let Some(buf) = self.edit.take() else { return };
        let Some(id) = self.current() else { return };
        let n = |s: &str, cur: usize| s.trim().parse::<usize>().unwrap_or(cur);
        let f = |s: &str, cur: f32| s.trim().parse::<f32>().unwrap_or(cur);
        match id {
            Id::Input => {
                self.cfg.input = buf.trim().to_string();
                self.rescan();
            }
            Id::Output => self.cfg.output = buf.trim().to_string(),
            Id::ExportDir => self.cfg.export_dir = buf.trim().to_string(),
            Id::DumpCorrection => self.cfg.dump_correction = buf.trim().to_string(),
            Id::Limit => self.cfg.limit = n(&buf, self.cfg.limit),
            Id::MaxStars => self.cfg.max_stars = n(&buf, self.cfg.max_stars).clamp(5, 2000),
            Id::FixedSearch => self.cfg.fixed_search = n(&buf, self.cfg.fixed_search).clamp(4, 1024),
            Id::ExportWindow => {
                // Even windows have no centre frame; the pipeline forces
                // them odd anyway, so say so here rather than there.
                self.cfg.export_window = (n(&buf, self.cfg.export_window).clamp(1, 199)) | 1;
            }
            Id::ScatterComp => self.cfg.scatter_comp = f(&buf, self.cfg.scatter_comp).clamp(0.0, 4.0),
            Id::SkyTolerance => self.cfg.sky_tolerance = f(&buf, self.cfg.sky_tolerance).clamp(1.0, 20.0),
            Id::MaxForeground => self.cfg.max_foreground = f(&buf, self.cfg.max_foreground).clamp(0.0, 1.0),
            _ => {}
        }
    }

    fn rescan(&mut self) {
        if self.scan.0 != self.cfg.input {
            self.scan = (self.cfg.input.clone(), count_arw(Path::new(&self.cfg.input)));
        }
    }

    /// ← / → on the focused field.
    fn bump(&mut self, dir: i32) {
        let Some(id) = self.current() else { return };
        let up = dir > 0;
        let c = &mut self.cfg;
        let step_u = |v: &mut usize, s: usize, lo: usize, hi: usize| {
            *v = if up { v.saturating_add(s) } else { v.saturating_sub(s) }.clamp(lo, hi);
        };
        let step_f = |v: &mut f32, s: f32, lo: f32, hi: f32| {
            *v = (*v + if up { s } else { -s }).clamp(lo, hi);
            // Keep the two decimals the field is displayed with, so
            // repeated steps do not accumulate a float tail.
            *v = (*v * 100.0).round() / 100.0;
        };
        match id {
            Id::Sequence => c.untracked = !c.untracked,
            Id::Stack => c.stack = !c.stack,
            Id::Stretch => c.stretch = !c.stretch,
            Id::Flatten => c.flatten = !c.flatten,
            Id::ResidualSurface => c.residual_surface = !c.residual_surface,
            Id::Export => {
                c.export = !c.export;
                // A directory is needed the moment the export is on, and
                // one named after the input is what the README's examples
                // use anyway.
                if c.export && c.export_dir.trim().is_empty() {
                    c.export_dir = format!("{}_clean", c.input.trim().trim_end_matches('/'));
                }
            }
            Id::ExportStabilize => c.export_stabilize = !c.export_stabilize,
            Id::ExportDeflicker => c.export_deflicker = !c.export_deflicker,
            Id::ExportTransients => c.export_transients = !c.export_transients,
            Id::ExportCloudGuard => c.export_cloud_guard = !c.export_cloud_guard,
            Id::FixedStabilize => c.fixed_stabilize = !c.fixed_stabilize,
            Id::FixedAnomaly => {
                c.fixed_anomaly = match (c.fixed_anomaly, up) {
                    (AnomalyArg::Coarse, true) => AnomalyArg::None,
                    (AnomalyArg::None, true) => AnomalyArg::Full,
                    (AnomalyArg::Full, true) => AnomalyArg::Coarse,
                    (AnomalyArg::Coarse, false) => AnomalyArg::Full,
                    (AnomalyArg::None, false) => AnomalyArg::Coarse,
                    (AnomalyArg::Full, false) => AnomalyArg::None,
                };
            }
            Id::Limit => step_u(&mut c.limit, 5, 0, 100_000),
            Id::MaxStars => step_u(&mut c.max_stars, 5, 5, 2000),
            Id::FixedSearch => step_u(&mut c.fixed_search, 8, 4, 1024),
            Id::ExportWindow => {
                let v = if up { c.export_window + 2 } else { c.export_window.saturating_sub(2) };
                c.export_window = v.clamp(1, 199) | 1;
            }
            Id::ScatterComp => step_f(&mut c.scatter_comp, 0.1, 0.0, 4.0),
            Id::SkyTolerance => step_f(&mut c.sky_tolerance, 0.1, 1.0, 20.0),
            Id::MaxForeground => {
                // One step below zero is "auto", which is a real setting
                // and not the same as 0.
                if !up && c.max_foreground <= 0.0 {
                    c.max_foreground = -1.0;
                } else if up && c.max_foreground < 0.0 {
                    c.max_foreground = 0.0;
                } else {
                    step_f(&mut c.max_foreground, 0.05, 0.0, 1.0);
                }
            }
            Id::Input | Id::Output | Id::ExportDir | Id::DumpCorrection => {}
        }
        // Hiding or revealing rows must not leave the cursor on a heading.
        if !matches!(rows(&self.cfg).get(self.cursor), Some(Row::Field(_))) {
            self.cursor = self.next_field(self.cursor, -1).unwrap_or(0);
        }
    }

    /// Checks the combination the pipeline would refuse anyway, so it is
    /// said here where it can still be corrected.
    fn validate(&self) -> Option<String> {
        let c = &self.cfg;
        if c.input.trim().is_empty() {
            return Some("no frames directory given".into());
        }
        match &self.scan.1 {
            Err(e) => return Some(e.clone()),
            Ok(0) => return Some(format!("no .ARW found in {}", c.input.trim())),
            Ok(_) => {}
        }
        if !c.stack && !c.export {
            return Some("with the stack off and no export the run would produce nothing".into());
        }
        if !c.stack && !c.untracked {
            return Some("the stack can only be skipped on an untracked sequence: on a tracked one it is what every exported frame is levelled against".into());
        }
        if c.stack && c.output.trim().is_empty() {
            return Some("no output DNG given".into());
        }
        if c.export && c.export_dir.trim().is_empty() {
            return Some("no export directory given".into());
        }
        if c.export && !c.flatten {
            return Some("the timelapse export needs the correction on".into());
        }
        None
    }

    fn start(&mut self) -> Step {
        self.commit();
        self.rescan();
        match self.validate() {
            Some(e) => {
                self.error = Some(e);
                Step::Continue
            }
            None => Step::Start,
        }
    }

    // -----------------------------------------------------------------
    // Drawing
    // -----------------------------------------------------------------

    fn draw(&mut self, f: &mut Frame) {
        let area = f.area();
        // On a short terminal the help text is what gives way: the command
        // line under it is the one panel that has to stay legible, since it
        // is the thing worth copying out of this screen.
        let help_h = if area.height < 22 { 3 } else { 5 };
        let [body, help_area, cmd, foot] = Layout::vertical([
            Constraint::Min(5),
            Constraint::Length(help_h),
            Constraint::Length(if self.error.is_some() { 4 } else { 3 }),
            Constraint::Length(1),
        ])
        .areas(area);
        self.draw_body(f, body);
        self.draw_help(f, help_area);
        self.draw_cmd(f, cmd);
        self.draw_footer(f, foot);
    }

    fn draw_body(&mut self, f: &mut Frame, area: Rect) {
        let frames = match &self.scan.1 {
            Ok(0) => Span::styled("  no .ARW here", Style::default().fg(Color::Red)),
            Ok(n) => Span::styled(format!("  {n} frames"), Style::default().fg(Color::Green)),
            Err(_) => Span::styled("  not a readable directory", Style::default().fg(Color::Red)),
        };
        let block = Block::bordered()
            .title(Line::from(vec![
                Span::styled(" apilaaa ", Style::default().add_modifier(Modifier::BOLD)),
                Span::styled(concat!("v", env!("CARGO_PKG_VERSION"), " — setup "), dim()),
            ]))
            .border_style(dim());
        let inner = block.inner(area);
        f.render_widget(block, area);

        let rows = rows(&self.cfg);
        let h = inner.height.max(1) as usize;
        // Scroll only as far as it takes: the cursor stays put on screen
        // while the list moves under it.
        if self.cursor < self.top {
            self.top = self.cursor;
        }
        if self.cursor >= self.top + h {
            self.top = self.cursor + 1 - h;
        }
        self.top = self.top.min(rows.len().saturating_sub(h));

        let label_w = 24usize;
        let mut lines: Vec<Line> = Vec::new();
        for (i, row) in rows.iter().enumerate().skip(self.top).take(h) {
            match row {
                Row::Section(name) => lines.push(Line::from(Span::styled(
                    format!(" {}", name.to_uppercase()),
                    Style::default().fg(Color::DarkGray).add_modifier(Modifier::BOLD),
                ))),
                Row::Field(id) => {
                    let sel = i == self.cursor;
                    let editing = sel && self.edit.is_some();
                    let (text, note) = if editing {
                        (self.edit.clone().unwrap(), String::new())
                    } else {
                        value(&self.cfg, *id, &self.scan)
                    };
                    let val_style = match kind(*id) {
                        Kind::Switch => {
                            if text.starts_with("on") {
                                Style::default().fg(Color::Green)
                            } else {
                                Style::default().fg(Color::DarkGray)
                            }
                        }
                        _ => Style::default().fg(Color::White),
                    };
                    let mut spans = vec![
                        Span::styled(
                            if sel { " ▸ " } else { "   " },
                            Style::default().fg(Color::Cyan),
                        ),
                        Span::styled(
                            format!("{:<w$}", label(*id), w = label_w),
                            if sel {
                                Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD)
                            } else {
                                Style::default()
                            },
                        ),
                        Span::styled(text, val_style),
                    ];
                    if editing {
                        spans.push(Span::styled("▏", Style::default().fg(Color::Cyan)));
                    }
                    if !note.is_empty() {
                        spans.push(Span::styled(note, dim()));
                    }
                    lines.push(Line::from(spans));
                }
            }
        }
        // The frame count belongs next to the directory it counted, and
        // that row is not always on screen.
        if let Some(Row::Field(Id::Input)) = rows.get(self.cursor) {
            if let Some(l) = lines.get_mut(self.cursor - self.top) {
                l.spans.push(frames);
            }
        }
        f.render_widget(Paragraph::new(lines), inner);
    }

    fn draw_help(&self, f: &mut Frame, area: Rect) {
        let block = Block::bordered().title(Span::styled(" what it does ", dim())).border_style(dim());
        let inner = block.inner(area);
        f.render_widget(block, area);
        let text = self.current().map(help).unwrap_or("");
        let lines: Vec<Line> = super::dash::wrap(text, inner.width.max(1) as usize)
            .into_iter()
            .map(|s| Line::from(Span::styled(s, Style::default().fg(Color::Gray))))
            .collect();
        f.render_widget(Paragraph::new(lines), inner);
    }

    fn draw_cmd(&self, f: &mut Frame, area: Rect) {
        let block = Block::bordered()
            .title(Span::styled(" the same run from the command line ", dim()))
            .border_style(dim());
        let inner = block.inner(area);
        f.render_widget(block, area);
        let mut lines: Vec<Line> = super::dash::wrap(&self.cfg.command_line(), inner.width.max(1) as usize)
            .into_iter()
            .map(|s| Line::from(Span::styled(s, Style::default().fg(Color::Cyan))))
            .collect();
        if let Some(e) = &self.error {
            lines.push(Line::from(Span::styled(
                format!("cannot start: {e}"),
                Style::default().fg(Color::Red),
            )));
        }
        f.render_widget(Paragraph::new(lines), inner);
    }

    fn draw_footer(&self, f: &mut Frame, area: Rect) {
        let editing = self.edit.is_some();
        let keys = if editing {
            "↑↓ field   type to edit   ⏎ accept   esc cancel"
        } else {
            "↑↓ field   ←→ change   type to edit   ⏎ start   esc quit"
        };
        f.render_widget(
            Paragraph::new(Line::from(vec![Span::raw(" "), Span::styled(keys, dim())])),
            area,
        );
    }
}

fn dim() -> Style {
    Style::default().fg(Color::DarkGray)
}

/// A field's value as it is drawn, plus a greyed note after it.
fn value(cfg: &Cfg, id: Id, scan: &(String, std::result::Result<usize, String>)) -> (String, String) {
    let on = |b: bool| String::from(if b { "on" } else { "off" });
    let _ = scan;
    match id {
        Id::Input => (cfg.input.clone(), String::new()),
        Id::Sequence => (
            String::from(if cfg.untracked { "untracked — fixed tripod" } else { "tracked — equatorial mount" }),
            String::new(),
        ),
        Id::Limit => (
            if cfg.limit == 0 { String::from("all") } else { cfg.limit.to_string() },
            String::new(),
        ),
        Id::Stack => (on(cfg.stack), String::new()),
        Id::Output => (cfg.output.clone(), String::new()),
        Id::Stretch => (on(cfg.stretch), String::new()),
        Id::FixedStabilize => (on(cfg.fixed_stabilize), String::new()),
        Id::FixedSearch => (format!("{} px", cfg.fixed_search), String::new()),
        Id::FixedAnomaly => (String::from(anomaly_name(cfg.fixed_anomaly)), String::new()),
        Id::Flatten => (on(cfg.flatten), String::new()),
        Id::ResidualSurface => (on(cfg.residual_surface), String::new()),
        Id::ScatterComp => (
            format!("{:.2}", cfg.scatter_comp),
            String::from(if cfg.scatter_comp == 0.0 { "  off" } else { "" }),
        ),
        Id::DumpCorrection => (
            if cfg.dump_correction.trim().is_empty() { String::from("off") } else { cfg.dump_correction.clone() },
            String::new(),
        ),
        Id::Export => (on(cfg.export), String::new()),
        Id::ExportDir => (cfg.export_dir.clone(), String::new()),
        Id::ExportWindow => (
            format!("{} frames", cfg.export_window),
            String::from(if cfg.export_window == 1 { "  no temporal reduction" } else { "" }),
        ),
        Id::ExportStabilize => (on(cfg.export_stabilize), String::new()),
        Id::ExportDeflicker => (on(cfg.export_deflicker), String::new()),
        Id::ExportTransients => (on(cfg.export_transients), String::new()),
        Id::ExportCloudGuard => (on(cfg.export_cloud_guard), String::new()),
        Id::SkyTolerance => (format!("x{:.2}", cfg.sky_tolerance), String::new()),
        Id::MaxForeground => (
            if cfg.max_foreground < 0.0 {
                String::from("auto")
            } else {
                format!("{:.0}%", 100.0 * cfg.max_foreground)
            },
            String::new(),
        ),
        Id::MaxStars => (cfg.max_stars.to_string(), String::new()),
    }
}

/// Counts the frames the run would pick up, with the same extension rule
/// the pipeline uses.
fn count_arw(dir: &Path) -> std::result::Result<usize, String> {
    if dir.as_os_str().is_empty() {
        return Err("no frames directory given".into());
    }
    let rd = std::fs::read_dir(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    Ok(rd
        .filter_map(|e| e.ok())
        .filter(|e| {
            e.path()
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("ARW"))
                .unwrap_or(false)
        })
        .count())
}
