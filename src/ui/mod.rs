//! Console front end.
//!
//! Two things live here. The **setup screen** (`wizard`) is what `apilaaa`
//! shows when it is started with no arguments at all and there is a
//! terminal on the other end: the same options the command line takes,
//! laid out as a form, with the equivalent command line written out at the
//! bottom so the next run can skip it. The **dashboard** (`dash`) then
//! replaces the scrolling log for the duration of that run, with one
//! progress bar per pass and the log kept underneath.
//!
//! Everything the pipeline reports goes through this module: `say!` is the
//! line channel that used to be `println!`, and the `task_*` calls drive
//! the bars. With arguments on the command line — or with no terminal, as
//! under a pipe, in a CI job or through `nohup` — the sink stays `Plain`
//! and every `say!` is a `println!` again, so scripted use is byte for byte
//! what it was.

pub mod dash;
pub mod wizard;

use std::collections::VecDeque;
use std::io::IsTerminal;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};
use std::thread::JoinHandle;
use std::time::Instant;

/// Log lines kept in memory for the dashboard's scrollback. Older ones are
/// dropped and counted, and the panel says how many are missing rather than
/// pretending it holds the whole run.
const LOG_CAP: usize = 8192;

/// One progress bar: a pass over the frames.
pub struct Task {
    /// Stable identifier the pipeline addresses the bar by.
    pub key: &'static str,
    pub label: String,
    pub done: u64,
    pub total: u64,
    pub finished: bool,
    /// Written after the bar once the pass is over ("298 frames, 41.2s").
    pub note: String,
}

/// Everything the dashboard draws. Written by the pipeline threads, read by
/// the render thread.
#[derive(Default)]
pub struct State {
    /// Run summary, drawn as `key value` pairs at the top.
    pub header: Vec<(String, String)>,
    pub tasks: Vec<Task>,
    pub log: VecDeque<String>,
    /// Log lines evicted by `LOG_CAP`.
    pub dropped: usize,
    /// Every line that looked like a warning, kept so the summary printed
    /// after the terminal is restored can repeat them: they are the part of
    /// the log a user actually needs once the screen is gone.
    pub warnings: Vec<String>,
    /// What the pipeline is doing right now, shown in the footer.
    pub phase: String,
    /// Whether the log panel sticks to the tail. Scrolling up turns it off,
    /// `End` turns it back on.
    pub follow: bool,
    /// Rows scrolled back from the tail when not following.
    pub scroll: usize,
}

pub struct Ui {
    pub state: Mutex<State>,
    pub started: Instant,
    /// Wall time of the run itself, in milliseconds, stamped the moment the
    /// pipeline returns. However long the finished screen is then left on
    /// the terminal is not part of it.
    pub finish_ms: AtomicU64,
    /// Set when the pipeline has returned an error, so the finished screen
    /// says so instead of saying "done" over a run that produced nothing.
    pub failed: AtomicBool,
    /// Set when the pipeline has returned. The render thread stops
    /// redrawing on its own and waits for the user to dismiss the screen.
    pub done: AtomicBool,
}

enum Sink {
    /// Plain `println!`, exactly as before this module existed.
    Plain,
    Tui(Arc<Ui>),
}

static SINK: OnceLock<Sink> = OnceLock::new();
static ABORT: AtomicBool = AtomicBool::new(false);
static RENDER: Mutex<Option<JoinHandle<()>>> = Mutex::new(None);

fn sink() -> &'static Sink {
    SINK.get_or_init(|| Sink::Plain)
}

/// Whether the setup screen should be offered: no arguments at all, and a
/// terminal on both stdin and stdout. `apilaaa > log` or `apilaaa | tee`
/// therefore still runs the defaults straight through.
pub fn should_prompt() -> bool {
    std::env::args_os().count() == 1
        && std::io::stdin().is_terminal()
        && std::io::stdout().is_terminal()
}

/// Takes over the terminal for the duration of the run. Called once, before
/// the pipeline starts; after it every `say!` goes to the log panel instead
/// of the screen.
pub fn start_dashboard(header: Vec<(String, String)>) {
    let ui = Arc::new(Ui {
        state: Mutex::new(State {
            header,
            follow: true,
            ..Default::default()
        }),
        started: Instant::now(),
        finish_ms: AtomicU64::new(0),
        failed: AtomicBool::new(false),
        done: AtomicBool::new(false),
    });
    if SINK.set(Sink::Tui(ui.clone())).is_err() {
        return;
    }
    let handle = std::thread::spawn(move || dash::run(ui));
    *RENDER.lock().unwrap() = Some(handle);
}

/// Gives the terminal back and prints the run's summary on the normal
/// screen: what each pass got through, every warning, and the total time.
/// The dashboard's log dies with the alternate screen, so this is what is
/// left in the scrollback afterwards. `error` is what the pipeline
/// returned, so the screen can say the run failed while it is still up.
pub fn shutdown(error: Option<String>) {
    let Sink::Tui(ui) = sink() else { return };
    if let Some(e) = &error {
        if !aborted() {
            say(format!("ERROR: {e}"));
            ui.failed.store(true, Ordering::SeqCst);
        }
    }
    let elapsed = ui.started.elapsed().as_secs_f32();
    ui.finish_ms.store((elapsed * 1000.0) as u64, Ordering::SeqCst);
    ui.done.store(true, Ordering::SeqCst);
    if let Some(h) = RENDER.lock().unwrap().take() {
        let _ = h.join();
    }
    let st = ui.state.lock().unwrap();
    for (k, v) in &st.header {
        println!("{k:<10} {v}");
    }
    for t in &st.tasks {
        let note = if t.note.is_empty() { String::new() } else { format!("  {}", t.note) };
        println!("{:<10} {}/{}{}", t.label, t.done, t.total, note);
    }
    if !st.warnings.is_empty() {
        println!("{} warning(s):", st.warnings.len());
        for w in st.warnings.iter().take(20) {
            println!("  {w}");
        }
        if st.warnings.len() > 20 {
            println!("  ... and {} more", st.warnings.len() - 20);
        }
    }
    println!("total time: {elapsed:.1}s");
}

/// A line of report. `say!` is the macro form and is what call sites use.
pub fn say(line: String) {
    match sink() {
        Sink::Plain => println!("{line}"),
        Sink::Tui(ui) => {
            let mut st = ui.state.lock().unwrap();
            if is_warning(&line) && st.warnings.len() < 500 {
                let w = line.trim().to_string();
                st.warnings.push(w);
            }
            st.log.push_back(line);
            if st.log.len() > LOG_CAP {
                st.log.pop_front();
                st.dropped += 1;
            }
        }
    }
}

/// Lines worth repeating once the screen is gone. The pipeline marks its
/// own in capitals already, so no call site has to say twice that it is
/// reporting a problem — and nothing is inferred from ordinary wording:
/// "drift not measurable" is the correct, expected reading on a sequence
/// framed on pure sky, not something to warn about.
fn is_warning(line: &str) -> bool {
    let l = line.trim_start();
    l.contains("WARNING")
        || l.contains("LOAD FAILED")
        || l.contains("NOT ALIGNED")
        || l.contains("different dimensions")
        || l.contains("implausible alignment")
}

/// Sets one field of the run summary at the top of the dashboard,
/// replacing it if it is already there.
pub fn head(key: &str, value: String) {
    let Sink::Tui(ui) = sink() else { return };
    let mut st = ui.state.lock().unwrap();
    match st.header.iter_mut().find(|(k, _)| k == key) {
        Some(slot) => slot.1 = value,
        None => st.header.push((key.to_string(), value)),
    }
}

/// What the pipeline is doing right now, for the footer.
pub fn phase(text: &str) {
    let Sink::Tui(ui) = sink() else { return };
    ui.state.lock().unwrap().phase = text.to_string();
}

/// Opens (or restarts) a progress bar. `total` of 0 draws an indeterminate
/// bar: the work is under way but its size is not known yet.
pub fn task_begin(key: &'static str, label: &str, total: u64) {
    let Sink::Tui(ui) = sink() else { return };
    let mut st = ui.state.lock().unwrap();
    match st.tasks.iter_mut().find(|t| t.key == key) {
        Some(t) => {
            t.label = label.to_string();
            t.done = 0;
            t.total = total;
            t.finished = false;
            t.note.clear();
        }
        None => st.tasks.push(Task {
            key,
            label: label.to_string(),
            done: 0,
            total,
            finished: false,
            note: String::new(),
        }),
    }
}

/// Advances a bar by `n`.
pub fn task_add(key: &'static str, n: u64) {
    with_task(key, |t| t.done += n);
}

/// Sets a bar's absolute position.
pub fn task_set(key: &'static str, done: u64) {
    with_task(key, |t| t.done = done);
}

/// Closes a bar. `note` is what stays next to it, and in the summary.
pub fn task_end(key: &'static str, note: String) {
    with_task(key, |t| {
        t.finished = true;
        if t.total > 0 {
            t.done = t.done.min(t.total);
        }
        t.note = note;
    });
}

fn with_task(key: &'static str, f: impl FnOnce(&mut Task)) {
    let Sink::Tui(ui) = sink() else { return };
    let mut st = ui.state.lock().unwrap();
    if let Some(t) = st.tasks.iter_mut().find(|t| t.key == key) {
        f(t);
    }
}

/// Whether the user has asked to stop. The pipeline checks it at frame
/// boundaries, so a run gives up between frames and never in the middle of
/// writing one.
pub fn aborted() -> bool {
    ABORT.load(Ordering::Relaxed)
}

pub fn request_abort() {
    ABORT.store(true, Ordering::Relaxed);
}

/// One line of report, formatted like `println!`. Replaces it everywhere in
/// the pipeline so the same call works with or without the dashboard.
#[macro_export]
macro_rules! say {
    ($($t:tt)*) => { $crate::ui::say(format!($($t)*)) };
}
