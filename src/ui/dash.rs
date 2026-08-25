//! The live dashboard: run summary at the top, one progress bar per pass,
//! and the log underneath. It owns the terminal from the moment
//! `ui::start_dashboard` is called until `ui::shutdown` joins it back.

use super::{State, Ui};
use ratatui::crossterm::event::{self, Event, KeyCode, KeyEventKind, KeyModifiers};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Paragraph};
use ratatui::Frame;
use std::sync::atomic::Ordering;
use std::sync::Arc;
use std::time::Duration;

/// Filled and empty cell of a progress bar. Half blocks would read better
/// on a wider terminal but not every font has them; these two are as safe
/// as box drawing gets.
const FILL: char = '█';
const VOID: char = '░';

/// Render loop. Runs on its own thread for the whole run: the pipeline
/// never blocks on drawing, and a pass that goes quiet for a minute still
/// ticks its clock.
pub fn run(ui: Arc<Ui>) {
    let mut term = ratatui::init();
    loop {
        let finished = ui.done.load(Ordering::SeqCst);
        let _ = term.draw(|f| draw(f, &ui, finished));
        // A run the user stopped needs no dismissing: they already asked
        // for the screen back.
        if finished && super::aborted() {
            break;
        }
        let wait = if finished { 200 } else { 80 };
        if event::poll(Duration::from_millis(wait)).unwrap_or(false) {
            if let Ok(ev) = event::read() {
                if handle(&ui, ev, finished) {
                    break;
                }
            }
        }
    }
    ratatui::restore();
}

/// Returns true when the screen should be given back.
fn handle(ui: &Ui, ev: Event, finished: bool) -> bool {
    let Event::Key(k) = ev else { return false };
    if k.kind == KeyEventKind::Release {
        return false;
    }
    let ctrl_c = k.modifiers.contains(KeyModifiers::CONTROL) && k.code == KeyCode::Char('c');
    let quit = ctrl_c || matches!(k.code, KeyCode::Char('q') | KeyCode::Esc);
    if finished {
        return quit || matches!(k.code, KeyCode::Enter | KeyCode::Char(' '));
    }
    if quit {
        super::request_abort();
        let mut st = ui.state.lock().unwrap();
        st.phase = String::from("stopping after the frame in flight...");
        return false;
    }
    let mut st = ui.state.lock().unwrap();
    match k.code {
        KeyCode::Up => scroll_back(&mut st, 1),
        KeyCode::Down => scroll_forward(&mut st, 1),
        KeyCode::PageUp => scroll_back(&mut st, 10),
        KeyCode::PageDown => scroll_forward(&mut st, 10),
        KeyCode::Home => scroll_back(&mut st, usize::MAX / 2),
        KeyCode::End => {
            st.follow = true;
            st.scroll = 0;
        }
        _ => {}
    }
    false
}

fn scroll_back(st: &mut State, n: usize) {
    st.follow = false;
    st.scroll = st.scroll.saturating_add(n);
}

fn scroll_forward(st: &mut State, n: usize) {
    st.scroll = st.scroll.saturating_sub(n);
    if st.scroll == 0 {
        st.follow = true;
    }
}

fn draw(f: &mut Frame, ui: &Ui, finished: bool) {
    let mut st = ui.state.lock().unwrap();
    let area = f.area();

    let head_rows = st.header.len().div_ceil(2).max(1) as u16;
    let task_rows = st.tasks.len().max(1) as u16;
    let [head, bars, log, foot] = Layout::vertical([
        Constraint::Length(head_rows + 2),
        Constraint::Length(task_rows + 2),
        Constraint::Min(3),
        Constraint::Length(1),
    ])
    .areas(area);

    draw_header(f, head, &st);
    draw_tasks(f, bars, &st, ui);
    draw_log(f, log, &mut st);
    draw_footer(f, foot, &st, ui, finished);
}

fn dim() -> Style {
    Style::default().fg(Color::DarkGray)
}

fn key_style() -> Style {
    Style::default().fg(Color::Cyan)
}

fn draw_header(f: &mut Frame, area: Rect, st: &State) {
    let title = Line::from(vec![
        Span::styled(" apilaaa ", Style::default().add_modifier(Modifier::BOLD)),
        Span::styled(concat!("v", env!("CARGO_PKG_VERSION"), " "), dim()),
    ]);
    let block = Block::bordered().title(title).border_style(dim());
    let inner = block.inner(area);
    f.render_widget(block, area);

    // Two pairs to a row: the summary is short labels and short values, and
    // one per row would push the log off a small terminal.
    let col = (inner.width as usize / 2).max(20);
    let mut lines: Vec<Line> = Vec::new();
    for pair in st.header.chunks(2) {
        let mut spans: Vec<Span> = Vec::new();
        for (i, (k, v)) in pair.iter().enumerate() {
            let value = clip(v, col.saturating_sub(10));
            spans.push(Span::styled(format!("{k:<8} "), key_style()));
            spans.push(Span::raw(value.clone()));
            if i + 1 < pair.len() {
                let used = 9 + width(&value);
                spans.push(Span::raw(" ".repeat(col.saturating_sub(used))));
            }
        }
        lines.push(Line::from(spans));
    }
    f.render_widget(Paragraph::new(lines), inner);
}

fn draw_tasks(f: &mut Frame, area: Rect, st: &State, ui: &Ui) {
    let block = Block::bordered().title(Span::styled(" progress ", dim())).border_style(dim());
    let inner = block.inner(area);
    f.render_widget(block, area);
    if st.tasks.is_empty() {
        f.render_widget(
            Paragraph::new(Line::from(Span::styled("starting...", dim()))),
            inner,
        );
        return;
    }

    let label_w = st.tasks.iter().map(|t| width(&t.label)).max().unwrap_or(6).clamp(6, 18);
    // " 312/426 73%" plus whatever the finished passes want to say after
    // the bar. The bar gives the room up, not the note: a bar three times
    // wider than it needs to be says no more than one that fits.
    let count_w = st
        .tasks
        .iter()
        .map(|t| width(&format!(" {}/{} 100%", t.done, t.total)))
        .max()
        .unwrap_or(12);
    let note_w = st.tasks.iter().map(|t| width(&t.note) + 2).max().unwrap_or(0);
    let bar_w = (inner.width as usize)
        .saturating_sub(label_w + 1 + count_w + note_w)
        .clamp(4, 56);

    let lines: Vec<Line> = st
        .tasks
        .iter()
        .map(|t| {
            let colour = if t.finished { Color::Green } else { Color::Cyan };
            let (filled, total_cells) = if t.total > 0 {
                let r = (t.done as f64 / t.total as f64).clamp(0.0, 1.0);
                ((r * bar_w as f64).round() as usize, bar_w)
            } else {
                (0, bar_w)
            };
            let mut spans = vec![
                Span::styled(format!("{:<w$} ", t.label, w = label_w), key_style()),
            ];
            if t.total > 0 {
                spans.push(Span::styled(FILL.to_string().repeat(filled), Style::default().fg(colour)));
                spans.push(Span::styled(VOID.to_string().repeat(total_cells - filled), dim()));
                let pct = 100.0 * t.done as f64 / t.total as f64;
                spans.push(Span::raw(format!(" {}/{} {:>3.0}%", t.done, t.total, pct)));
            } else {
                // Size not known yet: a block sliding along the track says
                // "working" without claiming a fraction that does not exist.
                let phase = (ui.started.elapsed().as_millis() / 90) as usize;
                let win = (bar_w / 5).max(2);
                let pos = phase % (bar_w + win);
                let mut track = String::new();
                for i in 0..bar_w {
                    track.push(if i + win > pos && i < pos { FILL } else { VOID });
                }
                spans.push(Span::styled(track, Style::default().fg(colour)));
                spans.push(Span::raw(format!(" {}", t.done)));
            }
            if !t.note.is_empty() {
                let room = (inner.width as usize).saturating_sub(label_w + 1 + bar_w + count_w + 2);
                spans.push(Span::styled(format!("  {}", clip(&t.note, room)), dim()));
            }
            Line::from(spans)
        })
        .collect();
    f.render_widget(Paragraph::new(lines), inner);
}

fn draw_log(f: &mut Frame, area: Rect, st: &mut State) {
    let mut title = vec![Span::styled(" log ", dim())];
    if st.dropped > 0 {
        title.push(Span::styled(format!("(+{} older dropped) ", st.dropped), dim()));
    }
    if !st.follow {
        title.push(Span::styled("PAUSED ", Style::default().fg(Color::Yellow)));
    }
    let block = Block::bordered().title(Line::from(title)).border_style(dim());
    let inner = block.inner(area);
    f.render_widget(block, area);

    let w = inner.width.max(1) as usize;
    let h = inner.height.max(1) as usize;
    // Only as much of the tail as can possibly fill the panel is wrapped:
    // the buffer holds thousands of lines and all of them are re-wrapped on
    // every frame otherwise.
    let take = h + st.scroll + 4;
    let mut rows: Vec<Line> = Vec::new();
    for line in st.log.iter().rev().take(take).collect::<Vec<_>>().into_iter().rev() {
        for piece in wrap(line, w) {
            rows.push(style_line(&piece));
        }
    }
    let max_scroll = rows.len().saturating_sub(h);
    if st.scroll > max_scroll {
        st.scroll = max_scroll;
        if st.scroll == 0 {
            st.follow = true;
        }
    }
    let end = rows.len().saturating_sub(st.scroll);
    let start = end.saturating_sub(h);
    f.render_widget(Paragraph::new(rows[start..end].to_vec()), inner);
}

/// The log is one stream of plain sentences; colour is only there to make
/// the few lines that report a problem findable in it.
fn style_line(s: &str) -> Line<'static> {
    let t = s.trim_start();
    let style = if t.starts_with("ERROR:") {
        Style::default().fg(Color::Red).add_modifier(Modifier::BOLD)
    } else if t.contains("WARNING") || t.contains("NOT ALIGNED") || t.contains("implausible") {
        Style::default().fg(Color::Yellow)
    } else if t.contains("LOAD FAILED") || t.contains("different dimensions") {
        Style::default().fg(Color::Red)
    } else if s.starts_with("  ") {
        Style::default().fg(Color::Gray)
    } else {
        Style::default()
    };
    Line::from(Span::styled(s.to_string(), style))
}

fn draw_footer(f: &mut Frame, area: Rect, st: &State, ui: &Ui, finished: bool) {
    // Once the run is over the clock stops: what it reports is the run,
    // not how long the finished screen was left up.
    let done_ms = ui.finish_ms.load(Ordering::SeqCst);
    let e = if done_ms > 0 { done_ms / 1000 } else { ui.started.elapsed().as_secs() };
    let clock = format!(" {:02}:{:02}:{:02} ", e / 3600, (e / 60) % 60, e % 60);
    let keys = if finished {
        "q close   ↑↓ PgUp/PgDn scroll"
    } else {
        "q stop   ↑↓ PgUp/PgDn scroll   End follow"
    };
    let failed = ui.failed.load(Ordering::SeqCst);
    let phase = match (finished, super::aborted(), failed) {
        (false, _, _) => st.phase.clone(),
        (true, true, _) => String::from("stopped"),
        (true, false, true) => String::from("failed — see the last line of the log"),
        (true, false, false) => String::from("done"),
    };
    let colour = match (finished, super::aborted() || failed) {
        (true, false) => Color::Green,
        (true, true) => Color::Red,
        _ => Color::Cyan,
    };
    let left = format!("{clock}{phase}");
    let pad = (area.width as usize)
        .saturating_sub(width(&left) + width(keys) + 1);
    f.render_widget(
        Paragraph::new(Line::from(vec![
            Span::styled(clock, Style::default().fg(colour).add_modifier(Modifier::BOLD)),
            Span::raw(clip(&phase, area.width as usize / 2)),
            Span::raw(" ".repeat(pad)),
            Span::styled(keys, dim()),
        ])),
        area,
    );
}

/// Display width. Everything the pipeline prints is Latin text with a few
/// box and degree signs in it, so counting characters is exact here and
/// avoids a dependency on a width table.
fn width(s: &str) -> usize {
    s.chars().count()
}

fn clip(s: &str, max: usize) -> String {
    if width(s) <= max {
        return s.to_string();
    }
    s.chars().take(max.saturating_sub(1)).collect::<String>() + "…"
}

/// Hard-wraps at the panel's width, breaking on a space when there is one
/// in the last third of the row so a filename is not split down the middle.
/// Continuation rows are indented, which is what keeps a wrapped line from
/// reading as two entries.
pub fn wrap(s: &str, w: usize) -> Vec<String> {
    if w < 8 {
        return vec![s.chars().take(w).collect()];
    }
    let chars: Vec<char> = s.chars().collect();
    if chars.len() <= w {
        return vec![s.to_string()];
    }
    let mut out = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let indent = if out.is_empty() { 0 } else { 4 };
        let room = w - indent;
        let end = (i + room).min(chars.len());
        let mut cut = end;
        if end < chars.len() {
            if let Some(p) = (i + room * 2 / 3..end).rev().find(|&k| chars[k] == ' ') {
                cut = p;
            }
        }
        let mut row = " ".repeat(indent);
        row.extend(&chars[i..cut]);
        out.push(row);
        i = if cut < end && chars[cut] == ' ' { cut + 1 } else { cut };
    }
    out
}
