# apilaaa

**A RAW astrophotography stacker that tells the sky apart from the camera.**

It aligns a night-sky sequence on its own stars, measures and removes the
optical defects of the camera system — lens halo, vignetting, flare wedges,
sensor banding, horizon glow — and writes the result as a linear DNG that any
raw developer can pick up. The very same model is then reused to export a
flicker-free, noise-reduced timelapse from the exact same exposures, so a single
night of frames yields both a deep still and a finished sequence.

| | |
|---|---|
| **Input** | A directory of Sony `.ARW` frames (one night, tracked or on a fixed tripod) |
| **Output** | A linear RGB DNG stack, an optional cleaned DNG sequence, an optional "what was removed" layer |
| **Language** | Rust 2021, no GPU, parallel over frames with Rayon |
| **Platform** | Linux (tested), CPU-bound, RAM-planned at startup |
| **Licence** | GPL-3.0-or-later |

## The idea

A tracked sky sequence contains two superimposed signals with opposite
behaviour: **the sky drifts across the sensor** between frames, whereas **the
lens halo, the vignetting, the flare wedges and the sensor banding stay fixed to
the sensor**. Because those defects are additive and stationary in sensor
coordinates, they can be estimated from the temporal median of the sequence and
subtracted *before* averaging — which is precisely the order of operations the
tool follows.

That single observation is what makes the difference between this and a plain
stacker. A plain stacker averages the defects in along with the sky and leaves
you to paint them out by hand later; here they are modelled explicitly, removed
in sensor space where they actually live, and the removed layer can be written
out as its own file so you can prove nothing else went with it.

For an **untracked** sequence (`--fixed-tripod`) the roles are swapped — the
landscape is what stays put and the sky is what moves — and the pipeline adapts:
star alignment is replaced by tripod-drift measurement, the per-frame anomaly is
stiffened so the drifting Milky Way is not mistaken for a defect, and the stack
becomes a star-trail image, which is what averaging a fixed sequence means.

## What you get

| Output | Flag | What it is |
|---|---|---|
| **The stack** | `-o stacked.dng` | Every selected frame aligned, cleaned and averaged. Linear RGB DNG, demosaic already done, white balance already baked in, per-channel stretch baked in (or not, with `--no-stretch`). |
| **The clean sequence** | `--export-clean DIR` | One `<name>_clean.dng` per frame: cleaned with the same model, stabilized into the same coordinate system, cropped identically, deflickered against the stack, temporally denoised, with meteors and satellites preserved. Ready to feed straight to a video encoder. |
| **The removed layer** | `--dump-correction layer.dng` | Exactly what was subtracted, plus the preserved sky pedestal, under the same crop and stretch as the stack. The audit file: it should contain defect geometry and nothing else. |

## Requirements

- A Rust toolchain (edition 2021).
- [`exiftool`](https://exiftool.org/) on `PATH`, optional. Every DNG already
  carries the EXIF that identifies the shot — camera, lens, focal length,
  aperture, exposure, ISO and capture time — written natively with no external
  tool, which is what a raw developer matches a lens profile on. exiftool is
  only needed to carry across everything *beyond* that: MakerNotes, XMP and
  ICC. If it is missing the run warns and the native EXIF stays.
- Linux is the tested platform. The RAM planner reads `/proc/meminfo` and falls
  back to a conservative 4 GiB assumption elsewhere.

## Build

```sh
cargo build --release
```

The release profile is tuned for throughput (`opt-level = 3`, thin LTO, a single
codegen unit), which matters because the whole pipeline is CPU-bound and
processes every frame twice.

## Quick start

Point the tool at a directory of Sony `.ARW` files and give it an output path:

```sh
./target/release/apilaaa --input res --output stacked.dng
```

Files are sorted by name, which for in-camera numbering is chronological order,
and the first one becomes the reference frame that defines the output geometry.
Everything downstream — alignment, cropping, the timelapse export — is expressed
in that reference's coordinate system.

Before committing a whole night, try the first few dozen frames:

```sh
./target/release/apilaaa -i res -o test.dng --limit 40
```

`--limit` cuts the run to the first N files, which is enough to see whether the
stars are detected, the alignment locks and the defect model looks sane.

## Tracked or fixed tripod?

This is the first decision of a run, and the only one the tool cannot make for
you: **was the camera following the sky, or was it not?**

| Your setup | What is fixed on the sensor | Run it as |
|---|---|---|
| The camera rode a star tracker or an equatorial mount, so the stars stayed put across the whole session | The **sky**; the landscape, if any, drifts through the frame | default, no flag |
| The camera sat on a plain tripod, head locked, nothing motorized | The **landscape**; the stars sweep across it | `--fixed-tripod` |

```sh
apilaaa -i res -o stacked.dng                  # tracked
apilaaa -i res -o trails.dng --fixed-tripod    # untracked
```

If you are unsure which you shot, open the first and the last frame of the night
side by side. Whichever of the two — stars or landscape — is in the same place in
both is what the sequence is fixed to, and that is the mode you want.

### What the flag changes

Both modes are the same two-pass pipeline; `--fixed-tripod` swaps the stages
that depend on which signal is the stationary one.

| Stage | Tracked (default) | `--fixed-tripod` |
|---|---|---|
| Registration | Triangle matching over the detected stars, RANSAC, least-squares similarity — rotation *and* translation | Two-scale normalized cross-correlation of the high-passed luminance, restricted to the landscape cells; translation only, sub-pixel, capped at `--fixed-search` px |
| Foreground mask | One per frame | A single consensus mask: a cell is landscape when more than half the frames see it as such |
| Per-frame anomaly | `Full` | `--fixed-anomaly`, `coarse` by default, so the drifting Milky Way is not read as a defect |
| `--stack-max-foreground` | Defaults to `0.0`: only clear-sky frames are averaged | Defaults to `1.0`: the landscape is part of the picture, never a reason to drop a frame |
| The stack (`-o`) | A deep still — pinpoint stars, the landscape smeared by its own motion | A star-trail image — fixed landscape, the stars drawn as arcs |
| `--export-clean` window | Combined on the stabilized frames directly: registration already put the sky in place | Warped onto the sky through the chain of star fits between consecutive frames, then combined |
| Export levels and stretch | Taken from the stack | Taken from one representative frame (the sky level closest to the session median), because a trail stack has no star cores to match against; the run names the frame it picked |

Everything else — the defect model, the frame selection by sky level, the
deflickering, the transient rescue, the dawn ramp — is identical in both modes,
which is why the flags in [Options](#options) apply to both unless the table
says otherwise.

### Picking the wrong one

Neither mistake damages anything — inputs are only ever read — but neither gives
a usable result. The second announces itself in the console; the first you see
in the stack:

- **Default on an untracked night.** The star alignment locks anyway, because
  stars are exactly what it matches: you get the sky registered and the
  landscape smeared across the frame, and the anomaly, at `Full` freedom over a
  Milky Way that is moving through the sensor, takes part of it with it. The
  symptom is a stack with round stars over a doubled or streaked horizon —
  every frame masks its landscape in a different place, so that region ends up
  a partial average rather than a picture.
- **`--fixed-tripod` on a tracked night.** The correlation goes looking for a
  stationary landscape and finds one that is moving, usually much further than
  the `--fixed-search` range, so the peak falls under the acceptance threshold
  (`0.25`) and the run prints `drift not measurable, identity` frame after
  frame, stacking everything unregistered. A framing with under 2 % landscape
  never gets that far and announces it at the start:
  `tripod drift: NOT measurable (landscape 0.4% of the frame)`.

A third case is legitimate rather than a mistake: an untracked sequence on a
tripod you are confident never moved at all. `--fixed-tripod --fixed-no-stabilize`
skips the drift measurement, keeps every frame in sensor coordinates and exports
uncropped.

## Recipes

### Deep stack from a tracked sequence

```sh
./target/release/apilaaa -i res -o stacked.dng
```

The default path: align on stars, select the clear-sky frames, fit the defect
model on their temporal median, clean and average. Open `stacked.dng` in
darktable, RawTherapee, Lightroom or anything else that reads LinearRaw DNG; the
stretch is already baked in, so it should look reasonable before you touch a
slider.

### Stack *and* clean timelapse in one run

```sh
./target/release/apilaaa -i res -o stacked.dng --export-clean res_clean
```

Both passes are shared, so the sequence costs the export time only, not a second
alignment. Every frame written by `--export-clean` uses the same stretch as the
stack, which means the still and the sequence are directly comparable — and
frames rejected from the stack (twilight, clouds, frames with the horizon in
them) are still exported, because a frame that is bad for averaging is usually
perfectly fine in a moving sequence.

There is no export-only mode, and that is deliberate: the stack is the reference
the deflickering matches every frame against, so it is always computed and
written. If the still is not what you came for, point `-o` at a throwaway path
and keep the sequence.

To turn the result into a video, develop the DNGs however you like and encode
them — see [Encoding the sequence](#encoding-the-sequence).

### Untracked sequence on a fixed tripod

```sh
./target/release/apilaaa -i res -o trails.dng --fixed-tripod --export-clean res_clean
```

One run gives both halves of what a fixed-tripod night is good for: a
star-trail still in `trails.dng` and a stabilized, cleaned sequence in
`res_clean` where the landscape stands still and the sky turns over it. The
residual drift — wind, legs settling into soft ground — is measured off the
landscape and taken out, which is what keeps the horizon from jittering through
the sequence; see [What the flag changes](#what-the-flag-changes) for the full
list of what the mode swaps.

Tune it with `--fixed-search` (how far the drift is looked for, in px),
`--fixed-anomaly coarse|none|full` (how much per-frame structure may be
subtracted before it starts eating the Milky Way) and `--fixed-no-stabilize`
(skip the drift measurement entirely).

If the trails come out as dashes rather than continuous arcs, that is the
interval between exposures, not the tool: the stack draws exactly what the
sensor saw, gaps included.

### Audit what was actually removed

```sh
./target/release/apilaaa -i res -o stacked.dng --dump-correction layer.dng
```

Then open `stacked.dng` and `layer.dng` side by side. See
[Auditing the correction](#auditing-the-correction) below for the full
procedure and what a healthy removed layer looks like.

### Leave the whole development to your raw developer

```sh
./target/release/apilaaa -i res -o stacked.dng --no-stretch
```

By default the output carries a per-channel black/white point rescale so the DNG
looks reasonable the moment it opens. The data stays linear either way — the
stretch is a rescale, not a curve — but `--no-stretch` keeps the native sensor
scale, which is what you want when you intend to *measure* the result: comparing
two runs numerically, reading the removed layer in sensor units, or feeding the
file to a photometry tool.

### Triage a night before committing to it

```sh
./target/release/apilaaa -i res -o test.dng --limit 40
./target/release/apilaaa -i res -o test.dng --limit 40 --export-clean test_clean --export-window 1
```

`--limit 40` stops after the first forty files, which is enough to see the star
count, whether the alignment locks, and what the defect model looks like.
Adding `--export-window 1` on top removes the temporal denoising, the most
expensive part of the export, so a preview sequence costs roughly one frame's
work per frame.

### Too few frames survive the selection

If the run stops with `no frame passes the stack selection`, or the frame count
that reaches the stack is far below what you shot, the two selection knobs are
what to reach for:

```sh
# a night that drifts in brightness (moonrise, twilight, thin cloud)
./target/release/apilaaa -i res -o stacked.dng --stack-sky-tolerance 2.5

# a tracked framing that always has some horizon or trees in it
./target/release/apilaaa -i res -o stacked.dng --stack-max-foreground 0.3
```

`--stack-sky-tolerance` widens the accepted sky-level band around the session
median; `--stack-max-foreground` admits frames that contain foreground and masks
those pixels out instead of dropping the frame. Both only change what is
*averaged* — `--export-clean` exports rejected frames regardless, because a frame
that is bad for averaging is usually fine in a moving sequence.

### The alignment cannot find its stars

```sh
./target/release/apilaaa -i res -o stacked.dng --max-stars 80
```

`too few stars in the reference` means the first file, not the sequence, is the
problem — it is the reference that defines the geometry. Raising `--max-stars`
gives the triangle matching more candidates on a sparse or short-exposure field;
lowering it speeds up matching on a dense one. If the reference frame itself is
unusable (cloud, a plane crossing it, a mistimed first exposure), move it out of
the directory rather than tuning around it.

### Keep or erase the things that move

Meteors, satellites and planes are preserved by default: the temporal window
would otherwise average them away, so pixels that deviate from the combination
by more than 4σ in elongated components are written back from the original
frame.

```sh
# a meteor shower — the default, spelled out
./target/release/apilaaa -i res -o stacked.dng --export-clean res_clean

# an airport approach path you would rather not see
./target/release/apilaaa -i res -o stacked.dng --export-clean res_clean --export-no-transients
```

### Export in sensor coordinates

```sh
./target/release/apilaaa -i res -o stacked.dng --export-clean res_raw --export-no-stabilize
```

Every frame keeps its own framing and full size, uncropped. Useful when you want
to stabilize elsewhere, or to see the defect correction on its own without the
alignment moving the picture underneath it. Temporal denoising is disabled with
it, because the sliding window needs the frames to be registered.

### When the model misbehaves

Three knobs, in order of severity:

```sh
apilaaa -i res -o stacked.dng --scatter-comp 0          # halo removed, contrast not rescaled
apilaaa -i res -o stacked.dng --no-residual-surface     # halo + gradient only
apilaaa -i res -o stacked.dng --no-flatten              # no correction at all
```

`--scatter-comp 0` is the mildest: the veil is still subtracted, but the local
contrast in the veiled area is left alone instead of being divided back up,
which is the right call when the compensation starts amplifying noise near the
corners. `--no-residual-surface` drops the surface, bands and spokes and keeps
the parametric halo and gradient. `--no-flatten` turns the correction off
entirely — and with it `--export-clean`, which refuses to run without a model.

## Options

### Core

| Option | Default | Purpose |
|---|---|---|
| `-i, --input <DIR>` | `res` | Directory holding the input `.ARW` files. |
| `-o, --output <DNG>` | `stacked.dng` | Output DNG for the stack. |
| `--max-stars <N>` | `40` | Stars used per frame for alignment. |
| `--limit <N>` | — | Process at most N frames; useful for quick tests. |
| `--no-stretch` | off | Keep the native sensor scale instead of baking a stretch. |

### Defect correction

| Option | Default | Purpose |
|---|---|---|
| `--no-flatten` | off | Disable the whole defect model; stack the frames as they are. |
| `--no-residual-surface` | off | Keep only the parametric halo + gradient model, skipping the residual surface, bands and spokes. |
| `--scatter-comp <BETA>` | `1.0` | Contrast compensation in the veiled areas; `0` turns it off. |
| `--dump-correction <DNG>` | — | Write the removed layer, with the same crop and stretch as the output. |

### Frame selection

| Option | Default | Purpose |
|---|---|---|
| `--stack-sky-tolerance <F>` | `1.6` | Frames whose sky level leaves `[median/F, median·F]` are excluded from the stack. |
| `--stack-max-foreground <FRAC>` | `0.0` (`1.0` under `--fixed-tripod`) | Frames with more than this fraction of foreground are excluded from the stack. |

### Timelapse export

| Option | Default | Purpose |
|---|---|---|
| `--export-clean <DIR>` | — | Export the stabilized, cleaned timelapse sequence. Requires the correction to be enabled. |
| `--export-window <N>` | `7` | Odd size of the temporal window used to reduce noise on export; `1` disables it. |
| `--export-no-stabilize` | off | Export in sensor coordinates, uncropped (this also disables temporal denoising). |
| `--export-no-deflicker` | off | Skip the per-channel level matching on export. |
| `--export-no-transients` | off | Let the temporal combination erase meteors, satellites and planes. |
| `--export-no-cloud-guard` | off | Let a cloud take part in the temporal window and in the frame's level statistics, as it did before the guard existed. |

### Untracked sequences

| Option | Default | Purpose |
|---|---|---|
| `--fixed-tripod` | off | The camera stayed fixed: measure tripod drift instead of star alignment, use a consensus landscape mask, combine the export window on the sky, and produce a star-trail stack. |
| `--fixed-no-stabilize` | off | Assume the camera was perfectly still; export in sensor coordinates, uncropped. |
| `--fixed-search <PX>` | `64` | Maximum tripod drift searched for, in sensor pixels. |
| `--fixed-anomaly <MODE>` | `coarse` | Freedom given to the per-frame anomaly: `coarse` (broad structure only), `none` (defect model alone), `full` (as in a tracked sequence — will eat part of a drifting Milky Way). |

## Command reference

Everything runnable in this repository, in one place: the build, the tool
itself, the audit scripts and the diagnostic invocations. `apilaaa` below is
shorthand for `./target/release/apilaaa`; the flags themselves are catalogued
under [Options](#options).

### Build and run

| Command | What it does |
|---|---|
| `cargo build --release` | Builds the optimized binary at `target/release/apilaaa`. A debug build works but is roughly an order of magnitude slower, and the pipeline reads every frame twice. |
| `cargo run --release -- -i res -o stacked.dng` | The same run without leaving Cargo; everything after `--` goes to the program. |
| `./target/release/apilaaa --help` | The generated flag list with its defaults. The [Options](#options) tables above are its annotated form. |
| `./target/release/apilaaa -i res -o stacked.dng 2>&1 \| tee run.log` | Keeps the narration. The console output is the run's record — see [Reading the console output](#reading-the-console-output). |

### Picking the run

| If you want to… | Command |
|---|---|
| Stack one tracked night | `apilaaa -i res -o stacked.dng` |
| Stack **and** get the timelapse | `apilaaa -i res -o stacked.dng --export-clean res_clean` |
| Work an untracked, fixed-tripod night | `apilaaa -i res -o trails.dng --fixed-tripod --export-clean res_clean` |
| Check the run before spending the night on it | `apilaaa -i res -o test.dng --limit 40` |
| Get a preview sequence fast | `apilaaa -i res -o test.dng --limit 40 --export-clean test_clean --export-window 1` |
| Keep the native sensor scale | `apilaaa -i res -o stacked.dng --no-stretch` |
| See what the correction removed | `apilaaa -i res -o stacked.dng --dump-correction layer.dng` |
| Loosen the frame selection | `apilaaa -i res -o stacked.dng --stack-sky-tolerance 2.5 --stack-max-foreground 0.3` |
| Back the correction off one stage | `apilaaa -i res -o stacked.dng --no-residual-surface` |
| Turn the correction off entirely | `apilaaa -i res -o stacked.dng --no-flatten` |

### Auditing an output

The DNGs this tool writes are ordinary linear DNGs, so an audit needs nothing
beyond the tool itself and whatever raw developer you already use. Run the same
session three ways, in linear scale so the differences stay in the units the
correction works in:

```sh
apilaaa -i res -o corrected.dng   --no-stretch --dump-correction layer.dng
apilaaa -i res -o uncorrected.dng --no-stretch --no-flatten
apilaaa -i res -o /tmp/p.dng      --no-stretch --no-residual-surface \
        --dump-correction layer_param.dng
```

`layer.dng` is everything the correction removed. `layer_param.dng` is what the
parametric stage alone removed, so `layer − layer_param` is the residual
surface's own contribution. Opening the four under one identical stretch is
what tells a defect from a piece of sky: anything in the removed layer that
looks like structure — a dark nebula, a dust lane — is the correction reaching
too far.

See [Auditing the correction](#auditing-the-correction) for how to read the
result.

### Encoding the sequence

`--export-clean` writes one `<name>_clean.dng` per frame. Develop them however
you like, then hand the result to a video encoder:

```sh
ffmpeg -framerate 24 -pattern_type glob -i 'res_clean/*.jpg' \
       -c:v libx264 -crf 16 -pix_fmt yuv420p timelapse.mp4
```

### Diagnostic runs

Three environment variables expose the internals without changing a pixel of
the output. Only their presence matters, except for `APILAAA_DEBUG_DIR`, whose
value is the destination directory — it is created if it does not exist.

```sh
APILAAA_DEBUG=1 ./target/release/apilaaa -i res -o stacked.dng
APILAAA_DEBUG_DIR=/tmp/dbg ./target/release/apilaaa -i res -o stacked.dng
APILAAA_DEBUG_CORR=1 ./target/release/apilaaa -i res -o stacked.dng --limit 5
```

`APILAAA_DEBUG_CORR` prints a numeric map per frame, so keep it under
`--limit`. What each variable dumps is described under
[Diagnostics](#diagnostics) below.

## How the flags interact

The flags are not fully independent, and the run tells you when one of them has
been overruled. The rules worth knowing before a long night:

- **`--export-clean` needs the correction.** Combined with `--no-flatten` the
  run stops immediately with `--export-clean requires the correction to be
  enabled (drop --no-flatten)`. There is no export-only mode either: the stack
  is the reference the deflickering matches against, so it is always computed.
- **`--dump-correction` with `--no-flatten` is a warning, not an error.** There
  is nothing to dump, so the run continues and prints
  `--dump-correction ignored (correction disabled)`.
- **`--export-window` is rounded up to the next odd number** and floored at
  `1`; `2` behaves as `3`. `1` disables the temporal noise reduction.
- **An unstabilized tracked export forces the window to `1`.** The sliding
  window combines registered frames, so `--export-no-stabilize` on a tracked
  sequence disables it whatever `--export-window` says. Under `--fixed-tripod`
  the window is combined on the sky instead and survives.
- **An untracked window needs stars.** On `--fixed-tripod` the window is
  aligned on the sky, which requires a usable star fit in the reference; if
  there are too few, the run warns that `--export-window will fall back to 1`.
- **The `--fixed-*` flags require `--fixed-tripod`.** `--fixed-search`,
  `--fixed-anomaly` and `--fixed-no-stabilize` are read only in untracked mode
  and are ignored, silently, without it.
- **`--stack-max-foreground` has two defaults**: `0` for a tracked sequence,
  `1` under `--fixed-tripod`, where the landscape is part of the picture. And
  it has a floor: when fewer than a fifth of the frames come in under the
  limit, frames with foreground are admitted with a mask up to `0.6` anyway,
  and the run says so.
- **`--limit` cuts the file listing first.** Files are sorted by name, the
  first one becomes the reference that defines the output geometry, and
  `--limit N` keeps the first N — so it always includes the reference.

## Reading the console output

The run narrates itself; every line is there to be checked rather than watched.

```
found 517 frames
reference: _DSC0001.ARW (6048 x 4024) in 0.28s
stars in reference: 40
  [1/517] reference: foreground 18.4%
pipeline: 12 workers, up to 90 frames in flight (total RAM 62.4 GiB, budget 80% = 50.0 GiB, 557 MiB/frame)
  [10/517] _DSC0010.ARW: 40 stars, 31 inliers, θ=-0.019°, t=(1.43,6.20), sky 0.2908, foreground 16.9% in 1.61s
  ...
aligned: 459  interpolated: 10  skipped: 48  time: 66.9s
stack selection: 153 of 469 aligned frames (median sky 0.2296, tolerance ×1.60); 293 with foreground (masked); ...
glare/gradient [stack]: ...
stacking 153 frames...
stacked: 153  total time: 165.1s
crop from drift: 6048×4024 → 5800×3848 (−176 px height, −248 px width)
DNG written: stacked.dng
exporting 451 frames to res_clean (window 7 frames, stabilize=true, deflicker=true, 5800×3848); 66 skipped (unaligned or sky saturated/above the white)
  [436/451] _DSC0436_clean.dng  window 7 frames  transients 37246 px  defects ×1.10/1.16/1.31  anomaly 19.5%  gain R/G/B 0.831/0.815/0.824  sky +10.61%
  [450/451] _DSC0450_clean.dng  window 5 frames  transients 0 px  defects ×1.00/1.00/1.00  anomaly 114.0%  gain R/G/B 0.930/0.932/0.676  sky +113.60%  dawn 93%
exported 451 frames in 414.2s
```

An untracked run replaces the star fit with the drift it measured, and names
the frame it takes the export's levels from:

```
tripod drift: measured against the reference's landscape (18.4% of the frame), up to 64 px
  [10/517] _DSC0010.ARW: 38 stars, drift (+1.43,-0.62) px, corr 0.81, sky 0.2908 in 1.61s
  [11/517] _DSC0011.ARW: 37 stars, drift not measurable, identity, sky 0.2911 in 1.58s
levels reference: _DSC0233.ARW (sky 0.2301, closest to the session median 0.2296) — the star-trail stack is not one
```

What to look at:

- **`inliers`** per frame. A healthy fit sits well above the minimum of six. A
  frame that suddenly drops is usually cloud.
- **`θ` and `t`** should evolve *smoothly*. They do here: the drift creeps from
  `t=(0.87,0.78)` to `t=(5.00,12.70)` over twenty frames. A jump means a hot
  pixel or a satellite trail was matched as a star — the plausibility check
  catches it and reports `interpolated`.
- **`stack selection`** tells you why frames were dropped, split by cause
  (foreground, bright sky, dark sky) and how many were saturated enough not to
  be worth exporting at all.
- **`defects ×`** on export is the per-channel gain applied to the defect model
  for that frame; it climbs towards dawn because stray light scales with the
  light entering the lens. **`anomaly`** is how much that frame differs from the
  session median, as a percentage of the sky level. **`dawn`** appears once the
  frame is being ramped towards its natural version.
- **`transients … px`** is how many pixels were rescued from the temporal
  combination — a meteor, a satellite, a plane.
- **`drift` and `corr`**, on an untracked run, are the equivalent of `θ`/`t`
  and `inliers`: the drift should stay small and evolve smoothly over the night
  — legs settling, wind — and the correlation should sit well clear of the
  `0.25` floor. A page of `drift not measurable, identity` means the landscape
  is not being found, which on a genuinely tracked sequence means the mode is
  wrong; see [Picking the wrong one](#picking-the-wrong-one).

## What it can achieve

The correction was built one stage at a time, and every stage was measured
against the one before it on real output DNGs: the previous stage on the left,
the new one on the right, **both stretched identically** (percentile 0.5–99.8)
so the comparison is honest rather than flattering. What each stage bought, and
what one of them cost, is the table below.

| # | Stage | What it does, and what the before/after proves |
|---|---|---|
| 01 | Parametric model | No correction → parametric model. `poly3(x, y)` for the light-pollution gradient, a free radial profile about a *searched* optical axis for the halo, pupil ring and vignetting, plus `cos2θ`/`cos4θ` harmonics for non-circular vignetting from a hood or filter holder. The bowl and the halo go; the Milky Way stays, because it is rejected as an outlier by the robust fit. |
| 02 | Residual edge surface | Adds the lower-envelope residual surface at the edges: horizontal bands and flare wedges that no radial model can describe. |
| 03 | Fine surface, 64 px nodes | Surface nodes tightened from 192 px to 64 px, solved matrix-free by conjugate gradient. The thin lines and wedges that survived at the edges disappear. |
| 04 | Row bands and flare spokes | Fixed 1D patterns whose geometry belongs to the hardware: sensor row banding (per-row median, high-passed in `y`) and flare spokes (per-angle median about the optical centre, high-passed in `θ`), both subtracted over the whole frame. |
| 05 | Audit of stage 04 | **A failure, kept on purpose.** The audit of stage 04 shows the edge surface eating a real, extensive dark patch — a dark nebula, top centre. This is exactly what the audit exists to catch. |
| 06 | Band-pass edge surface | The fix: the edge surface becomes band-pass (fine − Gaussian σ=150 px), so only fine structure is removed and extensive sky patches survive. The dark patch is gone from the removed layer. |
| 07 | Mirror-symmetric coarse component | From the coarse component, only the part **shared with its own left/right mirror** about the optical axis is subtracted. System defects are symmetric about the lens; the sky is not. Wide edge wedges removed, dark nebula intact. |
| 08 | Horizon glow | Horizon glow: a robust 1D profile along the direction of maximum residual gradient, searched over the full 0–360° so nothing is assumed about where the horizon sits in the framing. On export, every frame additionally gets its own `poly3` gradient for varying light pollution and twilight, pedestal preserved. |
| 09 | Stable timelapse | The sequence as a sequence: same similarity and same crop as the stack, a 7-frame sliding trimmed mean, per-channel deflickering against the stack. Frames 1, 61 and 121 come out identical in tone, brightness and framing. |
| 10 | Transients preserved | Transients survive the trimmed mean: `|frame − combination| > 4σ` over elongated connected components of ≥ 40 px, restored verbatim from the original frame. `_DSC0335` and `_DSC0406` keep their trails. |

The arc of that table is the point. Stages 01–04 remove more and more; stage 05
proves one of them went too far; stages 06–07 make the correction *narrower* on
purpose, trading a little defect removal for the certainty that no sky is being
taken. That trade is the whole design.

### A real session, measured

517 frames of 6048×4024 from a dusk-to-dawn night, on 12 worker threads with
62.4 GiB of RAM:

| Stage | Result |
|---|---|
| Alignment (pass 1) | 459 aligned, 10 recovered by interpolation, 48 skipped — 66.9 s |
| Stack selection | 153 clear-sky frames kept of 469; 293 dropped for foreground, 13 for a bright sky, 10 saturated beyond exporting |
| Stacking (pass 2) | 153 frames — 165.1 s total |
| Crop from tracking drift | 6048×4024 → 5800×3848 (−248 × −176 px) |
| Timelapse export | 451 frames — 414.2 s |

Two numbers are worth pulling out. **153 frames stack but 451 export**: rejection
from the average is not rejection from the sequence, so most of the night that a
stacker would throw away still becomes video. And a whole night of tracking
drift costs only 8 % of the frame area, because the crop is computed as the
largest rectangle every aligned frame actually covered rather than guessed at.

## How it works

### 1. Reference and star detection

The reference frame is demosaiced with a bilinear kernel (all four Bayer phases
are supported) and normalized from black to white level into `[0, 1]`. Its
luminance feeds a star detector that estimates the background by median and MAD,
keeps local maxima above `median + 5σ`, refines each position to sub-pixel
accuracy with a parabolic fit, and ranks the candidates by integrated flux. Only
the brightest `--max-stars` survive. A reference with fewer than six stars aborts
the run, because everything that follows depends on it.

### 2. Pass 1 — alignment, background maps and foreground masks

Every remaining frame is loaded in parallel and matched against the reference.
Matching is scale-invariant by construction: triangles built from the detected
stars vote for star correspondences through their normalized side ratios, RANSAC
then proposes a similarity (rotation plus translation) from pairs of matches, and
the surviving inliers are refined by least squares. A fit needs at least six
inliers to be accepted.

Because tracking drift is smooth, an accepted fit is still cross-examined against
the interpolation of its neighbours: if it displaces any sensor corner by more
than six pixels relative to that prediction — the signature of hot pixels or a
satellite trail being mistaken for stars — it is discarded and replaced by the
interpolated transform, and the frame is dropped from the stack. Frames that
never produced a usable fit are recovered the same way, which keeps them
available for the timelapse export even though they no longer contribute to the
stack.

The same pass also produces, for every frame, a low-resolution background map
(a clipped median over 32×32 px blocks, in white-balanced scale) and a foreground
mask derived from it. The mask is what separates trees, mountains and buildings
from sky: a smooth `poly3 + radial profiles` model is fitted to the frame's own
map, and cells falling below 55 % of that prediction are marked as foreground,
then filtered for tiny components and dilated by one cell. An opaque object
against the sky is nearly black; a dark nebula or real vignetting never drops
that far below a fit that already includes the radial profile.

Under `--fixed-tripod` this pass changes shape: the star fit is replaced by a
two-scale normalized cross-correlation of the high-passed luminance restricted
to the landscape cells (coarse search over the whole drift range, then a fine
pass with parabolic sub-pixel refinement), and the per-frame masks are collapsed
into a single consensus mask. Star fits between *consecutive* frames are still
computed when the export needs them, because the export's temporal window has to
be combined on the sky, not on the sensor.

### 3. Frame selection

With sky level and foreground fraction known for every frame, the stack membership
is decided. Frames are rejected when their sky level departs from the session
median by more than `--stack-sky-tolerance` — twilight, moonrise or passing
clouds — and when their foreground fraction exceeds `--stack-max-foreground`.
The default of `0.0` stacks only clear-sky frames, since frames showing the
horizon usually carry the horizon glow with them, and that glow varies over time
and leaves bands once averaged. If that rule leaves fewer than 20 % of the aligned
frames, the threshold is automatically relaxed to 60 % and the offending pixels
are simply masked out instead. Under `--fixed-tripod` the default becomes `1.0`:
the landscape is in every frame and is part of the picture, so it is never a
reason to drop one.

Rejection from the stack is not rejection from the export: excluded frames are
still written by `--export-clean`, the only exception being frames whose sky is
already saturated, where there is nothing left to clean.

### 4. The defect model

The maps of the selected frames are combined into a temporal median that ignores
each frame's foreground, so a tree occluding part of one frame does not bias the
cell where other frames see sky. Cells that were occluded in nearly every frame
are filled in from their neighbours and marked as having no real data, so they
influence no fit.

That median map is then decomposed in four successive stages, each one operating
on the residual left by the previous:

1. **Parametric model.** A single robust least-squares fit (IRLS) of
   `poly3(x, y) + Σ f_p(r)·cos(h_p θ)`. The cubic polynomial absorbs the
   non-linear light-pollution gradient; the free radial profile, centred on an
   optical axis that is searched for rather than assumed, captures the halo, the
   pupil ring and the vignetting bowl; the `cos2θ` and `cos4θ` harmonics account
   for non-circular vignetting from the lens hood or filter holder. Structure
   that does not fit this description — the Milky Way, nebulae — is rejected as
   an outlier, which is exactly what preserves it.
2. **Residual surface.** A smooth bilinear surface on 64 px nodes, solved by
   matrix-free conjugate gradient, is fitted to the residual as a *lower
   envelope*: asymmetric IRLS discards whatever sits clearly above it and follows
   what sits below. Only the surface's fine component is subtracted, plus the
   part of its coarse component that is shared with its own left/right mirror
   about the optical axis. System defects are symmetric; the sky is not. The
   surface is ramped in between `0.5·r_full` and `0.9·r_full`, so the centre is
   left to the parametric model alone.
3. **Fixed 1D patterns.** Two patterns whose geometry belongs to the hardware
   rather than to the sky: sensor row banding (a per-row median, high-passed in
   `y`) and flare spokes (a per-angle median around the optical centre,
   high-passed in `θ`). Both are subtracted over the whole frame.
4. **Horizon glow.** A robust 1D profile along the direction of maximum residual
   gradient, searched over the full 0–360° so that no assumption is made about
   where the horizon lies in the framing.

Throughout, the pedestal — the median sky level — is deliberately preserved.
Only spatial variation is removed, never the absolute level.

### 5. Pass 2 — stacking

Each selected frame is reloaded, white-balanced, and cleaned in sensor
coordinates with the global model plus a per-frame *temporal anomaly*. That
anomaly is the key to handling a sequence that changes over hours: everything
constant during the session cancels in `frame map − session median map`, so what
remains is what actually varies — horizon glow, twilight, halo amplitude
following sky brightness, faint clouds. A smooth surface is fitted to it and
subtracted alongside the model, together with a per-channel defect gain
constrained to `[1, 3]`, since stray light scales with the light entering the
lens and can never fall below the deep-night median. How much freedom that
surface gets is fixed at `Full` for a tracked sequence and chosen by
`--fixed-anomaly` for an untracked one.

Weather is not an anomaly of the instrument, and the surface's own amplitude is
what tells the two apart: over the test session a clear frame's surface spans a
few per cent of the sky and a clouded one's spans one to four times it. Past a
quarter of the sky the surface is therefore faded out, and by a full sky nothing
of it is left — the cloud stays in the picture as bright as it was. Subtracting
it instead would flatten the cloud to sky level and leave its stars standing on
a background that is no longer there, which is what gave a clouded frame a
harder, whiter star field than a clear one. The dome amplitude is bounded on the
same grounds (see `dome_bounds`), and neither guard needs a cloud detector.

Optionally, `--scatter-comp` compensates the contrast loss the veil implies:
where the defect is strong a fraction of the image light has been scattered into
it, so after subtracting the veil the deviation from the sky level is rescaled by
`1/(1 − β·veil/sky)`, capped at ×2.

Cleaned frames are then resampled into the reference system and accumulated
bilinearly, with foreground pixels contributing geometric coverage but no sample
and their edges ramped to avoid steps. Once every frame is in, the result is
cropped to the largest rectangle that every aligned frame covered, which removes
the borders left ragged by tracking drift.

### 6. Output

The stack is analysed per channel to derive a linear stretch: white at the 99.9th
percentile, black at `median − 3·MAD`. Computing both ends per channel is
deliberate — a global white point would compress the dominant channel against a
foreign one and introduce a tint. Those points are baked directly into the 16-bit
data rather than delegated to the developer, so the clipping is visible even in
raw processors that ignore `BlackLevel` on LinearRaw DNGs. `--no-stretch` keeps
the sensor-native scale instead.

The file is written as a linear RGB DNG (`PhotometricInterpretation = LinearRaw`,
DNG 1.4.0.0) with the demosaic step already done, `WhiteLevel = 65535`,
`BlackLevel = 0`, a D65 `ColorMatrix1` looked up per camera model (Sony A7 III,
A7 IV, A7R III/IIIA, A7R IV/IVA and A6400, with an identity fallback), and
`AsShotNeutral = [1, 1, 1]` because white balance is already baked in — which
prevents the developer from applying it a second time.

The photographic metadata of the source RAW is then written in natively: make
and model in IFD0, and lens make and model, focal length (and its 35 mm
equivalent), aperture, exposure time, ISO and capture time in a real Exif IFD,
which is where exiv2 — and therefore darktable and lensfun — look for them. That
is exactly the set a lens profile is matched on, so distortion, vignetting and
chromatic aberration corrections light up on a bare DNG with nothing installed.
The DNG is not rewritten to do it: TIFF places no constraint on where an IFD
lives, so the new IFD0 and the Exif IFD are appended to the end of the file and
the header is repointed at them, leaving the image strips and every value the
original IFD0 referenced exactly where they were.

Finally `exiftool`, if it is on `PATH`, copies the rest of the source metadata
across — MakerNotes, XMP, ICC — excluding the structural and DNG tags the writer
controls itself. Without it the DNG keeps the native EXIF and only loses those.

### 7. Timelapse export

`--export-clean` reuses everything above to write one cleaned DNG per frame, in
chronological order and with the same stretch as the stack, so the sequence and
the stack are directly comparable. Each frame is cleaned with the model and its
own anomaly, stabilized into the reference system with the very transform used
for stacking and cropped identically, then levelled: per channel, its sky level
and its star-core amplitude are matched to the stack's by a linear gain and
offset, which stops tone, brightness and contrast from drifting across the
session.

Neither of those two statistics is read the obvious way, because a cloud must
not move them. The sky level is the 20th percentile of the lightly blurred
channel rather than its median: a cloud only ever adds light, so the sky is read
from below it, and a frame the cloud covers for the most part still measures its
clear sky and not its cloud. The star-core amplitude is the 99.9th percentile of
the **high-pass** — the channel minus its own blur — rather than of the channel
itself: a cloud is smooth and contributes nothing to it, whereas a
light-polluted cloud is brighter than any star and used to take the plain
percentile over completely. The gain itself may compress the
contrast but never expand it: a frame whose star cores come out below the
reference's has almost always lost them to the atmosphere — cloud, haze, dew —
and stretching them back up is precisely what made the stars of a clouded frame
burn brighter than the stars of a clear one, and what pushed the cloud they sit
on past the white of the stretch on the way.

On an untracked sequence the levels reference is **not** the stack, and neither
is the stretch. There the stack is a star-trail image: a star lands on a given
pixel in only a few frames, so the mean divides its flux by the length of the
session and the stack's high percentile no longer stands for a star core.
Matching every exported frame to that would crush its contrast, frame after
frame. Both the reference and the stretch are taken instead from one
representative frame — the one whose sky level sits closest to the session
median — cleaned exactly as the export cleans it. The run names the frame it
picked.

Temporal noise is then reduced by a trimmed mean over a sliding window of
`--export-window` frames, discarding the per-pixel maximum and minimum — a choice
that reduces noise by roughly `√(N−2)` while rejecting trails and hot pixels for
free. Foreground pixels are excluded from that combination, so trees and horizon
stay sharp instead of being smeared with their neighbours. On an untracked
sequence the window is warped onto the sky through the chain of consecutive star
fits before combining, since there the sky is what moves between neighbours.

A cloud is the one thing in the frame that does not follow the sky: it crosses
it. Combining the window on the sky therefore hands back, at full brightness, a
star this frame has behind a cloud and its neighbours have in the clear — over a
cloud that the same combination has meanwhile averaged smooth. Before anything
else is restored, the export undoes its own combination wherever that happened.
What it measures is the loss of star signal, not brightness: the high-pass with
its noise floor subtracted, averaged over a neighbourhood, in the frame and in
the combination. Where the frame keeps as much as the combination the two agree;
where a cloud has taken the stars away, the frame's is only its noise and the
ratio collapses. The frame is faded back in over that ratio, with the weight
smoothed over 128 px so that it follows the body of a cloud rather than stamping
a disc around every dimmed star. Brightness is deliberately not part of the
test — the Milky Way is broad and bright too, and it is full of stars.
`--export-no-cloud-guard` turns both halves of this off.

Because a trimmed mean would also erase what appears in a single frame, genuine
transients are detected and restored afterwards. The frame's excess over the
combination is high-passed, thresholded at 4σ, and filtered to keep only
connected components spanning at least 40 px — a long exposure trail is always
elongated, whereas a star residue is not — after which those pixels are taken
verbatim from the original frame. Meteors, satellites and planes therefore
survive intact.

Dawn is handled as a continuous transition rather than a cutoff. As the sky level
climbs past 1.5× the session median, the export ramps from the night-sky treatment
towards a natural version carrying only the defect model, reaching it at 2.3×; the
sequence ends when a fully natural frame would clip to white anyway.

## Auditing the correction

The most important property of this tool is negative: **it must not remove sky.**
`--dump-correction` writes exactly what was subtracted from the stack, plus the
preserved pedestal, under the same crop and stretch as the output, so the two can
be compared side by side in a developer.

A healthy removed layer shows **defect geometry only** — radial lines, horizontal
bands, wedges towards the edges, the vignetting bowl. If a nebula-shaped patch
appears in it, a stage is taking sky, and stage 05 in the table above is what
that looks like when it happens.

The full procedure is the three runs spelled out under
[Auditing an output](#auditing-an-output); everything is done in linear scale,
with `--no-stretch`, so the numbers mean something.

Four things are worth putting side by side: the uncorrected stack, the corrected
stack, the total removed layer, and the **fine** removed layer, that is
`layer − layer_param`, exactly what the surface, bands and spokes took, read as
a percentage of the sky level with zero at mid grey. That last one is what
decides whether a change was an improvement.

Safety knobs when a scene defeats the model: `--no-residual-surface` restricts
the correction to halo and gradient, `--no-flatten` disables it entirely, and the
constants `SURF_R_START` / `SURF_R_END` / `SURF_COARSE_SIGMA_PX` in
`src/flatten.rs` control where and how gently the surface is allowed to act.

## Diagnostics

Three environment variables expose the internals without changing the output:

| Variable | Effect |
|---|---|
| `APILAAA_DEBUG` | Prints the per-channel radial profiles and lists every frame excluded from the stack with its sky level and foreground fraction. |
| `APILAAA_DEBUG_DIR=<dir>` | Dumps the measured background map and the evaluated model as CSV, the temporal medians of each chronological half with their mean transforms, and the transient masks as downsampled PGM. |
| `APILAAA_DEBUG_CORR` | Prints, per frame, coarse numeric maps of the anomaly, the model correction, the fitted surface and the resulting residual, as a percentage of the sky level. |

The two halves dumped by `APILAAA_DEBUG_DIR` are worth knowing about: if the
defect model is genuinely stationary, the first half of the night and the second
half must produce the same map. Where they differ is where the anomaly, not the
model, should be doing the work.

## Performance and memory

Both passes are parallel over frames, and the inner loops over rows are
parallelized with Rayon. Before starting, the tool sizes the pipeline against
the machine: it budgets 80 % of total RAM, subtracts the resident accumulator
(20 B per pixel) and divides the remainder by the cost of a frame in flight
(24 B per pixel) to decide how many workers to spawn and how deep the channel
between producers and consumer may be. The chosen figures are printed at startup,
so an unexpectedly slow run is easy to attribute.

Rough expectations from the measured session above (24 MP frames, 12 threads):
about 0.13 s per frame for pass 1, 0.6 s per stacked frame for pass 2, and
0.9 s per exported frame. Cutting `--export-window` to `1` removes most of the
export cost at the price of the noise reduction.

## Repository layout

| Path | Contents |
|---|---|
| `src/main.rs` | CLI, two-pass orchestration, frame selection, alignment plausibility checks, RAM planning. |
| `src/raw.rs` | RAW decoding, black/white normalization, bilinear demosaic, luminance. |
| `src/stars.rs` | Background estimation and sub-pixel star detection. |
| `src/align.rs` | Triangle matching, RANSAC and least-squares similarity fitting. |
| `src/fixed.rs` | Untracked sequences: tripod-drift measurement by cross-correlation, sky transform between consecutive frames. |
| `src/flatten.rs` | Background maps, foreground masks, consensus mask, the four-stage defect model, per-frame anomaly, correction grid. |
| `src/stack.rs` | Weighted accumulator, coverage tracking and valid-rectangle cropping. |
| `src/output.rs` | Stretch analysis, linear DNG writing, metadata copying through exiftool. |
| `src/exif.rs` | Native EXIF: reads the source RAW's identifying and photographic tags and writes them into the DNG by relocating its IFD0. |
| `src/timelapse.rs` | Clean sequence export: stabilization, deflickering, temporal denoising, transient preservation, dawn ramp. |
| `vendor/rawloader/` | `rawloader` 0.37.2 with `data/cameras/sony/a6400.toml` added, patched in from `Cargo.toml`; the upstream release refuses A6400 files outright. |

## Limits

- **Sony `.ARW` only.** The input listing accepts the `.ARW` extension and
  nothing else. Decoding goes through `rawloader`, which reads far more formats,
  so widening the filter is a one-line change — but the colour matrices are
  written for Sony bodies and another brand would stack with an identity matrix.
- **One night per run.** The defect model is fitted on the session's own median,
  so mixing sessions — or lenses, or focal lengths — breaks the assumption that
  the defects are stationary.
- **Similarity alignment only** (rotation + translation + implicit scale). Field
  rotation from an alt-az mount, or lens distortion across a very wide field, is
  not modelled.
- **Frames must be chronological by filename.** In-camera numbering is; renamed
  files may not be.
- **A framing that is all sky** gives `--fixed-tripod` nothing to correlate
  against; below 2 % landscape it declines to estimate the drift and leaves
  every frame at identity — which for a tripod that genuinely did not move is
  the right answer anyway.
- **The export window of an untracked sequence rests on a similarity.** The
  sky's apparent motion is a rotation on the celestial sphere; a similarity only
  approximates it, and the wider the field or the longer the interval, the worse
  the approximation. The run reports the median residual of the fit between
  consecutive frames and warns past a pixel; if the stars come out soft, lower
  `--export-window`.

## Licence

Copyright (C) 2026 Salvatore Josue Trimarchi Pinto.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.

The full text is in [`LICENSE`](LICENSE).
