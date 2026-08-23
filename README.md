# apilaaa

A RAW astrophotography stacker that aligns a night-sky sequence on its stars,
models and removes the optical defects of the camera system, and writes the
result as a linear DNG that any raw developer can pick up. The same model is
reused to export a flicker-free, noise-reduced timelapse sequence from the very
same exposures.

The guiding idea is that a tracked sky sequence contains two superimposed
signals with opposite behaviour: the sky drifts across the sensor between
frames, whereas the lens halo, the vignetting, the flare wedges and the sensor
banding stay fixed to the sensor. Because those defects are additive and
stationary in sensor coordinates, they can be estimated from the temporal
median of the sequence and subtracted before averaging — which is precisely the
order of operations the tool follows.

## Requirements

- A Rust toolchain (edition 2021).
- [`exiftool`](https://exiftool.org/) on `PATH`. It is not needed to produce the
  image data, only to copy EXIF, MakerNotes, XMP and ICC from the source RAW
  into the output DNG; if it is missing, the DNG is still valid and the run only
  warns.
- Linux is the tested platform. The RAM planner reads `/proc/meminfo` and falls
  back to a conservative 4 GiB assumption elsewhere.
- Optional, for the audit scripts in `comparaciones/`: Python 3 with NumPy and
  Pillow.

## Build

```sh
cargo build --release
```

The release profile is tuned for throughput (`opt-level = 3`, thin LTO, a single
codegen unit), which matters because the whole pipeline is CPU-bound and
processes every frame twice.

## Usage

Point the tool at a directory of Sony `.ARW` files and give it an output path:

```sh
./target/release/apilaaa --input res --output stacked.dng
```

Files are sorted by name, which for in-camera numbering is chronological order,
and the first one becomes the reference frame that defines the output geometry.
Everything downstream — alignment, cropping, the timelapse export — is expressed
in that reference's coordinate system.

To also obtain the cleaned timelapse sequence from the same run:

```sh
./target/release/apilaaa -i res -o stacked.dng --export-clean res_clean
```

And to audit what the correction actually removed, ask for the removed layer as
its own DNG:

```sh
./target/release/apilaaa -i res -o stacked.dng --dump-correction layer.dng
```

### Options

| Option | Default | Purpose |
|---|---|---|
| `-i, --input <DIR>` | `res` | Directory holding the input `.ARW` files. |
| `-o, --output <DNG>` | `stacked.dng` | Output DNG for the stack. |
| `--max-stars <N>` | `40` | Stars used per frame for alignment. |
| `--limit <N>` | — | Process at most N frames; useful for quick tests. |
| `--no-stretch` | off | Keep the native sensor scale instead of baking a stretch. |
| `--no-flatten` | off | Disable the whole defect model; stack the frames as they are. |
| `--no-residual-surface` | off | Keep only the parametric halo + gradient model, skipping the residual surface, bands and spokes. |
| `--dump-correction <DNG>` | — | Write the removed layer, with the same crop and stretch as the output. |
| `--export-clean <DIR>` | — | Export the stabilized, cleaned timelapse sequence. Requires the correction to be enabled. |
| `--export-window <N>` | `7` | Odd size of the temporal window used to reduce noise on export; `1` disables it. |
| `--export-no-stabilize` | off | Export in sensor coordinates, uncropped (this also disables temporal denoising). |
| `--export-no-deflicker` | off | Skip the per-channel level matching on export. |
| `--export-no-transients` | off | Let the temporal combination erase meteors, satellites and planes. |
| `--scatter-comp <BETA>` | `1.0` | Contrast compensation in the veiled areas; `0` turns it off. |
| `--stack-sky-tolerance <F>` | `1.6` | Frames whose sky level leaves `[median/F, median·F]` are excluded from the stack. |
| `--stack-max-foreground <FRAC>` | `0.0` | Frames with more than this fraction of foreground are excluded from the stack. |

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

### 3. Frame selection

With sky level and foreground fraction known for every frame, the stack membership
is decided. Frames are rejected when their sky level departs from the session
median by more than `--stack-sky-tolerance` — twilight, moonrise or passing
clouds — and when their foreground fraction exceeds `--stack-max-foreground`.
The default of `0.0` stacks only clear-sky frames, since frames showing the
horizon usually carry the horizon glow with them, and that glow varies over time
and leaves bands once averaged. If that rule leaves fewer than 20 % of the aligned
frames, the threshold is automatically relaxed to 60 % and the offending pixels
are simply masked out instead.

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
lens and can never fall below the deep-night median.

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
prevents the developer from applying it a second time. Finally, `exiftool` copies
the source metadata across, excluding the structural and DNG tags the writer
controls itself.

### 7. Timelapse export

`--export-clean` reuses everything above to write one cleaned DNG per frame, in
chronological order and with the same stretch as the stack, so the sequence and
the stack are directly comparable. Each frame is cleaned with the model and its
own anomaly, stabilized into the reference system with the very transform used
for stacking and cropped identically, then levelled: per channel, its sky median
and its 99.9th percentile are matched to the stack's by a linear gain and offset,
which stops tone, brightness and contrast from drifting across the session.

Temporal noise is then reduced by a trimmed mean over a sliding window of
`--export-window` frames, discarding the per-pixel maximum and minimum — a choice
that reduces noise by roughly `√(N−2)` while rejecting trails and hot pixels for
free. Foreground pixels are excluded from that combination, so trees and horizon
stay sharp instead of being smeared with their neighbours.

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

The most important property of this tool is negative: it must not remove sky.
`--dump-correction` writes exactly what was subtracted from the stack, plus the
preserved pedestal, under the same crop and stretch as the output, so the two can
be compared side by side in a developer. The removed layer should show defect
geometry only — radial lines, horizontal bands, wedges towards the edges. If a
nebula-shaped patch appears in it, a stage is taking sky.

`comparaciones/` documents that verification step by step, with a before/after
image per stage and the scripts that generate them:

```sh
cd comparaciones
python3 compare.py before.dng after.dng out.png 0.5 99.8 4 "label A" "label B"
python3 audit.py layer.dng layer_param.dng uncorrected.dng corrected.dng out.png label
```

Two safety knobs exist for when a scene defeats the model: `--no-residual-surface`
restricts the correction to halo and gradient, and `--no-flatten` disables it
entirely.

## Diagnostics

Three environment variables expose the internals without changing the output:

| Variable | Effect |
|---|---|
| `APILAAA_DEBUG` | Prints the per-channel radial profiles and lists every frame excluded from the stack with its sky level and foreground fraction. |
| `APILAAA_DEBUG_DIR=<dir>` | Dumps the measured background map and the evaluated model as CSV, the temporal medians of each chronological half with their mean transforms, and the transient masks as downsampled PGM. |
| `APILAAA_DEBUG_CORR` | Prints, per frame, coarse numeric maps of the anomaly, the model correction, the fitted surface and the resulting residual, as a percentage of the sky level. |

## Performance and memory

Both passes are parallel over frames, and the inner loops over rows are
parallelized with Rayon. Before starting, the tool sizes the pipeline against
the machine: it budgets 80 % of total RAM, subtracts the resident accumulator
(20 B per pixel) and divides the remainder by the cost of a frame in flight
(24 B per pixel) to decide how many workers to spawn and how deep the channel
between producers and consumer may be. The chosen figures are printed at startup,
so an unexpectedly slow run is easy to attribute.

## Repository layout

| Path | Contents |
|---|---|
| `src/main.rs` | CLI, two-pass orchestration, frame selection, alignment plausibility checks, RAM planning. |
| `src/raw.rs` | RAW decoding, black/white normalization, bilinear demosaic, luminance. |
| `src/stars.rs` | Background estimation and sub-pixel star detection. |
| `src/align.rs` | Triangle matching, RANSAC and least-squares similarity fitting. |
| `src/flatten.rs` | Background maps, foreground masks, the four-stage defect model, per-frame anomaly, correction grid. |
| `src/stack.rs` | Weighted accumulator, coverage tracking and valid-rectangle cropping. |
| `src/output.rs` | Stretch analysis, linear DNG writing, metadata copying. |
| `src/timelapse.rs` | Clean sequence export: stabilization, deflickering, temporal denoising, transient preservation. |
| `comparaciones/` | Stage-by-stage before/after evidence and the audit scripts. |
