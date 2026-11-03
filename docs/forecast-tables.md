# Forecast tables: wavesim as a per-spot transfer function for gogo

A scope, written 2026-11-03. Nothing here is built. It says what a table would hold, what it
would cost to compute, what wavesim has to gain first, which spots have a seabed good enough to
run it on, how gogo would use it, and how gogo's session ratings would become the first
real-world check of wavesim.

## In short

- **The table.** For one spot it maps the offshore swell (height, period, direction) and the
  tide to the wave at the spot: how tall it is when it breaks, where and in what depth, its surf
  similarity (spilling or plunging), and its peel angle. Four axes, about 10 × 6 × 6 × 3 points;
  gogo would read it the way it already reads `rays.json`, by interpolating in plain Python at
  score time.
- **The bathymetry and gogo's spots don't overlap well.** The good beds (NOAA 3 m, USACE 1 m)
  are in Hawaii, and gogo's 16 spots are all on the Lisbon–Ericeira–Peniche coast. For that
  coast, the 2011 DGT/APA coastal LiDAR, published by EMODnet at about 15 m, covers the surf
  zone of every gogo spot. It reaches about 0.5 to 1 km out and 5 to 12 m deep. The Coxos run
  (`out/coxos_2m13s_300`) used it. That is good enough for reefs and points, but not for sand,
  whose bars have moved since 2011.
- **Compute.** About 25 minutes per run on 3 efficiency-core threads. A 40-run pilot at Coxos
  takes about two nights on the laptop. A full table for one spot, cut down to the swells that
  actually happen there, is roughly 350 runs: about 19 nights, or about half a day on a rented
  16-vCPU machine for a few euros.
- **wavesim needs a few things first.** With today's side walls, one bed cannot serve every
  direction, so a direction needs its own bed (a re-crop). It also needs a mode that keeps only
  statistics, because one run's frames are 150 to 440 MB. It needs an in-run harmonic
  accumulator, so the peel angle can be computed without frames, a batch driver, wall time in
  `run.json`, and a defined lineup per spot. Spectral input is a separate research step.
  Spreading the forecast across the table's own axes gets most of the benefit for free.
- **gogo's gate.** gogo is in Stage 3, but its Stage 1 gate (about 100 real labels) is not met,
  and on 2026-10-14 it decided not to do break-scale physics until labels show a tide × size
  interaction that the low/mid/high classes miss. So wavesim can build the table now. In gogo it
  can at most be a candidate that serving never sets, and only once Max reopens that decision.
  It ships only on a backtest number from real same-day pairs.

## 1. The table

### Inputs and axes

The input is what gogo already stores for each hour (`HourForecast`): `swell_height_m`,
`swell_period_s` (the mean period; `swell_peak_period_s` where it is forecast), `swell_from_deg`,
and `sea_level_m`, read 32 minutes back (`tide.SEA_LEVEL_LEAD`). Like `rays.py`, the table
treats the forecast cell as open sea. Section 3 covers how that open-sea swell gets to
wavesim's wave-maker, which sits 1 to 3 km out in 23 to 35 m of water. Which period field feeds
the table should be a candidate option tested by the backtest, not a choice made here: the two
fields are a median 1.5 s apart on this coast (gogo S6).

The axes for Coxos, the first spot (gogo's spec: swell from 270–340°, 1.2–3.5 m, at least 10 s,
low and mid tide):

| Axis | Range | Step | Points | Why this resolution |
|---|---|---|---|---|
| Direction (open sea) | 260–350° | 10° | 10 | The spot's window and 10° past each edge. Direction matters most: on Pipeline, turning a 3.5 m, 16 s swell 20° moved the median break from 5.4 to 8.6 m deep and the rows plunging from 88 to 1 of 200. 10° is the coarsest step worth trying, and the pilot checks it at the 5° midpoints. |
| Period | 8–18 s | 2 s | 6 | Breaking height goes about as T^0.4 (Komar & Gaughan), which is smooth. Starts 2 s under the spot's minimum, so the veto edge is inside the table. |
| Height | 0.75–4.5 m | ×1.43 | 6: 0.75, 1.1, 1.6, 2.25, 3.2, 4.5 | Geometric, because the break changes in kind at the top: on Pipeline a 5 m swell moved to the outer reef 450–570 m further out. |
| Tide (sea level) | −1.2 to +1.2 m | 1.2 m | 3 | At Coxos, going from mid tide to −1.2 m moved the median break depth from 5.1 to 4.3 m and the plunging rows from 54 to 42 of 200 (the two existing runs), so the tide moves the break but not abruptly. At spring extremes beyond ±1.2 m the table holds its edge value and flags it. |

That is 1,080 points in full. gogo's `rays.json` uses 5° × 1 s, and wavesim cannot afford that;
section 2 shows how the full count is cut down to the swells that actually happen.

### Outputs at each point

A run covers 1.2 km of coast, but a spot is a few hundred metres of it. Each spot needs a
**lineup**: the rows (a `y` range) and the stretch of `x` where people take off, drawn once in
the bed's frame. Every output is summarised over those rows only.

| Output | From | Unit |
|---|---|---|
| Break height: median and 90th percentile over the lineup rows | `break_line[].height_m`: crest to trough, the tallest wave on the row | m |
| Where it breaks: distance offshore of the 0 m (MSL) contour, median; and the lat/lon of the median break point, for display | `break_line[].x_m` and the bed | m |
| Depth at the break, still water at that tide | `break_line[].depth_m` | m |
| Surf similarity at the break: median, 10th–90th percentile, and the share of rows that plunge | `break_line[].surf_similarity` (Battjes' ξ_b; wavesim's thresholds are 0.4 and 2.0) | – |
| Breaker type | From the interpolated ξ_b. Interpolate the number, never the label. | – |
| Peel angle: median and middle half; the share of rows under 10° (closing out) | `tools/direction.py` | ° |
| Wave direction at the break, relative to `+x` | `tools/direction.py` | ° |
| How uneven the crest is: the largest change in break height over 60 m along the shore | `tools/direction.py` | m |
| Quality flags | `run.json`: `unsettled_per_period` > 3%, `truncated`, `coarse_breaking_in_zone`, lineup rows breaking deeper than 10 m | – |
| Provenance | wavesim commit, the bed file's hash, the run options, the platform, the wall time | – |

The peel angle is the angle between the crest and the line of break points, as
`tools/direction.py` defines it: 0 when the wave breaks all at once along that line, larger the
more it peels. In the surfing-science literature (Walker 1974; Hutt, Black & Mead 2001), smaller
peel angles mean faster, harder waves, and near 0 a close-out.

### The file, and reading it each hour

Follow `rays.json`: a package-data file, `src/gogo/data/breaks.json`, with one entry per
tabulated spot (axes, outputs as arrays, provenance), read by a small module, `gogo/breaks.py`,
in plain Python. At about 15 numbers per point and 1,080 points, one spot is 100–200 KB of JSON.

- **Interpolation.** Multilinear in all four axes (the 16 corners around the hour). Direction
  wraps at north as in `rays.factor`. Interpolate the numbers (height, depth, ξ_b, peel) and
  derive the breaker type and the close-out flag afterwards.
- **Outside the table.** Hold the edge value and set a flag, never extrapolate. An hour outside
  a spot's table falls back to the incumbent score for that spot and says so in the reasons. A
  missing table must not become a veto, because a wrong veto leaves no trace (gogo's
  censoring trap).
- **Cost.** One lookup per tabulated spot per hour, which is microseconds. No change to the
  worker.

**I'd compute it from `forecast_snapshots`, not store it in them.** The brief put the table into
`forecast_snapshots` as a feature. gogo's own rules argue against writing it there. Snapshots
record what Open-Meteo said, deduplicated by payload. The table's output is a deterministic
function of fields already stored plus a table version, so it can be recomputed for any past
as-of, and by gogo's rule of asymmetric recoverability ("store the forecast fields now, score
them later") it can wait. Writing it into the payload would append a new row for every
(cell, hour) each time the table changes, which would break the dedup. It would tie the history
to one table version. And snapshots are per grid cell while the table is per spot. If a stored
copy is ever needed for speed, it belongs in its own derived table keyed on
`(spot, valid_at, fetched_at, table_version)`.

## 2. Compute

### What runs cost today

`run.json` does not record wall time. The times below come from the run directories' file
times: the bed is written at the start and the frames at the end. They match the README where it
quotes a figure (11.7, 25 and 57 minutes). All of them ran on 3 threads on the efficiency cores.

| Run | Swell | Fine grid (3 m cells) | Simulated | Wall time |
|---|---|---|---|---|
| `pipeline_swell_1.5m16s_h9` | 1.5 m, 16 s | 189 × 400 | 344 s | 10.2 min |
| `pipeline_new` | 2.5 m, 14 s | 179 × 400 | 328 s | 11.7 min |
| `coxos_2m13s_300_low` | 2 m, 13 s, tide −1.2 m | 438 × 400 | 282 s | 21.1 min |
| `coxos_2m13s_300` | 2 m, 13 s, mid tide | 423 × 400 | 271 s | 22.7 min |
| `pipeline_b110_3.5m16s` | 3.5 m, 16 s from 290° | 315 × 400 | 388 s | 24.6 min |
| `pipeline_swell_5m16s` | 5 m, 16 s | 494 × 400 | 401 s | 56.6 min |

Bigger and longer swells cost more three ways: the wave-maker moves deeper, so the domain
lengthens; the run lasts longer, because the swell travels further and its period is longer;
and taller waves shrink the time step. At Coxos, plan on 12–15 minutes for the small corners of
the table, 20–25 minutes for the middle, and 50–70 minutes for 4.5 m at 16–18 s. **A mean of
about 25 minutes per run** is the planning figure below.

### How many runs

| Plan | Points | Runs | Laptop: one run at a time on the efficiency cores, about 19 runs per 8-hour night | Rented 16 vCPU, about 30 runs an hour |
|---|---|---|---|---|
| Pilot, Coxos | 3 directions (290, 310, 330°) × 3 periods (10, 13, 16 s) × 2 heights (1.5, 3 m) × 2 tides (±1.2 m), plus 4 midpoints | 40 | 2 nights | 1–2 h |
| Coxos, full factorial | 10 × 6 × 6 × 3 | 1,080 | 57 nights | 1.5 days |
| **Coxos, cut to its climatology** | the points this coast's swell actually reaches, plus a ring around them | **~350** | **~19 nights** | **~12 h** |
| Five reef spots, each cut to its climatology | | ~1,750 | ~3 months | ~2.5 days |

**Cutting to the climatology.** gogo holds a year of ERA5 swell at Coxos's cell (`gogo
backfill`). Counting those hours per table point shows which points never occur, and that is
most of them: on this coast height and period correlate at 0.62 (swells under 1 m average 6.8 s,
over 2.5 m average 10.3 s; gogo S13a), so 4.5 m at 8 s and 0.75 m at 18 s are empty, and the
swell at the Ericeira cell is concentrated in the north-west (25% of it from 320–335°; S13c).
Keep every point that has hours, plus one point either side on every axis, so an interpolation
never reads a missing corner. The "about a third" in the table is a guess until those hours are
counted, and counting them is the first next step.

### Where to run them

- **The laptop, overnight, on the efficiency cores.** This is the setup the laptop compute rule
  allows: background priority, 3 threads, one run at a time, quiet. Running two at once would
  double the rate, but it would occupy every efficiency core, so it needs Max's OK first. A pilot
  and a two-spot table fit in a month of nights, which is slow but costs nothing.
- **A rented machine.** Hetzner bills its dedicated-vCPU cloud servers by the hour, and a
  16-vCPU one costs well under a euro an hour (check the current list). wavesim's threading
  scales as about threads^0.8 and tops out around 6, so throughput comes from running several
  independent runs at once, each on 4 threads. The README puts a performance core at 3.5 times
  an efficiency core. If a server thread is about as fast, a run takes about 8 minutes and four
  at once give about 30 runs an hour. That makes one spot's table half a day and a few euros.
  First do one calibration run there: `estimate_seconds` is calibrated on an M5, and the budget
  guard would mis-estimate on another CPU.
- **Not the gogo box.** It is a 2-vCPU product host, and its job is the hourly worker.

**Two practical points.** Disk: one 2D run writes 150–440 MB of frames (`out/` is 5.1 GB
today), so 350 runs would be about 90 GB. A table run must keep statistics only (section 3).
Reproducibility: results are byte-identical across thread counts, but not necessarily across
CPU architectures, because `sin` and `exp` come from the platform's libm. Record the platform in
the provenance, and run a few points on both machines to bound the difference.

## 3. What wavesim needs first

### Direction: can one bed serve every direction?

**Not with today's side walls.** The wave-maker supports oblique incidence. `Drive` in
[`crates/wavecore/src/forcing.rs`](../crates/wavecore/src/forcing.rs) sets
`k_y = k sin θ`, drives the phase `ωt − k_y y`, and divides the flux by `cos θ`. The Snell and
wave-action test checks it, in shallow water with periodic `y`. Three things stop it serving a
real bed:

1. `wavesim run` never uses it: it builds the maker with `angle: 0.0`
   ([`crates/wavesim/src/lib.rs:1241`](../crates/wavesim/src/lib.rs)).
2. The sides are walls, and a wall mirrors. A swell at θ reflects off each side wall as one at
   −θ, so the domain fills with the crossing pair, a pattern that stands along the shore. The
   wall is exact only for θ = 0 over a bed that is uniform along the shore, which is why it was
   chosen.
3. Periodic `y` lets an oblique wave through, but the bed must then be periodic too, and it must
   satisfy `k_y · L_y = 2πn`. On a 1.2 km box that quantises the direction to steps of about 9°
   for a 13 s swell in 25 m of water. A real reef is not periodic, and joining its two sides
   makes a seam that waves would refract on.

The options:

- **A. One bed per direction (recommended for now).** `tools/fetch_spot.py` crops and rotates a
  bed to any bearing, from a URL or a local GeoTIFF. That takes under a minute and about a
  megabyte. The 290° Pipeline run did exactly this. The cost is geometry: as the swell turns away
  from the shore normal, the shoreline crosses the box at an angle, so the box needs more water
  offshore and more land inland, and the walls mirror a shore that is not square to them. At 20°
  on Pipeline, a box twice as wide agreed within 9% (median). Beyond 20° this is untested, so
  repeat that double-width check at each spot's two extreme directions.
- **B. One bed, an oblique maker, and new side boundaries.** For example, relaxation strips along
  the sides, offshore of the surf zone, pulling towards the incoming oblique wave. This is real
  numerics work with known-answer tests to write first. One such test: an oblique swell over an
  alongshore-uniform beach must give rows that are exact phase-shifted copies of each other. The
  README already records that oblique waves with dispersion are untested, and that relaxing the
  fine grid's sides towards the coarse run failed (13–22% off near the break). B is worth doing
  only together with directional spreading (below), which needs it anyway.

**Which direction the table's axis means.** The open sea is tens of kilometres beyond the box,
and the shelf in between refracts the swell before it reaches the wave-maker. gogo's ray tracer
(`scripts/wave_rays.py`) already carries swell from open sea to 10 m off each spot. Traced
instead to the wave-maker's depth (23–35 m), and made to report the mean local direction there
as well as the height factor, it would map gogo's open-sea direction onto each run's bearing. At
Ericeira's open spots that correction is small: the rays give Coxos 0.82–1.04 of the open-sea
height across its window at 10–16 s. At sheltered spots, where the Peniche peninsula or Cabo
Raso stands tens of kilometres off, it is the larger part of the effect. This change belongs in
gogo's scripts, not in wavesim.

### Irregular and spectral input

The input is one monochromatic swell, and the README lists that as a limitation. There are three
levels of doing better:

1. **Monochromatic, with a stated convention (v1).** Run a regular wave of height Hs at the
   chosen period, and say so in the table's metadata. The error, mainly that every wave breaks
   in the same place, is part of what the labels then calibrate.
2. **Spread at read time (v1, free).** Weight the table across its own axes when reading it:
   over height by a Rayleigh distribution of individual waves around Hs, over direction by the
   cos^2s spread gogo's rays use (s = 20), over period by the ±12% the rays use. This gives a
   share of waves breaking here rather than further out, for example "the sets break on the outer
   reef", with no new runs. What it misses is everything waves do to each other: groups, surf
   beat, the change in setup, and a broken set wave changing the one behind it.
3. **True irregular input (a later research step).** The source is linear in its strength, so
   the maker could sum many `Drive`s, one per frequency (and direction) component, each with its
   own group speed and Gaussian roll-off, at random phases. That needs:
   - longer runs, 100–200 waves instead of the 4 a run records now, so about 6–10 times the
     cost, and several random seeds;
   - new break statistics: Hs from the variance in each cell, where Hs peaks, and the share of
     time each cell spends breaking (`stats.f32` already holds `break_fraction`);
   - attention to the long waves bound to the groups, for the same reason the maker sits where
     the second harmonic is under 10%;
   - option B above for directional spread, because with walls every oblique component is
     mirrored;
   - tests with known answers first, for example a JONSWAP sea keeping its spectrum over a flat
     bed, and Battjes–Janssen's height decay on a plane beach.

Recommendation: v1 is levels 1 and 2. Level 3 is worth it only if the labels show that v1's
errors depend on the sea state (for example, wrong on long-period groundswell, right on
windswell).

### Other things found

- **Peel needs frames today.** `tools/direction.py` takes the phase at the swell's frequency from
  the stored frames over the last periods. An accumulator in the run (the complex amplitude at ω
  in every cell over the statistics window, two more fields in `stats.f32`) would give
  direction and peel from the statistics alone. A `--no-frames` mode could then keep a table
  run at a few MB.
- **Wall time is not in `run.json`.** Record the estimate and the actual time. Budgeting a batch
  needs them, and they show when `estimate_seconds` drifts.
- **A batch driver.** Something like `wavesim table <spot> --plan plan.toml`. It would be
  resumable (skip a point whose inputs hash to a finished run), run one point at a time on the
  efficiency cores, keep no frames, and write one summary row per run. A reducer would then
  turn those rows into `breaks.json` with the provenance attached.
- **The lineup per spot.** The rows and `x` range, in `spots.toml`, beside the bed.
- **Coxos is not reproducible yet.** It is not in `spots.toml`. Its bed was cropped from
  `data/coxos/coxos_seabed_msl.tif`, a merge of the EMODnet HR LiDAR (~15 m) and the EMODnet DTM
  (~115 m), shifted to mean sea level. That merge is not a tool in `tools/` yet, so a re-crop
  per direction depends on it becoming one.
- **The datum.** The Coxos bed's datum is "lowest astronomical tide + 2.0 m, approximate". A
  0.3 m error there moves the tide axis by a quarter of a step. Check it against the
  Hidrográfico's mean sea level at Cascais and Peniche (gogo's `scripts/tide_check.py` already
  reads that source).
- **Some Coxos rows break deep.** Of the 200 trusted rows, 6 (y = 301.5–310.5 and 394.5–397.5 m)
  break in 10–12 m of water, 470–500 m along `x`, far outside the rest (median 5.1 m). That may
  be the outer reef, or something at the edge of the trusted width. Look at it before a lineup
  is drawn.
- **The bed is smoother than its cells suggest.** Resampled from ~15 m to 3 m, the extra
  resolution adds nothing: no ledge exists in the data. Separately, on 3 m cells the break
  height is about 9% low against converged 1.5 m cells (README, Limitations).
- **The bed may not reach deep enough.** The Coxos bed is 40.7 m deep at its deepest, and the
  4.5 m, 18 s corner may want the maker deeper than that (5 m at 16 s needed 42 m at Pipeline).
  Either fetch further offshore (the 115 m DTM is fine for a wave-maker) or cap the table.
- **Wind is not modelled.** gogo's wind term stays. Offshore wind holding up the face, which
  pushes towards plunging, is out of scope.

## 4. Which spots have bathymetry

| Area | Source | Cells | What it covers | gogo spots |
|---|---|---|---|---|
| Oahu North Shore | NOAA CUDEM (`spots.toml`: `pipeline`) | 3 m | Measured topo-bathy lidar over the surf zone | none |
| Oahu | 2013 USACE NCMP lidar DEM | 1 m | No gaps in 1–10 m along the lines checked | none |
| Lisbon–Peniche | EMODnet HR `590_HR_Lidar_Centro` (2011 DGT/APA coastal LiDAR) | ~15 m (0.00013°) | 38.50–40.43°N, 9.53–8.82°W; along each spot's facing it reaches 0.5–1 km out and 5–12 m deep | all 16, in the surf zone |
| Lisbon–Peniche | EMODnet DTM 2024 | ~115 m | The shelf, out to where the wave-maker goes | all 16, offshore only |
| Ericeira, Peniche, Cascais–Oeiras | DGT 2 m MDT, coastal zones, 2011 (600 m at sea, 400 m on land) | 2 m | By request only; a request is drafted in `data/dgt_pedido_mdt_2011.md` | would be reef-scale where it reaches |

Hawaii has the best data and none of gogo's spots. gogo is Lisbon–Ericeira–Peniche only, by its
frozen decisions. Pipeline stays useful for developing and verifying the method, not as a
forecast. [bathymetry.md](bathymetry.md) predates the EMODnet HR find, and its "no reef-scale
open data for Portugal" still holds: 15 m is not reef-scale. But 15 m does resolve the 100 m
scale over which a breaker's slope is measured, which the 115 m DTM cannot. In this file the
LiDAR's interpolation flag is unset (fill value) everywhere, so measured and interpolated cells
cannot be told apart.

gogo's spots, by bottom type. The type column is my reading, not something in `coast.yml`, so
Max should check it:

| Spot | Bottom (to confirm) | Table worth building? |
|---|---|---|
| Coxos | reef point | **yes**, first; already run |
| Ribeira d'Ilhas | reef point | **yes**, second; 1.5 km south of Coxos, so the two can share a source tile and a pilot's lessons |
| São Lourenço, Pedra Branca, Lagide, São Pedro do Estoril | reef | yes, after the pilot (Lagide's coordinates sit 0.79 km off the coastline; fix them first, gogo S13b) |
| Furnas, Empa, Consolação | reef or mixed? | check the bottom first |
| São Sebastião, Foz do Lizandro, Guincho, Carcavelos, Caparica, Supertubos, Baleal | sand | no: bars move with every swell and the survey is from 2011. For these the rays already carry the shelf, and nothing can see the bars (gogo S13f) |

**The Coxos run, briefly.** A 2 m, 13 s swell from 300°, with the bed cut along 120°. Over the
200 trusted rows:

- It breaks in 5.1 m of water (median; middle half 4.4–5.6 m), or 4.3 m at −1.2 m tide.
- It breaks at 4.2 m crest to trough (middle half 3.9–4.5 m).
- ξ_b is 0.31 (median), and 54 rows plunge (42 at low tide).
- Unlike at Pipeline, the reef turns the swell: at the break it travels +10° from `+x`
  (median), and 109 rows are more than 10° off.
- It peels at 28° (median; middle half 15–52°).

Two things to carry forward. The break height is 1.4 times Komar & Gaughan's 3.0 m for that
swell, while gogo's rays give no shelf focusing (0.97). So either the last 500 m focuses it, or
the tallest-wave convention and the slightly-too-tall inner surf zone (H/h near 1.0) inflate it.
That is the first number the labels can check, through the `size` faults. And a spot locals know
for heavy, hollow waves comes out mostly spilling. On Pipeline's 3 m bed the cause was a reef
smoother than the real one, and a 15 m bed is smoother still.

## 5. How gogo would use it

### Where gogo is, and what that allows

From gogo's README (Status, 2026-10-27) and `gogo-notes.md`:

- **Stage 2, the harness, is done.** That is as-of features, the event-grouped dataset, five
  metrics against five baselines, and `gogo backtest`. So "no score improvements until the
  harness exists" is satisfied.
- **Stage 3's rule then applies:** a change lands only if the backtest improves, or is neutral
  for a written reason. Every candidate so far (`size_period:typical`, `dir_taper`, `shadow`,
  `rays`, `rays:size`) is built and waiting, because **the Stage 1 gate, about 100 real labels
  across several swell events, is not met.** Without real pairs, every Δ reads n/a.
- **On 2026-10-14 gogo decided "not doing: break-scale physics".** That covers break depth with
  tide and breaker type, because no public data had the seabed at the break ("EMODnet at ~100
  m"). The decision is to revisit only with labels that show a tide × size interaction the
  low/mid/high classes miss. The 15 m LiDAR weakens the premise. It does not meet the condition.

So:

- **Now:** build the table in wavesim. That is wavesim's own work and touches nothing in gogo.
- **When Max reopens the 2026-10-14 decision** (the trigger it names, or the label gate being
  met): add the table to gogo as candidates only, following `rays`. Run `candidates.movement`,
  which needs no labels, to show how big the bet is and which days a `checked` pair at two
  spots would settle.
- **Never without a number:** making it the default and bumping `SCORE_VERSION`. That needs a
  paired Δ whose interval is clear of zero under `best_known`, with the lead policies no worse
  (gogo S13a).

### How it plugs in

- **Data and reader.** `src/gogo/data/breaks.json` and `gogo/breaks.py`, read at score time
  from the hour's own fields, as `rays.py` is. No migration and no worker change. Q1
  (`best_known`) feeds it ERA5 swell and the archive's sea level (from late 2022). Q2 feeds it
  forecast snapshots at the as-of. The table is the same in both, and only its inputs differ.
- **Separate candidates for separate claims**, so the backtest can tell which part helps:
  - `breaks:size`: the size gate reads the table's break height instead of the offshore height.
    The spot ranges were written in offshore height, so re-anchor them the way
    `size_period:typical` did: a typical day scores exactly as it does now.
  - `breaks:tide`: replace `tides: [low, mid]` with what the table says the tide does to this
    swell, such as break depth and ξ_b. This is the tide × size interaction the 2026-10-14
    decision asks labels to show.
  - `breaks:peel`: a close-out penalty where the peel angle is near 0 across the lineup.
- **Reasons stay one sentence.** For example: "1.6 m at 14 s from 305° breaks at 3.4 m in 4.5 m
  of water at mid tide, plunging, peeling about 30°."
- **Per spot.** Only tabulated spots change, which fits `pinned:<version>`: change one spot,
  hold the rest, and ask whether the change helped. A table version (the digest of its
  `breaks.json` entry) plays the role `spec_version` plays for specs, and impressions would
  carry it if it ever ships.
- **Sample size.** The headline metric counts same-day pairs of spots, and ties drop about half
  of them (gogo S10). Only pairs involving a tabulated spot can move, so a Coxos–Ribeira table
  needs labels at those two spots specifically, including `checked` rows, across several swells
  and tides.

## 6. Validation: the first real-world check of wavesim

So far wavesim is verified against exact solutions and laboratory beaches. Its only real-world
comparison, Pipeline, is against secondary sources that disagree with each other by 250 m. gogo's
labels would be the first real test, and an indirect one: a label is a blind 1–5 rating, a
would-return, optional fault codes (`size:±1`, `tide:±1`, …) and a time interval. It is not a
measured wave. What the labels can test is the **ordering** wavesim implies (which days, which
tides, which swell directions are better at a spot), not its metres.

1. **The ranking test.** Is pairwise accuracy on same-day pairs better with a `breaks`
   candidate than with the incumbent? gogo's harness answers this as it stands.
2. **Faults against the physics.**
   - `size` faults at Coxos under the incumbent should line up with where the table's break
     height differs most from the offshore height, and should fall under `breaks:size`. A
     per-spot scale factor fitted later (S18) says how far off wavesim's height is. One that
     holds steady across swells is a bias that can be corrected. One that moves with period or
     direction is a physics error in wavesim or in the bed.
   - `tide` faults should cluster where the table says the tide shuts the spot, for that swell
     size. A cluster there is the evidence the 2026-10-14 decision asks for. A cluster elsewhere
     says the model's tide response is wrong.
3. **Ratings against ξ_b and peel.** At a reef rated for advanced surfers, ratings should rise
   with ξ_b and with a peel neither near 0 (close-out) nor near 90°. A flat relation means the
   table's breaker type carries no information there. Given how smooth the bed is (section 4),
   that is a live possibility.
4. **Write the expectations down before looking**, as both projects do with their test bounds.
   For example: the `breaks:size` paired Δ is positive on pairs involving Coxos or Ribeira;
   `size` faults at Coxos fall by at least a third; the plunging share at Coxos is higher at low
   tide in the labels' best sessions.
5. **Optional, and a gogo product decision.** For tabulated spots only, the card could ask two
   absolute, blind questions: where it broke (inside, main peak, outside) and whether it closed
   out. They are observations, not comparisons, so they do not anchor the rating. They would
   lengthen a card whose gate is 15 seconds, though, and `kind=cam` sessions could answer the
   first instead.

Plan on at least 20–30 real labels at tabulated spots, across five or more swell events and
across tides, before reading anything. Ties halve the usable pairs. The labels cannot test
absolute height, tube size, or anything finer than the bed, which is 15 m in Portugal.

## Next steps

1. **Throughout, and the gate for everything on gogo's side:** keep logging real sessions at
   Coxos and Ribeira, including `checked` rows on the same morning.
2. Count the climatology: hours per table point from gogo's backfilled year at the Coxos and
   Ribeira cells (a read-only query). Fix the table's points and the run count from it.
3. Make Coxos reproducible: a `[coxos]` table in `spots.toml`, the LiDAR + DTM merge as a tool,
   and the datum checked against the Hidrográfico's mean sea level.
4. wavesim changes, each with a test: wall time in `run.json`, the in-run first-harmonic
   statistics so direction and peel need no frames, `--no-frames`, and the resumable
   `wavesim table` batch driver with its reducer.
5. Draw the Coxos lineup in the bed's frame, and decide on the six rows that break deep.
6. The pilot, 40 runs at Coxos, about two nights on the efficiency cores. Before running it, set
   bounds for interpolating to the midpoints, for example break depth within 0.5 m, break
   height within 10%, ξ_b within 0.05. Then run the double-width wall check at the window's
   extreme directions.
7. Choose the venue for the full tables: laptop nights or a rented machine, after one
   calibration run there and a cross-platform comparison of two points.
8. The full, climatology-cut tables for Coxos and Ribeira, reduced to `breaks.json` with
   provenance.
9. In gogo, and only after Max reopens the 2026-10-14 decision: rays traced to the wave-maker's
   depth, with the local direction; `breaks.py`; the three candidates; and the movement report.
10. Validation: the backtest on real pairs, against expectations written down before looking
    (section 6), and only then any question of shipping.
