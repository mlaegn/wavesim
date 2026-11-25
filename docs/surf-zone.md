# The surf zone: broken waves that stay too tall

A plan, written 2026-11-21. Step 1, the switch, is done (PR 5). Step 2 ran on 2026-11-23 and did
not pass on the 2D runs' 3 m cells: the fix leaves the break point where it was but hardly
changes the surf zone. On 2 m and 1.5 m cells it works, `H/h` about 1.0 down to 0.76–0.79
with the break unmoved: the cause is the grid (results under step 2). But Coxos in 2D on 2 m
cells, with and without the fix, shows the fix breaking the reef early: breaking starts 66 m
further out, in 6 m of water instead of 3.8, and the surf zone comes out at `H/h` 0.48, too
low. So the fix stays off; steps 3 and 4 do not go ahead as planned. The fallback tried next,
plain shallow water shoreward of a known break line (2026-11-24), passes the laboratory beach but
fails on the reef lines: the sheet's break line is the tallest wave, which on gentle reefs is
well seaward of where the wave breaks. Correcting the sheet's heights is what is left. It follows the 2026-11-17 audit, which named this the model's main accuracy gap
once its bugs were fixed.

## In short

- **The problem.** wavesim puts the break point in the right place, but once a wave has broken
  it hardly loses height on its way to the shore. In the inner surf zone it stays about as tall
  as the water is deep (`H/h` near 1.0). Regular laboratory waves sit at 0.58 to 0.61 there
  (Ting & Kirby's spilling breaker, which the model reproduces 1.66 times too tall), and
  random waves in the field at 0.5 to 0.68 by the Surfline camera fields.
- **What it spoils.** Heights shoreward of the break, which the spot sheets show as their height
  map, and how much the tide changes the size there (0.78 to 0.81 on the Coxos runs, against
  0.50 to 0.68 in the field). Not the break point, its depth or its height, and so not the
  takeoffs or the sheets' breaking zone either: since 2026-11-17 that zone runs from each row's
  break point to the shoreline, whatever the broken wave does in between. (The sheet preview's
  legend still calls it "breaks on most waves", the name from before; gogo's says "breaking and
  white water, to the shore".)
- **A fix exists.** `Breaking::viscosity` follows each breaking front from step to step and
  gives it Kennedy et al.'s (2000) eddy viscosity, the closure FUNWAVE uses (`wavecore::eddy`).
  On Ting & Kirby's beach it brings the inner surf zone to 1.04 times the laboratory's, no
  station off by more than 11%. It is off by default because it breaks the plunging laboratory
  wave 0.86 m early, against 0.5 m allowed. It has never run on a reef; `wavesim run` has a
  switch for it since step 1.
- **The plan.** Add the switch; try it on single cross-shore lines of Pipeline and Coxos, which
  is light work; only if the break points hold there, run Coxos in 2D; then decide whether it
  becomes the default. Each step has bounds written down before its first run, as the lab tests
  had.
- **The risk to watch.** Reef breaks plunge, and the plunging laboratory wave is exactly where
  the fix fails. If the reefs break early with it too, the fix stays off and the fallbacks below
  apply.

## Step 1: a switch on `wavesim run` (done, PR 5)

- `--breaking-viscosity` on `wavesim run`, off by default, building the solver with
  `Breaking::viscosity` set to `EddyViscosity::default()`. Refused without dispersion
  (`--dispersive false`), since plain shallow water has no breaking fronts to follow.
- `run.json` records `"breaking_viscosity"` (true or false) for both grids, and the `model`
  label says when it is on, so a sheet or a viewer can tell which model made a run.
- Refused together with `--handoff-y`, before anything is written: the slice brings its swell in
  without the viscosity, so a wave handed over from such a run would not match.
- Tests: the flag reaches the solver, and only when set; a run with it records it and stays
  finite; both refusals. What it does to breaking waves is the laboratory test's job, in the slow
  set that runs on GitHub every Monday.
- Not done: the time estimate a run checks against `--max-minutes` does not count the
  viscosity's extra sub-steps. Step 3 measures them; until then a run with it may stop at its
  limit sooner than estimated, marked `truncated`.

## Step 2: single lines across real reefs

A run along one cross-shore line, four cells wide, with walls at the sides, is a 2D run's
machinery on a thousandth of the cells: the same planning, wave-maker, statistics and break
line, through `wavesim::run` on one grid. It is enough to see whether the break point moves,
which is the failure the lab showed. A line cannot spread energy sideways, so its heights are
not the 2D run's (on Pipeline y = 511.5 m a line's crest came out three times the 2D run's);
what counts is the difference the viscosity makes on the same line.

- **Lines.** Pipeline y = 361.5, 421.5 and 511.5 m (`data/pipeline_line_*.json`, 1 m points,
  put on 3 m cells), and three rows of the Coxos bed for the 300° swell (`data/coxos/coxos.json`,
  `+x` along 120°, 3 m cells): y = 481.5, 601.5 and 721.5 m, rows 160, 200 and 240, whose break
  points in the 2D runs are a plunging one and two spilling ones at mid tide.
- **Waves.** Pipeline 2.5 m, 14 s and 16 s, at mean sea level; Coxos 2 m, 13 s, at mid tide
  (0 m) and low tide (−1.2 m): the sizes the sheets and slices already use. Twelve cases.
- **Each case twice**, with and without the viscosity, dispersive, on one grid, the run
  planning its own length. Results in `out/surf_lines/`.
- **What is measured,** on a middle row, from the run's own statistics: the break point (the
  run's `break_line`: where, how deep, how tall) and the inner surf zone's `H/h`, the median over
  the cells from the break point to where the still water is 1 m deep, counting only cells no
  deeper than the break (a reef's lagoon behind the break is not surf zone). A case with fewer
  than five such cells has no inner surf zone to judge, which a shore break may not.
- **Bounds, fixed on 2026-11-23 before the first run:**
  - the break point moves no more than 10 m, and its depth by no more than 10%; for scale,
    halving the cells moved Pipeline's by 3.7 m;
  - the height at the break point changes by less than 10%;
  - the inner surf zone's `H/h` falls into 0.55 to 0.75, the laboratory's regular-wave value
    with room for the reef.
- **Passing.** Step 2 passes if the break point and its height stay within their bounds in every
  case, and the inner surf zone comes into its range in at least three quarters of the cases
  that have one. A case that misses is reported, not dropped. If the break moves in any case,
  step 3 waits and the fallbacks come first.
- **Where it goes.** `crates/wavesim/examples/surf_lines.rs` runs the cases and prints the
  table, so the numbers can be rerun; the table goes into this doc.

### Step 2 results (2026-11-23)

All 24 runs took 15 minutes on one efficiency core (wavesim's estimate was 43), the Mac at
nominal throughout. Break: distance along the line, still-water depth, height there. Inner surf
zone: median `H/h` and, in brackets, how many cells it is taken over.

| Case | Break, without | Break, with | Moved | Depth | Height | Inner surf `H/h`, without → with | Verdict |
|---|---|---|---|---|---|---|---|
| Pipeline 361.5, 14 s | 1342 m, 8.4 m, 4.62 m | the same | 0 | 0% | 0% | 0.63 (77) → 0.64 | pass |
| Pipeline 361.5, 16 s | 1474 m, 8.3 m, 4.63 m | the same, 4.64 m | 0 | 0% | 0% | 0.66 (75) → 0.67 | pass |
| Pipeline 421.5, 14 s | 1424 m, 6.8 m, 4.64 m | the same, 4.63 m | 0 | 0% | 0% | 0.80 (66) → 0.80 | surf zone too tall |
| Pipeline 421.5, 16 s | 1646 m, 4.4 m, 4.61 m | the same, 4.52 m | 0 | 0% | −2% | 0.95 (41) → 0.92 | surf zone too tall |
| Pipeline 511.5, 14 s | 1028 m, 10.8 m, 4.16 m | the same | 0 | 0% | 0% | 0.36 (149) → 0.35 | below the range (see below) |
| Pipeline 511.5, 16 s | 1654 m, 4.0 m, 4.42 m | the same, 4.31 m | 0 | 0% | −2% | 0.96 (38) → 0.96 | surf zone too tall |
| Coxos 481.5, mid | 940 m, 3.4 m, 3.40 m | the same, 3.38 m | 0 | 0% | −1% | 1.30 (8) → 1.30 | surf zone too tall |
| Coxos 481.5, low | 962 m, 3.6 m, 3.35 m | the same | 0 | 0% | 0% | 1.01 (16) → 1.01 | surf zone too tall |
| Coxos 601.5, mid | 898 m, 4.8 m, 4.61 m | the same, 4.59 m | 0 | 0% | 0% | 0.99 (10) → 0.99 | surf zone too tall |
| Coxos 601.5, low | 970 m, 3.9 m, 4.29 m | the same, 4.27 m | 0 | 0% | 0% | 1.19 (8) → 1.19 | surf zone too tall |
| Coxos 721.5, mid | 898 m, 1.8 m, 3.86 m | the same, 3.89 m | 0 | 0% | +1% | none: breaks onto the shore | pass |
| Coxos 721.5, low | 962 m, 4.9 m, 3.67 m | the same, 3.68 m | 0 | 0% | 0% | 0.79 (12) → 0.79 | surf zone too tall |

**What it shows.**

- **The break point does not move,** in any of the twelve cases: not a cell, and its height by
  2% at most. The failure the laboratory showed, a plunging wave breaking early, does not appear
  on these reef lines. That half of step 2 passes.
- **But the fix does almost nothing to the surf zone.** The median `H/h` changes by 0.03 at
  most. By the rule fixed beforehand step 2 does not pass: the surf zone is in range in 2 of the
  11 cases that have one, and those two (Pipeline 361.5) were already in range without the fix.
- **It is not switched off by mistake.** With it, cells shoreward of the break read as breaking
  more of the time in every case (Coxos 721.5 at low tide 3% → 9%, at mid tide 31% → 42%;
  Pipeline 421.5 at 16 s 1.9% → 3.4%), and wave heights there change by up to 0.3 m in places:
  fronts are followed and given their viscosity. They are just too rare and too brief to matter.
  On the laboratory beach the switch reads breaking across most of the inner surf zone; on these
  lines, shoreward of the break, it does so 1% to 30% of the time, and a front can only be
  followed once the switch has started it.

**How far to read it.**

- **Coxos's inside is short and reflective.** Between the break and 1 m of water there are only
  8 to 16 cells (24 to 48 m), and the 2D runs measured 14% to 56% of the swell coming back off
  that shore. Heights there include the returning wave, which no dissipation of the incoming
  bores removes, so Coxos's `H/h` above 1 says as much about reflection as about the surf zone.
- **Pipeline 511.5 at 14 s** breaks in 10.8 m of water and the wave re-forms over the deeper reef
  behind: its 0.36 is not a surf zone, and the filter meant to keep a lagoon out let it in.
- **Pipeline 421.5 and 511.5 at 16 s** are the clean cases: a long, shallow inside (38 to 41
  cells), no lagoon, and `H/h` 0.92 to 0.96 with or without the fix. That is the problem this
  plan is about, standing.
- **Lines, not 2D,** as said above: the comparison is between the two runs of each line.

**So,** as tuned on the laboratory beach, the fix is safe on these reefs but does not help them.
Step 3 would only confirm that in 2D, so it is not worth the runs. The next move is a choice:

1. **Start fronts more easily on a reef.** The fix follows a front only once the switch has
   started it. Starting fronts by Kennedy's own rate-of-rise criterion as well (tried on the
   laboratory beach: 1.06, nearly as good, and the plunging wave broke early that way too)
   might catch the reef's broken bores. It would need the laboratory tests again, then these
   lines again, with the break point watched as closely: light work.
2. **Plain shallow water shoreward of the break** (fallback 2): the laboratory's inner surf zone
   comes out at 0.65 to 0.70 that way. Code first, then these lines.
3. **Correct the sheet, leave the model** (fallback 3), and say so on the sheet.
4. **Park it.** The sheets' takeoffs, zones and break-point numbers do not depend on it; only
   their height maps shoreward of the break read high, as the sheets already say.

Either way, nothing here holds up other work: the fix will not become the default soon, so
sheets made now, Ribeira's included, will not need redoing because of it.

### Finer cells (2026-11-23): it is the grid

Max asked what is best for the project rather than to park it, so the next cheap question was
whether the 3 m cells are why the fronts are so rarely followed. The two clean cases, Pipeline
421.5 and 511.5 at 16 s, ran again on 1.5 m cells (`--cell 1.5`, the line data has 1 m points):
four runs, 18 minutes on one efficiency core, nominal throughout.

| Case, 1.5 m cells | Break, without | Break, with | Inner surf `H/h`, without → with | Breaking inside, without → with |
|---|---|---|---|---|
| Pipeline 421.5, 16 s | 1633 m, 4.4 m, 5.11 m | see below | 1.10, 1.08, 0.97 → 0.83, 0.76, 0.70 (at 3.7, 2.8, 2.1 m deep) | 4.2% → 7.3% |
| Pipeline 511.5, 16 s | 1639 m, 4.1 m, 4.80 m | 1636 m, 4.3 m, 4.78 m | 1.05 (76) → 0.79 (78) | 6.5% → 29% |

- **On 1.5 m cells the fix works.** The surf zone comes down from about 1.05 to 0.70–0.83, just
  above the laboratory's 0.6, and the cells inside read as breaking far more often: the fronts
  are followed. On 3 m cells the same lines barely changed.
- **The break does not move.** At 511.5 it stays within 3 m. At 421.5 the table's "tallest
  point" jumped 381 m out to sea, but that is the measure, not the wave: breaking starts at
  about 1625 m, in about 5 m of water, with and without the fix, and nowhere seaward of it in
  either run (breaking share 0 to 1600 m). The fix lowered the peak at the break from 4.66 to
  4.45 m, below a 4.76 m bump at 1250 m that the swell reflected off the reef makes in both
  runs, so the tallest point moved there. `break_points` takes the shoreward of two peaks only
  within 3% of each other. If the fix ever becomes the default, the break point should be
  where breaking starts, not the tallest point, or a reef that reflects will mislead it.
- **Why the grid matters.** A front is followed while the surface rises faster than a
  threshold tied to the wave speed. On 3 m cells a broken front 1 to 2 m tall is spread over one
  or two cells, so its rise reads slow and the front is dropped. On 1.5 m cells it reads fast
  enough. On the laboratory beach the cells are 5 cm, finer again. So option 1 above, starting
  fronts by the same rate of rise, would meet the same smearing and is not the next step.
- **Without the fix, finer cells make it worse:** `H/h` 1.05 against 0.95 on 3 m cells. The
  coarse grid was taking out a little height by itself.

### 2 m cells (2026-11-23): enough

The same two cases on 2 m cells (four runs, 10 minutes, nominal throughout). Measured from where
breaking starts, the first cell from the sea that reads as breaking at least 0.5% of the time,
since the tallest point proved unreliable; the same measure for all three cell sizes:

| Case | Cells | Breaking starts, without → with | Inner surf `H/h`, without → with | Breaking inside, without → with |
|---|---|---|---|---|
| Pipeline 421.5, 16 s | 3 m | 1648 → 1648 m, 4.3 m deep, height −3% | 0.92 → 0.89 | 1.9% → 3.5% |
| | 2 m | 1625 → 1621 m, 4.9 → 5.1 m deep, height −4% | 1.00 → **0.76** | 3.0% → 24% |
| | 1.5 m | 1625 → 1624 m, 4.8 m deep, height −3% | 1.05 → 0.77 | 4.0% → 27% |
| Pipeline 511.5, 16 s | 3 m | 1654 → 1654 m, 4.0 m deep, height −2% | 0.94 → 0.93 | 5.0% → 6.3% |
| | 2 m | 1641 → 1635 m, 4.3 → 4.6 m deep, height −8% | 1.01 → **0.79** | 5.3% → 28% |
| | 1.5 m | 1637 → 1636 m, 4.2 → 4.3 m deep, height −2% | 1.05 → 0.79 | 6.5% → 29% |

- **2 m cells do nearly all that 1.5 m cells do.** The surf zone comes down to 0.76 and 0.79,
  the break moves 4 and 6 m, its depth by 4% and 7%, its height by at most 8%: inside the
  bounds for the break, a few hundredths above the 0.75 top of the surf-zone range.
- **The step is between 3 and 2 m,** for breaks in 4 to 5 m of water: cells up to about half
  the breaking depth work, cells of three quarters of it do not. Where between those the line
  falls, and whether it scales with depth on other reefs, these two cases cannot say.
- **The table's tallest point again jumped at 421.5** (54 m out to a reflected bump), while
  breaking starts within 4 m of where it did. The measure above is the one to keep.

**What it means.** The fix is right on these reefs once the cells are 2 m or finer, about half
the breaking depth; the 2D runs use 3 m. That leaves:

1. **The DGT 2 m seabed.** Runs on its own 2 m cells would bring the surf zone in line with the
   fix on, and the finer seabed is what the sheets need most anyway, for whether a wave pitches
   and how it peels. It comes when the data comes.
2. **2 m cells on today's 15 m-based beds.** Resampling the bed to 2 m adds no detail but gives
   the fix the cells it needs. About 3.4 times the work of a 3 m run (2.25 times the cells, 1.5
   times the steps): a Coxos run of half an hour becomes about an hour and three quarters, a
   sheet of 14 runs about a day. Two or three nights on the laptop, or a rented machine.
3. **A rule that knows the grid,** so 3 m cells would do: new tuning, re-checked on the
   laboratory beach and these lines at several cell sizes. Light work, uncertain.
4. **Correct the sheet's surf-zone heights** by the measured ratio, about 0.75 to 0.8.

### Coxos in 2D on 2 m cells (2026-11-23/24): the fix breaks the reef early

Two runs to see whether option 2 holds in 2D: Coxos 2 m, 13 s from 300°, mid tide, on the bed
resampled to 2 m (`tools/fetch_spot.py coxos --cell 2`), without and with the viscosity
(`out/coxos_2m13s_300_2m_plain`, `out/coxos_2m13s_300_2m_visc`). 72 and 76 minutes on four
efficiency-core threads, nominal throughout. Compared over the sheet's trusted rows, 300–1050 m,
on the same bed coordinates, with the existing 3 m run without it (`out/coxos_2m13s_300`).

`H/h` and the share of time breaking, by water depth, shoreward of the 3 m run's break point:

| Water depth | `H/h`: 3 m / 2 m / 2 m with fix | Breaking: 3 m / 2 m / 2 m with fix |
|---|---|---|
| 1–2 m | 0.94 / 1.02 / 0.65 | 1.4% / 2.0% / 22% |
| 2–3 m | 0.89 / 0.96 / 0.52 | 0.1% / 0.7% / 26% |
| 3–4 m | 0.88 / 0.93 / 0.45 | 0.0% / 1.0% / 30% |
| 4–6 m | 0.73 / 0.77 / 0.45 | 0.0% / 0.2% / 35% |

Medians over the rows:

| Run | Tallest point | Breaking starts |
|---|---|---|
| 3 m | 4.75 m deep, 4.07 m tall | 1.44 m deep |
| 2 m | 4.60 m deep, 4.34 m tall | 3.76 m deep, 3.46 m tall |
| 2 m with fix | 6.98 m deep, 3.63 m tall | 6.07 m deep, 2.73 m tall |

- **The cells alone change little.** On 2 m cells the tallest point stays where it was (a metre
  away, 0.15 m shallower) and a little taller, the surf zone a little taller too (0.84 → 0.92
  from each row's break point), as on the lines. Breaking starts further out than on 3 m cells because the
  switch now trips on the reef at all; on 3 m cells it hardly does before the shore.
- **The fix breaks the reef early.** Against the 2 m run without it, row by row: breaking starts
  a median 66 m further out (32–179 m, 10–90%), in 70% deeper water, more than 10 m out in every
  row. It starts where the wave is 2.7 m in 6.1 m of water, `H/h` 0.45, well short of the
  height a wave breaks at. Each row's tallest point is 18% lower, and it stays within 10 m in only
  34% of rows.
- **So the surf zone is overdone.** `H/h` 0.48 from each row's break point, 0.45–0.52 in 2–6 m of water,
  under the 0.55 floor. Only the last metre or two (0.65) lands in range, because the waves
  arrive there already small.
- **The lines did not show it.** On 2 m cells only the two Pipeline lines ran, where the fix
  moved the start 1–6 m. Why Coxos in 2D differs is open. Two guesses worth testing: Coxos
  plunges, as the laboratory wave the fix breaks early does; or in 2D a front started at one
  focusing point is carried along the crest into deeper water.
- **Cost:** 72 and 76 minutes (the fine grid 60 and 63 of them), 2.4–2.5 times the 3 m run's
  half hour; the viscosity added about 5%.

**What it means.** The fix as it stands fails at Coxos in 2D: it is the risk the plan named,
the reef breaking early. Option 2 above falls with it. The fix stays off, and the fallbacks
below are what is left; the first of them, two passes, can now start from the 2 m run's own
break line.

## The fallback tried next: shallow water inside the break (planned 2026-11-24)

Of the three fallbacks below, plain shallow water inside the surf zone came first:

- **It should work on the 3 m cells the sheets already use.** A shock takes the height out of
  a broken wave by an amount the jump across it sets, not the cells, while the viscosity
  needed its fronts resolved, on 2 m cells and at 2.5 times the work.
- **Nothing to tune,** and on the laboratory beach plain shallow water already gives `H/h`
  0.65 to 0.70 inside.
- **It cannot move the break.** It acts only shoreward of a break line already known, so
  the early break that sank the viscosity at Coxos cannot happen; only what the surf zone
  sends back seaward (reflection, setup) can reach the break.
- **The sheets' runs without it exist,** and their break lines with them; on Coxos the 2 m and
  3 m lines agree.

Viscosity only shoreward of the break needs 2 m cells and did badly in the one form tried on the
laboratory beach; correcting the sheet stays the fallback if this fails.

**How it works.** `wavesim run --shallow-inside <run>` takes each row's break point from an
earlier run of the same bed, swell and tide without a fix (its `break_line`, the tallest wave
before the shore) and holds every cell from there to the shore to plain shallow water, on both
grids of a nested run (`Solver::with_shallow_cells`). There the breaking switch reads 0: no
dispersive terms, the MC limiter, a broken wave goes on as a shock. Seaward the run is the
dispersive one it was. The breaking statistics still read the surface, so `break_fraction`
compares with runs without it. It refuses the viscosity, plain shallow water and hand-offs.

**Checks, in order, each only if the one before passes. Bounds fixed on 2026-11-24, before the
first run:**

1. **The laboratory beach** (slow tests in `tests/lab.rs`, three flume runs each: one to set the
   wave, one without, one with):
   - the spilling breaker under the surf-zone test's own bounds: model over laboratory within
     10% seaward of x = 5.5 m and within 15% over the inner surf zone by the median, no
     inner-surf station off by more than 30%;
   - both breakers' break point within 0.2 m (four cells) of the run without, and so within the
     0.5 m of the laboratory's that the break-point tests allow.
2. **The reef lines on 3 m cells:** the twelve step 2 cases with `--fix shallow`, each run
   without it first for its break line, under step 2's bounds and pass rule unchanged: the break
   point within 10 m, its depth within 10% and height within 10% in every case; inner-surf `H/h`
   in 0.55 to 0.75 in at least three quarters of the cases that have one. Pipeline 511.5 at
   14 s, where the wave re-forms over the deeper reef, is reported apart: a re-formed wave
   inside the line is held to shallow water too, which steepens it early. About 15 minutes on
   one efficiency core.
3. **The two Pipeline 16 s lines on 2 m cells,** the same bounds, and inner-surf `H/h` within
   0.08 of the 3 m result, the change the cells alone make without a fix: this is what says
   the fix does not need fine cells. About 10 minutes on one core.
4. **Only then, Coxos in 2D on 3 m cells,** mid tide, from `out/coxos_2m13s_300`'s break line:
   about half an hour on four efficiency-core threads, with Max's go-ahead. Checked row by row as
   the 2 m runs were (break line, where breaking starts, `H/h` by depth band).

### Results (2026-11-24): the laboratory passes, the reef lines do not

**1. The laboratory beach passes.** Both slow tests, 7.7 minutes on the efficiency cores after a
5-minute build, nominal throughout:

| | Without a fix | Viscosity | Shallow water inside | Bound |
|---|---|---|---|---|
| Spilling, inner surf zone, model over lab (median) | 1.66 | 1.04 | **1.13** | within 15% |
| Spilling, worst inner-surf station | | 11% | **24%** | under 30% |
| Spilling, seaward of 5.5 m (median) | | | **1.04** | within 10% |
| Spilling, break point | 6.66 m | | **6.61 m** | within 0.2 m of without |
| Plunging, break point (lab 7.795 m) | 7.48 m | 0.86 m early | **7.48 m** | within 0.2 m of without |

Inside, the model's `H/d` comes out 0.65 to 0.79 against the laboratory's 0.58 to 0.61: a little
taller than with the viscosity, and the plunging wave keeps its break.

**2. The reef lines on 3 m cells do not.** 24 runs, 18 minutes on one efficiency core, nominal
throughout (`out/surf_lines/<case>/shallow`, `summary_shallow.json`):

| Case | Break, without → with | Inner surf `H/h`, without → with | Verdict |
|---|---|---|---|
| Pipeline 361.5, 14 s | 1342 → 1336 m, 8.4 → 8.3 m deep, height −19% | 0.63 → 0.37 | break moved |
| Pipeline 361.5, 16 s | 1474 → 1468 m, 8.3 m deep, height −15% | 0.66 → 0.42 | break moved |
| Pipeline 421.5, 14 s | 1424 → 1420 m, 6.8 → 7.0 m deep, −7% | 0.80 → 0.49 | surf zone out |
| Pipeline 421.5, 16 s | 1646 → 1630 m, 4.4 → 5.0 m deep, −5% | 0.95 → 0.81 | break moved |
| Pipeline 511.5, 14 s (re-forms) | 1028 → 1012 m, 10.8 → 9.5 m deep, −14% | 0.36 → 0.32 | break moved |
| Pipeline 511.5, 16 s | 1654 → 1652 m, 4.0 → 4.3 m deep, +3% | 0.96 → 0.76 | surf zone out |
| Coxos 481.5, mid | 940 → 932 m, 3.4 → 3.9 m deep, −3% | 1.30 → 0.89 | break moved |
| Coxos 481.5, low | 962 → 952 m, 3.6 → 4.2 m deep, −3% | 1.01 → 0.81 | break moved |
| Coxos 601.5, mid | 898 → 896 m, 4.8 → 5.1 m deep, −5% | 0.99 → 0.85 | surf zone out |
| Coxos 601.5, low | 970 → 926 m, 3.9 → 5.4 m deep, −2% | 1.19 → 0.77 | break moved |
| Coxos 721.5, mid | 898 m, 1.8 m deep, −7% | none | pass |
| Coxos 721.5, low | 962 → 958 m, 4.9 → 5.0 m deep, 0% | 0.79 → 0.69 | pass |

The break moved past its bounds in 7 of the 12 cases, and the inner surf zone came into range in
1 of the 10 that have one and do not re-form. By the rule fixed beforehand that is a fail, so
steps 3 (2 m lines) and 4 (Coxos in 2D) did not run.

**Why.** Not the shallow water as such, which does take the height out (`H/h` fell in every case),
but where it starts:

- **The break line is not always where the wave breaks.** It is each row's tallest wave before
  the shore. On Pipeline 361.5 and 421.5 at 14 s that is in 7 to 8 m of water, 150 to 240 m
  seaward of where the switch first reads breaking (1580 and 1574 m, in about 2 m). Held from
  there, waves that have not broken steepen into shocks over 200 m and lose too much: 0.37 to
  0.49. On the lab beach the tallest wave and the break coincide, which is why it passed there.
- **Holding the peak's own cell flattens the peak.** In the run without, the wave still grows
  into that cell; in shallow water it starts losing height there. So the height at the break
  point falls (by 19% at Pipeline 361.5) and the tallest point steps one to five cells seaward:
  most of the "break moved" verdicts, by 3 to 15 m. Where breaking starts moved little in most
  cases, but in three of the Pipeline ones (361.5 both periods, 421.5 at 14 s) it moved 38 to
  78 m shoreward, the waves having been flattened before they could break.
- **On Coxos's short insides** (8 to 23 cells from the break to 1 m of water) `H/h` falls by
  0.1 to 0.4 but mostly stays above range, 0.77 to 0.89.

**What it means.** Shallow water held from the sheet's break line overdoes long, gentle surf zones
and underdoes short ones, and it moves the break point it was meant to leave alone. A start
further in (where the switch first reads breaking, or some distance shoreward of the tallest
wave) would answer the first two, but on 3 m cells the switch at Coxos reads breaking only near
the shore (0.6 to 1.7 m of water), so it would hardly act there, and choosing the distance is the
tuning these bounds were set to avoid. `--shallow-inside` stays in, off unless asked for, like
`--breaking-viscosity`. Of the fallbacks, correcting the sheet's surf-zone heights is what is
left without new physics.

## Step 3: Coxos in 2D (heavy, one run at a time, each approved)

Only if step 2 passes.

- **Runs.** Coxos 2 m, 13 s from 300°, mid and low tide, with the viscosity: two runs. A third
  without it, at mid tide, as a fresh baseline: the existing runs predate the 2026-11-17 fix to
  the fifth-order scheme at walls, which changes them slightly near the side walls. About half
  an hour each by the sheet's count (14 runs a night), on the efficiency cores, under the heat
  rules: nominal thermal state, nothing else heavy running, Max's go-ahead for each.
- **Checks, against the runs without it:**
  - the break line moves no more than the step 2 bounds, row by row, in the middle half of the
    width;
  - the inner surf zone's `H/h` falls into the step 2 range;
  - the tide response, by the same per-cell slope as the 2026-11-09 check, falls from 0.78–0.81
    towards 0.5–0.7;
  - the sheet both ways on the same map, for Max to look at: the height map shoreward of the
    break should come out lower, and the takeoffs and breaking zone should stay where they
    were. He knows where Coxos breaks and how big the white water is on the inside.
- **Cost of the extra work in each step,** measured on the first run, since the eddy viscosity
  adds explicit sub-steps.

## Step 4: decide

- **If it holds:** make it the default for dispersive runs. Sheets made before then keep their
  takeoffs, zones and break-point numbers, since those are what the bounds hold still; only
  their height maps and the tide's effect inside the break go stale. Bringing those in line
  means doing a sheet's runs again: about 28 for Coxos and about 14 for Ribeira, a night per
  spot on the laptop, or a rented machine. gogo's copy of each sheet follows, as for F2.
- **If the reefs break early with it:** keep it off, and pick from the fallbacks.
- **Either way,** the README's limits and the sheets' "how far to trust it" lines say what was
  found.

## Fallbacks

- **Viscosity only shoreward of the break.** Take the break line from the switch alone, then
  apply the viscosity only shoreward of it. That needs the break point before the run knows it:
  a first run without, or a second pass. Tried in one form on the lab beach (a front left to
  the switch until it had broken for `5 sqrt(h/g)`): the plunging break moved only to 7.08 m,
  and the spilling breaker broke 0.8 m late with its inner surf zone at 1.26. So not promising
  as tried, but untried in the two-pass form.
- **Plain shallow water inside the surf zone.** Without dispersion the inner surf zone comes out
  at `H/h` 0.65 to 0.70 on the lab beach, since the shock dissipates about the right amount once
  the switch leaves it alone. Switching dispersion off for good shoreward of the break line is
  close to the first fallback, with no viscosity to tune.
- **Leave the model, correct the sheet.** Scale the surf-zone heights on the sheets by the
  measured ratio and say so on the sheet. Cheap and honest about it, but it fixes the numbers,
  not the physics, and a single ratio cannot be right for every reef and tide.

## Not in this plan

- Random waves. The runs are monochromatic, so the comparison is with the laboratory's regular
  waves; the field's 0.5–0.68 is for the top third of random waves, which reads lower again.
- The break point itself, the flat-bed SGN operator and the wide dispersive stencil: other
  limits, listed in the README, with their own effects.
