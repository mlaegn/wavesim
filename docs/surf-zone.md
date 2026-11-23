# The surf zone: broken waves that stay too tall

A plan, written 2026-11-21. Step 1, the switch, is done (PR 5). Step 2 ran on 2026-11-23 and did
not pass: on reef lines the fix leaves the break point exactly where it was, but it hardly
changes the surf zone either (results under step 2). Steps 3 and 4 are on hold until the next
move is chosen. It follows the 2026-11-17 audit, which named this the model's main accuracy gap
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
