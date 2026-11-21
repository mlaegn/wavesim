# The surf zone: broken waves that stay too tall

A plan, written 2026-11-21, not started. It follows the 2026-11-17 audit, which named this the
model's main accuracy gap once its bugs were fixed.

## In short

- **The problem.** wavesim puts the break point in the right place, but once a wave has broken
  it hardly loses height on its way to the shore. In the inner surf zone it stays about as tall
  as the water is deep (`H/h` near 1.0). Regular laboratory waves sit at 0.58 to 0.61 there
  (Ting & Kirby's spilling breaker, which the model reproduces 1.66 times too tall), and
  random waves in the field at 0.5 to 0.68 by the Surfline camera fields.
- **What it spoils.** Heights shoreward of the break; the breaking zone on the spot sheets,
  which gogo shows; and how much the tide changes the size there (0.78 to 0.81 on the Coxos
  runs, against 0.50 to 0.68 in the field). The break point, its depth and its height are not
  affected.
- **A fix exists.** `Breaking::viscosity` follows each breaking front from step to step and
  gives it Kennedy et al.'s (2000) eddy viscosity, the closure FUNWAVE uses (`wavecore::eddy`).
  On Ting & Kirby's beach it brings the inner surf zone to 1.04 times the laboratory's, no
  station off by more than 11%. It is off by default because it breaks the plunging laboratory
  wave 0.86 m early, against 0.5 m allowed. It has never run on a reef, and `wavesim run` has no
  switch for it.
- **The plan.** Add the switch; try it on single cross-shore lines of Pipeline and Coxos, which
  is light work; only if the break points hold there, run Coxos in 2D; then decide whether it
  becomes the default. Each step has bounds written down before its first run, as the lab tests
  had.
- **The risk to watch.** Reef breaks plunge, and the plunging laboratory wave is exactly where
  the fix fails. If the reefs break early with it too, the fix stays off and the fallbacks below
  apply.

## Step 1: a switch on `wavesim run` (no runs)

- `--breaking-viscosity` on `wavesim run`, off by default, building the solver with
  `Breaking { viscosity: Some(EddyViscosity::default()), ..Default::default() }`. Only with
  dispersion: the flag with `--dispersive false` is an error, since plain shallow water has no
  switch to follow.
- `run.json` records it (`"breaking_viscosity": true`), so a sheet, a viewer and a hand-off can
  tell which model made a run. A hand-off from such a run is refused by the slice for now: the
  slice's swell along a line runs without it, and the two would disagree.
- Tests: the flag reaches the solver and the header; a short run on the synthetic beach still
  conserves water; the slow lab tests already cover the physics, and run on GitHub every Monday.
- Compute: the fast tests on one core, as for any change.

## Step 2: single lines across real reefs (light)

A run along one cross-shore line, four cells wide, costs about as much as a fast test: a few
minutes on one core, no Docker, no release build. That is enough to see whether the break point
moves, which is the failure the lab showed.

- **Lines.** Pipeline y = 361.5, 421.5 and 511.5 m (in `data/`, where earlier slices ran), and
  three rows of the Coxos bed for the 300° swell (`data/coxos/coxos.json`, `+x` along 120°),
  cut the same way.
- **Waves.** Pipeline 2.5 m, 14 s and 16 s; Coxos 2 m, 13 s, at mid and low tide: the sizes the
  sheets and slices already use.
- **Each line twice**, with and without the viscosity, on the cells the 2D runs use (3 m at
  Pipeline, the Coxos bed's own), recording the height per cell over the last periods as the
  runs' statistics do.
- **Bounds, to be fixed before the first run:**
  - the break point (the tallest point before the shore) moves no more than 10 m, and its depth
    by no more than 10%; for scale, halving the cells moved Pipeline's by 3.7 m;
  - the height at the break point changes by less than 10%;
  - the inner surf zone (between the break and 1 m of water) falls to `H/h` 0.55 to 0.75, the
    laboratory's regular-wave value with room for the reef.
- **Where it goes.** A small `#[ignore]`d test or example beside `tests/lab.rs`, so the numbers
  can be rerun, and a table in this doc.

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
  - the sheet's breaking zone, drawn both ways on the same map, for Max to look at: he knows
    where Coxos breaks and where the white water reaches.
- **Cost of the extra work in each step,** measured on the first run, since the eddy viscosity
  adds explicit sub-steps.

## Step 4: decide

- **If it holds:** make it the default for dispersive runs. Then every sheet's runs need doing
  again: about 28 for Coxos and about 14 for Ribeira, a night per spot on the laptop, or a
  rented machine if that is too much. gogo's copy of each sheet follows, as for F2.
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
  not the physics, and not the extent of the white water.

## Not in this plan

- Random waves. The runs are monochromatic, so the comparison is with the laboratory's regular
  waves; the field's 0.5–0.68 is for the top third of random waves, which reads lower again.
- The break point itself, the flat-bed SGN operator and the wide dispersive stencil: other
  limits, listed in the README, with their own effects.
