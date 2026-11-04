# Spot sheets: how each spot's wave changes with the swell's angle

A plan, written 2026-11-04. Nothing here is built yet. It replaces
[forecast-tables.md](forecast-tables.md) as the way wavesim reaches gogo. That doc is parked;
its findings still hold, and this one links to them rather than repeating them.

## In short

- **What it is.** For each reef or point spot, a short sheet gogo shows on the spot page. It
  covers a few swell directions across the spot's window, at one typical size and two tides,
  and says what the wave does at each: how big it breaks compared with the swell offshore, how
  far out and in how much water, whether it breaks onto the rock, and which way it peels. It is
  information for the surfer. gogo's score does not read it.
- **Why not the tables.** The tables were for the score, and gogo cannot measure a score change
  yet: half a dozen candidates are already waiting on labels, and a Coxos-only table would only
  move Coxos hours. They would also cost about 600 runs at Coxos and 960 at Ribeira (counted,
  forecast-tables.md section 2). A sheet is about 14 runs per spot: one night on the laptop.
- **How far to trust it today.** On the 15 m seabed, size, where it breaks, depth, how the tide
  moves it and whether it hits the rock are usable. Whether it pitches or crumbles (barrels) and
  how fast it peels depend on the bottom's exact shape at the break, which 15 m blurs. Those
  lines wait for the DGT 2 m data.
- **The check is Max.** He writes down what he knows about Coxos before seeing the runs, and
  the sheet goes into gogo only if it agrees.

## 1. What a sheet says

One row per direction. Coxos's window is 270–340°, so seven rows at 10° steps from 280 to 340°;
270°, the window's edge, is one more pair of runs if it turns out to matter. Each row at low and
mid tide, Coxos's two tides in `coast.yml`.

| Line on the sheet | From the run | On the 15 m bed today |
|---|---|---|
| Size at the break compared with offshore | `break_line[].height_m` over the offshore height, median over the lineup rows | **Usable.** Shown with the reflection off the rock beside it (Coxos 0.26, section 4 of forecast-tables.md), since the reported height includes it |
| How far out it breaks, and in how much water | `x_m` against the shoreline; `depth_m` | **Usable** |
| How the tide moves it | low against mid at the same direction | **Usable.** At 300°, going from mid to low tide moved the break from 4.8 to 4.2 m deep |
| Breaks onto the rock | the share of lineup rows with `shore_break` | **Usable.** 68 of 200 rows at 300° mid tide |
| Which way it peels, and whether it closes out | `tools/direction.py` | **Rough.** The direction of the peel is usable; the angle (28° at 300°) is not, at 15 m |
| Pitches or crumbles | `surf_similarity` | **Not yet.** In the 15 m bed the wave reaches the rock unbroken, so the type is unknown, not "crumbling" |

Every sheet carries the seabed it was run on, its resolution and date, and the words "model,
not observation".

**Which direction a row means.** The rows are labelled with the open-sea direction gogo
forecasts. The shelf between the open sea and the wave-maker (1–3 km out, 23–35 m deep) has
already bent the swell towards the shore. Over straight contours, a 13 s swell from 340°, 75°
off Coxos's shore normal (265°), arrives at the maker about 46° off; one from 280° about 11°
off. The runs use that bent direction. gogo's ray tracer could give it properly later
(forecast-tables.md section 3); Snell's law is enough to start.

## 2. The runs

- **Coxos:** 7 directions × 2 tides = 14 runs. 300° at both tides exists already
  (`out/coxos_2m13s_300`, `out/coxos_2m13s_300_low`), so 12 are new. The swell is 2 m, 13 s,
  the size those two runs used.
- **One width check** at 340°, the most oblique: the same run on a box twice as wide. The walls
  mirror a shore that is not square to them, and beyond 20° off the shore normal that has not
  been tested (forecast-tables.md section 3). If it disagrees by more than 10%, the oblique rows
  get a wider box.
- **Cost:** about 25 minutes per run at background priority on 3 efficiency-core threads, one at
  a time, so the 13 runs fit in one night. That is the setup the laptop compute rule allows, and
  it still needs Max's OK before it starts. Frames for 13 runs are 2–6 GB; nothing needs
  `--no-frames`.

A second size (3 m, 14 s) at three directions is the obvious addition once the first sheet has
been checked.

## 3. What wavesim needs first

Little, because the break reading is done (2026-11-04: shoreward-most near-peak, `reflection`,
`shore_break`):

1. **The Coxos seabed as a script.** `data/coxos/coxos_seabed_msl.tif` was merged by hand: the
   2011 LiDAR (`data/emodnet_hr_centro/`, ~15 m) over the EMODnet DTM (`data/emodnet_dtm_coxos/`,
   ~115 m) for the deep water, shifted from chart datum to mean sea level (LiDAR − 2.0 m, checked
   against Cascais). Nothing in `tools/` rebuilds it, and `data/` is not in git. Turn the merge
   into a tool and add `[coxos]` to `spots.toml`. Done when it rebuilds today's
   `data/coxos/coxos.json` bed, which faces 120° (swell from 300°), to the same numbers.
2. **A bearing override in `tools/fetch_spot.py`**, so one spot gives a bed per direction
   (`data/coxos_280/`, …). It already crops and turns a bed to any bearing from a local GeoTIFF;
   today the bearing comes only from `spots.toml`. One bed per direction is needed because
   `wavesim run` drives the swell square-on and the side walls only suit that.
3. **A reducer, `tools/spot_sheet.py`:** reads the runs of one spot and writes the sheet's rows
   as JSON, over the lineup rows.
4. **The lineup:** the rows (`y` range) and stretch of `x` where people take off at Coxos, drawn
   once on the bed. Max knows where that is; the reducer needs it as numbers.

No batch driver, wall time or `--no-frames` is needed for 14 runs.

## 4. How it shows in gogo

A section on the spot page (`web/src/screens/SpotScreen.tsx`), "How Coxos changes with the
swell's angle": the sheet's rows, the tide difference beside each, and the "model" label. The
data is a static file in gogo, the way `src/gogo/data/rays.json` is, copied from the reducer's
output. Nothing in the score reads it, so the 2026-10-14 decision against break-scale physics
in the score is untouched, and nothing waits on labels.

Later, cheaply: mark the row nearest today's forecast direction, so the sheet answers "what
does today's swell do here" without any new physics.

Whether the section goes on the page at all is Max's call, after he has seen Coxos's sheet.

## 5. Which spots

The reef and point spots, where the 2011 bottom is still the bottom: Coxos, then Ribeira d'Ilhas
(1.5 km south, the same source tiles), then São Lourenço, Pedra Branca, Lagide and São Pedro do
Estoril. The list and the bottom types are from section 4 of forecast-tables.md, still for Max
to confirm. Beach breaks get no sheet: their bars have moved since 2011, and the sheet would
describe sand that is no longer there.

## 6. The check

Before the runs, Max writes down three things he knows about Coxos, so the comparison is not
made after the fact:

- which direction works best, and which is too north or too west;
- where it breaks at low tide against mid tide;
- at what size or tide it starts hitting the rock.

If the sheet agrees, it goes into gogo, and Ribeira is next. If it disagrees, the likeliest
cause is the seabed: a shelf in front of the rock that 15 m does not resolve. Then the sheet
waits for the DGT 2 m data rather than going out wrong.

## Next steps

1. Keep logging sessions. That is gogo's own gate and unrelated to the sheets, but it is still
   what improves the score.
2. Max writes down what he knows about Coxos (section 6) and marks the lineup.
3. Script the Coxos seabed and add the bearing override (section 3, items 1–2).
4. The 13 runs at Coxos, one night, after Max's OK.
5. The reducer, and the sheet compared with step 2.
6. If it holds: the section on gogo's spot page, then Ribeira. If not: wait for DGT.
7. When the DGT 2 m data arrives: rebuild the beds, rerun, and add the pitch and peel lines.
