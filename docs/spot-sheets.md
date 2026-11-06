# Spot sheets: how each spot's wave changes with the swell's angle

A plan, written 2026-11-04; the Coxos runs and the reducer are done (2026-11-06, "Coxos, first
results" below). It replaces
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
- **It is a map.** The sheet shows on a satellite or aerial photo where the swell breaks at each
  direction and tide, and where it breaks first (the takeoffs), with a few lines of text beside
  it. Max asked for that on 2026-11-06, over a table of text.
- **The check is Max.** He looks at the Coxos map, and it goes into gogo only if it matches what
  he knows. He chose not to write his notes down first (2026-11-06), so the check is made after
  seeing the map, not before.

## 1. What a sheet says

One row per direction. Coxos's window is 270–340°, so seven rows at 10° steps from 280 to 340°;
270°, the window's edge, is one more pair of runs if it turns out to matter. Each row at low and
mid tide, Coxos's two tides in `coast.yml`.

| Line on the sheet | From the run | On the 15 m bed today |
|---|---|---|
| Size at the break compared with offshore | `break_line[].height_m` over the offshore height, median over the lineup rows | **Usable.** Shown with the reflection off the rock beside it (Coxos 0.26, section 4 of forecast-tables.md), since the reported height includes it |
| How far out it breaks, and in how much water | `x_m` against the shoreline; `depth_m` | **Usable** |
| How the tide moves it | low against mid at the same direction | **Usable.** On the existing runs (the 320° row), going from mid to low tide moved the break from 4.8 to 4.2 m deep |
| Breaks onto the rock | the share of lineup rows with `shore_break` | **Usable.** 68 of 200 rows on the existing mid-tide run |
| Which way it peels, and whether it closes out | `tools/direction.py` | **Rough.** The direction of the peel is usable; the angle (28° on the existing run) is not, at 15 m |
| Pitches or crumbles | `surf_similarity` | **Not yet.** In the 15 m bed the wave reaches the rock unbroken, so the type is unknown, not "crumbling" |

Every sheet carries the seabed it was run on, its resolution and date, and the words "model,
not observation".

**Which direction a row means.** The rows are labelled with the open-sea direction gogo
forecasts. The shelf between the open sea and the wave-maker (about 25 m deep for this swell)
has already bent the swell towards Coxos's shore normal (265°), and each run is driven square
at the bent direction. By Snell's law over straight contours, for 13 s:

| Open sea | 280° | 290° | 300° | 310° | 320° | 330° | 340° |
|---|---|---|---|---|---|---|---|
| At the wave-maker | 275° | 282° | 288° | 294° | 300° | 304° | 307° |

So the existing Coxos runs, driven at 300° at the maker, are the **320° row**, not the 300° one.
And the shelf squeezes the north end together: 330° and 340° arrive 3° apart, so those rows
should look alike, which is itself something the sheet can say. gogo's ray tracer could give
the bent direction properly later (forecast-tables.md section 3); Snell's law is enough to
start.

## 2. The runs

- **Coxos:** 7 directions × 2 tides = 14 runs. The 320° row at both tides exists already
  (`out/coxos_2m13s_300`, `out/coxos_2m13s_300_low`, named for their 300° at the maker), so 12
  are new. The beds are `fetch_spot.py coxos --from <maker direction>`: `coxos_275` to
  `coxos_307`, all seven built and checked on 2026-11-04 (shoreline inside the box on every row). The swell is 2 m, 13 s,
  the size those two runs used.
- **One width check** at 340° (307° at the maker), the most oblique: the same run on a box
  twice as wide (`--width 2400`). The walls
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

1. **The Coxos seabed as a script. Done 2026-11-04.** `tools/merge_bed.py coxos` lays the 2011
   LiDAR (`data/emodnet_hr_centro/`, ~15 m) over the EMODnet DTM (`data/emodnet_dtm_coxos/`,
   ~115 m, interpolated bilinearly in the LiDAR's gaps) and shifts both from chart datum to mean
   sea level (−2.0 m, checked against Cascais). `[coxos]` in `spots.toml` holds the sources and
   the frame. The bed cut from it matches the hand-made one the runs used to within 0.000002 m on
   every cell, and `tools/tests/test_merge_bed.py` checks that when `data/` is present. Over the
   whole merged file, 0.9% of pixels differ from the hand-made one, all outside the bed: some
   inland land that had been set to 8 m, and a few centimetres in deep water.
2. **A bed per direction. Done 2026-11-04.** `tools/fetch_spot.py coxos --from 288` turns the bed
   for a swell from 288° at the wave-maker and writes `coxos_288`; `--width` widens the box. One
   bed per direction is needed because `wavesim run` drives the swell square-on and the side
   walls only suit that.
3. **A reducer. Done 2026-11-06.** `tools/spot_sheet.py coxos` reads the runs listed under
   `[coxos.sheet]` in `spots.toml` and writes `out/sheets/coxos.json` (240 KB) and a preview map
   beside it. Per run: the breaking zone (cells breaking in at least 8% of the saved frames, about
   every wave), the break points, the takeoffs and a few numbers near the main one.
4. **The lineup: not needed.** The reducer finds the takeoffs itself (below), so Max does not have
   to mark one; he checks the ones it finds.

No batch driver, wall time or `--no-frames` is needed for 14 runs.

## Coxos, first results (2026-11-06)

All 13 runs finished in 5 h 50 min on 3 efficiency-core threads; the thermal state stayed nominal
throughout (701 readings).

- **The width check passes.** At 340° the 1200 m box and the 2400 m one put the break in the
  same place over the trusted rows (median 0 m apart) with heights 4–8% apart, inside the 10%
  limit. Past the north margin line, over the point (rows 900–1050 m), they agree within 16 m and
  1%, so the trusted band is 300–1050 m, not 300–900. Only the last 150 m by the north wall
  differ (19%).
- **It lines up with the photo.** The seabed's 0 m line follows the shoreline in the Esri imagery
  along the cliffs; in the sandy coves it sits up to about 50 m inland, which is the beach between
  the tide lines.
- **Low tide moves the break out** 9 to 39 m depending on direction (median over the rows both
  tides share). Runs crop the bed's offshore end differently, so this compares positions on the
  bed, not in each run's own frame; a first version that did not read 6 to 42 m the wrong way.
- **Takeoffs.** A takeoff is where the crest breaks before the reef on both sides: the crest leaves
  the wave-maker straight, its travel time to each row's break point is summed at the local phase
  speed, and a row is a takeoff if it is the earliest within 60 m and at least 2 s ahead of the
  latest within 150 m on both sides (rows onto the rock and the band's last 30 m never count). A
  coast that simply runs seaward to one end breaks first at that end, which is not a takeoff;
  an earlier version that took the single earliest row picked the band's end that way.
- **The reef's peaks.** Each run has two to five takeoffs, and which one is strongest changes from
  direction to direction. But they cluster: four places are a takeoff in half or more of the 14
  runs (39.0056 −9.4285, 39.0016 −9.4288, 39.0064 −9.4275, 39.0031 −9.4286). Those are features of
  the bed; the swell's direction switches them on and off. The map shows them as dashed rings.
- **Coxos is only part of it.** The runs show the whole coast; Max drew the spot itself on the
  map (2026-11-06), and `area` in `[coxos.sheet]` keeps only what lies inside: the rock ledge
  and the bay south-east of it, about 39.0028 to 38.9985. Takeoffs are still found over all the
  trusted rows, so the outline makes no edges of its own. Inside it, the takeoff sits on the
  outer ledge at about 39.0016 −9.4288 at every direction from 280° to 340° (39.0011 to 39.0022),
  a peak in 8 of the 14 runs, and the breaking zone runs from there south-east into the bay. Low
  tide moves the break 15 to 42 m out (0 m at 340°). The 2026-11-06 peaks above are the whole
  coast's.
- **What is left out.** The peel (which way and how fast) came out inconsistent between
  neighbouring directions and tides, sometimes negative, so the map does not show it; nor
  whether it barrels. Both wait for the 2 m seabed.

## 4. How it shows in gogo

A section on the spot page (`web/src/screens/SpotScreen.tsx`), "How Coxos changes with the
swell's angle": a map with a direction picker and a tide switch, the breaking zone, the break
points and the takeoffs drawn on it, a few lines of text beside it, and the "model" label.
`tools/sheet_preview.html` is that section as a standalone page. In gogo it needs Leaflet (about
40 KB) and an image source: Esri World Imagery in the preview; Portugal's DGT aerial photos
(25 cm, 2018; 10 cm along the coast, 2014–15) would be sharper, but no working tile service for
them was found on 2026-11-06, and their licence is still to be confirmed. The tile provider sees
which tiles a viewer loads. The data is a static file in gogo, the way `src/gogo/data/rays.json`
is, copied from the reducer's output. Nothing in the score reads it, so the 2026-10-14 decision against break-scale physics
in the score is untouched, and nothing waits on labels.

Later, cheaply: mark the row nearest today's forecast direction, so the sheet answers "what
does today's swell do here" without any new physics.

Built 2026-11-06 in gogo: `web/src/SpotSheet.tsx` (the section, on Leaflet), `web/src/sheets.ts`
(the types and the loader) and `web/src/sheets/coxos.json` (copied from `out/sheets/coxos.json`).
Leaflet and a sheet load only on a spot that has one. The photo is DGT's OrtoSat2023 (30 cm, 2023,
whole mainland; its WMS states no fees and no access constraints, though no licence is named),
with Esri to switch to; the 2023 image happens to show waves breaking along the Coxos reef. The
page does not yet open on today's forecast direction: the spot API does not return the swell's
direction or the tide per hour.

## 5. Which spots

The reef and point spots, where the 2011 bottom is still the bottom: Coxos, then Ribeira d'Ilhas
(1.5 km south, the same source tiles), then São Lourenço, Pedra Branca, Lagide and São Pedro do
Estoril. The list and the bottom types are from section 4 of forecast-tables.md, still for Max
to confirm. Beach breaks get no sheet: their bars have moved since 2011, and the sheet would
describe sand that is no longer there.

## 6. The check

Max looks at the Coxos map against what he knows (he chose on 2026-11-06 not to write it down
first, so this is a check made after seeing it):

- whether the takeoffs and the four peaks sit where people actually take off;
- which direction works best, and which is too north or too west;
- where it breaks at low tide against mid tide;
- at what size or tide it starts hitting the rock.

If the sheet agrees, it goes into gogo, and Ribeira is next. If it disagrees, the likeliest
cause is the seabed: a shelf in front of the rock that 15 m does not resolve. Then the sheet
waits for the DGT 2 m data rather than going out wrong.

## Next steps

1. Keep logging sessions. That is gogo's own gate and unrelated to the sheets, but it is still
   what improves the score.
2. Done 2026-11-04: the Coxos seabed as a script, and the bearing override (section 3).
3. Done 2026-11-05: the 13 runs at Coxos.
4. Done 2026-11-06: the reducer, the width check and the preview map.
5. Done 2026-11-06: Max checked the Coxos map, clipped to the spot: "fine enough for now".
6. If it holds: the image source (DGT or Esri), then the map section on gogo's spot page, then
   Ribeira. If not: wait for DGT.
7. When the DGT 2 m data arrives: rebuild the beds, rerun, and add the peel and pitch lines.
