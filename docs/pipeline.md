# Pipeline: what the model finds

Results of the 2D run and the slice solver on the Banzai Pipeline bed (`data/pipeline.json`, the 3 m NOAA CUDEM; see [bathymetry.md](bathymetry.md)), with the lines the slices break on taken from the 2013 USACE 1 m lidar. The short version is in the README, under [What it finds on Pipeline, in short](../README.md#break-one-wave-side-on), and so is how to break a wave on a real line. The first case below was run with:

```bash
uv run --python 3.12 --with rasterio --with numpy --with scipy --with pyproj \
    tools/fetch_line.py data/pipeline.json --y 511.5   # the 1 m lidar along y = 511.5 m, over 50 m
cargo run --release -p wavesim -- run data/pipeline.json --height 2.5 --period 14 \
    --threads 3 --max-minutes 15 --handoff-y 511.5 --handoff-depth 9 --out out/pipeline_handoff
cargo run --release -p wavesim -- slice --line data/pipeline_line_511.5_w50.json \
    --handoff-file out/pipeline_handoff/handoff_y511.5.json --period 14
```

## One line across Pipeline

On that line a 2.5 m, 14 s swell breaks at the edge of the reef flat, and the two models agree on where, within 10 m:

| | 2D run (where the wave is tallest) | Slice from the 2D run's wave (face vertical) |
|---|---|---|
| Where | 2224 m along | 2214 m along |
| Still-water depth | 2.7 m | 2.9 m |
| Wave | 2.3 m crest to trough | crest 1.73 m above still water |
| Lip and tube | | thrown 7.6 m, tube 0.34 m tall and 0.38 m wide |
| Computing | 12 min on 3 efficiency-core threads | 12 min on one efficiency core |

A modest plunge with a small tube: this stretch of the coast is one the reef spreads energy away from (the 2D run's biggest waves, about 5 m, break elsewhere).

On the line through those biggest waves, y = 301.5 m (handed over at 9 m, its crest 0.31 of the depth), the wave stands up with a face vertical 2159 m along, in 6.2 m of water, under a crest 4.6 m above still water (the 2D run: tallest at 2180 m, 4.9 m deep), and keeps a near-vertical face for about 30 m while only a small curl, under a metre across, forms at its top and falls onto the face: a heavy, crumbling break rather than a hollow tube, as the surf-similarity rule calls most of that stretch. At 0.6 m node spacing (13 minutes) the curl lands; at 0.4 m (19 minutes) the face stands vertical at the same moment and place and the curl grows the same way, 0.34, 0.50 and 0.65 m across, until its thin tip tangles and the run stops itself before it lands. For a break like this "face vertical" marks the first small overturn and "lip thrown" how far the steep face travelled before the curl fell (33 m); they were made for clean plunges. macOS reported no thermal warning, slowdown or power limit at any reading, 30 s apart, through both runs.

A bigger swell, 3.5 m and 16 s (the 2D run took 19 minutes on 3 threads), turns that heavy crumble into a plunge on the same line. Handed over at 9 m with its crest 3.8 m up (0.30 of the total depth, just where the main model's breaking switch starts to act), the slice's face stands vertical 2147.5 m along, in 6.7 m of water, under a crest 5.5 m above still water (the 2D run: tallest at 2156 m, 6.3 m deep, 6.2 m crest to trough). The lip is thrown 16 m and lands 2.6 m above still water, closing off a tube 1.0 m tall and 0.8 m wide: a real barrel, but a narrow one for a wave this size. 0.6 m node spacing, 8 minutes on one efficiency core; no thermal warning at any of 9 readings.

The same swell on y = 511.5 m, where the 2.5 m swell gave only a modest plunge, makes the hollowest wave of all. Handed over at 9 m with its crest 2.0 m up (0.18 of the total depth, the breaking switch off), the face stands vertical 2185 m along, in 3.9 m of water, under a crest 3.8 m above still water (the 2D run: tallest at 2170.5 m, 4.8 m deep, 4.1 m crest to trough, plunging by the surf-similarity rule). The lip is thrown 16 m and closes off a tube 2.5 m tall and 3.7 m wide, 4.5 m² of air: thirteen times the tube on the bigger wave at 301.5. Here the reef rises more steeply before the break (1 in 20 against 1 in 27 at 301.5, over the half wavelength seaward), and the wave stands up in water about as deep as its crest is high (at 301.5 the crest is 0.8 of the depth). At 0.4 m node spacing instead of 0.6 the face stands vertical at the same moment and place (within 0.01 s and 0.1 m) under the same crest, and the tube comes out 5% bigger (2.6 by 3.8 m, 4.7 m²). 15 and 25 minutes on one efficiency core, after a 21-minute 2D run on 3 threads; no thermal warning at any of 62 readings.

## Seven lines along the reef

Along the reef, one 2D run (18 minutes) saved seven lines across the middle half of the width, and each was sliced at 0.6 m (6 to 15 minutes each; no thermal warning at any of 73 readings). The 2D run repeats exactly: its saves at 301.5 and 511.5 are the same bytes as before.

| y (m) | ξ₀ | Face vertical | Depth | Crest | Lip thrown | Tube (tall × wide) | 2D run: tallest |
|---|---|---|---|---|---|---|---|
| 301.5 | 0.40 | 2147.5 m | 6.7 m | 5.5 m | 16 m | 1.0 × 0.8 m | 2155.5 m, 6.3 m deep |
| 421.5 | 0.26 | 2112.7 m | 6.6 m | 4.9 m | 15 m | 0.7 × 0.6 m | 2116.5 m, 6.4 m deep |
| 511.5 | 0.54 | 2185.2 m | 3.9 m | 3.8 m | 16 m | **2.5 × 3.7 m** | 2170.5 m, 4.8 m deep |
| 571.5 | 0.52 | 2181.4 m | 3.8 m | 2.2 m | 7 m | 0.7 × 0.7 m | 2188.5 m, 3.7 m deep |
| 661.5 | 0.40 | 2198.4 m | 4.2 m | 2.6 m | 7 m | 0.7 × 0.8 m | 2197.5 m, 4.0 m deep |
| 751.5 | 0.50 | 2147.9 m | 5.3 m | 4.6 m | 14 m | 0.9 × 0.7 m | 2167.5 m, 5.1 m deep |
| 841.5 | 0.47 | 2190.5 m | 5.5 m | 4.0 m | 13 m | 0.8 × 0.6 m | 2194.5 m, 5.4 m deep |

Every line plunges, and the slice and the 2D run agree on where within 20 m (mostly within 8). Only 511.5 makes a big tube; everywhere else the lip closes a tube under a metre across, whether the wave is 2 m or 5.5 m. On 511.5 the reef falls from 9 to 5 m deep over 50 m just before the break, about 1 in 13, against 1 in 20 to 1 in 25 on the other lines, and a fair-sized wave (4.1 m crest to trough in the 2D run) arrives over it: a likely reason, not a tested one. Lines 60 m apart cannot say how wide that stretch is.

ξ₀ is Battjes' (1974) surf similarity from the deep-water swell, `tan β / sqrt(H0 / L0)` with `L0 = g T² / 2π` (400 m here) and `tan β` the reef's mean rise over the half wavelength seaward of the 2D run's break: below 0.5 a wave should spill, from 0.5 to 3.3 plunge. All seven lines sit near that border, and the slices bear it out in kind: every one overturns, but only just, except 511.5. The number does not tell the barrel apart, though: 571.5 has nearly the same ξ₀ as 511.5 and closes a tube a quarter as tall, because the reef brings it a wave with a crest little more than half as high. ξ₀ sorts the kind of break, not the size of the tube.

## The barrel is one line wide

The barrel is one line wide. Slices 30 m either side of it, from another run of the same swell (21 minutes, then 14 and 12), close tubes only 0.5 m across: at y = 481.5 the face stands vertical 2152.9 m along in 6.0 m of water under a 4.2 m crest, at 541.5 2171.3 m along in 4.4 m under 3.4 m. The 2D run, on 3 m cells, changes smoothly across those 60 m (its tallest wave 4.3 to 3.8 m, all at 2170 to 2174 m), so the difference comes from the 1 m line. On 511.5 the reef falls away to 9.7 m just seaward and rises 1 in 13 to the break; on 481.5 and 541.5 it is 7.6 m deep there and rises 1 in 25. A slice treats its line's bed as if the reef were the same all along the shore, so a feature this narrow counts in full in a slice and only in part in the sea, where a wave crest spans far more than 60 m. The barrel at 511.5 therefore shows what that profile does, not what the reef there does.

`tools/fetch_line.py` (50 m unless `--width` says otherwise) makes the line the average of parallel lines a metre apart across a strip 50 m wide, so that a narrow feature counts in proportion. On that bed the barrel is gone: the tube is 0.7 m tall and 0.6 m wide, like its neighbours'. Averaged across 20 m it is the same (0.7 by 0.7 m), so the pocket behind it is narrower than 20 m.

| y (m) | Bed | Face vertical | Depth | Crest | Tube (tall × wide) |
|---|---|---|---|---|---|
| 481.5 | one line / 50 m strip | 2152.9 / 2157.2 m | 6.0 / 6.0 m | 4.2 / 4.3 m | 0.5 × 0.6 / 0.6 × 0.6 m |
| 511.5 | one line | 2185.2 m | 3.9 m | 3.8 m | 2.5 × 3.7 m |
| 511.5 | 20 m / 50 m strip | 2192.7 / 2192.5 m | 3.8 / 3.9 m | 3.6 / 3.5 m | 0.7 × 0.7 / 0.7 × 0.6 m |
| 541.5 | one line / 50 m strip | 2171.3 / 2171.1 m | 4.4 / 4.5 m | 3.4 / 3.4 m | 0.5 × 0.4 / 0.7 × 0.6 m |

Where the wave breaks, how deep and how tall hardly move with the averaging: by at most 8 m, 0.1 m and 0.3 m. Only the tube on 511.5 does. So on this swell the slices give heavy plunges with tubes under a metre across all along this stretch. That is much smaller than the barrels Pipeline is known for, and a likely reason, untested, is what a slice leaves out: the real wave is shaped across the reef as well as along it, and peels along the shore.

## Heavy and mellow swells

Four swells on two of those lines, 50 m strips, each from its own 2D run (9 to 35 minutes on 3 threads), all handed over at 9 m, and sliced at 0.6 m (2 to 24 minutes each; no thermal warning at any reading):

| Swell | y (m) | ξ₀ | Hand-off (crest / total depth) | Face vertical | Depth | Crest | Lip thrown | Tube (tall × wide) |
|---|---|---|---|---|---|---|---|---|
| 1.5 m, 16 s | 421.5 | 0.43 | 9 m (0.08) | 2233.9 m | 2.3 m | 2.2 m | | the lip did not land |
| 1.5 m, 16 s | 511.5 | 1.13 | 9 m (0.03) | 2299.6 m | 0.8 m | 0.9 m | 2 m | 1.0 × 1.1 m, against the slice's shallow wall: not a clean break |
| 2.5 m, 12 s | 421.5 | 0.40 | 9 m (0.16) | 2170.6 m | 4.3 m | 3.3 m | 10 m | 0.6 × 0.6 m |
| 2.5 m, 12 s | 511.5 | 0.39 | 9 m (0.09) | 2209.6 m | 3.3 m | 2.3 m | | the lip did not land |
| 2.5 m, 16 s | 421.5 | 0.54 | 9 m (0.20) | 2183.9 m | 3.7 m | 4.2 m | 11 m | 0.5 × 0.4 m |
| 2.5 m, 16 s | 511.5 | 0.33 | 9 m (0.09) | 2231.6 m | 3.0 m | 2.1 m | | the lip did not land |
| 3.5 m, 16 s | 421.5 | 0.21 | 9 m (0.29) | 2116.1 m | 6.5 m | 5.0 m | 15 m | 1.1 × 0.8 m |
| 3.5 m, 16 s | 511.5 | 0.51 | 9 m (0.19) | 2192.5 m | 3.9 m | 3.5 m | 11 m | 0.7 × 0.6 m |

Here ξ₀'s slope is taken over the half wavelength seaward of the slice's vertical face. The bigger the swell, the further out and the deeper it breaks: on 421.5 from 2234 m in 2.3 m of water to 2116 m in 6.5 m; on 511.5 less so. The same 2.5 m swell is far bigger on the reef at 16 s than at 12 s (a 4.2 m crest against 3.3 m on 421.5), and 421.5 gathers more of every swell than 511.5. Every lip that lands closes a tube about a metre across or less, so on this stretch the size of the wave varies far more than its shape. Three lips do not land: the face stands vertical, a small curl forms at the top, and it is still in the air when the slice stops it. All three have ξ₀ between 0.33 and 0.43, under the 0.5 where a wave should spill, and a small curl falling from the top is how a spilling break starts (as on the plane beaches below); but ξ₀ does not sort them cleanly, for a lip lands at 0.40, and the 3.5 m swell on 421.5 has the lowest ξ₀ of all and the biggest tube. The 1.5 m swell is smallest on 511.5, its crest only 0.26 m up at the 9 m hand-off, and it reaches the shallow end of the line before it stands up, so its lip lands against the slice's wall.

Handing over at 5 m on 421.5 first gave a 2.6 by 3.0 m tube for the 2.5 m, 16 s swell, from a crest already 0.37 of the total depth beneath it at the hand-off. At 0.4 m node spacing it came out the same (2.7 by 3.0 m), but handed over at 6 m, the crest 0.33 of the total depth, the face stands vertical 4 m further out under the same crest and the tube is 0.6 by 0.5 m. Past 0.30 the 2D run's breaking switch has begun to turn its dispersion down, so the wave is saved already part broken, and the tube is the first thing to show it; where, how deep and how tall barely move. A slice now refuses such a wave and says so. This one could not be handed over cleanly at first: at 6.5 m its crest is 0.31, and at 7 m the hand-off point was too near the edge of the run's fine grid, which started in 6.2 m of water. The 2D run now starts its fine grid as far offshore as a hand-off point needs (here 111 m further, which leaves its waves as they were: the same break points, heights within 2 cm), and handed over at 9 m, at 0.20, the wave closes a tube 0.5 by 0.4 m, like the rest. Every hand-off in the table is now under 0.30.

The hand-off depth moves the break by tens of metres. Most of the smaller swells were first handed over at 4 and 5 m, when a fine grid still started in shallow water. Against those, a 9 m hand-off moves where the face stands vertical by 7 to 32 m, mostly shoreward, and, but for the 1.5 m swell on 511.5, raises the crest then by 0.3 to 0.7 m. On 421.5 the 1.5 m, 16 s swell stood up at 2221.3 m under a 1.7 m crest and closed a tube 0.6 by 0.7 m; from 9 m it stands up 13 m further in under 2.2 m and its lip does not land. The 2.5 m, 12 s swell stands up 16 m further out under 3.3 m rather than 2.6 m, and closes the same tube, 0.6 by 0.6 m against 0.5 by 0.6. On 511.5 the 2.5 m, 16 s swell, from 5 m, closed a tube 0.5 by 1.1 m 27 m further out; the 1.5 m swell, from 4 m, stood up at 2268.0 m in 2.5 m of water and did not land either. Where the 2D run's tallest wave marks the break, the two hand-offs fall either side of it: 6 m from each for the 1.5 m swell on 421.5, 11 and 16 m for the 2.5 m, 16 s swell on 511.5. (For the 12 s swell the 2D run's tallest wave is far offshore, over the outer reef, so it marks nothing here.) A slice from deeper carries its one wave further over still water, without the backwash of the wave ahead and without the reef spreading it sideways, and that adds up; from shallower it gets a wave nearer breaking. Neither is plainly the better, and the table uses 9 m throughout. The tube, the point of it, does not depend on it: from either depth every lip that lands closes a tube about a metre across or less.

## Will it barrel? The surf-similarity rule on plane beaches

ξ₀ says a wave spills below 0.5 and plunges from 0.5 to 3.3. To see what that means for the tube, one swell, 2 m and 12 s (`L0` 225 m), runs onto plane beaches of five slopes, 600 m of water 25 m deep and then the slope, brought in along the line alone (exact on a beach the same all along the shore) and handed over at 6 m, or 8 m on the two steepest, where at 6 m the crest was already 0.31 of the total depth on one of them. Each took 3 seconds to 5 minutes on one efficiency core; no thermal warning at any reading.

| Slope | ξ₀ | Face vertical in | Crest | Lip thrown | Tube (tall × wide) | Air in it |
|---|---|---|---|---|---|---|
| 0.028 (1 in 36) | 0.30 | 4.6 m | 3.1 m | 8 m | 0.4 × 0.6 m | 0.10 m² |
| 0.05 (1 in 20) | 0.53 | 4.0 m | 2.9 m | 8 m | 0.7 × 0.7 m | 0.18 m² |
| 0.10 (1 in 10) | 1.06 | 3.3 m | 2.1 m | 9 m | 1.8 × 2.8 m | 2.78 m² |
| 0.151 (1 in 6.6) | 1.60 | 3.0 m | 1.4 m | 7 m | 1.5 × 2.7 m | 2.62 m² |
| 0.208 (1 in 4.8) | 2.21 | 3.1 m | 2.8 m | 3 m | 1.2 × 2.8 m | 2.52 m² |

Every one overturns, a slice computing the flow as frictionless, but the size of what it closes changes by nearly thirty times, and not gradually: under ξ₀ 0.5 a curl less than a metre across falls from the top of the wave, which is what a spilling break is at its start, and by 1.0 the lip throws over a tube nearly 2 m tall and 3 m wide. Steeper still the tube stays about that size while the lip is thrown less far: the wave stands up later, nearer the shore, and pitches down more than out. At 0.4 m node spacing instead of 0.6 the tube on the 1 in 10 beach is the same (2.85 m²), and on the 1 in 20 beach it stays small (0.07 m²). Handed over at 8 m instead of 6 the two steepest stand up in the same place, within 2 m, and close tubes about a third larger.

So the rule holds here in the sense that matters for a barrel: a hollow tube needs ξ₀ near 1 or more. Every Pipeline line above has ξ₀ between 0.17 and 0.55, on averaged strips, and their tubes are the size this table gives there. Either the real reef is steeper where the wave stands up than a 50 m strip and a half wavelength average show, or the barrel comes from what a slice leaves out; this does not say which.

The first of those does not hold for this seabed. Right where the waves stand up, the slope over the last 10 and 20 m on 421.5 and 511.5 is 0.03 to 0.06, the same as over the half wavelength, on a single 1 m line and on strips 10, 20 and 50 m wide alike; the narrower strips find no hidden step. The steepest 10 m anywhere in that half wavelength rises 0.05 to 0.10, which would give ξ₀ 0.65 to 1.08, but a stretch a tenth of a wavelength long is not a beach the wave meets as a whole. Across the whole trusted width, 200 rows of the 3 m bed, the 50 m seaward of where the reef is 4 m deep rises at most 1 in 18 (ξ₀ 0.6 for the 3.5 m swell), and a 9 m stretch steep enough for ξ₀ 1 turns up on one row. The lidar has no gaps there, so these are measured slopes, not filled ones. The steep parts of the profile lie shoreward: from 6 to 3 m the reef rises about 1 in 20 to 30, from 3 to 2 m it is nearly flat (1 in 40 to 120), and from 2 m to the shore it rises 1 in 4 to 14. On this survey, then, Pipeline's waves meet a moderate slope where they stand up, and the slices' small tubes are what such a slope gives; the much bigger barrels of the real place come from something a slice leaves out (the shape of the reef across the wave, which focuses it and makes it peel), from the reef being steeper than a 2013 survey shows, or from the sea being bigger or longer than these swells.

## A bigger sea: 5 m, 16 s

A bigger sea, 5 m and 16 s, needs its wave-maker in 42 m of water, beyond the end of the bed; it ran on a copy fetched 3 km offshore instead of 2.4 km (`offshore_m = 3000`; over the old bed it is the same to 2e-7 m, and its x is the old x plus 600 m, given below in the old frame). The 2D run took 57 minutes on 3 threads (estimated 47), handing over at 12 m. That swell does not reach the inner reef whole: it breaks on the outer reef, 450 to 570 m further out than the 3.5 m swell, in 9 to 11 m of water, the 2D run's tallest waves 8.5 to 8.8 m crest to trough. On 421.5 the slice (handed over at 0.28 of the total depth) stands up 1547.8 m along, in 9.4 m of water, under a 5.3 m crest, 23 m from the 2D run's tallest wave; the lip is thrown 19 m and closes a tube 1.0 by 1.2 m. On 511.5 the outer reef rises to 10.7 m, and the water deepens again behind it, so the 12 m hand-off point lies 57 m shoreward of where the 2D run's wave broke: the slice gets a wave the 2D run has already broken once, which stands up again at 1865 m, in 10 m of water, and reaches the slice's far wall before its lip lands. That one is no result, and a run now says so: of every line saved on Pipeline so far, it is the only one whose waves broke seaward of the hand-off point. A bigger sea, then, moves the break out to an outer reef rather than hollowing the inner one, and there too the tube is about a metre across.

## Is it the direction?

A slice cuts the reef along the bed's `+x`, 130°, and the real wave may cross it at an angle, meeting a different slope along its own path. The 2D runs say where each wave goes (`tools/direction.py out/<run>` works it out from what a run saved): over its last four periods the swell's phase at its own frequency changes from cell to cell, and its gradient points the way the wave travels. At the break point of each of the 200 rows across the trusted width, averaged over the half wavelength seaward of it, for the 3.5 m and the 2.5 m, 16 s swells:

| | 3.5 m, 16 s | 2.5 m, 16 s |
|---|---|---|
| Wave direction from `+x`: median, middle half | −0.7°, −3.2° to +2.3° | +0.9°, −2.3° to +3.4° |
| Rows more than 10° off | 24 of 200 | 21 of 200 |
| Slope along the wave over slope along `x`: median, middle half | 1.00, 0.98 to 1.01 | 1.00, 0.99 to 1.02 |
| Steepest slope along the wave at any break | 0.054 | 0.055 |
| Peel angle (crest against the line of break points, fitted over 63 m): median, middle half | 23°, 7° to 46° | 19°, 13° to 36° |
| Height at the break, and its largest change over 60 m along the shore | 3.1 to 6.2 m, 1.8 m | 2.4 to 5.3 m, 2.1 m |

The waves arrive square to the slices, so a slice along each wave's own path would meet the same slope, within 2% on most rows, and no steeper than 1 in 18 anywhere: the direction is not what keeps the tubes small. What a slice does leave out shows in the last two rows: the waves break in turn along the shore, not all at once (their crests lie at about 20° to the line they break along, so they peel), and their height changes by up to 2 m over 60 m of shore, so a crest is far from uniform along its length, which a slice assumes it is. Whether either makes a bigger tube needs a model that is three-dimensional at the break.

## From another direction

All the swells above come from 310°, square to the shore (the bed cut along 130°). The same 3.5 m, 16 s swell from 290°, 20° further west, ran on a bed cut along 110° (fetched 3 km offshore and 200 m inland, since the shore now crosses the box at 20°; the 2D run took 25 minutes on 3 threads):

| | From 310° | From 290° |
|---|---|---|
| Height in 15 to 25 m of water | 3.7 to 4.1 m | 3.6 to 4.2 m |
| Height in 8 m of water | 4.4 m | 6.5 m |
| Where it breaks: depth, median | 5.4 m | 8.6 m |
| Height at the break, median | 4.9 m | 6.8 m |
| Surf similarity at the break (the 2D run's own, from the height there), median and 10th to 90th percentile | 0.36, 0.19 to 0.51 | 0.10, 0.00 to 0.33 |
| Rows that plunge by it, of 200 | 88 | 1 |
| Peel angle: median, middle half | 23°, 7° to 46° | 53°, 29° to 72° |

Offshore the two are the same; the difference builds over the outer reef, 8 to 12 m deep, where the bed is uneven (on half its cells the slope faces between −54° and +66° of `+x`), and the 290° swell meets more of it head on and is focused: by 8 m deep it is half as tall again. It breaks there, further out and in deeper water, on the outer reef's gentle slope, and in the 2D run almost every row spills. Like the 5 m swell, then, a more westerly swell moves the break out rather than hollowing the inner reef. The side walls, which mirror a bed whose shore now crosses them at an angle, are not the cause: on a box twice as wide (70 minutes), over the same cells, the heights agree within 9% (median), the wide box's the taller, and its breaks are as deep (8.4 m) and as gentle (0.13). Only where along a row the tallest wave stands moves, by 50 m (median), because over the outer reef the waves stay near their tallest for a long way. No slice was run: handed over at 9 m the wave was already 0.39 of the total depth beneath it, part broken, and both slices refused it; at 12 m the fine grid would start 700 m further out and the 2D run take about two hours.

## Against the real place

Against the real place, the check can only be coarse. Published positions for Pipeline disagree with each other by up to 250 m along the shore: in this frame they fall at y = 326 m (Wikipedia), 410 m (Mondo Surf, Pipeline) and 576 m (Mondo Surf, Ehukai). The barrel line, 511.5, lies among them, as does the heavy break at 301.5, so the comparison neither confirms nor rules out that the model's barrel is the real one. Across the shore, the 3.5 m swell breaks on 511.5 about 115 m out from where the water is half a metre deep, in 3.9 m of water. Pipeline's First Reef is usually described as about 60 m out and 0.6 to 2.4 m deep. The 2.5 m swell on the same line breaks closer to that, 85 m out in 2.9 m of water. That fits bigger swells breaking further out, but the descriptions are secondary sources, not a survey.
