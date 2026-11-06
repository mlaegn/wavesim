# A 3D box around the break (scope)

This is a plan, not a feature. It scopes handing the 2D run's swell to a three-dimensional
solver over a small box of reef, so that a wave can overturn and peel at the same time. Nothing
here is built yet. The numbers marked *estimate* still have to be measured.

## The question it answers

The slices break one line of reef at a time. On Pipeline, every line closes a tube about a metre
across or less, except the single 1 m lidar line at y = 511.5 m, which goes away once the bed is
averaged across 20 m. The plane beaches say why: a hollow tube needs ξ₀ near 1, and Pipeline's
lines give 0.17 to 0.55. The README leaves two explanations open:

1. the real reef is steeper where the wave stands up than the survey shows, or
2. the barrel comes from something a slice leaves out. The waves peel, with their crests about
   20° to the line of break points, and their height changes by up to 2 m over 60 m of shore.

A slice cannot test the second. A 3D box can, and on the same 1 m lidar it holds the first fixed.
That gives one question:

> **On the same bed and the same swell, is the tube in 3D bigger than the slices' tube?**

If it is, the peel is what makes Pipeline hollow, and the slices were never going to show it. If
it is not, explanation 2 is ruled out and the answer lies in the bed.

## Where the box sits

```mermaid
graph LR
    C[coarse 2D run<br/>6 m cells, shelf] --> F[fine 2D run<br/>3 m cells, reef]
    F -->|hand-off line| S[waveslice<br/>2DV, one line]
    F -->|box edges, every step| B[3D box<br/>Basilisk, VOF on an octree]
    L[1 m lidar patch] --> B
    L --> S
    B --> M[tube per metre of shore<br/>same measure as a slice]
    B --> V[viewer: box.html]
```

The box is a third level of nesting, set up like the coarse-to-fine one that already exists. The
fine run already feeds a relaxation zone from the coarse one. Here, the fine run would record
the sea along three edges of the box for the 3D solver to use as inflow.

## Which solver, and not something written here

| Option | What it would take | Verdict |
|---|---|---|
| Write a 3D two-phase Navier–Stokes VOF solver in Rust | Years of work, and a solver that has to be validated from scratch | No |
| Extend `waveslice` to a 3D boundary element method (Grilli et al. 2001) | The 2D dense solve already takes 8 to 27 min per break and grows with the cube of the node count. In 3D the node count is roughly squared, so it would need a fast multipole method first | No: a research project in itself |
| OpenFOAM (`interFoam`, VOF, ESI release v2512) | Glue code only: export inflow, mesh the bed, measure the result | The first choice, dropped after V1: its tube stayed at a third to a half of waveslice's and the measured ones at every cell size, solver and time step tried |
| Basilisk (octree VOF, Popinet) | Adaptive mesh, so cells only where the surface is, which is far cheaper for a lip. Real seabeds are less well supported, and the community is smaller | **Yes, after V1** (2026-11-05): it gives waveslice's tube, which OpenFOAM did not, at a seventh of the cells. See the V1 results |
| SPH (DualSPHysics) | Fast on a CUDA GPU, which this machine does not have | No |

This keeps to the project's rule of verification first. The new code is the coupling, and the
coupling is what gets tested. OpenFOAM's VOF has been validated against breaking waves by
others. What has to be shown here is that it agrees with `waveslice` where both apply. OpenFOAM
did not (V1, below), and Basilisk did, so the box is planned in Basilisk. The build plan below
was rewritten for it on 2026-11-06. The mesh and cost sections still hold the OpenFOAM
estimates until the 3D strip test replaces them.

**Physics the box adds that a slice lacks:** the second horizontal direction (peel, the crest
varying along its length, flow along the shore), air, viscosity, and life after the lip lands
(splash-up, the tube collapsing, foam). **Physics it gives up:** nothing that a slice up to
landing gets right. VOF smears the surface over one or two cells, though, so a tube needs about
eight cells across before it can be measured. That is the cost driver below.

## The first case: 3.5 m, 16 s on y = 511.5 m

From `out/pipeline_big_map` (the hand-off line at y = 511.5 m):

| Along `x` | Still-water depth |
|---|---|
| 1774.5 m | 12 m |
| 2095.5 m | 9 m (the slices' hand-off) |
| 2143.5 m | 6 m |
| 2185.5 m | 4 m (the slice's face goes vertical at 2185–2193 m) |
| 2305.5 m | 0.5 m |

At 16 s, the wavelength is 147 m in 9 m of water and 101 m in 4 m.

**Box, in the bed's frame:**

- **`x`:** 2040 to 2320 m (280 m). The inflow edge sits about half a wavelength outside the 9 m
  point, and the box runs to the beach.
- **`y`:** 431.5 to 591.5 m (160 m). This takes in the lines at 481.5, 511.5 and 541.5 m, with
  40 m of buffer to each side. With a 20° peel, a crest breaks along about 55 m of shore for each
  20 m it moves forward.
- **`z`:** −11 to +9 m above still water. That covers the 10.5 m deep inflow and a 5.5 m crest,
  plus room for the splash.

**Cells and cost** (*estimate*, from V0's 3D strip, 2026-11-06; the OpenFOAM tiers this
replaces assumed 0.5 m and 0.25 m cells, which V1 showed are far too coarse):

The strip needed 165 cells per m² of sea surface at 0.25 m on the surface. The finer levels are
assumed to scale with the surface's cell count, so 660 per m² at 0.126 m and 2,650 at 0.063 m.
That assumption is untested in 3D; a strip at 0.126 m (2 to 3 h) would test it. A step costs
16.6 µs per cell on one core, and the time step halves with the cell size: about 0.0025 s at
0.126 m and 0.001 s at 0.063 m during a break. Only the 100 m of the box around the break
(x 2140–2240 m, 16,000 m²) needs the fine cells; the other 28,800 m² stays at 0.25 m.

| Tier | Surface cells at the break | Cells | Tube it can measure | Core-hours for one wave (15 s) | 48 cores |
|---|---|---|---|---|---|
| A | 0.126 m | about 15 M | about 1 m or more (V1: 23% small on a 3 m tube) | about 400 | about 10 h |
| B | 0.063 m | about 47 M | about 0.5 m or more (V1: 19% small) | about 3,200 | about 3 days |

Both are beyond the laptop: one core would take weeks even for tier A, and its memory (about
3 KB per cell in 3D) would need 45 GB. Both need a rented many-core machine. Basilisk runs
there with MPI from the same source. The 18 s before the first crest reaches the break can run
at 0.25 m, about 120 core-hours. Tier A answers the yes-or-no question: the slices' tubes are
0.7 m across, and a 2 m tube would be news. Tier B puts a size on the answer.

Three ways to bring this down, each testable in 2D first, in minutes:
- refine only where the surface is steep or overturning, not along all of it: at any moment
  most of the 100 m is not breaking. This is likely the biggest saving, a factor of 3 or more;
- stop refining the bed to the finest level, which the 2D and strip runs did;
- a box narrower along the shore, if V4 shows the peel settles in less than 160 m.

These are guesses, good to a factor of about three.

## What has to be built

### In wavesim (Rust)

1. **`wavesim run --box X0,X1,Y0,Y1`** records the fine grid along the box's offshore edge and its
   two sides at every solver step, from the start of the settled sea to the end of the run. It
   writes `box/boundary.json` and `box/boundary.f32` in a new format, `wavesim-boundary` version
   1, documented in [data-format.md](data-format.md). The fields per edge point are `eta`, `ubar`,
   `vbar`, the still depth and the bed slope. The run refuses a box that does not lie inside the
   fine grid's trusted width (outside `side_margin_cells`), or whose offshore edge lies seaward of
   the fine grid's relaxation zone. The data is small: about 600 edge points, five fields, and a
   few thousand steps come to tens of MB.
2. **The velocity profile in depth, in 2D.** `handoff.rs` already uses the SGN profile at the
   surface (`u = ubar − H²/3 ubar_xx`). The box needs it at every height `ζ` above the bed:
   `u(ζ) = ū + (H²/6 − ζ²/2) ∇(∇·ū)` and `w(ζ) = −ζ ∇·ū − ū·∇h`. This belongs in `wavecore` as a
   pure function and is tested against the surface rule (`ζ = H`) and against linear theory.
3. **One measure for slices and box.** `cases/v1-steep/measure.py` measures face vertical,
   landing and tube from a VOF surface in the way `waveslice::Tank::tube` does. For V1 it was
   checked against waveslice's own numbers only on the cases it compared. Before P1, either it
   is made to measure waveslice's surfaces too and agree with `Tank::tube` on them, or
   `wavesim measure` takes a surface and runs `Tank::tube` itself. Either way, the 3D tube and
   the slice tube are measured the same way.

### In tools (Python, uv)

4. **`tools/fetch_patch.py`** crops the 2013 USACE 1 m lidar over the box in the bed's frame. It
   reuses `fetch_spot.py` with the box's limits and 1 m cells, and writes the bed as a grid of
   heights (`patch.json` plus `patch.f32`), which the case reads and interpolates bilinearly.
   Basilisk needs no STL or bed-following mesh: the bed is an embedded boundary cut through the
   octree's cells. Outside the lidar's coverage it falls back to the 3 m bed and records where
   it did so.
5. **`cases/box/box.c`**, a Basilisk case grown from `cases/v1-basilisk/v1.c`, plus
   **`tools/box_case.py`**, which writes its inputs from the patch and the boundary record:
   - The set-up V1 passed with: two-phase Navier–Stokes, momentum-conserving VOF, reduced
     gravity, surface tension, the bed embedded (`solid (cs, fs, y - zb(x, z))`, with `y` up
     in Basilisk), laminar, CFL 0.5. The octree refines on the surface, the bed and the
     velocity, down to the cell size V1 passed at on the surface (0.063 m; 258 m / 2^12 there),
     with no cell larger than 0.5 m in the water.
   - The inflow as a **relaxation zone** along the offshore edge, about half a wavelength wide:
     every step, inside the zone, `f` and `u` are blended toward the 2D run's sea (`f` from
     `eta`, `u` from item 2), with a weight rising from 0 at the zone's inner edge to 1 at the
     wall. Basilisk has nothing like OpenFOAM's `timeVaryingMappedFixedValue`, and a zone is
     the better choice anyway: it also absorbs what the beach reflects, which a fixed inflow
     does not (see the risks below). The zone is one event of a few dozen lines that reads the
     boundary record and interpolates it in time and along the edge.
   - The sides: either the same zones, fed from the 2D run, or free-slip walls. V3 and V4
     decide.
   - An open top (`u.n` Neumann, `p` zero), at a cell edge of the octree's coarse levels.
   - Outputs: `f` sampled on cross-sections every 2 m along the shore, every 0.025 s during the
     break, in the format `measure.py` already reads (`case.json` plus one `.bin` per frame);
     and the surface's facets (`output_facets`) every 0.1 s for the viewer.
   - One core on the laptop under the usual guard. Basilisk runs on many cores with MPI and the
     same source (`qcc -D_MPI=1`), for a rented machine if tier B needs one.
6. **`tools/box_tubes.py`** reads each cross-section's frames, keeps the free surface and any
   air it encloses, and measures them as `cases/v1-steep/measure.py` does (the same code, moved
   to `tools/`). The output is the tube per metre of shore, plus two measures only 3D has:
   **the peel speed**, how fast the landing point moves along the shore, and **the tube's
   length**, how far along the shore it stays closed at once.

### In the viewer

7. **`viewer/box.html`** plays the surface frames over the terrain, using the same camera and
   colouring as the side view. This needs a compact mesh format, `wavesim-surface` version 1,
   with indexed triangles per frame and decimated: Basilisk's facets are one polygon per surface
   cell, so a frame is hundreds of thousands of them. This is last in the list. The numbers come first.

## Validation, in order: each step must pass before the next

| Step | Case | Compare with | Passes when |
|---|---|---|---|
| **V0** | Basilisk built from source in Docker, then a 3D strip of V1 | — | it runs on one core, and the 3D cost per cell-step and the cells a break needs are measured. That sets the cost table above. (For OpenFOAM, V0 was `damBreak`: 1.7 µs per cell-step in 2D.) |
| **V1** | 2D (one cell thick): solitary wave on the `steep` seabed | `waveslice`: face vertical 15.97 s at 222.3 m, crest 3.27 m, lip 10.72 m, tube 1.61 × 3.14 m | the face goes vertical within 2 m, and the crest is within 5%. The tube is within 20% at some cell size. **That cell size becomes the 3D box's cell size.** Passed in Basilisk at 0.063 m; OpenFOAM failed at every size |
| **V2** | 2D: plane beach 1 in 10, 2 m and 12 s, inflow from a `wavesim run --box` on that beach | the slice table: tube 1.8 × 2.8 m | the same tolerances. This tests the export and the velocity profile |
| **V3** | 3D box over the same plane beach, uniform along the shore | V2 | the tube is the same on every cross-section, and the same as V2 within 10%. This tests the sides and the 3D mesh |
| **V4** | 3D, the same beach with the swell 20° oblique | the 2D run's own peel angle | the peel measured from landing points matches the angle the 2D run gives where waves break, within a few degrees |
| **P1** | Pipeline box, tier A | the slices at 481.5, 511.5 (50 m strip) and 541.5: tubes 0.6 to 0.7 m | — this is the answer |
| **P2** | Pipeline box, tier B | P1 | a tube the same size within 20% between tiers, or else the result is not converged and the report says so |

V1 to V3 are cheap: 2D cases take minutes to two hours on one core, and V3 can use a short box. They fit on the laptop
under the usual guard, a time limit plus thermal readings. One difference from wavesim's runs:
inside OrbStack's VM, a run cannot be pinned to the efficiency cores. The VM should be limited
to 4 CPUs, and the thermal readings taken as before.

### V1 result, 2026-11-04: OpenFOAM converges, but on the wrong tube

`cases/v1-steep/`: `make_case.py` builds the case, `run.sh` runs it on one core under the guard,
and `measure.py` measures it as `wavesim slice` does. The thermal state stayed nominal
throughout all six runs. `interFoam` cost 1.7 µs per cell-step on one core in 2D.

Four runs went from the solitary wave at 0 s. Two started from the 0.1 m
run at 14.5 s, when it still lay on waveslice's surface. One used `interIsoFoam` (isoAdvector,
plicRDF), which keeps the surface within one cell where `interFoam`'s MULES spreads it over
two. The other used 0.05 m cells, in a tank cut at x = 170 m. The cut keeps all the water
(246.20 m² either way). What the new wall sends shoreward cannot reach the break before about
21 s.

| | waveslice | 0.4 m | 0.2 m | 0.2 m, Co 0.05 | 0.1 m | 0.1 m, isoAdvector | 0.05 m |
|---|---|---|---|---|---|---|---|
| Cells, wall time | — | 22,575, 2 min | 90,300, 14 min | 90,300, 1 h 59 min | 361,200, 1 h 46 min | from 14.5 s, 42 min | 492,800, from 14.5 s, 1 h 51 min |
| Face vertical | 15.97 s, 222.26 m | 15.75 s, 220.56 m | 15.88 s, 221.48 m | 15.90 s, 221.56 m | 15.82 s, 221.14 m | 15.78 s, 220.70 m | 15.75 s, 220.46 m |
| Crest then | 3.27 m | 3.03 m | 3.14 m | 3.17 m | 3.14 m | 3.14 m | 3.15 m |
| Lip lands | 16.93 s | 16.80 s | 16.80 s | 16.77 s | 16.62 s | 16.55 s | 16.52 s |
| Lip thrown | 10.72 m | 11.13 m | 9.97 m | 9.27 m | 8.22 m | 7.77 m | 7.67 m |
| Tube, height × width | 1.61 × 3.14 m | 0.98 × 3.25 m | 1.49 × 2.80 m | 1.50 × 2.42 m | 1.75 × 1.85 m | 1.64 × 1.53 m | 1.54 × 1.41 m |
| Tube area | 3.16 m² | 0.97 m² | 1.60 m² | 1.68 m² | 1.37 m² | 1.20 m² | 1.14 m² |

**Against the gate:** the face passes at every size, within 0.8–1.8 m of waveslice's. The crest
passes from 0.2 m down: it is 4% low. The tube fails at every size. OpenFOAM's tube is as tall
as waveslice's, but half as wide, with a third of the area, and its lip is thrown 7.7 m rather
than 10.7.

**OpenFOAM has converged in cell size.** From 0.1 to 0.05 m the area changes by 17%, the
height by 12% and the throw by 7%. A sharper surface on the same 0.1 m cells moves it the same
way, and by less than the cells do. Frame by frame, the three finest runs lie on one another
until the lip lands, and differ only in the splash after that. Both codes agree to 15.6 s.
After that, waveslice's crest goes on rising, to 3.3 m, and throws the lip from the top.
OpenFOAM's crest stops at 3.15 m, and its jet leaves the face lower down and falls more
steeply.

**The measurements side with waveslice.** The overturn has been measured at field scale, on the
wave in the video that started this: Feddersen, Fincham, Pieszka and others put lidar on 22
Surf Ranch waves in 2019 (Feddersen et al. 2023). They measured the overturn area at the moment
the jet lands, as A/H², with H the wave height just before breaking. It came out at 0.2 to 0.42,
smaller with onshore wind and larger with offshore wind. Their 2-phase Basilisk model of the
same wave with no wind (Feddersen et al. 2024) gives 0.35. Potential flow gives 0.38
(Pick & Feddersen 2026). That model is on a slope of 0.0693, almost our 1 in 15. Its grid had
converged at Δx/h₀ = 6 × 10⁻³, which for our 8 m depth is 0.05 m. Earlier, a VOF model without
air had matched the BEM of Grilli et al. (1997) (Guignard et al. 2001).

Ours, with H the crest height at the curl: waveslice 3.16 / 3.27² = **0.30**, inside the
measured range. OpenFOAM at 0.05 m: 1.14 / 3.15² = **0.12**, below all of them.

**Not the time step.** Every run above but one used `maxCo 0.5` with the damBreak tutorial's
settings. Larsen, Fuhrman & Roenby (2019) showed that this makes interFoam's crest velocities too
high, and its air move fast just above the surface. They recommend Co ≈ 0.05. Ours showed both
signs: the face went vertical 0.2 s early and up to 1.8 m offshore, and the air within 1 m above
the crest moved onshore at 15–21 m/s, about twice the wave's speed of 10 m/s, from as early as
10 s. Feddersen's measurements show that an onshore wind shrinks the overturn.

So the 0.2 m case was run again from 0 s at Co 0.05 (2026-11-05, 1 h 59 min on one core). The tube
did not move: 1.68 m² against 1.60 m² at Co 0.5, so A/H² 0.17 against 0.16, and waveslice's 0.30.
The face goes vertical at the same time and place, within 0.1 m, and the jet is thrown 0.7 m
shorter. Frame by frame, the two runs and the 0.05 m run lie on one another until the lip lands.
The fast air stays too: its mean onshore speed within 1 m of the crest falls from 7–11 m/s to
6–9 m/s, and it still peaks at 16–22 m/s.

Six OpenFOAM runs now agree with one another: cells from 0.4 to 0.05 m, MULES and isoAdvector,
Co 0.5 and 0.05. None of them comes near waveslice or the measurements. Whatever lowers the
tube is in the method itself. The likeliest candidate is how interFoam couples the air and the
water across the surface, where the density jumps by a factor of 1000: the fast air layer is a
sign of that, but V1 has not shown that it is the cause. Basilisk can carry momentum with the
volume fraction itself (`navier-stokes/conserving.h`), which is made for large density ratios,
and Feddersen et al. (2024) reproduced the Surf Ranch with it as a code.

### V1 in Basilisk, 2026-11-05: it passes, and it gives waveslice's tube

`cases/v1-basilisk/`: `v1.c` is the case, `run.sh` compiles and runs it on one core under the
same guard, and `measure.py` reads its frames. Basilisk is built from its source tarball in
`out/basilisk/`, inside the OpenFOAM image, which has the compiler. The set-up follows
Feddersen et al. (2024):
- two-phase Navier-Stokes with momentum-conserving VOF (`navier-stokes/conserving.h`);
- reduced gravity;
- the bed as an embedded boundary;
- a quadtree refined on the surface, the bed and the velocity, with no cell larger than 0.5 m.

Same tank, wave, fluids and free-slip walls as the OpenFOAM case. The top is at +8.1 m rather
than +6 m, because Basilisk cuts the square domain only along cell edges. Water is conserved to
0.01%, and to 14.5 s the wave lies on OpenFOAM's and waveslice's (crest 3.0 m at 206.5 m). The
thermal state stayed nominal throughout.

| | waveslice | Basilisk 0.25 m | 0.126 m | 0.063 m | OpenFOAM 0.05 m |
|---|---|---|---|---|---|
| Cells (most at once), wall time | — | 23,845, 1.5 min | 40,129, 4 min | 76,261, 18 min | 492,800, from 14.5 s, 1 h 51 min |
| Face vertical | 15.97 s, 222.26 m | 16.23 s, 224.56 m | 16.10 s, 223.48 m | 16.10 s, 223.51 m | 15.75 s, 220.46 m |
| Crest then | 3.27 m | 3.31 m | 3.23 m | 3.22 m | 3.15 m |
| Lip lands | 16.93 s | 17.07 s | 16.98 s | 16.95 s | 16.52 s |
| Lip thrown | 10.72 m | 9.25 m | 9.75 m | 9.35 m | 7.67 m |
| Tube, height × width | 1.61 × 3.14 m | 1.21 × 2.84 m | 1.42 × 3.00 m | 1.49 × 2.83 m | 1.54 × 1.41 m |
| Tube area | 3.16 m² | 1.89 m² | 2.42 m² | 2.55 m² | 1.14 m² |
| A/H² | 0.30 | 0.17 | 0.23 | 0.25 | 0.12 |

**Against the gate, at 0.063 m:** the face goes vertical 1.25 m shoreward of waveslice's (within
2 m), the crest is 1.7% low (within 5%), and the tube is 19% small (within 20%). V1 passes. At
0.126 m the tube is 23% small, just outside the gate.

**It has converged.** From 0.126 to 0.063 m the area grows by 5%, after 28% the step before, so
finer cells would add about 0.05 m² more (*Richardson estimate*). Frame by frame, the two finest
runs lie on each other, and both lie on waveslice until the lip lands, the overturn included. The
only visible difference is the last metre of the lip, which in Basilisk is thinner and falls a
little sooner, 9.35 m out where waveslice's goes 10.72 m. That is where most of the 19% is. A
VOF lip a few cells thick can break up or fall, while a potential-flow lip has no thickness
limit, so some of the gap is expected.

**The air stays calm.** Within 1 m above the crest, Basilisk's air moves onshore at -0.5 to
+1.2 m/s on average from 14.5 to 16 s, with 95% of it below 4-9 m/s. OpenFOAM's moved at
7-11 m/s on average, at either Courant number. The same air-water coupling that Basilisk's
conserving scheme is built for is the likeliest reason for OpenFOAM's small tube: an onshore
layer of air that fast acts as a strong onshore wind, and Feddersen measured that onshore wind
shrinks the overturn.

**Cost:** 4.5 µs per cell-step on one core, against interFoam's 1.7. But the quadtree needed a
seventh of OpenFOAM's cells at a finer size, so the whole tank from 0 s took 18 min, where
OpenFOAM's last 4 s at 0.05 m took 1 h 51 min.

**For the 3D box:** Basilisk replaces OpenFOAM, at 0.063 m on the surface (258 m / 2^12, the
level V1 passed at). The time step stays at Basilisk's own CFL 0.5. The box needs one thing V1
did not: the inflow and the sides fed from a 2D run, in place of a solitary wave in a closed
tank (V2).

### V0 in Basilisk, 2026-11-06: a 3D strip of V1

`cases/v1-basilisk/v1-3d.c`, run by `run.sh 10 3600 --3d`: the V1 tank 4.03 m wide between
free-slip walls, started at 14.5 s from the 2D run's state extruded across, and run through the
break to 17.2 s. It took 46 min on one core, 406 steps, at most 419,000 cells and 1.2 GB. The
thermal state stayed nominal.

**It breaks as the 2D run does.** On the middle cross-section, at 0.25 m: the face goes
vertical at 16.23 s and 224.55 m (2D: 16.23 s, 224.56 m), the crest is 3.30 m (3.31), and the lip
lands at 17.07 s (17.07). The tube is 2.14 m² against 1.89 m², 13% larger. The 1 m cells in the
water (below) account for a few per cent of that in 2D. The rest is within what 0.25 m resolves
of a lip: there the 2D tube was still changing by 28% per level.

**What it costs:** 16.6 µs per cell-step, against 4.5 in 2D. Of the 419,000 cells, 172,000 lay
below +8.125 m, the top of the 2D tank: 165 per m² of sea surface. The cost table above is built
from these numbers.

**A Basilisk limitation found on the way:** the embedded bed crashes (a floating-point
exception in the first mesh adaptation) when the domain is masked in two directions, here the
top and the strip's width. Either mask alone works. So the strip has the air up to the octree's
top at +250 m, and no cell larger than 1 m, rather than 0.5 m, to keep that air affordable. In
2D at 0.126 m that change moved the tube from 2.42 to 2.48 m². The real box needs a top and two
sides, so this has to be solved first: report it upstream, find which of Basilisk's restriction
steps fails, or bound the box by the embedded solid instead of masks.

References for this section:

- Feddersen, F., Fincham, A. M., Brodie, K. L., Young, A. P., Spydell, M. S., Grimes, D. J.,
  Pieszka, M. & Hanson, K. (2023). Cross-shore wind-induced changes to field-scale overturning
  wave shape. *J. Fluid Mech.* 958, A4. doi:10.1017/jfm.2023.40
- Feddersen, F., Hanson, K., Mostert, W. & Fincham, A. (2024). Modelling wind-induced changes to
  overturning wave shape. *J. Fluid Mech.* 1000, A48. doi:10.1017/jfm.2024.1012
- Pick & Feddersen (2026). Scaling the shape of shoaling and overturning solitary waves.
  *J. Fluid Mech.* 1040, A8.
- Larsen, B. E., Fuhrman, D. R. & Roenby, J. (2019). Performance of interFoam on the simulation
  of progressive waves. *Coastal Engineering Journal* 61(3), 380–400. arXiv:1804.01158
- Grilli, S. T., Svendsen, I. A. & Subramanya, R. (1997). Breaking criterion and characteristics
  for solitary waves on slopes. *J. Waterway, Port, Coastal, Ocean Eng.* 123(3), 102–112.

## Risks and open choices

- **Reflection at the inflow.** The plan already uses a relaxation zone, which absorbs what the
  beach sends back. The reflected fraction should still be measured in V2, and the zone made
  longer along `x` if it matters.
- **Masks and the embedded bed.** Basilisk's embedded bed fails in a domain masked in two
  directions (V0). The box needs that, or another way to bound it.
- **The sides.** Feeding the sides from the 2D run is consistent seaward of the break. Shoreward,
  the 2D run's bores will differ from the 3D overturning. The plan keeps 40 m between each side
  and the lines measured. An alternative is slip walls on the sides, which mirror the bed as the
  2D run's walls do. V3 and V4 should decide.
- **The bed outside the lidar.** If the lidar has gaps inside the box, the patch mixes 1 m and 3 m
  data. `fetch_patch.py` records where, and the measurements report it per cross-section.
- **The resolution needed for the tube.** V1 passed at 0.063 m. That puts even tier A beyond the
  laptop (see the cost table), so the box is a rented-machine job unless the refinement can be
  limited to where the wave breaks.
- **Time to first answer.** V0 to P1 is roughly three weeks of evenings, most of it in items 1, 5
  and 6 (*estimate*). The viewer can wait.

## What it will not do

- It will not replace the 2D run. The box is a few hundred metres across and costs hours. The
  2D run covers the reef and the whole shelf in minutes.
- It will not make the swell input or the bed more certain. A 3D answer on a 3 m bed, or for a
  swell no buoy saw, is still only as good as those inputs.
- It will not model foam or air bubbles beyond what VOF resolves. Air entrained after landing is
  outside what the comparison uses.

## References

- Grilli, S. T., Guyenne, P. & Dias, F. (2001). A fully non-linear model for three-dimensional
  overturning waves over an arbitrary bottom. *Int. J. Numer. Meth. Fluids* 35.
- Larsen, B. E. & Fuhrman, D. R. (2018). On the over-production of turbulence beneath surface waves
  in Reynolds-averaged Navier–Stokes models. *J. Fluid Mech.* 853.
- Mostert, W., Popinet, S. & Deike, L. (2022). High-resolution direct simulation of deep water
  breaking waves: transition to turbulence, bubbles and droplets production. *J. Fluid Mech.* 942.
- Fincham, A. & Pieszka, M. CFD of the Surf Ranch wave, *Proc. ICCE*:
  https://icce-ojs-tamu.tdl.org/icce/article/view/12736
- OpenFOAM (ESI) Docker images: https://hub.docker.com/r/opencfd/openfoam-default/tags
