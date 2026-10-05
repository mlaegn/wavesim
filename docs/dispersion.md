# Dispersion and breaking

How the solver goes beyond shallow water, what each choice is based on, what was tried and dropped, and where it stops being trustworthy.

## Why shallow water is not enough

In the shallow-water equations every wave in a given depth travels at `sqrt(g h)`, whatever its length. A tall wave's crest then runs ahead of its trough, so the front steepens into a shock **wherever the wave is**: a 2.5 m swell in 8 m of water turns into bores within tens of metres, long before it reaches a reef, and the model has lost over 80% of its energy by 500 m.

Real water disperses: short components lag, which spreads a steep front out. A swell crosses an ocean and only breaks when its height approaches the depth. The solver adds that back with the Serre–Green–Naghdi (SGN) equations, which are fully nonlinear (unlike Boussinesq-type models, they stay sensible when `h` changes a lot under a wave) and have an exact solitary-wave solution to test against.

## The model

For a flat bed, a vertical velocity that grows linearly from the bed (`w = -z ∇·u`) gives a non-hydrostatic pressure, and depth-averaging the horizontal momentum gives

```text
∂(hu)/∂t + ∇·(h u⊗u) + g h ∇η = (1/3) ∇( h³ Γ )
Γ = ∂(∇·u)/∂t + u·∇(∇·u) − (∇·u)²
```

The linear dispersion relation is `ω² = g h k² / (1 + (kh)²/3)`, accurate to 3% in `ω²` at `kh = 1.26` (the full theory is `ω² = g k tanh(kh)`), and the group speed is the phase speed divided by `1 + (kh)²/3`.

**Not included:** the bed-slope terms of the dispersive operator. The operator uses the local total depth as if the bed were flat. The omitted terms are of order `slope × (kh)²`, under 1% on a 1:30 beach for `kh < 0.7`, but they are not zero on the steep ledge of a reef.

## What is solved each step

`Γ` contains the time derivative of the velocity, so with `w = ∂u/∂t` every evaluation of the time derivative needs

```text
h w − ∇( c ∇·w ) = r,        c = φ h³ / 3
r = (shallow-water momentum rate) − u (mass rate) + ∇( c N ),     N = u·∇(∇·u) − (∇·u)²
```

then `∂(hu)/∂t = h w + u ∂h/∂t`. The operator is symmetric and positive definite (checked numerically, with walls and with periodic wrapping), so the system is solved by conjugate gradients with a Jacobi preconditioner, starting from the previous solution. A relative tolerance of 1e-4 gives the same waves as 1e-6 to three digits with half the iterations; it takes 7 to 11 iterations per solve. The solves are evaluated in both Runge–Kutta stages.

The hyperbolic part is the existing finite-volume solver, so mass and momentum conservation, well-balancing over a rough bed and wet/dry handling are unchanged. The mass source of the wave-maker enters `∂h/∂t`, and so the SGN momentum correction accounts for it.

## The breaking switch

`φ` is 1 where the wave is smooth and falls to 0 where it is steep or tall for its depth, which leaves plain shallow water there: the shock then dissipates the wave, which is how breaking is represented. Two criteria, combined by taking the larger:

| Criterion | Fades from | Fully off at |
|---|---|---|
| `η / h` (surface above still water over total depth) | 0.30 | 0.55 |
| `|∇η|` (surface slope) | 0.40 | 0.65 |

`φ` sits inside the gradient, so the operator stays symmetric. Water shallower than 5 cm gets `φ = 0`.

A run counts a cell as breaking where the switch reads 0.8, which the slope criterion reaches at a surface angle of 30°, the onset criterion of the hybrid breaking models this follows (Tissier et al. 2012, Kazolea et al. 2014). It used to fade from a slope of 0.20 (11°) and be off by 24°: on Ting & Kirby's laboratory beach that broke the spilling and the plunging breaker 1.0 and 1.3 m early, and on a 1:30 beach at `H/h` 0.7 where Goda's formula suggests 0.9. With the 30° onset the 1:30 beach breaks at `H/h` 0.81.

**Against the laboratory.** Ting & Kirby (1994) recorded where regular waves broke on a 1:35 beach from a flume 0.40 m deep: a spilling breaker (0.125 m, 2 s) at x = 6.40 m in 0.196 m of water, and a plunging one (0.128 m, 5 s) at x = 7.795 m in 0.156 m. The test makes the waves with a relaxation zone pulling towards a first-order cnoidal wave, tuned once so the incident height matches the measured one, as the laboratory tuned its paddle. On 5 cm cells:

| | Tallest (the break point) | Switch first reads breaking | Laboratory |
|---|---|---|---|
| Spilling | x = 6.76 m, 0.187 m deep, 0.179 m tall | x = 6.91 m | x = 6.40 m, 0.196 m |
| Plunging | x = 7.57 m, 0.164 m deep, 0.192 m tall | x = 7.22 m | x = 7.795 m, 0.156 m |

On 10, 5 and 2.5 cm cells the tallest point is at 6.47, 6.76, 6.63 m (spilling) and 7.57, 7.57, 7.46 m (plunging): it hardly moves. The switch's first reading moves seaward as the cells shrink (7.47, 6.91, 5.48 m and 7.57, 7.22, 4.81 m), because it measures a slope over two cells; on 2.5 cm cells, an eighth of the breaking depth, ripples a few cells long that these equations carry wrongly set it off far too early. So the break point is where the wave is tallest, and the seaward edge of a break map is uncertain by about that much.

**On the Pipeline reef.** One cross-shore line of the Pipeline bed (y = 361.5 m, through the part of the reef where the 2.5 m, 14 s swell breaks well out), run in 1D on 3, 1.5 and 0.75 m cells from the same profile:

| Cells | Tallest (the break point) | Height there | Switch first reads breaking |
|---|---|---|---|
| 3 m | x = 2191.5 m, 4.05 m deep | 4.38 m | x = 2203.5 m |
| 1.5 m | x = 2195.2 m, 3.83 m deep | 4.78 m | x = 2193.8 m |
| 0.75 m | x = 2196.4 m, 3.77 m deep | 4.81 m | x = 980.6 m, in 23 m of water |

The break point moves 3.7 m and then 1.2 m, so where it breaks is settled on 3 m cells. The height at breaking grows 9% and then 0.6%: 3 m cells under-predict it by about 9%, 1.5 m cells are converged. On 0.75 m cells, thirty times finer than the 23 m of water offshore, the ripples of the laboratory flume appear again and set the switch off a kilometre out to sea; that is why a nested run keeps 6 m cells offshore. (In 1D the line has no refraction, so its heights are not the 2D run's.)

## Reconstruction: what was tried

Dispersion exposed a problem that had nothing to do with dispersion: at the resolution of the Pipeline bed (3 m cells), the MC limiter damps a smooth wave noticeably, because it clips every wave crest and trough to first order. Energy height of an 8 s, 1.5 m swell on a 1:30 beach (requested 1.5 m offshore):

| Scheme | Cells | At 7.3 m depth | At 3.5 m | At 1.5 m |
|---|---|---|---|---|
| MC limiter | 3.0 m | 1.36 | 1.12 | 0.64 |
| MC limiter | 1.5 m | 1.50 | 1.52 | 1.02 |
| WENO3 | 3.0 m | 1.21 | 0.83 | 0.46 |
| WENO3 | 1.5 m | 1.47 | 1.40 | 0.89 |
| **Third order, blended with MC by `φ`** | 3.0 m | 1.45 | 1.37 | 0.93 |
| **Third order, blended with MC by `φ`** | 1.5 m | 1.52 | 1.60 | 1.09 |

- **WENO3 was worse than MC** at these resolutions and was removed. Three-cell WENO drifts to its lower-order stencil too readily at 20 to 40 cells per wavelength.
- **Third order** (the unlimited upwind-biased scheme `(2·low + 5·c − high)/6` where `φ` says the wave is smooth, and the MC limiter where it is not) made a 3 m grid behave roughly like a 1.5 m grid with MC. It was the default until the swell had to cross the Pipeline shelf: over that kilometre, 6 m cells lost 14% of the height they kept on 3 m cells, and 18–24% by the reef.
- **`Order::Fifth`** replaced it: the unlimited fifth-order upwind-biased scheme `(2 v₀ − 13 v₁ + 47 v₂ + 27 v₃ − 3 v₄)/60` where the five cells around a face are smooth and wet, the MC limiter where they are not. On 6 m cells across the Pipeline shelf it comes within 1.5% of the third-order run on 3 m cells, and within 3% on the reef; over 1 km of flat 12 m water at 24 cells per wavelength it keeps 100.7% of a 14 s swell's height, where third order keeps 94.2%. It costs no more per step, and it is what lets a nested run carry the swell across the shelf on cells twice as large.
- **SSP-RK3** in place of RK2 made no measurable difference (differences under 0.01 m at both resolutions) and costs half as much again, so it was dropped.

The table is for the schemes before fifth order. Dissipation does not vanish at 3 m cells with third order: the offshore energy height is 1.45 against 1.5 requested.

## Wave-maker calibration, and a mistake that is worth recording

A source of strength `Q` per unit length radiates an amplitude `Q / (2 c_g cos θ)` on each side, where **`c_g` is the group speed**, not the phase speed. The first dispersive version used the phase speed and radiated waves 11% too tall at `kh = 0.59` — exactly `1 + (kh)²/3`. In shallow water the two speeds coincide, which is why the shallow-water tests never showed it. The source also has a finite width, so a Gaussian of width `σ` radiates a wave of wavenumber `k_x` with its strength scaled by `exp(−k_x²σ²/2)`; the strength is boosted by the inverse. That accounts for about 1% of the amplitude deficit the wave-maker had before dispersion was added (1.4–1.6% before, 0.4–0.6% after).

## A failure on the real bed, and the mask that fixed it

The synthetic tests all passed, then the first run on the Pipeline bed blew up at 52 s (88 s on the coarser 6 m grid). The time step had been collapsing for several seconds before, which made the run look merely slow. Stepping the 6 m case and reporting the fastest cell showed a **thin-film runaway** at the edge of the run-up: a 1 to 2 mm film on a bed 0.9 m above sea level, whose speed grew by a factor of about 1.3 every half second until it reached hundreds of kilometres per second.

The dispersive terms were already switched off in water thinner than 5 cm (`c = 0`), but the *linear system was still coupled* across those cells. A film's acceleration is `w = r / h`, which is enormous for `h` of a millimetre; it entered its neighbours' divergence, and the neighbours' dispersive force fed back into the film's equation, again divided by that depth. The system `A w = r` is positive definite and has a bounded solution at every step; the loop was in the time evolution.

The fix removes cells shallower than 5 cm from the dispersive system entirely. They are not unknowns, they contribute zero velocity to the divergence and to the nonlinear term, and their own momentum rate passes through as plain shallow water, as it did before dispersion existed. The restricted operator `M A M` is still symmetric and positive definite (unit-tested with a patchy mask, including that masked cells stay exactly zero). A white-box regression test checks that changing a film's speed from 1 to 1000 m/s changes nothing about its neighbours, and fails when the masking is removed. On the 6 m bed the same case now runs the full 200 s.

## Verification

| Check | Result | Reference |
|---|---|---|
| Standing wave in a closed basin, frequency | 1.4278 rad/s (MC), 1.4261 (third order) | 1.4247 SGN theory; full linear theory 1.4478; shallow water 1.7629 |
| Exact SGN solitary wave, `a/h = 0.2`, 20 s | amplitude 0.2003 / 0.2001, crest 88.65 / 88.63 m, shape error 0.7% / 0.2% (MC / third) | exact 0.2, 88.62 m; shallow water: 0.136, 57% |
| Wave-maker, `kh = 0.59` | amplitude 0.0498 and 0.0499, phase lag 4.737 rad | 0.05; 4.741 (shallow water would give 4.486) |
| 2.5 m, 14 s swell in 8 m, energy kept over 400 m | 93% | shallow water keeps 14% |
| 1:30 beach, 8 s, 1.5 m: shoaling from 6 m to 3.5 m | ×1.14 | Green's law ×1.14 |
| Same beach: breaking | peak `H/h` 0.68 at 3.0 m depth; surf zone `H/h` up to 0.96 in 1.0–2.5 m | textbook range 0.7–1.2 (see the caveat above) |
| Still water over a rough bed, with and without dry bumps | momentum below 1e-9 | zero |
| Linear solves | 4–11 iterations each, none unconverged | |

## Known limits

- **Flat-bed dispersive operator**, as above.
- **`kh` range.** SGN is accurate for `kh` up to about 1 to 2. The wave-maker refuses a wave too short for the model at its depth (`ω² h / g ≥ 3`).
- **Breaking is a switch, not a model.** There is no distinction between spilling and plunging, no roller and no turbulence, and the thresholds are tuned to textbook indices.
- **A sinusoid is not a swell.** At Ursell number `H λ² / h³` above about 25 (a 2.5 m, 14 s swell in 8 m is about 37), a sinusoidal wave-maker produces a wave that evolves into sharp crests and flat troughs, conserving energy but not crest-to-trough height. Compare energy height `2√2 σ` where height matters, not the requested crest-to-trough value.
- **Oblique waves with dispersion are untested.** The periodic-boundary refraction tests are shallow-water only.
- **Cost.** About 78 ms per step on 200,000 cells, against 25 ms for shallow water.
