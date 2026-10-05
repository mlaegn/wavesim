# Data formats

Two formats connect the pieces: a **bed file**, written by the fetch tool and read by the solver, and a **run**, written by the solver and meant for a viewer. Both are a small JSON header next to raw little-endian `float32` data, so anything that can read JSON and bytes can read them (Rust, Python, JavaScript).

Both have `format` and `version` fields. Readers must reject a `format` they do not know and a `version` they do not support.

## Coordinate frame

Grids are metric and rotated so the interesting direction is `+x`:

- `+x` points along a compass bearing (`x_bearing_deg`, degrees clockwise from north). It is the direction the waves travel, towards the beach.
- `+y` is 90° counter-clockwise from `+x`. The frame is right-handed.
- Cell `(i, j)` covers `[i·dx, (i+1)·dx] × [j·dy, (j+1)·dy]`; the origin is the lower-left corner of the domain.
- `frame.origin_easting` / `origin_northing` place that corner in `frame.crs` (a projected CRS in metres), so a viewer can put the grid back on a map.

Because the swell direction is fixed by the frame, a different swell direction means fetching the bed again with a different `bearing`.

## Bed file (`wavesim-bed`, version 1)

`data/<name>.json` plus `data/<name>.f32`.

```json
{
  "format": "wavesim-bed",
  "version": 1,
  "name": "pipeline",
  "description": "Banzai Pipeline, Oahu North Shore, Hawaii",
  "nx": 500, "ny": 400, "dx": 3.0, "dy": 3.0,
  "data": "pipeline.f32",
  "dtype": "f32-le",
  "layout": "row-major, y outer (increasing), x inner (increasing)",
  "units": "metres above mean sea level; negative is underwater",
  "frame": {
    "crs": "EPSG:32604",
    "origin_easting": 0.0, "origin_northing": 0.0,
    "x_bearing_deg": 130.0,
    "note": "..."
  },
  "stats": { "min": -14.6, "max": 10.0, "cells_clipped_above": 23847, "cells_clipped_below": 0 },
  "source": { "name": "...", "url": "...", "datum": "...", "attribution": "...", "retrieved": "2026-09-30" }
}
```

The data file holds exactly `nx · ny` values: element `j·nx + i` is the bed elevation of cell `(i, j)`. There must be no gaps; the reader rejects NaN and infinite values and any file of the wrong length.

`stats.cells_clipped_above` / `below` count cells the fetch tool flattened to `clip_above` / `clip_below` from `spots.toml`. Land higher than the ceiling is flattened because waves never reach it; deep water is flattened to keep the time step large.

## Run (`wavesim-run`, version 1)

A directory written by `wavesim run`:

```text
out/pipeline/
├── run.json     # header
├── bed.f32      # copy of the bed, same layout as a bed file's data
└── frames.f32   # frame_count frames, back to back
```

`run.json` carries `nx`, `ny`, `dx`, `dy`, `frame` (copied from the bed), `source`, `frame_count`, `times` (seconds, one per frame), `fields` and `waves` (the settings used: wave height, period, tide, wave-maker position, sponge widths, friction).

`fields` lists what each frame contains, in order. `wavesim run` writes `["eta", "breaking"]`:

- **`eta`** is the free-surface elevation in metres above mean sea level. Where a cell is dry, `eta` equals the bed elevation, so `depth = max(eta − bed, 0)` and a cell is wet where `eta > bed`.
- **`breaking`** is how close each cell is to breaking, from 0 to 1. It is the switch that turns the dispersive terms off: 0 for smooth water, rising to 1 where the surface is steep or tall for its depth, and 0 where there is too little water (under 5 cm) to carry a wave. A value above about 0.8 means the model treats the wave as breaking there. Mid values (0.3 to 0.6) mean *tall or steep for this depth, not yet breaking*. With `--dispersive false` it is the same criteria applied to plain shallow water.

A frame is the blocks of its fields, one after another, each `nx · ny` values in the same layout as the bed; field `f` of frame `k` starts at value `(k · len(fields) + f) · nx · ny`. A reader must look fields up by name, not assume `eta` comes first.

`waves` also records `maker_x_m` (the wave-maker line), `maker_sigma_m` (the source's Gaussian width), `maker_depth_m` and `maker_depth_min_m` (the depth there along the middle row, and the shallowest along the line), `maker_second_harmonic` (the second harmonic a sinusoid lacks there, relative to the wave), `sponge_offshore_cells` (the absorbing strip behind the maker), `cropped_offshore_m` (how much of the bed file's offshore end the run left out; `nx`, the bed and `frame.origin_*` describe the cropped grid), `side_margin_cells` (rows along each side that a viewer should trim: the sides are walls, which mirror the bed, and the middle half of the width is what to look at), and `near_field_end_m`, the distance along `x` past which the wave-maker's own bump has died away (six source widths beyond the line). **Seaward of `near_field_end_m` the surface is the source and the sponge, not sea**; a viewer should start there. Runs written before the sides were walls record `sponge_side_cells` instead of `side_margin_cells`, the width of absorbing strips along them, which a viewer should trim the same way.

**Statistics.** Frames sample a breaking wave only now and then (frames 2 s apart catch about 40% of the cells that break), so a run also records, at every step over the last four wave periods, `stats.f32`: one block per field in the order `run.json`'s `stats.fields` gives, laid out like a frame. `break_fraction` is the fraction of that time each cell read as breaking (the `breaking` field at 0.8 or more); `eta_max` and `eta_min` are the highest and lowest the surface got there; and `wave_height` is the crest-to-trough height of each wave, period by period, averaged, which unlike `eta_max - eta_min` is not inflated by a slow change of the mean level. `stats.from_s` and `stats.to_s` give the window. A run that stopped before the window has no `stats`.

Runs with `wave_height` also record in `waves`: `break_line`, one entry per row, `[x_m, depth_m, height_m]` where the wave is tallest between the wave-maker's bump (or the relaxation zone) and the shore, in water at least 0.3 m deep, or `null`; `waves_recorded`, the number of waves averaged; and `unsettled_per_period`, how much the per-wave height still changed per period over those waves (the median over the sea in the middle half of the width), which should be under 0.03.

**Nested runs.** `wavesim run` writes the reef, on the bed's own cells, as the run in the output directory, and the coarse run that carried the swell across the shelf, with the wave-maker, as a run of its own in `coarse/` inside it. The reef's run has no wave-maker; its `waves` record `nested_in` (`"coarse"`), `coarse_cell_m`, `fine_start_s` (its first frame is then, not at 0), `edge_depth_m` (the depth where the reef grid begins), `zone_offshore_cells` (the relaxation zone along its offshore edge, which is the coarse run, not this one; `near_field_end_m` is its end, and `sponge_offshore_cells` its width, so a viewer starts past it and marks it as a model zone) and `coarse_breaking_in_zone` (the largest fraction of the time the coarse run spent breaking under that zone, in the middle half of the width, which should be 0). A reader must not assume `times` starts at 0.

## Reading a run

Python:

```python
import json, numpy as np
h = json.load(open("out/pipeline/run.json"))
frames = np.memmap("out/pipeline/frames.f32", dtype="<f4", mode="r",
                   shape=(h["frame_count"], len(h["fields"]), h["ny"], h["nx"]))
eta = frames[:, h["fields"].index("eta")]
breaking = frames[:, h["fields"].index("breaking")]
bed = np.fromfile("out/pipeline/bed.f32", dtype="<f4").reshape(h["ny"], h["nx"])
depth = np.maximum(eta[10] - bed, 0)
```

JavaScript:

```js
const h = await (await fetch("run.json")).json();
const frames = new Float32Array(await (await fetch("frames.f32")).arrayBuffer());
const field = (name, k) => {
  const f = h.fields.indexOf(name), n = h.nx * h.ny;
  return frames.subarray((k * h.fields.length + f) * n, (k * h.fields.length + f + 1) * n);
};
```

`Float32Array` uses the machine's byte order, which is little-endian on every current browser and desktop CPU.

## Spot registry (`spots.toml`)

One table per spot, read by `tools/fetch_spot.py`:

| Key | Meaning |
|---|---|
| `url` | GeoTIFF to crop; must be tiled (range requests) and in lon/lat degrees |
| `lat`, `lon` | Anchor point, normally on the shoreline |
| `utm_crs` | Projected CRS in metres used for the grid |
| `bearing` | Compass bearing `+x` points along |
| `offshore_m`, `onshore_m` | Domain extent seaward and landward of the anchor, along `+x` |
| `width_m` | Domain extent along `y`, centred on the anchor |
| `cell_m` | Cell size |
| `clip_above`, `clip_below` | Elevation limits, in metres (optional) |
| `datum`, `source`, `attribution`, `description` | Recorded in the header |
