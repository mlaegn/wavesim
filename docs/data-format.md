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

`waves` also records `maker_x_m` (the wave-maker line), `maker_sigma_m` (the source's Gaussian width), `sponge_offshore_cells` and `sponge_side_cells` (the absorbing strips), and `near_field_end_m`, the distance along `x` past which the wave-maker's own bump has died away (six source widths beyond the line). **Seaward of `near_field_end_m` and inside the side strips the surface is the source and the sponges, not sea**; a viewer should start at `near_field_end_m` and trim the sides. The Pipeline viewer does.

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
| `clip_above`, `clip_below` | Elevation limits, in metres |
| `datum`, `source`, `attribution`, `description` | Recorded in the header |
