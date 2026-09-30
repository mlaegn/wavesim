#!/usr/bin/env python3
"""Crop a real seabed from a remote GeoTIFF into a wavesim bed file.

    uv run --with rasterio --with numpy --with scipy --with pyproj \
        tools/fetch_spot.py pipeline

Reads only the window it needs (HTTP range requests on a tiled GeoTIFF), resamples it
onto a metric grid rotated so +x points along the wave direction, and writes
`data/<name>.json` (header) and `data/<name>.f32` (float32 little-endian elevation in
metres relative to mean sea level, row-major, y outer, x inner, y increasing).
See docs/data-format.md.
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import sys
import tomllib
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.windows import Window
from scipy.ndimage import map_coordinates

FORMAT = "wavesim-bed"
VERSION = 1


def frame_vectors(bearing_deg: float):
    """Unit vectors (east, north) of the grid's +x and +y axes."""
    b = math.radians(bearing_deg)
    ex = (math.sin(b), math.cos(b))  # +x points along the compass bearing
    ey = (-math.cos(b), math.sin(b))  # +y is 90 degrees counter-clockwise from +x
    return ex, ey


def sample_bed(
    ds,
    *,
    lat: float,
    lon: float,
    utm_crs: str,
    bearing: float,
    offshore_m: float,
    onshore_m: float,
    width_m: float,
    cell_m: float,
    clip_above: float = math.inf,
    clip_below: float = -math.inf,
):
    """Resample the raster `ds` (lon/lat degrees) onto the rotated metric grid.

    Returns `(z, info)` with `z` of shape (ny, nx), row j = y index.
    """
    if abs(ds.transform.a) > 0.01:
        raise ValueError("expected a raster in geographic (lon/lat) coordinates")

    to_utm = Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True)
    to_ll = Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True)
    e0, n0 = to_utm.transform(lon, lat)
    ex, ey = frame_vectors(bearing)

    nx = round((offshore_m + onshore_m) / cell_m)
    ny = round(width_m / cell_m)
    # Lower-left corner of the domain, i.e. the point (x, y) = (0, 0).
    origin_e = e0 - offshore_m * ex[0] - (width_m / 2) * ey[0]
    origin_n = n0 - offshore_m * ex[1] - (width_m / 2) * ey[1]

    x = (np.arange(nx) + 0.5) * cell_m
    y = (np.arange(ny) + 0.5) * cell_m
    xx, yy = np.meshgrid(x, y)
    easting = origin_e + xx * ex[0] + yy * ey[0]
    northing = origin_n + xx * ex[1] + yy * ey[1]
    lons, lats = to_ll.transform(easting, northing)

    # Pixel coordinates of every sample, then read just the window that covers them.
    cols, rows = (~ds.transform) @ (lons, lats)
    c0 = max(int(math.floor(cols.min())) - 2, 0)
    r0 = max(int(math.floor(rows.min())) - 2, 0)
    c1 = min(int(math.ceil(cols.max())) + 3, ds.width)
    r1 = min(int(math.ceil(rows.max())) + 3, ds.height)
    if cols.min() < 0 or rows.min() < 0 or cols.max() > ds.width or rows.max() > ds.height:
        raise ValueError("the requested domain reaches outside this raster")
    window = Window(c0, r0, c1 - c0, r1 - r0)
    data = ds.read(1, window=window, masked=True).astype("float64").filled(np.nan)

    # A pixel centre sits half a pixel in from the corner the affine refers to.
    z = map_coordinates(
        data, [rows - 0.5 - r0, cols - 0.5 - c0], order=1, mode="nearest", cval=np.nan
    )
    missing = int(np.isnan(z).sum())
    if missing:
        raise ValueError(f"{missing} of {z.size} cells have no data in the source raster")

    above = int((z > clip_above).sum())
    below = int((z < clip_below).sum())
    z = np.clip(z, clip_below, clip_above)

    info = {
        "nx": nx,
        "ny": ny,
        "origin_easting": origin_e,
        "origin_northing": origin_n,
        "read_window_pixels": [window.width, window.height],
        "cells_clipped_above": above,
        "cells_clipped_below": below,
        "min": float(z.min()),
        "max": float(z.max()),
    }
    return z, info


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name", help="spot name in spots.toml")
    ap.add_argument("--spots", default="spots.toml", type=Path)
    ap.add_argument("--out", default="data", type=Path)
    ap.add_argument("--cell", type=float, help="override cell size in metres")
    args = ap.parse_args(argv)

    spots = tomllib.loads(args.spots.read_text())
    if args.name not in spots:
        print(f"unknown spot {args.name!r}; known: {', '.join(spots)}", file=sys.stderr)
        return 2
    s = dict(spots[args.name])
    if args.cell:
        s["cell_m"] = args.cell

    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"):
        with rasterio.open("/vsicurl/" + s["url"]) as ds:
            z, info = sample_bed(
                ds,
                lat=s["lat"],
                lon=s["lon"],
                utm_crs=s["utm_crs"],
                bearing=s["bearing"],
                offshore_m=s["offshore_m"],
                onshore_m=s["onshore_m"],
                width_m=s["width_m"],
                cell_m=s["cell_m"],
                clip_above=s.get("clip_above", math.inf),
                clip_below=s.get("clip_below", -math.inf),
            )

    stem = args.name if not args.cell else f"{args.name}_{s['cell_m']:g}m"
    args.out.mkdir(parents=True, exist_ok=True)
    data_path = args.out / f"{stem}.f32"
    data_path.write_bytes(z.astype("<f4").tobytes())
    header = {
        "format": FORMAT,
        "version": VERSION,
        "name": args.name,
        "description": s.get("description", ""),
        "nx": info["nx"],
        "ny": info["ny"],
        "dx": s["cell_m"],
        "dy": s["cell_m"],
        "data": data_path.name,
        "dtype": "f32-le",
        "layout": "row-major, y outer (increasing), x inner (increasing)",
        "units": "metres above mean sea level; negative is underwater",
        "frame": {
            "crs": s["utm_crs"],
            "origin_easting": info["origin_easting"],
            "origin_northing": info["origin_northing"],
            "x_bearing_deg": s["bearing"],
            "note": "+x points along the bearing (the direction waves travel, towards the "
            "beach); +y is 90 degrees counter-clockwise from +x; the origin is the "
            "domain's lower-left corner",
        },
        "stats": {k: info[k] for k in ("min", "max", "cells_clipped_above", "cells_clipped_below")},
        "source": {
            "name": s.get("source", ""),
            "url": s["url"],
            "datum": s.get("datum", ""),
            "attribution": s.get("attribution", ""),
            "retrieved": datetime.date.today().isoformat(),
            "read_window_pixels": info["read_window_pixels"],
        },
    }
    header_path = args.out / f"{stem}.json"
    header_path.write_text(json.dumps(header, indent=2) + "\n")
    print(
        f"{header_path}: {info['nx']} x {info['ny']} cells of {s['cell_m']} m, "
        f"elevation {info['min']:.1f}..{info['max']:.1f} m, "
        f"clipped {info['cells_clipped_above']} above / {info['cells_clipped_below']} below"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
