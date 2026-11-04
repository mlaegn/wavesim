#!/usr/bin/env python3
"""Merge a fine seabed over a coarse one into a single GeoTIFF that fetch_spot.py can crop.

    uv run --with rasterio --with numpy --with scipy tools/merge_bed.py coxos

For spots where no single raster reaches from the surf zone out to the wave-maker, such as
Portugal's: the 2011 coastal LiDAR (~15 m) covers the surf zone only, and the EMODnet DTM
(~115 m) covers the shelf. The output is on the fine raster's grid, over the fine pixels whose
centres fall inside the coarse raster (or inside `bounds`). Each pixel is the fine value where
the fine raster has one, and the coarse raster interpolated bilinearly where it has none; then
`shift_m` moves both onto mean sea level. The spot's `[name.merge]` table in spots.toml names
the two sources and the shift, and the spot's `url` is where the merged file goes.
"""

from __future__ import annotations

import argparse
import math
import sys
import tomllib
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import Affine
from rasterio.windows import Window
from scipy.ndimage import map_coordinates


def merge(fine, coarse, *, shift_m: float = 0.0, bounds=None):
    """The merged elevation over `bounds` (west, south, east, north) and its transform.

    `fine` and `coarse` are open rasters in geographic (lon/lat) coordinates, north up.
    """
    for ds in (fine, coarse):
        if abs(ds.transform.a) > 0.01 or ds.transform.e >= 0:
            raise ValueError("expected north-up rasters in geographic (lon/lat) coordinates")
    if bounds is None:
        bounds = tuple(coarse.bounds)
    west, south, east, north = bounds

    # Fine pixels whose centres lie inside the bounds.
    t = fine.transform
    c0 = max(math.ceil((west - t.c) / t.a - 0.5), 0)
    c1 = min(math.floor((east - t.c) / t.a - 0.5) + 1, fine.width)
    r0 = max(math.ceil((north - t.f) / t.e - 0.5), 0)
    r1 = min(math.floor((south - t.f) / t.e - 0.5) + 1, fine.height)
    if c1 <= c0 or r1 <= r0:
        raise ValueError("the fine raster has no pixels inside the bounds")
    window = Window(c0, r0, c1 - c0, r1 - r0)
    z = fine.read(1, window=window, masked=True).astype("float64").filled(np.nan)
    transform = fine.window_transform(window)

    # The coarse raster at every fine pixel centre, where the fine one has nothing.
    gap = np.isnan(z)
    rows, cols = np.nonzero(gap)
    lon = transform.c + (cols + 0.5) * transform.a
    lat = transform.f + (rows + 0.5) * transform.e
    pc, pr = (~coarse.transform) @ (lon, lat)
    if pc.size and (pc.min() < 0 or pr.min() < 0 or pc.max() > coarse.width or pr.max() > coarse.height):
        raise ValueError("the coarse raster does not cover every gap in the fine one")
    data = coarse.read(1, masked=True).astype("float64").filled(np.nan)
    z[gap] = map_coordinates(data, [pr - 0.5, pc - 0.5], order=1, mode="nearest")
    missing = int(np.isnan(z).sum())
    if missing:
        raise ValueError(f"{missing} of {z.size} cells have no data in either raster")
    return z + shift_m, transform


def write_geotiff(path: Path, z, transform: Affine) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        path, "w", driver="GTiff", height=z.shape[0], width=z.shape[1], count=1,
        dtype="float32", crs="EPSG:4326", transform=transform, compress="deflate",
    ) as dst:
        dst.write(z.astype("float32"), 1)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name", help="spot name in spots.toml, with a [name.merge] table")
    ap.add_argument("--spots", default="spots.toml", type=Path)
    ap.add_argument("--out", type=Path, help="write here instead of the spot's url")
    args = ap.parse_args(argv)

    spots = tomllib.loads(args.spots.read_text())
    if args.name not in spots or "merge" not in spots[args.name]:
        known = [k for k, v in spots.items() if "merge" in v]
        print(f"no [{args.name}.merge] in {args.spots}; spots with one: {', '.join(known)}", file=sys.stderr)
        return 2
    m = spots[args.name]["merge"]
    out = args.out or Path(spots[args.name]["url"])

    with rasterio.open(m["fine"]) as fine, rasterio.open(m["coarse"]) as coarse:
        z, transform = merge(fine, coarse, shift_m=m.get("shift_m", 0.0), bounds=m.get("bounds"))
    write_geotiff(out, z, transform)
    print(f"{out}: {z.shape[1]} x {z.shape[0]} pixels, elevation {z.min():.1f}..{z.max():.1f} m")
    return 0


if __name__ == "__main__":
    sys.exit(main())
