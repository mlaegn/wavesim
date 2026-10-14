#!/usr/bin/env python3
"""Sample one cross-shore line of a bed file's frame from a finer remote DEM.

    uv run --python 3.12 --with rasterio --with numpy --with scipy --with pyproj \\
        tools/fetch_line.py data/pipeline.json --y 511.5

Takes the line at `y` (metres, in the bed's frame) from `x = 0` to the bed's far end and samples
the DEM along it every `--step` metres, by bilinear interpolation, reading only the window the
line crosses. Also samples the bed file itself along the same line, for comparison. Writes
`data/<bed>_line_<y>.json`: a `wavesim-line` (see docs/data-format.md).

The default DEM is the 2013 USACE NCMP topobathy lidar of Oahu, 1 m (see docs/bathymetry.md),
whose gaps are masked, not filled: a gap is written as null.
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import sys
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.windows import Window
from scipy.ndimage import map_coordinates

FORMAT = "wavesim-line"
VERSION = 1
LIDAR_1M = (
    "https://noaa-nos-coastal-lidar-pds.s3.amazonaws.com/dem/USACE_Oahu_HI_LMSL_DEM_2013_9365/"
    "USACE_Oahu_HI_LMSL_DEM_2013_m9365_EPSG-6322.vrt"
)


def line_points(frame: dict, y: float, xs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Easting and northing of the points `(x, y)` of a bed frame."""
    b = math.radians(frame["x_bearing_deg"])
    ex = (math.sin(b), math.cos(b))
    ey = (-math.cos(b), math.sin(b))
    easting = frame["origin_easting"] + xs * ex[0] + y * ey[0]
    northing = frame["origin_northing"] + xs * ex[1] + y * ey[1]
    return easting, northing


def sample(ds, lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
    """Bilinear samples of the lon/lat raster `ds`, NaN where it has no data."""
    if abs(ds.transform.a) > 0.01:
        raise ValueError("expected a raster in geographic (lon/lat) coordinates")
    cols, rows = (~ds.transform) @ (lons, lats)
    c0 = max(int(math.floor(cols.min())) - 2, 0)
    r0 = max(int(math.floor(rows.min())) - 2, 0)
    c1 = min(int(math.ceil(cols.max())) + 3, ds.width)
    r1 = min(int(math.ceil(rows.max())) + 3, ds.height)
    window = Window(c0, r0, c1 - c0, r1 - r0)
    data = ds.read(1, window=window, masked=True).astype("float64").filled(np.nan)
    # A pixel centre sits half a pixel in from the corner the affine refers to. Bilinear
    # interpolation next to a gap gives NaN, so a gap is never filled in.
    return map_coordinates(data, [rows - 0.5 - r0, cols - 0.5 - c0], order=1, mode="nearest")


def bed_along(header: dict, bed_path: Path, y: float, xs: np.ndarray) -> np.ndarray:
    """The bed file's elevation along the line, straight between its cell centres in x, from
    the row nearest `y`."""
    nx, ny, dx, dy = header["nx"], header["ny"], header["dx"], header["dy"]
    z = np.fromfile(bed_path, dtype="<f4").reshape(ny, nx).astype("float64")
    row = min(max(int(round(y / dy - 0.5)), 0), ny - 1)
    centres = (np.arange(nx) + 0.5) * dx
    return np.interp(xs, centres, z[row])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bed", type=Path, help="bed file header, e.g. data/pipeline.json")
    ap.add_argument("--y", type=float, required=True, help="the line, metres along the shore in the bed's frame")
    ap.add_argument("--step", type=float, default=1.0, help="metres between samples")
    ap.add_argument("--url", default=LIDAR_1M, help="DEM to sample (lon/lat GeoTIFF or VRT)")
    ap.add_argument("--out", type=Path, help="output file (default data/<bed>_line_<y>.json)")
    args = ap.parse_args(argv)

    header = json.loads(args.bed.read_text())
    frame = header["frame"]
    length = header["nx"] * header["dx"]
    if not 0 <= args.y <= header["ny"] * header["dy"]:
        print(f"y = {args.y} is outside the bed (0 to {header['ny'] * header['dy']})", file=sys.stderr)
        return 1
    xs = np.arange(0.0, length + 1e-9, args.step)
    easting, northing = line_points(frame, args.y, xs)
    to_ll = Transformer.from_crs(frame["crs"], "EPSG:4326", always_xy=True)
    lons, lats = to_ll.transform(easting, northing)
    with rasterio.open(args.url) as ds:
        z = sample(ds, lons, lats)
    near = bed_along(header, args.bed.parent / header["data"], args.y, xs)

    missing = int(np.isnan(z).sum())
    wet = (near < -1.0) & (near > -10.0)
    both = wet & ~np.isnan(z)
    out = args.out or args.bed.parent / f"{header['name']}_line_{args.y:g}.json"
    doc = {
        "format": FORMAT,
        "version": VERSION,
        "name": f"{header['name']} line y = {args.y:g} m",
        "bed": args.bed.name,
        "y": args.y,
        "step": args.step,
        "x": [round(float(v), 3) for v in xs],
        "z": [None if math.isnan(v) else round(float(v), 3) for v in z],
        "z_bed": [round(float(v), 3) for v in near],
        "units": "metres above local mean sea level; negative is underwater; x along the bed's frame",
        "source": {
            "url": args.url,
            "retrieved": datetime.date.today().isoformat(),
            "missing": missing,
        },
    }
    out.write_text(json.dumps(doc) + "\n")
    diff = float(np.mean(np.abs(z[both] - near[both]))) if both.any() else float("nan")
    print(
        f"wrote {out}: {len(xs)} points, {missing} without data; "
        f"in 1 to 10 m of water the bed file differs by {diff:.2f} m on average"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
