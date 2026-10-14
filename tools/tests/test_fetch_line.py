"""A line must land on the bed frame exactly where the fetch tool laid that frame out.

The fetch tool puts +x along the bearing and +y 90 degrees counter-clockwise from it, from the
frame's origin; the line tool must reach the same points, and sample a known plane exactly.
"""

import math
import sys
from pathlib import Path

import numpy as np
from rasterio.io import MemoryFile
from rasterio.transform import from_origin

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fetch_line  # noqa: E402
import fetch_spot  # noqa: E402


def test_the_line_follows_the_fetch_tools_frame():
    for bearing in (0.0, 90.0, 130.0, 245.0):
        frame = {"origin_easting": 1000.0, "origin_northing": 2000.0, "x_bearing_deg": bearing}
        xs = np.array([0.0, 10.0, 250.0])
        e, n = fetch_line.line_points(frame, 40.0, xs)
        ex, ey = fetch_spot.frame_vectors(bearing)
        for k, x in enumerate(xs):
            assert math.isclose(e[k], 1000.0 + x * ex[0] + 40.0 * ey[0], abs_tol=1e-9)
            assert math.isclose(n[k], 2000.0 + x * ex[1] + 40.0 * ey[1], abs_tol=1e-9)


def test_a_plane_is_sampled_exactly_and_a_gap_is_not_filled():
    # z = 2 lon + 3 lat over a small lon/lat raster, with one pixel masked.
    res, cols, rows = 0.01, 40, 30
    west, north = -158.2, 21.8
    transform = from_origin(west, north, res, res)
    lon = west + (np.arange(cols) + 0.5) * res
    lat = north - (np.arange(rows) + 0.5) * res
    z = (2.0 * lon[None, :] + 3.0 * lat[:, None]).astype("float32")
    nodata = -9999.0
    z[10, 25] = nodata
    profile = dict(driver="GTiff", width=cols, height=rows, count=1, dtype="float32",
                   crs="EPSG:4326", transform=transform, nodata=nodata)
    with MemoryFile() as mem:
        with mem.open(**profile) as ds:
            ds.write(z, 1)
        with mem.open() as ds:
            plons = np.array([west + 5.3 * res, west + 20.7 * res])
            plats = np.array([north - 4.2 * res, north - 18.9 * res])
            got = fetch_line.sample(ds, plons, plats)
            assert np.allclose(got, 2.0 * plons + 3.0 * plats, atol=1e-4)
            # Next to the masked pixel the sample is missing, not interpolated across it.
            near_gap = fetch_line.sample(ds, np.array([west + 25.2 * res]), np.array([north - 10.4 * res]))
            assert np.isnan(near_gap[0])
