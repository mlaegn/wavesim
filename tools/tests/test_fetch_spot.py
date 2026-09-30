"""The rotated resampling must put +x along the bearing and +y to its left.

A plane with a known slope direction is written to an in-memory GeoTIFF; the sampled
bed must then slope only along the intended axis, by exactly the expected amount.
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest
import rasterio
from pyproj import Transformer
from rasterio.io import MemoryFile
from rasterio.transform import from_origin

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fetch_spot  # noqa: E402

LAT, LON, UTM = 21.665, -158.051, "EPSG:32604"
RES = 3.0864e-5  # degrees, about 3.2 m
SIZE = 600  # about 1.9 km, comfortably larger than the test domain


def plane_dataset(bearing_of_slope: float, slope: float):
    """A GeoTIFF whose elevation rises `slope` m/m along the given compass bearing."""
    west, north = LON - SIZE * RES / 2, LAT + SIZE * RES / 2
    transform = from_origin(west, north, RES, RES)
    cols, rows = np.meshgrid(np.arange(SIZE) + 0.5, np.arange(SIZE) + 0.5)
    lon = west + cols * RES
    lat = north - rows * RES
    to_utm = Transformer.from_crs("EPSG:4326", UTM, always_xy=True)
    e, n = to_utm.transform(lon, lat)
    e0, n0 = to_utm.transform(LON, LAT)
    b = math.radians(bearing_of_slope)
    along = (e - e0) * math.sin(b) + (n - n0) * math.cos(b)
    z = (-10.0 + slope * along).astype("float32")

    mem = MemoryFile()
    with mem.open(
        driver="GTiff", height=SIZE, width=SIZE, count=1, dtype="float32",
        crs="EPSG:4326", transform=transform, nodata=-9999.0,
    ) as dst:
        dst.write(z, 1)
    return mem


def sample(mem, bearing):
    with mem.open() as ds:
        return fetch_spot.sample_bed(
            ds, lat=LAT, lon=LON, utm_crs=UTM, bearing=bearing,
            offshore_m=300.0, onshore_m=100.0, width_m=300.0, cell_m=3.0,
        )


@pytest.mark.parametrize("bearing", [0.0, 45.0, 130.0, 250.0])
def test_plus_x_points_along_the_bearing(bearing):
    mem = plane_dataset(bearing_of_slope=bearing, slope=0.05)
    z, info = sample(mem, bearing)
    assert z.shape == (info["ny"], info["nx"]) == (100, 133)
    dz_dx = np.diff(z, axis=1) / 3.0
    dz_dy = np.diff(z, axis=0) / 3.0
    assert np.allclose(dz_dx, 0.05, atol=2e-4), "slope along +x should be the full 0.05"
    assert np.allclose(dz_dy, 0.0, atol=2e-4), "no slope across the bearing"
    # The anchor point sits offshore_m along +x, at the vertical centre: elevation -10 there.
    assert z[50, 100] == pytest.approx(-10.0, abs=0.1)


@pytest.mark.parametrize("bearing", [0.0, 130.0])
def test_plus_y_is_ninety_degrees_counter_clockwise(bearing):
    # A plane rising towards the compass direction (bearing - 90 degrees) rises along +y.
    left = (bearing - 90.0) % 360.0
    mem = plane_dataset(bearing_of_slope=left, slope=0.03)
    z, _ = sample(mem, bearing)
    assert np.allclose(np.diff(z, axis=0) / 3.0, 0.03, atol=2e-4)
    assert np.allclose(np.diff(z, axis=1) / 3.0, 0.0, atol=2e-4)


def test_clipping_flattens_and_counts():
    mem = plane_dataset(bearing_of_slope=130.0, slope=0.05)
    with mem.open() as ds:
        z, info = fetch_spot.sample_bed(
            ds, lat=LAT, lon=LON, utm_crs=UTM, bearing=130.0,
            offshore_m=300.0, onshore_m=100.0, width_m=300.0, cell_m=3.0,
            clip_above=-8.0, clip_below=-15.0,
        )
    assert z.max() == pytest.approx(-8.0)
    assert z.min() == pytest.approx(-15.0)
    assert info["cells_clipped_above"] > 0 and info["cells_clipped_below"] > 0


def test_domain_outside_the_raster_is_an_error():
    mem = plane_dataset(bearing_of_slope=130.0, slope=0.05)
    with mem.open() as ds, pytest.raises(ValueError, match="outside"):
        fetch_spot.sample_bed(
            ds, lat=LAT, lon=LON, utm_crs=UTM, bearing=130.0,
            offshore_m=5000.0, onshore_m=100.0, width_m=300.0, cell_m=3.0,
        )
