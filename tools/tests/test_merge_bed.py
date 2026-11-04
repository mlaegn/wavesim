"""The merge must keep the fine raster where it has data and fill only its gaps from the coarse one.

The coarse raster is a plane, which bilinear interpolation reproduces exactly, so every filled
cell has a known answer. The fine raster is a different plane with a hole in it, so a cell
taken from the wrong source is off by metres.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fetch_spot  # noqa: E402
import merge_bed  # noqa: E402

FINE_RES, COARSE_RES = 0.0001, 0.001  # degrees, about 10 m and 100 m
WEST, NORTH = -9.43, 39.01


def coarse_plane(lon, lat):
    return -20.0 + 1000.0 * (lon - WEST) - 2000.0 * (NORTH - lat)


def write(path, z, res, west, north, nodata=None):
    with rasterio.open(
        path, "w", driver="GTiff", height=z.shape[0], width=z.shape[1], count=1, dtype="float32",
        crs="EPSG:4326", transform=from_origin(west, north, res, res), nodata=nodata,
    ) as dst:
        dst.write(z.astype("float32"), 1)
    return path


@pytest.fixture
def sources(tmp_path):
    # Coarse: 20 x 20 pixels from (WEST, NORTH). Fine: 150 x 150 pixels starting 2.5 coarse
    # pixels in, so the coarse raster covers all of it with room to spare.
    rows, cols = np.mgrid[0:20, 0:20] + 0.5
    coarse = coarse_plane(WEST + cols * COARSE_RES, NORTH - rows * COARSE_RES)
    fw, fn = WEST + 0.0025, NORTH - 0.0025
    fine = np.full((150, 150), -5.0)  # flat, far from the coarse plane
    fine[40:80, 60:120] = np.nan  # a hole the coarse raster must fill
    fine[:, 140:] = np.nan  # and a strip along the east edge
    return (
        write(tmp_path / "fine.tif", fine, FINE_RES, fw, fn, nodata=float("nan")),
        write(tmp_path / "coarse.tif", coarse, COARSE_RES, WEST, NORTH),
        (fw, fn),
    )


def test_fine_data_is_kept_and_only_its_gaps_are_filled(sources):
    fine_path, coarse_path, (fw, fn) = sources
    with rasterio.open(fine_path) as fine, rasterio.open(coarse_path) as coarse:
        z, transform = merge_bed.merge(fine, coarse, shift_m=-2.0)
    assert z.shape == (150, 150), "the coarse raster covers all of the fine one"
    assert transform.almost_equals(from_origin(fw, fn, FINE_RES, FINE_RES))

    rows, cols = np.mgrid[0:150, 0:150] + 0.5
    expected = coarse_plane(fw + cols * FINE_RES, fn - rows * FINE_RES) - 2.0
    gap = np.zeros((150, 150), bool)
    gap[40:80, 60:120] = gap[:, 140:] = True
    assert np.allclose(z[~gap], -7.0), "fine cells keep their value, shifted"
    assert np.allclose(z[gap], expected[gap], atol=1e-3), "gaps are the coarse plane, shifted"


def test_bounds_select_the_fine_pixels_whose_centres_lie_inside(sources):
    fine_path, coarse_path, (fw, fn) = sources
    bounds = (fw + 0.00101, fn - 0.0050, fw + 0.0030, fn - 0.00049)
    with rasterio.open(fine_path) as fine, rasterio.open(coarse_path) as coarse:
        z, transform = merge_bed.merge(fine, coarse, bounds=bounds)
    # Columns 10..29 (centres fw + 0.00105 .. fw + 0.00295), rows 5..49.
    assert z.shape == (45, 20)
    assert transform.c == pytest.approx(fw + 10 * FINE_RES)
    assert transform.f == pytest.approx(fn - 5 * FINE_RES)


def test_a_gap_the_coarse_raster_does_not_cover_is_an_error(sources, tmp_path):
    fine_path, _, (fw, fn) = sources
    small = write(tmp_path / "small.tif", np.zeros((5, 5)), COARSE_RES, fw, fn)
    with rasterio.open(fine_path) as fine, rasterio.open(small) as coarse:
        with pytest.raises(ValueError, match="does not cover"):
            merge_bed.merge(fine, coarse, bounds=tuple(fine.bounds))


def test_a_spot_with_a_merge_table_builds_a_bed_for_any_swell_direction(sources, tmp_path):
    fine_path, coarse_path, (fw, fn) = sources
    merged = tmp_path / "merged.tif"
    spots = tmp_path / "spots.toml"
    spots.write_text(f"""
[here]
url = "{merged}"
lat = {fn - 0.0075}
lon = {fw + 0.0075}
utm_crs = "EPSG:32629"
bearing = 120.0
offshore_m = 300.0
onshore_m = 100.0
width_m = 300.0
cell_m = 3.0

[here.merge]
fine = "{fine_path}"
coarse = "{coarse_path}"
shift_m = -2.0
""")
    assert merge_bed.main(["here", "--spots", str(spots)]) == 0
    out = tmp_path / "beds"
    assert fetch_spot.main(["here", "--spots", str(spots), "--out", str(out), "--from", "290"]) == 0
    header = json.loads((out / "here_290.json").read_text())
    assert header["frame"]["x_bearing_deg"] == 110.0, "a swell from 290 deg travels towards 110 deg"
    assert (header["nx"], header["ny"]) == (133, 100)

    assert fetch_spot.main(["here", "--spots", str(spots), "--out", str(out), "--width", "600"]) == 0
    header = json.loads((out / "here_w600.json").read_text())
    assert header["frame"]["x_bearing_deg"] == 120.0 and header["ny"] == 200


REPO = Path(__file__).resolve().parents[2]
COXOS = REPO / "data" / "coxos" / "coxos.f32"


@pytest.mark.skipif(not COXOS.exists(), reason="needs the Coxos sources and bed in data/, which git does not hold")
def test_the_coxos_bed_is_rebuilt_from_its_two_sources(tmp_path, monkeypatch):
    # The bed the Coxos runs used was cut from a seabed merged by hand on 2026-11-03.
    monkeypatch.chdir(REPO)
    spots = (REPO / "spots.toml").read_text().replace(
        'url = "data/coxos/coxos_seabed_msl.tif"', f'url = "{tmp_path / "merged.tif"}"'
    )
    (tmp_path / "spots.toml").write_text(spots)
    assert merge_bed.main(["coxos", "--spots", str(tmp_path / "spots.toml")]) == 0
    assert fetch_spot.main(["coxos", "--spots", str(tmp_path / "spots.toml"), "--out", str(tmp_path)]) == 0
    rebuilt = np.fromfile(tmp_path / "coxos.f32", "<f4")
    assert np.abs(rebuilt - np.fromfile(COXOS, "<f4")).max() < 1e-5
