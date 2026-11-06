"""The spot sheet's pieces: map coordinates, takeoffs, the tide comparison, and ring simplifying."""

import math
import sys
from pathlib import Path

import numpy as np
import pytest
from pyproj import Transformer

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import spot_sheet  # noqa: E402

FRAME = {"crs": "EPSG:32629", "origin_easting": 462000.0, "origin_northing": 4317000.0, "x_bearing_deg": 120.0}


def test_frame_puts_x_along_the_bearing_and_y_to_its_left():
    f = spot_sheet.Frame(FRAME)
    to_utm = Transformer.from_crs("EPSG:4326", FRAME["crs"], always_xy=True)
    e, n = to_utm.transform(*f.lonlat(100.0, 0.0))
    assert e - 462000.0 == pytest.approx(100 * math.sin(math.radians(120)), abs=1e-3)
    assert n - 4317000.0 == pytest.approx(100 * math.cos(math.radians(120)), abs=1e-3)
    e, n = to_utm.transform(*f.lonlat(0.0, 100.0))  # +y is 90 deg counter-clockwise: bearing 30
    assert (e - 462000.0, n - 4317000.0) == pytest.approx((50.0, 86.6025), abs=1e-3)


def rows_for(ts):
    return [{"rock": False} for _ in ts]


def test_a_reef_that_breaks_before_both_sides_is_a_takeoff():
    y = np.arange(200) * 3.0
    ts = 80 + 0.02 * y - 6 * np.exp(-(((y - 300) / 40) ** 2))  # a dip at y = 300 on a sloping coast
    found = spot_sheet.takeoffs(rows_for(ts), ts, 3.0)
    assert len(found) == 1 and abs(found[0]["k"] - 100) <= 2  # the slope moves the bottom a row
    assert found[0]["stands_out_s"] > 2


def test_a_coast_that_runs_seaward_to_one_end_has_no_takeoff():
    ts = 80 - 0.02 * np.arange(200) * 3.0  # breaks ever earlier towards the end of the band
    assert spot_sheet.takeoffs(rows_for(ts), ts, 3.0) == []


def test_a_row_onto_the_rock_is_never_a_takeoff():
    y = np.arange(200) * 3.0
    ts = 80 - 6 * np.exp(-(((y - 300) / 40) ** 2))
    rows = rows_for(ts)
    for k in range(90, 111):
        rows[k]["rock"] = True
    assert all(abs(p["k"] - 100) > 10 for p in spot_sheet.takeoffs(rows, ts, 3.0))


def test_the_tide_shift_compares_positions_on_the_bed_not_in_each_run():
    # The low-tide run is cropped 45 m further offshore, so its x reads 45 m more for the same place.
    mid = {"_rows": [{"j": j, "x_bed": 3300.0} for j in range(10)]}
    low = {"_rows": [{"j": j, "x_bed": 3270.0} for j in range(10)]}
    assert spot_sheet.tide_shift(mid, low) == 30.0


def test_simplify_keeps_a_squares_corners_only():
    side = np.linspace(0, 30, 11)
    ring = np.concatenate([
        np.c_[side, 0 * side], np.c_[0 * side + 30, side], np.c_[side[::-1], 0 * side + 30], np.c_[0 * side, side[::-1]],
    ])
    r = spot_sheet.simplify(ring, 3.0)
    assert {tuple(p) for p in r} == {(0.0, 0.0), (30.0, 0.0), (30.0, 30.0), (0.0, 30.0)}


def test_inside_reads_the_area_as_lat_lon_corners():
    area = [[39.0, -9.43], [39.0, -9.42], [39.01, -9.42], [39.01, -9.43]]
    lon = np.array([-9.425, -9.435, -9.425])
    lat = np.array([39.005, 39.005, 39.02])
    assert spot_sheet.inside(lon, lat, area).tolist() == [True, False, False]


def test_frame_xy_inverts_lonlat():
    f = spot_sheet.Frame(FRAME)
    x, y = f.xy(*f.lonlat(np.array([10.0, 250.0]), np.array([600.0, -30.0])))
    assert np.allclose(x, [10.0, 250.0], atol=1e-4) and np.allclose(y, [600.0, -30.0], atol=1e-4)


def test_png_grey_round_trips():
    import zlib

    img = np.arange(12, dtype=np.uint8).reshape(3, 4) * 20
    png = spot_sheet.png_grey(img)
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert int.from_bytes(png[16:20], "big") == 4 and int.from_bytes(png[20:24], "big") == 3
    idat = png.index(b"IDAT")
    n = int.from_bytes(png[idat - 4 : idat], "big")
    raw = zlib.decompress(png[idat + 4 : idat + 4 + n])
    assert [raw[i * 5 + 1 : i * 5 + 5] for i in range(3)] == [bytes(img[i]) for i in range(3)]
