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


def test_a_break_point_off_its_neighbours_or_too_deep_to_break_is_stray():
    rows = [{"x": 500.0, "height": 3.0, "depth": 4.0} for _ in range(60)]
    rows[20]["x"] = 400.0  # one row read 100 m seaward
    for k in range(40, 46):  # a run of rows on the ledge in 9 m of water with a 2 m wave
        rows[k].update(height=2.0, depth=9.0)
    spot_sheet.mark_strays(rows)
    assert [k for k, r in enumerate(rows) if r["stray"]] == [20, *range(40, 46)]


def test_the_zone_runs_from_each_break_point_to_the_shore():
    # A beach rising 1 in 30 from 6 m deep at x = 0, at a tide 0.5 m up: the shoreline is the
    # first cell whose bed reaches 0.5 m, at x = 195 m (cell 65 of 3 m cells).
    nx, ny, dx = 100, 30, 3.0
    x = (np.arange(nx) + 0.5) * dx
    bed = np.tile(-6.0 + x / 30.0, (ny, 1))
    run = {"header": {"nx": nx, "dx": dx, "waves": {"tide_m": 0.5}}, "bed": bed}
    shore = int(np.nonzero(bed[0] >= 0.5)[0][0])
    rows = [{"j": j, "x": 120.0 + j, "height": 2.0, "depth": 3.0} for j in range(ny)]
    rows[12]["x"] = 30.0  # one row misread 100 m seaward: it starts with its neighbours
    del rows[20]  # a row without a break point has no zone
    spot_sheet.mark_strays(rows)
    zone = spot_sheet.zone_mask(run, rows)
    for r in rows:
        start = int(round((r["x_kept"] if r["stray"] else r["x"]) / dx - 0.5))
        expect = np.zeros(nx)
        expect[start:shore] = 1.0
        assert np.array_equal(zone[r["j"]], expect), r["j"]
    assert rows[12]["stray"] and not zone[12, :35].any()
    assert not zone[20].any()
    # Never on land, never seaward of a break point.
    assert not zone[:, shore:].any()


def beach_run(height_over_depth: float, tide: float = 0.0):
    """A beach rising 1 in 30 from 6 m deep at x = 0 on 1 m cells, four rows, whose wave
    height is `height_over_depth` times the depth everywhere wet, and 2.5 m seaward of x = 60 m
    whatever the depth; a break point at x = 120.5 m, 2 m deep."""
    nx, ny, dx = 200, 4, 1.0
    x = (np.arange(nx) + 0.5) * dx
    bed = np.tile(-6.0 + x / 30.0, (ny, 1))
    depth = np.maximum(tide - bed, 0.0)
    height = np.where(x < 60.0, 2.5, height_over_depth * depth)
    run = {"header": {"nx": nx, "dx": dx, "waves": {"tide_m": tide}}, "bed": bed, "height": height}
    rows = [{"j": j, "x": 120.5, "stray": False} for j in range(ny)]
    return run, rows, x, depth


def test_the_surf_cap_falls_from_the_break_to_the_laboratorys_share_of_the_depth():
    run, rows, _, depth = beach_run(1.0)
    cap = spot_sheet.surf_cap(run, rows)
    capped = run["height"] * cap
    i0 = 120
    d0 = depth[0, i0]
    assert d0 == pytest.approx(2.0, abs=0.02)
    # Seaward of the break, and the break point's own cell: exactly as the run had it.
    assert np.array_equal(cap[:, :i0 + 1], np.ones((4, i0 + 1)))
    # It settles at SURF_HD 5.6 break depths in and holds there, to the shore.
    settled = i0 + int(math.ceil(spot_sheet.SURF_TRANSITION * d0))
    shore = int(np.nonzero(depth[0] == 0)[0][0])
    ratio = capped[0, settled:shore] / depth[0, settled:shore]
    assert np.allclose(ratio, spot_sheet.SURF_HD)
    # On the way it falls steadily from the break's own 1.0.
    between = capped[0, i0:settled] / depth[0, i0:settled]
    assert between[0] == pytest.approx(1.0) and np.all(np.diff(between) < 0)
    # Dry land and the cells past the shoreline are left alone.
    assert np.all(cap[:, shore:] == 1.0)


def test_a_surf_zone_already_under_the_cap_is_left_alone():
    run, rows, _, _ = beach_run(0.5)
    assert np.array_equal(spot_sheet.surf_cap(run, rows), np.ones(run["bed"].shape))


def test_a_stray_row_is_capped_from_its_kept_break_point_and_a_row_without_one_not_at_all():
    run, rows, _, _ = beach_run(1.0)
    rows[1].update(stray=True, x=30.0, x_kept=120.5)
    rows[2].update(stray=True, x=30.0, x_kept=None)
    del rows[3]
    cap = spot_sheet.surf_cap(run, rows)
    assert np.array_equal(cap[1], cap[0])
    assert np.array_equal(cap[2], np.ones(200)) and np.array_equal(cap[3], np.ones(200))


def test_a_long_run_of_stray_rows_does_not_drag_the_zone_out_to_sea():
    # 60 rows breaking at x = 120 m, but rows 20-49 have their tallest wave far out in deep water
    # (the swell never peaks inshore there). Their neighbours' plain median is stray too, so the
    # zone must not start from it: the run's edges start with the kept rows, its middle has none.
    nx, ny, dx = 100, 60, 3.0
    x = (np.arange(nx) + 0.5) * dx
    bed = np.tile(-6.0 + x / 30.0, (ny, 1))
    run = {"header": {"nx": nx, "dx": dx, "waves": {"tide_m": 0.0}}, "bed": bed}
    rows = [{"j": j, "x": 120.0, "height": 2.0, "depth": 3.0} for j in range(ny)]
    for j in range(20, 50):
        rows[j].update(x=15.0, height=3.5, depth=20.0)
    spot_sheet.mark_strays(rows)
    zone = spot_sheet.zone_mask(run, rows)
    assert all(rows[j]["stray"] for j in range(20, 50))
    assert not zone[:, : int(120.0 / dx) - 1].any()  # nothing seaward of the kept break points
    assert zone[20:30, 40].all() and zone[40:50, 40].all()  # within STRAY_ROWS of a kept row
    assert not zone[31:39].any()  # no kept row within reach: no zone


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


def test_harmonics_are_as_many_as_the_frames_phases_allow():
    # 2 s frames in a 13 s sea land on 13 phases; in a 12 s sea on 6; 1 s frames on 12.
    assert spot_sheet.harmonics(np.arange(0, 52, 2.0), 13.0)[1] == 6
    assert spot_sheet.harmonics(np.arange(0, 52, 2.0), 12.0)[1] == 2
    assert spot_sheet.harmonics(np.arange(0, 52, 1.0), 12.0)[1] == 5


def test_a_periodic_sea_is_rebuilt_at_any_phase():
    t = 220.0 + np.arange(26) * 2.0
    w = 2 * np.pi / 13.0
    sea = lambda tt: 0.1 + np.cos(w * tt) + 0.4 * np.sin(2 * w * tt) - 0.2 * np.cos(5 * w * tt)  # noqa: E731
    pinv, k = spot_sheet.harmonics(t, 13.0)
    coef = pinv @ sea(t)[:, None]
    for ph in (0.0, 0.13, 0.5, 0.77):
        assert spot_sheet.at_phase(coef, k, ph)[0] == pytest.approx(sea(ph * 13.0), abs=1e-9)
