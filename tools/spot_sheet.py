#!/usr/bin/env python3
"""Boil a spot's finished runs down to its spot sheet: where the swell breaks for each direction.

    uv run --python 3.12 --with numpy --with scipy --with pyproj --with contourpy \
        tools/spot_sheet.py coxos

The spot's `[name.sheet]` table in spots.toml lists the runs, by swell (`[[name.sheet.swells]]`,
each a height and period) and within a swell one per open-sea direction and tide, and the band of
rows to trust. For each run the sheet gives, in latitude and longitude:

- the breaking zone: cells where the surface is breaking in at least `zone_fraction` of the
  saved frames of the settled sea (every wave breaks there, roughly);
- the break points: the run's own break line (one point per row: where the wave stops growing
  and breaks, or reaches the rock unbroken), with the height and depth there;
- the takeoffs: where the crest breaks before the reef on both sides. The crest leaves the
  wave-maker straight, so the time it takes to reach each row's break point, summed cell by cell
  at the local phase speed, says which rows break first. (Which way it then peels is left out:
  on the 15 m bed it is not consistent from one direction to the next.)
- a few numbers near the takeoff: height at the break over the offshore height, depth, how far
  from the shoreline, and the share of rows where the wave reaches the rock unbroken;
- the wave-height map: the run's wave height over the recorded waves, over the offshore
  height, smoothed over the wet cells (sigma HEIGHT_SMOOTH_M) so that the stripes the swell
  reflected off the rock lays over the incoming one, half a wavelength apart, do not read as
  peaks, and resampled onto a latitude-longitude grid of about one cell, as an 8-bit grey PNG
  (0 = no data; v = 1 + round(ratio / HEIGHT_TOP * 254)) that the page colours itself.

A break point more than STRAY_M seaward or shoreward of the median of its neighbours (STRAY_ROWS
either side) is a misreading of a single row on the 15 m bed, and is left out.

Only the rows inside `trust_y_m` count: the rest is the run's side margins (README, "side
walls"). An optional `area`, a polygon of [lat, lon] corners, keeps only what lies inside it: the
surf spot itself, where the runs show the whole coast. Writes `out/sheets/<name>.json`, small enough to ship with gogo.
"""

from __future__ import annotations

import argparse
import base64
import datetime
import json
import math
import sys
import tomllib
import zlib
from pathlib import Path

import numpy as np
from contourpy import FillType, contour_generator
from pyproj import Transformer
from scipy.ndimage import gaussian_filter, map_coordinates, median_filter, uniform_filter

G = 9.81
RIDE_M = 150.0  # how far along the shore from the takeoff the numbers are read
STRAY_M, STRAY_ROWS = 40.0, 10  # a break point this far off its neighbours' median is dropped
HEIGHT_TOP = 3.0  # the height map's grey levels span 0 to this many times the offshore height
HEIGHT_SMOOTH_M = 15.0  # Gaussian sigma: the swell reflected off the rock stripes the raw height


def phase_speed(depth: np.ndarray, period: float) -> np.ndarray:
    """Linear phase speed at each depth (Fenton and McKee's explicit wavelength, within 2%)."""
    l0 = G * period**2 / (2 * math.pi)
    h = np.maximum(depth, 0.1)
    return l0 * np.tanh((2 * math.pi * h / l0) ** 0.75) ** (2 / 3) / period


def load_run(run: Path) -> dict:
    h = json.loads((run / "run.json").read_text())
    nx, ny = h["nx"], h["ny"]
    bed = np.fromfile(run / h["bed"], dtype="<f4").reshape(ny, nx).astype("float64")
    times = np.array(h["times"])
    frames = np.memmap(run / h["frames"], dtype="<f4", mode="r", shape=(len(times), len(h["fields"]), ny, nx))
    use = times >= h["stats"]["from_s"] - 1e-9
    breaking = np.asarray(frames[use, h["fields"].index("breaking")]) > 0
    stats = np.fromfile(run / h["stats"]["file"], dtype="<f4").reshape(len(h["stats"]["fields"]), ny, nx)
    height = stats[h["stats"]["fields"].index("wave_height")].astype("float64")
    out = {"header": h, "bed": bed, "break_share": breaking.mean(axis=0), "height": height, "coarse": None}
    if (run / "coarse" / "run.json").exists():
        c = json.loads((run / "coarse" / "run.json").read_text())
        cbed = np.fromfile(run / "coarse" / c["bed"], dtype="<f4").reshape(c["ny"], c["nx"]).astype("float64")
        out["coarse"] = {"header": c, "bed": cbed}
    return out


def crest_times(run: dict) -> np.ndarray:
    """Seconds for a straight crest to reach each fine cell, from the wave-maker (any constant
    offset is shared by every row, so only differences mean anything)."""
    h = run["header"]
    w = h["waves"]
    tide, period, dx = w["tide_m"], w["wave_period_s"], h["dx"]
    t = np.cumsum(dx / phase_speed(tide - run["bed"], period), axis=1)
    if run["coarse"] is not None:
        c = run["coarse"]["header"]
        cw, cdx = c["waves"], c["dx"]
        start = cw.get("maker_x_m", cw.get("sponge_offshore_cells", 0) * cdx)
        end = w["cropped_offshore_m"] - cw["cropped_offshore_m"]  # the fine grid's x = 0
        i0, i1 = int(start // cdx), int(end // cdx)
        ct = (cdx / phase_speed(tide - run["coarse"]["bed"][:, i0:i1], period)).sum(axis=1)
        rows = np.minimum(((np.arange(h["ny"]) + 0.5) * h["dy"] / c["dy"]).astype(int), c["ny"] - 1)
        t += ct[rows, None]
    return t


class Frame:
    """Run-frame metres (x along the swell, y 90 degrees counter-clockwise) to longitude, latitude."""

    def __init__(self, frame: dict):
        b = math.radians(frame["x_bearing_deg"])
        self.e0, self.n0 = frame["origin_easting"], frame["origin_northing"]
        self.ex, self.ey = (math.sin(b), math.cos(b)), (-math.cos(b), math.sin(b))
        self.to_ll = Transformer.from_crs(frame["crs"], "EPSG:4326", always_xy=True)
        self.to_utm = Transformer.from_crs("EPSG:4326", frame["crs"], always_xy=True)

    def xy(self, lon, lat):
        """Longitude, latitude to run-frame metres: the inverse of lonlat."""
        e, n = self.to_utm.transform(np.asarray(lon, float), np.asarray(lat, float))
        de, dn = e - self.e0, n - self.n0
        return de * self.ex[0] + dn * self.ex[1], de * self.ey[0] + dn * self.ey[1]

    def lonlat(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        e = self.e0 + x * self.ex[0] + y * self.ey[0]
        n = self.n0 + x * self.ex[1] + y * self.ey[1]
        return self.to_ll.transform(e, n)


def simplify(ring: np.ndarray, tol: float) -> np.ndarray:
    """Douglas-Peucker: drop the vertices of a closed ring that lie within `tol` of the line kept."""
    keep = np.zeros(len(ring), bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(ring) - 1)]
    while stack:
        a, b = stack.pop()
        if b - a < 2:
            continue
        p, q = ring[a], ring[b]
        d = ring[a + 1 : b] - p
        seg = q - p
        n = math.hypot(*seg)
        dist = np.abs(d[:, 0] * seg[1] - d[:, 1] * seg[0]) / n if n > 0 else np.hypot(d[:, 0], d[:, 1])
        k = int(np.argmax(dist))
        if dist[k] > tol:
            keep[a + 1 + k] = True
            stack += [(a, a + 1 + k), (a + 1 + k, b)]
    return ring[keep]


def ring_area(r: np.ndarray) -> float:
    return 0.5 * abs(float(np.dot(r[:, 0], np.roll(r[:, 1], 1)) - np.dot(r[:, 1], np.roll(r[:, 0], 1))))


def zone_polygons(share: np.ndarray, rows: slice, level: float, dx: float, frame: Frame) -> list:
    """Polygons (outer ring, then holes; [lon, lat] pairs) where `share` >= level, inside `rows`.

    Rings are simplified to within a cell, and pieces under 4 cells dropped: the bed is 15 m data
    on 3 m cells, so nothing finer means anything.
    """
    z = uniform_filter(share, size=3)  # one-cell speckle is not a zone
    z[: rows.start] = 0.0
    z[rows.stop :] = 0.0
    ny, nx = z.shape
    gen = contour_generator(x=(np.arange(nx) + 0.5) * dx, y=(np.arange(ny) + 0.5) * dx, z=z, fill_type=FillType.OuterOffset)
    points, offsets = gen.filled(level, np.inf)
    polys = []
    for pts, offs in zip(points, offsets):
        rings = []
        for k, (a, b) in enumerate(zip(offs[:-1], offs[1:])):
            r = simplify(pts[a:b], dx)
            if len(r) < 4 or ring_area(r) < 4 * dx * dx:
                if k == 0:
                    break  # the outer ring is too small: drop the piece, holes and all
                continue
            lon, lat = frame.lonlat(r[:, 0], r[:, 1])
            rings.append([[round(float(o), 5), round(float(t), 5)] for o, t in zip(lon, lat)])
        if rings:
            polys.append(rings)
    return polys


def inside(lon, lat, area) -> np.ndarray:
    """Whether each point lies inside the polygon `area` ([lat, lon] corners), by ray casting."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    out = np.zeros(lon.shape, bool)
    n = len(area)
    for k in range(n):
        (a_lat, a_lon), (b_lat, b_lon) = area[k], area[(k + 1) % n]
        crosses = (a_lat > lat) != (b_lat > lat)
        with np.errstate(divide="ignore", invalid="ignore"):
            at = a_lon + (lat - a_lat) * (b_lon - a_lon) / (b_lat - a_lat)
        out ^= crosses & (lon < at)
    return out


def read_run(run: dict, trust_y: tuple[float, float], zone_fraction: float, area=None) -> dict:
    """One run's entry in the sheet. Takeoffs are found over all the trusted rows, so the
    spot's `area` does not make edges of its own; then only what lies inside it is kept."""
    h = run["header"]
    w = h["waves"]
    ny, dx, dy = h["ny"], h["dx"], h["dy"]
    tide, h0 = w["tide_m"], w["wave_height_m"]
    frame = Frame(h["frame"])
    j0 = max(int(math.ceil(trust_y[0] / dy - 0.5)), 0)
    j1 = min(int(math.floor(trust_y[1] / dy - 0.5)) + 1, ny)

    t = crest_times(run)
    bed = run["bed"]
    rows = []
    for j in range(j0, j1):
        p = w["break_line"][j]
        if not p:
            continue
        i = min(max(int(round(p["x_m"] / dx - 0.5)), 0), h["nx"] - 1)
        dry = np.nonzero(bed[j, i:] >= tide)[0]
        shore_x = (i + dry[0] + 0.5) * dx if dry.size else None
        rows.append(
            {
                "j": j,
                "y": (j + 0.5) * dy,
                "x": p["x_m"],
                "x_bed": p["x_m"] + w.get("cropped_offshore_m", 0.0),  # runs crop differently: compare these
                "t": t[j, i],
                "height": p["height_m"],
                "depth": p["depth_m"],
                "rock": bool(p["shore_break"]),
                "reflection": p.get("reflection"),
                "from_shore": None if shore_x is None else shore_x - p["x_m"],
            }
        )
    if not rows:
        raise ValueError("no break points inside the trusted rows")
    lon, lat = frame.lonlat([r["x"] for r in rows], [r["y"] for r in rows])
    for r, o, a, ok in zip(rows, lon, lat, inside(lon, lat, area) if area else [True] * len(rows)):
        r["lon"], r["lat"], r["in_area"] = float(o), float(a), bool(ok)
    xs = median_filter(np.array([r["x"] for r in rows]), size=2 * STRAY_ROWS + 1, mode="nearest")
    for r, m in zip(rows, xs):
        r["stray"] = abs(r["x"] - m) > STRAY_M

    ts = median_filter(np.array([r["t"] for r in rows]), size=max(3, int(30 / dy) | 1), mode="nearest")
    peaks = [p for p in takeoffs(rows, ts, dy) if rows[p["k"]]["in_area"]]
    shared = common(run, rows, frame, j0, j1, zone_fraction, area)
    if not peaks:
        return {"tide_m": tide, "takeoff": None, "takeoffs": [], "near_takeoff": None, "_rows": rows, **shared}
    k0 = peaks[0]["k"]
    peak = rows[k0]
    ride = [r for r in rows if abs(r["y"] - peak["y"]) <= RIDE_M and r["in_area"] and not r["stray"]]

    def med(key, rs):
        v = [r[key] for r in rs if r[key] is not None]
        return round(float(np.median(v)), 2) if v else None

    def where(k, extra=None):
        lon, lat = frame.lonlat(rows[k]["x"], rows[k]["y"])
        return {"lon": round(float(lon), 6), "lat": round(float(lat), 6), "depth_m": round(rows[k]["depth"], 1), **(extra or {})}

    return {
        "tide_m": tide,
        "takeoff": where(k0),
        "takeoffs": [where(p["k"], {"stands_out_s": p["stands_out_s"]}) for p in peaks],
        "near_takeoff": {
            "height_over_offshore": None if med("height", ride) is None else round(med("height", ride) / h0, 2),
            "reflection": med("reflection", ride),
            "depth_m": med("depth", ride),
            "from_shore_m": med("from_shore", ride),
            "rock_share": round(float(np.mean([r["rock"] for r in ride])), 2),
        },
        "_rows": rows,
        **shared,
    }


def common(run: dict, rows: list, frame: Frame, j0: int, j1: int, zone_fraction: float, area=None) -> dict:
    """What every run's entry carries, takeoff or not: the zone and the break points, inside the
    spot's area."""
    share = run["break_share"]
    if area:
        h = run["header"]
        xs, ys = np.meshgrid((np.arange(h["nx"]) + 0.5) * h["dx"], (np.arange(h["ny"]) + 0.5) * h["dy"])
        share = np.where(inside(*frame.lonlat(xs, ys), area), share, 0.0)
    kept = [r for r in rows if r["in_area"] and not r["stray"]]
    return {
        "rock_share": round(float(np.mean([r["rock"] for r in kept])), 2) if kept else None,
        "breaks": [
            [round(r["lon"], 6), round(r["lat"], 6), round(r["height"], 2), round(r["depth"], 1), r["rock"]]
            for r in kept
        ],
        "zone": zone_polygons(share, slice(j0, j1), zone_fraction, run["header"]["dx"], frame),
        "height_map": height_map(run, frame, j0, j1, area),
    }


def png_grey(img: np.ndarray) -> bytes:
    """An 8-bit greyscale PNG, row 0 at the top. No dependency beyond zlib."""
    def chunk(kind: bytes, data: bytes) -> bytes:
        body = kind + data
        return len(data).to_bytes(4, "big") + body + zlib.crc32(body).to_bytes(4, "big")

    h, w = img.shape
    raw = b"".join(b"\x00" + img[i].astype(np.uint8).tobytes() for i in range(h))
    ihdr = w.to_bytes(4, "big") + h.to_bytes(4, "big") + bytes([8, 0, 0, 0, 0])
    return b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b"")


def height_map(run: dict, frame: Frame, j0: int, j1: int, area=None, pad_m: float = 30.0) -> dict | None:
    """The wave height over the offshore height on a latitude-longitude grid of about one cell,
    over the area (or the trusted rows), wet cells only: {bounds: [[s, w], [n, e]], png: base64}.
    """
    h = run["header"]
    dx, dy, ny, nx = h["dx"], h["dy"], h["ny"], h["nx"]
    tide, h0 = h["waves"]["tide_m"], h["waves"]["wave_height_m"]
    if area:
        lats, lons = [a for a, _ in area], [o for _, o in area]
    else:
        xs, ys = np.meshgrid([0.0, nx * dx], [j0 * dy, j1 * dy])
        lons, lats = frame.lonlat(xs, ys)
        lons, lats = list(np.ravel(lons)), list(np.ravel(lats))
    lat0 = (min(lats) + max(lats)) / 2
    dlat = dx / 111_320.0
    dlon = dx / (111_320.0 * math.cos(math.radians(lat0)))
    pad_lat, pad_lon = pad_m / 111_320.0, pad_m / (111_320.0 * math.cos(math.radians(lat0)))
    south, north = min(lats) - pad_lat, max(lats) + pad_lat
    west, east = min(lons) - pad_lon, max(lons) + pad_lon
    rows, cols = int(math.ceil((north - south) / dlat)), int(math.ceil((east - west) / dlon))
    lat = north - (np.arange(rows) + 0.5) * dlat  # row 0 at the top
    lon = west + (np.arange(cols) + 0.5) * dlon
    LON, LAT = np.meshgrid(lon, lat)
    x, y = frame.xy(LON, LAT)
    i, j = x / dx - 0.5, y / dy - 0.5
    ok = (i >= 0) & (i <= nx - 1) & (j >= j0) & (j <= j1 - 1)
    wet_cells = (run["bed"] < tide).astype("float64")
    s = HEIGHT_SMOOTH_M / dx
    smooth = gaussian_filter(run["height"] * wet_cells, s) / np.maximum(gaussian_filter(wet_cells, s), 1e-9)
    ratio = map_coordinates(smooth, [j, i], order=1, mode="nearest") / h0
    wet = map_coordinates(run["bed"], [j, i], order=1, mode="nearest") < tide
    ok &= wet & np.isfinite(ratio)
    if area:
        ok &= inside(LON, LAT, area)
    if not ok.any():
        return None
    grey = np.where(ok, 1 + np.round(np.clip(ratio, 0, HEIGHT_TOP) / HEIGHT_TOP * 254), 0)
    return {
        "bounds": [[round(south, 6), round(west, 6)], [round(south + rows * dlat, 6), round(west + cols * dlon, 6)]],
        "top": HEIGHT_TOP,
        "png": base64.b64encode(png_grey(grey)).decode(),
    }


def takeoffs(rows: list, ts: np.ndarray, dy: float, window_m: float = 60.0, reach_m: float = 150.0,
             edge_m: float = 30.0, min_s: float = 2.0) -> list:
    """Where the crest breaks before the shore on either side of it: the takeoffs, best first.

    A takeoff is a row whose (smoothed) crest time is the earliest within `window_m` either side,
    and earlier by at least `min_s` than the latest row within `reach_m` on *both* sides. A
    coast that simply runs seaward towards one end breaks first at that end, which is not a
    takeoff, so a row within `edge_m` of the trusted band's ends is never one. Nor is a row
    where the wave reaches the rock unbroken.
    """
    n = len(rows)
    w, reach, edge = int(window_m / dy), int(reach_m / dy), int(edge_m / dy)
    found = []
    for k in range(edge, n - edge):
        if rows[k]["rock"] or ts[k] > ts[max(0, k - w) : k + w + 1].min():
            continue
        if found and k - found[-1]["k"] <= w:  # the same flat bottom
            continue
        left, right = ts[max(0, k - reach) : k], ts[k + 1 : k + 1 + reach]
        if not left.size or not right.size:
            continue
        stands_out = min(left.max(), right.max()) - ts[k]
        if stands_out >= min_s:
            found.append({"k": k, "stands_out_s": round(float(stands_out), 1)})
    return sorted(found, key=lambda p: -p["stands_out_s"])


def tide_shift(mid: dict, low: dict) -> float | None:
    """How far seaward the break moves from mid to low tide: median over rows both share (m).

    Each run crops the bed's offshore end differently, so this compares positions on the bed.
    """
    a = {r["j"]: r["x_bed"] for r in mid["_rows"] if r.get("in_area", True)}
    d = [a[r["j"]] - r["x_bed"] for r in low["_rows"] if r["j"] in a]
    return round(float(np.median(d)), 1) if d else None


def recurring_peaks(swells: list, min_s: float = 5.0, radius_m: float = 50.0) -> list:
    """Takeoffs that recur across the runs: the bed's own peaks, which the swell's size and
    direction and the tide only switch on and off. Greedy clusters of every takeoff standing out
    by `min_s`, over every swell."""
    hits = []
    for sw in swells:
        for r in sw["rows"]:
            for tide_name, t in r["tides"].items():
                for p in t["takeoffs"]:
                    if p["stands_out_s"] >= min_s:
                        hits.append((p["lat"], p["lon"], p["stands_out_s"], f"{sw['name']} {r['from_deg']}", tide_name))
    hits.sort(key=lambda h: -h[2])
    clusters = []
    for lat, lon, s, deg, tide_name in hits:
        for c in clusters:
            dn = (lat - c["lat"]) * 111_320.0
            de = (lon - c["lon"]) * 111_320.0 * math.cos(math.radians(lat))
            if math.hypot(dn, de) <= radius_m:
                c["runs"].append([deg, tide_name])
                break
        else:
            clusters.append({"lat": lat, "lon": lon, "runs": [[deg, tide_name]]})
    runs = sum(len(r["tides"]) for sw in swells for r in sw["rows"])
    for c in clusters:
        c["share"] = round(len({tuple(x) for x in c["runs"]}) / runs, 2)
    return sorted(clusters, key=lambda c: -c["share"])


def read_swell(cfg: dict, trust, level: float, area) -> dict:
    """One swell's rows: a direction each, with its runs at each tide. Every run must have been
    made with the swell this table names."""
    rows = []
    for entry in cfg["rows"]:
        tides = {}
        for tide_name in ("mid", "low"):
            if tide_name not in entry:
                continue
            run = load_run(Path(entry[tide_name]))
            w = run["header"]["waves"]
            if (w["wave_height_m"], w["wave_period_s"]) != (cfg["height_m"], cfg["period_s"]):
                raise ValueError(
                    f"{entry[tide_name]} is a {w['wave_height_m']} m, {w['wave_period_s']} s run, "
                    f"listed under the {cfg['height_m']} m, {cfg['period_s']} s swell"
                )
            tides[tide_name] = read_run(run, trust, level, area)
            tides[tide_name]["run"] = entry[tide_name]
        shift = tide_shift(tides["mid"], tides["low"]) if {"mid", "low"} <= tides.keys() else None
        for v in tides.values():
            v.pop("_rows")
        rows.append({"from_deg": entry["from"], "at_maker_deg": entry["maker"], "low_tide_seaward_m": shift, "tides": tides})
        mid = tides.get("mid") or {}
        print(
            f"{cfg['height_m']} m {cfg['period_s']} s from {entry['from']:3d} (maker {entry['maker']}): "
            f"takeoff {mid.get('takeoff')} low tide moves it {shift} m seaward"
        )
    name = f"{cfg['height_m']:g} m, {cfg['period_s']:g} s"
    return {"name": name, "height_m": cfg["height_m"], "period_s": cfg["period_s"], "rows": rows}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name", help="spot name in spots.toml, with a [name.sheet] table")
    ap.add_argument("--spots", default="spots.toml", type=Path)
    ap.add_argument("--out", type=Path, help="default: out/sheets/<name>.json")
    args = ap.parse_args(argv)

    spots = tomllib.loads(args.spots.read_text())
    if args.name not in spots or "sheet" not in spots[args.name]:
        print(f"no [{args.name}.sheet] in {args.spots}", file=sys.stderr)
        return 2
    spot, cfg = spots[args.name], spots[args.name]["sheet"]
    trust = tuple(cfg["trust_y_m"])
    level = cfg.get("zone_fraction", 0.08)
    area = cfg.get("area")

    swells = [read_swell(s, trust, level, area) for s in cfg["swells"]]

    sheet = {
        "format": "wavesim-spot-sheet",
        "version": 1,
        "spot": args.name,
        "description": spot.get("description"),
        "note": "model, not observation",
        "seabed": {"source": spot.get("source"), "cell_m": spot.get("cell_m"), "attribution": spot.get("attribution")},
        "zone_fraction": level,
        "trust_y_m": list(trust),
        "area": area,
        "generated": datetime.date.today().isoformat(),
        "peaks": recurring_peaks(swells),
        "swells": swells,
    }
    for p in sheet["peaks"]:
        print(f"peak at {p['lat']:.5f}, {p['lon']:.5f}: in {p['share']:.0%} of the runs")
    out = args.out or Path("out/sheets") / f"{args.name}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(sheet, separators=(",", ":"))
    out.write_text(text)
    print(f"{out}: {len(swells)} swell(s), {sum(len(s['rows']) for s in swells)} directions, {out.stat().st_size / 1024:.0f} KB")
    page = (Path(__file__).parent / "sheet_preview.html").read_text().replace("/*SHEET*/null", text.replace("</", "<\\/"))
    out.with_suffix(".html").write_text(page)
    print(f"{out.with_suffix('.html')}: the preview map")
    return 0


if __name__ == "__main__":
    sys.exit(main())
