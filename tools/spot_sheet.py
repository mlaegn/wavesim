#!/usr/bin/env python3
"""Boil a spot's finished runs down to its spot sheet: where the swell breaks for each direction.

    uv run --python 3.12 --with numpy --with scipy --with pyproj --with contourpy \
        tools/spot_sheet.py coxos

The spot's `[name.sheet]` table in spots.toml lists the runs, one per open-sea direction and
tide, and the band of rows to trust. For each run the sheet gives, in latitude and longitude:

- the breaking zone: cells where the surface is breaking in at least `zone_fraction` of the
  saved frames of the settled sea (every wave breaks there, roughly);
- the break points: the run's own break line (one point per row: where the wave stops growing
  and breaks, or reaches the rock unbroken), with the height and depth there;
- the takeoff: where the crest breaks first. The crest leaves the wave-maker straight, so the
  time it takes to reach each row's break point, summed cell by cell at the local phase speed,
  says which row breaks first and which way the break runs from there: the peel. Towards -y is
  to a surfer's right as they face the shore, so a peel that way is a right-hander;
- a few numbers near the takeoff: height at the break over the offshore height, depth, how far
  from the shoreline, and the share of rows where the wave reaches the rock unbroken.

Only the rows inside `trust_y_m` count: the rest is the run's side margins (README, "side
walls"). Writes `out/sheets/<name>.json`, small enough to ship with gogo.
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
from contourpy import FillType, contour_generator
from pyproj import Transformer
from scipy.ndimage import median_filter, uniform_filter

G = 9.81
RIDE_M = 150.0  # how far along the shore from the takeoff the numbers are read


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
    out = {"header": h, "bed": bed, "break_share": breaking.mean(axis=0), "coarse": None}
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


def read_run(run: dict, trust_y: tuple[float, float], zone_fraction: float) -> dict:
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

    ts = median_filter(np.array([r["t"] for r in rows]), size=max(3, int(30 / dy) | 1), mode="nearest")
    peaks = takeoffs(rows, ts, dy)
    if not peaks:
        return {"tide_m": tide, "takeoff": None, "takeoffs": [], "peel": None, "near_takeoff": None,
                "_rows": rows, **common(run, rows, frame, j0, j1, zone_fraction)}
    k0 = peaks[0]["k"]
    peak = rows[k0]

    # The peel: how the break runs away from the takeoff on each side, over RIDE_M. Towards -y
    # (lower rows) is to the surfer's right.
    sides = {}
    for side, sel in (("right", lambda r: peak["y"] - RIDE_M <= r["y"] < peak["y"]),
                      ("left", lambda r: peak["y"] < r["y"] <= peak["y"] + RIDE_M)):
        ks = [k for k, r in enumerate(rows) if sel(r)]
        if len(ks) < 5:
            sides[side] = None
            continue
        ys = np.array([rows[k]["y"] for k in ks])
        later = np.polyfit(np.abs(ys - peak["y"]), ts[ks] - ts[k0], 1)[0]  # s per m along the shore
        sides[side] = {
            "metres": float(abs(ys - peak["y"]).max()),
            "seconds_per_100m": round(float(later * 100), 1),
            "rock_share": round(float(np.mean([rows[k]["rock"] for k in ks])), 2),
        }
    ride = [r for r in rows if abs(r["y"] - peak["y"]) <= RIDE_M]

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
        "peel": sides,
        "near_takeoff": {
            "height_over_offshore": None if med("height", ride) is None else round(med("height", ride) / h0, 2),
            "reflection": med("reflection", ride),
            "depth_m": med("depth", ride),
            "from_shore_m": med("from_shore", ride),
            "rock_share": round(float(np.mean([r["rock"] for r in ride])), 2),
        },
        "_rows": rows,
        **common(run, rows, frame, j0, j1, zone_fraction),
    }


def common(run: dict, rows: list, frame: Frame, j0: int, j1: int, zone_fraction: float) -> dict:
    """What every run's entry carries, takeoff or not: the zone and the break points."""
    lon, lat = frame.lonlat([r["x"] for r in rows], [r["y"] for r in rows])
    return {
        "rock_share": round(float(np.mean([r["rock"] for r in rows])), 2),
        "breaks": [
            [round(float(o), 6), round(float(a), 6), round(r["height"], 2), round(r["depth"], 1), r["rock"]]
            for o, a, r in zip(lon, lat, rows)
        ],
        "zone": zone_polygons(run["break_share"], slice(j0, j1), zone_fraction, run["header"]["dx"], frame),
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
    a = {r["j"]: r["x_bed"] for r in mid["_rows"]}
    d = [a[r["j"]] - r["x_bed"] for r in low["_rows"] if r["j"] in a]
    return round(float(np.median(d)), 1) if d else None


def recurring_peaks(rows: list, min_s: float = 5.0, radius_m: float = 50.0) -> list:
    """Takeoffs that recur across the runs: the bed's own peaks, which the swell's direction and
    the tide only switch on and off. Greedy clusters of every takeoff standing out by `min_s`."""
    hits = []
    for r in rows:
        for tide_name, t in r["tides"].items():
            for p in t["takeoffs"]:
                if p["stands_out_s"] >= min_s:
                    hits.append((p["lat"], p["lon"], p["stands_out_s"], r["from_deg"], tide_name))
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
    runs = sum(len(r["tides"]) for r in rows)
    for c in clusters:
        c["share"] = round(len({tuple(x) for x in c["runs"]}) / runs, 2)
    return sorted(clusters, key=lambda c: -c["share"])


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

    rows, swell = [], None
    for entry in cfg["rows"]:
        tides = {}
        for tide_name in ("mid", "low"):
            if tide_name not in entry:
                continue
            run = load_run(Path(entry[tide_name]))
            w = run["header"]["waves"]
            swell = swell or {"height_m": w["wave_height_m"], "period_s": w["wave_period_s"]}
            tides[tide_name] = read_run(run, trust, level)
            tides[tide_name]["run"] = entry[tide_name]
        shift = tide_shift(tides["mid"], tides["low"]) if {"mid", "low"} <= tides.keys() else None
        for v in tides.values():
            v.pop("_rows")
        rows.append({"from_deg": entry["from"], "at_maker_deg": entry["maker"], "low_tide_seaward_m": shift, "tides": tides})
        mid = tides.get("mid", {})
        print(
            f"{entry['from']:3d} (maker {entry['maker']}): takeoff {mid.get('takeoff')} "
            f"peel {mid.get('peel')} near {mid.get('near_takeoff')} low tide moves it {shift} m seaward"
        )

    sheet = {
        "format": "wavesim-spot-sheet",
        "version": 1,
        "spot": args.name,
        "description": spot.get("description"),
        "note": "model, not observation",
        "seabed": {"source": spot.get("source"), "cell_m": spot.get("cell_m"), "attribution": spot.get("attribution")},
        "swell": swell,
        "zone_fraction": level,
        "trust_y_m": list(trust),
        "generated": datetime.date.today().isoformat(),
        "peaks": recurring_peaks(rows),
        "rows": rows,
    }
    for p in sheet["peaks"]:
        print(f"peak at {p['lat']:.5f}, {p['lon']:.5f}: in {p['share']:.0%} of the runs")
    out = args.out or Path("out/sheets") / f"{args.name}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(sheet, separators=(",", ":"))
    out.write_text(text)
    print(f"{out}: {len(rows)} directions, {out.stat().st_size / 1024:.0f} KB")
    page = (Path(__file__).parent / "sheet_preview.html").read_text().replace("/*SHEET*/null", text.replace("</", "<\\/"))
    out.with_suffix(".html").write_text(page)
    print(f"{out.with_suffix('.html')}: the preview map")
    return 0


if __name__ == "__main__":
    sys.exit(main())
