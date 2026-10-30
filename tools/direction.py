#!/usr/bin/env python3
"""Which way the waves of a finished 2D run travel where they break, and what slope they meet.

    uv run --python 3.12 --with numpy tools/direction.py out/pipeline_big_map

A slice cuts the reef along the bed's +x and treats the reef as the same all along the shore. This
asks how far that holds, from what the run saved, without simulating anything. Over the waves the
run recorded (its last periods) the surface's phase at the swell's own frequency changes from cell
to cell, and its gradient is the local wavenumber: it points the way the wave travels. At each
row's break point (from the run's break line), averaged over the half wavelength seaward of it, the
tool gives the wave's direction from +x, the bed's slope along x (what a slice meets) and along the
wave's own path, and the peel angle: the angle between the crest and the line of break points along
the shore, fitted over `--smooth` metres, 0 where the crest breaks all at once along that line.
Only the middle half of the width counts; the outer quarters are the run's side margins.

On Pipeline the waves arrive within 3 degrees of +x on half the rows, so a slice along their own
path meets the same slope within 2%; they peel, at about 20 degrees.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

G = 9.81


def load(run: Path) -> dict:
    """The run's header, its frames' surface, its bed and its frame times."""
    h = json.loads((run / "run.json").read_text())
    nx, ny = h["nx"], h["ny"]
    times = np.array(h["times"])
    frames = np.fromfile(run / h["frames"], dtype="<f4").reshape(len(times), len(h["fields"]), ny, nx)
    eta = frames[:, h["fields"].index("eta")].astype("float64")
    bed = np.fromfile(run / h["bed"], dtype="<f4").reshape(ny, nx).astype("float64")
    return {"header": h, "times": times, "eta": eta, "bed": bed}


def wavenumber(eta: np.ndarray, times: np.ndarray, period: float, dx: float) -> tuple[np.ndarray, np.ndarray]:
    """The local wavenumber (kx, ky) of the part of `eta` (time, y, x) at the swell's frequency.

    The complex amplitude A = sum of eta exp(-i w t) carries the wave's phase as exp(-i k.x), so
    the wave travels down the phase gradient: k = -Im(conj(A) grad A) / |A|^2, which needs no
    unwrapping.
    """
    w = 2 * math.pi / period
    a = np.tensordot(np.exp(-1j * w * times), eta, axes=1)
    gy, gx = np.gradient(a, dx)
    mag2 = np.abs(a) ** 2 + 1e-30
    return -(np.conj(a) * gx).imag / mag2, -(np.conj(a) * gy).imag / mag2


def at_breaks(run: dict, smooth: float = 63.0) -> list[dict]:
    """Direction, slopes and peel angle at each break point of the middle half of the width."""
    h = run["header"]
    nx, ny, dx = h["nx"], h["ny"], h["dx"]
    period = h["waves"]["wave_period_s"]
    stats = h["stats"]
    use = (run["times"] >= stats["from_s"] - 1e-9) & (run["times"] <= stats["to_s"] + 1e-9)
    kx, ky = wavenumber(run["eta"][use], run["times"][use], period, dx)
    # +x is shoreward, so a bed rising towards the shore has a positive slope.
    by, bx = np.gradient(run["bed"], dx)
    rows = []
    for j in range(ny // 4, ny - ny // 4):
        p = h["waves"]["break_line"][j]
        if not p:
            continue
        i = min(max(int(round(p["x_m"] / dx - 0.5)), 0), nx - 1)
        half = 0.5 * period * math.sqrt(G * p["depth_m"])
        seg = slice(max(i - max(2, int(half / dx)), 0), i + 1)
        kxm, kym = kx[j, seg].mean(), ky[j, seg].mean()
        ux, uy = kxm / math.hypot(kxm, kym), kym / math.hypot(kxm, kym)
        rows.append(
            {
                "y": (j + 0.5) * h["dy"],
                "x": p["x_m"] + h["waves"].get("cropped_offshore_m", 0.0),
                "depth": p["depth_m"],
                "height": p["height_m"],
                "direction": math.degrees(math.atan2(kym, kxm)),
                "slope_x": bx[j, seg].mean(),
                "slope_wave": (bx[j, seg] * ux + by[j, seg] * uy).mean(),
            }
        )
    # The line of break points: at each row, the slope of a straight line fitted to the break
    # points within half of `smooth` either side (fewer at the ends, so a straight line stays
    # straight there), and the crest, square to the wavenumber.
    ys = np.array([r["y"] for r in rows])
    xb = np.array([r["x"] for r in rows])
    for r in rows:
        near = np.abs(ys - r["y"]) <= smooth / 2 + 1e-9
        g = np.polyfit(ys[near], xb[near], 1)[0] if near.sum() > 1 else 0.0
        line = math.degrees(math.atan2(1.0, g))
        r["peel"] = abs((r["direction"] + 90.0 - line + 90.0) % 180.0 - 90.0)
    return rows


def summary(rows: list[dict]) -> dict:
    """The medians and ranges the README quotes."""
    d = np.array([r["direction"] for r in rows])
    sx = np.array([r["slope_x"] for r in rows])
    sw = np.array([r["slope_wave"] for r in rows])
    peel = np.array([r["peel"] for r in rows])
    hb = np.array([r["height"] for r in rows])
    rising = sx > 0.01
    ratio = sw[rising] / sx[rising]
    span = max(1, int(round(60.0 / (rows[1]["y"] - rows[0]["y"])))) if len(rows) > 1 else 1
    return {
        "rows": len(rows),
        "direction_median": float(np.median(d)),
        "direction_middle_half": [float(np.percentile(d, 25)), float(np.percentile(d, 75))],
        "rows_over_10_deg": int((np.abs(d) > 10).sum()),
        "slope_ratio_median": float(np.median(ratio)) if ratio.size else None,
        "slope_ratio_middle_half": [float(np.percentile(ratio, 25)), float(np.percentile(ratio, 75))]
        if ratio.size
        else None,
        "steepest_along_wave": float(sw.max()),
        "peel_median": float(np.median(peel)),
        "peel_middle_half": [float(np.percentile(peel, 25)), float(np.percentile(peel, 75))],
        "height_range": [float(hb.min()), float(hb.max())],
        "height_change_60m": float(max((abs(hb[k + span] - hb[k]) for k in range(len(hb) - span)), default=0.0)),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path, help="a finished run's directory (for a nested run, its fine grid's)")
    ap.add_argument("--smooth", type=float, default=63.0, help="metres along the shore to fit the break line over")
    ap.add_argument("--every", type=int, default=10, help="print every this many rows")
    ap.add_argument("--y", type=float, action="append", default=[], help="also print the row nearest this y")
    args = ap.parse_args(argv)
    try:
        run = load(args.run)
    except (OSError, KeyError, ValueError) as e:
        print(f"{args.run}: not a finished run with frames and a break line ({e})", file=sys.stderr)
        return 1
    rows = at_breaks(run, args.smooth)
    if not rows:
        print(f"{args.run}: no break points in the middle half of the width", file=sys.stderr)
        return 1
    print("     y |  break x  depth height | direction | slope along x  along wave | peel")
    shown = rows[:: max(1, args.every)] + [min(rows, key=lambda r: abs(r["y"] - y)) for y in args.y]
    for r in shown:
        print(
            f"{r['y']:6.1f} | {r['x']:7.1f} {r['depth']:5.1f} {r['height']:5.1f}  | {r['direction']:+6.1f}    "
            f"| {r['slope_x']:7.3f}     {r['slope_wave']:7.3f}   | {r['peel']:4.0f}"
        )
    s = summary(rows)
    print(
        f"direction from +x: median {s['direction_median']:+.1f} deg, middle half "
        f"{s['direction_middle_half'][0]:+.1f} to {s['direction_middle_half'][1]:+.1f}; "
        f"{s['rows_over_10_deg']} of {s['rows']} rows more than 10 deg off"
    )
    if s["slope_ratio_median"] is not None:
        print(
            f"slope along the wave over slope along x: median {s['slope_ratio_median']:.2f}, middle half "
            f"{s['slope_ratio_middle_half'][0]:.2f} to {s['slope_ratio_middle_half'][1]:.2f}; "
            f"steepest along the wave {s['steepest_along_wave']:.3f}"
        )
    print(
        f"peel angle: median {s['peel_median']:.0f} deg, middle half {s['peel_middle_half'][0]:.0f} to "
        f"{s['peel_middle_half'][1]:.0f}"
    )
    print(
        f"height at the break {s['height_range'][0]:.1f} to {s['height_range'][1]:.1f} m; "
        f"largest change over 60 m along the shore {s['height_change_60m']:.1f} m"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
