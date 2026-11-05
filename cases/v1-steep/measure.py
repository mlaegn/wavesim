"""Measure an OpenFOAM run of V1 the way `wavesim slice` measures its own, and compare.

The surface is the `alpha.water = 0.5` contour of each saved frame. On it, as in waveslice:

- **face vertical** (the slice's `curl`): the first frame in which the front face, followed from
  the crest down to half its height, runs back on itself along x by at least 10 cm, and goes
  on doing so for the next 0.1 s: the face has passed vertical. Its x is where the overhang
  starts, its crest the highest point of the surface then, its depth the still water under it.
- **landing**: the first frame with a pocket of air of at least 0.1 m² closed off by water,
  near and shoreward of the face: the lip has reached the water ahead. Its x is the pocket's
  shoreward end, the throw that less the face's x, and the tube the pocket's area, width and
  height.

The limits are in metres, the same at every cell size, so that finer cells do not lower them.
waveslice's surface is smooth, and any fold-back there is the face. A VOF surface also has
spray thrown off the crest, which folds back by up to 0.3 m for a frame or two, and bubbles of
up to 0.03 m², and the finer the cells the more of both it resolves.

    measure.py out/foam-v1-steep-0.4 [--slice out/slice-steep]
    measure.py out/basilisk-v1-L12          # a Basilisk run of cases/v1-basilisk/

A Basilisk run writes f on a uniform grid over the break (its `case.json` says where), with the
bed counted as water, and is contoured the same way.

Writes `measure.json` and `compare.png` into the case.
"""

import argparse
import json
from pathlib import Path

import contourpy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.path import Path as Polygon  # noqa: E402


def read_scalars(path: Path) -> np.ndarray:
    text = path.read_text()
    body = text[text.index("internalField"):]
    head, rest = body.split("\n", 1)
    if "nonuniform" not in head:
        raise ValueError(f"{path}: expected a nonuniform field")
    n_line, rest = rest.lstrip().split("\n", 1)
    n = int(n_line)
    start, end = rest.index("("), rest.index(")")
    vals = np.array(rest[start + 1:end].split(), dtype=float)
    assert len(vals) == n, (path, len(vals), n)
    return vals


def read_vectors(path: Path) -> np.ndarray:
    text = path.read_text()
    body = text[text.index("internalField"):]
    n = int(body.split("\n")[1])
    start = body.index("(", body.index("\n", body.index("\n") + 1))
    vals = body[start + 1:].replace("(", " ").replace(")", " ").split()[: 3 * n]
    return np.array(vals, dtype=float).reshape(n, 3)


class Grid:
    """The case's cells as a structured (z, x) array: blockMesh numbers each block's cells x
    fastest, then z, block after block along x."""

    def __init__(self, case: Path):
        self.case = case
        self.meta = json.loads((case / "case.json").read_text())
        self.nz = self.meta["nz"]
        self.order = []
        offset = 0
        for n in self.meta["blocks_nx"]:
            idx = offset + np.arange(n * self.nz).reshape(self.nz, n)
            self.order.append(idx)
            offset += n * self.nz
        self.order = np.hstack(self.order)
        c = read_vectors(case / "0" / "C")
        self.x, self.z = c[:, 0][self.order], c[:, 2][self.order]
        assert np.all(np.diff(self.x, axis=1) > 0) and np.all(np.diff(self.z, axis=0) > 0), \
            "cells are not in the expected order"
        self.dx = self.meta["dx"]

    def read(self, path):
        return read_scalars(path)[self.order]

    def frames(self, start):
        return frames(self.case, start)

    def still_depth(self, x):
        xs, ds = zip(*self.meta["profile"])
        return float(np.interp(x, xs, ds))

    def surfaces(self, a):
        lines = contourpy.contour_generator(self.x, self.z, a, line_type="Separate").lines(0.5)
        return [ln for ln in lines if len(ln) >= 3], a


class BasiliskGrid(Grid):
    """A Basilisk run's samples: f at the centres of a uniform grid of the finest cells over
    the window x0-x1, z0-z1, as float32 rows of constant z from the bottom up."""

    def __init__(self, case: Path):
        self.case = case
        self.meta = json.loads((case / "case.json").read_text())
        self.dx = self.meta["dx"]
        x0, _, z0, _ = self.meta["window"]
        nx, nz = self.meta["nx"], self.meta["nz"]
        self.x, self.z = np.meshgrid(x0 + (np.arange(nx) + 0.5) * self.dx,
                                     z0 + (np.arange(nz) + 0.5) * self.dx)

    def read(self, path):
        return np.fromfile(path, dtype="<f4").reshape(self.x.shape).astype(float)

    def frames(self, start):
        out = [(float(p.stem[2:]), p) for p in (self.case / "frames").glob("f-*.bin")]
        return sorted(f for f in out if f[0] >= start - 1e-9)


def load_grid(case: Path):
    meta = json.loads((case / "case.json").read_text())
    return BasiliskGrid(case) if meta.get("code") == "basilisk" else Grid(case)


def polygon_area(p):
    x, z = p[:, 0], p[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(z, -1)) - np.dot(z, np.roll(x, -1)))


def closed(p, tol):
    return np.hypot(*(p[0] - p[-1])) < tol


FOLD = 0.1  # m, the least fold-back that counts as an overhang
LASTING = 0.1  # s, for how long it must last
POCKET = 0.1  # m², the least pocket of air that counts as a tube


def air_pockets(lines, grid, a):
    """Closed contours with air inside: the cells whose centres fall inside are mostly air."""
    out = []
    for p in lines:
        if not closed(p, 1e-6) or polygon_area(p) < POCKET:
            continue
        inside = Polygon(p).contains_points(np.c_[grid.x.ravel(), grid.z.ravel()])
        if inside.any() and a.ravel()[inside].mean() < 0.5:
            out.append(p)
    return out


def main_surface(lines):
    """The free surface: the open contour spanning the most of the tank along x."""
    open_lines = [p for p in lines if not closed(p, 1e-6)]
    p = max(open_lines, key=lambda q: np.ptp(q[:, 0]))
    return p if p[0, 0] < p[-1, 0] else p[::-1]


def overhang(surface):
    """How far the front face runs back along x, followed from the crest down to half its
    height: the largest drop below the furthest x reached so far, and the point where that
    furthest x was reached. Only the face: the ripples behind the crest, and those in the
    shallow water ahead, fold back by a quarter of a 0.1 m cell or more, and in the 0.2 and
    0.1 m runs they did so at 13 s, long before the wave reached the reef."""
    surface = surface[int(np.argmax(surface[:, 1])):]
    low = np.nonzero(surface[:, 1] < surface[0, 1] / 2)[0]
    surface = surface[: low[0] if len(low) else None]
    x = surface[:, 0]
    reach = np.maximum.accumulate(x)
    back = reach - x
    i = int(np.argmax(back))
    j = int(np.argmax(x[: i + 1])) if i > 0 else 0
    return back[i], surface[j]


def frames(case: Path, start: float):
    out = []
    for d in case.iterdir():
        try:
            t = float(d.name)
        except ValueError:
            continue
        if t >= start - 1e-9 and (d / "alpha.water").exists():
            out.append((t, d / "alpha.water"))
    return sorted(out)


def read_slice(slice_dir: Path):
    s = json.loads((slice_dir / "slice.json").read_text())
    data = np.fromfile(slice_dir / s["frames"], dtype="<f4")
    n = s["nodes"]
    fr = data.reshape(s["frame_count"], 3, n)
    return s, np.array(s["times"]), fr


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("case", type=Path)
    p.add_argument("--slice", type=Path, default=Path("out/slice-steep"))
    a = p.parse_args()
    grid = load_grid(a.case)
    curl = landing = None
    kept = {}
    all_frames = grid.frames(grid.meta["frames_from"])
    faces = []
    for t, path in all_frames:
        surface = main_surface(grid.surfaces(grid.read(path))[0])
        faces.append((t, surface, *overhang(surface)))
    for i, (t, surface, back, at) in enumerate(faces):
        ahead = [f for f in faces[i:] if f[0] <= t + LASTING + 1e-9]
        if ahead[-1][0] >= t + LASTING - 1e-9 and all(f[2] >= FOLD for f in ahead):
            curl = {"time": t, "x": float(at[0]), "crest": float(surface[:, 1].max()),
                    "depth": grid.still_depth(at[0]), "overhang": float(back)}
            kept["curl"] = (t, surface, [])
            break
    for t, path in all_frames if curl else []:
        if t < curl["time"]:
            continue
        lines, alpha = grid.surfaces(grid.read(path))
        surface = main_surface(lines)
        if landing is None:
            pockets = [q for q in air_pockets(lines, grid, alpha)
                       if q[:, 0].max() > curl["x"] - 2 * grid.dx and q[:, 0].min() < curl["x"] + 30]
            if pockets:
                q = max(pockets, key=polygon_area)
                landing = {"time": t, "x": float(q[:, 0].max()), "throw": float(q[:, 0].max() - curl["x"]),
                           "tube_area": polygon_area(q), "tube_width": float(np.ptp(q[:, 0])),
                           "tube_height": float(np.ptp(q[:, 1]))}
                kept["landing"] = (t, surface, [q])
                break
    s, times, sf = read_slice(a.slice)
    result = {"case": str(a.case), "dx": grid.dx, "cells": int(grid.x.size), "curl": curl,
              "landing": landing, "slice": {"curl": s["curl"], "landing": s["landing"]}}
    (a.case / "measure.json").write_text(json.dumps(result, indent=1))

    rows = [
        ("face vertical: time (s)", curl and curl["time"], s["curl"]["time"]),
        ("face vertical: x (m)", curl and curl["x"], s["curl"]["x"]),
        ("face vertical: still depth (m)", curl and curl["depth"], s["curl"]["depth"]),
        ("crest then (m)", curl and curl["crest"], s["curl"]["crest"]),
        ("lip lands: time (s)", landing and landing["time"], s["landing"]["time"]),
        ("lip thrown (m)", landing and landing["throw"], s["landing"]["throw"]),
        ("tube height (m)", landing and landing["tube_height"], s["landing"]["tube_height"]),
        ("tube width (m)", landing and landing["tube_width"], s["landing"]["tube_width"]),
        ("tube area (m²)", landing and landing["tube_area"], s["landing"]["tube_area"]),
    ]
    print(f"{a.case.name}: {grid.x.size} cells, {grid.dx:.3g} m")
    print(f"{'':32} {'OpenFOAM':>10} {'waveslice':>10} {'diff':>8}")
    for name, f, w in rows:
        fs = f"{f:10.2f}" if f is not None else f"{'—':>10}"
        d = f"{100 * (f - w) / w:+7.1f}%" if f is not None and w else ""
        print(f"{name:32} {fs} {w:10.2f} {d}")

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    bed = np.array(grid.meta["profile"])
    for ax, key in zip(axes, ("curl", "landing")):
        ax.fill_between(bed[:, 0], -bed[:, 1], -9, color="0.8")
        t_slice = s[key]["time"]
        k = int(np.argmin(abs(times - t_slice)))
        ax.plot(sf[k, 0], sf[k, 1], color="C1", lw=1.5, label=f"waveslice {times[k]:.2f} s")
        if key in kept:
            t, surf, pockets = kept[key]
            ax.plot(surf[:, 0], surf[:, 1], color="C0", lw=1, label=f"OpenFOAM {t:.2f} s")
            for q in pockets:
                ax.fill(q[:, 0], q[:, 1], color="C0", alpha=0.3)
        ax.set_xlim(195, 245)
        ax.set_ylim(-2, 5)
        ax.set_aspect("equal")
        ax.set_title({"curl": "face vertical", "landing": "lip lands"}[key])
        ax.legend(loc="upper left")
    axes[1].set_xlabel("x (m)")
    fig.tight_layout()
    fig.savefig(a.case / "compare.png", dpi=120)


if __name__ == "__main__":
    main()
