"""The figures of docs/v1-breaking-wave.md, from the V1 runs in out/.

    cd cases/v1-steep && uv run -q --python 3.12 --with numpy --with matplotlib --with contourpy \\
        python ../v1-basilisk/figures.py ../../img/v1

Needs out/slice-steep (wavesim slice steep), out/foam-v1-steep-0.05-x170 (cases/v1-steep) and
out/basilisk-v1-L12 (cases/v1-basilisk). Run from cases/v1-steep, where measure.py is.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path.cwd()))
import measure as m  # noqa: E402

OUT = Path("../../out")
INK, MUTED, GRID, BED = "#0b0b0b", "#52514e", "#e4e3df", "#d9d8d3"
WAVESLICE, OPENFOAM, BASILISK = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({
    "font.size": 9, "text.color": INK, "axes.labelcolor": MUTED, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white", "savefig.facecolor": "white",
})


def surface_at(grid, t):
    frames = dict(grid.frames(t - 0.02))
    tt = min(frames, key=lambda q: abs(q - t))
    lines, _ = grid.surfaces(grid.read(frames[tt]))
    return [ln for ln in lines if np.ptp(ln[:, 0]) > 0.05]


def sequence(path):
    """The surface at four moments: the three models on one another."""
    _, times, sf = m.read_slice(OUT / "slice-steep")
    foam = m.load_grid(OUT / "foam-v1-steep-0.05-x170")
    bas = m.load_grid(OUT / "basilisk-v1-L12")
    bed = np.array(bas.meta["profile"])
    ts = [15.6, 16.2, 16.5, 16.8]
    fig, axes = plt.subplots(len(ts), 1, figsize=(8, 7.2), sharex=True)
    for ax, t in zip(axes, ts):
        ax.fill_between(bed[:, 0], -bed[:, 1], -9, color=BED, lw=0)
        k = int(np.argmin(abs(times - t)))
        ax.plot(sf[k, 0], sf[k, 1], color=WAVESLICE, lw=2.4, label="waveslice (potential flow)")
        for i, ln in enumerate(surface_at(foam, t)):
            ax.plot(ln[:, 0], ln[:, 1], color=OPENFOAM, lw=1.2,
                    label="OpenFOAM interFoam, 0.05 m" if i == 0 else None)
        for i, ln in enumerate(surface_at(bas, t)):
            ax.plot(ln[:, 0], ln[:, 1], color=BASILISK, lw=1.2,
                    label="Basilisk, 0.063 m" if i == 0 else None)
        ax.set_xlim(214, 236)
        ax.set_ylim(-0.6, 3.6)
        ax.set_aspect("equal")
        ax.set_yticks([0, 2])
        ax.text(0.005, 0.92, f"t = {t:.1f} s", transform=ax.transAxes, va="top", color=INK)
    axes[-1].set_xlabel("distance along the tank (m)")
    axes[1].set_ylabel("height above still water (m)")
    axes[0].legend(loc="lower right", bbox_to_anchor=(1, 1.02), ncol=3, frameon=False,
                   fontsize=8, handlelength=1.6)
    fig.tight_layout(h_pad=0.4)
    fig.savefig(path, dpi=150)


def overturn(path):
    """A/H² against cell size: OpenFOAM converges below the measured range, Basilisk into it."""
    foam = {0.4: (0.97, 3.03), 0.2: (1.60, 3.14), 0.1: (1.37, 3.14), 0.05: (1.14, 3.15)}
    foam_var = {"Courant 0.05": (0.2, 1.68, 3.17), "isoAdvector": (0.1, 1.20, 3.14)}
    bas = {0.25: (1.89, 3.31), 0.126: (2.42, 3.23), 0.063: (2.55, 3.22)}
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.axhspan(0.20, 0.42, color="#efeeea", lw=0)
    ax.text(0.47, 0.413, "measured: Surf Ranch lidar, 22 waves (Feddersen et al. 2023)", ha="left",
            va="top", color=MUTED, fontsize=8)
    ax.axhline(0.30, color=WAVESLICE, lw=2, zorder=2)
    ax.text(0.47, 0.305, "waveslice 0.30", ha="left", va="bottom", color=INK, fontsize=8)
    x = np.array(sorted(foam, reverse=True))
    y = np.array([foam[k][0] / foam[k][1] ** 2 for k in x])
    ax.plot(x, y, color=OPENFOAM, lw=2, marker="o", ms=7, mec="white", mew=1.5, zorder=3)
    ax.text(x[-1] * 0.92, y[-1], f"OpenFOAM {y[-1]:.2f}", ha="left", va="center", fontsize=8)
    for lab, (dx, area, h) in foam_var.items():
        ax.plot(dx, area / h ** 2, ls="none", marker="o", ms=7, mfc="white", mec=OPENFOAM,
                mew=1.5, zorder=3)
        ax.text(dx, area / h ** 2 - 0.014, lab, ha="center", va="top", color=MUTED,
                fontsize=7.5)
    x = np.array(sorted(bas, reverse=True))
    y = np.array([bas[k][0] / bas[k][1] ** 2 for k in x])
    ax.plot(x, y, color=BASILISK, lw=2, marker="o", ms=7, mec="white", mew=1.5, zorder=3)
    ax.text(x[-1] * 0.92, y[-1], f"Basilisk {y[-1]:.2f}", ha="left", va="center", fontsize=8)
    ax.set_xscale("log")
    ax.set_xlim(0.5, 0.032)
    ax.set_xticks([0.4, 0.2, 0.1, 0.05])
    ax.set_xticklabels(["0.4", "0.2", "0.1", "0.05"])
    ax.minorticks_off()
    ax.set_ylim(0, 0.45)
    ax.set_xlabel("cell size at the surface (m), finer to the right")
    ax.set_ylabel("overturn area / H²")
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def air(path, t=15.5):
    """Horizontal air speed around the crest: OpenFOAM's air rushes onshore, Basilisk's does not."""
    foam = m.load_grid(OUT / "foam-v1-steep-0.05-x170")
    fa = foam.read(OUT / "foam-v1-steep-0.05-x170" / f"{t:g}" / "alpha.water")
    fu = m.read_vectors(OUT / "foam-v1-steep-0.05-x170" / f"{t:g}" / "U")[foam.order][..., 0]
    bas = m.load_grid(OUT / "basilisk-v1-L12")
    ba = bas.read(OUT / "basilisk-v1-L12" / "frames" / f"f-{t:.3f}.bin")
    bu = np.fromfile(OUT / "basilisk-v1-L12" / "frames" / f"u-{t:.3f}.bin",
                     dtype="<f4").reshape(2, *bas.x.shape)[0].astype(float)
    fig, axes = plt.subplots(2, 1, figsize=(8, 3.9), sharex=True)
    for ax, (g, a, u, lab) in zip(axes, [(foam, fa, fu, "OpenFOAM interFoam, 0.05 m"),
                                        (bas, ba, bu, "Basilisk, 0.063 m")]):
        sel = (g.x[0] > 198) & (g.x[0] < 232)
        x, z = g.x[:, sel], g.z[:, sel]
        ua = np.ma.masked_where(a[:, sel] > 0.5, u[:, sel])
        mesh = ax.pcolormesh(x, z, ua, cmap="RdBu_r", vmin=-15, vmax=15, shading="nearest",
                             rasterized=True)
        ax.contourf(x, z, a[:, sel], levels=[0.5, 2], colors=["#c9d8ea"])
        ax.contour(x, z, a[:, sel], levels=[0.5], colors=[INK], linewidths=0.8)
        ax.set_ylim(-1, 6)
        ax.set_aspect("equal")
        ax.text(0.005, 0.95, lab, transform=ax.transAxes, va="top", color=INK)
    axes[-1].set_xlabel("distance along the tank (m)")
    for ax in axes:
        ax.set_ylabel("height (m)")
    cb = fig.colorbar(mesh, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label("air speed along the tank (m/s), onshore > 0")
    fig.suptitle(f"t = {t} s, as the face steepens; the wave itself travels at about 10 m/s",
                 x=0.07, ha="left", fontsize=9, color=MUTED)
    fig.savefig(path, dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../../img/v1")
    out.mkdir(parents=True, exist_ok=True)
    sequence(out / "sequence.png")
    overturn(out / "overturn.png")
    air(out / "air.png")
