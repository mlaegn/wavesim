"""The peel figures of docs/3d-box.md, from the whole-slope peel run in out/.

    cd cases/v1-steep && uv run -q --python 3.12 --with numpy --with matplotlib \\
        --with contourpy --with pillow python ../v1-basilisk/peel_figures.py ../../img/peel

Needs out/basilisk-v1-peel20-full-L10 (`run.sh 10 18000 --peel 20 full`). Run from
cases/v1-steep, where measure.py is. Writes sections.png (three cross-sections over one another)
and peel.gif (the surface in 3D every 0.1 s, mirrored at the strip's wall).
"""
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

sys.path.insert(0, str(Path.cwd()))
import measure as m  # noqa: E402

RUN = Path("../../out/basilisk-v1-peel20-full-L10")
INK, MUTED = "#0b0b0b", "#52514e"
SECTIONS = [("sec-1.00", "z = 1 m (the bed's peak)", "#2a78d6"),
            ("sec-8.06", "z = 8 m", "#1baf7a"),
            ("sec-15.12", "z = 15 m", "#eb6834")]
plt.rcParams.update({
    "font.size": 9, "text.color": INK, "axes.labelcolor": MUTED, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white", "savefig.facecolor": "white",
})


def sections(path):
    """The three cross-sections at five moments: the break reaches each in turn."""
    ts = [16.3, 16.6, 16.9, 17.2, 17.5]
    grids = [(m.load_grid(RUN / s), lab, c) for s, lab, c in SECTIONS]
    fig, axes = plt.subplots(len(ts), 1, figsize=(8, 8.5), sharex=True)
    for ax, t in zip(axes, ts):
        for g, lab, c in grids:
            bed = np.array(g.meta["profile"])
            ax.plot(bed[:, 0], -bed[:, 1], color=c, lw=0.8, ls=":")
            frames = dict(g.frames(t - 0.02))
            tt = min(frames, key=lambda q: abs(q - t))
            lines, _ = g.surfaces(g.read(frames[tt]))
            for i, ln in enumerate([ln for ln in lines if np.ptp(ln[:, 0]) > 0.05]):
                ax.plot(ln[:, 0], ln[:, 1], color=c, lw=1.3,
                        label=lab if (i == 0 and t == ts[0]) else None)
        ax.set_xlim(218, 242)
        ax.set_ylim(-1.2, 3.8)
        ax.set_aspect("equal")
        ax.set_yticks([0, 2])
        ax.text(0.01, 0.92, f"t = {t:.1f} s", transform=ax.transAxes, va="top")
    axes[0].legend(loc="lower right", bbox_to_anchor=(1, 1.0), ncol=3, frameon=False, fontsize=8)
    axes[2].set_ylabel("height above still water (m)")
    axes[-1].set_xlabel("distance along the tank (m); dotted: the bed under each section")
    fig.tight_layout(h_pad=0.3)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def facets(path, x0, x1):
    """Basilisk's output_facets, as triangles, cropped to x0-x1 and mirrored at z = 0."""
    tris = []
    for block in path.read_text().split("\n\n"):
        v = [ln.split() for ln in block.strip().split("\n") if ln.strip()]
        if len(v) < 3:
            continue
        p = np.array(v, float)
        if not (x0 < p[:, 0].mean() < x1 and p[:, 1].mean() > -0.5):
            continue
        for k in range(1, len(p) - 1):
            tris.append([p[0], p[k], p[k + 1]])
    t = np.array(tris)
    return np.concatenate([t, t * [1, 1, -1]])


def animation(path, x0=219., x1=243.):
    """The surface in 3D through the break, seen from offshore and above."""
    from PIL import Image

    light = np.array([0.5, 0.75, 0.45])
    light /= np.linalg.norm(light)
    deep, crest = np.array([0.05, 0.30, 0.55]), np.array([0.55, 0.85, 0.95])
    images = []
    for f in sorted(RUN.glob("facets-*.txt")):
        t = float(re.search(r"facets-([\d.]+)\.txt", f.name).group(1))
        tris = facets(f, x0, x1)
        n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
        n /= np.linalg.norm(n, axis=1)[:, None] + 1e-12
        shade = 0.3 + 0.7 * np.abs(n @ light)
        h = np.clip(tris[:, :, 1].mean(1) / 3.4, 0, 1)[:, None]
        cols = np.clip((deep * (1 - h) + crest * h) * shade[:, None], 0, 1)
        fig = plt.figure(figsize=(9, 5.4))
        ax = fig.add_subplot(projection="3d")
        # matplotlib's axes: along the tank, alongshore, height
        ax.add_collection3d(Poly3DCollection(tris[:, :, [0, 2, 1]], facecolors=cols,
                                             edgecolors="none"))
        ax.set_xlim(x0, x1)
        ax.set_ylim(-16, 16)
        ax.set_zlim(-0.5, 3.6)
        ax.set_box_aspect((x1 - x0, 32, 7))
        ax.view_init(elev=30, azim=-50)
        ax.set_axis_off()
        fig.text(0.03, 0.93, f"t = {t:.1f} s", fontsize=11)
        fig.text(0.03, 0.88, "depth lines at 20° to the crest; mirrored at the strip's wall, so "
                 "the bed's peak runs down the middle", fontsize=8.5, color=MUTED)
        fig.subplots_adjust(0, 0, 1, 1)
        fig.canvas.draw()
        rgb = np.asarray(fig.canvas.buffer_rgba())[..., :3]
        img = Image.fromarray(rgb)
        img = img.resize((img.width * 4 // 5, img.height * 4 // 5), Image.LANCZOS)
        images.append(img.convert("P", palette=Image.ADAPTIVE, colors=64))  # small for the repo
        plt.close(fig)
    images[0].save(path, save_all=True, append_images=images[1:], duration=350, loop=0,
                   optimize=True)


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../../img/peel")
    out.mkdir(parents=True, exist_ok=True)
    sections(out / "sections.png")
    animation(out / "peel.gif")
