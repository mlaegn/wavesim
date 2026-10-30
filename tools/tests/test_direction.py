"""The direction tool must read a known wave's direction, slopes and peel angle back out of a run.

A plane wave crossing a plane beach at 10 degrees to +x, written as a run would write it: the
tool must find the angle, the slope along x and along the wave (cos 10 degrees of it), and a peel
of 10 degrees where the waves break along a line straight along the shore, 0 where they break along
a line parallel to their crests, all at once.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import direction  # noqa: E402

NX, NY, D = 200, 40, 3.0
PERIOD, SLOPE, ANGLE = 12.0, 0.04, 10.0


def write_run(dir: Path, break_x) -> Path:
    x = (np.arange(NX) + 0.5) * D
    y = (np.arange(NY) + 0.5) * D
    bed = np.broadcast_to(-(16.0 - SLOPE * x), (NY, NX)).astype("<f4")
    times = np.arange(0.0, 4 * PERIOD, 1.5)
    k = 2 * math.pi / 60.0
    kx, ky = k * math.cos(math.radians(ANGLE)), k * math.sin(math.radians(ANGLE))
    w = 2 * math.pi / PERIOD
    xx, yy = np.meshgrid(x, y)
    frames = np.zeros((len(times), 2, NY, NX), dtype="<f4")
    for n, t in enumerate(times):
        frames[n, 0] = 0.5 * np.cos(kx * xx + ky * yy - w * t)
    bed.tofile(dir / "bed.f32")
    frames.tofile(dir / "frames.f32")
    breaks = [
        {"x_m": break_x(yj), "depth_m": 16.0 - SLOPE * break_x(yj), "height_m": 1.0} for yj in y
    ]
    header = {
        "nx": NX,
        "ny": NY,
        "dx": D,
        "dy": D,
        "bed": "bed.f32",
        "frames": "frames.f32",
        "fields": ["eta", "breaking"],
        "times": times.tolist(),
        "stats": {"from_s": 0.0, "to_s": float(times[-1])},
        "waves": {"wave_period_s": PERIOD, "cropped_offshore_m": 0.0, "break_line": breaks},
    }
    (dir / "run.json").write_text(json.dumps(header))
    return dir


def test_a_plane_wave_at_an_angle_is_read_back(tmp_path):
    rows = direction.at_breaks(direction.load(write_run(tmp_path, lambda y: 300.0)))
    assert len(rows) == NY // 2
    for r in rows:
        assert math.isclose(r["direction"], ANGLE, abs_tol=0.2), r
        assert math.isclose(r["slope_x"], SLOPE, rel_tol=1e-3), r
        assert math.isclose(r["slope_wave"], SLOPE * math.cos(math.radians(ANGLE)), rel_tol=1e-2), r
        assert math.isclose(r["peel"], ANGLE, abs_tol=0.2), r


def test_waves_breaking_along_their_own_crests_do_not_peel(tmp_path):
    # Crests run at 100 degrees from +x; so does a break line whose x falls as y grows by tan 10.
    tilt = math.tan(math.radians(ANGLE))
    rows = direction.at_breaks(direction.load(write_run(tmp_path, lambda y: 330.0 - tilt * y)))
    for r in rows:
        assert r["peel"] < 0.5, r
