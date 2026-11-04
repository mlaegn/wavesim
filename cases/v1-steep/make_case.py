"""V1 of docs/3d-box.md: the `steep` slice case, set up for OpenFOAM's interFoam.

The same tank as `wavesim slice steep`: water 8 m deep, flat for 120 m, a plane slope of 1 in 15
up to a shelf 0.4 m deep, the shelf 24 m long, walls at both ends. The same wave: a solitary wave
2.4 m high with its crest 60 m from the left wall, from the Serre-Green-Naghdi theory that
`waveslice::Solitary` uses, with that theory's velocity profile in depth through the whole water
column (`waveslice` needs it only on the surface).

The mesh is one cell thick (2D in x and z) and follows the bed: every column has the same number
of cells between the bed and the top at +6 m, so cells are `dx` tall offshore and shorter where
the water is shallow. Walls and bed are free-slip, the model laminar, to compare with potential
flow up to the moment the lip lands.

    make_case.py mesh --dx 0.4 --out CASE   # dictionaries and fields with no wave yet
    make_case.py init --out CASE            # the wave, after blockMesh and writeCellCentres
"""

import argparse
import json
from pathlib import Path

import numpy as np

G = 9.81
DEPTH = 8.0
HEIGHT = 2.4
CREST = 60.0
TOP = 6.0
# The steep profile of crates/wavesim/src/slice.rs for an 8 m depth, as (x, still depth).
PROFILE = [(0.0, 8.0), (120.0, 8.0), (234.0, 0.4), (258.0, 0.4)]
END_TIME = 18.5
# Frames of alpha around the break: the slice's face goes vertical at 15.97 s, its lip lands
# at 16.93 s.
FRAMES_FROM, FRAME_EVERY = 13.0, 0.025

HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""


def write(case: Path, rel: str, body: str, cls: str = "dictionary") -> None:
    path = case / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HEADER.format(cls=cls, obj=path.name) + body)


def bed_z(x):
    xs, ds = zip(*PROFILE)
    return -np.interp(x, xs, ds)


def blocks(dx: float):
    """The three blocks along x (flat, slope, shelf): their x limits and cell counts, and the
    cells per column, the same in every block so their faces match."""
    xs = [p[0] for p in PROFILE]
    nx = [max(1, round((b - a) / dx)) for a, b in zip(xs, xs[1:])]
    nz = round((TOP + DEPTH) / dx)
    return xs, nx, nz


def mesh(case: Path, dx: float) -> None:
    xs, nx, nz = blocks(dx)
    t = dx
    vertices = []
    for x in xs:
        zb = float(bed_z(x))
        vertices += [(x, 0, zb), (x, 0, TOP), (x, t, zb), (x, t, TOP)]
    hexes, faces = [], {"leftWall": [], "rightWall": [], "bed": [], "atmosphere": [], "frontAndBack": []}
    for i, n in enumerate(nx):
        b0, t0, b0y, t0y = 4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3
        b1, t1, b1y, t1y = b0 + 4, t0 + 4, b0y + 4, t0y + 4
        v = [b0, b1, b1y, b0y, t0, t1, t1y, t0y]
        hexes.append(f"    hex ({' '.join(map(str, v))}) ({n} 1 {nz}) simpleGrading (1 1 1)")
        q = lambda *k: "(" + " ".join(str(v[j]) for j in k) + ")"
        if i == 0:
            faces["leftWall"].append(q(0, 4, 7, 3))
        if i == len(nx) - 1:
            faces["rightWall"].append(q(1, 2, 6, 5))
        faces["frontAndBack"] += [q(0, 1, 5, 4), q(3, 7, 6, 2)]
        faces["bed"].append(q(0, 3, 2, 1))
        faces["atmosphere"].append(q(4, 5, 6, 7))
    kinds = {"leftWall": "wall", "rightWall": "wall", "bed": "wall", "atmosphere": "patch", "frontAndBack": "empty"}
    boundary = "\n".join(
        f"    {name}\n    {{\n        type {kinds[name]};\n        faces\n        (\n"
        + "\n".join(f"            {f}" for f in fs)
        + "\n        );\n    }"
        for name, fs in faces.items()
    )
    write(case, "system/blockMeshDict", f"""
scale 1;
vertices
(
{chr(10).join(f'    ({x:g} {y:g} {z:g})' for x, y, z in vertices)}
);
blocks
(
{chr(10).join(hexes)}
);
boundary
(
{boundary}
);
""")
    write(case, "system/controlDict", f"""
application     interFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {END_TIME};
deltaT          0.001;
writeControl    adjustable;
writeInterval   0.5;
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable no;
adjustTimeStep  yes;
maxCo           0.5;
maxAlphaCo      0.5;
maxDeltaT       0.05;

functions
{{
    frames
    {{
        type            writeObjects;
        libs            (utilityFunctionObjects);
        objects         (alpha.water);
        writeOption     anyWrite;
        timeStart       {FRAMES_FROM};
        writeControl    adjustableRunTime;
        writeInterval   {FRAME_EVERY};
    }}
}}
""")
    write(case, "system/fvSchemes", """
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    div(rhoPhi,U)   Gauss linearUpwind grad(U);
    div(phi,alpha)  Gauss vanLeer;
    div(phirb,alpha) Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
""")
    write(case, "system/fvSolution", """
solvers
{
    "alpha.water.*"
    {
        nAlphaCorr      2;
        nAlphaSubCycles 1;
        cAlpha          1;
        MULESCorr       yes;
        nLimiterIter    5;
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0;
    }
    "pcorr.*"   { solver PCG; preconditioner DIC; tolerance 1e-5; relTol 0; }
    p_rgh       { solver GAMG; smoother DIC; tolerance 1e-7; relTol 0.05; }
    p_rghFinal  { $p_rgh; relTol 0; }
    "U.*"       { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-6; relTol 0; }
}
PIMPLE
{
    momentumPredictor   no;
    nOuterCorrectors    1;
    nCorrectors         3;
    nNonOrthogonalCorrectors 1;
}
relaxationFactors { equations { ".*" 1; } }
""")
    write(case, "constant/g", f"dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 {-G});\n", cls="uniformDimensionedVectorField")
    write(case, "constant/transportProperties", """
phases (water air);
water { transportModel Newtonian; nu 1e-06; rho 1000; }
air   { transportModel Newtonian; nu 1.48e-05; rho 1; }
sigma 0.07;
""")
    write(case, "constant/turbulenceProperties", "simulationType laminar;\n")
    walls = ("leftWall", "rightWall", "bed")

    def field(cls, name, dims, internal, bc):
        body = f"dimensions {dims};\ninternalField {internal};\nboundaryField\n{{\n"
        for w in walls:
            body += f"    {w} {{ {bc['wall']} }}\n"
        body += f"    atmosphere {{ {bc['atmosphere']} }}\n    frontAndBack {{ type empty; }}\n}}\n"
        write(case, f"0/{name}", body, cls=cls)

    field("volVectorField", "U", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)",
          {"wall": "type slip;", "atmosphere": "type pressureInletOutletVelocity; value uniform (0 0 0);"})
    field("volScalarField", "p_rgh", "[1 -1 -2 0 0 0 0]", "uniform 0",
          {"wall": "type fixedFluxPressure; value uniform 0;", "atmosphere": "type totalPressure; p0 uniform 0;"})
    field("volScalarField", "alpha.water", "[0 0 0 0 0 0 0]", "uniform 0",
          {"wall": "type zeroGradient;", "atmosphere": "type inletOutlet; inletValue uniform 0; value uniform 0;"})
    (case / "case.json").write_text(json.dumps({
        "case": "v1-steep", "dx": dx, "blocks_x": xs, "blocks_nx": nx, "nz": nz, "top": TOP,
        "profile": PROFILE, "depth": DEPTH, "height": HEIGHT, "crest": CREST, "g": G,
        "end_time": END_TIME, "frames_from": FRAMES_FROM, "frame_every": FRAME_EVERY,
    }, indent=1))


def solitary(x):
    """The SGN solitary wave of waveslice::Solitary: surface, depth-averaged velocity and its
    first two derivatives along x (analytically, rather than waveslice's finite differences)."""
    h, a = DEPTH, HEIGHT
    kappa = np.sqrt(3 * a / (4 * h * h * (h + a)))
    c = np.sqrt(G * (h + a))
    s = 1 / np.cosh(kappa * (x - CREST))
    th = np.tanh(kappa * (x - CREST))
    eta = a * s**2
    eta_x = -2 * kappa * eta * th
    eta_xx = 2 * kappa**2 * eta * (2 * th**2 - s**2)
    big = h + eta
    ubar = c * eta / big
    # ubar = c (1 - h / H): derivatives through H = h + eta.
    ubar_x = c * h * eta_x / big**2
    ubar_xx = c * h * (eta_xx / big**2 - 2 * eta_x**2 / big**3)
    return eta, big, ubar, ubar_x, ubar_xx


def read_vectors(path: Path) -> np.ndarray:
    """An ascii volVectorField's internal values (writeCellCentres writes ascii here)."""
    text = path.read_text()
    body = text[text.index("internalField"):]
    n = int(body.split("\n")[1])
    start = body.index("(", body.index("\n", body.index("\n") + 1))
    vals = body[start + 1:].replace("(", " ").replace(")", " ").split()[: 3 * n]
    return np.array(vals, dtype=float).reshape(n, 3)


def replace_internal(path: Path, values: str) -> None:
    text = path.read_text()
    start = text.index("internalField")
    end = text.index(";", start)
    path.write_text(text[:start] + "internalField " + values + text[end:])


def init(case: Path) -> None:
    meta = json.loads((case / "case.json").read_text())
    centres = read_vectors(case / "0" / "C")
    x, z = centres[:, 0], centres[:, 2]
    zb = bed_z(x)
    dz = (TOP - zb) / meta["nz"]
    eta, big, ubar, ubar_x, ubar_xx = solitary(x)
    # Volume fraction of water in each cell, along its vertical through the centre.
    alpha = np.clip((eta - (z - dz / 2)) / dz, 0.0, 1.0)
    zeta = np.clip(z - zb, 0.0, big)  # height above the bed, held at the surface in the air
    u = ubar + (big**2 / 6 - zeta**2 / 2) * ubar_xx
    w = -zeta * ubar_x  # the bed is flat under the wave at the start
    # Air starts at rest; the water in a cell carries the velocity in proportion.
    u, w = alpha * u, alpha * w
    n = len(x)
    replace_internal(case / "0" / "alpha.water",
                     f"nonuniform List<scalar>\n{n}\n(\n" + "\n".join(f"{v:.8g}" for v in alpha) + "\n)\n")
    replace_internal(case / "0" / "U",
                     f"nonuniform List<vector>\n{n}\n(\n"
                     + "\n".join(f"({a:.8g} 0 {b:.8g})" for a, b in zip(u, w)) + "\n)\n")
    volume = float(np.sum(alpha * dz) * meta["dx"])  # m² per metre of width
    print(f"{n} cells; wave crest {eta.max():.3f} m at {x[np.argmax(eta)]:.1f} m; "
          f"max u {u.max():.3f} m/s; water area {volume:.1f} m² (still water alone "
          f"{np.trapezoid(-bed_z(np.linspace(0, 258, 2581)), dx=0.1):.1f} m²)")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("step", choices=["mesh", "init"])
    p.add_argument("--dx", type=float, default=0.4)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    if a.step == "mesh":
        mesh(a.out, a.dx)
    else:
        init(a.out)


if __name__ == "__main__":
    main()
