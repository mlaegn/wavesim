/**
# V1 as a 3D strip: what a break costs in 3D

V0 of docs/3d-box.md. The tank of `v1.c`, 4.03 m wide (258 m / 64, a cell edge of the octree's
coarse level) between free-slip walls, so the wave is the same all along the strip. It starts
from `v1.c`'s whole-tank state at 14.5 s, extruded across the strip, and runs through the break
to 17.2 s, past the lip's landing (16.95 s at 0.063 m in 2D). The set-up is `v1.c`'s.

The point is the cost: the cells an octree needs for a break in 3D and the time per cell-step.
The tube on the middle cross-section must also come out as in 2D, since nothing varies across.

Two differences from `v1.c`, both forced by Basilisk (2026-11): the embedded bed fails (SIGFPE in
the first `adapt_wavelet`) when the domain is masked in two directions, here the top and the
strip's width. So only the width is masked, and the air goes up to the cube's top at +250 m,
open there. And to keep that air affordable, no cell is larger than 1 m (level 8) rather than
0.5 m: in 2D at level 11 that moved the tube from 2.42 to 2.48 m² and the face by 0.08 s. The
cells above +8.125 m (v1.c's top) are logged apart, since a box would not have them.

    v1-3d LEVEL     # needs state.bin: v1.c's state-14.500.bin at the same LEVEL

Writes `case.json` and `frames/f-<t>.bin` (the middle cross-section, z = 2.02 m, in `v1.c`'s
format, so `cases/v1-steep/measure.py` reads it), and `log`: one line every 0.05 s with the time,
step, cells, cells below +8.125 m, water volume per metre of width and wall-clock seconds.
*/

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include <time.h>

#define G_ 9.81
#define DEPTH 8.
#define TOP 8.125 // v1.c's top: cells above it are counted apart
#define WIDTH 4.03125 // 258/64
#define LENGTH 258.
#define START 14.5
#define END_TIME 17.2
#define FRAME_EVERY 0.025
#define WX0 150.
#define WX1 258.
#define WZ0 -2.
#define WZ1 6.

int LEVEL = 10, MINLEVEL = 8;
double uemax = 0.05;
float * S = NULL;
int SNX, SNZ;
double SDX;
time_t wall0;

static const double PX[] = {0., 120., 234., 258.};
static const double PD[] = {8., 8., 0.4, 0.4};

double bed_z (double x)
{
  for (int i = 0; i < 3; i++)
    if (x <= PX[i + 1])
      return - (PD[i] + (PD[i + 1] - PD[i])*(x - PX[i])/(PX[i + 1] - PX[i]));
  return - PD[3];
}

/** Field k (0 f, 1 u, 2 w) of the 2D state, bilinear between its cell centres. */
double state (int k, double x, double y)
{
  double fx = x/SDX - 0.5, fy = (y + DEPTH)/SDX - 0.5;
  int i = clamp ((int) floor (fx), 0, SNX - 2), j = clamp ((int) floor (fy), 0, SNZ - 2);
  double a = clamp (fx - i, 0., 1.), b = clamp (fy - j, 0., 1.);
  float * F = S + (long) k*SNX*SNZ;
  return (1. - a)*(1. - b)*F[j*SNX + i] + a*(1. - b)*F[j*SNX + i + 1]
    + (1. - a)*b*F[(j + 1)*SNX + i] + a*b*F[(j + 1)*SNX + i + 1];
}

u.n[embed] = dirichlet (0.);
u.t[embed] = neumann (0.);
u.r[embed] = neumann (0.);
u.n[top] = neumann (0.);
p[top] = dirichlet (0.);
pf[top] = dirichlet (0.);

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  MINLEVEL = min (LEVEL, 8);
  SDX = LENGTH/(1 << LEVEL);
  SNX = round (LENGTH/SDX), SNZ = round ((TOP + DEPTH)/SDX);
  long n = 3L*SNX*SNZ;
  S = malloc (n*sizeof (float));
  FILE * fp = fopen ("state.bin", "rb");
  if (!fp || fread (S, sizeof (float), n, fp) != n) {
    fprintf (stderr, "state.bin: expected %ld floats (%d x %d x 3)\n", n, SNX, SNZ);
    return 1;
  }
  fclose (fp);
  size (LENGTH);
  origin (0., - DEPTH, 0.);
  rho1 = 1000., rho2 = 1.;
  mu1 = 1e-3, mu2 = 1.48e-5;
  f.sigma = 0.07;
  G.y = - G_;
  N = 64;
  CFL = 0.5;
  DT = 0.01;
  wall0 = time (NULL);
  run();
}

void set_fields ()
{
  solid (cs, fs, y - bed_z (x) - 1e-3); // 1 mm up: a bed exactly on the domain edge trips embed-tree.h in 3D
  foreach() {
    f[] = clamp (state (0, x, y), 0., 1.);
    u.x[] = cs[]*state (1, x, y);
    u.y[] = cs[]*state (2, x, y);
    u.z[] = 0.;
  }
}

event init (i = 0)
{
  mask (z > WIDTH ? front : none);
  refine (level < MINLEVEL);
  do
    set_fields();
  while (adapt_wavelet ({cs, f, u.x, u.y, u.z}, (double[]){0.01, 0.01, uemax, uemax, uemax},
			LEVEL, MINLEVEL).nf);
  set_fields();

  double dx = LENGTH/(1 << LEVEL);
  FILE * fp = fopen ("case.json", "w");
  fprintf (fp, "{\"code\": \"basilisk\", \"dims\": 3, \"level\": %d, \"dx\": %.10g, "
	   "\"width\": %g, \"section_z\": %g, "
	   "\"window\": [%g, %g, %g, %g], \"nx\": %d, \"nz\": %d, "
	   "\"profile\": [[0, 8], [120, 8], [234, 0.4], [258, 0.4]], "
	   "\"start\": %g, \"end_time\": %g, \"frames_from\": %g, \"frame_every\": %g}\n",
	   LEVEL, dx, WIDTH, WIDTH/2., WX0, WX1, WZ0, WZ1,
	   (int) round ((WX1 - WX0)/dx), (int) round ((WZ1 - WZ0)/dx),
	   START, END_TIME, START, FRAME_EVERY);
  fclose (fp);
  system ("mkdir -p frames");
}

event adapt (i++)
{
  adapt_wavelet ({cs, f, u.x, u.y, u.z}, (double[]){0.01, 0.01, uemax, uemax, uemax},
		 LEVEL, MINLEVEL);
}

event logfile (t += 0.05; t <= END_TIME - START)
{
  double vol = 0.;
  long cells = 0, below = 0;
  foreach (reduction(+:vol) reduction(+:cells) reduction(+:below)) {
    vol += f[]*dv();
    cells++;
    if (y - Delta/2. < TOP)
      below++;
  }
  FILE * fp = fopen ("log", i == 0 ? "w" : "a");
  fprintf (fp, "t %.4f i %d dt %.3g cells %ld below %ld water %.3f wall %ld\n",
	   START + t, i, dt, cells, below, vol/WIDTH, (long) (time (NULL) - wall0));
  fclose (fp);
}

/** The middle cross-section, as `v1.c` samples its window. */
event frames (t = 0; t <= END_TIME - START + 1e-6; t += FRAME_EVERY)
{
  double dx = LENGTH/(1 << LEVEL);
  int nx = round ((WX1 - WX0)/dx), nz = round ((WZ1 - WZ0)/dx);
  char path[80];
  sprintf (path, "frames/f-%.3f.bin", START + t);
  FILE * fp = fopen (path, "wb");
  float * row = malloc (nx*sizeof (float));
  for (int k = 0; k < nz; k++) {
    double zp = WZ0 + (k + 0.5)*dx;
    for (int j = 0; j < nx; j++) {
      double xp = WX0 + (j + 0.5)*dx;
      if (zp < bed_z (xp))
	row[j] = 1.;
      else {
	double v = interpolate (f, xp, zp, WIDTH/2., linear = true);
	row[j] = v == nodata ? 0. : v;
      }
    }
    fwrite (row, sizeof (float), nx, fp);
  }
  free (row);
  fclose (fp);
}
