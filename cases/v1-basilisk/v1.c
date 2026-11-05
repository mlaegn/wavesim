/**
# V1 of docs/3d-box.md in Basilisk: the `steep` slice case

The same tank and wave as `cases/v1-steep/` (OpenFOAM) and `wavesim slice steep`: water 8 m deep,
flat for 120 m, a plane slope of 1 in 15 up to a shelf 0.4 m deep and 24 m long, walls at both
ends, air up to +8.125 m with an open top (OpenFOAM's top was at +6 m; Basilisk's mask cuts the
square domain along cell edges, and the nearest one is 2.1 m higher). A solitary wave 2.4 m high, crest 60 m from the left wall,
from the Serre-Green-Naghdi theory, with that theory's velocity profile through the water column.
Air and water as in the OpenFOAM case (rho 1000 and 1, nu 1e-6 and 1.48e-5, sigma 0.07), laminar,
free-slip walls and bed.

The set-up follows Feddersen et al. (2024): two-phase Navier-Stokes with momentum-conserving VOF
(`navier-stokes/conserving.h`), reduced gravity, the bed as an embedded boundary, and a quadtree
that refines where the surface, the bed and the flow need it, down to 258 m / 2^LEVEL. No
cell is larger than 0.5 m (level 9): with larger ones, the water column under the wave was a
single cell deep, the bed's kink at the slope's toe added 20 m² of water, and the wave slowed.

    v1 LEVEL     # 10: 0.25 m, 11: 0.126 m, 12: 0.063 m, 13: 0.031 m

Writes, in the current folder:

- `case.json`, what `cases/v1-steep/measure.py` needs to read the frames;
- `frames/f-<t>.bin`: f (water fraction) sampled on a uniform grid at the finest cell size over
  x 150-258 m, z -2-6 m, float32, z-major, every 0.025 s from 13 s; the bed counts as water there,
  so that the only 0.5 contour is the free surface;
- `frames/u-<t>.bin`: u and w on the same grid at a few times, for the air-speed check;
- `log`: one line every 0.05 s (time, step, cells, water area, highest surface).
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"

#define G_ 9.81
#define DEPTH 8.
#define HEIGHT 2.4
#define CREST 60.
#define TOP 8.125 // -8 + 3 x 258/64
#define LENGTH 258.
#define END_TIME 18.5
#define FRAMES_FROM 13.
#define FRAME_EVERY 0.025
#define WX0 150.
#define WX1 258.
#define WZ0 -2.
#define WZ1 6.

int LEVEL = 10, MINLEVEL = 9;
double uemax = 0.05; // m/s, the velocity error that refines the mesh

/** The steep profile, as (x, still depth). */
static const double PX[] = {0., 120., 234., 258.};
static const double PD[] = {8., 8., 0.4, 0.4};

double bed_z (double x)
{
  for (int i = 0; i < 3; i++)
    if (x <= PX[i + 1])
      return - (PD[i] + (PD[i + 1] - PD[i])*(x - PX[i])/(PX[i + 1] - PX[i]));
  return - PD[3];
}

/** The SGN solitary wave, as in `make_case.py`'s `solitary()`. */
void solitary (double x, double * eta, double * big, double * ubar,
	       double * ubar_x, double * ubar_xx)
{
  double h = DEPTH, a = HEIGHT;
  double kappa = sqrt (3.*a/(4.*h*h*(h + a)));
  double c = sqrt (G_*(h + a));
  double s = 1./cosh (kappa*(x - CREST)), th = tanh (kappa*(x - CREST));
  *eta = a*s*s;
  double eta_x = -2.*kappa*(*eta)*th;
  double eta_xx = 2.*kappa*kappa*(*eta)*(2.*th*th - s*s);
  *big = h + *eta;
  *ubar = c*(*eta)/(*big);
  *ubar_x = c*h*eta_x/sq(*big);
  *ubar_xx = c*h*(eta_xx/sq(*big) - 2.*sq(eta_x)/cube(*big));
}

/** Walls left and right and at the bottom (free-slip, Basilisk's default), the bed free-slip,
an open top. */
u.n[embed] = dirichlet (0.);
u.t[embed] = neumann (0.);
u.n[top] = neumann (0.);
p[top] = dirichlet (0.);
pf[top] = dirichlet (0.);

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  MINLEVEL = min (LEVEL, 9);
  size (LENGTH);
  origin (0., - DEPTH);
  rho1 = 1000., rho2 = 1.;
  mu1 = 1e-3, mu2 = 1.48e-5;
  f.sigma = 0.07;
  G.y = - G_;
  N = 64;
  CFL = 0.5;
  DT = 0.01;
  run();
}

double surface (double x)
{
  double eta, big, ub, ubx, ubxx;
  solitary (x, &eta, &big, &ub, &ubx, &ubxx);
  return eta;
}

void set_fields ()
{
  solid (cs, fs, y - bed_z (x));
  fraction (f, surface (x) - y);
  foreach() {
    double eta, big, ub, ubx, ubxx;
    solitary (x, &eta, &big, &ub, &ubx, &ubxx);
    double zeta = clamp (y - bed_z (x), 0., big);
    u.x[] = f[]*cs[]*(ub + (sq(big)/6. - sq(zeta)/2.)*ubxx);
    u.y[] = f[]*cs[]*(- zeta*ubx);
  }
}

event init (i = 0)
{
  mask (y > TOP ? top : none); // at N = 64: cells [4.09, 8.125] stay, [8.125, 12.16] go
  refine (level < MINLEVEL);
  do
    set_fields();
  while (adapt_wavelet ({cs, f, u.x, u.y}, (double[]){0.01, 0.01, uemax, uemax},
			LEVEL, MINLEVEL).nf);
  set_fields();

  double dx = LENGTH/(1 << LEVEL);
  FILE * fp = fopen ("case.json", "w");
  fprintf (fp, "{\"code\": \"basilisk\", \"level\": %d, \"dx\": %.10g, "
	   "\"window\": [%g, %g, %g, %g], \"nx\": %d, \"nz\": %d, "
	   "\"profile\": [[0, 8], [120, 8], [234, 0.4], [258, 0.4]], "
	   "\"end_time\": %g, \"frames_from\": %g, \"frame_every\": %g}\n",
	   LEVEL, dx, WX0, WX1, WZ0, WZ1,
	   (int) round ((WX1 - WX0)/dx), (int) round ((WZ1 - WZ0)/dx),
	   END_TIME, FRAMES_FROM, FRAME_EVERY);
  fclose (fp);
  system ("mkdir -p frames");
}

/** The adaptive mesh follows the surface, the bed and the flow. */
event adapt (i++)
{
  adapt_wavelet ({cs, f, u.x, u.y}, (double[]){0.01, 0.01, uemax, uemax}, LEVEL, MINLEVEL);
}

event logfile (t += 0.05; t <= END_TIME)
{
  double area = 0., top = - HUGE;
  long cells = 0;
  foreach (reduction(+:area) reduction(+:cells) reduction(max:top)) {
    area += f[]*dv();
    cells++;
    if (f[] > 0.01 && f[] < 0.99 && y > top)
      top = y;
  }
  FILE * fp = fopen ("log", i == 0 ? "w" : "a");
  fprintf (fp, "t %.4f i %d dt %.3g cells %ld water %.3f top %.3f\n",
	   t, i, dt, cells, area, top);
  fclose (fp);
}

/** Samples on the window's uniform grid, at cell centres of the finest level. The bed counts
as water. */
void sample (const char * name, scalar * list, double t)
{
  double dx = LENGTH/(1 << LEVEL);
  int nx = round ((WX1 - WX0)/dx), nz = round ((WZ1 - WZ0)/dx);
  char path[80];
  sprintf (path, "frames/%s-%.3f.bin", name, t);
  FILE * fp = fopen (path, "wb");
  float * row = malloc (nx*sizeof (float));
  for (scalar s in list)
    for (int k = 0; k < nz; k++) {
      double zp = WZ0 + (k + 0.5)*dx;
      for (int j = 0; j < nx; j++) {
	double xp = WX0 + (j + 0.5)*dx;
	if (zp < bed_z (xp))
	  row[j] = s.i == f.i ? 1. : 0.;
	else {
	  double v = interpolate (s, xp, zp, linear = true);
	  row[j] = v == nodata ? 0. : v;
	}
      }
      fwrite (row, sizeof (float), nx, fp);
    }
  free (row);
  fclose (fp);
}

event frames (t = FRAMES_FROM; t <= END_TIME + 1e-6; t += FRAME_EVERY)
{
  sample ("f", {f}, t);
}

event velocities (t = {10., 14.5, 15.5, 16.})
{
  sample ("u", {u.x, u.y}, t);
}
