/**
# A made-up peel: V1's wave on V1's slope, with the depth lines turned

Does a peeling break make a bigger tube than the same wave breaking straight on? This is the
question of the 3D box (docs/3d-box.md), asked without Pipeline's reef or wavesim's swell: the
wave, slope, fluids and set-up of `v1.c`, in a strip 16.1 m wide whose bed is V1's profile
shifted shoreward by z tan(ANGLE). The depth lines then lie at ANGLE to the crest, so the wave
breaks first at z = 0 and the break runs across the strip. The strip's sides are free-slip
walls, mirrors: at z = 0 the bed is a peak, the break peeling away from it both ways.

The shift starts shoreward of the wave as it is at 12 s (its front 5 cm high at about 205 m):
it ramps from nothing at 206 m to full at 216 m. So the 2D run's state at 12 s, extruded across
the strip, is still right for the new bed, and the run starts from it. ANGLE 0 is the straight
control: V1 in 3D.

Set-up as `v1-3d.c` (V0), with three differences, all for cost:

- the domain is a cube 129 m on a side, from x = 129 m to the shore wall at 258 m, so that its
  cells are the 2D run's at one level more: level 9 here is level 10 there (0.252 m), level 10
  is level 11 (0.126 m);
- no cell is larger than 1 m (level 7) below +4 m, as in V0; above, the air is left as coarse
  as the mesh allows, unless the flow refines it;
- the side mask is the only mask, as V0 found it must be: the top is the cube's, at +121 m.

    v1-peel LEVEL ANGLE    # needs state.bin: v1.c's state-12.000.bin at level LEVEL + 1

Writes `case.json`, and for each of three cross-sections, at z = 1, 8.06 and 15.1 m,
`sec-<z>/case.json` and `sec-<z>/frames/f-<t>.bin` in `v1.c`'s format, so that
`cases/v1-steep/measure.py out/<run>/sec-<z>` measures it as it measures a 2D run. Also
`facets-<t>.txt`, the surface as triangles at a few moments, and `log`: one line every 0.05 s
with the time, step, cells, water volume per metre of width and wall-clock seconds.
*/

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "fractions.h"
#include <time.h>

#define G_ 9.81
#define DEPTH 8.
#define TOP2D 8.125 // v1.c's top: the 2D state reaches up to here
#define LENGTH 258. // the 2D tank
#define X0 129.
#define SIZE 129.
#define WIDTH 16.125 // SIZE/8, a cell edge of the octree's level 3
#define FILL 4. // no cell larger than 1 m below this height
#define RAMP0 206.
#define RAMP1 216.
#define START 12.
#define END_TIME 18.
#define FRAMES_FROM 15.
#define FRAME_EVERY 0.025
#define WX0 150.
#define WX1 258.
#define WZ0 -2.
#define WZ1 6.

int LEVEL = 9, MINLEVEL = 7;
double ANGLE = 20.; // degrees
double uemax = 0.05;
float * S = NULL;
int SNX, SNZ;
double SDX;
time_t wall0;
const double SECTIONS[] = {1., WIDTH/2., WIDTH - 1.};

static const double PX[] = {0., 120., 234., 258.};
static const double PD[] = {8., 8., 0.4, 0.4};

double profile_z (double x)
{
  for (int i = 0; i < 3; i++)
    if (x <= PX[i + 1])
      return - (PD[i] + (PD[i + 1] - PD[i])*(x - PX[i])/(PX[i + 1] - PX[i]));
  return - PD[3];
}

/** V1's profile, shifted shoreward by z tan(ANGLE) beyond the ramp. */
double bed_z (double x, double z)
{
  double r = clamp ((x - RAMP0)/(RAMP1 - RAMP0), 0., 1.);
  r = r*r*(3. - 2.*r);
  return profile_z (x - r*z*tan (ANGLE*pi/180.));
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
  if (argc > 2)
    ANGLE = atof (argv[2]);
  MINLEVEL = min (LEVEL, 7);
  SDX = LENGTH/(1 << (LEVEL + 1));
  SNX = round (LENGTH/SDX), SNZ = round ((TOP2D + DEPTH)/SDX);
  long n = 3L*SNX*SNZ;
  S = malloc (n*sizeof (float));
  FILE * fp = fopen ("state.bin", "rb");
  if (!fp || fread (S, sizeof (float), n, fp) != n) {
    fprintf (stderr, "state.bin: expected %ld floats (%d x %d x 3)\n", n, SNX, SNZ);
    return 1;
  }
  fclose (fp);
  size (SIZE);
  origin (X0, - DEPTH, 0.);
  rho1 = 1000., rho2 = 1.;
  mu1 = 1e-3, mu2 = 1.48e-5;
  f.sigma = 0.07;
  G.y = - G_;
  N = 8;
  CFL = 0.5;
  DT = 0.01;
  wall0 = time (NULL);
  run();
}

void set_fields ()
{
  solid (cs, fs, y - bed_z (x, z) - 1e-3); // 1 mm up: a bed exactly on the domain edge trips embed-tree.h in 3D
  foreach() {
    f[] = clamp (state (0, x, y), 0., 1.);
    u.x[] = cs[]*state (1, x, y);
    u.y[] = cs[]*state (2, x, y);
    u.z[] = 0.;
  }
}

void write_case (const char * path, double zs)
{
  double dx = SIZE/(1 << LEVEL);
  FILE * fp = fopen (path, "w");
  fprintf (fp, "{\"code\": \"basilisk\", \"dims\": 3, \"level\": %d, \"dx\": %.10g, "
	   "\"angle\": %g, \"width\": %g, \"section_z\": %g, "
	   "\"window\": [%g, %g, %g, %g], \"nx\": %d, \"nz\": %d, \"profile\": [",
	   LEVEL, dx, ANGLE, WIDTH, zs, WX0, WX1, WZ0, WZ1,
	   (int) round ((WX1 - WX0)/dx), (int) round ((WZ1 - WZ0)/dx));
  // the bed under this section, as (x, still depth) every metre
  for (int i = 0; i <= 258; i++)
    fprintf (fp, "%s[%d, %.4f]", i ? ", " : "", i, - bed_z (i, zs < 0. ? 0. : zs));
  fprintf (fp, "], \"start\": %g, \"end_time\": %g, \"frames_from\": %g, \"frame_every\": %g}\n",
	   START, END_TIME, FRAMES_FROM, FRAME_EVERY);
  fclose (fp);
}

event init (i = 0)
{
  mask (z > WIDTH ? front : none);
  refine (level < MINLEVEL && y - Delta/2. < FILL);
  fprintf (stderr, "init: %ld cells after the fill, %lds\n", grid->tn, (long) (time (NULL) - wall0));
  astats s;
  int pass = 0;
  do { // it can end up flipping a few cells back and forth, hence the cap
    set_fields();
    s = adapt_wavelet ({cs, f, u.x, u.y, u.z}, (double[]){0.01, 0.01, uemax, uemax, uemax},
		       LEVEL, MINLEVEL);
    fprintf (stderr, "init: %d refined, %d coarsened, %ld cells, %lds\n",
	     s.nf, s.nc, grid->tn, (long) (time (NULL) - wall0));
  } while (s.nf > 100 && ++pass < 30);
  set_fields();

  /** Where the cells are, by level and by region. */
  long by_level[20] = {0}, air = 0, offshore = 0, water = 0, bed = 0;
  foreach (serial) {
    by_level[level]++;
    if (y - Delta/2. > FILL) air++;
    else if (x < RAMP0 - 20.) offshore++;
    else if (cs[] > 0. && cs[] < 1.) bed++;
    else water++;
  }
  for (int l = 0; l < 20; l++)
    if (by_level[l])
      fprintf (stderr, "init: level %d (%.3g m): %ld cells\n", l, SIZE/(1 << l), by_level[l]);
  fprintf (stderr, "init: above %g m %ld, x < %g m %ld, at the bed %ld, the rest %ld\n",
	   FILL, air, RAMP0 - 20., offshore, bed, water);
  if (getenv ("INIT_ONLY"))
    exit (0);

  write_case ("case.json", -1.);
  for (int k = 0; k < 3; k++) {
    char cmd[80], path[80];
    sprintf (cmd, "mkdir -p sec-%.2f/frames", SECTIONS[k]);
    system (cmd);
    sprintf (path, "sec-%.2f/case.json", SECTIONS[k]);
    write_case (path, SECTIONS[k]);
  }
}

/** Basilisk traps division by zero and invalid operations. At 20° something traps at 14.25 s
(SIGFPE, exit 136) without a NaN ever reaching a field: with the trap off, the run goes on and
conserves its water. So the trap is off, and instead the fields are checked after each step; a NaN
in one prints the first cells that have it and stops the run (exit 5). */
event nancheck (i++)
{
  if (i == 0)
    disable_fpe (FE_DIVBYZERO|FE_INVALID);
  long bad = 0;
  foreach (reduction(+:bad))
    if (isnan (f[]) || isnan (u.x[]) || isnan (u.y[]) || isnan (u.z[]) || isnan (p[])
	|| isnan (cs[]) || isnan (fs.x[]) || isnan (fs.y[]) || isnan (fs.z[]))
      bad++;
  if (!bad)
    return 0;
  fprintf (stderr, "nancheck: t %.4f i %d: %ld cells with a NaN\n", START + t, i, bad);
  int shown = 0;
  foreach (serial)
    if (shown < 30 && (isnan (f[]) || isnan (u.x[]) || isnan (u.y[]) || isnan (u.z[])
		       || isnan (p[]) || isnan (cs[]) || isnan (fs.x[]) || isnan (fs.y[])
		       || isnan (fs.z[]))) {
      fprintf (stderr, "  x %.3f y %.3f z %.3f level %d cs %g fs %g %g %g f %g u %g %g %g p %g "
	       "bed %.3f\n", x, y, z, level, cs[], fs.x[], fs.y[], fs.z[], f[], u.x[], u.y[],
	       u.z[], p[], bed_z (x, z));
      shown++;
    }
  exit (5);
}

event adapt (i++)
{
  adapt_wavelet ({cs, f, u.x, u.y, u.z}, (double[]){0.01, 0.01, uemax, uemax, uemax},
		 LEVEL, MINLEVEL);
}

event logfile (t += 0.05; t <= END_TIME - START)
{
  double vol = 0.;
  long cells = 0;
  foreach (reduction(+:vol) reduction(+:cells)) {
    vol += f[]*cs[]*dv();
    cells++;
  }
  FILE * fp = fopen ("log", i == 0 ? "w" : "a");
  fprintf (fp, "t %.4f i %d dt %.3g cells %ld water %.3f wall %ld\n",
	   START + t, i, dt, cells, vol/WIDTH, (long) (time (NULL) - wall0));
  fclose (fp);
}

/** Each cross-section, as `v1.c` samples its window; below the bed counts as water. */
event frames (t = FRAMES_FROM - START; t <= END_TIME - START + 1e-6; t += FRAME_EVERY)
{
  double dx = SIZE/(1 << LEVEL);
  int nx = round ((WX1 - WX0)/dx), nz = round ((WZ1 - WZ0)/dx);
  float * row = malloc (nx*sizeof (float));
  for (int s = 0; s < 3; s++) {
    double zs = SECTIONS[s];
    char path[80];
    sprintf (path, "sec-%.2f/frames/f-%.3f.bin", zs, START + t);
    FILE * fp = fopen (path, "wb");
    for (int k = 0; k < nz; k++) {
      double yp = WZ0 + (k + 0.5)*dx;
      for (int j = 0; j < nx; j++) {
	double xp = WX0 + (j + 0.5)*dx;
	if (yp < bed_z (xp, zs))
	  row[j] = 1.;
	else {
	  double v = interpolate (f, xp, yp, zs, linear = true);
	  row[j] = v == nodata ? 0. : v;
	}
      }
      fwrite (row, sizeof (float), nx, fp);
    }
    fclose (fp);
  }
  free (row);
}

/** The surface in 3D at a few moments of the break, for a picture of the peel. */
event facets (t = {16. - START, 16.4 - START, 16.8 - START, 17.2 - START, 17.6 - START})
{
  char path[80];
  sprintf (path, "facets-%.3f.txt", START + t);
  FILE * fp = fopen (path, "w");
  output_facets (f, fp);
  fclose (fp);
}
