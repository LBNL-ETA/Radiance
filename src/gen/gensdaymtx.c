#ifndef lint
static const char RCSid[] = "$Id: gensdaymtx.c,v 1.4 2024/08/08 02:00:20 greg Exp $";
#endif
#include "atmos.h"
#include "copyright.h"
#include "data.h"
#include "platform.h"
#include "rtio.h"
#include <ctype.h>
#include <stdlib.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

char *progname;

double altitude;   /* Solar altitude (radians) */
double azimuth;    /* Solar azimuth (radians) */
int julian_date;   /* Julian date */
double sun_zenith; /* Sun zenith angle (radians) */
int input = 0;     /* Input type */
int output = 0;    /* Output type */
FVECT sundir;

const double ARCTIC_LAT = 67.;
const double TROPIC_LAT = 23.;
const int SUMMER_START = 4;
const int SUMMER_END = 9;
const double GNORM = 0.777778;

const double D65EFF = 203.; /* standard illuminant D65 */

/* Mean normalized relative daylight spectra where CCT = 6415K for overcast */
const double D6415[NSSAMP] = {0.63231, 1.06171, 1.00779, 1.36423, 1.34133,
                              1.27258, 1.26276, 1.26352, 1.22201, 1.13246,
                              1.0434,  1.05547, 0.98212, 0.94445, 0.9722,
                              0.82387, 0.87853, 0.82559, 0.75111, 0.78925};
/* Degrees into radians */
#define DegToRad(deg) ((deg) * (PI / 180.))

/* Radiuans into degrees */
#define RadToDeg(rad) ((rad) * (180. / PI))

#ifndef NSUNPATCH
#define NSUNPATCH 4 /* max. # patches to spread sun into */
#endif

#define SUN_ANG_DEG 0.533 /* sun full-angle in degrees */

int nsuns = NSUNPATCH;    /* number of sun patches to use */
double fixed_sun_sa = -1; /* fixed solid angle per sun? */

int verbose = 0; /* progress reports to stderr? */

int outfmt = 'a'; /* output format */

int rhsubdiv = 1; /* Reinhart sky subdivisions */

COLOR skycolor = {.96, 1.004, 1.118}; /* sky coloration */
COLOR suncolor = {1., 1., 1.};        /* sun color */
double grefl = .2;                    /* ground reflectance */

int nskypatch;  /* number of Reinhart patches */
float *rh_palt; /* sky patch altitudes (radians) */
float *rh_pazi; /* sky patch azimuths (radians) */
float *rh_dom;  /* sky patch solid angle (sr) */

double sun_ct;

#define vector(v, alt, azi)                                                    \
  ((v)[1] = cos(alt), (v)[0] = (v)[1] * sin(azi), (v)[1] *= cos(azi),          \
   (v)[2] = sin(alt))

#define rh_vector(v, i) vector(v, rh_palt[i], rh_pazi[i])

#define rh_cos(i) tsin(rh_palt[i])

#define solar_minute(jd, hr) ((24 * 60) * ((jd) - 1) + (int)((hr) * 60. + .5))

inline void vectorize(double altitude, double azimuth, FVECT v) {
  v[1] = cos(altitude);
  v[0] = (v)[1] * sin(azimuth);
  v[1] *= cos(azimuth);
  v[2] = sin(altitude);
}

static int make_directory(const char *path) {
#ifdef _WIN32
  if (CreateDirectory(path, NULL) || GetLastError() == ERROR_ALREADY_EXISTS) {
    return 1;
  }
  return 0;
#else
  if (mkdir(path, 0777) == 0 || errno == EEXIST) {
    return 1;
  }
  return 0;
#endif
}

static const char *getfmtname(int fmt) {
  switch (fmt) {
  case 'a':
    return ("ascii");
  case 'f':
    return ("float");
  case 'd':
    return ("double");
  }
  return ("unknown");
}

static inline double wmean2(const double a, const double b, const double x) {
  return a * (1 - x) + b * x;
}

static inline double wmean(const double a, const double x, const double b,
                           const double y) {
  return (a * x + b * y) / (a + b);
}

static double get_zenith_brightness(const double sundir[3]) {
  double zenithbr;
  if (sundir[2] < 0) {
    zenithbr = 0;
  } else {
    zenithbr = (8.6 * sundir[2] + .123) * 1000.0 / D65EFF;
  }
  return zenithbr;
}

/* from gensky.c */
static double get_overcast_brightness(const double dz, const double zenithbr) {
  double groundbr = zenithbr * GNORM;
  return wmean(pow(dz + 1.01, 10), zenithbr * (1 + 2 * dz) / 3,
               pow(dz + 1.01, -10), groundbr);
}

int rh_init(void) {
#define NROW 7
  static const int tnaz[NROW] = {30, 30, 24, 24, 18, 12, 6};
  const double alpha = (PI / 2.) / (NROW * rhsubdiv + .5);
  int p, i, j;
  /* allocate patch angle arrays */
  nskypatch = 0;
  for (p = 0; p < NROW; p++)
    nskypatch += tnaz[p];
  nskypatch *= rhsubdiv * rhsubdiv;
  nskypatch += 2;
  rh_palt = (float *)malloc(sizeof(float) * nskypatch);
  rh_pazi = (float *)malloc(sizeof(float) * nskypatch);
  rh_dom = (float *)malloc(sizeof(float) * nskypatch);
  if ((rh_palt == NULL) | (rh_pazi == NULL) | (rh_dom == NULL)) {
    fprintf(stderr, "%s: out of memory in rh_init()\n", progname);
    exit(1);
  }
  rh_palt[0] = -PI / 2.; /* ground & zenith patches */
  rh_pazi[0] = 0.;
  rh_dom[0] = 2. * PI;
  rh_palt[nskypatch - 1] = PI / 2.;
  rh_pazi[nskypatch - 1] = 0.;
  rh_dom[nskypatch - 1] = 2. * PI * (1. - cos(alpha * .5));
  p = 1; /* "normal" patches */
  for (i = 0; i < NROW * rhsubdiv; i++) {
    const float ralt = alpha * (i + .5);
    const int ninrow = tnaz[i / rhsubdiv] * rhsubdiv;
    const float dom =
        2. * PI * (sin(alpha * (i + 1)) - sin(alpha * i)) / (double)ninrow;
    for (j = 0; j < ninrow; j++) {
      rh_palt[p] = ralt;
      rh_pazi[p] = 2. * PI * j / (double)ninrow;
      rh_dom[p++] = dom;
    }
  }
  return nskypatch;
#undef NROW
}

/* Resize daylight matrix (GW) */
float *resize_dmatrix(float *mtx_data, int nsteps, int npatch) {
  if (mtx_data == NULL)
    mtx_data = (float *)malloc(sizeof(float) * NSSAMP * nsteps * npatch);
  else
    mtx_data =
        (float *)realloc(mtx_data, sizeof(float) * NSSAMP * nsteps * npatch);
  if (mtx_data == NULL) {
    fprintf(stderr, "%s: out of memory in resize_dmatrix(%d,%d)\n", progname,
            nsteps, npatch);
    exit(1);
  }
  return (mtx_data);
}

static Atmosphere init_atmos(const double aod, const double grefl) {
  Atmosphere atmos = {.ozone_density = {.layers =
                                            {
                                                {.width = 25000.0,
                                                 .exp_term = 0.0,
                                                 .exp_scale = 0.0,
                                                 .linear_term = 1.0 / 15000.0,
                                                 .constant_term = -2.0 / 3.0},
                                                {.width = AH,
                                                 .exp_term = 0.0,
                                                 .exp_scale = 0.0,
                                                 .linear_term = -1.0 / 15000.0,
                                                 .constant_term = 8.0 / 3.0},
                                            }},
                      .rayleigh_density = {.layers =
                                               {
                                                   {.width = AH,
                                                    .exp_term = 1.0,
                                                    .exp_scale = -1.0 / HR_MS,
                                                    .linear_term = 0.0,
                                                    .constant_term = 0.0},
                                               }},
                      .beta_r0 = BR0_MS,
                      .beta_scale = aod / AOD0_CA,
                      .beta_m = NULL,
                      .grefl = grefl};
  return atmos;
}

static DpPaths get_dppaths(const char *dir, const double aod, const char *mname,
                           const char *tag) {
  DpPaths paths;

  snprintf(paths.tau, PATH_MAX, "%s%ctau_%s_%s_%.2f.dat", dir, DIRSEP, tag,
           mname, aod);
  snprintf(paths.scat, PATH_MAX, "%s%cscat_%s_%s_%.2f.dat", dir, DIRSEP, tag,
           mname, aod);
  snprintf(paths.scat1m, PATH_MAX, "%s%cscat1m_%s_%s_%.2f.dat", dir, DIRSEP,
           tag, mname, aod);
  snprintf(paths.irrad, PATH_MAX, "%s%cirrad_%s_%s_%.2f.dat", dir, DIRSEP, tag,
           mname, aod);

  return paths;
}
static void set_rayleigh_density_profile(Atmosphere *atmos, char *tag,
                                         const int is_summer,
                                         const double s_latitude) {
  /* Set rayleigh density profile */
  if (fabs(s_latitude * 180.0 / PI) > ARCTIC_LAT) {
    tag[0] = 's';
    if (is_summer) {
      tag[1] = 's';
      atmos->rayleigh_density.layers[0].exp_scale = -1.0 / HR_SS;
      atmos->beta_r0 = BR0_SS;
    } else {
      tag[1] = 'w';
      atmos->rayleigh_density.layers[0].exp_scale = -1.0 / HR_SW;
      atmos->beta_r0 = BR0_SW;
    }
  } else if (fabs(s_latitude * 180.0 / PI) > TROPIC_LAT) {
    tag[0] = 'm';
    if (is_summer) {
      tag[1] = 's';
      atmos->rayleigh_density.layers[0].exp_scale = -1.0 / HR_MS;
      atmos->beta_r0 = BR0_MS;
    } else {
      tag[1] = 'w';
      atmos->rayleigh_density.layers[0].exp_scale = -1.0 / HR_MW;
      atmos->beta_r0 = BR0_MW;
    }
  } else {
    tag[0] = 't';
    tag[1] = 'r';
    atmos->rayleigh_density.layers[0].exp_scale = -1.0 / HR_T;
    atmos->beta_r0 = BR0_T;
  }
  tag[2] = '\0';
}
/* Add in solar direct to nearest sky patches (GW) */
void add_direct(DATARRAY *tau, DATARRAY *scat, DATARRAY *scat1m,
                DATARRAY *irrad, double ccover, float *parr) {
  FVECT svec;
  double near_dprod[NSUNPATCH];
  int near_patch[NSUNPATCH];
  double wta[NSUNPATCH], wtot;
  int i, j, p;

  /* identify nsuns closest patches */
  for (i = nsuns; i--;)
    near_dprod[i] = -1.;
  vectorize(altitude, azimuth, svec);
  for (p = 1; p < nskypatch; p++) {
    FVECT pvec;
    double dprod;
    vectorize(rh_palt[p], rh_pazi[p], pvec);
    dprod = DOT(pvec, svec);
    for (i = 0; i < nsuns; i++)
      if (dprod > near_dprod[i]) {
        for (j = nsuns; --j > i;) {
          near_dprod[j] = near_dprod[j - 1];
          near_patch[j] = near_patch[j - 1];
        }
        near_dprod[i] = dprod;
        near_patch[i] = p;
        break;
      }
  }
  /* Get solar radiance */
  double sun_radiance[NSSAMP] = {0};
  get_solar_radiance(tau, scat, scat1m, sundir, ER, sun_ct, sun_radiance);
  if (ccover > 0) {
    double zenithbr = get_zenith_brightness(sundir);
    double skybr = get_overcast_brightness(sundir[2], zenithbr);
    int l;
    for (l = 0; l < NSSAMP; ++l) {
      sun_radiance[l] =
          wmean2(sun_radiance[l], D6415[l] * skybr / WVLSPAN, ccover);
    }
  }
  /* weight by proximity */
  wtot = 0;
  for (i = nsuns; i--;)
    wtot += wta[i] = 1. / (1.002 - near_dprod[i]);
  /* add to nearest patch radiances */
  for (i = nsuns; i--;) {
    float *pdest = parr + NSSAMP * near_patch[i];
    int k;
    for (k = 0; k < NSSAMP; k++) {
      *pdest++ = sun_radiance[k] * wta[i] / wtot;
    }
  }
}

void calc_sky_patch_radiance(DATARRAY *scat, DATARRAY *scat1m, double ccover,
                             float *parr) {
  int i;
  double mu_sky; /* Sun-sky point azimuthal angle */
  double sspa;   /* Sun-sky point angle */
  FVECT view_point = {0, 0, ER};
  for (i = 1; i < nskypatch; i++) {
    FVECT rdir_sky;
    int k;
    vectorize(rh_palt[i], rh_pazi[i], rdir_sky);
    mu_sky = fdot(view_point, rdir_sky) / ER;
    sspa = fdot(rdir_sky, sundir);
    SCOLOR sky_radiance = {0};

    get_sky_radiance(scat, scat1m, ER, mu_sky, sun_ct, sspa, sky_radiance);
    for (k = 0; k < NSSAMP; ++k) {
      sky_radiance[k] *= WVLSPAN;
    }

    if (ccover > 0) {
      double zenithbr = get_zenith_brightness(sundir);
      double grndbr = zenithbr * GNORM;
      double skybr = get_overcast_brightness(rdir_sky[2], zenithbr);
      int k;
      for (k = 0; k < NSSAMP; ++k) {
        sky_radiance[k] = wmean2(sky_radiance[k], skybr * D6415[k], ccover);
      }
    }

    for (k = 0; k < NSSAMP; ++k) {
      parr[NSSAMP * i + k] = sky_radiance[k];
    }
  }
}

/* Return maximum of two doubles */
static inline double dmax(double a, double b) { return (a > b) ? a : b; }

/* Compute sky patch radiance values (modified by GW) */
void compute_sky(DATARRAY *tau, DATARRAY *scat, DATARRAY *scat1m,
                 DATARRAY *irrad, double ccover, float *parr) {
  int index; /* Category index */
  int i;
  float sun_zenith;
  SCOLOR sky_radiance = {0};
  SCOLOR ground_radiance = {0};
  SCOLR sky_sclr = {0};
  SCOLR ground_sclr = {0};
  FVECT view_point = {0, 0, ER};
  const double radius = VLEN(view_point);
  const double sun_ct = fdot(view_point, sundir) / radius;
  const FVECT rdir_grnd = {0, 0, -1};
  const double mu_grnd = fdot(view_point, rdir_grnd) / radius;
  const double nu_grnd = fdot(rdir_grnd, sundir);
  int j;

  /* Calculate sun zenith angle (don't let it dip below horizon) */
  /* Also limit minimum angle to keep circumsolar off zenith */
  if (altitude <= 0.0)
    sun_zenith = DegToRad(90.0);
  else if (altitude >= DegToRad(87.0))
    sun_zenith = DegToRad(3.0);
  else
    sun_zenith = DegToRad(90.0) - altitude;

  /* Compute ground radiance (include solar contribution if any) */
  get_ground_radiance(tau, scat, scat1m, irrad, view_point, rdir_grnd, radius,
                      mu_grnd, sun_ct, nu_grnd, grefl, sundir, parr);
  for (j = 0; j < NSSAMP; j++) {
    parr[j] *= WVLSPAN;
  }
  calc_sky_patch_radiance(scat, scat1m, ccover, parr);
}

int main(int argc, char *argv[]) {

  char buf[256];
  int doheader = 1; /* output header? */
  double rotation = 0.0;
  double elevation = 0;
  int leap_day = 0;       /* add leap day? */
  int sun_hours_only = 0; /* only output sun hours? */
  float *mtx_data = NULL;
  int ntsteps = 0;      /* number of time steps */
  int tstorage = 0;     /* number of allocated time steps */
  int nstored = 0;      /* number of time steps in matrix */
  int last_monthly = 0; /* month of last report */
  int mo, da;
  double hr, aod, cc;
  double dni, dhi;
  int mtx_offset = 0;
  int i, j;
  char lstag[3];
  char *mie_path = getpath("mie_ca.dat", getrlibpath(), R_OK);
  char *ddir = ".";
  char mie_name[20] = "mie_ca";
  int num_threads = 1;
  int sorder = 4;
  int solar_only = 0;
  int sky_only = 0;
  FVECT view_point = {0, 0, ER};

  progname = argv[0];

  for (i = 1; i < argc && argv[i][0] == '-'; i++) {
    switch (argv[i][1]) {
    case 'd': /* solar (direct) only */
      solar_only = 1;
      break;
    case 's': /* sky only (no direct) */
      sky_only = 1;
      break;
    case 'g':
      grefl = atof(argv[++i]);
      break;
    case 'm':
      rhsubdiv = atoi(argv[++i]);
      break;
    case 'n':
      num_threads = atoi(argv[++i]);
      break;
    case 'r': /* rotate distribution */
      if (argv[i][2] && argv[i][2] != 'z')
        goto userr;
      rotation = atof(argv[++i]);
      break;
    case 'u': /* solar hours only */
      sun_hours_only = 1;
      break;
    case 'p':
      ddir = argv[++i];
      break;
    case 'v': /* verbose progress reports */
      verbose++;
      break;
    case 'h': /* turn off header */
      doheader = 0;
      break;
    case '5': /* 5-phase calculation */
      nsuns = 1;
      fixed_sun_sa = PI / 360. * atof(argv[++i]);
      if (fixed_sun_sa <= 0) {
        fprintf(stderr,
                "%s: missing solar disk size argument for '-5' option\n",
                progname);
        exit(1);
      }
      fixed_sun_sa *= fixed_sun_sa * PI;
      break;
    case 'o': /* output format */
      switch (argv[i][2]) {
      case 'f':
      case 'd':
      case 'a':
        outfmt = argv[i][2];
        break;
      default:
        goto userr;
      }
      break;
    default:
      goto userr;
    }
  }
  if (i < argc - 1)
    goto userr;
  if (i == argc - 1 && freopen(argv[i], "r", stdin) == NULL) {
    fprintf(stderr, "%s: cannot open '%s' for input\n", progname, argv[i]);
    exit(1);
  }
  if (verbose) {
    if (i == argc - 1)
      fprintf(stderr, "%s: reading weather tape '%s'\n", progname, argv[i]);
    else
      fprintf(stderr, "%s: reading weather tape from <stdin>\n", progname);
  }
  /* read weather tape header */
  if (scanf("place %[^\r\n] ", buf) != 1)
    goto fmterr;
  if (scanf("latitude %lf\n", &s_latitude) != 1)
    goto fmterr;
  if (scanf("longitude %lf\n", &s_longitude) != 1)
    goto fmterr;
  if (scanf("time_zone %lf\n", &s_meridian) != 1)
    goto fmterr;
  if (scanf("site_elevation %lf\n", &elevation) != 1)
    goto fmterr;
  if (scanf("weather_data_file_units %d\n", &input) != 1)
    goto fmterr;

  rh_init();
  if (verbose) {
    fprintf(stderr, "%s: location '%s'\n", progname, buf);
    fprintf(stderr, "%s: (lat,long)=(%.1f,%.1f) degrees north, west\n",
            progname, s_latitude, s_longitude);
    if (rotation != 0)
      fprintf(stderr, "%s: rotating output %.0f degrees\n", progname, rotation);
  }

  s_latitude = DegToRad(s_latitude);
  s_longitude = DegToRad(s_longitude);
  s_meridian = DegToRad(s_meridian);
  /* initial allocation */
  mtx_data = resize_dmatrix(mtx_data, tstorage = 2, nskypatch);

  /* Load mie density data */
  DATARRAY *mie_dp = getdata(mie_path);
  if (mie_dp == NULL) {
    fprintf(stderr, "Error reading mie data\n");
    return 0;
  }

  while (scanf("%d %d %lf %lf %lf %lf %lf\n", &mo, &da, &hr, &dni, &dhi, &aod,
               &cc) == 7) {
    if (aod == 0.0) {
      aod = AOD0_CA;
      fprintf(stderr, "aod is zero, using default value %.3f\n", AOD0_CA);
    }
    double sda, sta;
    int sun_in_sky;
    /* compute solar position */
    if ((mo == 2) & (da == 29)) {
      julian_date = 60;
      leap_day = 1;
    } else
      julian_date = jdate(mo, da) + leap_day;
    sda = sdec(julian_date);
    sta = stadj(julian_date);
    altitude = salt(sda, hr + sta);
    sun_in_sky = (altitude > -DegToRad(SUN_ANG_DEG / 2.));

    azimuth = sazi(sda, hr + sta) + PI - DegToRad(rotation);

    vectorize(altitude, azimuth, sundir);
    if (sun_hours_only && sundir[2] <= 0.) {
      continue; /* skipping nighttime points */
    }
    sun_ct = fdot(view_point, sundir) / ER;

    mtx_offset = NSSAMP * nskypatch * nstored;
    nstored += 1;
    printf("mtx_offset = %d nstored = %d nskypatch = %d\n", mtx_offset, nstored,
           nskypatch);
    /* make space for next row */
    if (nstored > tstorage) {
      printf("make space for next row nstored = %d tstorage = %d\n", nstored,
             tstorage);
      tstorage += (tstorage >> 1) + nstored + 7;
      mtx_data = resize_dmatrix(mtx_data, tstorage, nskypatch);
    }
    ntsteps++; /* keep count of time steps */
               /* compute sky patch values */
    Atmosphere clear_atmos = init_atmos(aod, grefl);
    int is_summer = (mo >= SUMMER_START && mo <= SUMMER_END);
    if (s_latitude < 0) {
      is_summer = !is_summer;
    }
    set_rayleigh_density_profile(&clear_atmos, lstag, is_summer, s_latitude);

    clear_atmos.beta_m = mie_dp;

    char gsdir[PATH_MAX];
    size_t siz = strlen(ddir);
    if (ISDIRSEP(ddir[siz - 1]))
      ddir[siz - 1] = '\0';
    snprintf(gsdir, PATH_MAX, "%s%catmos_data", ddir, DIRSEP);
    if (!make_directory(gsdir)) {
      fprintf(stderr, "Failed creating atmos_data directory");
      exit(1);
    }
    DpPaths clear_paths = get_dppaths(gsdir, aod, mie_name, lstag);

    if (getpath(clear_paths.tau, ".", R_OK) == NULL ||
        getpath(clear_paths.scat, ".", R_OK) == NULL ||
        getpath(clear_paths.scat1m, ".", R_OK) == NULL ||
        getpath(clear_paths.irrad, ".", R_OK) == NULL) {
      printf("# Pre-computing...\n");
      if (!precompute(sorder, clear_paths, &clear_atmos, num_threads)) {
        fprintf(stderr, "Pre-compute failed\n");
        return 0;
      }
    }

    DATARRAY *tau_clear_dp = getdata(clear_paths.tau);
    DATARRAY *irrad_clear_dp = getdata(clear_paths.irrad);
    DATARRAY *scat_clear_dp = getdata(clear_paths.scat);
    DATARRAY *scat1m_clear_dp = getdata(clear_paths.scat1m);

    if (!solar_only)
      compute_sky(tau_clear_dp, scat_clear_dp, scat1m_clear_dp, irrad_clear_dp,
                  cc, mtx_data + mtx_offset);
    if (!sky_only)
      add_direct(tau_clear_dp, scat_clear_dp, scat1m_clear_dp, irrad_clear_dp,
                 cc, mtx_data + mtx_offset);
    /* monthly reporting */
    if (verbose && mo != last_monthly)
      fprintf(stderr, "%s: stepping through month %d...\n", progname,
              last_monthly = mo);
  }
  freedata(mie_dp);
  if (!ntsteps) {
    fprintf(stderr, "%s: no valid time steps on input\n", progname);
    exit(1);
  }
  /* check for junk at end */
  while ((i = fgetc(stdin)) != EOF)
    if (!isspace(i)) {
      fprintf(stderr, "%s: warning - unexpected data past EOT: ", progname);
      buf[0] = i;
      buf[1] = '\0';
      fgets(buf + 1, sizeof(buf) - 1, stdin);
      fputs(buf, stderr);
      fputc('\n', stderr);
      break;
    }
  /* write out matrix */
  if (outfmt != 'a')
    SET_FILE_BINARY(stdout);
#ifdef getc_unlocked
  flockfile(stdout);
#endif
  if (verbose)
    fprintf(stderr, "%s: writing %smatrix with %d time steps...\n", progname,
            outfmt == 'a' ? "" : "binary ", nstored);
  if (doheader) {
    newheader("RADIANCE", stdout);
    printargs(argc, argv, stdout);
    printf("LATLONG= %.8f %.8f\n", RadToDeg(s_latitude),
           -RadToDeg(s_longitude));
    printf("NROWS=%d\n", nskypatch);
    printf("NCOLS=%d\n", nstored);
    printf("NCOMP=%d\n", NSSAMP);
    if ((outfmt == 'f') | (outfmt == 'd'))
      fputendian(stdout);
    fputformat((char *)getfmtname(outfmt), stdout);
    putchar('\n');
  }
  /* patches are rows (outer sort) */
  for (i = 0; i < nskypatch; i++) {
    mtx_offset = NSSAMP * i;
    switch (outfmt) {
    case 'a':
      for (j = 0; j < nstored; j++) {
	int	k;
        for (k = 0; k < NSSAMP; k++) {
          printf("%.3g ", mtx_data[mtx_offset + k]);
        }
        printf("\n");
        mtx_offset += NSSAMP * nskypatch;
      }
      if (nstored > 1)
        fputc('\n', stdout);
      break;
    case 'f':
      for (j = 0; j < nstored; j++) {
        putbinary(mtx_data + mtx_offset, sizeof(float), NSSAMP, stdout);
        mtx_offset += NSSAMP * nskypatch;
      }
      break;
    case 'd':
      for (j = 0; j < nstored; j++) {
        double ment[NSSAMP];
        for (j = 0; j < NSSAMP; j++)
          ment[j] = mtx_data[mtx_offset + j];
        putbinary(ment, sizeof(double), NSSAMP, stdout);
        mtx_offset += NSSAMP * nskypatch;
      }
      break;
    }
    if (ferror(stdout))
      goto writerr;
  }
alldone:
  if (fflush(NULL) == EOF)
    goto writerr;
  if (verbose)
    fprintf(stderr, "%s: done.\n", progname);
  exit(0);
userr:
  fprintf(stderr,
          "Usage: %s [-v][-h][-A][-d|-s|-n][-u][-D file [-M modfile]][-r "
          "deg][-m N][-g r g b][-c r g b][-o{f|d}][-O{0|1}] [tape.wea]\n",
          progname);
  exit(1);
fmterr:
  fprintf(stderr, "%s: weather tape format error in header\n", progname);
  exit(1);
writerr:
  fprintf(stderr, "%s: write error on output\n", progname);
  exit(1);
}
