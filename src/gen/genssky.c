#include "color.h"
#ifndef lint
static const char RCSid[] =
    "$Id: genssky.c,v 2.7 2025/04/10 23:30:58 greg Exp $";
#endif
/* Main function for generating spectral sky */
/* Cloudy sky computed as weight average of clear and cie overcast sky */

#include "atmos.h"
#include "copyright.h"
#include "resolu.h"
#include "rtio.h"
#include <ctype.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

char *progname;

const double ARCTIC_LAT = 67.;
const double TROPIC_LAT = 23.;
const int SUMMER_START = 4;
const int SUMMER_END = 9;
const double GNORM = 0.777778;

const double D65EFF = 203.; /* standard illuminant D65 */

/* Mean normalized relative daylight spectra where CCT = 6415K for overcast; */
const double D6415[NSSAMP] = {0.63231, 1.06171, 1.00779, 1.36423, 1.34133,
                              1.27258, 1.26276, 1.26352, 1.22201, 1.13246,
                              1.0434,  1.05547, 0.98212, 0.94445, 0.9722,
                              0.82387, 0.87853, 0.82559, 0.75111, 0.78925};

/* European and North American zones */
struct {
  char zname[8]; /* time zone name (all caps) */
  float zmer;    /* standard meridian */
} tzone[] = {{"YST", 135},   {"YDT", 120},   {"PST", 120},  {"PDT", 105},
             {"MST", 105},   {"MDT", 90},    {"CST", 90},   {"CDT", 75},
             {"EST", 75},    {"EDT", 60},    {"AST", 60},   {"ADT", 45},
             {"NST", 52.5},  {"NDT", 37.5},  {"GMT", 0},    {"BST", -15},
             {"CET", -15},   {"CEST", -30},  {"EET", -30},  {"EEST", -45},
             {"AST", -45},   {"ADT", -60},   {"GST", -60},  {"GDT", -75},
             {"IST", -82.5}, {"IDT", -97.5}, {"JST", -135}, {"NDT", -150},
             {"NZST", -180}, {"NZDT", -195}, {"", 0}};

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

inline static float deg2rad(float deg) { return deg * (PI / 180.); }

static int cvthour(char *hs, int *tsolar, double *hour) {
  char *cp = hs;
  int i, j;

  if ((*tsolar = *cp == '+'))
    cp++; /* solar time? */
  while (isdigit(*cp))
    cp++;
  if (*cp == ':')
    *hour = atoi(hs) + atoi(++cp) / 60.0;
  else {
    *hour = atof(hs);
    if (*cp == '.')
      cp++;
  }
  while (isdigit(*cp))
    cp++;
  if (!*cp)
    return (0);
  if (*tsolar || !isalpha(*cp)) {
    fprintf(stderr, "%s: bad time format: %s\n", progname, hs);
    exit(1);
  }
  i = 0;
  do {
    for (j = 0; cp[j]; j++)
      if (toupper(cp[j]) != tzone[i].zname[j])
        break;
    if (!cp[j] && !tzone[i].zname[j]) {
      s_meridian = tzone[i].zmer * (PI / 180);
      return (1);
    }
  } while (tzone[i++].zname[0]);

  fprintf(stderr, "%s: unknown time zone: %s\n", progname, cp);
  fprintf(stderr, "Known time zones:\n\t%s", tzone[0].zname);
  for (i = 1; tzone[i].zname[0]; i++)
    fprintf(stderr, " %s", tzone[i].zname);
  putc('\n', stderr);
  exit(1);
}

static void basename(const char *path, char *output, size_t outsize) {
  const char *last_slash = strrchr(path, '/');
  const char *last_backslash = strrchr(path, '\\');
  const char *filename = path;
  const char *last_dot;

  if (last_slash && last_backslash) {
    filename =
        (last_slash > last_backslash) ? last_slash + 1 : last_backslash + 1;
  } else if (last_slash) {
    filename = last_slash + 1;
  } else if (last_backslash) {
    filename = last_backslash + 1;
  }

  last_dot = strrchr(filename, '.');
  if (last_dot) {
    size_t length = last_dot - filename;
    if (length < outsize) {
      strncpy(output, filename, length);
      output[length] = '\0';
    } else {
      strncpy(output, filename, outsize - 1);
      output[outsize - 1] = '\0';
    }
  }
}

static char *join_paths(const char *path1, const char *path2) {
  size_t len1 = strlen(path1);
  size_t len2 = strlen(path2);
  int need_separator = (path1[len1 - 1] != DIRSEP);

  char *result = malloc(len1 + len2 + (need_separator ? 2 : 1));
  if (!result)
    return NULL;

  strcpy(result, path1);
  if (need_separator) {
    result[len1] = DIRSEP;
    len1++;
  }
  strcpy(result + len1, path2);

  return result;
}

static inline double wmean2(const double a, const double b, const double x) {
  return a * (1 - x) + b * x;
}

static inline double wmean(const double a, const double x, const double b,
                           const double y) {
  return (a * x + b * y) / (a + b);
}

static double get_overcast_zenith_brightness(const double sundir[3]) {
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

static void write_header(const int argc, char **argv, const double cloud_cover,
                         const double grefl, const int res) {
  int i;
  printf("# ");
  for (i = 0; i < argc; i++) {
    printf("%s ", argv[i]);
  }
  printf("\n");
  printf(
      "#Cloud cover: %g\n#Ground reflectance: %g\n#Sky map resolution: %d\n\n",
      cloud_cover, grefl, res);
}

static void write_rad(const double *sun_radiance, const double intensity,
                      const FVECT sundir, const char *ddir,
                      const char *skyfile) {
  if (sundir[2] > 0) {
    printf("void spectrum sunrad\n0\n0\n22 380 780 ");
    int i;
    for (i = 0; i < NSSAMP; ++i) {
      printf("%.3f ", sun_radiance[i]);
    }
    printf("\n\nsunrad light solar\n0\n0\n3 %.1f %.1f %.1f\n\n", intensity,
           intensity, intensity);
    printf("solar source sun\n0\n0\n4 %f %f %f 0.533\n\n", sundir[0], sundir[1],
           sundir[2]);
  }
  printf("void specpict skyfunc\n5 noop %s . 'Atan2(Dy,Dx)/PI+1' "
         "'1-Acos(Dz)/PI'\n0\n0\n\n",
         skyfile);
}

static void write_hsr_header(FILE *fp, RESOLU *res) {
  float wvsplit[4] = {380, 480, 588, 780};
  newheader("RADIANCE", fp);
  fputncomp(NSSAMP, fp);
  fputwlsplit(wvsplit, fp);
  fputformat(SPECFMT, fp);
  fputc('\n', fp);
  fputsresolu(res, fp);
}

static inline float frac(float x) { return x - floor(x); }

int gen_spect_sky(DATARRAY *tau_clear, DATARRAY *scat_clear,
                  DATARRAY *scat1m_clear, DATARRAY *irrad_clear,
                  const double cloud_cover, const FVECT sundir,
                  const double grefl, const int res, const char *outname,
                  const char *ddir, const double dirnorm, const double difhor) {
  char skyfile[PATH_MAX];
  if (!snprintf(skyfile, sizeof(skyfile), "%s%c%s_sky.hsr", ddir, DIRSEP,
                outname)) {
    fprintf(stderr, "Error setting sky file name\n");
    return 0;
  };
  int xres = res;
  int yres = xres / 2;
  RESOLU rs = {PIXSTANDARD, xres, yres};
  FILE *skyfp = fopen(skyfile, "w");
  write_hsr_header(skyfp, &rs);

  CNDX[3] = NSSAMP;

  FVECT view_point = {0, 0, ER + 10};
  const double radius = VLEN(view_point);
  const double sun_ct = fdot(view_point, sundir) / radius;

  double overcast_zenithbr = get_overcast_zenith_brightness(sundir);
  double overcast_grndbr = overcast_zenithbr * GNORM;

  double dif_ratio = 1;
  if (difhor > 0) {
    DATARRAY *indirect_irradiance_clear = get_indirect_irradiance(irrad_clear, radius, sun_ct);
    double overcast_ghi = overcast_zenithbr * 7.0 * PI / 9.0;
    double diffuse_irradiance = 0;
    int l;
    for (l = 0; l < NSSAMP; ++l) {
      diffuse_irradiance += indirect_irradiance_clear->arr.d[l] * 20;  /* 20nm interval */
    }
    free(indirect_irradiance_clear);
    diffuse_irradiance = wmean2(diffuse_irradiance, overcast_ghi, cloud_cover);
    if (diffuse_irradiance > 0) {
        dif_ratio = difhor / WHTEFFICACY / diffuse_irradiance / 1.15;       /* fudge */
    }
  }
  int i, j, k;
  for (j = 0; j < yres; ++j) {
    for (i = 0; i < xres; ++i) {
      SCOLOR radiance = {0};
      SCOLR sky_sclr = {0};

      float px = i / (xres - 1.0);
      float py = j / (yres - 1.0);
      float lambda = ((1 - py) * PI) - (PI / 2.0);
      float phi = (px * 2.0 * PI) - PI;

      FVECT rdir = {cos(lambda) * cos(phi), cos(lambda) * sin(phi),
                    sin(lambda)};

      const double mu = fdot(view_point, rdir) / radius;
      const double nu = fdot(rdir, sundir);

      /* hit ground */
      if (rdir[2] < 0) {
        get_ground_radiance(tau_clear, scat_clear, scat1m_clear, irrad_clear,
                            view_point, rdir, radius, mu, sun_ct, nu, grefl,
                            sundir, radiance);
      } else {
        get_sky_radiance(scat_clear, scat1m_clear, radius, mu, sun_ct, nu,
                         radiance);
      }

      for (k = 0; k < NSSAMP; ++k) {
        radiance[k] *= WVLSPAN;
      }

      if (cloud_cover > 0) {
        double skybr = get_overcast_brightness(rdir[2], overcast_zenithbr);
        if (rdir[2] < 0) {
          for (k = 0; k < NSSAMP; ++k) {
            radiance[k] = wmean2(radiance[k], overcast_grndbr * D6415[k], cloud_cover);
          }
        } else {
          for (k = 0; k < NSSAMP; ++k) {
            radiance[k] = wmean2(radiance[k], skybr * D6415[k], cloud_cover);
          }
        }
      }

      for (k = 0; k < NSSAMP; ++k) {
        radiance[k] *= dif_ratio;
      }

      scolor2scolr(sky_sclr, radiance, NSSAMP);
      putbinary(sky_sclr, LSCOLR, 1, skyfp);
    }
  }
  fclose(skyfp);

  /* Get solar radiance */
  double sun_radiance[NSSAMP] = {0};
  get_solar_radiance(tau_clear, scat_clear, scat1m_clear, sundir, radius,
                     sun_ct, sun_radiance);
  if (cloud_cover > 0) {
    double skybr = get_overcast_brightness(sundir[2], overcast_zenithbr);
    int i;
    for (i = 0; i < NSSAMP; ++i) {
      sun_radiance[i] =
          wmean2(sun_radiance[i], D6415[i] * skybr / WVLSPAN, cloud_cover);
    }
  }

  /* Normalize */
  double sum = 0.0;
  for (i = 0; i < NSSAMP; ++i) {
    sum += sun_radiance[i];
  }
  double mean = sum / NSSAMP;
  for (i = 0; i < NSSAMP; ++i) {
    sun_radiance[i] /= mean;
  }
  double intensity = mean * WVLSPAN;
  if (dirnorm > 0) {
    intensity = dirnorm / SOLOMG / WHTEFFICACY;
  }

  write_rad(sun_radiance, intensity, sundir, ddir, skyfile);
  return 1;
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

int main(int argc, char *argv[]) {
  progname = argv[0];
  int month, day;
  double hour;
  FVECT sundir;
  int num_threads = 1;
  int sorder = 4;
  int year = 0;
  int tsolar = 0;
  int got_meridian = 0;
  double grefl = 0.2;
  double ccover = 0.0;
  int res = 64;
  double aod = AOD0_CA;
  char *outname = "out";
  char *mie_path = getpath("mie_ca.dat", getrlibpath(), R_OK);
  char mie_name[20] = "mie_ca";
  char lstag[3];
  char *ddir = ".";
  int i;
  double dirnorm = 0; /* direct normal illuminance */
  double difhor = 0;  /* diffuse horizontal illuminance */

  if (argc == 2 && !strcmp(argv[1], "-defaults")) {
    printf("-i %d\t\t\t\t#scattering order\n", sorder);
    printf("-g %f\t\t\t#ground reflectance\n", grefl);
    printf("-c %f\t\t\t#cloud cover\n", ccover);
    printf("-r %d\t\t\t\t#image resolution\n", res);
    printf("-d %f\t\t\t#broadband aerosol optical depth\n", AOD0_CA);
    printf("-f %s\t\t\t\t#output name (-f)\n", outname);
    printf("-p %s\t\t\t\t#atmos data directory\n", ddir);
    exit(0);
  }

  if (argc < 4) {
    fprintf(stderr,
            "Usage: %s month day hour -y year -a lat -o lon -m tz -d aod -r "
            "res -n nproc -c ccover -l mie -L dirnorm_illum difhor_illum "
	    "-g grefl -f outpath\n",
            argv[0]);
    return 0;
  }

  month = atoi(argv[1]);
  if (month < 1 || month > 12) {
    fprintf(stderr, "bad month");
    exit(1);
  }
  day = atoi(argv[2]);
  if (day < 1 || day > 31) {
    fprintf(stderr, "bad month");
    exit(1);
  }
  got_meridian = cvthour(argv[3], &tsolar, &hour);

  if (!compute_sundir(year, month, day, hour, tsolar, sundir)) {
    fprintf(stderr, "Cannot compute solar angle\n");
    exit(1);
  }

  for (i = 4; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'a':
        s_latitude = atof(argv[++i]) * (PI / 180.0);
        break;
      case 'c':
        ccover = atof(argv[++i]);
        break;
      case 'd':
        aod = atof(argv[++i]);
        break;
      case 'f':
        outname = argv[++i];
        break;
      case 'g':
        grefl = atof(argv[++i]);
        break;
      case 'i':
        sorder = atoi(argv[++i]);
        break;
      case 'l':
        mie_path = argv[++i];
        basename(mie_path, mie_name, sizeof(mie_name));
        break;
      case 'm':
        if (got_meridian) {
          ++i;
          break;
        }
        s_meridian = atof(argv[++i]) * (PI / 180.0);
        break;
      case 'n':
        num_threads = atoi(argv[++i]);
        break;
      case 'o':
        s_longitude = atof(argv[++i]) * (PI / 180.0);
        break;
      case 'L':
        dirnorm = atof(argv[++i]);
        difhor = atof(argv[++i]);
        break;
      case 'p':
        ddir = argv[++i];
        break;
      case 'r':
        res = atoi(argv[++i]);
        break;
      case 'y':
        year = atoi(argv[++i]);
        break;
      default:
        fprintf(stderr, "Unknown option %s\n", argv[i]);
        exit(1);
      }
    }
  }
  if (year && (year < 1950) | (year > 2050))
    fprintf(stderr, "%s: warning - year should be in range 1950-2050\n",
            progname);
  if (month && !tsolar && fabs(s_meridian - s_longitude) > 45 * PI / 180)
    fprintf(stderr,
            "%s: warning - %.1f hours btwn. standard meridian and longitude\n",
            progname, (s_longitude - s_meridian) * 12 / PI);

  Atmosphere clear_atmos = init_atmos(aod, grefl);

  int is_summer = (month >= SUMMER_START && month <= SUMMER_END);
  if (s_latitude < 0) {
    is_summer = !is_summer;
  }
  set_rayleigh_density_profile(&clear_atmos, lstag, is_summer, s_latitude);

  /* Load mie density data */
  DATARRAY *mie_dp = getdata(mie_path);
  if (mie_dp == NULL) {
    fprintf(stderr, "Error reading mie data\n");
    return 0;
  }
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

  write_header(argc, argv, ccover, grefl, res);

  if (!gen_spect_sky(tau_clear_dp, scat_clear_dp, scat1m_clear_dp,
                     irrad_clear_dp, ccover, sundir, grefl, res, outname, ddir,
                     dirnorm, difhor)) {
    fprintf(stderr, "gen_spect_sky failed\n");
    exit(1);
  }

  freedata(mie_dp);
  freedata(tau_clear_dp);
  freedata(scat_clear_dp);
  freedata(irrad_clear_dp);
  freedata(scat1m_clear_dp);

  return 1;
}
