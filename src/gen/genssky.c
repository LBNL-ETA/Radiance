// Main function for generating spectral sky
// Cloudy sky computed as weight average of clear and cie overcast sky

#include "copyright.h"
#include "atmos.h"
#include "resolu.h"
#include "view.h"


char *progname;

const double ARCTIC_LAT = 67.;
const double TROPIC_LAT = 23.;
const int SUMMER_START = 4;
const int SUMMER_END = 9;
const double GNORM = 0.777778;

const double D65EFF = 203.; /* standard illuminant D65 */

// Mean normalized relative daylight spectra where CCT = 6415K for overcast;
const double D6415[NSSAMP] = {0.63231, 1.06171, 1.00779, 1.36423, 1.34133,
                              1.27258, 1.26276, 1.26352, 1.22201, 1.13246,
                              1.0434,  1.05547, 0.98212, 0.94445, 0.9722,
                              0.82387, 0.87853, 0.82559, 0.75111, 0.78925};

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

// from gensky.c
static double get_overcast_brightness(const double dz, const double zenithbr) {
  double groundbr = zenithbr * GNORM;
  return wmean(pow(dz + 1.01, 10), zenithbr * (1 + 2 * dz) / 3,
               pow(dz + 1.01, -10), groundbr);
}

static void write_rad_file(FILE *fp, const double *sun_radiance,
                           const FVECT sundir, const char skyfile[PATH_MAX],
                           const char grndfile[PATH_MAX]) {
  if (sundir[2] > 0) {
    fprintf(fp, "void spectrum sunrad\n0\n0\n22 380 780 ");
    for (int i = 0; i < NSSAMP; ++i) {
      fprintf(fp, "%.1f ", sun_radiance[i] * WVLSPAN);
    }
    fprintf(fp, "\n\nsunrad light solar\n0\n0\n3 1 1 1\n\n");
    fprintf(fp, "solar source sun\n0\n0\n4 %f %f %f 0.533\n\n", sundir[0],
            sundir[1], sundir[2]);
  }
  fprintf(fp,
          "void specpict skyfunc\n8 noop %s fisheye.cal fish_u fish_v -rx 90 "
          "-mx\n0\n0\n\n",
          skyfile);
  fprintf(fp, "skyfunc glow sky_glow\n0\n0\n4 1 1 1 0\n\n");
  fprintf(fp, "sky_glow source sky\n0\n0\n4 0 0 1 180\n\n");

  fprintf(fp,
          "void specpict grndmap\n8 noop %s fisheye.cal fish_u fish_v -rx -90 "
          "-my\n0\n0\n\n",
          grndfile);
  fprintf(fp, "grndmap glow ground_glow\n0\n0\n4 1 1 1 0\n\n");
  fprintf(fp, "ground_glow source ground_source\n0\n0\n4 0 0 -1 180\n\n");
}

static void write_hsr_header(FILE *fp, RESOLU *res) {
  float wvsplit[4] = {380, 480, 588,
                      780}; // RGB wavelength limits+partitions (nm)
  newheader("RADIANCE", fp);
  fputncomp(NSSAMP, fp);
  fputwlsplit(wvsplit, fp);
  fputformat(SPECFMT, fp);
  fputc('\n', fp);
  fputsresolu(res, fp);
}

int gen_spect_sky(DATARRAY *tau_clear, DATARRAY *scat_clear,
                  DATARRAY *scat1m_clear, DATARRAY *irrad_clear,
                  const double cloud_cover, const FVECT sundir,
                  const double grefl, const int res, const char *outname) {

  char radfile[PATH_MAX];
  char skyfile[PATH_MAX];
  char grndfile[PATH_MAX];
  if (!snprintf(radfile, sizeof(radfile), "%s.rad", outname)) {
    fprintf(stderr, "Error setting rad file name\n");
    return 0;
  };
  if (!snprintf(skyfile, sizeof(skyfile), "%s_sky.hsr", outname)) {
    fprintf(stderr, "Error setting sky file name\n");
    return 0;
  };
  if (!snprintf(grndfile, sizeof(grndfile), "%s_ground.hsr", outname)) {
    fprintf(stderr, "Error setting ground file name\n");
    return 0;
  }
  RESOLU rs = {PIXSTANDARD, res, res};
  FILE *skyfp = fopen(skyfile, "w");
  FILE *grndfp = fopen(grndfile, "w");
  write_hsr_header(grndfp, &rs);
  write_hsr_header(skyfp, &rs);
  VIEW skyview = {VT_ANG, {0., 0., 0.}, {0., 0., 1.}, {0., 1., 0.}, 1.,
                  180.,   180.,         0.,           0.,           0.,
                  0.,     {0., 0., 0.}, {0., 0., 0.}, 0.,           0.};
  VIEW grndview = {
      VT_ANG, {0., 0., 0.}, {0., 0., -1.}, {0., 1., 0.}, 1., 180., 180., 0., 0.,
      0.,     0.,           {0., 0., 0.},  {0., 0., 0.}, 0., 0.};
  setview(&skyview);
  setview(&grndview);

  CNDX[3] = NSSAMP;

  FVECT view_point = {0, 0, ER};
  const double radius = VLEN(view_point);
  const double sun_ct = fdot(view_point, sundir) / radius;
  for (unsigned int j = 0; j < res; ++j) {
    for (unsigned int i = 0; i < res; ++i) {
      RREAL loc[2];
      FVECT rorg = {0};
      FVECT rdir_sky = {0};
      FVECT rdir_grnd = {0};
      SCOLOR sky_radiance = {0};
      SCOLOR ground_radiance = {0};
      SCOLR sky_sclr = {0};
      SCOLR ground_sclr = {0};

      pix2loc(loc, &rs, i, j);
      viewray(rorg, rdir_sky, &skyview, loc[0], loc[1]);
      viewray(rorg, rdir_grnd, &grndview, loc[0], loc[1]);

      const double mu_sky = fdot(view_point, rdir_sky) / radius;
      const double nu_sky = fdot(rdir_sky, sundir);

      const double mu_grnd = fdot(view_point, rdir_grnd) / radius;
      const double nu_grnd = fdot(rdir_grnd, sundir);

      get_sky_radiance(scat_clear, scat1m_clear, radius, mu_sky, sun_ct, nu_sky,
                       sky_radiance);
      get_ground_radiance(tau_clear, scat_clear, scat1m_clear, irrad_clear,
                          view_point, rdir_grnd, radius, mu_grnd, sun_ct,
                          nu_grnd, grefl, sundir, ground_radiance);

      for (int k = 0; k < NSSAMP; ++k) {
        sky_radiance[k] *= WVLSPAN;
        ground_radiance[k] *= WVLSPAN;
      }

      if (cloud_cover > 0) {
        double zenithbr = get_zenith_brightness(sundir);
        double grndbr = zenithbr * GNORM;
        double skybr = get_overcast_brightness(rdir_sky[2], zenithbr);
        for (int k = 0; k < NSSAMP; ++k) {
          sky_radiance[k] =
              wmean2(sky_radiance[k], skybr * D6415[k], cloud_cover);
          ground_radiance[k] =
              wmean2(ground_radiance[k], grndbr * D6415[k], cloud_cover);
        }
      }

      scolor2scolr(sky_sclr, sky_radiance, 20);
      putbinary(sky_sclr, LSCOLR, 1, skyfp);

      scolor2scolr(ground_sclr, ground_radiance, 20);
      putbinary(ground_sclr, LSCOLR, 1, grndfp);
    }
  }
  fclose(skyfp);
  fclose(grndfp);

  // Get solar radiance
  double sun_radiance[NSSAMP] = {0};
  get_solar_radiance(tau_clear, scat_clear, scat1m_clear, sundir, radius,
                     sun_ct, sun_radiance);
  if (cloud_cover > 0) {
    double zenithbr = get_zenith_brightness(sundir);
    double skybr = get_overcast_brightness(sundir[2], zenithbr);
    for (int i = 0; i < NSSAMP; ++i) {
      sun_radiance[i] =
          wmean2(sun_radiance[i], D6415[i] * skybr / WVLSPAN, cloud_cover);
    }
  }

  FILE *rfp = fopen(radfile, "w");
  write_rad_file(rfp, sun_radiance, sundir, skyfile, grndfile);
  fclose(rfp);
  return 1;
}

static DpPaths get_dppaths(const double aod, const char *tag) {
  DpPaths paths;

  snprintf(paths.tau, PATH_MAX, "tau_%s_%.2f.dat", tag, aod);
  snprintf(paths.scat, PATH_MAX, "scat_%s_%.2f.dat", tag, aod);
  snprintf(paths.scat1m, PATH_MAX, "scat1m_%s_%.2f.dat", tag, aod);
  snprintf(paths.irrad, PATH_MAX, "irrad_%s_%.2f.dat", tag, aod);

  return paths;
}

static void set_rayleigh_density_profile(Atmosphere *atmos, char *tag, const int is_summer,
                                         const double s_latitude) {
  // Set rayleigh density profile
  if (fabs(s_latitude*180.0 / PI) > ARCTIC_LAT) {
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
  } else if (fabs(s_latitude*180.0/PI) > TROPIC_LAT) {
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
  Atmosphere atmos = {
      .ozone_density = {.layers =
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
      .grefl = grefl
  };
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
  double grefl = 0.2;
  double ccover = 0.0;
  int res = 128;
  double aod = AOD0_CA;
  char *outname = "out";
  char *mie_path = getpath("mie_ca.dat", getrlibpath(), R_OK);
  char lstag[3];

  if (argc < 4) {
    fprintf(stderr, "Usage: %s month day hour -y year -a lat -o lon -m tz -d aod -r res -n nproc -c ccover -l mie -g grefl -f outpath\n",
            argv[0]);
    return 0;
  }

  month = atoi(argv[1]);
  day = atoi(argv[2]);
  hour = atof(argv[3]);

  if (!compute_sundir(year, month, day, hour, tsolar, sundir)) {
    fprintf(stderr, "Cannot compute solar angle\n");
    exit(1);
  }

  for (int i = 4; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'a':
        s_latitude = atof(argv[++i]) * (PI / 180.0);
        break;
      case 'g':
        grefl = atof(argv[++i]);
        break;
      case 'c':
        ccover = atof(argv[++i]);
        break;
      case 'd':
        aod = atof(argv[++i]);
        break;
      case 'i':
        sorder = atoi(argv[++i]);
        break;
      case 'l':
        mie_path = argv[++i];
        break;
      case 'm':
        s_meridian = atof(argv[++i]) * (PI / 180.0);
        break;
      case 'o':
        s_longitude = atof(argv[++i]) * (PI / 180.0);
        break;
      case 'n':
        num_threads = atoi(argv[++i]);
        break;
      case 'y':
        year = atoi(argv[++i]);
        break;
      case 'f':
        outname = argv[++i];
        break;
      case 'r':
        res = atoi(argv[++i]);
        break;
      default:
        fprintf(stderr, "Unknown option %s\n", argv[i]);
        exit(1);
      }
    }
  }

  Atmosphere clear_atmos = init_atmos(aod, grefl);

  int is_summer = (month >= SUMMER_START && month <= SUMMER_END);
  if (s_latitude < 0) {
    is_summer = !is_summer;
  }
  set_rayleigh_density_profile(&clear_atmos, lstag, is_summer, s_latitude);

  // Load mie density data
  DATARRAY *mie_dp = getdata(mie_path);
  if (mie_dp == NULL) {
    fprintf(stderr, "Error reading mie data\n");
    return 0;
  }
  clear_atmos.beta_m = mie_dp;

  DpPaths clear_paths = get_dppaths(aod, lstag);

  if (getpath(clear_paths.tau, ".", R_OK) == NULL ||
      getpath(clear_paths.scat, ".", R_OK) == NULL ||
      getpath(clear_paths.scat1m, ".", R_OK) == NULL ||
      getpath(clear_paths.irrad, ".", R_OK) == NULL) {
    printf("# Precomputing...\n");
    if (!precompute(sorder, clear_paths, &clear_atmos, num_threads)) {
      fprintf(stderr, "Precompute failed\n");
      return 0;
    }
  }

  DATARRAY *tau_clear_dp = getdata(clear_paths.tau);
  DATARRAY *irrad_clear_dp = getdata(clear_paths.irrad);
  DATARRAY *scat_clear_dp = getdata(clear_paths.scat);
  DATARRAY *scat1m_clear_dp = getdata(clear_paths.scat1m);

  if (!gen_spect_sky(tau_clear_dp, scat_clear_dp, scat1m_clear_dp,
                     irrad_clear_dp, ccover, sundir, grefl, res, outname)) {
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
