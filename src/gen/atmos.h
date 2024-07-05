#ifndef ATMOS_H
#define ATMOS_H

#include "color.h"
#include "data.h"
#include "fvect.h"
#include "paths.h"
#include "rtio.h"
#include "rtmath.h"
#include "sun.h"

#define NSSAMP 20

typedef struct {
  double width;
  double exp_term;
  double exp_scale;
  double linear_term;
  double constant_term;
} DensityProfileLayer;

typedef struct {
  DensityProfileLayer layers[2];
} DensityProfile;

typedef struct {
  DensityProfile rayleigh_density;
  DensityProfile ozone_density;
  const float *beta_r0;
  float beta_scale;
  DATARRAY *beta_m;
  const double grefl;
} Atmosphere;

typedef struct {
  char tau[PATH_MAX];
  char scat[PATH_MAX];
  char scat1m[PATH_MAX];
  char irrad[PATH_MAX];
} DpPaths;

extern const double ER;
extern const double AH;
extern const double HR_MS;
extern const double HR_MW;
extern const double HR_SS;
extern const double HR_SW;
extern const double HR_T;
extern const int WVLSPAN;
extern const float EXTSOL[NSSAMP];
extern const float BR0_MS[NSSAMP];
extern const float BR0_MW[NSSAMP];
extern const float BR0_SS[NSSAMP];
extern const float BR0_SW[NSSAMP];
extern const float BR0_T[NSSAMP];
extern const float BCLOUD;
extern const double AOD0_CA;
extern const double SOLOMG;

extern void get_rmumusnu(FVECT vpt, FVECT vdir, FVECT sundir, double *r,
                         double *mu, double *mu_s, double *nu);

extern void get_transmittance_to_sun(DATARRAY *tau_dp, const double r,
                                     const double mu_s, double *result);

extern void get_sky_transmittance(DATARRAY *tau, double r, double mu, float *result);

extern void get_sky_radiance(DATARRAY *scat, DATARRAY *scat1m, const double radius, 
                             const double mu, const double mu_s, const double nu, float *result);

extern void get_solar_radiance(DATARRAY *tau, DATARRAY *scat, DATARRAY *scat1m, const FVECT sundir, const double radius, const double sun_ct, double *sun_radiance); 

extern void get_ground_radiance(DATARRAY *tau, DATARRAY *scat, DATARRAY *scat1m, DATARRAY *irrad, 
                                const FVECT view_point, const FVECT view_direction, const double radius, const double mu, const double sun_ct, const double nu, 
                                const double grefl, const FVECT sundir, float *result);

extern void add_cloud_radiance(DATARRAY *scat, double nu, double pt[4],
                               float *result);

extern int compute_sundir(const int year, const int month, const int day,
                          const double hour, const int tsolar,
                          double sundir[3]);

extern int precompute(const int sorder, const DpPaths dppaths, const Atmosphere *atmos,
                      int num_threads);

#endif // ATMOS_H
