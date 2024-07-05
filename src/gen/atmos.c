#include "atmos.h"
#include "data.h"

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#else
#include <pthread.h>
#endif

#if defined(_WIN32) || defined(_WIN64)
#define THREAD HANDLE
#define CREATE_THREAD(thread, func, context)                                   \
  ((*(thread) = CreateThread(NULL, 0, (func),          \
                             (context), 0, NULL)) != NULL)
#define THREAD_RETURN DWORD WINAPI
#define THREAD_JOIN(thread) WaitForSingleObject(thread, INFINITE)
#define THREAD_CLOSE(thread) CloseHandle(thread)
#else
#define THREAD pthread_t
#define CREATE_THREAD(thread, func, context)                                   \
  (pthread_create((thread), NULL, (func), (context)) == 0)
#define THREAD_RETURN void *
#define THREAD_JOIN(thread) pthread_join(thread, NULL)
#endif

typedef struct {
  int start_l;
  int end_l;
  const Atmosphere *atmos;
  DATARRAY *tau_dp;
  DATARRAY *delta_rayleigh_scattering_dp;
  DATARRAY *delta_mie_scattering_dp;
  DATARRAY *scattering_dp;
} Scat1Tdat;

typedef struct {
  int start_l;
  int end_l;
  DATARRAY *tau_dp;
  DATARRAY *delta_scattering_density_dp;
  DATARRAY *delta_multiple_scattering_dp;
  DATARRAY *scattering_dp;
} ScatNTdat;

typedef struct {
  int start_l;
  int end_l;
  const Atmosphere *atmos;
  DATARRAY *tau_dp;
  DATARRAY *delta_rayleigh_scattering_dp;
  DATARRAY *delta_mie_scattering_dp;
  DATARRAY *delta_multiple_scattering_dp;
  DATARRAY *delta_irradiance_dp;
  int scattering_order;
  DATARRAY *delta_scattering_density_dp;
} ScatDenseTdat;

const double ER = 6360000.0; // Earth radius in meters
const double AR = 6420000; // Atmosphere radius in meters
const double AH = 60000;   // Atmosphere thickness in meters

const float START_WVL = 390.0;
const float END_WVL = 770.0;
const double SUNRAD = 0.004625;

// The cosine of the maximum Sun zenith angle
const double MU_S_MIN = -0.2;

// Scale heights (m)
// Rayleigh scattering
const double HR_MS = 8820; // Midlatitude summer
const double HR_MW = 8100; // Midlatitude winter
const double HR_SS = 8550; // Subarctic summer
const double HR_SW = 7600; // Subarctic winter
const double HR_T = 9050;  // Tropical

const double AOD0_CA = 0.115; // Broadband AOD for continental average

const double SOLOMG = 6.7967e-5; // .533 apex angle
const int WVLSPAN = 400;         // 380nm to 780nm

// Aerosol optical depth at 550nm for continental clean
const double AOD_CC_550 = 0.05352;

/** 380nm to 780nm at 20nm intervals **/
// Extraterrestrial solar W/m^2/nm
const float EXTSOL[NSSAMP] = {1.11099345, 1.75363199, 1.67126921, 1.97235667,
                              2.01863015, 1.93612482, 1.87310758, 1.88792588,
                              1.86274213, 1.84433248, 1.80370942, 1.73020456,
                              1.67036654, 1.57482149, 1.53646383, 1.45515319,
                              1.38207435, 1.3265141,  1.26648111, 1.20592043};

// Rayleight scattering coefficients at sea level in m^-1
const float BR0_MS[NSSAMP] = {
    4.63889967e-05, 3.76703293e-05, 3.09193474e-05, 2.56223081e-05,
    2.14198405e-05, 1.80483398e-05, 1.53177709e-05, 1.30867298e-05,
    1.12487269e-05, 9.72370830e-06, 8.44932642e-06, 7.37769845e-06,
    6.47116653e-06, 5.70005327e-06, 5.04091603e-06, 4.47430240e-06,
    3.98549839e-06, 3.56178512e-06, 3.19293761e-06, 2.87072599e-06};

const float BR0_MW[NSSAMP] = {
    5.03857564e-05, 4.09159104e-05, 3.35832809e-05, 2.78298619e-05,
    2.32653203e-05, 1.96033395e-05, 1.66375116e-05, 1.42142496e-05,
    1.22178890e-05, 1.05614786e-05, 9.17729921e-06, 8.01334246e-06,
    7.02870602e-06, 6.19115557e-06, 5.47522872e-06, 4.85979708e-06,
    4.32887894e-06, 3.86865960e-06, 3.46803311e-06, 3.11806054e-06};

const float BR0_SS[NSSAMP] = {
    4.73789047e-05, 3.84741872e-05, 3.15791442e-05, 2.61690698e-05,
    2.18769246e-05, 1.84334785e-05, 1.56446412e-05, 1.33659912e-05,
    1.14887667e-05, 9.93120528e-06, 8.62962901e-06, 7.53513326e-06,
    6.60925660e-06, 5.82168833e-06, 5.14848557e-06, 4.56978082e-06,
    4.07054608e-06, 3.63779107e-06, 3.26107261e-06, 2.93198523e-06};

const float BR0_SW[NSSAMP] = {
    5.30623659e-05, 4.30894595e-05, 3.53673035e-05, 2.93082494e-05,
    2.45012287e-05, 2.06447149e-05, 1.75213353e-05, 1.49693438e-05,
    1.28669320e-05, 1.11225291e-05, 9.66481886e-06, 8.43902999e-06,
    7.40208736e-06, 6.52004427e-06, 5.76608570e-06, 5.11796089e-06,
    4.55883914e-06, 4.07417187e-06, 3.65226316e-06, 3.28369923e-06};
const float BR0_T[NSSAMP] = {
    4.55376661e-05, 3.69790036e-05, 3.03519157e-05, 2.51520876e-05,
    2.10267438e-05, 1.77171168e-05, 1.50366593e-05, 1.28465622e-05,
    1.10422904e-05, 9.54525887e-06, 8.29426444e-06, 7.24230298e-06,
    6.35240772e-06, 5.59544593e-06, 4.94840517e-06, 4.39219003e-06,
    3.91235655e-06, 3.49641927e-06, 3.13434084e-06, 2.81804244e-06};

// Ozone extinction coefficients at sea level
const float AO0[NSSAMP] = {
    3.5661864e-09, 1.4848362e-08, 4.5415674e-08, 1.2446184e-07, 2.6461576e-07,
    4.8451984e-07, 8.609148e-07,  1.480537e-06,  1.8809e-06,    2.5107328e-06,
    2.5263174e-06, 2.313507e-06,  1.727741e-06,  1.2027012e-06, 7.915902e-07,
    5.0639202e-07, 3.5285684e-07, 2.23021e-07,   1.7395638e-07, 1.5052574e-07};

const double MIE_G = 0.7;  // Mie eccentricity

const int TRANSMITTANCE_TEXTURE_WIDTH = 256;
const int TRANSMITTANCE_TEXTURE_HEIGHT = 64;

// Resolution halved by original impl
const int SCATTERING_TEXTURE_R_SIZE = 16;
const int SCATTERING_TEXTURE_MU_SIZE = 64;
const int SCATTERING_TEXTURE_MU_S_SIZE = 16;
const int SCATTERING_TEXTURE_NU_SIZE = 4;

const int IRRADIANCE_TEXTURE_WIDTH = 64;
const int IRRADIANCE_TEXTURE_HEIGHT = 16;

static inline double clamp_cosine(double mu) {
  return fmax(-1.0, fmin(1.0, mu));
}

static inline double clamp_distance(double d) { return fmax(d, 0.0); }

static inline double clamp_radius(double r) { return fmax(ER, fmin(r, AR)); }

static inline double safe_sqrt(double a) { return sqrt(fmax(a, 0.0)); }

static inline double distance_to_space(const double radius, const double mu) {
  const double descriminant = radius * radius * (mu * mu - 1.0) + AR * AR;
  return clamp_distance(-radius * mu + safe_sqrt(descriminant));
}

static inline double distance_to_earth(const double radius, const double mu) {
  const double descriminant = radius * radius * (mu * mu - 1.0) + ER * ER;
  return clamp_distance(-radius * mu - safe_sqrt(descriminant));
}

static inline int ray_hits_ground(double radius, double mu) {
  return mu < 0.0 && radius * radius * (mu * mu - 1.0) + ER * ER >= 0.0;
}

static inline double get_layer_density(const DensityProfileLayer *layer,
                                const double altitude) {
  double density = layer->exp_term * exp(layer->exp_scale * altitude) +
                   layer->linear_term * altitude + layer->constant_term;
  return fmax(0.0, fmin(1.0, density));
}

static inline double get_profile_density(const DensityProfile *profile,
                                  const double altitude) {
  if (altitude < profile->layers[0].width) {
    return get_layer_density(&(profile->layers[0]), altitude);
  } else {
    return get_layer_density(&(profile->layers[1]), altitude);
  }
}

static int compute_optical_length_to_space_mie(DATARRAY *mdp, const double radius,
                                               const double mu, double *result) {
  // Number of intervals for the numerical integration.
  const int nsamp = 500;
  // The integration step, i.e. the length of each integration interval.
  const double dx = distance_to_space(radius, mu) / nsamp;
  // Integration loop.
  for (int i = 0; i <= nsamp; ++i) {
    double d_i = i * dx;
    // Distance between the current sample point and the planet center.
    double r_i = sqrt(d_i * d_i + 2.0 * radius * mu * d_i + radius * radius);
    double pt[1] = {r_i - ER};
    DATARRAY *mie = datavector(mdp, pt);
    double weight_i = i == 0 || i == nsamp ? 0.5 : 1.0;
    for (int j = 0; j < NSSAMP; ++j) {
      result[j] += mie->arr.d[j] * weight_i * dx;
    }
    free(mie);
  }
  return 1;
}

static double compute_optical_length_to_space(const DensityProfile *profile,
                                              const double radius, const double mu) {
  // Number of intervals for the numerical integration.
  const int SAMPLE_COUNT = 500;
  // The integration step, i.e. the length of each integration interval.
  const double dx = distance_to_space(radius, mu) / SAMPLE_COUNT;
  // Integration loop.
  double result = 0.0;
  for (int i = 0; i <= SAMPLE_COUNT; ++i) {
    double d_i = i * dx;
    // Distance between the current sample point and the planet center.
    double r_i = sqrt(d_i * d_i + 2.0 * radius * mu * d_i + radius * radius);
    // Number density at the current sample point (divided by the number density
    // at the bottom of the atmosphere, yielding a dimensionless number).
    double y_i = get_profile_density(profile, r_i - ER);
    // Sample weight (from the trapezoidal rule).
    double weight_i = i == 0 || i == SAMPLE_COUNT ? 0.5 : 1.0;
    result += y_i * weight_i * dx;
  }
  return result;
}

static void compute_transmittance_to_space(const Atmosphere *atmos,
                                           const double radius, const double mu,
                                           float *result) {
  const double taur =
      compute_optical_length_to_space(&atmos->rayleigh_density, radius, mu);
  double taum[NSSAMP] = {0};
  compute_optical_length_to_space_mie(atmos->beta_m, radius, mu, taum);
  const double tauo =
      compute_optical_length_to_space(&atmos->ozone_density, radius, mu);
  for (int i = 0; i < NSSAMP; ++i) {
    result[i] = exp(-(taur * atmos->beta_r0[i] + taum[i] * atmos->beta_scale +
                      tauo * AO0[i]));
  }
}

static inline double get_texture_coord_from_unit_range(const double x,
                                                       const int texture_size) {
  return 0.5 / (double)texture_size + x * (1.0 - 1.0 / (double)texture_size);
}

static inline double get_unit_range_from_texture_coord(const double u,
                                                       const int texture_size) {
  return (u - 0.5 / (double)texture_size) / (1.0 - 1.0 / (double)texture_size);
}

static void to_transmittance_uv(const double radius, const double mu, double *u,
                                double *v) {
  // Distance to top atmosphere boundary for a horizontal ray at ground level.
  double H = sqrt(AR * AR - ER * ER);
  // Distance to the horizon.
  double rho = safe_sqrt(radius * radius - ER * ER);
  // Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
  // and maximum values over all mu - obtained for (r,1) and (r,mu_horizon).
  double d = distance_to_space(radius, mu);
  double d_min = AR - radius;
  double d_max = rho + H;
  double x_mu = (d - d_min) / (d_max - d_min);
  double x_r = rho / H;
  *u = fmax(0.0, fmin(1.0, x_mu));
  *v = fmax(0.0, fmin(1.0, x_r));
}

static void from_transmittance_uv(const double u, const double v, double *radius,
                                  double *mu) {
  double x_mu = u;
  double x_r = v;
  // Distance to top atmosphere boundary for a horizontal ray at ground level.
  double H = sqrt(AR * AR - ER * ER);
  // Distance to the horizon, from which we can compute r:
  double rho = H * x_r;
  *radius = sqrt(rho * rho + ER * ER);
  // Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
  // and maximum values over all mu - obtained for (r,1) and (r,mu_horizon) -
  // from which we can recover mu:
  double d_min = AR - *radius;
  double d_max = rho + H;
  double d = d_min + x_mu * (d_max - d_min);
  *mu = d == 0.0 ? 1.0 : (H * H - rho * rho - d * d) / (2.0 * *radius * d);
  *mu = clamp_cosine(*mu);
}

static DATARRAY *get_transmittance_to_space(DATARRAY *dp, const double radius,
                                            const double mu) {
  double u, v;
  to_transmittance_uv(radius, mu, &u, &v);

  double pt[2] = {v, u};
  return datavector(dp, pt);
}

static void get_transmittance(DATARRAY *tau_dp, const double radius, const double mu, const double d,
                              const int intersects_ground, double *result) {

  DATARRAY *result1;
  DATARRAY *result2;

  const double r_d = clamp_radius(sqrt(d * d + 2.0 * radius * mu * d + radius * radius));
  const double mu_d = clamp_cosine((radius * mu + d) / r_d);

  if (intersects_ground) {
    result1 = get_transmittance_to_space(tau_dp, r_d, -mu_d);
    result2 = get_transmittance_to_space(tau_dp, radius, -mu);
  } else {
    result1 = get_transmittance_to_space(tau_dp, radius, mu);
    result2 = get_transmittance_to_space(tau_dp, r_d, mu_d);
  }
  for (int i = 0; i < NSSAMP; ++i) {
    result[i] = fmin(result1->arr.d[i] / result2->arr.d[i], 1.0);
  }
  free(result1);
  free(result2);
}

void get_transmittance_to_sun(DATARRAY *tau_dp, const double r,
                                     const double mu_s, double *result) {
  const double sin_theta_h = ER / r;
  const double cos_theta_h = -sqrt(fmax(1.0 - sin_theta_h * sin_theta_h, 0.0));
  DATARRAY *tau_sun = get_transmittance_to_space(tau_dp, r, mu_s);
  // Smooth step
  const double edge0 = -sin_theta_h * SUNRAD;
  const double edge1 = edge0 * -1;
  double x = mu_s - cos_theta_h;
  x = (x - edge0) / (edge1 - edge0);
  x = x < 0.0 ? 0.0 : x > 1.0 ? 1.0 : x;
  const double st = x * x * (3.0 - 2.0 * x);
  for (int i = 0; i < NSSAMP; ++i) {
    result[i] = tau_sun->arr.d[i] * st;
  }
  free(tau_sun);
}

static void compute_single_scattering_integrand(
    const Atmosphere *atmos, DATARRAY *tau_dp, const double radius, const double mu,
    const double mu_s, const double nu, const double distance,
    const int r_mu_intersects_ground, double *rayleigh, double *mie) {

  const double radius_d = clamp_radius(sqrt(distance * distance + 2.0 * radius * mu * distance + radius * radius));
  const double mu_s_d = clamp_cosine((radius * mu_s + distance * nu) / radius_d);
  // double transmittance[NSSAMP];
  double tau_r_mu[NSSAMP] = {0};
  double tau_sun[NSSAMP] = {0};
  get_transmittance(tau_dp, radius, mu, distance, r_mu_intersects_ground, tau_r_mu);
  get_transmittance_to_sun(tau_dp, radius_d, mu_s_d, tau_sun);
  const double rayleigh_profile_density =
      get_profile_density(&atmos->rayleigh_density, radius_d - ER);
  DATARRAY *mie_scat;
  double pt[1] = {radius_d - ER};
  mie_scat = datavector(atmos->beta_m, pt);
  for (int i = 0; i < NSSAMP; ++i) {
    double transmittance = tau_r_mu[i] * tau_sun[i];
    rayleigh[i] = transmittance * rayleigh_profile_density;
    mie[i] = transmittance * mie_scat->arr.d[i];
  }
  free(mie_scat);
}

static double
distance_to_nearst_atmosphere_boundary(const double radius, const double mu,
                                       const int hits_ground) {
  if (hits_ground) {
    return distance_to_earth(radius, mu);
  } else {
    return distance_to_space(radius, mu);
  }
}

static void compute_single_scattering(const Atmosphere *atmos, DATARRAY *tau_dp,
                                      const double radius, const double mu,
                                      const double mu_s, const double nu,
                                      const int hits_ground,
                                      float *rayleigh, float *mie) {

  const int nsamp = 50;
  double distance_to_boundary;

  if (hits_ground) {
    distance_to_boundary = distance_to_earth(radius, mu);
  } else {
    distance_to_boundary = distance_to_space(radius, mu);
  }
  const double dx = distance_to_boundary / (double)nsamp;

  double rayleigh_sum[NSSAMP] = {0};
  double mie_sum[NSSAMP] = {0};

  for (int i = 0; i <= nsamp; ++i) {
    const double d_i = i * dx;
    double rayleigh_i[NSSAMP] = {0};
    double mie_i[NSSAMP] = {0};
    compute_single_scattering_integrand(atmos, tau_dp, radius, mu, mu_s, nu, d_i,
                                        hits_ground, rayleigh_i,
                                        mie_i);
    double weight_i = (i == 0 || i == nsamp) ? 0.5 : 1.0;
    for (int j = 0; j < NSSAMP; ++j) {
      rayleigh_sum[j] += rayleigh_i[j] * weight_i;
      mie_sum[j] += mie_i[j] * weight_i;
    }
  }
  for (int i = 0; i < NSSAMP; ++i) {
    rayleigh[i] = rayleigh_sum[i] * dx * EXTSOL[i] * atmos->beta_r0[i];
    mie[i] = dx * EXTSOL[i] * (mie_sum[i] * atmos->beta_scale);
  }
}

inline static double rayleigh_phase_function(const double nu) {
  const double k = 3.0 / (16.0 * PI);
  return k * (1.0 + nu * nu);
}

inline static double mie_phase_function(const double g, const double nu) {
  const double k = 3.0 / (8.0 * PI) * (1.0 - g * g) / (2.0 + g * g);
  return k * (1.0 + nu * nu) / pow(1.0 + g * g - 2.0 * g * nu, 1.5);
}

inline static double cloud_phase_function(const double g, const double nu) {
  const double alpha = 0.5;
  const double g1 = -0.5;
  const double hg0 = mie_phase_function(g, nu);
  const double hg1 = mie_phase_function(g1, nu);
  return alpha * hg0 + (1.0 - alpha) * hg1;
}

static void to_scattering_uvwz(const double radius, const double mu, const double mu_s, const double nu,
                        const int r_mu_hits_ground, double *u, double *v,
                        double *w, double *z) {

  double H = sqrt(AR * AR - ER * ER);
  double rho = safe_sqrt(radius * radius - ER * ER);
  double u_r = rho / H;

  double r_mu = radius * mu;
  double discriminant = r_mu * r_mu - radius * radius + ER * ER;

  double u_mu;
  if (r_mu_hits_ground) {
    double d = -r_mu - safe_sqrt(discriminant);
    double d_min = radius - ER;
    double d_max = rho;
    u_mu = 0.5 - 0.5 * get_texture_coord_from_unit_range(
                           d_max == d_min ? 0.0 : (d - d_min) / (d_max - d_min),
                           SCATTERING_TEXTURE_MU_SIZE / 2);
  } else {
    double d = -r_mu + sqrt(discriminant + H * H);
    double d_min = AR - radius;
    double d_max = rho + H;
    u_mu = 0.5 + 0.5 * get_texture_coord_from_unit_range(
                           (d - d_min) / (d_max - d_min),
                           SCATTERING_TEXTURE_MU_SIZE / 2);
  }

  double d = distance_to_space(ER, mu_s);
  double d_min = AH;
  double d_max = H;
  double a = (d - d_min) / (d_max - d_min);
  double D = distance_to_space(ER, MU_S_MIN);
  double A = (D - d_min) / (d_max - d_min);
  double u_mu_s = get_texture_coord_from_unit_range(
      fmax(1.0 - a / A, 0.0) / (1.0 + a), SCATTERING_TEXTURE_MU_S_SIZE);
  double u_nu = (nu + 1.0) / 2;
  *u = u_nu;
  *v = u_mu_s;
  *w = u_mu;
  *z = u_r;
}

static void from_scattering_uvwz(const double x, const double y, const double z, const double w,
                                 double *radius, double *mu, double *mu_s,
                                 double *nu, int *r_mu_hits_ground) {


  const double H = sqrt(AR * AR - ER * ER);
  const double rho = H * w;
  *radius = sqrt(rho * rho + ER * ER);

  if (z < 0.5) {
    // Distance to the ground for the ray (r,mu), and its minimum and maximum
    // values over all mu - obtained for (r,-1) and (r,mu_horizon) - from which
    // we can recover mu:
    const double d_min = *radius - ER;
    const double d_max = rho;
    const double d = d_min + (d_max - d_min) *
                           get_unit_range_from_texture_coord(
                               1.0 - 2.0 * z, SCATTERING_TEXTURE_MU_SIZE / 2);
    *mu = d == 0.0 ? -1.0 : clamp_cosine(-(rho * rho + d * d) / (2.0 * *radius * d));
    *r_mu_hits_ground = 1;
  } else {
    // Distance to the top atmosphere boundary for the ray (r,mu), and its
    // minimum and maximum values over all mu - obtained for (r,1) and
    // (r,mu_horizon) - from which we can recover mu:
    const double d_min = AR - *radius;
    const double d_max = rho + H;
    const double d = d_min + (d_max - d_min) *
                           get_unit_range_from_texture_coord(
                               2.0 * z - 1.0, SCATTERING_TEXTURE_MU_SIZE / 2);
    *mu = d == 0.0
              ? 1.0
              : clamp_cosine((H * H - rho * rho - d * d) / (2.0 * (*radius) * d));
    *r_mu_hits_ground = 0;
  }

  const double x_mu_s = y;
  const double d_min = AH;
  const double d_max = H;
  const double D = distance_to_space(ER, MU_S_MIN);
  const double A = (D - d_min) / (d_max - d_min);
  const double a = (A - x_mu_s * A) / (1.0 + x_mu_s * A);
  const double d = d_min + fmin(a, A) * (d_max - d_min);
  *mu_s = d == 0.0 ? 1.0 : clamp_cosine((H * H - d * d) / (2.0 * ER * d));
  *nu = clamp_cosine(x * 2.0 - 1.0);
}

static DATARRAY *interpolate_scattering(DATARRAY *dp, const double radius,
                                        const double mu, const double mu_s,
                                        const double nu,
                                        const int r_mu_hits_ground) {
  double u, v, w, z;
  to_scattering_uvwz(radius, mu, mu_s, nu, r_mu_hits_ground, &u, &v, &w,
                     &z);
  double pt[4] = {z, w, v, u};
  return datavector(dp, pt);
}

static void get_scattering(DATARRAY *scat1r, DATARRAY *scat1m, DATARRAY *scat,
                           const double radius, const double mu, const double mu_s,
                           const double nu,
                           const int ray_r_mu_intersects_ground,
                           const int scattering_order, double *result) {
  if (scattering_order == 1) {
    DATARRAY *rayleigh = interpolate_scattering(scat1r, radius, mu, mu_s, nu,
                                                ray_r_mu_intersects_ground);
    DATARRAY *mie = interpolate_scattering(scat1m, radius, mu, mu_s, nu,
                                           ray_r_mu_intersects_ground);
    const double rayleigh_phase = rayleigh_phase_function(nu);
    const double mie_phase = mie_phase_function(MIE_G, nu);
    for (int i = 0; i < NSSAMP; ++i) {
      result[i] =
          rayleigh->arr.d[i] * rayleigh_phase + mie->arr.d[i] * mie_phase;
    }
    free(rayleigh);
    free(mie);
  } else {
    DATARRAY *scattering = interpolate_scattering(scat, radius, mu, mu_s, nu,
                                                  ray_r_mu_intersects_ground);
    for (int i = 0; i < NSSAMP; ++i) {
      result[i] = scattering->arr.d[i];
    }
    free(scattering);
  }
}

static void to_irradiance_uv(const double radius, const double mu_s, double *u,
                             double *v) {
  const double x_r = (radius - ER) / AH;
  const double x_mu_s = mu_s * 0.5 + 0.5;
  *u = x_mu_s;
  *v = x_r;
}

static void from_irradiance_uv(const double u, const double v, double *radius,
                               double *mu_s) {
  *radius = v * AH + ER;
  *mu_s = clamp_cosine(u * 2.0 - 1.0);
}

DATARRAY *get_indirect_irradiance(DATARRAY *dp, const double radius,
                                const double mu_s) {
  double u, v;
  to_irradiance_uv(radius, mu_s, &u, &v);
  double pt[2] = {v, u};
  return datavector(dp, pt);
}

static void
compute_scattering_density(const Atmosphere *atmos, DATARRAY *tau_dp,
                           DATARRAY *scat1r, DATARRAY *scat1m, DATARRAY *scat,
                           DATARRAY *irrad_dp, const double radius, const double mu,
                           const double mu_s, const double nu,
                           const int scattering_order, double *result) {
  // Compute unit direction vectors for the zenith, the view direction omega and
  // and the sun direction omega_s, such that the cosine of the view-zenith
  // angle is mu, the cosine of the sun-zenith angle is mu_s, and the cosine of
  // the view-sun angle is nu. The goal is to simplify computations below.
  const double height = radius - ER;
  const double epsilon = 1e-6;
  const FVECT zenith_direction = {0.0, 0.0, 1.0};
  double omega_0 = 1.0 - mu * mu;
  omega_0 = fabs(omega_0) <= epsilon ? 0.0 : sqrt(omega_0);
  const FVECT omega = {omega_0, 0.0, mu};
  const double sun_dir_x = omega[0] == 0.0 ? 0.0 : (nu - mu * mu_s) / omega[0];
  const double sun_dir_y =
      sqrt(fmax(1.0 - sun_dir_x * sun_dir_x - mu_s * mu_s, 0.0));
  const FVECT omega_s = {sun_dir_x, sun_dir_y, mu_s};

  const int SAMPLE_COUNT = 16;
  const double dphi = PI / SAMPLE_COUNT;
  const double dtheta = PI / SAMPLE_COUNT;

  // Nested loops for the integral over all the incident directions omega_i.
  for (int l = 0; l < SAMPLE_COUNT; ++l) {
    const double theta = (l + 0.5) * dtheta;
    const double cos_theta = cos(theta);
    const double sin_theta = sin(theta);
    const int r_theta_hits_ground = ray_hits_ground(radius, cos_theta);

    // The distance and transmittance to the ground only depend on theta, so we
    // can compute them in the outer loop for efficiency.
    double distance_to_ground = 0.0;
    double transmittance_to_ground[NSSAMP] = {0};
    double ground_albedo = 0;
    if (r_theta_hits_ground) {
      distance_to_ground = distance_to_earth(radius, cos_theta);
      // assert(distance_to_ground >= 0.0 && distance_to_ground <= HMAX);
      get_transmittance(tau_dp, radius, cos_theta, distance_to_ground, r_theta_hits_ground,
                        transmittance_to_ground);
      ground_albedo = atmos->grefl;
    }

    for (int m = 0; m < 2 * SAMPLE_COUNT; ++m) {
      const double phi = (m + 0.5) * dphi;
      const FVECT omega_i = {cos(phi) * sin_theta, sin(phi) * sin_theta,
                             cos_theta};
      const double domega_i = dtheta * dphi * sin(theta);

      // The radiance L_i arriving from direction omega_i after n-1 bounces is
      // the sum of a term given by the precomputed scattering texture for the
      // (n-1)-th order:
      double nu1 = fdot(omega_s, omega_i);
      if (nu1 <= -1.0) {
        nu1 += 0.1;
      } else if (nu1 >= 1.0) {
        nu1 -= 0.1;
      }
      double incident_radiance[NSSAMP] = {0};
      get_scattering(scat1r, scat1m, scat, radius, omega_i[2], mu_s, nu1,
                     r_theta_hits_ground, scattering_order - 1,
                     incident_radiance);

      // and of the contribution from the light paths with n-1 bounces and whose
      // last bounce is on the ground. This contribution is the product of the
      // transmittance to the ground, the ground albedo, the ground BRDF, and
      // the irradiance received on the ground after n-2 bounces.
      FVECT ground_normal = {
          zenith_direction[0] * radius + omega_i[0] * distance_to_ground,
          zenith_direction[1] * radius + omega_i[1] * distance_to_ground,
          zenith_direction[2] * radius + omega_i[2] * distance_to_ground};
      normalize(ground_normal);
      DATARRAY *ground_irradiance =
          get_indirect_irradiance(irrad_dp, ER, fdot(ground_normal, omega_s));
      for (int i = 0; i < NSSAMP; ++i) {
        incident_radiance[i] += transmittance_to_ground[i] * ground_albedo *
                                (1.0 / PI) * ground_irradiance->arr.d[i];
      }
      free(ground_irradiance);

      // The radiance finally scattered from direction omega_i towards direction
      // -omega is the product of the incident radiance, the scattering
      // coefficient, and the phase function for directions omega and omega_i
      // (all this summed over all particle types, i.e. Rayleigh and Mie).
      const double nu2 = fdot(omega, omega_i);
      const double rayleigh_density =
          get_profile_density(&atmos->rayleigh_density, height);
      DATARRAY *mie_scat;
      double pt[1] = {height};
      mie_scat = datavector(atmos->beta_m, pt);
      const double rayleigh_phase = rayleigh_phase_function(nu2);
      const double mie_phase = mie_phase_function(MIE_G, nu2);
      for (int j = 0; j < NSSAMP; ++j) {
        result[j] += incident_radiance[j] * domega_i *
                     (atmos->beta_r0[j] * rayleigh_density * rayleigh_phase +
                      mie_scat->arr.d[j] * mie_phase * atmos->beta_scale);
      }
      free(mie_scat);
    }
  }
}

static void compute_multi_scattering(DATARRAY *tau_dp,
                                        DATARRAY *scattering_density_dp,
                                        const double radius, const double mu,
                                        const double mu_s, const double nu,
                                        const int r_mu_hits_ground,
                                        double *result) {


  // Numerical integration sample count.
  const int nsamp = 50;
  const double dx = distance_to_nearst_atmosphere_boundary(
                        radius, mu, r_mu_hits_ground) /
                    (double)nsamp;

  for (int i = 0; i <= nsamp; ++i) {
    const double d_i = i * dx;

    // The r, mu and mu_s parameters at the current integration point (see the
    // single scattering section for a detailed explanation).
    const double r_i =
        clamp_radius(sqrt(d_i * d_i + 2.0 * radius * mu * d_i + radius * radius));
    const double mu_i = clamp_cosine((radius * mu + d_i) / r_i);
    const double mu_s_i = clamp_cosine((radius * mu_s + d_i * nu) / r_i);

    // The Rayleigh and Mie multiple scattering at the current sample point.
    double transmittance[NSSAMP] = {0};
    DATARRAY *rayleigh_mie_s =
        interpolate_scattering(scattering_density_dp, r_i, mu_i, mu_s_i, nu,
                               r_mu_hits_ground);
    get_transmittance(tau_dp, radius, mu, d_i, r_mu_hits_ground,
                      transmittance);
    // Sample weight (from the trapezoidal rule).
    const double weight_i = (i == 0 || i == nsamp) ? 0.5 : 1.0;
    for (int j = 0; j < NSSAMP; ++j) {
      result[j] += rayleigh_mie_s->arr.d[j] * transmittance[j] * dx * weight_i;
    }
    free(rayleigh_mie_s);
  }
}

static void compute_direct_irradiance(DATARRAY *tau_dp, const double r,
                                      const double mu_s, float *result) {

  // Approximate the average of the cosine factor mu_s over the visible fraction
  // of the sun disc
  const double average_cosine_factor =
      mu_s < -SUNRAD
          ? 0.0
          : (mu_s > SUNRAD
                 ? mu_s
                 : (mu_s + SUNRAD) * (mu_s + SUNRAD) / (4.0 * SUNRAD));

  DATARRAY *transmittance = get_transmittance_to_space(tau_dp, r, mu_s);
  for (int i = 0; i < NSSAMP; ++i) {
    result[i] = EXTSOL[i] * transmittance->arr.d[i] * average_cosine_factor;
  }
  free(transmittance);
}

static void compute_indirect_irradiance(DATARRAY *scat1r, DATARRAY *scat1m,
                                        DATARRAY *scat, const double radius,
                                        const double mu_s,
                                        const int scattering_order,
                                        double *result) {

  const int SAMPLE_COUNT = 32;
  const double dphi = PI / SAMPLE_COUNT;
  const double dtheta = PI / SAMPLE_COUNT;

  const FVECT omega_s = {sqrt(1.0 - mu_s * mu_s), 0.0, mu_s};
  for (int j = 0; j < SAMPLE_COUNT / 2; ++j) {
    const double theta = (j + 0.5) * dtheta;
    for (int i = 0; i < 2 * SAMPLE_COUNT; ++i) {
      const double phi = (i + 0.5) * dphi;
      const FVECT omega = {sin(theta) * cos(phi), sin(theta) * sin(phi),
                           cos(theta)};
      const double domega = dtheta * dphi * sin(theta);
      const double nu = fdot(omega, omega_s);
      double result1[NSSAMP] = {0};
      get_scattering(scat1r, scat1m, scat, radius, omega[2], mu_s, nu, 0,
                     scattering_order, result1);
      for (int k = 0; k < NSSAMP; ++k) {
        result[k] += result1[k] * omega[2] * domega;
      }
    }
  }
}

DATARRAY *allocate_3d_datarray(const char *name, const int ri, const int rj,
                               const int rk) {
  DATARRAY *dp;
  if ((dp = (DATARRAY *)malloc(sizeof(DATARRAY))) == NULL)
    goto memerr;
  int asize = ri * rj * rk;
  dp->name = savestr(name);
  dp->type = DATATY;
  dp->nd = 3;
  dp->dim[0].org = 0;
  dp->dim[0].siz = 1;
  dp->dim[0].ne = ri;
  dp->dim[0].p = NULL;
  dp->dim[1].org = 0;
  dp->dim[1].siz = 1;
  dp->dim[1].ne = rj;
  dp->dim[1].p = NULL;
  dp->dim[2].org = 0;
  dp->dim[2].siz = rk - 1;
  dp->dim[2].ne = rk;
  dp->dim[2].p = NULL;
  // dp->arr.p =
  // malloc(sizeof(DATATYPE)*dp->dim[0].ne*dp->dim[1].ne*dp->dim[2].ne);
  if ((dp->arr.d = (DATATYPE *)malloc(asize * sizeof(DATATYPE))) == NULL)
    goto memerr;
  return (dp);

memerr:
  fprintf(stderr, "Memory allocation error in allocate_3d_datarray\n");
  return (NULL);
}

DATARRAY *allocate_5d_datarray(const char *name, const int ri, const int rj,
                               const int rk, const int rl, const int rm) {
  DATARRAY *dp;
  int asize = ri * rj * rk * rl * rm;
  if ((dp = (DATARRAY *)malloc(sizeof(DATARRAY))) == NULL)
    goto memerr;
  dp->name = savestr(name);
  dp->type = DATATY;
  dp->nd = 5;
  dp->dim[0].org = 0;
  dp->dim[0].siz = 1;
  dp->dim[0].ne = ri;
  dp->dim[0].p = NULL;
  dp->dim[1].org = 0;
  dp->dim[1].siz = 1;
  dp->dim[1].ne = rj;
  dp->dim[1].p = NULL;
  dp->dim[2].org = 0;
  dp->dim[2].siz = 1;
  dp->dim[2].ne = rk;
  dp->dim[2].p = NULL;
  dp->dim[3].org = 0;
  dp->dim[3].siz = 1;
  dp->dim[3].ne = rl;
  dp->dim[3].p = NULL;
  dp->dim[4].org = 0;
  dp->dim[4].siz = rm - 1;
  dp->dim[4].ne = rm;
  dp->dim[4].p = NULL;
  if ((dp->arr.d = (DATATYPE *)malloc(asize * sizeof(DATATYPE))) == NULL)
    goto memerr;
  return (dp);
memerr:
  fprintf(stderr, "Memory allocation error in allocate_5d_datarray\n");
  return (NULL);
}

/* save to a .dat file */
void savedata(DATARRAY *dp) {
  if (dp == NULL)
    return;
  FILE *fp = fopen(dp->name, "w");
  if (fp == NULL) {
    fprintf(stderr, "Error opening file %s\n", dp->name);
    return;
  }
  int i, j;
  int nvals = 1;
  for (i = 0; i < dp->nd; i++) {
    nvals *= dp->dim[i].ne;
  }
  fprintf(fp, "%d\n", dp->nd);
  for (i = 0; i < dp->nd; i++) {
    fprintf(fp, "%f %f %d\n", dp->dim[i].org, dp->dim[i].siz, dp->dim[i].ne);
  }
  for (i = 0; i < nvals; i++) {
    fprintf(fp, "%f\n", dp->arr.d[i]);
  }
  fclose(fp);
}

void increment_dp(DATARRAY *dp1, DATARRAY *dp2) {
  if (dp1 == NULL || dp2 == NULL)
    perror("null pointer in increment_dp\n");

  if (dp1->nd != dp2->nd)
    perror("dimension mismatch in increment_dp\n");

  int i;
  int nvals1 = 1;
  int nvals2 = 1;
  for (i = 0; i < dp1->nd; i++) {
    nvals1 *= dp1->dim[i].ne;
    nvals2 *= dp2->dim[i].ne;
  }
  if (nvals1 != nvals2)
    perror("size mismatch in increment_dp\n");
  for (i = 0; i < nvals1; i++) {
    dp1->arr.d[i] += dp2->arr.d[i];
  }
}

THREAD_RETURN compute_scattering_density_thread(void *arg) {
  ScatDenseTdat *tdata = (ScatDenseTdat *)arg;
  for (unsigned int l = tdata->start_l; l < tdata->end_l; ++l) {
    for (unsigned int k = 0; k < SCATTERING_TEXTURE_MU_SIZE; ++k) {
      for (unsigned int j = 0; j < SCATTERING_TEXTURE_MU_S_SIZE; ++j) {
        for (unsigned int i = 0; i < SCATTERING_TEXTURE_NU_SIZE; ++i) {
          double scattering_density[NSSAMP] = {0};
          double r, mu, mu_s, nu;
          int ray_r_mu_intersects_ground;
          double xr = (double)i / SCATTERING_TEXTURE_NU_SIZE;
          double yr = (double)j / SCATTERING_TEXTURE_MU_S_SIZE;
          double zr = (double)k / SCATTERING_TEXTURE_MU_SIZE;
          double wr = (double)l / SCATTERING_TEXTURE_R_SIZE;
          from_scattering_uvwz(xr, yr, zr, wr, &r, &mu, &mu_s, &nu,
                               &ray_r_mu_intersects_ground);
          nu =
              fmax(mu * mu_s - sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)),
                   fmin(mu * mu_s + sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)),
                        nu));
          compute_scattering_density(
              tdata->atmos, tdata->tau_dp, tdata->delta_rayleigh_scattering_dp,
              tdata->delta_mie_scattering_dp,
              tdata->delta_multiple_scattering_dp, tdata->delta_irradiance_dp,
              r, mu, mu_s, nu, tdata->scattering_order, scattering_density);
          for (int m = 0; m < NSSAMP; ++m) {
            int idx = l * SCATTERING_TEXTURE_MU_SIZE *
                          SCATTERING_TEXTURE_MU_S_SIZE *
                          SCATTERING_TEXTURE_NU_SIZE * NSSAMP +
                      k * SCATTERING_TEXTURE_MU_S_SIZE *
                          SCATTERING_TEXTURE_NU_SIZE * NSSAMP +
                      j * SCATTERING_TEXTURE_NU_SIZE * NSSAMP + i * NSSAMP + m;
            tdata->delta_scattering_density_dp->arr.d[idx] =
                scattering_density[m];
          }
        }
      }
    }
  }
  return 0;
}

THREAD_RETURN compute_single_scattering_thread(void *arg) {
  Scat1Tdat *tdata = (Scat1Tdat *)arg;
  for (unsigned int l = tdata->start_l; l < tdata->end_l; ++l) {
    for (int k = 0; k < SCATTERING_TEXTURE_MU_SIZE; ++k) {
      for (int j = 0; j < SCATTERING_TEXTURE_MU_S_SIZE; ++j) {
        for (int i = 0; i < SCATTERING_TEXTURE_NU_SIZE; ++i) {
          float rayleigh[NSSAMP] = {0};
          float mie[NSSAMP] = {0};
          double r, mu, mu_s, nu;
          int r_mu_hits_ground;
          const double xr = (double)i / SCATTERING_TEXTURE_NU_SIZE;
          const double yr = (double)j / SCATTERING_TEXTURE_MU_S_SIZE;
          const double zr = (double)k / SCATTERING_TEXTURE_MU_SIZE;
          const double wr = (double)l / SCATTERING_TEXTURE_R_SIZE;
          from_scattering_uvwz(xr, yr, zr, wr, &r, &mu, &mu_s, &nu,
                               &r_mu_hits_ground);
          const double mu_mu_s = mu * mu_s;
          const double sqrt_term = sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s));
          nu = fmax(mu_mu_s - sqrt_term, fmin(mu_mu_s + sqrt_term, nu));
          compute_single_scattering(tdata->atmos, tdata->tau_dp, r, mu, mu_s,
                                    nu, r_mu_hits_ground, rayleigh,
                                    mie);
          for (int m = 0; m < NSSAMP; ++m) {
            int idx = l * SCATTERING_TEXTURE_MU_SIZE *
                          SCATTERING_TEXTURE_MU_S_SIZE *
                          SCATTERING_TEXTURE_NU_SIZE * NSSAMP +
                      k * SCATTERING_TEXTURE_MU_S_SIZE *
                          SCATTERING_TEXTURE_NU_SIZE * NSSAMP +
                      j * SCATTERING_TEXTURE_NU_SIZE * NSSAMP + i * NSSAMP + m;
            tdata->delta_rayleigh_scattering_dp->arr.d[idx] = rayleigh[m];
            tdata->delta_mie_scattering_dp->arr.d[idx] = mie[m];
            tdata->scattering_dp->arr.d[idx] = rayleigh[m];
          }
        }
      }
    }
  }
  return 0;
}

THREAD_RETURN compute_multi_scattering_thread(void *arg) {
  ScatNTdat *tdata = (ScatNTdat *)arg;
  for (unsigned int l = tdata->start_l; l < tdata->end_l; ++l) {
    for (unsigned int k = 0; k < SCATTERING_TEXTURE_MU_SIZE; ++k) {
      for (unsigned int j = 0; j < SCATTERING_TEXTURE_MU_S_SIZE; ++j) {
        for (unsigned int i = 0; i < SCATTERING_TEXTURE_NU_SIZE; ++i) {
          double delta_multiple_scattering[NSSAMP] = {0};
          double r, mu, mu_s, nu;
          int r_mu_hits_ground;
          const double xr = (double)i / SCATTERING_TEXTURE_NU_SIZE;
          const double yr = (double)j / SCATTERING_TEXTURE_MU_S_SIZE;
          const double zr = (double)k / SCATTERING_TEXTURE_MU_SIZE;
          const double wr = (double)l / SCATTERING_TEXTURE_R_SIZE;
          from_scattering_uvwz(xr, yr, zr, wr, &r, &mu, &mu_s, &nu,
                               &r_mu_hits_ground);
          const double mu_mu_s = mu * mu_s;
          const double sqrt_term = sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s));
          nu = fmax(mu_mu_s - sqrt_term, fmin(mu_mu_s + sqrt_term, nu));
          compute_multi_scattering(
              tdata->tau_dp, tdata->delta_scattering_density_dp, r, mu, mu_s,
              nu, r_mu_hits_ground, delta_multiple_scattering);
          const double rayleigh_phase = rayleigh_phase_function(nu);
          for (int m = 0; m < NSSAMP; ++m) {
            int idx = l * SCATTERING_TEXTURE_MU_SIZE *
                          SCATTERING_TEXTURE_MU_S_SIZE *
                          SCATTERING_TEXTURE_NU_SIZE * NSSAMP +
                      k * SCATTERING_TEXTURE_MU_S_SIZE *
                          SCATTERING_TEXTURE_NU_SIZE * NSSAMP +
                      j * SCATTERING_TEXTURE_NU_SIZE * NSSAMP + i * NSSAMP + m;
            tdata->delta_multiple_scattering_dp->arr.d[idx] =
                delta_multiple_scattering[m];
            tdata->scattering_dp->arr.d[idx] +=
                delta_multiple_scattering[m] * (1.0 / rayleigh_phase);
          }
        }
      }
    }
  }
  return 0;
}

int precompute(const int sorder, const DpPaths dppaths, const Atmosphere *atmos,
               int num_threads) {
  num_threads = (num_threads < SCATTERING_TEXTURE_R_SIZE)
                    ? num_threads
                    : SCATTERING_TEXTURE_R_SIZE;
  const int tchunk = SCATTERING_TEXTURE_R_SIZE / num_threads;
  const int tremainder = SCATTERING_TEXTURE_R_SIZE % num_threads;

  if (sorder < 2) {
    fprintf(stderr, "scattering order must be at least 2\n");
    return 0;
  }

  DATARRAY *tau_dp =
      allocate_3d_datarray(dppaths.tau, TRANSMITTANCE_TEXTURE_HEIGHT,
                           TRANSMITTANCE_TEXTURE_WIDTH, NSSAMP);
  DATARRAY *delta_irradiance_dp = allocate_3d_datarray(
      "ssky_delta_irradiance.dat", IRRADIANCE_TEXTURE_HEIGHT,
      IRRADIANCE_TEXTURE_WIDTH, NSSAMP);
  DATARRAY *irradiance_dp =
      allocate_3d_datarray(dppaths.irrad, IRRADIANCE_TEXTURE_HEIGHT,
                           IRRADIANCE_TEXTURE_WIDTH, NSSAMP);

  DATARRAY *delta_rayleigh_scattering_dp = allocate_5d_datarray(
      "ssky_delta_rayleigh_scattering_dp.dat", SCATTERING_TEXTURE_R_SIZE,
      SCATTERING_TEXTURE_MU_SIZE, SCATTERING_TEXTURE_MU_S_SIZE,
      SCATTERING_TEXTURE_NU_SIZE, NSSAMP);
  DATARRAY *delta_mie_scattering_dp = allocate_5d_datarray(
      dppaths.scat1m, SCATTERING_TEXTURE_R_SIZE, SCATTERING_TEXTURE_MU_SIZE,
      SCATTERING_TEXTURE_MU_S_SIZE, SCATTERING_TEXTURE_NU_SIZE, NSSAMP);
  DATARRAY *scattering_dp = allocate_5d_datarray(
      dppaths.scat, SCATTERING_TEXTURE_R_SIZE, SCATTERING_TEXTURE_MU_SIZE,
      SCATTERING_TEXTURE_MU_S_SIZE, SCATTERING_TEXTURE_NU_SIZE, NSSAMP);
  DATARRAY *delta_multiple_scattering_dp = allocate_5d_datarray(
      "ssky_delta_multiple_scattering.dat", SCATTERING_TEXTURE_R_SIZE,
      SCATTERING_TEXTURE_MU_SIZE, SCATTERING_TEXTURE_MU_S_SIZE,
      SCATTERING_TEXTURE_NU_SIZE, NSSAMP);
  DATARRAY *delta_scattering_density_dp = allocate_5d_datarray(
      "ssky_delta_scattering_density.dat", SCATTERING_TEXTURE_R_SIZE,
      SCATTERING_TEXTURE_MU_SIZE, SCATTERING_TEXTURE_MU_S_SIZE,
      SCATTERING_TEXTURE_NU_SIZE, NSSAMP);

  int idx = 0;
  // Computing transmittance...
  for (int j = 0; j < TRANSMITTANCE_TEXTURE_HEIGHT; ++j) {
    for (int i = 0; i < TRANSMITTANCE_TEXTURE_WIDTH; ++i) {
      double r, mu;
      float tau[NSSAMP] = {0};
      const float xr = (float)i / TRANSMITTANCE_TEXTURE_WIDTH;
      const float yr = (float)j / TRANSMITTANCE_TEXTURE_HEIGHT;
      from_transmittance_uv(xr, yr, &r, &mu);
      compute_transmittance_to_space(atmos, r, mu, tau);
      for (int k = 0; k < NSSAMP; ++k) {
        tau_dp->arr.d[idx] = tau[k];
        idx += 1;
      }
    }
  }
  savedata(tau_dp);

  // Computing direct irradiance...
  idx = 0;
  for (int j = 0; j < IRRADIANCE_TEXTURE_HEIGHT; ++j) {
    for (int i = 0; i < IRRADIANCE_TEXTURE_WIDTH; ++i) {
      double r, mu_s;
      float result[NSSAMP] = {0};
      const double xr = (double)i / IRRADIANCE_TEXTURE_WIDTH;
      const double yr = (double)j / IRRADIANCE_TEXTURE_HEIGHT;
      from_irradiance_uv(xr, yr, &r, &mu_s);
      compute_direct_irradiance(tau_dp, r, mu_s, result);
      for (int k = 0; k < NSSAMP; ++k) {
        delta_irradiance_dp->arr.d[idx] = result[k];
        irradiance_dp->arr.d[idx] = 0.0;
        idx += 1;
      }
    }
  }

  // Computing single scattering...
  THREAD *threads = (THREAD *)malloc(num_threads * sizeof(THREAD));
  Scat1Tdat *tdata = (Scat1Tdat *)malloc(num_threads * sizeof(Scat1Tdat));
  for (int i = 0; i < num_threads; i++) {
    tdata[i].atmos = atmos;
    tdata[i].tau_dp = tau_dp;
    tdata[i].delta_rayleigh_scattering_dp = delta_rayleigh_scattering_dp;
    tdata[i].delta_mie_scattering_dp = delta_mie_scattering_dp;
    tdata[i].scattering_dp = scattering_dp;
    tdata[i].start_l = i * tchunk;
    tdata[i].end_l = (i + 1) * tchunk;
    if (i == num_threads - 1)
      tdata[i].end_l += tremainder;
    CREATE_THREAD(&threads[i], compute_single_scattering_thread,
                  (void *)&tdata[i]);
  }

  for (int i = 0; i < num_threads; i++) {
    THREAD_JOIN(threads[i]);
#if defined(_WIN32) || defined(_WIN64)
    THREAD_CLOSE(threads[i]);
#endif
  }
  free(tdata);
  free(threads);
  savedata(delta_mie_scattering_dp);


  // Compute the 2nd, 3rd and 4th order of scattering, in sequence.
  for (int scattering_order = 2; scattering_order <= sorder; ++scattering_order) {
    // Compute the scattering density, and store it in
    // delta_scattering_density_texture.
    THREAD *threads = (THREAD *)malloc(num_threads * sizeof(THREAD));
    ScatDenseTdat *tdata = (ScatDenseTdat *)malloc(num_threads * sizeof(ScatDenseTdat));
    for (int i = 0; i < num_threads; i++) {
      tdata[i].atmos = atmos;
      tdata[i].tau_dp = tau_dp;
      tdata[i].delta_rayleigh_scattering_dp = delta_rayleigh_scattering_dp;
      tdata[i].delta_mie_scattering_dp = delta_mie_scattering_dp;
      tdata[i].delta_multiple_scattering_dp = delta_multiple_scattering_dp;
      tdata[i].delta_irradiance_dp = delta_irradiance_dp;
      tdata[i].scattering_order = scattering_order;
      tdata[i].delta_scattering_density_dp = delta_scattering_density_dp;
      tdata[i].start_l = i * tchunk;
      tdata[i].end_l = (i + 1) * tchunk;
      if (i == num_threads - 1)
        tdata[i].end_l += tremainder;
      CREATE_THREAD(&threads[i], compute_scattering_density_thread,
                    (void *)&tdata[i]);
    }

    for (int i = 0; i < num_threads; i++) {
      THREAD_JOIN(threads[i]);
#if defined(_WIN32) || defined(_WIN64)
      THREAD_CLOSE(threads[i]);
#endif
    }
    free(tdata);
    free(threads);

    // Compute the indirect irradiance, store it in delta_irradiance_texture
    // and accumulate it in irradiance_texture_.
    idx = 0;
    for (unsigned int j = 0; j < IRRADIANCE_TEXTURE_HEIGHT; ++j) {
      for (unsigned int i = 0; i < IRRADIANCE_TEXTURE_WIDTH; ++i) {
        double delta_irradiance[NSSAMP] = {0};
        double r, mu_s;
        const double xr = (double)i / IRRADIANCE_TEXTURE_WIDTH;
        const double yr = (double)j / IRRADIANCE_TEXTURE_HEIGHT;
        from_irradiance_uv(xr, yr, &r, &mu_s);
        compute_indirect_irradiance(delta_rayleigh_scattering_dp,
                                    delta_mie_scattering_dp,
                                    delta_multiple_scattering_dp, r, mu_s,
                                    scattering_order - 1, delta_irradiance);
        for (int k = 0; k < NSSAMP; ++k) {
          delta_irradiance_dp->arr.d[idx] = delta_irradiance[k];
          idx += 1;
        }
      }
    }

    increment_dp(irradiance_dp, delta_irradiance_dp);

    // Computing multiple scattering...
    THREAD *threads2 = (THREAD *)malloc(num_threads * sizeof(THREAD));
    ScatNTdat *tdata2 = (ScatNTdat *)malloc(num_threads * sizeof(ScatNTdat));
    for (int i = 0; i < num_threads; i++) {
      tdata2[i].tau_dp = tau_dp;
      tdata2[i].delta_multiple_scattering_dp = delta_multiple_scattering_dp;
      tdata2[i].delta_scattering_density_dp = delta_scattering_density_dp;
      tdata2[i].scattering_dp = scattering_dp;
      tdata2[i].start_l = i * tchunk;
      tdata2[i].end_l = (i + 1) * tchunk;
      if (i == num_threads - 1)
        tdata2[i].end_l += tremainder;
      CREATE_THREAD(&threads2[i], compute_multi_scattering_thread,
                    (void *)&tdata2[i]);
    }

    for (int i = 0; i < num_threads; i++) {
      THREAD_JOIN(threads2[i]);
#if defined(_WIN32) || defined(_WIN64)
      THREAD_CLOSE(threads2[i]);
#endif
    }
    free(tdata2);
    free(threads2);
  }
  savedata(scattering_dp);
  savedata(irradiance_dp);

  freedata(tau_dp);
  freedata(delta_rayleigh_scattering_dp);
  freedata(delta_mie_scattering_dp);
  freedata(scattering_dp);
  freedata(delta_multiple_scattering_dp);
  freedata(delta_scattering_density_dp);
  freedata(delta_irradiance_dp);
  freedata(irradiance_dp);
  return 1;
}

int compute_sundir(const int year, const int month, const int day,
                   const double hour, const int tsolar, double sundir[3]) {
  double sd, st = hour;

  if (year) {
    /* Michalsky algorithm? */
    double mjd = mjdate(year, month, day, hour);
    if (tsolar)
      sd = msdec(mjd, NULL);
    else
      sd = msdec(mjd, &st);
  } else {
    int jd = jdate(month, day); /* Julian date */
    sd = sdec(jd);              /* solar declination */
    if (!tsolar)                /* get solar time? */
      st = hour + stadj(jd);
  }
  const double altitude = salt(sd, st);
  const double azimuth = sazi(sd, st);
  sundir[0] = -sin(azimuth) * cos(altitude);
  sundir[1] = -cos(azimuth) * cos(altitude);
  sundir[2] = sin(altitude);
  return 1;
}

void get_sky_transmittance(DATARRAY *tau, const double radius, const double mu, float *result) {
  DATARRAY *trans = get_transmittance_to_space(tau, radius, mu);
  for (int i = 0; i < NSSAMP; ++i) {
    result[i] = trans->arr.d[i];
  }
  free(trans);
}

void get_sky_radiance(DATARRAY *scat, DATARRAY *scat1m, const double radius, 
                      const double mu, const double mu_s, const double nu, float *result) {
  double u, v, w, z;
  to_scattering_uvwz(radius, mu, mu_s, nu, 0, &u, &v, &w, &z);
  double pt[4] = {z, w, v, u};

  DATARRAY *scattering = datavector(scat, pt);
  DATARRAY *single_mie_scattering = datavector(scat1m, pt);

  const double rayleigh_phase = rayleigh_phase_function(nu);
  const double mie_phase = mie_phase_function(MIE_G, nu);

  for (int i = 0; i < NSSAMP; ++i) {
    result[i] = scattering->arr.d[i] * rayleigh_phase +
                single_mie_scattering->arr.d[i] * mie_phase;
  }
  free(single_mie_scattering);
  free(scattering);
}

static void get_irradiance(DATARRAY *tau, DATARRAY *irrad, const double radius, 
                             const FVECT point, const FVECT normal, 
                             const FVECT sundir, double *result) {
  double mu_s = fdot(point, sundir) / radius;

  double point_trans_sun[NSSAMP] = {0};
  get_transmittance_to_sun(tau, radius, mu_s, point_trans_sun);

  DATARRAY *indirect_irradiance = get_indirect_irradiance(irrad, radius, mu_s);
  const double sun_ct = fmax(fdot(normal, sundir), 0.0);

  for (int i = 0; i < NSSAMP; ++i) {
    result[i] = indirect_irradiance->arr.d[i] + point_trans_sun[i] * EXTSOL[i] * sun_ct;
  }

  free(indirect_irradiance);
}


void get_ground_radiance(DATARRAY *tau, DATARRAY *scat, DATARRAY *scat1m, DATARRAY *irrad, 
                         const FVECT view_point, const FVECT view_direction, 
                         const double radius, const double mu, const double sun_ct, 
                         const double nu, const double grefl, const FVECT sundir, 
                         float *ground_radiance) {

  const int intersect = ray_hits_ground(radius, mu);

  if (intersect) {
    FVECT point, normal;
    const double distance = distance_to_earth(radius, mu);

    // direct + indirect irradiance
    VSUM(point, view_point, view_direction, distance);
    VCOPY(normal, point);
    normalize(normal);
    double irradiance[NSSAMP] = {0};
    get_irradiance(tau, irrad, ER, point, normal, sundir, irradiance);

    // transmittance between view point and ground point
    double trans[NSSAMP] = {0};
    get_transmittance(tau, radius, mu, distance, intersect, trans);
    if (trans[0] == 0) {printf("trans 0\n");}

    // inscattering
    float inscatter[NSSAMP] = {0}; 
    get_sky_radiance(scat, scat1m, radius, mu, sun_ct, nu, inscatter);

    for (int i = 0; i < NSSAMP; ++i) {
      ground_radiance[i] = inscatter[i] + irradiance[i] * trans[i] * grefl / PI;
    }
  }
}

void get_solar_radiance(DATARRAY *tau, DATARRAY *scat, DATARRAY *scat1m, const FVECT sundir, const double radius, const double sun_ct, double *sun_radiance) {
  double trans_sun[NSSAMP] = {0};
  float sky_radiance_sun[NSSAMP] = {0};
  double u, v, w, z;
  const double nu = fdot(sundir, sundir);
  get_transmittance_to_sun(tau, radius, sun_ct, trans_sun);
  get_sky_radiance(scat, scat1m, radius, sun_ct, sun_ct, nu, sky_radiance_sun);
  for (int i = 0; i < NSSAMP; ++i) {
    sun_radiance[i] = sky_radiance_sun[i] + trans_sun[i] * EXTSOL[i] / SOLOMG;
  }
}
