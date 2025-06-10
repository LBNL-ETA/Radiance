/*
	Glazing system multi-layer optics calculation.
	Source equations: LBNL WINDOW technical documention.
	System calcuation: Section 7.6.1
	Anuglar calculation Section 7.7
	T. Wang
*/

#include "rtmath.h"
#include "color.h"
#include "data.h"

enum {
	NTHETA = 59,
	MAX_NAME = 256,
};

const double DBL_EPSILON = 1E-9;
const double LARGE_DOUBE = 3.40282347E+38;
const double RAD_PER_DEG = 0.0174532925199;
const double COEFFS_TAU_CLEAR[] = {-0.0015, 3.355, -3.84, 1.46, 0.0288};
const double COEFFS_TAU_BRONZE[] = {0.002, 2.813, -2.341, -0.05725, 0.599};
const double COEFFS_RHO_CLEAR[] = {0.999, -0.563, 2.043, -2.532, 1.054};
const double COEFFS_RHO_BRONZE[] = {0.997, -1.868, 6.513, -7.862, 3.225};

const double THETAS[NTHETA] = {
	0., 5., 10., 15., 20., 25., 30., 35., 40.,
	41., 42., 43., 44., 45., 46., 47., 48., 49, 50.,
	51., 52., 53., 54., 55., 56., 57., 58., 59, 60.,
	61., 62., 63., 64., 65., 66., 67., 68., 69, 70.,
	71., 72., 73., 74., 75., 76., 77., 78., 79, 80.,
	81., 82., 83., 84., 85., 86., 87., 88., 89, 90.,
};


char *progname;

typedef struct {
	char *filename;
	int is_mono;
	double thickness_m;

	/* Interpolated data at standard wavelengths */
	double *t_0;
	double *rf_0;
	double *rb_0;

	/* Calculated angular data */
	double *t_lambda_theta;
	double *rf_lambda_theta;
	double *rb_lambda_theta;

} GlazingLayer;


static inline double polynomial_5(double x, const double coeffs[5]) {
	return x * ( x * ( x * (coeffs[4] * x + coeffs[3]) + coeffs[2]) + coeffs[1]) + coeffs[0];
}


/* Calculate rho0 (unpolarized reflectance of a single surface) */
static void rho0_calc(const double *t_0, const double *r_0, const int nwvl, double *rho0) {
	for (int i = 0; i < nwvl; ++i) {
		const double beta = t_0[i] * t_0[i] - r_0[i] * r_0[i] + 2.0 * r_0[i] + 1.0;
		double denominator = 2.0 * (2.0 - r_0[i]);
		if (fabs(denominator) < DBL_EPSILON) 
			denominator = DBL_EPSILON;
		double discriminant = beta * beta - 2.0 * denominator * r_0[i];
		if (discriminant < 0.0) 
			discriminant = 0.0;
		rho0[i] = (beta - sqrt(discriminant)) / denominator;
		if (rho0[i] < 0.0) 
			rho0[i] = 0.0;
		if (rho0[i] > 1.0) 
			rho0[i] = 1.0;
	}
}

/* refraction index calc */
static void n_calc(const double *rho_0, const int nwvl, double *n)
{
	for (int i=0; i < nwvl; ++i) {
		const double sqrt_rho = sqrt(rho_0[i]);
		n[i] = ((1.0 - sqrt_rho) < DBL_EPSILON) ? 1.0 : (1.0 + sqrt_rho) / (1.0 - sqrt_rho);
	}
}

/* Calculate absoprtion coefficient */
static void alpha_calc(const double *t0, const double *r0, const double *rho0,
					const double thickness_m, const int nwvl, double *abs_coeff)
{
	for (int i = 0; i < nwvl; ++i) {
		const double numerator = r0[i] - rho0[i];
		const double denominator = rho0[i] * t0[i];
		if (denominator > DBL_EPSILON && numerator > DBL_EPSILON) {
			const double log_val = log(numerator / denominator);
			abs_coeff[i] = -log_val / thickness_m;
		} else {
			abs_coeff[i] = 1.0;
		}
		if (abs_coeff[i] < 0.0 || isnan(abs_coeff[i]) || isinf(abs_coeff[i])) {
			abs_coeff[i] = 0.0;
		}
	}
}


void angular_monolithic(GlazingLayer *layer, const int nwvl) {
	double *rho_0 = malloc(nwvl * sizeof(double));
	double *refrac = malloc(nwvl * sizeof(double));
	double *abs_coeff = malloc(nwvl * sizeof(double));

	rho0_calc(layer->t_0, layer->rf_0, nwvl, rho_0);
	n_calc(rho_0, nwvl, refrac);
	alpha_calc(layer->t_0, layer->rf_0, rho_0, layer->thickness_m, nwvl, abs_coeff);

	for (int itheta = 0; itheta < NTHETA; ++itheta) {
		double phi = THETAS[itheta] * RAD_PER_DEG;
		double cos_phi = cos(phi);
		double sin_phi = sin(phi);

		for (int iwvl = 0; iwvl < nwvl; ++iwvl) {
			double sin_phi_prime_sq = sin_phi * sin_phi / (refrac[iwvl] * refrac[iwvl]);

			double cos_phi_prime = (sin_phi_prime_sq >= 1.0) ? 0.0 : sqrt(1.0 - sin_phi_prime_sq);

			double rho_lambda;
			if (cos_phi_prime <= DBL_EPSILON) {
				rho_lambda = 1.0;
			} else {
				double nf_cos_phi = refrac[iwvl] * cos_phi;
				double nf_cos_phi_p = refrac[iwvl] * cos_phi_prime;
				double par_num = nf_cos_phi - cos_phi_prime;
				double par_den = nf_cos_phi + cos_phi_prime;
				double per_num = nf_cos_phi_p - cos_phi;
				double per_den = nf_cos_phi_p + cos_phi;
				if (fabs(par_den) < DBL_EPSILON || fabs(per_den) < DBL_EPSILON) {
					rho_lambda = 1.0;
				} else {
					double par_comp = pow(par_num / par_den, 2.0);
					double per_comp = pow(per_num / per_den, 2.0);
					rho_lambda = 0.5 * (par_comp + per_comp);
				}
			}
			if (rho_lambda < 0.0)
				rho_lambda = 0.0;
			if (rho_lambda > 1.0)
				rho_lambda = 1.0;

			double exp_arg = (cos_phi_prime < DBL_EPSILON) ? -LARGE_DOUBE : (-abs_coeff[iwvl] * layer->thickness_m / cos_phi_prime);
			/* Clamp exponent argument to prevent overflow in exp() */
			if (exp_arg < -700.0) 
				exp_arg = -700.0;

			double T_internal = exp(exp_arg);

			double denominator = 1.0 - rho_lambda * rho_lambda * T_internal * T_internal;

			double t_lambda, r_lambda;

			if (fabs(denominator) < DBL_EPSILON || rho_lambda >= 1.0) {
				t_lambda = 0.0;
				r_lambda = 1.0;
			} else {
				double tau_lambda = 1.0 - rho_lambda;
				t_lambda = tau_lambda * tau_lambda * T_internal / denominator;
				r_lambda = rho_lambda * (1.0 + t_lambda * T_internal);

			}
			if (t_lambda < 0.0)
				t_lambda = 0.0;
			if (t_lambda > 1.0)
				t_lambda = 1.0;

			const int flat_idx = itheta * nwvl + iwvl;
			layer->t_lambda_theta[flat_idx] = t_lambda;
			layer->rf_lambda_theta[flat_idx] = r_lambda;
			layer->rb_lambda_theta[flat_idx] = r_lambda;
		}
	}
	free(rho_0);
	free(refrac);
	free(abs_coeff);
}


void angular_coated(GlazingLayer *layer, const int nwvl) {
	for (int itheta = 0; itheta < NTHETA; ++itheta) {
		const double phi = THETAS[itheta] * RAD_PER_DEG;
		const double cos_phi = cos(phi);

		for (int iwvl = 0; iwvl < nwvl; ++iwvl) {
			const double t_0_lambda = layer->t_0[iwvl];
			const double *coeffs_tau = (t_0_lambda > 0.645) ? COEFFS_TAU_CLEAR : COEFFS_TAU_BRONZE;
			const double *coeffs_rho = (t_0_lambda > 0.645) ? COEFFS_RHO_CLEAR : COEFFS_RHO_BRONZE;

			const double tau_bar = polynomial_5(cos_phi, coeffs_tau);
			const double rho_bar_term = polynomial_5(cos_phi, coeffs_rho);
			const double rho_bar = rho_bar_term - tau_bar;

			const int flat_idx = itheta * nwvl + iwvl;
			layer->t_lambda_theta[flat_idx] = t_0_lambda * tau_bar;
			layer->rf_lambda_theta[flat_idx] = layer->rf_0[iwvl] * (1.0 - rho_bar) + rho_bar;
			layer->rb_lambda_theta[flat_idx] = layer->rb_0[iwvl] * (1.0 - rho_bar) + rho_bar;
		}
	}
}


void multi_layer_calc(
	GlazingLayer *layers,
	int nlayers,
	int nwvl,
	double *total_t,
	double *total_rf,
	double *total_rb)
{
	if (nlayers <= 0) 
		return;

	size_t total_size = (size_t)NTHETA * nwvl * sizeof(double);

	memcpy(total_t, layers[0].t_lambda_theta, total_size);
	memcpy(total_rf, layers[0].rf_lambda_theta, total_size);
	memcpy(total_rb, layers[0].rb_lambda_theta, total_size);

	if (nlayers == 1) 
		return;

	double *prev_t = malloc(total_size);
	double *prev_rf = malloc(total_size);
	double *prev_rb = malloc(total_size);

	if (!prev_t || !prev_rf || !prev_rb)
		perror("Failed to allocate temporary storage in multi_layer");

	for (int j = 1; j < nlayers; ++j) {
		memcpy(prev_t, total_t, total_size);
		memcpy(prev_rf, total_rf, total_size);
		memcpy(prev_rb, total_rb, total_size);

		double *t_j = layers[j].t_lambda_theta;
		double *rf_j = layers[j].rf_lambda_theta;
		double *rb_j = layers[j].rb_lambda_theta;

		for (int itheta = 0; itheta < NTHETA; ++itheta) {
			for (int iwvl = 0; iwvl < nwvl; ++iwvl) {
				const int idx = itheta * nwvl + iwvl;
				double denominator = 1.0 - rf_j[idx] * prev_rb[idx];
				if (fabs(denominator) < DBL_EPSILON) {
					denominator = DBL_EPSILON;
				}
				total_t[idx] = prev_t[idx] * t_j[idx] / denominator;
				total_rf[idx] = prev_rf[idx] + prev_t[idx] * prev_t[idx] * rf_j[idx] / denominator;
				total_rb[idx] = rb_j[idx] + t_j[idx] * t_j[idx] * prev_rb[idx] / denominator;
			}
		}
	}
	free(prev_t);
	free(prev_rf);
	free(prev_rb);
}


int add_layer(GlazingLayer **layers_array, int *count, int *capacity,
			  const char *filename, double thickness_m, int is_monolithic) {
	if (*count >= *capacity) {
		*capacity = (*capacity == 0) ? 4 : *capacity * 2;
		GlazingLayer *new_layers = realloc(*layers_array, *capacity * sizeof(GlazingLayer));
		if (!new_layers) {
			perror("Failed to reallocate memory for layers");
			return 0;
		}
		*layers_array = new_layers;
	}

	memset(&((*layers_array)[*count]), 0, sizeof(GlazingLayer));

	#ifdef _WIN32
		(*layers_array)[*count].filename = _strdup(filename);
	#else
		(*layers_array)[*count].filename = strdup(filename);
	#endif

	if (!(*layers_array)[*count].filename) {
		 perror("Failed to duplicate filename");
		 return 0;
	}
	(*layers_array)[*count].is_mono = is_monolithic;
	(*layers_array)[*count].thickness_m = thickness_m;
	(*layers_array)[*count].t_0 = NULL;
	(*layers_array)[*count].rf_0 = NULL;
	(*layers_array)[*count].rb_0 = NULL;
	(*layers_array)[*count].t_lambda_theta = NULL;
	(*layers_array)[*count].rf_lambda_theta = NULL;
	(*layers_array)[*count].rb_lambda_theta = NULL;

	(*count)++;
	return 1;
}


int interpolate(GlazingLayer *layer, const double wvl_start, const double wvl_end, const double wvl_interval, const int nwvl) {
	DATARRAY *dp = getdata(layer->filename);
	if (!dp) {
		fprintf(stderr, "Error: Cannot open file '%s'\n", layer->filename);
		return 0;
	}

	layer->t_0 = malloc(nwvl * sizeof(double));
	layer->rf_0 = malloc(nwvl * sizeof(double));
	layer->rb_0 = malloc(nwvl * sizeof(double));

	double wvl = wvl_start;
	int i = 0;
	while (wvl <= wvl_end) {
		double t_pt[2] = {2., wvl};
		layer->t_0[i] = datavalue(dp, t_pt);
		double rf_pt[2] = {0., wvl};
		layer->rf_0[i] = datavalue(dp, rf_pt);
		double rb_pt[2] = {1., wvl};
		layer->rb_0[i] = datavalue(dp, rb_pt);
		wvl = wvl + wvl_interval;
		i = i + 1;

	}
	freedata(dp);
	return 1;
}


int write_output_file(
		const char *tfilename, 
		const char *rfilename, 
		const double *tdata, 
		const double *rfdata, 
		const double *rbdata, 
		const double wvl_start, 
		const double wvl_end, 
		const int nwvl,
		int argc, char *argv[])
{
	FILE *tfp = fopen(tfilename, "w");
	FILE *rfp = fopen(rfilename, "w");
	if (!tfp || !rfp) {
		fprintf(stderr, "Error: Cannot open output files\n");
		return(0);
	}

	fprintf(tfp, "# ");
	for (int i = 0; i< argc; ++i) {
		fprintf(tfp, "%s ", argv[i]);
	}
	fprintf(tfp, "\n2\n0 0 %d", NTHETA);
	for (int i = 0; i < NTHETA; i++) {
		fprintf(tfp, " %d", (int)THETAS[i]);
	}
	fprintf(tfp, "\n%.0f %.0f %d\n", wvl_start, wvl_end, nwvl);

	for (int itheta = 0; itheta < NTHETA; ++itheta) {
		for (int iwvl = 0; iwvl < nwvl; ++iwvl) {
			fprintf(tfp, "%.6f\n", tdata[itheta*nwvl+iwvl]);
		}
	}
	fclose(tfp);

	fprintf(rfp, "# ");
	for (int i = 0; i< argc; ++i) {
		fprintf(rfp, "%s ", argv[i]);
	}
	fprintf(rfp, "\n2\n0 0 %d", NTHETA * 2 - 1);
	for (int i = 0; i < NTHETA; i++) {
		fprintf(rfp, " %d", (int)THETAS[i]);
	}
	for (int i = NTHETA-2; i >= 0; --i) {
		fprintf(rfp, " %d", (int)(180 - THETAS[i]));
	}
	fprintf(rfp, "\n%.0f %.0f %d\n", wvl_start, wvl_end, nwvl);

	for (int itheta = 0; itheta < NTHETA; ++itheta) {
		for (int iwvl = 0; iwvl < nwvl; ++iwvl) {
			fprintf(rfp, "%.6f\n", rfdata[itheta * nwvl + iwvl]);
		}
	}
	for (int itheta = NTHETA - 2; itheta >= 0 ; --itheta) {
		for (int iwvl = 0; iwvl < nwvl; ++iwvl) {
			fprintf(rfp, "%.6f\n", rbdata[itheta * nwvl + iwvl]);
		}
	}

	fclose(rfp);
	return(1);
}


void print_usage() {
	fprintf(stderr, "Usage: %s [-m monolithic_layer.dat thickness | -c coated_layer.dat] ...\n", progname);
	fprintf(stderr, "Calculate multi-layer glazing optics from spectral data files.\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -m <filename> thickness    Specify an uncoated (monolithic) glazing layer .dat file and its thickness (meter).\n");
	fprintf(stderr, "  -c <filename>    Specify a coated or laminate glazing layer .dat file.\n");
	fprintf(stderr, "  -s start_wvl end_wvl interval    Specify wavelength range and interval.\n");
	fprintf(stderr, "  -p prefix Specify prefix name to the output files.\n");
	fprintf(stderr, "  -h, --help        Show this help message.\n");
	fprintf(stderr, "Layer order is determined by the sequence of options on the command line.\n");
	fprintf(stderr, "Output Files: prefix_t.dat, prefix_r.dat\n");
}


static void free_layer_data(GlazingLayer *layer) {
    if (!layer) return;
    
    free(layer->filename);
    free(layer->t_0);
    free(layer->rf_0);
    free(layer->rb_0);
    free(layer->t_lambda_theta);
    free(layer->rf_lambda_theta);
    free(layer->rb_lambda_theta);
    
    // Set pointers to NULL to prevent double-free errors
    layer->filename = NULL;
    layer->t_0 = NULL;
    layer->rf_0 = NULL;
    layer->rb_0 = NULL;
    layer->t_lambda_theta = NULL;
    layer->rf_lambda_theta = NULL;
    layer->rb_lambda_theta = NULL;
}


void cleanup_layers(GlazingLayer *layers, int num_layers) {
	 if (!layers) return;
	 for (int i = 0; i < num_layers; ++i) {
		if (layers[i].filename)
			free_layer_data(&layers[i]);
	}
	free(layers);
}


int main(int argc, char *argv[])
{
	GlazingLayer *layers = NULL;
	int num_layers = 0;
	int layer_capacity = 0;
	int success = 1;
	double wvl_start_nm = 380.;
	double wvl_end_nm = 780.;
	double wvl_interval_nm = 5.;
	int nwvl = 81;
	double thickness_m = 0.003;
	char *filename;
	char *prefix = "unnamed";
	char file_t[MAX_NAME];
	char file_r[MAX_NAME];
	progname = argv[0];
	if (argc <= 1) {
	  print_usage();
	  return EXIT_FAILURE;
	}

	for (int i=1; i < argc; ++i) {
		switch (argv[i][1]) {
		case 'm':
			filename = argv[++i];
			thickness_m = atof(argv[++i]);
			if (!add_layer(&layers, &num_layers, &layer_capacity, filename, thickness_m, 1)) {
				cleanup_layers(layers, num_layers);
				return EXIT_FAILURE;
			}
			break;
		case 'c':
			filename = argv[++i];
			if (!add_layer(&layers, &num_layers, &layer_capacity, filename, thickness_m, 0)) {
				 cleanup_layers(layers, num_layers);
				 return EXIT_FAILURE;
			}
			break;
		case 'p':
			prefix = argv[++i];
			break;
		case 's':
			wvl_start_nm = atof(argv[++i]);
			wvl_end_nm = atof(argv[++i]);
			wvl_interval_nm = atof(argv[++i]);
			if (wvl_start_nm > wvl_end_nm) {
				fprintf(stderr, "Starting wavelength > End wavelength\n");
				cleanup_layers(layers, num_layers);
				return EXIT_FAILURE;
			}
			if (((int)(wvl_end_nm - wvl_start_nm) % (int)wvl_interval_nm) > 0) {
				fprintf(stderr, 
						"Error: Wavelength range (%f to %f nm) must be evenly divisible by the interval (%f nm).\n", 
						wvl_start_nm, wvl_end_nm, wvl_interval_nm);
				cleanup_layers(layers, num_layers);
				return EXIT_FAILURE;
			}
			break;
		case 'h':
			print_usage();
			return 1;
		default:
			break;
		}
	}
	snprintf(file_t, MAX_NAME, "%s_t.dat", prefix);
	snprintf(file_r, MAX_NAME, "%s_r.dat", prefix);

	if (fabs(wvl_start_nm - wvl_end_nm) < DBL_EPSILON && wvl_interval_nm > DBL_EPSILON) {
		 nwvl = 1;
	} else {
		 nwvl = (int)((wvl_end_nm - wvl_start_nm) / wvl_interval_nm + 1.5);
	}

	if (num_layers == 0) {
		fprintf(stderr, "Error: No layers specified.\n");
		print_usage();
		return EXIT_FAILURE;
	}

	for (int i = 0; i < num_layers; ++i) {
		if (!interpolate(&layers[i], wvl_start_nm, wvl_end_nm, wvl_interval_nm, nwvl)) {
			fprintf(stderr, "Error: Failed to parse or interpolate data for layer %d.\n", i + 1);
			success = 0;
			break;
		}
		layers[i].t_lambda_theta = malloc(NTHETA * nwvl * sizeof(double));
		layers[i].rf_lambda_theta = malloc(NTHETA * nwvl * sizeof(double));
		layers[i].rb_lambda_theta = malloc(NTHETA * nwvl * sizeof(double));

		if (layers[i].is_mono) {
			angular_monolithic(&layers[i], nwvl);
		} else {
			angular_coated(&layers[i], nwvl);
		}
	}

	if (!success) {
		cleanup_layers(layers, num_layers);
		return EXIT_FAILURE;
	}


	size_t total_flat_size = (size_t)NTHETA * nwvl * sizeof(double);
	double *total_t = malloc(total_flat_size);
	double *total_rf = malloc(total_flat_size);
	double *total_rb = malloc(total_flat_size);

	if (!total_t || !total_rf || !total_rb) {
		perror("Failed to allocate memory for final results");
		cleanup_layers(layers, num_layers);
		free(total_t);
		free(total_rf);
		free(total_rb);
		return EXIT_FAILURE;
	}

	multi_layer_calc(layers, num_layers, nwvl, total_t, total_rf, total_rb);

	success &= write_output_file(file_t, file_r, total_t, total_rf, total_rb, wvl_start_nm, wvl_end_nm, nwvl, argc, argv);

	if (success) {
		printf("# ");
		for (int i = 0;i < argc; ++i) {
			printf("%s ", argv[i]);
		}
		printf("\nvoid specdata refl_spec_%s\n", prefix);
		printf("4 noop %s . 'Acos(Rdot)/DEGREE'\n0\n0\n\n", file_r);
		printf("void specdata trans_spec_%s\n", prefix);
		printf("4 noop %s . 'Acos(abs(Rdot))/DEGREE'\n0\n0\n\n", file_t);
		printf("void WGMDfunc glaze_mat_%s\n13\n\trefl_spec_%s 1 0 0\n\ttrans_spec_%s 1 0 0\n\tvoid\n\t0 0 1 .\n0\n9 0 0 0 0 0 0 0 0 0\n", prefix, prefix, prefix);
	} else {
		fprintf(stderr, "Error: Failed to write one or more output files.\n");
	}

	cleanup_layers(layers, num_layers);
	free(total_t);
	free(total_rf);
	free(total_rb);

	return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
