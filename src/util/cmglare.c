/*
 *  cmglare.c - routines for calculating glare autonomy.
 *
 *  N. Jones
 */

/*
 * Copyright (c) 2017-2019 Nathaniel Jones
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "rtmath.h"
#include "cmglare.h"

#define LUMINOUS_EFFICACY	179	/* lumens per Watt */
#define LUMINANCE_THRESHOLD	100 /* minimum threshold that will be interpreted as luminance rather than ratio between Ev and glare source */

#define ANGLE(u,v)	acos(DOT((u),(v)))	// TODO need to normalize?

typedef struct reinhart_sky {
	//int mf;					/* Linear divisions per Tregenza patch. */
	int rings;					/* Number of rings of sky patches. */
	//int patches;				/* Number of sky patches. */
	double ringElevationAngle;	/* Angle in radians between rings. */
	int *patchesPerRow;			/* Sky patches per row. */
	int *firstPatchIndex;		/* Index of first sky patch in each row. */
	double *solidAngle;			/* Solid angle of each patch in each row. */
} ReinhartSky;

static const int tnaz[] = { 30, 30, 24, 24, 18, 12, 6 };	/* Number of patches per row */

extern char* progname;

static ReinhartSky* make_sky(const CMATRIX *smx)
{
	int mf, count, i, j;
	double remaining = PI / 2;

	/* Check sky partitionss */
	switch (smx->nrows) { // nrows is Reinhart sky subdivisoins plus one for miss
	case 2: mf = 0; break;
	case 146: mf = 1; break;
	case 578: mf = 2; break;
	case 1298: mf = 3; break;
	case 2306: mf = 4; break;
	case 3602: mf = 5; break;
	case 5186: mf = 6; break;
	case 7058: mf = 7; break;
	case 9218: mf = 8; break;
	case 11666: mf = 9; break;
	case 14402: mf = 10; break;
	case 17426: mf = 11; break;
	case 20738: mf = 12; break;
	default:
		fprintf(stderr,
			"%s: unknown number of sky patches %d\n",
			progname, smx->nrows);
		return NULL;
	}

	/* Allocate sky */
	ReinhartSky *sky = (ReinhartSky*)malloc(sizeof(ReinhartSky));
	if (!sky) goto memerr;

	/* Calculate patches per row */
	sky->rings = 7 * mf + 1;
	sky->patchesPerRow = (int*)malloc(sky->rings * sizeof(int));
	sky->firstPatchIndex = (int*)malloc(sky->rings * sizeof(int));
	if (!sky->patchesPerRow || !sky->firstPatchIndex) goto memerr;
	count = 1; // The below horizon patch
	for (i = 0; i < 7; i++)
		for (j = 0; j < mf; j++) {
			sky->firstPatchIndex[i * mf + j] = count;
			count += sky->patchesPerRow[i * mf + j] = tnaz[i] * mf;
		}
	sky->firstPatchIndex[7 * mf] = count;
	sky->patchesPerRow[7 * mf] = 1;
	//sky->patches = count + 1;

	/* Calculate solid angle of patches in each row */
	sky->solidAngle = (double*)malloc(sky->rings * sizeof(double));
	if (!sky->solidAngle) goto memerr;
	sky->ringElevationAngle = PI / (2 * sky->rings - 1);
	sky->solidAngle[0] = 2 * PI;
	for (i = 1; i < sky->rings; i++) {
		remaining -= sky->ringElevationAngle;
		sky->solidAngle[i] = 2 * PI * (1 - cos(remaining)); // solid angle of cap
		sky->solidAngle[i - 1] -= sky->solidAngle[i];
		sky->solidAngle[i - 1] /= sky->patchesPerRow[i - 1];
	}

	return sky;

memerr:
	fprintf(stderr,
		"%s: out of memory\n",
		progname);
	return NULL;
}

static void free_sky(ReinhartSky *sky)
{
	if (sky) {
		free(sky->patchesPerRow);
		free(sky->firstPatchIndex);
		free(sky->solidAngle);
		free(sky);
	}
}

static void get_patch_direction(const ReinhartSky *sky, const int patch, FVECT target)
{
	//if (patch >= sky->patches || patch < 0) throw new RuntimeException("Illegal patch " + patch);
	if (!patch) {
		target[0] = target[1] = 0;
		target[2] = -1;
		return; // Ignore below horizon?
	}
	int row = sky->rings - 1;
	while (patch < sky->firstPatchIndex[row]) row--;
	const double alt = PI / 2 - sky->ringElevationAngle * (row - (sky->rings - 1));
	const double azi = 2 * PI * (patch - sky->firstPatchIndex[row]) / sky->patchesPerRow[row];
	const double cos_alt = cos(alt);
	target[0] = cos_alt * -sin(azi);
	target[1] = cos_alt * -cos(azi);
	target[2] = sin(alt);
}

static double get_patch_solid_angle(const ReinhartSky *sky, const int patch, const double cos_theta)
{
	//if (patch >= sky->patches || patch < 0) throw new RuntimeException("Illegal patch " + patch);
	if (!patch) return 2 * (PI - 2 * acos(cos_theta)); // Solid angle overlap between visible hemisphere and ground hemisphere
	int row = sky->rings - 1;
	while (patch < sky->firstPatchIndex[row]) row--;
	return sky->solidAngle[row];
}

static double get_guth(const FVECT dir, const FVECT forward, const FVECT up)
{
	double posindex;
	FVECT hv, temp;
	VCROSS(hv, forward, up);
	normalize(hv);
	VCROSS(temp, forward, hv);
	double phi = ANGLE(dir, temp) - PI / 2;

	/* Guth model, equation from IES lighting handbook */
	if (phi >= 0) {
		double sigma = ANGLE(dir, forward);
		VSUM(hv, forward, dir, 1 / DOT(dir, forward));
		normalize(hv);
		double tau = ANGLE(hv, up);
		tau *= 180.0 / PI;
		sigma *= 180.0 / PI;

		if (sigma <= 0)
			sigma = -sigma;

		posindex = exp((35.2 - 0.31889 * tau - 1.22 * exp(-2.0 * tau / 9.0)) / 1000.0 * sigma + (21.0 + 0.26667 * tau - 0.002963 * tau * tau) / 100000.0 * sigma * sigma);
	}
	/* below line of sight, using Iwata model */
	else {
		double teta = PI / 2 - ANGLE(dir, hv);

		if (teta == 0.0)
			teta = FTINY;

		double fact = 0.8;
		double d = 1 / tan(phi);
		double s = tan(teta) / tan(phi);
		double r = sqrt((1 + s * s) / (d * d));
		if (r > 0.6)
			fact = 1.2;
		if (r > 3)
			r = 3.0;

		posindex = 1.0 + fact * r;
	}
	if (posindex > 16)
		posindex = 16.0;

	return posindex;
}

float* cm_glare(const CMATRIX *dcmx, const CMATRIX *evmx, const CMATRIX *smx, const int *occupied, const double dgp_limit, const double dgp_threshold, const FVECT *views, const FVECT dir, const FVECT up)
{
	int p, t, c;
	int hourly_output = dgp_limit < 0;
	float *dgp_list;
	ReinhartSky *sky;
	FVECT vdir;

	/* Check consistancy */
	if ((dcmx->nrows != evmx->nrows) | (dcmx->ncols != smx->nrows) | (evmx->ncols != smx->ncols)) {
		fprintf(stderr,
			"%s: inconsistant matrix dimensions: dc(%d, %d) ev(%d, %d) s(%d, %d)\n",
			progname, dcmx->nrows, dcmx->ncols, evmx->nrows, evmx->ncols, smx->nrows, smx->ncols);
		return NULL;
	}

	/* Create output buffer */
	dgp_list = (float*)malloc(evmx->nrows * (hourly_output ? evmx->ncols : 1) * sizeof(float));
	if (!dgp_list) {
		fprintf(stderr,
			"%s: out of memory in cm_glare()\n",
			progname);
		return NULL;
	}

	/* Create sky */
	sky = make_sky(smx);
	if (!sky) return NULL;

	/* Calculate glare limit */
	double ev_max = -1;
	if (!hourly_output) {
		ev_max = (dgp_limit - 0.16) / 5.87e-5;
		if (ev_max < 0) ev_max = 0;
	}

	/* For each position and direction */
	if (!views) VCOPY(vdir, dir);
	for (p = 0; p < evmx->nrows; p++) {
		/* For each time step */
		int occupied_hours = 0;
		int glare_hours = 0;
		if (views) VCOPY(vdir, views[p]);
		for (t = 0; t < evmx->ncols; t++) {
			if (!occupied[t]) {
				/* Not occupied */
				if (hourly_output) dgp_list[p * evmx->ncols + t] = 0.0f;
			}
			else {
				/* Occupied */
				double illum = LUMINOUS_EFFICACY * bright(cm_lval(evmx, p, t));
				occupied_hours++;

				if (illum <= FTINY) {
					/* No light, so no glare */
					if (hourly_output) dgp_list[p * evmx->ncols + t] = 0.0f;
				}
				else if ((illum >= ev_max) & (!hourly_output)) {
					/* Guarangeed glare */
					glare_hours++;
				}
				else {
					/* Calculate enhanced simplified daylight glare probability */
					double sum = 0.0;
					FVECT patch_normal;
					for (c = 0; c < smx->nrows; c++) {
						const double dc = bright(cm_lval(dcmx, p, c));
						if (dc > 0) {
							get_patch_direction(sky, c, patch_normal);
							if (!c) {
								/* Direction toward the center of the visible ground */
								VADD(patch_normal, patch_normal, vdir);
								if (normalize(patch_normal) == 0) continue;
							}
							const double cos_theta = DOT(vdir, patch_normal);
							if (cos_theta <= FTINY) continue;
							const double s = LUMINOUS_EFFICACY * bright(cm_lval(smx, c, t));
							const double omega = get_patch_solid_angle(sky, c, cos_theta);
							const double patch_luminance = (dc * s) / (omega * cos_theta);

							double min_patch_luminance = dgp_threshold;
							if (dgp_threshold < LUMINANCE_THRESHOLD)
								min_patch_luminance *= illum / PI; // TODO should use average luminance, not illuminance
							if (patch_luminance < min_patch_luminance) continue;
							const double P = get_guth(patch_normal, vdir, up);
							sum += (patch_luminance * patch_luminance * omega) / (P * P);
						}
					}

					//double dgp = 5.87e-5 * illum + 0.092 * log10(1 + dgp / pow(illum, 1.87)) + 0.159;
					double eDGPs = 5.87e-5 * illum + 0.0918 * log10(1 + sum / pow(illum, 1.87)) + 0.16;
					if (illum < 1000) /* low light correction */
						eDGPs *= exp(0.024 * illum - 4) / (1 + exp(0.024 * illum - 4));
					//eDGPs /= 1.1 - 0.5 * age / 100.0; /* age correction */
					if (eDGPs > 1.0) eDGPs = 1.0;

					if (hourly_output)
						dgp_list[p * evmx->ncols + t] = (float)eDGPs;
					else if (eDGPs >= dgp_limit)
						glare_hours++;
				}
			}
		}
		if (!hourly_output) {
			/* Save glare autonomy */
			dgp_list[p] = (float)(occupied_hours - glare_hours) / occupied_hours;
		}
	}

	free_sky(sky);

	return dgp_list;
}

static int getvec(FVECT vec, const int dtype, FILE *fp)		/* get a vector from fp */
{
	static float  vf[3];
	static double  vd[3];
	char  buf[32];
	int  i;

	switch (dtype) {
	case DTascii:
		for (i = 0; i < 3; i++) {
			if (fgetword(buf, sizeof(buf), fp) == NULL ||
				!isflt(buf))
				return(-1);
			vec[i] = atof(buf);
		}
		break;
	case DTfloat:
		if (getbinary(vf, sizeof(float), 3, fp) != 3)
			return(-1);
		VCOPY(vec, vf);
		break;
	case DTdouble:
		if (getbinary(vd, sizeof(double), 3, fp) != 3)
			return(-1);
		VCOPY(vec, vd);
		break;
	default:
		fprintf(stderr,
			"%s: botched input format\n",
			progname);
		return(-1);
	}
	return(0);
}

int cm_load_schedule(const int count, int* schedule, FILE *fp)
{
	char buf[512];
	char *cp;
	char *comma;
	double val;
	int i = 0;

	while (fgetline(buf, sizeof(buf), fp) != NULL) {
		if (buf[0] == '#') continue; // Comment line
		comma = NULL;
		for (cp = buf; *cp; cp++) {
			/* If there are multiple commas, assume the value is after the last comma */
			if (*cp == ',') {
				comma = cp; /* Record position of last comma */
			}
		}
		if (comma)
			val = atof(comma + 1);
		else
			val = atof(buf);

		if (i < count) {
			/* Add the value to the schedule */
			schedule[i++] = (val > 0);
		}
		else {
			fprintf(stderr,
				"%s: too many schedule entries\n",
				progname);
			return(1);
		}
	}
	fclose(fp);

	if (i < count) {
		fprintf(stderr,
			"%s: too few schedule entries\n",
			progname);
		return(-1);
	}

	return 0;
}

FVECT* cm_load_views(const int nrows, const int dtype, FILE *fp)
{
	int i;
	double d;
	FVECT orig;

	FVECT *views = (FVECT*)malloc(nrows * sizeof(FVECT));
	if (!views) {
		fprintf(stderr,
			"%s: out of memory in cm_load_views()\n",
			progname);
		return NULL;
	}

	for (i = 0; i < nrows; i++) {
		if (getvec(orig, dtype, fp) | getvec(views[i], dtype, fp)) {
			fprintf(stderr,
				"%s: unexpected end of input, missing %d entries\n",
				progname, i);
			return NULL;
		}

		d = normalize(views[i]);
		if (d == 0.0) {				/* zero ==> flush */
			fprintf(stderr,
				"%s: zero length direction detected\n",
				progname);
			return NULL;
		}
	}

	return views;
}

int cm_write_glare(const float *mp, const int nrows, const int ncols, const int dtype, FILE *fp)
{
	static const char	tabEOL[2] = { '\t', '\n' };
	int			r, c;
	double	dc[1];

	switch (dtype) {
	case DTascii:
		for (r = 0; r < nrows; r++)
			for (c = 0; c < ncols; c++, mp++)
				fprintf(fp, "%.6e%c",
				mp[0],
				tabEOL[c >= ncols - 1]);
		break;
	case DTfloat:
		r = ncols*nrows;
		while (r > 0) {
			c = putbinary(mp, sizeof(float), r, fp);
			if (c <= 0)
				return(0);
			mp += c;
			r -= c;
		}
		break;
	case DTdouble:
		r = ncols*nrows;
		while (r--) {
			dc[0] = mp[0];
			if (putbinary(dc, sizeof(double), 1, fp) != 1)
				return(0);
			mp++;
		}
		break;
	default:
		fputs("Unsupported data type in cm_write_glare()!\n", stderr);
		return(0);
	}
	return(fflush(fp) == 0);
}