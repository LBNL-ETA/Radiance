#ifndef lint
static const char RCSid[] = "$Id$";
#endif
/*
 * Program to convert Alstan's spectral database to Radiance materials.
 *
 * Each material is output to a separate file, and materials are sorted
 * into separate subfolders of the current directory based on category.
 *
 *  G. Ward	Feb. 2026
 */

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>
#include <ctype.h>
#include "readCSV.h"

#define SPECMULT	0.01

/* Spectral data holder */
typedef struct {
	double	lambda_min, lambda_max;
	int	nmeas;
	float	curv[1];	/* extends struct */
} SPECTRUM;

#define alloc_spectrum(n)	(SPECTRUM *)malloc(sizeof(SPECTRUM) + \
					sizeof(float)*((n)-1))

/* 1 means force overwrite, 0 means append file, -1 means silently skip */
static int	overwrite = 0;

/* Our list of created subfolders */
static CSVREC	*subfldr = NULL;

static const char *
get_subfolder(const char *dir_name)
{
	int	i;

	if (!dir_name || !*dir_name | !strcmp(dir_name, "."))
		return ".";

	if (!subfldr)
		subfldr = alloc_csvrec(16);

	for (i = 0; i < subfldr->nf && subfldr->f[i]; i++)
		if (!strcasecmp(subfldr->f[i], dir_name))
			return subfldr->f[i];

	if (i >= subfldr->nf)
		subfldr = realloc_csvrec(subfldr, i + i/2 + 16);

	if ((mkdir(dir_name, 0777) < 0 && errno != EEXIST) ||
			!set_csvfield(subfldr, i, dir_name)) {
		perror(dir_name);
		exit(1);
	}
	return subfldr->f[i];
}

/* Lookup for fields we may need */
static struct {
	short	Name, ObjectType, RadianceType, MeasurementCredit,
		Comments, Specularity, Roughness, MeasurementType,
		SCIMeasures, SCEMeasures;
}	fndx;

static int
find_field(const CSVREC *rp, const char *fn)
{
	int	i = rp->nf;

	while (i--)
		if (!strcasecmp(rp->f[i], fn))
			return i;

	fprintf(stderr, "Cannot find requied field '%s' in header\n", fn);
	exit(1);
}

static void
index_fields(const CSVREC *rp)
{
	fndx.Name = find_field(rp, "Name");
	fndx.ObjectType = find_field(rp, "ObjectType");
	fndx.RadianceType = find_field(rp, "RadianceType");
	fndx.MeasurementCredit = find_field(rp, "MeasurementCredit");
	fndx.Comments = find_field(rp, "Comments");
	fndx.Specularity = find_field(rp, "Specularity");
	fndx.Roughness = find_field(rp, "Roughness");
	fndx.MeasurementType = find_field(rp, "MeasurementType");
	fndx.SCIMeasures = find_field(rp, "SCIMeasures");
	fndx.SCEMeasures = find_field(rp, "SCEMeasures");
}

static int
mk_name(char *nm, const char *orig)
{
	if (!orig || !*orig) {
		fputs("Missing material name!\n", stderr);
		return 0;
	}
	while (*orig) {
		if ((*orig > ' ') & (*orig != '/'))
			*nm++ = tolower(*orig);
		else if (orig[1] > ' ')
			*nm++ = '_';
		++orig;
	}
	*nm = '\0';
	return 1;
}

static SPECTRUM *
get_spectrum(const char *sp)
{
	SPECTRUM	shead, *spec;
	const char	*cp;
	double		wl;
	int		i;

	if (!sp || !*sp)
		return NULL;
	shead.lambda_min = 100000;
	shead.lambda_max = 0;
	shead.nmeas = 0;
	cp = sp;
	if (*cp++ != '{')
		goto badspec;
	do {
		while (isspace(*cp))
			cp++;
		wl = atof(cp);
		if (wl < 100)
			goto badspec;
		if (wl > 10000)
			goto badspec;
		if (wl < shead.lambda_min)
			shead.lambda_min = wl;
		if (wl > shead.lambda_max)
			shead.lambda_max = wl;
		shead.nmeas++;
		while (*cp != ':')
			if (!*cp++)
				goto badspec;
		while (isspace(*++cp))
			;
		while (strchr("9876543210.+-eE", *cp))
			if (!*cp++)
				goto badspec;
		while (isspace(*cp))
			cp++;
		if ((*cp != ',') & (*cp != '}'))
			goto badspec;
	} while (*cp++ != '}');

	if ((*cp != '\0') | (shead.lambda_min >= shead.lambda_max)
			| (shead.nmeas < 3))
		goto badspec;

	spec = alloc_spectrum(shead.nmeas);
	if (!spec) {
		perror("malloc");
		exit(1);
	}
	*spec = shead;		/* go back for values */
	cp = sp+1;
	do {
		while (isspace(*cp))
			cp++;
		wl = atof(cp);
		i = (wl - shead.lambda_min) /
				(shead.lambda_max - shead.lambda_min)
				* (shead.nmeas - 1) + .5;
		while (*cp != ':')
			cp++;
		while (isspace(*++cp))
			;
		spec->curv[i] = SPECMULT * atof(cp);
		while ((*cp != ',') & (*cp != '}'))
			cp++;
	} while (*cp++ != '}');

	return spec;
badspec:
	fprintf(stderr, "Bad spectrum: %s\n", sp);
	return NULL;
}

#define	free_spectrum	free

static double
max_spectrum(const SPECTRUM *spec, double *wlp)
{
	int	maxi = 0;
	int	i = spec->nmeas;

	while (--i > 0)
		if (spec->curv[i] > spec->curv[maxi])
			maxi = i;
	if (wlp)
		*wlp = (double)maxi / (spec->nmeas - 1.0) *
			 (spec->lambda_max - spec->lambda_min)
			 + spec->lambda_min;

	return spec->curv[maxi];
}

static double
avg_spectrum(const SPECTRUM *spec)
{
	double	sum = 0;
	int	i = spec->nmeas;

	while (i--)
		sum += spec->curv[i];

	return sum / (double)spec->nmeas;
}

enum {sBad, sLambertian, sWhite, sMetallic, sOther};

#define MAX2(a,b)	((a)>(b) ? (a) : (b))

static int
classify_reflection(const SPECTRUM *sci, const SPECTRUM *sce)
{
	const double	nf = 1./(double)sci->nmeas;
	double		dsum=0, d2sum=0, rsum=0, r2sum=0, scisum=0;
	int		i = sci->nmeas;

	if (sci == sce)
		return sLambertian;	/* only option! */

	while (i--) {
		double	d = sci->curv[i] - sce->curv[i];
		dsum += d;
		d2sum += d*d;
		d /= MAX2(sci->curv[i], 1e-5);
		rsum += d;
		r2sum += d*d;
		scisum += sci->curv[i];
	}
	if (dsum*nf <= -0.010)
		return sBad;		/* sci <==> sce ??? */
	if (dsum >= 0.95*scisum)
		return sMetallic;	/* treat as pure specular */
	if ((rsum*nf <= 0.010) | (dsum*nf <= 0.003))
		return sLambertian;	/* treat as pure diffuse */
	if (dsum*nf <= 0.020 || (dsum*nf <= 0.075) &
				(d2sum/(dsum*dsum*nf) - 1. < 0.15*0.15))
		return sWhite;		/* use white highlight */
	if (r2sum/(rsum*rsum*nf) - 1. < 0.15*0.15)
		return sMetallic;	/* highlight matches color */

	return sOther;			/* otherwise need both spectra */
}

static void
put_spectrum(const SPECTRUM *spec, double mf, const char *mname, FILE *fp)
{
	const int	vals_row = 6;
	int		i;

	fprintf(fp, "void spectrum %s\n0\n0\n%d  %.0f %.0f", mname,
			spec->nmeas+2, spec->lambda_min, spec->lambda_max);
	for (i = 0; i < spec->nmeas; i++) {
		if (!(i % vals_row)) fputs("\n       ", fp);
		fprintf(fp, " %.4f", mf*spec->curv[i]);
	}
	fputs("\n\n", fp);
}

static void
create_radspec(const CSVREC *rp)
{
	static const char	diffuse_spec_name[] = "_diffuseSpectrum*";
	static const char	specular_spec_name[] = "_specularSpectrum*";
	char			lc_name[80];
	char			outfile[256];
	const char		*folder;
	FILE			*fout;
	const char		*fs;
	SPECTRUM		*sci, *sce;
	double			sci_max, sce_max;
	double			sci_avg, sce_avg;
	double			rough;
	double			d;
	int			i;

	if (!mk_name(lc_name, get_csvfield(rp, fndx.Name)))
		return;
	folder = get_subfolder(get_csvfield(rp, fndx.ObjectType));
	sprintf(outfile, "%s/%s.mat", folder, lc_name);
	errno = 0;
	if (overwrite <= 0 && access(outfile, F_OK) == 0) {
		if (overwrite < 0)
			return;		/* silently skip */
		fprintf(stderr, "%s: Warning - appending existing file\n",
				outfile);
	}
	if ((errno != 0) & (errno != ENOENT))
		goto syserror;

	sci = get_spectrum(get_csvfield(rp, fndx.SCIMeasures));
	sce = get_spectrum(get_csvfield(rp, fndx.SCEMeasures));
	if (!sci) sci = sce;
	else if (!sce) sce = sci;
	if (!sci) {
		fprintf(stderr, "%s: Need SCI/SCE Measures - skipping...\n",
				rp->f[fndx.Name]);
		return;
	}
	if ((sci->lambda_min != sce->lambda_min) |
			(sci->lambda_max != sce->lambda_max) |
			(sci->nmeas != sce->nmeas)) {
		fprintf(stderr, "%s: Mismatched measures for SCI/SCE - skipping...\n",
				rp->f[fndx.Name]);
		return;
	}
	fout = fopen(outfile, overwrite ? "w" : "a");
	if (!fout)
		goto syserror;
	fputs("# Spectral material translated by csv2radspec\n", fout);
	fprintf(fout, "# Material name: %s\n", rp->f[fndx.Name]);
	if ((fs = get_csvfield(rp, fndx.ObjectType)) && *fs)
		fprintf(fout, "# Object type: %s\n", fs);
	if ((fs = get_csvfield(rp, fndx.MeasurementType)) && *fs)
		fprintf(fout, "# Measurement type: %s\n", fs);
	if ((fs = get_csvfield(rp, fndx.MeasurementCredit)) && *fs)
		fprintf(fout, "# Measured by: %s\n", fs);
	if ((fs = get_csvfield(rp, fndx.Comments)) && *fs)
		fprintf(fout, "# %s\n", fs);
	fputc('\n', fout);
	if ((fs = get_csvfield(rp, fndx.Roughness)) && *fs)
		rough = atof(fs);
	else
		rough = 0;
	sci_max = max_spectrum(sci, NULL);
	sce_max = max_spectrum(sce, NULL);
	sci_avg = avg_spectrum(sci);
	sce_avg = avg_spectrum(sce);
#if 0
	printf("%s: SCI (max,avg)=(%.3f,%.3f), SCE (max,avg)=(%.3f,%.3f)\n",
		rp->f[fndx.Name], sci_max, sci_avg, sce_max, sce_avg);
#endif
	if ((sci_max*1.01 <= sce_max) | (sci_avg <= sce_avg))
		fprintf(stderr, "%s: Warning - SCE is greater than SCI\n",
				rp->f[fndx.Name]);

	switch (classify_reflection(sci, sce)) {
	case sLambertian:		/* pure diffuse */
		put_spectrum(sci, 1./sci_max, diffuse_spec_name, fout);
		fprintf(fout, "%s plastic %s\n0\n0\n5 ",
				diffuse_spec_name, lc_name);
		fprintf(fout, " %.4f %.4f %.4f  0  %.4f\n\n",
				sci_max, sci_max, sci_max, rough);
		break;
	case sWhite:			/* plastic (white highlight) */
		put_spectrum(sce, 1./sce_max, diffuse_spec_name, fout);
		fprintf(fout, "%s plastic %s\n0\n0\n5 ",
				diffuse_spec_name, lc_name);
		d = sce_max/(1. + sce_avg - sci_avg);
		fprintf(fout, " %.4f %.4f %.4f  %.3f  %.4f\n\n",
				d, d, d, sci_avg - sce_avg, rough);
		break;
	case sMetallic:			/* metal */
		put_spectrum(sci, 1./sci_max, specular_spec_name, fout);
		fprintf(fout, "%s metal %s\n0\n0\n5 ",
				specular_spec_name, lc_name);
		fprintf(fout, " %.4f %.4f %.4f  %.3f  %.4f\n\n",
				sci_max, sci_max, sci_max,
				1. - sce_avg/sci_avg, rough);
		break;
	case sOther:			/* WGMDfunc */
		i = sci->nmeas;		/* update sci = sci - sce */
		sci_max = 0;
		while (i--)
			if ((sci->curv[i] -= sce->curv[i]) > sci_max)
				sci_max = sci->curv[i];
		sci_avg -= sce_avg;
		put_spectrum(sce, 1./sce_max, diffuse_spec_name, fout);
		put_spectrum(sci, 1./sci_max, specular_spec_name, fout);
		fprintf(fout, "%s WGMDfunc %s\n13\n",
				diffuse_spec_name, lc_name);
		fprintf(fout, "\t%s  %.3f  %.4f %.4f\n",
				specular_spec_name, sci_max, rough, rough);
		fprintf(fout, "\tvoid  0  0 0\n\tvoid\n\t1 1 1  .\n0\n9");
		fprintf(fout, "\t%.4f %.4f %.4f\n", sce_max, sce_max, sce_max);
		fprintf(fout, "\t%.4f %.4f %.4f\n", sce_max, sce_max, sce_max);
		fprintf(fout, "\t0 0 0\n\n");
		break;
	default:
		fputs("# Bad SCE/SCI combo, aborting this material\n\n", fout);
		fprintf(stderr, "%s: Bad SCE/SCI - skipping...\n",
				rp->f[fndx.Name]);
		break;
	}
	free_spectrum(sci);
	if (sce != sci) free_spectrum(sce);
	if (fclose(fout) == 0)
		return;
syserror:
	perror(outfile);
}

int
main(int argc, char *argv[])
{
	char	*inpname = "<stdin>";
	CSVREC	*rcp;

	if (argc > 1 && argv[1][0] == '-') {
		if (argv[1][1] == 'f')
			overwrite = 1;
		else if (argv[1][1] == 'n')
			overwrite = -1;
		else
			goto userror;
	}
	if (argc > (overwrite!=0)+2)
		goto userror;
	if (argc == (overwrite!=0)+2 &&
			!freopen(inpname=argv[(overwrite!=0)+1], "r", stdin))
		goto inperror;
					/* get header and index */
	rcp = read_csvrec(stdin, NULL);
	if (!rcp)
		goto inperror;
	index_fields(rcp);
	free_csv(rcp);
					/* read each record */
	while ((rcp = read_csvrec(stdin, NULL))) {
		create_radspec(rcp);
		free_csv(rcp);
	}
	if (feof(stdin))
		return 0;
inperror:
	perror(inpname);
	return 1;
userror:
	fprintf(stderr, "Usage: %s [-f|-n][input.csv]\n", argv[0]);
	return 1;
}
