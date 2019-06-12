#ifndef lint
static const char RCSid[] = "$Id: pabopto2xyz.c,v 2.5 2019/06/12 00:09:18 greg Exp $";
#endif
/*
 * Combine PAB-Opto data files for color (CIE-XYZ) interpolation.
 *
 *	G. Ward
 */

#define _USE_MATH_DEFINES
#include "rtio.h"
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "interp2d.h"

/* #define DEBUG	1 */

#ifndef MAX_INPUTS
#define MAX_INPUTS	512		/* can handle this many files per chan */
#endif

#define ANGRES		256		/* resolution per 180 degrees */

/* Input file header information */
typedef struct {
	const char	*fname;		/* input file path */
	double		theta, phi;	/* input angles (degrees) */
	double		up_phi;		/* azimuth for "up" direction */
	int		igp[2];		/* input grid position */
	int		isDSF;		/* data is DSF (rather than BSDF)? */
	int		ndata;		/* number of data points (or 0) */
	long		dstart;		/* data start offset in file */
} PGINPUT;

/* 2-D interpolant for a particular incident direction */
typedef struct {
	INTERP2		*ip2;		/* 2-D interpolatiing function */
	float		*va;		/* corresponding value array */
} PGINTERP;

extern void		SDdisk2square(double sq[2], double diskx, double disky);

char			*progname;	/* global argv[0] */

int			incident_side = 0;	/* signum(intheta) */
int			exiting_side = 0;	/* signum(outtheta) */

int			nprocs = 1;	/* number of processes to run */

int			nchild = 0;	/* number of children (-1 in child) */

char			*basename = "pabopto_xyz";

					/* color conversion matrix */
double			sens2xyz[3][3] = {
				8.7539e-01, 3.5402e-02, 9.0381e-02,
				1.3463e+00, -3.9082e-01, 7.2721e-02,
				2.7306e-01, -2.1111e-01, 1.1014e+00,
			};

#if defined(_WIN32) || defined(_WIN64)

#define await_children(n)	(void)(n)
#define run_subprocess()	0
#define end_subprocess()	(void)0

#else	/* ! _WIN32 */

#include <unistd.h>
#include <sys/wait.h>

/* Wait for the specified number of child processes to complete */
static void
await_children(int n)
{
	int	exit_status = 0;

	if (n > nchild)
		n = nchild;
	while (n-- > 0) {
		int	status;
		if (wait(&status) < 0) {
			fprintf(stderr, "%s: missing child(ren)!\n", progname);
			nchild = 0;
			break;
		}
		--nchild;
		if (status) {			/* something wrong */
			if ((status = WEXITSTATUS(status)))
				exit_status = status;
			else
				exit_status += !exit_status;
			fprintf(stderr, "%s: subprocess died\n", progname);
			n = nchild;		/* wait for the rest */
		}
	}
	if (exit_status)
		exit(exit_status);
}

/* Start child process if multiprocessing selected */
static pid_t
run_subprocess(void)
{
	int	status;
	pid_t	pid;

	if (nprocs <= 1)			/* any children requested? */
		return(0);
	await_children(nchild + 1 - nprocs);	/* free up child process */
	if ((pid = fork())) {
		if (pid < 0) {
			fprintf(stderr, "%s: cannot fork subprocess\n",
					progname);
			await_children(nchild);
			exit(1);
		}
		++nchild;			/* subprocess started */
		return(pid);
	}
	nchild = -1;
	return(0);				/* put child to work */
}

/* If we are in subprocess, call exit */
#define	end_subprocess()	if (nchild < 0) _exit(0); else

#endif	/* ! _WIN32 */

/* Compute square location from normalized input/output vector */
static void
sq_from_ang(double sq[2], double theta, double phi)
{
	double	vec[3];
	double	norm;
					/* uniform hemispherical projection */
	vec[2] = sin(M_PI/180.*theta);
	vec[0] = cos(M_PI/180.*phi)*vec[2];
	vec[1] = sin(M_PI/180.*phi)*vec[2];
	vec[2] = sqrt(1. - vec[2]*vec[2]);
	norm = 1./sqrt(1. + vec[2]);

	SDdisk2square(sq, vec[0]*norm, vec[1]*norm);
}

/* Compute quantized grid position from normalized input/output vector */
static void
pos_from_ang(int gp[2], double theta, double phi)
{
	double	sq[2];

	sq_from_ang(sq, theta, phi);
	gp[0] = (int)(sq[0]*ANGRES);
	gp[1] = (int)(sq[1]*ANGRES);
}

/* Compare incident angles */
static int
cmp_indir(const void *p1, const void *p2)
{
	const PGINPUT	*inp1 = (const PGINPUT *)p1;
	const PGINPUT	*inp2 = (const PGINPUT *)p2;
	int		ydif;

	if (!inp1->dstart)		/* put terminator last */
		return(inp2->dstart > 0);
	if (!inp2->dstart)
		return(-1);

	ydif = inp1->igp[1] - inp2->igp[1];
	if (ydif)
		return(ydif);

	return(inp1->igp[0] - inp2->igp[0]);
}

/* Prepare a PAB-Opto input file by reading its header */
static int
init_pabopto_inp(PGINPUT *pgi, const char *fname)
{
	FILE	*fp = fopen(fname, "r");
	char	buf[2048];
	int	c;
	
	if (fp == NULL) {
		fputs(fname, stderr);
		fputs(": cannot open\n", stderr);
		return(0);
	}
	pgi->fname = fname;
	pgi->isDSF = -1;
	pgi->ndata = 0;
	pgi->up_phi = 0;
	pgi->theta = pgi->phi = -10001.;
				/* read header information */
	while ((c = getc(fp)) == '#' || c == EOF) {
		char	typ[64];
		if (fgets(buf, sizeof(buf), fp) == NULL) {
			fputs(fname, stderr);
			fputs(": unexpected EOF\n", stderr);
			fclose(fp);
			return(0);
		}
		if (sscanf(buf, "colorimetry: %s", typ) == 1) {
			if (!strcasecmp(typ, "CIE-XYZ")) {
				fputs(fname, stderr);
				fputs(": already in XYZ color space!\n", stderr);
				return(0);
			}
			continue;
		}
		if (sscanf(buf, "datapoints_in_file:%d", &pgi->ndata) == 1)
			continue;
		if (sscanf(buf, "format: theta phi %s", typ) == 1) {
			if (!strcasecmp(typ, "DSF"))
				pgi->isDSF = 1;
			else if (!strcasecmp(typ, "BSDF") ||
					!strcasecmp(typ, "BRDF") ||
					!strcasecmp(typ, "BTDF"))
				pgi->isDSF = 0;
			continue;
		}
		if (sscanf(buf, "upphi %lf", &pgi->up_phi) == 1)
			continue;
		if (sscanf(buf, "intheta %lf", &pgi->theta) == 1)
			continue;
		if (sscanf(buf, "inphi %lf", &pgi->phi) == 1)
			continue;
		if (sscanf(buf, "incident_angle %lf %lf",
				&pgi->theta, &pgi->phi) == 2)
			continue;
	}
	pgi->dstart = ftell(fp) - 1;
	fclose(fp);
	if (pgi->isDSF < 0) {
		fputs(fname, stderr);
		fputs(": unknown format\n", stderr);
		return(0);
	}
	if ((pgi->theta < -10000.) | (pgi->phi < -10000.)) {
		fputs(fname, stderr);
		fputs(": unknown incident angle\n", stderr);
		return(0);
	}
	if (!incident_side) {
		incident_side = (pgi->theta < 90.) ? 1 : -1;
	} else if ((incident_side > 0) ^ (pgi->theta < 90.)) {
		fputs(fname, stderr);
		fputs(": incident on opposite side of surface\n", stderr);
		return(0);
	}
				/* convert angle to grid position */
	pos_from_ang(pgi->igp, pgi->theta, pgi->phi);
	return(1);
}

/* Load interpolating function for a particular incident direction */
static int
load_interp(PGINTERP *pgint, const int igp[2], const PGINPUT *pginp)
{
	int	nv = 0;
	int	nread = 0;
	pgint->ip2 = NULL;
	pgint->va = NULL;
					/* load input file(s) */
	while (pginp->dstart && (pginp->igp[0] == igp[0]) &
				(pginp->igp[1] == igp[1])) {
		FILE	*fp = fopen(pginp->fname, "r");
		double	theta, phi;
		float	val;

		if (fp == NULL || fseek(fp, pginp->dstart, SEEK_SET) < 0) {
			fputs(pginp->fname, stderr);
			fputs(": cannot re-open input file!\n", stderr);
			return(0);
		}
#ifdef DEBUG
		fprintf(stderr, "Loading channel from '%s'\n", pginp->fname);
#endif
		while (fscanf(fp, "%lf %lf %f", &theta, &phi, &val) == 3) {
			double	sq[2];
			if (nread >= nv) {
				if (pginp->ndata > 0)
					nv += pginp->ndata;
				else
					nv += (nv>>1) + 1024;
				pgint->ip2 = interp2_realloc(pgint->ip2, nv);
				pgint->va = (float *)realloc(pgint->va, sizeof(float)*nv);
				if ((pgint->ip2 == NULL) | (pgint->va == NULL))
					goto memerr;
			}
			if (!exiting_side) {
				exiting_side = (theta < 90.) ? 1 : -1;
			} else if ((exiting_side > 0) ^ (theta < 90.)) {
				fputs(pginp->fname, stderr);
				fputs(": exiting measurement on wrong side\n", stderr);
				return(0);
			}
			sq_from_ang(sq, theta, phi);
			pgint->ip2->spt[nread][0] = sq[0]*ANGRES;
			pgint->ip2->spt[nread][1] = sq[1]*ANGRES;
			pgint->va[nread++] = val;
		}
		fclose(fp);
		++pginp;		/* advance to next input */
	}
	if (nv > nread) {		/* fix array sizes */
		pgint->ip2 = interp2_realloc(pgint->ip2, nread);
		pgint->va = (float *)realloc(pgint->va, nread);
		if ((pgint->ip2 == NULL) | (pgint->va == NULL))
			goto memerr;
	}
	return(1);
memerr:
	fputs(progname, stderr);
	fputs(": Out of memory in load_interp()!\n", stderr);
	return(0);
}

/* Interpolate value at the given position */
static double
interp2val(const PGINTERP *pgint, double px, double py)
{
#define	NSMP	36
	float	wt[NSMP];
	int	si[NSMP];
	int	n = interp2_topsamp(wt, si, NSMP, pgint->ip2, px, py);
	double	v = .0;

	while (n-- > 0)
		v += wt[n]*pgint->va[si[n]];

	return(v);
#undef NSMP
}

static void
free_interp(PGINTERP *pgint)
{
	interp2_free(pgint->ip2);
	pgint->ip2 = NULL;
	free(pgint->va);
	pgint->va = NULL;
}

/* Interpolate the given incident direction to an output file */
static int
interp_xyz(const PGINPUT *inp0, int nf, const PGINTERP *int1, const PGINTERP *int2)
{
	int	dtotal = 0;
	char	outfname[256];
	FILE	*ofp;
	int	i;

	if (nf <= 0)
		return(0);
				/* create output file & write header */
	sprintf(outfname, "%s_t%03.0fp%03.0f.txt", basename, inp0->theta, inp0->phi);
	if ((ofp = fopen(outfname, "w")) == NULL) {
		fputs(outfname, stderr);
		fputs(": cannot open for writing\n", stderr);
		return(0);
	}
	fprintf(ofp, "#data written using %s\n", progname);
	fprintf(ofp, "#incident_angle %g %g\n", inp0->theta, inp0->phi);
	fprintf(ofp, "#intheta %g\n", inp0->theta);
	fprintf(ofp, "#inphi %g\n", inp0->phi);
	fprintf(ofp, "#upphi %g\n", inp0->up_phi);
	fprintf(ofp, "#colorimetry: CIE-XYZ\n");
	for (i = nf; i--; )
		dtotal += inp0[i].ndata;
	if (dtotal > 0)
		fprintf(ofp, "#datapoints_in_file: %d\n", dtotal);
	fprintf(ofp, "#format: theta phi %s\n", inp0->isDSF ? "DSF" : "BSDF");
	dtotal = 0;
	while (nf-- > 0) {	/* load channel 1 file(s) */
		FILE	*ifp = fopen(inp0->fname, "r");
		double	theta, phi, val[3], sq[2];

		if (ifp == NULL || fseek(ifp, inp0->dstart, SEEK_SET) < 0) {
			fputs(inp0->fname, stderr);
			fputs(": cannot re-open input file!\n", stderr);
			return(0);
		}
#ifdef DEBUG
		fprintf(stderr, "Adding points from '%s' to '%s'\n",
				inp0->fname, outfname);
#endif
		while (fscanf(ifp, "%lf %lf %lf", &theta, &phi, &val[0]) == 3) {
			if (!exiting_side) {
				exiting_side = (theta < 90.) ? 1 : -1;
			} else if ((exiting_side > 0) ^ (theta < 90.)) {
				fputs(inp0->fname, stderr);
				fputs(": exiting measurement on wrong side\n", stderr);
				return(0);
			}
			sq_from_ang(sq, theta, phi);
			val[1] = interp2val(int1, sq[0]*ANGRES, sq[1]*ANGRES);
			val[2] = interp2val(int2, sq[0]*ANGRES, sq[1]*ANGRES);
			fprintf(ofp, "%.3f %.3f", theta, phi);
			for (i = 0; i < 3; i++)
				fprintf(ofp, " %.4f",
					sens2xyz[i][0]*val[0] +
					sens2xyz[i][1]*val[1] +
					sens2xyz[i][2]*val[2]);
			fputc('\n', ofp);
			++dtotal;
		}
		fclose(ifp);
		++inp0;		/* advance to next input file */
	}
#ifdef DEBUG
	fprintf(stderr, "Wrote %d values to '%s'\n", dtotal, outfname);
#endif
	return(fclose(ofp) == 0);
}

/* Do next incident direction and advance counters */
static int
advance_incidence(PGINPUT *slist[3], int ndx[3])
{
	const int	i0 = ndx[0];
	PGINTERP	pgi1, pgi2;
	int		c;
				/* match up directions */
	while ((c = cmp_indir(slist[1]+ndx[1], slist[0]+ndx[0])) < 0)
		ndx[1]++;
	if (c) {
		fputs(slist[0][ndx[0]].fname, stderr);
		fputs(": warning - missing channel 2 at this incidence\n", stderr);
		do
			++ndx[0];
		while (!cmp_indir(slist[0]+ndx[0], slist[0]+i0));
		return(0);
	}
	while ((c = cmp_indir(slist[2]+ndx[2], slist[0]+ndx[0])) < 0)
		ndx[2]++;
	if (c) {
		fputs(slist[0][ndx[0]].fname, stderr);
		fputs(": warning - missing channel 3 at this incidence\n", stderr);
		do
			++ndx[0];
		while (!cmp_indir(slist[0]+ndx[0], slist[0]+i0));
		return(0);
	}
	do			/* advance first channel index */
		++ndx[0];
	while (!cmp_indir(slist[0]+ndx[0], slist[0]+i0));

	if (run_subprocess())	/* fork here if multiprocessing */
		return(1);
				/* interpolate channels 2 & 3 */
	load_interp(&pgi1, slist[0][i0].igp, slist[1]+ndx[1]);
	load_interp(&pgi2, slist[0][i0].igp, slist[2]+ndx[2]);
	if (!interp_xyz(slist[0]+i0, ndx[0]-i0, &pgi1, &pgi2))
		return(-1);
	end_subprocess();
	free_interp(&pgi1);
	free_interp(&pgi2);
	return(1);
}

/* Read in single-channel PAB-Opto BSDF files and output XYZ versions */
int
main(int argc, char *argv[])
{
	PGINPUT	*slist[3];
	int	i, j;
	int	ndx[3];
						/* get options */
	progname = argv[0];
	for (i = 1; i < argc-3 && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'n':
			nprocs = atoi(argv[++i]);
			break;
		case 'm':
			if (i >= argc-(3+9))
				goto userr;
			for (j = 0; j < 3; j++) {
				sens2xyz[j][0] = atof(argv[++i]);
				sens2xyz[j][1] = atof(argv[++i]);
				sens2xyz[j][2] = atof(argv[++i]);
			}
			break;
		case 'o':
			basename = argv[++i];
			break;
		default:
			goto userr;
		}
	if (i != argc-3)
		goto userr;
	for (j = 0; j < 3; j++) {		/* prep input channels */
		char	*flist[MAX_INPUTS];
		int	k, n;
		n = wordfile(flist, MAX_INPUTS, argv[i+j]);
		if (n <= 0) {
			fputs(argv[i+j], stderr);
			fputs(": cannot load input file names\n", stderr);
			return(1);
		}
		slist[j] = (PGINPUT *)malloc(sizeof(PGINPUT)*(n+1));
		if (slist[j] == NULL) {
			fputs(argv[0], stderr);
			fputs(": out of memory!\n", stderr);
			return(1);
		}
#ifdef DEBUG
		fprintf(stderr, "Checking %d files from '%s'\n", n, argv[i+j]);
#endif
		for (k = 0; k < n; k++)
			if (!init_pabopto_inp(slist[j]+k, flist[k]))
				return(1);
		memset(slist[j]+n, 0, sizeof(PGINPUT));
						/* sort by incident direction */
		qsort(slist[j], n, sizeof(PGINPUT), cmp_indir);
	}
	ndx[0]=ndx[1]=ndx[2]=0;			/* compile measurements */
	while (slist[0][ndx[0]].dstart)
		if (advance_incidence(slist, ndx) < 0)
			return(1);
	await_children(nchild);			/* wait for all to finish */
	return(0);
userr:
	fputs("Usage: ", stderr);
	fputs(argv[0], stderr);
	fputs(" [-m X1 X2 X3 Y1 Y2 Y3 Z1 Z2 Z3][-o basename][-n nprocs]", stderr);
	fputs(" s1files.txt s2files.txt s3files.txt\n", stderr);
	return(1);
}
