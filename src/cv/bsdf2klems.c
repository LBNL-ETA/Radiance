#ifndef lint
static const char RCSid[] = "$Id: bsdf2klems.c,v 2.20 2016/02/03 01:57:06 greg Exp $";
#endif
/*
 * Load measured BSDF interpolant and write out as XML file with Klems matrix.
 *
 *	G. Ward
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "random.h"
#include "platform.h"
#include "paths.h"
#include "rtio.h"
#include "calcomp.h"
#include "bsdfrep.h"
#include "bsdf_m.h"
				/* tristimulus components */
enum {CIE_X, CIE_Y, CIE_Z};
				/* assumed maximum # Klems patches */
#define MAXPATCHES	145
				/* global argv[0] */
char			*progname;
				/* selected basis function name */
static const char	klems_full[] = "LBNL/Klems Full";
static const char	klems_half[] = "LBNL/Klems Half";
static const char	klems_quarter[] = "LBNL/Klems Quarter";
static const char	*kbasis = klems_full;
				/* number of BSDF samples per patch */
static int		npsamps = 256;
				/* limit on number of RBF lobes */
static int		lobe_lim = 15000;
				/* progress bar length */
static int		do_prog = 79;

#define MAXCARG		512	/* wrapBSDF command */
static char		*wrapBSDF[MAXCARG] = {"wrapBSDF", "-W", "-UU"};
static int		wbsdfac = 3;

/* Add argument to wrapBSDF, allocating space if isstatic */
static void
add_wbsdf(const char *arg, int isstatic)
{
	if (arg == NULL)
		return;
	if (wbsdfac >= MAXCARG-1) {
		fputs(progname, stderr);
		fputs(": too many command arguments to wrapBSDF\n", stderr);
		exit(1);
	}
	if (!*arg)
		arg = "";
	else if (!isstatic)
		arg = savqstr((char *)arg);

	wrapBSDF[wbsdfac++] = (char *)arg;
}

/* Start new progress bar */
#define prog_start(s)	if (do_prog) fprintf(stderr, "%s: %s...\n", progname, s); else

/* Draw progress bar of the appropriate length */
static void
prog_show(double frac)
{
	static unsigned	call_cnt = 0;
	static char	lastc[] = "-\\|/";
	char		pbar[256];
	int		nchars;

	if (do_prog <= 1) return;
	if (do_prog > sizeof(pbar)-2)
		do_prog = sizeof(pbar)-2;
	if (frac < 0) frac = 0;
	else if (frac >= 1) frac = .9999;
	nchars = do_prog*frac;
	pbar[0] = '\r';
	memset(pbar+1, '*', nchars);
	pbar[nchars+1] = lastc[call_cnt++ & 3];
	memset(pbar+2+nchars, '-', do_prog-nchars-1);
	pbar[do_prog+1] = '\0';
	fputs(pbar, stderr);
}

/* Finish progress bar */
static void
prog_done(void)
{
	int	n = do_prog;

	if (n <= 1) return;
	fputc('\r', stderr);
	while (n--)
		fputc(' ', stderr);
	fputc('\r', stderr);
}

/* Return angle basis corresponding to the given name */
static ANGLE_BASIS *
get_basis(const char *bn)
{
	int	n = nabases;
	
	while (n-- > 0)
		if (!strcasecmp(bn, abase_list[n].name))
			return &abase_list[n];
	return NULL;
}

/* Copy geometry string to file for wrapBSDF */
static char *
save_geom(const char *mgf)
{
	char	*tfname = mktemp(savqstr(TEMPLATE));
	int	fd = open(tfname, O_CREAT|O_WRONLY, 0600);

	if (fd < 0)
		return(NULL);
	write(fd, mgf, strlen(mgf));
	close(fd);
	add_wbsdf("-g", 1);
	add_wbsdf(tfname, 1);
	return(tfname);
}

/* Open XYZ component file for output and add appropriate arguments */
static FILE *
open_component_file(int c)
{
	static const char	sname[3][6] = {"CIE-X", "CIE-Y", "CIE-Z"};
	static const char	cname[4][4] = {"-rf", "-tf", "-tb", "-rb"};
	char			*tfname = mktemp(savqstr(TEMPLATE));
	FILE			*fp = fopen(tfname, "w");

	if (fp == NULL) {
		fprintf(stderr, "%s: cannot open '%s' for writing\n",
				progname, tfname);
		exit(1);
	}
	add_wbsdf("-s", 1); add_wbsdf(sname[c], 1);
	add_wbsdf(cname[(input_orient>0)<<1 | (output_orient>0)], 1);
	add_wbsdf(tfname, 1);
	return(fp);
}

/* Load and resample XML BSDF description using Klems basis */
static void
eval_bsdf(const char *fname)
{
	ANGLE_BASIS	*abp = get_basis(kbasis);
	FILE		*cfp[3];
	SDData		bsd;
	SDError		ec;
	FVECT		vin, vout;
	SDValue		sdv;
	double		sum, xsum, ysum;
	int		i, j, n;

	initurand(npsamps);
	SDclearBSDF(&bsd, fname);		/* load BSDF file */
	if ((ec = SDloadFile(&bsd, fname)) != SDEnone)
		goto err;
	if (bsd.mgf != NULL)			/* save geometry */
		save_geom(bsd.mgf);
	if (bsd.matn[0])			/* save identifier(s) */
		strcpy(bsdf_name, bsd.matn);
	if (bsd.makr[0])
		strcpy(bsdf_manuf, bsd.makr);
	if (bsd.dim[2] > 0) {			/* save dimension(s) */
		char	buf[64];
		if ((bsd.dim[0] > 0) & (bsd.dim[1] > 0))
			sprintf(buf, "w=%g;h=%g;t=%g",
					bsd.dim[0], bsd.dim[1], bsd.dim[2]);
		else
			sprintf(buf, "t=%g", bsd.dim[2]);
		add_wbsdf("-f", 1);
		add_wbsdf(buf, 0);
	}
						/* front reflection */
	if (bsd.rf != NULL || bsd.rLambFront.cieY > .002) {
	    input_orient = 1; output_orient = 1;
	    cfp[CIE_Y] = open_component_file(CIE_Y);
	    if (bsd.rf != NULL && bsd.rf->comp[0].cspec[2].flags) {
		rbf_colorimetry = RBCtristimulus;
		cfp[CIE_X] = open_component_file(CIE_X);
		cfp[CIE_Z] = open_component_file(CIE_Z);
	    } else
		rbf_colorimetry = RBCphotopic;
	    for (j = 0; j < abp->nangles; j++) {
	        for (i = 0; i < abp->nangles; i++) {
		    sum = 0;			/* average over patches */
		    xsum = ysum = 0;
		    for (n = npsamps; n-- > 0; ) {
			fo_getvec(vout, j+(n+frandom())/npsamps, abp);
			fi_getvec(vin, i+urand(n), abp);
			ec = SDevalBSDF(&sdv, vout, vin, &bsd);
			if (ec != SDEnone)
				goto err;
			sum += sdv.cieY;
			if (rbf_colorimetry == RBCtristimulus) {
				c_ccvt(&sdv.spec, C_CSXY);
				xsum += sdv.cieY*sdv.spec.cx;
				ysum += sdv.cieY*sdv.spec.cy;
			}
		    }
		    fprintf(cfp[CIE_Y], "\t%.3e\n", sum/npsamps);
		    if (rbf_colorimetry == RBCtristimulus) {
			fprintf(cfp[CIE_X], "\t%3e\n", xsum*sum/(npsamps*ysum));
			fprintf(cfp[CIE_Z], "\t%3e\n",
				(sum - xsum - ysum)*sum/(npsamps*ysum));
		    }
		}
		fputc('\n', cfp[CIE_Y]);	/* extra space between rows */
		if (rbf_colorimetry == RBCtristimulus) {
			fputc('\n', cfp[CIE_X]);
			fputc('\n', cfp[CIE_Z]);
		}
	    }
	    if (fclose(cfp[CIE_Y])) {
		fprintf(stderr, "%s: error writing Y output\n", progname);
		exit(1);
	    }
	    if (rbf_colorimetry == RBCtristimulus &&
			(fclose(cfp[CIE_X]) || fclose(cfp[CIE_Z]))) {
		fprintf(stderr, "%s: error writing X/Z output\n", progname);
		exit(1);
	    }
	}
						/* back reflection */
	if (bsd.rb != NULL || bsd.rLambBack.cieY > .002) {
	    input_orient = -1; output_orient = -1;
	    cfp[CIE_Y] = open_component_file(CIE_Y);
	    if (bsd.rb != NULL && bsd.rb->comp[0].cspec[2].flags) {
		rbf_colorimetry = RBCtristimulus;
		cfp[CIE_X] = open_component_file(CIE_X);
		cfp[CIE_Z] = open_component_file(CIE_Z);
	    } else
		rbf_colorimetry = RBCphotopic;
	    for (j = 0; j < abp->nangles; j++) {
	        for (i = 0; i < abp->nangles; i++) {
		    sum = 0;			/* average over patches */
		    xsum = ysum = 0;
		    for (n = npsamps; n-- > 0; ) {
			bo_getvec(vout, j+(n+frandom())/npsamps, abp);
			bi_getvec(vin, i+urand(n), abp);
			ec = SDevalBSDF(&sdv, vout, vin, &bsd);
			if (ec != SDEnone)
				goto err;
			sum += sdv.cieY;
			if (rbf_colorimetry == RBCtristimulus) {
				c_ccvt(&sdv.spec, C_CSXY);
				xsum += sdv.cieY*sdv.spec.cx;
				ysum += sdv.cieY*sdv.spec.cy;
			}
		    }
		    fprintf(cfp[CIE_Y], "\t%.3e\n", sum/npsamps);
		    if (rbf_colorimetry == RBCtristimulus) {
			fprintf(cfp[CIE_X], "\t%3e\n", xsum*sum/(npsamps*ysum));
			fprintf(cfp[CIE_Z], "\t%3e\n",
				(sum - xsum - ysum)*sum/(npsamps*ysum));
		    }
		}
		if (rbf_colorimetry == RBCtristimulus) {
			fputc('\n', cfp[CIE_X]);
			fputc('\n', cfp[CIE_Z]);
		}
	    }
	    if (fclose(cfp[CIE_Y])) {
		fprintf(stderr, "%s: error writing Y output\n", progname);
		exit(1);
	    }
	    if (rbf_colorimetry == RBCtristimulus &&
			(fclose(cfp[CIE_X]) || fclose(cfp[CIE_Z]))) {
		fprintf(stderr, "%s: error writing X/Z output\n", progname);
		exit(1);
	    }
	}
						/* front transmission */
	if (bsd.tf != NULL || bsd.tLamb.cieY > .002) {
	    input_orient = 1; output_orient = -1;
	    cfp[CIE_Y] = open_component_file(CIE_Y);
	    if (bsd.tf != NULL && bsd.tf->comp[0].cspec[2].flags) {
		rbf_colorimetry = RBCtristimulus;
		cfp[CIE_X] = open_component_file(CIE_X);
		cfp[CIE_Z] = open_component_file(CIE_Z);
	    } else
		rbf_colorimetry = RBCphotopic;
	    for (j = 0; j < abp->nangles; j++) {
	        for (i = 0; i < abp->nangles; i++) {
		    sum = 0;			/* average over patches */
		    xsum = ysum = 0;
		    for (n = npsamps; n-- > 0; ) {
			bo_getvec(vout, j+(n+frandom())/npsamps, abp);
			fi_getvec(vin, i+urand(n), abp);
			ec = SDevalBSDF(&sdv, vout, vin, &bsd);
			if (ec != SDEnone)
				goto err;
			sum += sdv.cieY;
			if (rbf_colorimetry == RBCtristimulus) {
				c_ccvt(&sdv.spec, C_CSXY);
				xsum += sdv.cieY*sdv.spec.cx;
				ysum += sdv.cieY*sdv.spec.cy;
			}
		    }
		    fprintf(cfp[CIE_Y], "\t%.3e\n", sum/npsamps);
		    if (rbf_colorimetry == RBCtristimulus) {
			fprintf(cfp[CIE_X], "\t%3e\n", xsum*sum/(npsamps*ysum));
			fprintf(cfp[CIE_Z], "\t%3e\n",
				(sum - xsum - ysum)*sum/(npsamps*ysum));
		    }
		}
		if (rbf_colorimetry == RBCtristimulus) {
			fputc('\n', cfp[CIE_X]);
			fputc('\n', cfp[CIE_Z]);
		}
	    }
	    if (fclose(cfp[CIE_Y])) {
		fprintf(stderr, "%s: error writing Y output\n", progname);
		exit(1);
	    }
	    if (rbf_colorimetry == RBCtristimulus &&
			(fclose(cfp[CIE_X]) || fclose(cfp[CIE_Z]))) {
		fprintf(stderr, "%s: error writing X/Z output\n", progname);
		exit(1);
	    }
	}
						/* back transmission */
	if ((bsd.tb != NULL) | (bsd.tf != NULL)) {
	    input_orient = -1; output_orient = 1;
	    cfp[CIE_Y] = open_component_file(CIE_Y);
	    if (bsd.tb != NULL && bsd.tb->comp[0].cspec[2].flags) {
		rbf_colorimetry = RBCtristimulus;
		cfp[CIE_X] = open_component_file(CIE_X);
		cfp[CIE_Z] = open_component_file(CIE_Z);
	    } else
		rbf_colorimetry = RBCphotopic;
	    for (j = 0; j < abp->nangles; j++) {
	        for (i = 0; i < abp->nangles; i++) {
		    sum = 0;		/* average over patches */
		    xsum = ysum = 0;
		    for (n = npsamps; n-- > 0; ) {
			fo_getvec(vout, j+(n+frandom())/npsamps, abp);
			bi_getvec(vin, i+urand(n), abp);
			ec = SDevalBSDF(&sdv, vout, vin, &bsd);
			if (ec != SDEnone)
				goto err;
			sum += sdv.cieY;
			if (rbf_colorimetry == RBCtristimulus) {
				c_ccvt(&sdv.spec, C_CSXY);
				xsum += sdv.cieY*sdv.spec.cx;
				ysum += sdv.cieY*sdv.spec.cy;
			}
		    }
		    fprintf(cfp[CIE_Y], "\t%.3e\n", sum/npsamps);
		    if (rbf_colorimetry == RBCtristimulus) {
			fprintf(cfp[CIE_X], "\t%3e\n", xsum*sum/(npsamps*ysum));
			fprintf(cfp[CIE_Z], "\t%3e\n",
				(sum - xsum - ysum)*sum/(npsamps*ysum));
		    }
		}
		if (rbf_colorimetry == RBCtristimulus) {
			fputc('\n', cfp[CIE_X]);
			fputc('\n', cfp[CIE_Z]);
		}
	    }
	    if (fclose(cfp[CIE_Y])) {
		fprintf(stderr, "%s: error writing Y output\n", progname);
		exit(1);
	    }
	    if (rbf_colorimetry == RBCtristimulus &&
			(fclose(cfp[CIE_X]) || fclose(cfp[CIE_Z]))) {
		fprintf(stderr, "%s: error writing X/Z output\n", progname);
		exit(1);
	    }
	}
	SDfreeBSDF(&bsd);			/* all done */
	return;
err:
	SDreportError(ec, stderr);
	exit(1);
}

/* Interpolate and output a BSDF function using Klems basis */
static void
eval_function(char *funame)
{
	ANGLE_BASIS	*abp = get_basis(kbasis);
	int		assignD = (fundefined(funame) < 6);
	FILE		*ofp = open_component_file(CIE_Y);
	double		iovec[6];
	double		sum;
	int		i, j, n;

	initurand(npsamps);
	for (j = 0; j < abp->nangles; j++) {	/* run through directions */
	    for (i = 0; i < abp->nangles; i++) {
		sum = 0;
		for (n = npsamps; n--; ) {	/* average over patches */
		    if (output_orient > 0)
			fo_getvec(iovec+3, j+(n+frandom())/npsamps, abp);
		    else
			bo_getvec(iovec+3, j+(n+frandom())/npsamps, abp);

		    if (input_orient > 0)
			fi_getvec(iovec, i+urand(n), abp);
		    else
			bi_getvec(iovec, i+urand(n), abp);

		    if (assignD) {
			varset("Dx", '=', -iovec[3]);
			varset("Dy", '=', -iovec[4]);
			varset("Dz", '=', -iovec[5]);
			++eclock;
		    }
		    sum += funvalue(funame, 6, iovec);
		}
		fprintf(ofp, "\t%.3e\n", sum/npsamps);
	    }
	    fputc('\n', ofp);
	    prog_show((j+1.)/abp->nangles);
	}
	prog_done();
	if (fclose(ofp)) {
		fprintf(stderr, "%s: error writing Y output\n", progname);
		exit(1);
	}
}

/* Interpolate and output a radial basis function BSDF representation */
static void
eval_rbf(void)
{
	ANGLE_BASIS	*abp = get_basis(kbasis);
	float		(*XZarr)[2] = NULL;
	float		bsdfarr[MAXPATCHES*MAXPATCHES];
	FILE		*cfp[3];
	FVECT		vin, vout;
	double		sum, xsum, ysum;
	int		i, j, n;
						/* sanity check */
	if (abp->nangles > MAXPATCHES) {
		fprintf(stderr, "%s: too many patches!\n", progname);
		exit(1);
	}
	if (rbf_colorimetry == RBCtristimulus)
		XZarr = (float (*)[2])malloc(sizeof(float)*2*abp->nangles*abp->nangles);
	for (i = 0; i < abp->nangles; i++) {
	    RBFNODE	*rbf;
	    if (input_orient > 0)		/* use incident patch center */
		fi_getvec(vin, i+.5*(i>0), abp);
	    else
		bi_getvec(vin, i+.5*(i>0), abp);

	    rbf = advect_rbf(vin, lobe_lim);	/* compute radial basis func */

	    for (j = 0; j < abp->nangles; j++) {
	        sum = 0;			/* sample over exiting patch */
		xsum = ysum = 0;
		for (n = npsamps; n--; ) {
		    SDValue	sdv;
		    if (output_orient > 0)
			fo_getvec(vout, j+(n+frandom())/npsamps, abp);
		    else
			bo_getvec(vout, j+(n+frandom())/npsamps, abp);

		    eval_rbfcol(&sdv, rbf, vout);
		    sum += sdv.cieY;
		    if (XZarr != NULL) {
			c_ccvt(&sdv.spec, C_CSXY);
			xsum += sdv.cieY*sdv.spec.cx;
			ysum += sdv.cieY*sdv.spec.cy;
		    }
		}
		n = j*abp->nangles + i;
		bsdfarr[n] = sum / (double)npsamps;
		if (XZarr != NULL) {
		    XZarr[n][0] = xsum*sum/(npsamps*ysum);
		    XZarr[n][1] = (sum - xsum - ysum)*sum/(npsamps*ysum);
		}
	    }
	    if (rbf != NULL)
		free(rbf);
	    prog_show((i+1.)/abp->nangles);
	}
						/* write out our matrix */
	cfp[CIE_Y] = open_component_file(CIE_Y);
	n = 0;
	for (j = 0; j < abp->nangles; j++) {
	    for (i = 0; i < abp->nangles; i++, n++)
		fprintf(cfp[CIE_Y], "\t%.3e\n", bsdfarr[n]);
	    fputc('\n', cfp[CIE_Y]);
	}
	prog_done();
	if (fclose(cfp[CIE_Y])) {
		fprintf(stderr, "%s: error writing Y output\n", progname);
		exit(1);
	}
	if (XZarr == NULL)			/* no color? */
		return;
	cfp[CIE_X] = open_component_file(CIE_X);
	cfp[CIE_Z] = open_component_file(CIE_Z);
	n = 0;
	for (j = 0; j < abp->nangles; j++) {
	    for (i = 0; i < abp->nangles; i++, n++) {
		fprintf(cfp[CIE_X], "\t%.3e\n", XZarr[n][0]);
		fprintf(cfp[CIE_Z], "\t%.3e\n", XZarr[n][1]);
	    }
	    fputc('\n', cfp[CIE_X]);
	    fputc('\n', cfp[CIE_Z]);
	}
	free(XZarr);
	if (fclose(cfp[CIE_X]) || fclose(cfp[CIE_Z])) {
		fprintf(stderr, "%s: error writing X/Z output\n", progname);
		exit(1);
	}
}

#ifdef _WIN32
/* Execute wrapBSDF command (may never return) */
static int
wrap_up(void)
{
	char	cmd[8192];

	if (bsdf_manuf[0]) {
		add_wbsdf("-f", 1);
		strcpy(cmd, "m=");
		strcpy(cmd+2, bsdf_manuf);
		add_wbsdf(cmd, 0);
	}
	if (bsdf_name[0]) {
		add_wbsdf("-f", 1);
		strcpy(cmd, "n=");
		strcpy(cmd+2, bsdf_name);
		add_wbsdf(cmd, 0);
	}
	if (!convert_commandline(cmd, sizeof(cmd), wrapBSDF)) {
		fputs(progname, stderr);
		fputs(": command line too long in wrap_up()\n", stderr);
		return(1);
	}
	return(system(cmd));
}
#else
/* Execute wrapBSDF command (may never return) */
static int
wrap_up(void)
{
	char	buf[256];
	char	*compath = getpath((char *)wrapBSDF[0], getenv("PATH"), X_OK);

	if (compath == NULL) {
		fprintf(stderr, "%s: cannot locate %s\n", progname, wrapBSDF[0]);
		return(1);
	}
	if (bsdf_manuf[0]) {
		add_wbsdf("-f", 1);
		strcpy(buf, "m=");
		strcpy(buf+2, bsdf_manuf);
		add_wbsdf(buf, 0);
	}
	if (bsdf_name[0]) {
		add_wbsdf("-f", 1);
		strcpy(buf, "n=");
		strcpy(buf+2, bsdf_name);
		add_wbsdf(buf, 0);
	}
	execv(compath, wrapBSDF);	/* successful call never returns */
	perror(compath);
	return(1);
}
#endif

/* Read in BSDF and interpolate as Klems matrix representation */
int
main(int argc, char *argv[])
{
	int	dofwd = 0, dobwd = 1;
	char	buf[2048];
	char	*cp;
	int	i, na;

	progname = argv[0];
	esupport |= E_VARIABLE|E_FUNCTION|E_RCONST;
	esupport &= ~(E_INCHAN|E_OUTCHAN);
	scompile("PI:3.14159265358979323846", NULL, 0);
	biggerlib();
	for (i = 1; i < argc && (argv[i][0] == '-') | (argv[i][0] == '+'); i++)
		switch (argv[i][1]) {		/* get options */
		case 'n':
			npsamps = atoi(argv[++i]);
			if (npsamps <= 0)
				goto userr;
			break;
		case 'e':
			scompile(argv[++i], NULL, 0);
			single_plane_incident = 0;
			break;
		case 'f':
			if (!argv[i][2]) {
				if (strchr(argv[++i], '=') != NULL) {
					add_wbsdf("-f", 1);
					add_wbsdf(argv[i], 1);
				} else {
					fcompile(argv[i]);
					single_plane_incident = 0;
				}
			} else
				dofwd = (argv[i][0] == '+');
			break;
		case 'b':
			dobwd = (argv[i][0] == '+');
			break;
		case 'h':
			kbasis = klems_half;
			add_wbsdf("-a", 1);
			add_wbsdf("kh", 1);
			break;
		case 'q':
			kbasis = klems_quarter;
			add_wbsdf("-a", 1);
			add_wbsdf("kq", 1);
			break;
		case 'l':
			lobe_lim = atoi(argv[++i]);
			break;
		case 'p':
			do_prog = atoi(argv[i]+2);
			break;
		case 'C':
			add_wbsdf(argv[i], 1);
			add_wbsdf(argv[++i], 1);
			break;
		default:
			goto userr;
		}
	if (kbasis == klems_full) {		/* default (full) basis? */
		add_wbsdf("-a", 1);
		add_wbsdf("kf", 1);
	}
	strcpy(buf, "File produced by: ");
	if (convert_commandline(buf+18, sizeof(buf)-18, argv) != NULL) {
		add_wbsdf("-C", 1); add_wbsdf(buf, 0);
	}
	if (single_plane_incident >= 0) {	/* function-based BSDF? */
		if (i != argc-1 || fundefined(argv[i]) < 3) {
			fprintf(stderr,
	"%s: need single function with 6 arguments: bsdf(ix,iy,iz,ox,oy,oz)\n",
					progname);
			fprintf(stderr, "\tor 3 arguments using Dx,Dy,Dz: bsdf(ix,iy,iz)\n");
			goto userr;
		}
		++eclock;
		if (dofwd) {
			input_orient = -1;
			output_orient = -1;
			prog_start("Evaluating outside reflectance");
			eval_function(argv[i]);
			output_orient = 1;
			prog_start("Evaluating outside->inside transmission");
			eval_function(argv[i]);
		}
		if (dobwd) {
			input_orient = 1;
			output_orient = 1;
			prog_start("Evaluating inside reflectance");
			eval_function(argv[i]);
			output_orient = -1;
			prog_start("Evaluating inside->outside transmission");
			eval_function(argv[i]);
		}
		return(wrap_up());
	}
						/* XML input? */
	if (i == argc-1 && (cp = argv[i]+strlen(argv[i])-4) > argv[i] &&
				!strcasecmp(cp, ".xml")) {
		eval_bsdf(argv[i]);		/* load & resample BSDF */
		return(wrap_up());
	}
	if (i < argc) {				/* open input files if given */
		int	nbsdf = 0;
		for ( ; i < argc; i++) {	/* interpolate each component */
			char	pbuf[256];
			FILE	*fpin = fopen(argv[i], "rb");
			if (fpin == NULL) {
				fprintf(stderr, "%s: cannot open BSDF interpolant '%s'\n",
						progname, argv[i]);
				return(1);
			}
			if (!load_bsdf_rep(fpin))
				return(1);
			fclose(fpin);
			sprintf(pbuf, "Interpolating component '%s'", argv[i]);
			prog_start(pbuf);
			eval_rbf();
		}
		return(wrap_up());
	}
	SET_FILE_BINARY(stdin);			/* load from stdin */
	if (!load_bsdf_rep(stdin))
		return(1);
	prog_start("Interpolating from standard input");
	eval_rbf();				/* resample dist. */
	return(wrap_up());
userr:
	fprintf(stderr,
	"Usage: %s [-n spp][-h|-q][-l maxlobes] [bsdf.sir ..] > bsdf.xml\n", progname);
	fprintf(stderr,
	"   or: %s [-n spp][-h|-q] bsdf_in.xml > bsdf_out.xml\n", progname);
	fprintf(stderr,
	"   or: %s [-n spp][-h|-q][{+|-}for[ward]][{+|-}b[ackward]][-e expr][-f file] bsdf_func > bsdf.xml\n",
				progname);
	return(1);
}
