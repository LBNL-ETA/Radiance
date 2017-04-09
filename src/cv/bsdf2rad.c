#ifndef lint
static const char RCSid[] = "$Id: bsdf2rad.c,v 2.18 2017/04/09 22:51:19 greg Exp $";
#endif
/*
 *  Plot 3-D BSDF output based on scattering interpolant or XML representation
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "paths.h"
#include "rtmath.h"
#include "resolu.h"
#include "bsdfrep.h"

#define NINCIDENT	37		/* number of samples/hemisphere */

#define	GRIDSTEP	2		/* our grid step size */
#define SAMPRES		(GRIDRES/GRIDSTEP)

int	front_comp = 0;			/* front component flags (SDsamp*) */
int	back_comp = 0;			/* back component flags */
double	overall_min = 1./PI;		/* overall minimum BSDF value */
double	min_log10;			/* smallest log10 value for plotting */
double	overall_max = .0;		/* overall maximum BSDF value */

char	ourTempDir[TEMPLEN] = "";	/* our temporary directory */

const FVECT	Xaxis = {1., 0., 0.};
const FVECT	Yaxis = {0., 1., 0.};
const FVECT	Zaxis = {0., 0., 1.};

const char	frpref[] = "frefl";
const char	ftpref[] = "ftrans";
const char	brpref[] = "brefl";
const char	btpref[] = "btrans";
const char	dsuffix[] = ".txt";

const char	sph_mat[] = "BSDFmat";
const double	sph_rad = 10.;
const double	sph_xoffset = 15.;

#define bsdf_rad	(sph_rad*.25)
#define arrow_rad	(bsdf_rad*.015)

#define	FEQ(a,b)	((a)-(b) <= 1e-7 && (b)-(a) <= 1e-7)

#define	set_minlog()	(min_log10 = log10(overall_min + 1e-5) - .1)

char	*progname;

/* Get Fibonacci sphere vector (0 to NINCIDENT-1) */
static void
get_ivector(FVECT iv, int i)
{
	const double	phistep = PI*(3. - 2.236067978);
	double		r;

	iv[2] = 1. - (i+.5)*(1./NINCIDENT);
	r = sqrt(1. - iv[2]*iv[2]);
	iv[0] = r * cos((i+1.)*phistep);
	iv[1] = r * sin((i+1.)*phistep);
}

/* Get temporary file name */
static char *
tfile_name(const char *prefix, const char *suffix, int i)
{
	static char	buf[128];

	if (!ourTempDir[0]) {		/* create temporary directory */
		mktemp(strcpy(ourTempDir,TEMPLATE));
		if (mkdir(ourTempDir, 0777) < 0) {
			perror("mkdir");
			exit(1);
		}
	}
	if (!prefix) prefix = "T";
	if (!suffix) suffix = "";
	sprintf(buf, "%s/%s%03d%s", ourTempDir, prefix, i, suffix);
	return(buf);
}

/* Remove temporary directory & contents */
static void
cleanup_tmp(void)
{
	char	buf[128];

	if (!ourTempDir[0])
		return;
#if defined(_WIN32) || defined(_WIN64)
	sprintf(buf, "RMDIR %s /S /Q", ourTempDir);
#else
	sprintf(buf, "rm -rf %s", ourTempDir);
#endif
	system(buf);
}

/* Run the specified command, returning 1 if OK */
static int
run_cmd(const char *cmd)
{
	fflush(stdout);
	if (system(cmd)) {
		fprintf(stderr, "%s: error running: %s\n", progname, cmd);
		return(0);
	}
	return(1);
}

/* Plot surface points for the given BSDF incident angle */
static int
plotBSDF(const char *fname, const FVECT ivec, int dfl, const SDData *sd)
{
	FILE		*fp = fopen(fname, "w");
	int		i, j;

	if (fp == NULL) {
		fprintf(stderr, "%s: cannot open '%s' for writing\n",
				progname, fname);
		return(0);
	}
	if (ivec[2] > 0) {
		input_orient = 1;
		output_orient = dfl&SDsampR ? 1 : -1;
	} else {
		input_orient = -1;
		output_orient = dfl&SDsampR ? -1 : 1;
	}
	for (i = SAMPRES; i--; )
	    for (j = 0; j < SAMPRES; j++) {
		FVECT	ovec;
		SDValue	sval;
		double	bsdf;
		ovec_from_pos(ovec, i*GRIDSTEP, j*GRIDSTEP);
		if (SDreportError(SDevalBSDF(&sval, ovec,
						ivec, sd), stderr))
			return(0);
		if (sval.cieY > overall_max)
			overall_max = sval.cieY;
		bsdf = (sval.cieY < overall_min) ? overall_min : sval.cieY;
		bsdf = log10(bsdf) - min_log10;
		fprintf(fp, "%.5f %.5f %.5f\n",
				ovec[0]*bsdf, ovec[1]*bsdf, ovec[2]*bsdf);
	    }
	if (fclose(fp) == EOF) {
		fprintf(stderr, "%s: error writing data to '%s'\n",
				progname, fname);
		return(0);
	}
	return(1);
}

/* Build BSDF values from loaded XML file */
static int
build_wBSDF(const SDData *sd)
{
	FVECT	ivec;
	int	i;

	if (front_comp & SDsampR)
		for (i = 0; i < NINCIDENT; i++) {
			get_ivector(ivec, i);
			if (!plotBSDF(tfile_name(frpref, dsuffix, i),
					ivec, SDsampR, sd))
				return(0);
		}
	if (front_comp & SDsampT)
		for (i = 0; i < NINCIDENT; i++) {
			get_ivector(ivec, i);
			if (!plotBSDF(tfile_name(ftpref, dsuffix, i),
					ivec, SDsampT, sd))
				return(0);
		}
	if (back_comp & SDsampR)
		for (i = 0; i < NINCIDENT; i++) {
			get_ivector(ivec, i);
			ivec[0] = -ivec[0]; ivec[2] = -ivec[2];
			if (!plotBSDF(tfile_name(brpref, dsuffix, i),
					ivec, SDsampR, sd))
				return(0);
		}
	if (back_comp & SDsampT)
		for (i = 0; i < NINCIDENT; i++) {
			get_ivector(ivec, i);
			ivec[0] = -ivec[0]; ivec[2] = -ivec[2];
			if (!plotBSDF(tfile_name(btpref, dsuffix, i),
					ivec, SDsampT, sd))
				return(0);
		}
	return(1);
}

/* Plot surface points using radial basis function */
static int
plotRBF(const char *fname, const RBFNODE *rbf)
{
	FILE		*fp = fopen(fname, "w");
	int		i, j;

	if (fp == NULL) {
		fprintf(stderr, "%s: cannot open '%s' for writing\n",
				progname, fname);
		return(0);
	}
	for (i = SAMPRES; i--; )
	    for (j = 0; j < SAMPRES; j++) {
		FVECT	ovec;
		double	bsdf;
		ovec_from_pos(ovec, i*GRIDSTEP, j*GRIDSTEP);
		bsdf = eval_rbfrep(rbf, ovec);
		if (bsdf > overall_max)
			overall_max = bsdf;
		else if (bsdf < overall_min)
			bsdf = overall_min;
		bsdf = log10(bsdf) - min_log10;
		fprintf(fp, "%.5f %.5f %.5f\n",
				ovec[0]*bsdf, ovec[1]*bsdf, ovec[2]*bsdf);
	    }
	if (fclose(fp) == EOF) {
		fprintf(stderr, "%s: error writing data to '%s'\n",
				progname, fname);
		return(0);
	}
	return(1);
}

/* Build BSDF values from scattering interpolant representation */
static int
build_wRBF(void)
{
	const char	*pref;
	int		i;

	if (input_orient > 0) {
		if (output_orient > 0)
			pref = frpref;
		else
			pref = ftpref;
	} else if (output_orient < 0)
		pref = brpref;
	else
		pref = btpref;

	for (i = 0; i < NINCIDENT; i++) {
		FVECT	ivec;
		RBFNODE	*rbf;
		get_ivector(ivec, i);
		if (input_orient < 0) {
			ivec[0] = -ivec[0]; ivec[1] = -ivec[1]; ivec[2] = -ivec[2];
		}
		rbf = advect_rbf(ivec, 15000);
		if (!plotRBF(tfile_name(pref, dsuffix, i), rbf))
			return(0);
		if (rbf) free(rbf);
	}
	return(1);				/* next call frees */
}

/* Put out mirror arrow for the given incident vector */
static void
put_mirror_arrow(const FVECT ivec, int inc_side)
{
	const double	arrow_len = 1.2*bsdf_rad;
	const double	tip_len = 0.2*bsdf_rad;
	FVECT		origin, refl;
	int		i;

	for (i = 3; i--; ) origin[i] = ivec[i]*sph_rad;
	origin[0] -= inc_side*sph_xoffset;

	refl[0] = 2.*ivec[2]*ivec[0];
	refl[1] = 2.*ivec[2]*ivec[1];
	refl[2] = 2.*ivec[2]*ivec[2] - 1.;

	printf("\n# Mirror arrow\n");
	printf("\narrow_mat cylinder inc_dir\n0\n0\n7");
	printf("\n\t%f %f %f\n\t%f %f %f\n\t%f\n",
			origin[0], origin[1], origin[2]+arrow_len,
			origin[0], origin[1], origin[2],
			arrow_rad);
	printf("\narrow_mat cylinder mir_dir\n0\n0\n7");
	printf("\n\t%f %f %f\n\t%f %f %f\n\t%f\n",
			origin[0], origin[1], origin[2],
			origin[0] + arrow_len*refl[0],
			origin[1] + arrow_len*refl[1],
			origin[2] + arrow_len*refl[2],
			arrow_rad);
	printf("\narrow_mat cone mir_tip\n0\n0\n8");
	printf("\n\t%f %f %f\n\t%f %f %f\n\t%f 0\n",
			origin[0] + (arrow_len-.5*tip_len)*refl[0],
			origin[1] + (arrow_len-.5*tip_len)*refl[1],
			origin[2] + (arrow_len-.5*tip_len)*refl[2],
			origin[0] + (arrow_len+.5*tip_len)*refl[0],
			origin[1] + (arrow_len+.5*tip_len)*refl[1],
			origin[2] + (arrow_len+.5*tip_len)*refl[2],
			2.*arrow_rad);
}

/* Put out transmitted direction arrow for the given incident vector */
static void
put_trans_arrow(const FVECT ivec, int inc_side)
{
	const double	arrow_len = 1.2*bsdf_rad;
	const double	tip_len = 0.2*bsdf_rad;
	FVECT		origin;
	int		i;

	for (i = 3; i--; ) origin[i] = ivec[i]*sph_rad;
	origin[0] -= inc_side*sph_xoffset;

	printf("\n# Transmission arrow\n");
	printf("\narrow_mat cylinder trans_dir\n0\n0\n7");
	printf("\n\t%f %f %f\n\t%f %f %f\n\t%f\n",
			origin[0], origin[1], origin[2],
			origin[0], origin[1], origin[2]-arrow_len,
			arrow_rad);
	printf("\narrow_mat cone trans_tip\n0\n0\n8");
	printf("\n\t%f %f %f\n\t%f %f %f\n\t%f 0\n",
			origin[0], origin[1], origin[2]-arrow_len+.5*tip_len,
			origin[0], origin[1], origin[2]-arrow_len-.5*tip_len,
			2.*arrow_rad);	
}

/* Compute rotation (x,y,z) => (xp,yp,zp) */
static int
addrot(char *xf, const FVECT xp, const FVECT yp, const FVECT zp)
{
	int	n = 0;
	double	theta;

	if (yp[2]*yp[2] + zp[2]*zp[2] < 2.*FTINY*FTINY) {
		/* Special case for X' along Z-axis */
		theta = -atan2(yp[0], yp[1]);
		sprintf(xf, " -ry %f -rz %f",
				xp[2] < 0.0 ? 90.0 : -90.0,
				theta*(180./PI));
		return(4);
	}
	theta = atan2(yp[2], zp[2]);
	if (!FEQ(theta,0.0)) {
		sprintf(xf, " -rx %f", theta*(180./PI));
		while (*xf) ++xf;
		n += 2;
	}
	theta = Asin(-xp[2]);
	if (!FEQ(theta,0.0)) {
		sprintf(xf, " -ry %f", theta*(180./PI));
		while (*xf) ++xf;
		n += 2;
	}
	theta = atan2(xp[1], xp[0]);
	if (!FEQ(theta,0.0)) {
		sprintf(xf, " -rz %f", theta*(180./PI));
		/* while (*xf) ++xf; */
		n += 2;
	}
	return(n);
}

/* Put out BSDF surfaces */
static int
put_BSDFs(void)
{
	const double	scalef = bsdf_rad/(log10(overall_max) - min_log10);
	FVECT		ivec;
	RREAL		vMtx[3][3];
	char		*fname;
	char		cmdbuf[256];
	char		xfargs[128];
	int		nxfa;
	int		i;

	printf("\n# Gensurf output corresponding to %d incident directions\n",
			NINCIDENT);

	printf("\nvoid glow arrow_glow\n0\n0\n4 1 0 1 0\n");
	printf("\nvoid mixfunc arrow_mat\n4 arrow_glow void .5 .\n0\n0\n");

	if (front_comp & SDsampR)
		for (i = 0; i < NINCIDENT; i++) {
			get_ivector(ivec, i);
			put_mirror_arrow(ivec, 1);
			sprintf(xfargs, "-s %f -t %f %f %f", bsdf_rad,
					ivec[0]*sph_rad - sph_xoffset,
					ivec[1]*sph_rad, ivec[2]*sph_rad);
			nxfa = 6;
			printf("\nvoid colorfunc scale_pat\n");
			printf("%d bsdf_red bsdf_grn bsdf_blu bsdf2rad.cal\n\t%s\n0\n0\n",
					4+nxfa, xfargs);
			printf("\nscale_pat glow scale_mat\n0\n0\n4 1 1 1 0\n");
			SDcompXform(vMtx, ivec, Yaxis);
			nxfa = addrot(xfargs, vMtx[0], vMtx[1], vMtx[2]);
			sprintf(xfargs+strlen(xfargs), " -s %f -t %f %f %f",
					scalef, ivec[0]*sph_rad - sph_xoffset,
					ivec[1]*sph_rad, ivec[2]*sph_rad);
			nxfa += 6;
			fname = tfile_name(frpref, dsuffix, i);
			sprintf(cmdbuf, "gensurf scale_mat %s%d %s %s %s %d %d | xform -mx -my %s",
					frpref, i+1, fname, fname, fname, SAMPRES-1, SAMPRES-1,
					xfargs);
			if (!run_cmd(cmdbuf))
				return(0);
		}
	if (front_comp & SDsampT)
		for (i = 0; i < NINCIDENT; i++) {
			get_ivector(ivec, i);
			put_trans_arrow(ivec, 1);
			sprintf(xfargs, "-s %f -t %f %f %f", bsdf_rad,
					ivec[0]*sph_rad - sph_xoffset,
					ivec[1]*sph_rad, ivec[2]*sph_rad);
			nxfa = 6;
			printf("\nvoid colorfunc scale_pat\n");
			printf("%d bsdf_red bsdf_grn bsdf_blu bsdf2rad.cal\n\t%s\n0\n0\n",
					4+nxfa, xfargs);
			printf("\nscale_pat glow scale_mat\n0\n0\n4 1 1 1 0\n");
			SDcompXform(vMtx, ivec, Yaxis);
			nxfa = addrot(xfargs, vMtx[0], vMtx[1], vMtx[2]);
			sprintf(xfargs+strlen(xfargs), " -s %f -t %f %f %f",
					scalef, ivec[0]*sph_rad - sph_xoffset,
					ivec[1]*sph_rad, ivec[2]*sph_rad);
			nxfa += 6;
			fname = tfile_name(ftpref, dsuffix, i);
			sprintf(cmdbuf, "gensurf scale_mat %s%d %s %s %s %d %d | xform -I -mx -my %s",
					ftpref, i+1, fname, fname, fname, SAMPRES-1, SAMPRES-1,
					xfargs);
			if (!run_cmd(cmdbuf))
				return(0);
		}
	if (back_comp & SDsampR)
		for (i = 0; i < NINCIDENT; i++) {
			get_ivector(ivec, i);
			put_mirror_arrow(ivec, -1);
			fname = tfile_name(brpref, dsuffix, i);
			sprintf(xfargs, "-s %f -t %f %f %f", bsdf_rad,
					ivec[0]*sph_rad + sph_xoffset,
					ivec[1]*sph_rad, ivec[2]*sph_rad);
			nxfa = 6;
			printf("\nvoid colorfunc scale_pat\n");
			printf("%d bsdf_red bsdf_grn bsdf_blu bsdf2rad.cal\n\t%s\n0\n0\n",
					4+nxfa, xfargs);
			printf("\nscale_pat glow scale_mat\n0\n0\n4 1 1 1 0\n");
			SDcompXform(vMtx, ivec, Yaxis);
			nxfa = addrot(xfargs, vMtx[0], vMtx[1], vMtx[2]);
			sprintf(xfargs+strlen(xfargs), " -s %f -t %f %f %f",
					scalef, ivec[0]*sph_rad + sph_xoffset,
					ivec[1]*sph_rad, ivec[2]*sph_rad);
			nxfa += 6;
			fname = tfile_name(brpref, dsuffix, i);
			sprintf(cmdbuf, "gensurf scale_mat %s%d %s %s %s %d %d | xform -I -ry 180 -mx -my %s",
					brpref, i+1, fname, fname, fname, SAMPRES-1, SAMPRES-1,
					xfargs);
			if (!run_cmd(cmdbuf))
				return(0);
		}
	if (back_comp & SDsampT)
		for (i = 0; i < NINCIDENT; i++) {
			get_ivector(ivec, i);
			put_trans_arrow(ivec, -1);
			fname = tfile_name(btpref, dsuffix, i);
			sprintf(xfargs, "-s %f -t %f %f %f", bsdf_rad,
					ivec[0]*sph_rad + sph_xoffset,
					ivec[1]*sph_rad, ivec[2]*sph_rad);
			nxfa = 6;
			printf("\nvoid colorfunc scale_pat\n");
			printf("%d bsdf_red bsdf_grn bsdf_blu bsdf2rad.cal\n\t%s\n0\n0\n",
					4+nxfa, xfargs);
			printf("\nscale_pat glow scale_mat\n0\n0\n4 1 1 1 0\n");
			SDcompXform(vMtx, ivec, Yaxis);
			nxfa = addrot(xfargs, vMtx[0], vMtx[1], vMtx[2]);
			sprintf(xfargs+strlen(xfargs), " -s %f -t %f %f %f",
					scalef, ivec[0]*sph_rad + sph_xoffset,
					ivec[1]*sph_rad, ivec[2]*sph_rad);
			nxfa += 6;
			fname = tfile_name(btpref, dsuffix, i);
			sprintf(cmdbuf, "gensurf scale_mat %s%d %s %s %s %d %d | xform -ry 180 -mx -my %s",
					btpref, i+1, fname, fname, fname, SAMPRES-1, SAMPRES-1,
					xfargs);
			if (!run_cmd(cmdbuf))
				return(0);
		}
	return(1);
}

/* Put our hemisphere material */
static void
put_matBSDF(const char *XMLfile)
{
	const char	*curdir = "./";

	if (!XMLfile) {			/* simple material */
		printf("\n# Simplified material because we have no XML input\n");
		printf("\nvoid brightfunc latlong\n2 latlong bsdf2rad.cal\n0\n0\n");
		if ((front_comp|back_comp) & SDsampT)
			printf("\nlatlong trans %s\n0\n0\n7 .75 .75 .75 0 0 .5 .8\n",
					sph_mat);
		else
			printf("\nlatlong plastic %s\n0\n0\n5 .5 .5 .5 0 0\n",
					sph_mat);
		return;
	}
	switch (XMLfile[0]) {		/* avoid RAYPATH search */
	case '.':
	CASEDIRSEP:
		curdir = "";
		break;
	case '\0':
		fprintf(stderr, "%s: empty file name in put_matBSDF\n", progname);
		exit(1);
		break;
	}
	printf("\n# Actual BSDF material for rendering the hemispheres\n");
	printf("\nvoid BSDF BSDFmat\n6 0 \"%s%s\" 0 1 0 .\n0\n0\n",
			curdir, XMLfile);
	printf("\nvoid plastic black\n0\n0\n5 0 0 0 0 0\n");
	printf("\nvoid mixfunc %s\n4 BSDFmat black latlong bsdf2rad.cal\n0\n0\n",
			sph_mat);
}

/* Put out overhead parallel light source */
static void
put_source(void)
{
	printf("\n# Overhead parallel light source\n");
	printf("\nvoid light bright\n0\n0\n3 2000 2000 2000\n");
	printf("\nbright source light\n0\n0\n4 0 0 1 2\n");
	printf("\n# Material used for labels\n");
	printf("\nvoid trans vellum\n0\n0\n7 1 1 1 0 0 .5 0\n");
}

/* Put out hemisphere(s) */
static void
put_hemispheres(void)
{
	printf("\n# Hemisphere(s) for showing BSDF appearance (if XML file)\n");
	printf("\nvoid antimatter anti_sph\n2 void %s\n0\n0\n", sph_mat);
	if (front_comp) {
		printf("\n%s sphere Front\n0\n0\n4 %f 0 0 %f\n",
				sph_mat, -sph_xoffset, sph_rad);
		printf("\n!genbox anti_sph sph_eraser %f %f %f | xform -t %f %f %f\n",
				2.02*sph_rad, 2.02*sph_rad, 1.02*sph_rad,
				-1.01*sph_rad - sph_xoffset, -1.01*sph_rad, -1.01*sph_rad);
		printf("\nvoid brighttext front_text\n3 helvet.fnt . FRONT\n0\n");
		printf("12\n\t%f %f 0\n\t%f 0 0\n\t0 %f 0\n\t.01 1 -.1\n",
				-.22*sph_rad - sph_xoffset, -1.4*sph_rad,
				.35/5.*sph_rad, -1.6*.35/5.*sph_rad);
		printf("\nfront_text alias front_label_mat vellum\n");
		printf("\nfront_label_mat polygon front_label\n0\n0\n12");
		printf("\n\t%f %f 0\n\t%f %f 0\n\t%f %f 0\n\t%f %f 0\n",
				-.25*sph_rad - sph_xoffset, -1.3*sph_rad,
				-.25*sph_rad - sph_xoffset, (-1.4-1.6*.35/5.-.1)*sph_rad,
				.25*sph_rad - sph_xoffset, (-1.4-1.6*.35/5.-.1)*sph_rad,
				.25*sph_rad - sph_xoffset, -1.3*sph_rad );
	}
	if (back_comp) {
		printf("\n%s bubble Back\n0\n0\n4 %f 0 0 %f\n",
				sph_mat, sph_xoffset, sph_rad);
		printf("\n!genbox anti_sph sph_eraser %f %f %f | xform -t %f %f %f\n",
				2.02*sph_rad, 2.02*sph_rad, 1.02*sph_rad,
				-1.01*sph_rad + sph_xoffset, -1.01*sph_rad, -1.01*sph_rad);
		printf("\nvoid brighttext back_text\n3 helvet.fnt . BACK\n0\n");
		printf("12\n\t%f %f 0\n\t%f 0 0\n\t0 %f 0\n\t.01 1 -.1\n",
				-.22*sph_rad + sph_xoffset, -1.4*sph_rad,
				.35/4.*sph_rad, -1.6*.35/4.*sph_rad);
		printf("\nback_text alias back_label_mat vellum\n");
		printf("\nback_label_mat polygon back_label\n0\n0\n12");
		printf("\n\t%f %f 0\n\t%f %f 0\n\t%f %f 0\n\t%f %f 0\n",
				-.25*sph_rad + sph_xoffset, -1.3*sph_rad,
				-.25*sph_rad + sph_xoffset, (-1.4-1.6*.35/4.-.1)*sph_rad,
				.25*sph_rad + sph_xoffset, (-1.4-1.6*.35/4.-.1)*sph_rad,
				.25*sph_rad + sph_xoffset, -1.3*sph_rad );
	}
}

/* Put out falsecolor scale and name label */
static void
put_scale(void)
{
	const double	max_log10 = log10(overall_max);
	const double	leg_width = 2.*.75*(sph_xoffset - sph_rad);
	const double	leg_height = 2.*sph_rad;
	const int	text_lines = 6;
	const int	text_digits = 8;
	char		fmt[16];
	int		i;

	printf("\n# BSDF legend with falsecolor scale\n");
	printf("\nvoid colorfunc lscale\n10 sca_red(Py) sca_grn(Py) sca_blu(Py)");
	printf("\n\tbsdf2rad.cal -s %f -t 0 %f 0\n0\n0\n", leg_height, -.5*leg_height);
	sprintf(fmt, "%%.%df", text_digits-3);
	for (i = 0; i < text_lines; i++) {
		char	vbuf[16];
		sprintf(vbuf, fmt, pow(10., (i+.5)/text_lines*(max_log10-min_log10)+min_log10));
		printf("\nlscale brighttext lscale\n");
		printf("3 helvet.fnt . %s\n0\n12\n", vbuf);
		printf("\t%f %f 0\n", -.45*leg_width, ((i+.9)/text_lines-.5)*leg_height);
		printf("\t%f 0 0\n", .8*leg_width/strlen(vbuf));
		printf("\t0 %f 0\n", -.9/text_lines*leg_height);
		printf("\t.01 1 -.1\n");
	}
	printf("\nlscale alias legend_mat vellum\n");
	printf("\nlegend_mat polygon legend\n0\n0\n12");
	printf("\n\t%f %f 0\n\t%f %f 0\n\t%f %f 0\n\t%f %f 0\n",
			-.5*leg_width, .5*leg_height,
			-.5*leg_width, -.5*leg_height,
			.5*leg_width, -.5*leg_height,
			.5*leg_width, .5*leg_height);
	printf("\nvoid brighttext BSDFtitle\n3 helvet.fnt . BSDF\n0\n12\n");
	printf("\t%f %f 0\n", -.25*leg_width, .7*leg_height);
	printf("\t%f 0 0\n", .4/4.*leg_width);
	printf("\t0 %f 0\n", -.1*leg_height);
	printf("\t.01 1 -.1\n");
	printf("\nBSDFtitle alias title_mat vellum\n");
	printf("\ntitle_mat polygon title\n0\n0\n12");
	printf("\n\t%f %f 0\n\t%f %f 0\n\t%f %f 0\n\t%f %f 0\n",
			-.3*leg_width, .75*leg_height,
			-.3*leg_width, .55*leg_height,
			.3*leg_width, .55*leg_height,
			.3*leg_width, .75*leg_height);
	if (!bsdf_name[0])
		return;
	printf("\nvoid brighttext BSDFname\n3 helvet.fnt . \"%s\"\n0\n12\n", bsdf_name);
	printf("\t%f %f 0\n", -.95*leg_width, -.6*leg_height);
	printf("\t%f 0 0\n", 1.8/strlen(bsdf_name)*leg_width);
	printf("\t0 %f 0\n", -.1*leg_height);
	printf("\t.01 1 -.1\n");
	printf("\nBSDFname alias name_mat vellum\n");
	printf("\nname_mat polygon name\n0\n0\n12");
	printf("\n\t%f %f 0\n\t%f %f 0\n\t%f %f 0\n\t%f %f 0\n",
			-leg_width, -.55*leg_height,
			-leg_width, -.75*leg_height,
			leg_width, -.75*leg_height,
			leg_width, -.55*leg_height);
}

/* Convert MGF to Radiance in output */
static void
convert_mgf(const char *mgfdata)
{
	int	len = strlen(mgfdata);
	char	mgfn[128];
	char	radfn[128];
	char	cmdbuf[256];
	float	xmin, xmax, ymin, ymax, zmin, zmax;
	double	max_dim;
	int	fd;
	FILE	*fp;

	if (!len) return;
	strcpy(mgfn, tfile_name("geom", ".mgf", 0));
	fd = open(mgfn, O_WRONLY|O_CREAT, 0666);
	if (fd < 0 || write(fd, mgfdata, len) != len) {
		fprintf(stderr, "%s: cannot write file '%s'\n",
				progname, mgfn);
		return;
	}
	close(fd);
	strcpy(radfn, tfile_name("geom", ".rad", 0));
	sprintf(cmdbuf, "mgf2rad %s > %s", mgfn, radfn);
	if (!run_cmd(cmdbuf))
		return;
	sprintf(cmdbuf, "getbbox -w -h %s", radfn);
	if ((fp = popen(cmdbuf, "r")) == NULL ||
			fscanf(fp, "%f %f %f %f %f %f",
				&xmin, &xmax, &ymin, &ymax, &zmin, &zmax) != 6
			|| pclose(fp) < 0) {
		fprintf(stderr, "%s: error reading from command: %s\n",
				progname, cmdbuf);
		return;
	}
	max_dim = ymax - ymin;
	if (xmax - xmin > max_dim)
		max_dim = xmax - xmin;
	if (front_comp) {
		printf("\n# BSDF system geometry (front view)\n");
		sprintf(cmdbuf, "xform -t %f %f %f -s %f -t %f %f 0 %s",
				-.5*(xmin+xmax), -.5*(ymin+ymax), -zmax,
				1.5*sph_rad/max_dim,
				-sph_xoffset, -2.5*sph_rad,
				radfn);
		if (!run_cmd(cmdbuf))
			return;
	}
	if (back_comp) {
		printf("\n# BSDF system geometry (back view)\n");
		sprintf(cmdbuf, "xform -t %f %f %f -s %f -ry 180 -t %f %f 0 %s",
				-.5*(xmin+xmax), -.5*(ymin+ymax), -zmin,
				1.5*sph_rad/max_dim,
				sph_xoffset, -2.5*sph_rad,
				radfn);
		if (!run_cmd(cmdbuf))
			return;
	}
}

/* Check RBF input header line & get minimum BSDF value */
static int
rbf_headline(char *s, void *p)
{
	char	fmt[64];

	if (formatval(fmt, s)) {
		if (strcmp(fmt, BSDFREP_FMT))
			return(-1);
		return(0);
	}
	if (!strncmp(s, "IO_SIDES=", 9)) {
		sscanf(s+9, "%d %d", &input_orient, &output_orient);
		if (input_orient == output_orient) {
			if (input_orient > 0)
				front_comp |= SDsampR;
			else
				back_comp |= SDsampR;
		} else if (input_orient > 0)
			front_comp |= SDsampT;
		else
			back_comp |= SDsampT;
		return(0);
	}
	if (!strncmp(s, "BSDFMIN=", 8)) {
		sscanf(s+8, "%lf", &bsdf_min);
		if (bsdf_min < overall_min)
			overall_min = bsdf_min;
		return(0);
	}
	return(0);
}

/* Produce a Radiance model plotting the given BSDF representation */
int
main(int argc, char *argv[])
{
	int	inpXML = -1;
	SDData	myBSDF;
	int	n;
						/* check arguments */
	progname = argv[0];
	if (argc > 1 && (n = strlen(argv[1])-4) > 0) {
		if (!strcasecmp(argv[1]+n, ".xml"))
			inpXML = 1;
		else if (!strcasecmp(argv[1]+n, ".sir"))
			inpXML = 0;
	}
	if (inpXML < 0 || inpXML & (argc > 2)) {
		fprintf(stderr, "Usage: %s bsdf.xml > output.rad\n", progname);
		fprintf(stderr, "   Or: %s hemi1.sir hemi2.sir .. > output.rad\n", progname);
		return(1);
	}
	fputs("# ", stdout);			/* copy our command */
	printargs(argc, argv, stdout);
						/* evaluate BSDF */
	if (inpXML) {
		SDclearBSDF(&myBSDF, argv[1]);
		if (SDreportError(SDloadFile(&myBSDF, argv[1]), stderr))
			return(1);
		if (myBSDF.rf != NULL) front_comp |= SDsampR;
		if (myBSDF.tf != NULL) front_comp |= SDsampT;
		if (myBSDF.rb != NULL) back_comp |= SDsampR;
		if (myBSDF.tb != NULL) back_comp |= SDsampT;
		if (!front_comp & !back_comp) {
			fprintf(stderr, "%s: nothing to plot in '%s'\n",
					progname, argv[1]);
			return(1);
		}
		if (front_comp & SDsampR && myBSDF.rLambFront.cieY < overall_min*PI)
			overall_min = myBSDF.rLambFront.cieY/PI;
		if (back_comp & SDsampR && myBSDF.rLambBack.cieY < overall_min*PI)
			overall_min = myBSDF.rLambBack.cieY/PI;
		if ((front_comp|back_comp) & SDsampT &&
				myBSDF.tLamb.cieY < overall_min*PI)
			overall_min = myBSDF.tLamb.cieY/PI;
		set_minlog();
		if (!build_wBSDF(&myBSDF))
			return(1);
		if (myBSDF.matn[0])
			strcpy(bsdf_name, myBSDF.matn);
		else
			strcpy(bsdf_name, myBSDF.name);
		strcpy(bsdf_manuf, myBSDF.makr);
		put_matBSDF(argv[1]);
	} else {
		FILE	*fp;
		for (n = 1; n < argc; n++) {
			fp = fopen(argv[n], "rb");
			if (fp == NULL) {
				fprintf(stderr, "%s: cannot open BSDF interpolant '%s'\n",
						progname, argv[n]);
				return(1);
			}
			if (getheader(fp, rbf_headline, NULL) < 0) {
				fprintf(stderr, "%s: bad BSDF interpolant '%s'\n",
						progname, argv[n]);
				return(1);
			}
			fclose(fp);
		}
		set_minlog();
		for (n = 1; n < argc; n++) {
			fp = fopen(argv[n], "rb");
			if (!load_bsdf_rep(fp))
				return(1);
			fclose(fp);
			if (!build_wRBF())
				return(1);
		}
		put_matBSDF(NULL);
	}
	put_source();			/* before hemispheres & labels */
	put_hemispheres();
	put_scale();
	if (inpXML && myBSDF.mgf)
		convert_mgf(myBSDF.mgf);
	if (!put_BSDFs())
		return(1);
	cleanup_tmp();
	return(0);
}
