#ifndef lint
static const char	RCSid[] = "$Id: tmapcolrs.c,v 3.30 2011/06/28 21:11:04 greg Exp $";
#endif
/*
 * Routines for tone mapping on Radiance RGBE and XYZE pictures.
 *
 * Externals declared in tonemap.h
 */

#include "copyright.h"

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	<time.h>

#ifdef PCOND
#include	"rtprocess.h"
#endif
#include	"tmprivat.h"
#include	"resolu.h"

#define GAMTSZ	4096

typedef struct {
	uby8		gamb[GAMTSZ];	/* gamma lookup table */
	int		clfb[3];	/* encoded tm->clf */
	int32		cmatb[3][3];	/* encoded color transform */
	TMbright	inpsfb;		/* encoded tm->inpsf */
} COLRDATA;

static MEM_PTR	colrInit(TMstruct *);
static void	colrNewSpace(TMstruct *);
static gethfunc headline;

static struct tmPackage	colrPkg = {	/* our package functions */
	colrInit, colrNewSpace, free
};
static int	colrReg = -1;		/* our package registration number */

static TMbright	logi[256];


int
tmCvColrs(				/* convert RGBE/XYZE colors */
TMstruct	*tms,
TMbright	*ls,
uby8	*cs,
COLR	*scan,
int	len
)
{
	static const char	funcName[] = "tmCvColrs";
	int	cmon[4];
	COLRDATA	*cd;
	int	i, j, bi;
	int32	li, vl;

	if (tms == NULL)
		returnErr(TM_E_TMINVAL);
	if ((ls == NULL) | (scan == NULL) | (len < 0))
		returnErr(TM_E_ILLEGAL);
	if (colrReg < 0) {			/* build tables if necessary */
		colrReg = tmRegPkg(&colrPkg);
		if (colrReg < 0)
			returnErr(TM_E_CODERR1);
		for (i = 256; i--; )
			logi[i] = TM_BRTSCALE*log((i+.5)/256.) - .5;
		tmMkMesofact();
	}
	if ((cd = (COLRDATA *)tmPkgData(tms,colrReg)) == NULL)
		returnErr(TM_E_NOMEM);
	for (i = len; i--; ) {
		if (tmNeedMatrix(tms)) {		/* apply color xform */
			for (j = 3; j--; ) {
				vl =	cd->cmatb[j][RED]*(int32)scan[i][RED] +
					cd->cmatb[j][GRN]*(int32)scan[i][GRN] +
					cd->cmatb[j][BLU]*(int32)scan[i][BLU] ;
				if (vl < 0) cmon[j] = vl/0x10000;
				else cmon[j] = vl>>16;
			}
			cmon[EXP] = scan[i][EXP];
		} else
			copycolr(cmon, scan[i]);
							/* world luminance */
		li =	cd->clfb[RED]*(int32)cmon[RED] +
			cd->clfb[GRN]*(int32)cmon[GRN] +
			cd->clfb[BLU]*(int32)cmon[BLU] ;
		if (li >= 1L<<(12+8)) li = 255;
		else li >>= 12;
		bi = BRT2SCALE(cmon[EXP]-COLXS) + cd->inpsfb;
		if (li > 0)
			bi += logi[li];
		else {
			bi += logi[0];
			li = 1;				/* avoid /0 */
		}
		ls[i] = bi;
		if (cs == TM_NOCHROM)			/* no color? */
			continue;
							/* mesopic adj. */
		if (tms->flags & TM_F_MESOPIC && bi < BMESUPPER) {
			int	pf, sli = normscot(cmon);
			if (bi < BMESLOWER) {
				cmon[RED] = cmon[GRN] = cmon[BLU] = sli;
			} else {
				if (tms->flags & TM_F_BW)
					cmon[RED] = cmon[GRN] = cmon[BLU] = li;
				pf = tmMesofact[bi-BMESLOWER];
				sli *= 256 - pf;
				for (j = 3; j--; ) {
					cmon[j] = sli + pf*cmon[j];
					if (cmon[j] <= 0) cmon[j] = 0;
					else cmon[j] >>= 8;
				}
			}
		} else if (tms->flags & TM_F_BW) {
			cmon[RED] = cmon[GRN] = cmon[BLU] = li;
		} else {
			for (j = 3; j--; )
				if (cmon[j] < 0) cmon[j] = 0;
		}
		bi = ( (uint32)GAMTSZ*cd->clfb[RED]*cmon[RED]/li ) >> 12;
		cs[3*i  ] = bi>=GAMTSZ ? 255 : cd->gamb[bi];
		bi = ( (uint32)GAMTSZ*cd->clfb[GRN]*cmon[GRN]/li ) >> 12;
		cs[3*i+1] = bi>=GAMTSZ ? 255 : cd->gamb[bi];
		bi = ( (uint32)GAMTSZ*cd->clfb[BLU]*cmon[BLU]/li ) >> 12;
		cs[3*i+2] = bi>=GAMTSZ ? 255 : cd->gamb[bi];
	}
	returnOK;
}


#define	FMTRGB		1	/* Input is RGBE */
#define	FMTCIE		2	/* Input is CIE XYZE */
#define FMTUNK		3	/* Input format is unspecified */
#define FMTBAD		4	/* Input is not a recognized format */

static struct radhead {
	int	format;		/* FMTRGB, FMTCIE, FMTUNK, FMTBAD */
	double	expos;		/* input exposure value */
	RGBPRIMP	primp;	/* input primaries */
	RGBPRIMS	mypri;	/* custom primaries */
} rhdefault = {FMTUNK, 1., stdprims, STDPRIMS};


static int
headline(			/* grok a header line */
	char	*s,
	void	*vrh
)
{
	char	fmt[32];
	struct radhead	*rh = vrh;

	if (formatval(fmt, s)) {
		if (!strcmp(fmt, COLRFMT))
			rh->format = FMTRGB;
		else if (!strcmp(fmt, CIEFMT))
			rh->format = FMTCIE;
		else
			rh->format = FMTBAD;
		return(0);
	}
	if (isexpos(s)) {
		rh->expos *= exposval(s);
		return(0);
	}
	if (isprims(s)) {
		primsval(rh->mypri, s);
		rh->primp = rh->mypri;
		return(0);
	}
	return(0);
}


int
tmLoadPicture(				/* convert Radiance picture */
TMstruct	*tms,
TMbright	**lpp,
uby8	**cpp,
int	*xp,
int	*yp,
char	*fname,
FILE	*fp
)
{
	char	*funcName = fname==NULL ? "tmLoadPicture" : fname;
	FILE	*inpf;
	struct radhead	info;
	int	err;
	COLR	*scanin = NULL;
	int	i;
						/* check arguments */
	if (tms == NULL)
		returnErr(TM_E_TMINVAL);
	if ((lpp == NULL) | (xp == NULL) | (yp == NULL) |
			((fname == NULL) & (fp == TM_GETFILE)))
		returnErr(TM_E_ILLEGAL);
	*xp = *yp = 0;				/* error precaution */
	if ((inpf = fp) == TM_GETFILE && (inpf = fopen(fname, "rb")) == NULL)
		returnErr(TM_E_BADFILE);
	*lpp = NULL;
	if (cpp != TM_NOCHROMP) *cpp = NULL;
	info = rhdefault;			/* get our header */
	getheader(inpf, headline, &info);
	if ((info.format == FMTBAD) | (info.expos <= 0.) ||
			fgetresolu(xp, yp, inpf) < 0) {
		err = TM_E_BADFILE; goto done;
	}
	if (info.format == FMTUNK)		/* assume RGBE format */
		info.format = FMTRGB;
	if (info.format == FMTRGB)
		info.expos /= WHTEFFICACY;
	else if (info.format == FMTCIE)
		info.primp = TM_XYZPRIM;
						/* prepare library */
	if ((err = tmSetSpace(tms, info.primp, 1./info.expos, NULL)) != TM_E_OK)
		goto done;
	err = TM_E_NOMEM;			/* allocate arrays */
	*lpp = (TMbright *)malloc(sizeof(TMbright) * *xp * *yp);
	if (*lpp == NULL)
		goto done;
	if (cpp != TM_NOCHROMP) {
		*cpp = (uby8 *)malloc(3*sizeof(uby8) * *xp * *yp);
		if (*cpp == NULL)
			goto done;
	}
	scanin = (COLR *)malloc(sizeof(COLR) * *xp);
	if (scanin == NULL)
		goto done;
	err = TM_E_BADFILE;			/* read & convert scanlines */
	for (i = 0; i < *yp; i++) {
		if (freadcolrs(scanin, *xp, inpf) < 0) {
			err = TM_E_BADFILE; break;
		}
		err = tmCvColrs(tms, *lpp + (i * *xp),
			cpp==TM_NOCHROMP ? TM_NOCHROM : *cpp + (i * 3 * *xp),
				scanin, *xp);
		if (err != TM_E_OK)
			break;
	}
done:						/* clean up */
	if (fp == NULL)
		fclose(inpf);
	if (scanin != NULL)
		free((MEM_PTR)scanin);
	if (err != TM_E_OK) {
		if (*lpp != NULL)
			free((MEM_PTR)*lpp);
		if (cpp != TM_NOCHROMP && *cpp != NULL)
			free((MEM_PTR)*cpp);
		returnErr(err);
	}
	returnOK;
}


#ifdef PCOND
static int					/* run pcond to map picture */
dopcond(psp, xp, yp, flags, monpri, gamval, Lddyn, Ldmax, fname)
uby8	**psp;
int	*xp, *yp;
int	flags;
RGBPRIMP	monpri;
double	gamval, Lddyn, Ldmax;
char	*fname;
{
	char	*funcName = fname;
	TMstruct	*tms = NULL;
	char	cmdbuf[1024];
	FILE	*infp;
	COLR	*scan;
	uby8	*rp;
	int	y;
	int	x;
					/* set up gamma correction */
	if (setcolrcor(pow, 1./gamval) < 0)
		returnErr(TM_E_NOMEM);
					/* create command */
	strcpy(cmdbuf, PCOND);
	if (flags & TM_F_HCONTR)
		strcat(cmdbuf, " -s");
	if (flags & TM_F_MESOPIC)
		strcat(cmdbuf, " -c");
	if (flags & TM_F_LINEAR)
		strcat(cmdbuf, " -l");
	if (flags & TM_F_ACUITY)
		strcat(cmdbuf, " -a");
	if (flags & TM_F_VEIL)
		strcat(cmdbuf, " -v");
	if (flags & TM_F_CWEIGHT)
		strcat(cmdbuf, " -w");
	if (monpri != stdprims)
		sprintf(cmdbuf+strlen(cmdbuf), " -p %f %f %f %f %f %f %f %f",
				monpri[RED][CIEX], monpri[RED][CIEY],
				monpri[GRN][CIEX], monpri[GRN][CIEY],
				monpri[BLU][CIEX], monpri[BLU][CIEY],
				monpri[WHT][CIEX], monpri[WHT][CIEY]);
	sprintf(cmdbuf+strlen(cmdbuf), " -d %f -u %f %s", Lddyn, Ldmax, fname);
					/* start pcond */
	if ((infp = popen(cmdbuf, "r")) == NULL)
		returnErr(TM_E_BADFILE);
					/* check picture format and size */
	if (checkheader(infp, COLRFMT, NULL) < 0 ||
			fgetresolu(xp, yp, infp) < 0) {
		pclose(infp);
		returnErr(TM_E_BADFILE);
	}
					/* allocate arrays */
	scan = (COLR *)malloc(sizeof(COLR) * *xp);
	if (flags & TM_F_BW)
		rp = (uby8 *)malloc(sizeof(uby8) * *xp * *yp);
	else
		rp = (uby8 *)malloc(3*sizeof(uby8) * *xp * *yp);
	if (((*psp = rp) == NULL) | (scan == NULL)) {
		pclose(infp);
		returnErr(TM_E_NOMEM);
	}
					/* read and gamma map file */
	for (y = 0; y < *yp; y++) {
		if (freadcolrs(scan, *xp, infp) < 0) {
			pclose(infp);
			free((MEM_PTR)scan);
			free((MEM_PTR)*psp);
			*psp = NULL;
			returnErr(TM_E_BADFILE);
		}
		colrs_gambs(scan, *xp);
		if (flags & TM_F_BW)
			for (x = 0; x < *xp; x++)
				*rp++ = normbright(scan[x]);
		else
			for (x = 0; x < *xp; x++) {
				*rp++ = scan[x][RED];
				*rp++ = scan[x][GRN];
				*rp++ = scan[x][BLU];
			}
	}
	free((MEM_PTR)scan);
	pclose(infp);
	returnOK;
}
#endif


int					/* map a Radiance picture */
tmMapPicture(psp, xp, yp, flags, monpri, gamval, Lddyn, Ldmax, fname, fp)
uby8	**psp;
int	*xp, *yp;
int	flags;
RGBPRIMP	monpri;
double	gamval, Lddyn, Ldmax;
char	*fname;
FILE	*fp;
{
	char	*funcName = fname==NULL ? "tmMapPicture" : fname;
	TMstruct	*tms = NULL;
	uby8	*cp;
	TMbright	*lp;
	int	err;
						/* check arguments */
	if ((psp == NULL) | (xp == NULL) | (yp == NULL) | (monpri == NULL) |
			((fname == NULL) & (fp == TM_GETFILE)))
		returnErr(TM_E_ILLEGAL);
						/* set defaults */
	if (gamval < MINGAM) gamval = DEFGAM;
	if (Lddyn < MINLDDYN) Lddyn = DEFLDDYN;
	if (Ldmax < MINLDMAX) Ldmax = DEFLDMAX;
	if (flags & TM_F_BW) monpri = stdprims;
#ifdef PCOND
						/* check for pcond run */
	if (fp == TM_GETFILE && flags & TM_F_UNIMPL)
		return( dopcond(psp, xp, yp, flags,
				monpri, gamval, Lddyn, Ldmax, fname) );
#endif
						/* initialize tone mapping */
	if ((tms = tmInit(flags, monpri, gamval)) == NULL)
		returnErr(TM_E_NOMEM);
						/* load & convert picture */
	err = tmLoadPicture(tms, &lp, (flags&TM_F_BW) ? TM_NOCHROMP : &cp,
			xp, yp, fname, fp);
	if (err != TM_E_OK) {
		tmDone(tms);
		return(err);
	}
						/* allocate space for result */
	if (flags & TM_F_BW) {
		*psp = (uby8 *)malloc(sizeof(uby8) * *xp * *yp);
		if (*psp == NULL) {
			free((MEM_PTR)lp);
			tmDone(tms);
			returnErr(TM_E_NOMEM);
		}
		cp = TM_NOCHROM;
	} else
		*psp = cp;
						/* compute color mapping */
	err = tmAddHisto(tms, lp, *xp * *yp, 1);
	if (err != TM_E_OK)
		goto done;
	err = tmComputeMapping(tms, gamval, Lddyn, Ldmax);
	if (err != TM_E_OK)
		goto done;
						/* map colors */
	err = tmMapPixels(tms, *psp, lp, cp, *xp * *yp);

done:						/* clean up */
	free((MEM_PTR)lp);
	tmDone(tms);
	if (err != TM_E_OK) {			/* free memory on error */
		free((MEM_PTR)*psp);
		*psp = NULL;
		returnErr(err);
	}
	returnOK;
}


static void
colrNewSpace(tms)		/* color space changed for tone mapping */
TMstruct	*tms;
{
	COLRDATA	*cd;
	double	d;
	int	i, j;

	cd = (COLRDATA *)tms->pd[colrReg];
	for (i = 3; i--; )
		cd->clfb[i] = 0x1000*tms->clf[i] + .5;
	cd->inpsfb = tmCvLuminance(tms->inpsf);
	for (i = 3; i--; )
		for (j = 3; j--; ) {
			d = tms->cmat[i][j] / tms->inpsf;
			cd->cmatb[i][j] = 0x10000*d + (d<0. ? -.5 : .5);
		}
}


static MEM_PTR
colrInit(tms)			/* initialize private data for tone mapping */
TMstruct	*tms;
{
	COLRDATA	*cd;
	int	i;
					/* allocate our data */
	cd = (COLRDATA *)malloc(sizeof(COLRDATA));
	if (cd == NULL)
		return(NULL);
	tms->pd[colrReg] = (MEM_PTR)cd;
					/* compute gamma table */
	for (i = GAMTSZ; i--; )
		cd->gamb[i] = 256.*pow((i+.5)/GAMTSZ, 1./tms->mongam);
					/* compute color and scale factors */
	colrNewSpace(tms);
	return((MEM_PTR)cd);
}
