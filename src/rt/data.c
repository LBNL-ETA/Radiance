#ifndef lint
static const char	RCSid[] = "$Id: data.c,v 2.37 2024/03/12 16:54:51 greg Exp $";
#endif
/*
 *  data.c - routines dealing with interpolated data.
 */

#include "copyright.h"

#include  <time.h>

#include  "platform.h"
#include  "paths.h"
#include  "standard.h"
#include  "color.h"
#include  "resolu.h"
#include  "view.h"
#include  "data.h"

				/* picture memory usage before warning */
#ifndef PSIZWARN
#ifdef SMLMEM
#define PSIZWARN	3000000
#else
#define PSIZWARN	50000000
#endif
#endif

#ifndef TABSIZ
#define TABSIZ		997		/* table size (prime) */
#endif

#define hash(s)		(shash(s)%TABSIZ)


static DATARRAY	 *dtab[TABSIZ];		/* data array list */

static gethfunc headaspect;


DATARRAY *
getdata(				/* get data array dname */
	char  *dname
)
{
	char  *dfname;
	FILE  *fp;
	int  asize=0;
	int  i, j;
	DATARRAY  *dp;
						/* look for array in list */
	for (dp = dtab[hash(dname)]; dp != NULL; dp = dp->next)
		if (!strcmp(dname, dp->name))
			return(dp);		/* found! */
	/*
	 *	If we haven't loaded the data already, we will look
	 *  for it in the directories specified by the library path.
	 *
	 *	The file has the following format:
	 *
	 *		N
	 *		beg0	end0	n0
	 *		beg1	end1	n1
	 *		. . .
	 *		begN	endN	nN
	 *		data, later dimensions changing faster
	 *		. . .
	 *
	 *	For irregularly spaced points, the following can be
	 *  substituted for begi endi ni:
	 *
	 *		0 0 ni p0i p1i .. pni
	 */

	if ((dfname = getpath(dname, getrlibpath(), R_OK)) == NULL) {
		sprintf(errmsg, "cannot find data file \"%s\"", dname);
		error(SYSTEM, errmsg);
	}
	if ((fp = fopen(dfname, "r")) == NULL) {
		sprintf(errmsg, "cannot open data file \"%s\"", dfname);
		error(SYSTEM, errmsg);
	}
							/* get dimensions */
	if (fgetval(fp, 'i', &asize) <= 0)
		goto scanerr;
	if ((asize <= 0) | (asize > MAXDDIM)) {
		sprintf(errmsg, "bad number of dimensions for \"%s\"", dname);
		error(USER, errmsg);
	}
	if ((dp = (DATARRAY *)malloc(sizeof(DATARRAY))) == NULL)
		goto memerr;
	dp->name = savestr(dname);
	dp->type = DATATY;
	dp->nd = asize;
	asize = 1;
	for (i = 0; i < dp->nd; i++) {
		if (fgetval(fp, DATATY, &dp->dim[i].org) <= 0)
			goto scanerr;
		if (fgetval(fp, DATATY, &dp->dim[i].siz) <= 0)
			goto scanerr;
		if (fgetval(fp, 'i', &dp->dim[i].ne) <= 0)
			goto scanerr;
		if (dp->dim[i].ne < 2)
			goto scanerr;
		asize *= dp->dim[i].ne;
		if ((dp->dim[i].siz -= dp->dim[i].org) == 0) {
			dp->dim[i].p = (DATATYPE *)
					malloc(dp->dim[i].ne*sizeof(DATATYPE));
			if (dp->dim[i].p == NULL)
				goto memerr;
			for (j = 0; j < dp->dim[i].ne; j++)
				if (fgetval(fp, DATATY, &dp->dim[i].p[j]) <= 0)
					goto scanerr;
			for (j = 1; j < dp->dim[i].ne-1; j++)
				if ((dp->dim[i].p[j-1] < dp->dim[i].p[j]) !=
					(dp->dim[i].p[j] < dp->dim[i].p[j+1]))
					goto scanerr;
			dp->dim[i].org = dp->dim[i].p[0];
			dp->dim[i].siz = dp->dim[i].p[dp->dim[i].ne-1]
						- dp->dim[i].p[0];
		} else
			dp->dim[i].p = NULL;
	}
	if ((dp->arr.d = (DATATYPE *)malloc(asize*sizeof(DATATYPE))) == NULL)
		goto memerr;
	
	for (i = 0; i < asize; i++)
		if (fgetval(fp, DATATY, &dp->arr.d[i]) <= 0)
			goto scanerr;
	fclose(fp);
	i = hash(dname);
	dp->next = dtab[i];
	return(dtab[i] = dp);
memerr:
	error(SYSTEM, "out of memory in getdata");
scanerr:
	sprintf(errmsg, "%s in data file \"%s\"",
			feof(fp) ? "unexpected EOF" : "bad format", dfname);
	error(USER, errmsg);
	return NULL;	/* pro forma return */
}


static int
headaspect(			/* check string for aspect ratio */
	char  *s,
	void  *iap
)
{
	char	fmt[MAXFMTLEN];

	if (isaspect(s))
		*(double*)iap *= aspectval(s);
	else if (formatval(fmt, s) && strcmp(fmt, COLRFMT))
		*(double*)iap = 0.0;
	return(0);
}

DATARRAY *
getpict(				/* get picture pname */
	char  *pname
)
{
	double  inpaspect;
	char  *pfname;
	FILE  *fp;
	COLR  *scanin;
	int  sl, ns;
	RESOLU	inpres;
	RREAL  loc[2];
	int  y;
	int  x, i;
	DATARRAY  *pp;
						/* look for array in list */
	for (pp = dtab[hash(pname)]; pp != NULL; pp = pp->next)
		if (!strcmp(pname, pp->name))
			return(pp);		/* found! */

	if ((pfname = getpath(pname, getrlibpath(), R_OK)) == NULL) {
		sprintf(errmsg, "cannot find picture file \"%s\"", pname);
		error(SYSTEM, errmsg);
	}
	if ((pp = (DATARRAY *)malloc(3*sizeof(DATARRAY))) == NULL)
		goto memerr;

	pp[0].name = savestr(pname);

	if ((fp = fopen(pfname, "rb")) == NULL) {
		sprintf(errmsg, "cannot open picture file \"%s\"", pfname);
		error(SYSTEM, errmsg);
	}
						/* get dimensions */
	inpaspect = 1.0;
	getheader(fp, headaspect, &inpaspect);
	if (inpaspect <= FTINY || !fgetsresolu(&inpres, fp))
		goto readerr;
	pp[0].nd = 2;
	pp[0].dim[0].ne = inpres.yr;
	pp[0].dim[1].ne = inpres.xr;
	pp[0].dim[0].org =
	pp[0].dim[1].org = 0.0;
	if (inpres.xr <= inpres.yr*inpaspect) {
		pp[0].dim[0].siz = inpaspect *
					(double)inpres.yr/inpres.xr;
		pp[0].dim[1].siz = 1.0;
	} else {
		pp[0].dim[0].siz = 1.0;
		pp[0].dim[1].siz = (double)inpres.xr/inpres.yr /
					inpaspect;
	}
	pp[0].dim[0].p = pp[0].dim[1].p = NULL;
	sl = scanlen(&inpres);				/* allocate array */
	ns = numscans(&inpres);
	i = ns*sl*sizeof(COLR);
#if PSIZWARN
	if (i > PSIZWARN) {				/* memory warning */
		sprintf(errmsg, "picture file \"%s\" using %.1f MB of memory",
				pname, i*(1.0/(1024*1024)));
		error(WARNING, errmsg);
	}
#endif
	if ((pp[0].arr.c = (COLR *)malloc(i)) == NULL)
		goto memerr;
							/* load picture */
	if ((scanin = (COLR *)malloc(sl*sizeof(COLR))) == NULL)
		goto memerr;
	for (y = 0; y < ns; y++) {
		if (freadcolrs(scanin, sl, fp) < 0)
			goto readerr;
		for (x = 0; x < sl; x++) {
			pix2loc(loc, &inpres, x, y);
			i = (int)(loc[1]*inpres.yr)*inpres.xr +
					(int)(loc[0]*inpres.xr);
			copycolr(pp[0].arr.c[i], scanin[x]);
		}
	}
	free(scanin);
	fclose(fp);
	i = hash(pname);
	pp[0].next = dtab[i];		/* link into picture list */
	pp[1] = pp[0];
	pp[2] = pp[0];
	pp[0].type = RED;		/* differentiate RGB records */
	pp[1].type = GRN;
	pp[2].type = BLU;
	return(dtab[i] = pp);
memerr:
	error(SYSTEM, "out of memory in getpict");
readerr:
	sprintf(errmsg, "bad picture file \"%s\"", pfname);
	error(USER, errmsg);
	return NULL;	/* pro forma return */
}


/* header info type for hyperspectral image */
typedef struct {
	float	wlpart[4];	/* wavelength partitions */
	int	nc;		/* number of components */
	double	inpaspect;	/* pixel aspect ratio */
} SPECINFO;

static int
specheadline(				/* get info for spectral image */
	char  *s,
	void  *cdp
)
{
	SPECINFO	*sip = (SPECINFO *)cdp;
	char		fmt[MAXFMTLEN];

	if (isaspect(s))
		sip->inpaspect *= aspectval(s);
	else if (isncomp(s))
		sip->nc = ncompval(s);
	else if (iswlsplit(s))
		wlsplitval(sip->wlpart, s);
	else if (formatval(fmt, s) && strcmp(fmt, SPECFMT))
		return(-1);
	return(0);
}

DATARRAY *
getspec(		/* load hyperspectral image as data */
	char *sname
)
{
	SPECINFO	si;
	char		*pfname;
	FILE		*fp;
	int		sl, ns;
	int		y, i;
	DATARRAY	*pp;
						/* look for array in list */
	for (pp = dtab[hash(sname)]; pp != NULL; pp = pp->next)
		if (!strcmp(sname, pp->name))
			return(pp);		/* found! */

	if ((pfname = getpath(sname, getrlibpath(), R_OK)) == NULL) {
		sprintf(errmsg, "cannot find hyperspectral image \"%s\"", sname);
		error(SYSTEM, errmsg);
	}
	if ((fp = fopen(pfname, "rb")) == NULL) {
		sprintf(errmsg, "cannot open hyperspectral image \"%s\"", pfname);
		error(SYSTEM, errmsg);
	}
	si.wlpart[3] = 0;
	si.nc = 0;
	si.inpaspect = 1.0;
	if (getheader(fp, specheadline, &si) < 0 ||
			(si.nc <= 3) | (si.nc > MAXCSAMP) | (si.wlpart[3] < 1) ||
			!fscnresolu(&sl, &ns, fp))
		goto readerr;

	if ((pp = (DATARRAY *)malloc(sizeof(DATARRAY))) == NULL)
		goto memerr;

	pp->name = savestr(sname);
	pp->type = SPECTY;
	pp->nd = 3;
	pp->dim[0].ne = ns;
	pp->dim[1].ne = sl;
	pp->dim[0].org =
	pp->dim[1].org = 0.0;
	if (sl <= ns*si.inpaspect) {
		pp->dim[0].siz = si.inpaspect * (double)ns/sl;
		pp->dim[1].siz = 1.0;
	} else {
		pp->dim[0].siz = 1.0;
		pp->dim[1].siz = (double)sl/ns / si.inpaspect;
	}
	pp->dim[2].ne = si.nc;
	pp->dim[2].siz = si.wlpart[3] - si.wlpart[0];
	pp->dim[2].org = si.wlpart[0] + 0.5*pp->dim[2].siz/si.nc;
	pp->dim[2].siz *= (si.nc - 1.0)/si.nc;
	pp->dim[0].p = pp->dim[1].p = pp->dim[2].p = NULL;
	i = ns*sl*(si.nc+1);
#if PSIZWARN
	if (i > PSIZWARN) {			/* memory warning */
		sprintf(errmsg, "hyperspectral image \"%s\" using %.1f MB of memory",
				sname, i*(1.0/(1024*1024)));
		error(WARNING, errmsg);
	}
#endif
	if ((pp->arr.s = (uby8 *)malloc(i)) == NULL)
		goto memerr;
	for (y = 0; y < ns; y++)		/* read each scanline */
		if (freadscolrs(pp->arr.s + y*sl*(si.nc+1), si.nc, sl, fp) < 0)
			goto readerr;
	fclose(fp);
	i = hash(sname);			/* insert in hash table */
	pp->next = dtab[i];
	return(dtab[i] = pp);
memerr:
	error(SYSTEM, "out of memory in getspec");
readerr:
	sprintf(errmsg, "bad hyperspectral image \"%s\"", pfname);
	error(USER, errmsg);
	return NULL;	/* pro forma return */
}


void
freedata(			/* release data array reference */
	DATARRAY  *dta
)
{
	DATARRAY  head;
	int  hval, nents;
	DATARRAY  *dpl, *dp;
	int  i;

	if (dta == NULL) {			/* free all if NULL */
		hval = 0; nents = TABSIZ;
	} else {
		hval = hash(dta->name); nents = 1;
		if (!*dta->name) {		/* not a data file? */
			dta->next = dtab[hval];
			dtab[hval] = dta;	/* ...fake position */
		}
	}
	while (nents--) {
		head.next = dtab[hval];
		dpl = &head;
		while ((dp = dpl->next) != NULL)
			if ((dta == NULL) | (dta == dp)) {
				dpl->next = dp->next;
				free(dp->arr.p);
				for (i = 0; i < dp->nd; i++)
					if (dp->dim[i].p != NULL)
						free(dp->dim[i].p);
				freestr(dp->name);
				free(dp);
			} else
				dpl = dp;
		dtab[hval++] = head.next;
	}
}


/* internal call to interpolate data value or vector */
static double
data_interp(DATARRAY *dp, double *pt, double coef, DATATYPE *rvec)
{
	DATARRAY	sd;
	int		stride, i;
	double		x, c0, c1, y0, y1;
					/* set up dimensions for recursion */
	if (dp->nd > 1) {
		sd.name = dp->name;
		sd.type = dp->type;
		sd.nd = dp->nd - 1;
		memcpy(sd.dim, dp->dim+1, sd.nd*sizeof(struct dadim));
		stride = sd.dim[i = sd.nd-1].ne + (sd.type==SPECTY);
		while (i-- > 0)
			stride *= sd.dim[i].ne;
	}
					/* get independent variable */
	if (dp->dim[0].p == NULL) {		/* evenly spaced points */
		x = (pt[0] - dp->dim[0].org)/dp->dim[0].siz;
		x *= (double)(dp->dim[0].ne - 1);
		i = x;
		if (i < 0)
			i = 0;
		else if (i > dp->dim[0].ne - 2)
			i = dp->dim[0].ne - 2;
	} else {				/* unevenly spaced points */
		int	lower, upper;
		if (dp->dim[0].siz > 0.0) {
			lower = 0;
			upper = dp->dim[0].ne;
		} else {
			lower = dp->dim[0].ne;
			upper = 0;
		}
		do {
			i = (lower + upper) >> 1;
			if (pt[0] >= dp->dim[0].p[i])
				lower = i;
			else
				upper = i;
		} while (i != (lower + upper) >> 1);

		if (i > dp->dim[0].ne - 2)
			i = dp->dim[0].ne - 2;

		x = i + (pt[0] - dp->dim[0].p[i]) /
				(dp->dim[0].p[i+1] - dp->dim[0].p[i]);
	}
	/*
	 * Compute interpolation coefficients:
	 * extrapolate as far as one division, then
	 * taper off harmonically to zero.
	 */
	if (x > i+2) {
		c0 = 1./(i-1 - x);
		c1 = -2.*c0;
	} else if (x < i-1) {
		c1 = 1./(i - x);
		c0 = -2.*c1;
	} else {
		c0 = i+1 - x;
		c1 = x - i;
	}
	c0 *= coef;
	c1 *= coef;
					/* check if vector interp */
	if ((dp->nd == 2) & (rvec != NULL)) {
		if (dp->type == DATATY) {
			sd.arr.d = dp->arr.d + i*stride;
			for (i = sd.dim[0].ne; i--; )
				rvec[i] += c0*sd.arr.d[i]
					+ c1*sd.arr.d[i+stride];
			return(0.);
		}
		if (dp->type == SPECTY) {
			double	f;
			sd.arr.s = dp->arr.s + i*stride;
			f = ldexp(1.0, (int)sd.arr.s[sd.dim[0].ne]
					- (COLXS+8));
			for (i = sd.dim[0].ne; i--; )
				rvec[i] += c0*f*(sd.arr.s[i] + 0.5);
			sd.arr.s += stride;
			f = ldexp(1.0, (int)sd.arr.s[sd.dim[0].ne]
					- (COLXS+8));
			for (i = sd.dim[0].ne; i--; )
				rvec[i] += c1*f*(sd.arr.s[i] + 0.5);
			return(0.);
		}
		sd.arr.c = dp->arr.c + i*stride;
		for (i = sd.dim[0].ne; i--; )
			rvec[i] += c0*colrval(sd.arr.c[i],sd.type)
				+ c1*colrval(sd.arr.c[i+stride],sd.type);
		return(0.);
	}
					/* get dependent variable */
	if (dp->nd > 1) {
		if (dp->type == DATATY) {
			sd.arr.d = dp->arr.d + i*stride;
			y0 = data_interp(&sd, pt+1, c0, rvec);
			sd.arr.d += stride;
		} else if (dp->type == SPECTY) {
			sd.arr.s = dp->arr.s + i*stride;
			y0 = data_interp(&sd, pt+1, c0, rvec);
			sd.arr.s += stride;
		} else {
			sd.arr.c = dp->arr.c + i*stride;
			y0 = data_interp(&sd, pt+1, c0, rvec);
			sd.arr.c += stride;
		}
		y1 = data_interp(&sd, pt+1, c1, rvec);
	} else {			/* end of recursion */
		if (dp->type == DATATY) {
			y0 = dp->arr.d[i];
			y1 = dp->arr.d[i+1];
		} else if (dp->type == SPECTY) {
			if (dp->arr.s[dp->dim[0].ne]) {
				double	f = ldexp(1.0, -(COLXS+8) +
						(int)dp->arr.s[dp->dim[0].ne]);
				y0 = (dp->arr.s[i] + 0.5)*f;
				y1 = (dp->arr.s[i+1] + 0.5)*f;
			} else
				y0 = y1 = 0.0;
		} else {
			y0 = colrval(dp->arr.c[i],dp->type);
			y1 = colrval(dp->arr.c[i+1],dp->type);
		}
		y0 *= c0;
		y1 *= c1;
	}
	return(y0 + y1);	/* coefficients already applied */
}


double
datavalue(		/* interpolate data value at a point */
	DATARRAY  *dp,
	double	*pt
)
{
	return(data_interp(dp, pt, 1., NULL));
}


/* Interpolate final vector corresponding to last dimension in data array */
DATARRAY *
datavector(DATARRAY *dp, double *pt)
{
	DATARRAY	*newdp;

	if (dp->nd < 2)
		error(INTERNAL, "datavector() called with 1-D array");
					/* create vector array */
	newdp = (DATARRAY *)malloc(sizeof(DATARRAY) -
				(MAXDDIM-1)*sizeof(struct dadim) +
				sizeof(DATATYPE)*dp->dim[dp->nd-1].ne);
	if (newdp == NULL)
		error(SYSTEM, "out of memory in datavector");
	newdp->next = NULL;
	newdp->name = dp->name;
	newdp->type = DATATY;
	newdp->nd = 1;			/* vector data goes here */
	newdp->dim[0] = dp->dim[dp->nd-1];
	newdp->arr.d = (DATATYPE *)(newdp->dim + 1);
	memset(newdp->arr.d, 0, sizeof(DATATYPE)*newdp->dim[0].ne);

	(void)data_interp(dp, pt, 1., newdp->arr.d);

	return(newdp);			/* will be free'd using free() */
}


#if 0
double
datavalue(		/* interpolate data value at a point */
	DATARRAY  *dp,
	double	*pt
)
{
	DATARRAY  sd;
	int  asize;
	int  lower, upper;
	int  i;
	double	x, y0, y1;
					/* set up dimensions for recursion */
	if (dp->nd > 1) {
		sd.name = dp->name;
		sd.type = dp->type;
		sd.nd = dp->nd - 1;
		asize = 1;
		for (i = 0; i < sd.nd; i++) {
			sd.dim[i].org = dp->dim[i+1].org;
			sd.dim[i].siz = dp->dim[i+1].siz;
			sd.dim[i].p = dp->dim[i+1].p;
			asize *= (sd.dim[i].ne = dp->dim[i+1].ne) +
				((sd.type==SPECTY) & (i==sd.nd-1));
		}
	}
					/* get independent variable */
	if (dp->dim[0].p == NULL) {		/* evenly spaced points */
		x = (pt[0] - dp->dim[0].org)/dp->dim[0].siz;
		x *= (double)(dp->dim[0].ne - 1);
		i = x;
		if (i < 0)
			i = 0;
		else if (i > dp->dim[0].ne - 2)
			i = dp->dim[0].ne - 2;
	} else {				/* unevenly spaced points */
		if (dp->dim[0].siz > 0.0) {
			lower = 0;
			upper = dp->dim[0].ne;
		} else {
			lower = dp->dim[0].ne;
			upper = 0;
		}
		do {
			i = (lower + upper) >> 1;
			if (pt[0] >= dp->dim[0].p[i])
				lower = i;
			else
				upper = i;
		} while (i != (lower + upper) >> 1);

		if (i > dp->dim[0].ne - 2)
			i = dp->dim[0].ne - 2;

		x = i + (pt[0] - dp->dim[0].p[i]) /
				(dp->dim[0].p[i+1] - dp->dim[0].p[i]);
	}
					/* get dependent variable */
	if (dp->nd > 1) {
		if (dp->type == DATATY) {
			sd.arr.d = dp->arr.d + i*asize;
			y0 = datavalue(&sd, pt+1);
			sd.arr.d += asize;
			y1 = datavalue(&sd, pt+1);
		} else if (dp->type == SPECTY) {
			sd.arr.s = dp->arr.s + i*asize;
			y0 = datavalue(&sd, pt+1);
			sd.arr.s += asize;
			y1 = datavalue(&sd, pt+1);
		} else {
			sd.arr.c = dp->arr.c + i*asize;
			y0 = datavalue(&sd, pt+1);
			sd.arr.c += asize;
			y1 = datavalue(&sd, pt+1);
		}
	} else {
		if (dp->type == DATATY) {
			y0 = dp->arr.d[i];
			y1 = dp->arr.d[i+1];
		} else if (dp->type == SPECTY) {
			if (dp->arr.s[dp->dim[0].ne]) {
				double	f = ldexp(1.0, -(COLXS+8) +
						(int)dp->arr.s[dp->dim[0].ne]);
				y0 = (dp->arr.s[i] + 0.5)*f;
				y1 = (dp->arr.s[i+1] + 0.5)*f;
			} else
				y0 = y1 = 0.0;
		} else {
			y0 = colrval(dp->arr.c[i],dp->type);
			y1 = colrval(dp->arr.c[i+1],dp->type);
		}
	}
	/*
	 * Extrapolate as far as one division, then
	 * taper off harmonically to zero.
	 */
	if (x > i+2)
		return( (2*y1-y0)/(x-(i-1)) );

	if (x < i-1)
		return( (2*y0-y1)/(i-x) );

	return( y0*((i+1)-x) + y1*(x-i) );
}
#endif
