#ifndef lint
static const char	RCSid[] = "$Id: data.c,v 2.17 2003/05/13 17:58:33 greg Exp $";
#endif
/*
 *  data.c - routines dealing with interpolated data.
 */

#include "copyright.h"

#include  "standard.h"

#include  "color.h"

#include  "resolu.h"

#include  "data.h"

				/* picture memory usage before warning */
#ifndef PSIZWARN
#ifdef BIGMEM
#define PSIZWARN	5000000
#else
#define PSIZWARN	1500000
#endif
#endif

#ifndef TABSIZ
#define TABSIZ		97		/* table size (prime) */
#endif

#define hash(s)		(shash(s)%TABSIZ)


static DATARRAY	 *dtab[TABSIZ];		/* data array list */


DATARRAY *
getdata(dname)				/* get data array dname */
char  *dname;
{
	char  *dfname;
	FILE  *fp;
	int  asize;
	register int  i, j;
	register DATARRAY  *dp;
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
		error(USER, errmsg);
	}
	if ((fp = fopen(dfname, "r")) == NULL) {
		sprintf(errmsg, "cannot open data file \"%s\"", dfname);
		error(SYSTEM, errmsg);
	}
							/* get dimensions */
	if (fgetval(fp, 'i', (char *)&asize) <= 0)
		goto scanerr;
	if (asize <= 0 | asize > MAXDDIM) {
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
		if (fgetval(fp, DATATY, (char *)&dp->dim[i].org) <= 0)
			goto scanerr;
		if (fgetval(fp, DATATY, (char *)&dp->dim[i].siz) <= 0)
			goto scanerr;
		if (fgetval(fp, 'i', (char *)&dp->dim[i].ne) <= 0)
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
				if (fgetval(fp, DATATY,
						(char *)&dp->dim[i].p[j]) <= 0)
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
		if (fgetval(fp, DATATY, (char *)&dp->arr.d[i]) <= 0)
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
}


static int
headaspect(s, iap)			/* check string for aspect ratio */
char  *s;
double  *iap;
{
	char	fmt[32];

	if (isaspect(s))
		*iap *= aspectval(s);
	else if (formatval(fmt, s) && !globmatch(PICFMT, fmt))
		*iap = 0.0;
	return(0);
}


DATARRAY *
getpict(pname)				/* get picture pname */
char  *pname;
{
	double  inpaspect;
	char  *pfname;
	FILE  *fp;
	COLR  *scanin;
	int  sl, ns;
	RESOLU	inpres;
	FLOAT  loc[2];
	int  y;
	register int  x, i;
	register DATARRAY  *pp;
						/* look for array in list */
	for (pp = dtab[hash(pname)]; pp != NULL; pp = pp->next)
		if (!strcmp(pname, pp->name))
			return(pp);		/* found! */

	if ((pfname = getpath(pname, getrlibpath(), R_OK)) == NULL) {
		sprintf(errmsg, "cannot find picture file \"%s\"", pname);
		error(USER, errmsg);
	}
	if ((pp = (DATARRAY *)malloc(3*sizeof(DATARRAY))) == NULL)
		goto memerr;

	pp[0].name = savestr(pname);

	if ((fp = fopen(pfname, "r")) == NULL) {
		sprintf(errmsg, "cannot open picture file \"%s\"", pfname);
		error(SYSTEM, errmsg);
	}
#ifdef MSDOS
	setmode(fileno(fp), O_BINARY);
#endif
						/* get dimensions */
	inpaspect = 1.0;
	getheader(fp, headaspect, (char *)&inpaspect);
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
		sprintf(errmsg, "picture file \"%s\" using %d bytes of memory",
				pname, i);
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
	free((void *)scanin);
	fclose(fp);
	i = hash(pname);
	pp[0].next = dtab[i];		/* link into picture list */
	copystruct(&pp[1], &pp[0]);
	copystruct(&pp[2], &pp[0]);
	pp[0].type = RED;		/* differentiate RGB records */
	pp[1].type = GRN;
	pp[2].type = BLU;
	return(dtab[i] = pp);

memerr:
	error(SYSTEM, "out of memory in getpict");
readerr:
	sprintf(errmsg, "bad picture file \"%s\"", pfname);
	error(USER, errmsg);
}


void
freedata(dta)			/* release data array reference */
DATARRAY  *dta;
{
	DATARRAY  head;
	int  hval, nents;
	register DATARRAY  *dpl, *dp;
	register int  i;

	if (dta == NULL) {			/* free all if NULL */
		hval = 0; nents = TABSIZ;
	} else {
		hval = hash(dta->name); nents = 1;
	}
	while (nents--) {
		head.next = dtab[hval];
		dpl = &head;
		while ((dp = dpl->next) != NULL)
			if ((dta == NULL | dta == dp)) {
				dpl->next = dp->next;
				if (dp->type == DATATY)
					free((void *)dp->arr.d);
				else
					free((void *)dp->arr.c);
				for (i = 0; i < dp->nd; i++)
					if (dp->dim[i].p != NULL)
						free((void *)dp->dim[i].p);
				freestr(dp->name);
				free((void *)dp);
			} else
				dpl = dp;
		dtab[hval++] = head.next;
	}
}


double
datavalue(dp, pt)		/* interpolate data value at a point */
register DATARRAY  *dp;
double	*pt;
{
	DATARRAY  sd;
	int  asize;
	int  lower, upper;
	register int  i;
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
			asize *= sd.dim[i].ne = dp->dim[i+1].ne;
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
			sd.arr.d = dp->arr.d + (i+1)*asize;
			y1 = datavalue(&sd, pt+1);
		} else {
			sd.arr.c = dp->arr.c + i*asize;
			y0 = datavalue(&sd, pt+1);
			sd.arr.c = dp->arr.c + (i+1)*asize;
			y1 = datavalue(&sd, pt+1);
		}
	} else {
		if (dp->type == DATATY) {
			y0 = dp->arr.d[i];
			y1 = dp->arr.d[i+1];
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
