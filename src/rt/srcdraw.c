#ifndef lint
static const char	RCSid[] = "$Id: srcdraw.c,v 2.5 2003/02/22 02:07:29 greg Exp $";
#endif
/*
 * Draw small sources into image in case we missed them.
 *
 *  External symbols declared in ray.h
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

#include  "ray.h"

#include  "view.h"

#include  "source.h"


#define	CLIP_ABOVE	1
#define CLIP_BELOW	2
#define CLIP_RIGHT	3
#define CLIP_LEFT	4

#define MAXVERT		10

typedef struct splist {
	struct splist	*next;			/* next source in list */
	int	sn;				/* source number */
	short	nv;				/* number of vertices */
	FLOAT	vl[3][2];			/* vertex array (last) */
} SPLIST;				/* source polygon list */

extern VIEW	ourview;		/* our view parameters */
extern int	hres, vres;		/* our image resolution */
static SPLIST	*sphead = NULL;		/* our list of source polys */


static int
inregion(p, cv, crit)			/* check if vertex is in region */
FLOAT	p[2];
double	cv;
int	crit;
{
	switch (crit) {
	case CLIP_ABOVE:
		return(p[1] < cv);
	case CLIP_BELOW:
		return(p[1] >= cv);
	case CLIP_RIGHT:
		return(p[0] < cv);
	case CLIP_LEFT:
		return(p[0] >= cv);
	}
	return(-1);
}


static
clipregion(a, b, cv, crit, r)		/* find intersection with boundary */
register FLOAT	a[2], b[2];
double	cv;
int	crit;
FLOAT	r[2];	/* return value */
{
	switch (crit) {
	case CLIP_ABOVE:
	case CLIP_BELOW:
		r[1] = cv;
		r[0] = a[0] + (cv-a[1])/(b[1]-a[1])*(b[0]-a[0]);
		return;
	case CLIP_RIGHT:
	case CLIP_LEFT:
		r[0] = cv;
		r[1] = a[1] + (cv-a[0])/(b[0]-a[0])*(b[1]-a[1]);
		return;
	}
}


static int
hp_clip_poly(vl, nv, cv, crit, vlo)	/* clip polygon to half-plane */
FLOAT	vl[][2];
int	nv;
double	cv;
int	crit;
FLOAT	vlo[][2];	/* return value */
{
	FLOAT	*s, *p;
	register int	j, nvo;

	s = vl[nv-1];
	nvo = 0;
	for (j = 0; j < nv; j++) {
		p = vl[j];
		if (inregion(p, cv, crit)) {
			if (!inregion(s, cv, crit))
				clipregion(s, p, cv, crit, vlo[nvo++]);
			vlo[nvo][0] = p[0]; vlo[nvo++][1] = p[1];
		} else if (inregion(s, cv, crit))
			clipregion(s, p, cv, crit, vlo[nvo++]);
		s = p;
	}
	return(nvo);
}


static int
box_clip_poly(vl, nv, xl, xr, yb, ya, vlo)	/* clip polygon to box */
FLOAT	vl[MAXVERT][2];
int	nv;
double	xl, xr, yb, ya;
FLOAT	vlo[MAXVERT][2];	/* return value */
{
	FLOAT	vlt[MAXVERT][2];
	int	nvt, nvo;

	nvt = hp_clip_poly(vl, nv, yb, CLIP_BELOW, vlt);
	nvo = hp_clip_poly(vlt, nvt, ya, CLIP_ABOVE, vlo);
	nvt = hp_clip_poly(vlo, nvo, xl, CLIP_LEFT, vlt);
	nvo = hp_clip_poly(vlt, nvt, xr, CLIP_RIGHT, vlo);

	return(nvo);
}


static double
minw2(vl, nv, ar2)			/* compute square of minimum width */
FLOAT	vl[][2];
int	nv;
double	ar2;
{
	double	d2, w2, w2min, w2max;
	register FLOAT	*p0, *p1, *p2;
	int	i, j;
				/* find minimum for all widths */
	w2min = FHUGE;
	p0 = vl[nv-1];
	for (i = 0; i < nv; i++) {	/* for each edge */
		p1 = vl[i];
		d2 = (p1[0]-p0[0])*(p1[0]-p0[0]) +
				(p1[1]-p0[1])*(p1[1]-p0[1])*ar2;
		w2max = 0.;		/* find maximum for this side */
		for (j = 1; j < nv-1; j++) {
			p2 = vl[(i+j)%nv];
			w2 = (p1[0]-p0[0])*(p2[1]-p0[1]) -
					(p1[1]-p0[1])*(p2[0]-p0[0]);
			w2 = w2*w2*ar2/d2;	/* triangle height squared */
			if (w2 > w2max)
				w2max = w2;
		}
		if (w2max < w2min)	/* global min. based on local max.'s */
			w2min = w2max;
		p0 = p1;
	}
	return(w2min);
}


static
convex_center(vl, nv, cv)		/* compute center of convex polygon */
register FLOAT	vl[][2];
int	nv;
FLOAT	cv[2];		/* return value */
{
	register int	i;
					/* simple average (suboptimal) */
	cv[0] = cv[1] = 0.;
	for (i = 0; i < nv; i++) {
		cv[0] += vl[i][0];
		cv[1] += vl[i][1];
	}
	cv[0] /= (double)nv;
	cv[1] /= (double)nv;
}


static double
poly_area(vl, nv)			/* compute area of polygon */
register FLOAT	vl[][2];
int	nv;
{
	double	a;
	FLOAT	v0[2], v1[2];
	register int	i;

	a = 0.;
	v0[0] = vl[1][0] - vl[0][0];
	v0[1] = vl[1][1] - vl[0][1];
	for (i = 2; i < nv; i++) {
		v1[0] = vl[i][0] - vl[0][0];
		v1[1] = vl[i][1] - vl[0][1];
		a += v0[0]*v1[1] - v0[1]*v1[0];
		v0[0] = v1[0]; v0[1] = v1[1];
	}
	return(a * (a >= 0. ? .5 : -.5));
}


static int
convex_hull(vl, nv, vlo)		/* compute polygon's convex hull */
FLOAT	vl[][2];
int	nv;
FLOAT	vlo[][2];	/* return value */
{
	int	nvo, nvt;
	FLOAT	vlt[MAXVERT][2];
	double	voa, vta;
	register int	i, j;
					/* start with original polygon */
	for (i = nvo = nv; i--; ) {
		vlo[i][0] = vl[i][0]; vlo[i][1] = vl[i][1];
	}
	voa = poly_area(vlo, nvo);	/* compute its area */
	for (i = 0; i < nvo; i++) {		/* for each output vertex */
		for (j = 0; j < i; j++) {
			vlt[j][0] = vlo[j][0]; vlt[j][1] = vlo[j][1];
		}
		nvt = nvo - 1;			/* make poly w/o vertex */
		for (j = i; j < nvt; j++) {
			vlt[j][0] = vlo[j+1][0]; vlt[j][1] = vlo[j+1][1];
		}
		vta = poly_area(vlt, nvt);
		if (vta >= voa) {		/* is simpler poly bigger? */
			voa = vta;			/* then use it */
			for (j = nvo = nvt; j--; ) {
				vlo[j][0] = vlt[j][0]; vlo[j][1] = vlt[j][1];
			}
			i--;				/* next adjust */
		}
	}
	return(nvo);
}


static
spinsert(sn, vl, nv)			/* insert new source polygon */
int	sn;
FLOAT	vl[][2];
int	nv;
{
	register SPLIST	*spn;
	register int	i;

	if (nv < 3)
		return;
	if (nv > 3)
		spn = (SPLIST *)malloc(sizeof(SPLIST)+sizeof(FLOAT)*2*(nv-3));
	else
		spn = (SPLIST *)malloc(sizeof(SPLIST));
	if (spn == NULL)
		error(SYSTEM, "out of memory in spinsert");
	spn->sn = sn;
	for (i = spn->nv = nv; i--; ) {
		spn->vl[i][0] = vl[i][0]; spn->vl[i][1] = vl[i][1];
	}
	spn->next = sphead;		/* push onto global list */
	sphead = spn;
}


int
sourcepoly(sn, sp)			/* compute image polygon for source */
int	sn;
FLOAT	sp[MAXVERT][2];
{
	static char	cubeord[8][6] = {{1,3,2,6,4,5},{0,4,5,7,3,2},
					 {0,1,3,7,6,4},{0,1,5,7,6,2},
					 {0,2,6,7,5,1},{0,4,6,7,3,1},
					 {0,2,3,7,5,4},{1,5,4,6,2,3}};
	register SRCREC	*s = source + sn;
	FVECT	ap, ip;
	FLOAT	pt[6][2];
	int	dir;
	register int	i, j;

	if (s->sflags & (SDISTANT|SFLAT)) {
		if (s->sflags & SDISTANT && ourview.type == VT_PAR)
			return(0);		/* all or nothing case */
		if (s->sflags & SFLAT) {
			for (i = 0; i < 3; i++)
				ap[i] = s->sloc[i] - ourview.vp[i];
			if (DOT(ap, s->snorm) >= 0.)
				return(0);	/* source faces away */
		}
		for (j = 0; j < 4; j++) {	/* four corners */
			for (i = 0; i < 3; i++) {
				ap[i] = s->sloc[i];
				if (j==1|j==2) ap[i] += s->ss[SU][i];
				else ap[i] -= s->ss[SU][i];
				if (j==2|j==3) ap[i] += s->ss[SV][i];
				else ap[i] -= s->ss[SV][i];
				if (s->sflags & SDISTANT) {
					ap[i] *= 1. + ourview.vfore;
					ap[i] += ourview.vp[i];
				}
			}
			viewloc(ip, &ourview, ap);	/* find image point */
			if (ip[2] <= 0.)
				return(0);		/* in front of view */
			sp[j][0] = ip[0]; sp[j][1] = ip[1];
		}
		return(4);
	}
					/* identify furthest corner */
	for (i = 0; i < 3; i++)
		ap[i] = s->sloc[i] - ourview.vp[i];
	dir =	(DOT(ap,s->ss[SU])>0.) |
		(DOT(ap,s->ss[SV])>0.)<<1 |
		(DOT(ap,s->ss[SW])>0.)<<2 ;
					/* order vertices based on this */
	for (j = 0; j < 6; j++) {
		for (i = 0; i < 3; i++) {
			ap[i] = s->sloc[i];
			if (cubeord[dir][j] & 1) ap[i] += s->ss[SU][i];
			else ap[i] -= s->ss[SU][i];
			if (cubeord[dir][j] & 2) ap[i] += s->ss[SV][i];
			else ap[i] -= s->ss[SV][i];
			if (cubeord[dir][j] & 4) ap[i] += s->ss[SW][i];
			else ap[i] -= s->ss[SW][i];
		}
		viewloc(ip, &ourview, ap);	/* find image point */
		if (ip[2] <= 0.)
			return(0);		/* in front of view */
		pt[j][0] = ip[0]; pt[j][1] = ip[1];
	}
	return(convex_hull(pt, 6, sp));		/* make sure it's convex */
}


			/* initialize by finding sources smaller than rad */
init_drawsources(rad)
int	rad;				/* source sample size */
{
	FLOAT	spoly[MAXVERT][2];
	int	nsv;
	register SPLIST	*sp;
	register int	i;
					/* free old source list if one */
	for (sp = sphead; sp != NULL; sp = sphead) {
		sphead = sp->next;
		free((void *)sp);
	}
					/* loop through all sources */
	for (i = nsources; i--; ) {
					/* compute image polygon for source */
		if (!(nsv = sourcepoly(i, spoly)))
			continue;
					/* clip to image boundaries */
		if (!(nsv = box_clip_poly(spoly, nsv, 0., 1., 0., 1., spoly)))
			continue;
					/* big enough for standard sampling? */
		if (minw2(spoly, nsv, ourview.vn2/ourview.hn2) >
				(double)rad*rad/hres/hres)
			continue;
					/* OK, add to our list */
		spinsert(i, spoly, nsv);
	}
}

void			/* add sources smaller than rad to computed subimage */
drawsources(pic, zbf, x0, xsiz, y0, ysiz)
COLOR	*pic[];				/* subimage pixel value array */
float	*zbf[];				/* subimage distance array (opt.) */
int	x0, xsiz, y0, ysiz;		/* origin and size of subimage */
{
	FLOAT	spoly[MAXVERT][2], ppoly[MAXVERT][2];
	int	nsv, npv;
	int	xmin, xmax, ymin, ymax, x, y;
	FLOAT	cxy[2];
	double	w;
	RAY	sr;
	register SPLIST	*sp;
	register int	i;
					/* check each source in our list */
	for (sp = sphead; sp != NULL; sp = sp->next) {
					/* clip source poly to subimage */
		nsv = box_clip_poly(sp->vl, sp->nv,
				(double)x0/hres, (double)(x0+xsiz)/hres,
				(double)y0/vres, (double)(y0+ysiz)/vres, spoly);
		if (!nsv)
			continue;
					/* find common subimage (BBox) */
		xmin = x0 + xsiz; xmax = x0;
		ymin = y0 + ysiz; ymax = y0;
		for (i = 0; i < nsv; i++) {
			if ((double)xmin/hres > spoly[i][0])
				xmin = spoly[i][0]*hres + FTINY;
			if ((double)xmax/hres < spoly[i][0])
				xmax = spoly[i][0]*hres - FTINY;
			if ((double)ymin/vres > spoly[i][1])
				ymin = spoly[i][1]*vres + FTINY;
			if ((double)ymax/vres < spoly[i][1])
				ymax = spoly[i][1]*vres - FTINY;
		}
					/* evaluate each pixel in BBox */
		for (y = ymin; y <= ymax; y++)
			for (x = xmin; x <= xmax; x++) {
							/* subarea for pixel */
				npv = box_clip_poly(spoly, nsv,
						(double)x/hres, (x+1.)/hres,
						(double)y/vres, (y+1.)/vres,
						ppoly);
				if (!npv)
					continue;	/* no overlap */
				convex_center(ppoly, npv, cxy);
				if ((sr.rmax = viewray(sr.rorg,sr.rdir,&ourview,
						cxy[0],cxy[1])) < -FTINY)
					continue;	/* not in view */
				if (source[sp->sn].sflags & SSPOT &&
						spotout(&sr, source[sp->sn].sl.s))
					continue;	/* outside spot */
				rayorigin(&sr, NULL, SHADOW, 1.0);
				sr.rsrc = sp->sn;
				rayvalue(&sr);		/* compute value */
				if (bright(sr.rcol) <= FTINY)
					continue;	/* missed/blocked */
							/* modify pixel */
				if (zbf[y-y0] != NULL &&
						sr.rt < 0.999*zbf[y-y0][x-x0])
					zbf[y-y0][x-x0] = sr.rt;
				else if (!bigdiff(sr.rcol, pic[y-y0][x-x0],
						0.001))	/* source sample */
					setcolor(pic[y-y0][x-x0], 0., 0., 0.);
				w = poly_area(ppoly, npv) * hres * vres;
				scalecolor(sr.rcol, w);
				scalecolor(pic[y-y0][x-x0], 1.-w);
				addcolor(pic[y-y0][x-x0], sr.rcol);
			}
	}
}
