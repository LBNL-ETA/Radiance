#ifndef lint
static const char	RCSid[] = "$Id: holo.c,v 3.18 2003/06/26 00:58:10 schorsch Exp $";
#endif
/*
 * Routines for converting holodeck coordinates, etc.
 *
 *	10/22/97	GWLarson
 */

#include "holo.h"

float	hd_depthmap[DCINF-DCLIN];

int	hdwg0[6] = {1,1,2,2,0,0};
int	hdwg1[6] = {2,2,0,0,1,1};

static double	logstep;


hdcompgrid(hp)			/* compute derived grid vector and index */
register HOLO	*hp;
{
	double	d;
	register int	i, j;
				/* initialize depth map */
	if (hd_depthmap[0] < 1.) {
		d = 1. + .5/DCLIN;
		for (i = 0; i < DCINF-DCLIN; i++) {
			hd_depthmap[i] = d;
			d *= 1. + 1./DCLIN;
		}
		logstep = log(1. + 1./DCLIN);
	}
				/* compute grid coordinate vectors */
	for (i = 0; i < 3; i++) {
		fcross(hp->wg[i], hp->xv[(i+1)%3], hp->xv[(i+2)%3]);
		d = DOT(hp->wg[i],hp->xv[i]);
		if (d <= FTINY & d >= -FTINY)
			error(USER, "degenerate holodeck section");
		d = hp->grid[i] / d;
		hp->wg[i][0] *= d; hp->wg[i][1] *= d; hp->wg[i][2] *= d;
	}
				/* compute linear depth range */
	hp->tlin = VLEN(hp->xv[0]) + VLEN(hp->xv[1]) + VLEN(hp->xv[2]);
				/* compute wall super-indices from grid */
	hp->wi[0] = 1;		/**** index values begin at 1 ****/
	for (i = 1; i < 6; i++) {
		hp->wi[i] = 0;
		for (j = i; j < 6; j++)
			hp->wi[i] += hp->grid[hdwg0[j]] * hp->grid[hdwg1[j]];
		hp->wi[i] *= hp->grid[hdwg0[i-1]] * hp->grid[hdwg1[i-1]];
		hp->wi[i] += hp->wi[i-1];
	}
}


hdbcoord(gc, hp, i)		/* compute beam coordinates from index */
GCOORD	gc[2];		/* returned */
register HOLO	*hp;
register int	i;
{
	register int	j, n;
	int	n2, reverse;
	GCOORD	g2[2];
					/* check range */
	if (i < 1 | i > nbeams(hp))
		return(0);
	if (reverse = i >= hp->wi[5])
		i -= hp->wi[5] - 1;
	for (j = 0; j < 5; j++)		/* find w0 */
		if (hp->wi[j+1] > i)
			break;
	i -= hp->wi[gc[0].w=j];
					/* find w1 */
	n2 = hp->grid[hdwg0[j]] * hp->grid[hdwg1[j]];
	while (++j < 5)	{
		n = n2 * hp->grid[hdwg0[j]] * hp->grid[hdwg1[j]];
		if (n > i)
			break;
		i -= n;
	}
	gc[1].w = j;
					/* find position on w0 */
	n2 = hp->grid[hdwg0[j]] * hp->grid[hdwg1[j]];
	n = i / n2;
	gc[0].i[1] = n / hp->grid[hdwg0[gc[0].w]];
	gc[0].i[0] = n - gc[0].i[1]*hp->grid[hdwg0[gc[0].w]];
	i -= n*n2;
					/* find position on w1 */
	gc[1].i[1] = i / hp->grid[hdwg0[gc[1].w]];
	gc[1].i[0] = i - gc[1].i[1]*hp->grid[hdwg0[gc[1].w]];
	if (reverse) {
		copystruct(g2, gc+1);
		copystruct(gc+1, gc);
		copystruct(gc, g2);
	}
	return(1);			/* we're done */
}


int
hdbindex(hp, gc)		/* compute index from beam coordinates */
register HOLO	*hp;
register GCOORD	gc[2];
{
	GCOORD	g2[2];
	int	reverse;
	register int	i, j;
					/* check ordering and limits */
	if (reverse = gc[0].w > gc[1].w) {
		copystruct(g2, gc+1);
		copystruct(g2+1, gc);
		gc = g2;
	} else if (gc[0].w == gc[1].w)
		return(0);
	if (gc[0].w < 0 | gc[1].w > 5)
		return(0);
	i = 0;				/* compute index */
	for (j = gc[0].w+1; j < gc[1].w; j++)
		i += hp->grid[hdwg0[j]] * hp->grid[hdwg1[j]];
	i *= hp->grid[hdwg0[gc[0].w]] * hp->grid[hdwg1[gc[0].w]];
	i += hp->wi[gc[0].w];
	i += (hp->grid[hdwg0[gc[0].w]]*gc[0].i[1] + gc[0].i[0]) *
			hp->grid[hdwg0[gc[1].w]] * hp->grid[hdwg1[gc[1].w]] ;
	i += hp->grid[hdwg0[gc[1].w]]*gc[1].i[1] + gc[1].i[0];
	if (reverse)
		i += hp->wi[5] - 1;
	return(i);
}


hdcell(cp, hp, gc)		/* compute cell coordinates */
register FVECT	cp[4];	/* returned (may be passed as FVECT cp[2][2]) */
register HOLO	*hp;
register GCOORD	*gc;
{
	register RREAL	*v;
	double	d;
					/* compute common component */
	VCOPY(cp[0], hp->orig);
	if (gc->w & 1) {
		v = hp->xv[gc->w>>1];
		cp[0][0] += v[0]; cp[0][1] += v[1]; cp[0][2] += v[2];
	}
	v = hp->xv[hdwg0[gc->w]];
	d = (double)gc->i[0] / hp->grid[hdwg0[gc->w]];
	VSUM(cp[0], cp[0], v, d);
	v = hp->xv[hdwg1[gc->w]];
	d = (double)gc->i[1] / hp->grid[hdwg1[gc->w]];
	VSUM(cp[0], cp[0], v, d);
					/* compute x1 sums */
	v = hp->xv[hdwg0[gc->w]];
	d = 1.0 / hp->grid[hdwg0[gc->w]];
	VSUM(cp[1], cp[0], v, d);
	VSUM(cp[3], cp[0], v, d);
					/* compute y1 sums */
	v = hp->xv[hdwg1[gc->w]];
	d = 1.0 / hp->grid[hdwg1[gc->w]];
	VSUM(cp[2], cp[0], v, d);
	VSUM(cp[3], cp[3], v, d);
}


hdlseg(lseg, hp, gc)			/* compute line segment for beam */
register int	lseg[2][3];
register HOLO	*hp;
GCOORD	gc[2];
{
	register int	k;

	for (k = 0; k < 2; k++) {		/* compute end points */
		lseg[k][gc[k].w>>1] = gc[k].w&1	? hp->grid[gc[k].w>>1]-1 : 0 ;
		lseg[k][hdwg0[gc[k].w]] = gc[k].i[0];
		lseg[k][hdwg1[gc[k].w]] = gc[k].i[1];
	}
	return(1);
}


unsigned
hdcode(hp, d)			/* compute depth code for d */
HOLO	*hp;
double	d;
{
	double	tl = hp->tlin;
	register long	c;

	if (d <= 0.)
		return(0);
	if (d >= .99*FHUGE)
		return(DCINF);
	if (d < tl)
		return((unsigned)(d*DCLIN/tl));
	c = (long)(log(d/tl)/logstep) + DCLIN;
	return(c > DCINF ? (unsigned)DCINF : (unsigned)c);
}


hdgrid(gp, hp, wp)		/* compute grid coordinates */
FVECT	gp;		/* returned */
register HOLO	*hp;
FVECT	wp;
{
	FVECT	vt;

	VSUB(vt, wp, hp->orig);
	gp[0] = DOT(vt, hp->wg[0]);
	gp[1] = DOT(vt, hp->wg[1]);
	gp[2] = DOT(vt, hp->wg[2]);
}


hdworld(wp, hp, gp)		/* compute world coordinates */
register FVECT	wp;
register HOLO	*hp;
FVECT	gp;
{
	register double	d;

	d = gp[0]/hp->grid[0];
	VSUM(wp, hp->orig, hp->xv[0], d);

	d = gp[1]/hp->grid[1];
	VSUM(wp, wp, hp->xv[1], d);

	d = gp[2]/hp->grid[2];
	VSUM(wp, wp, hp->xv[2], d);
}


double
hdray(ro, rd, hp, gc, r)	/* compute ray within a beam */
FVECT	ro, rd;		/* returned */
HOLO	*hp;
GCOORD	gc[2];
BYTE	r[2][2];
{
	FVECT	cp[4], p[2];
	register int	i, j;
	double	d0, d1;
					/* compute entry and exit points */
	for (i = 0; i < 2; i++) {
		hdcell(cp, hp, gc+i);
		d0 = (1./256.)*(r[i][0]+.5);
		d1 = (1./256.)*(r[i][1]+.5);
		for (j = 0; j < 3; j++)
			p[i][j] = (1.-d0-d1)*cp[0][j] +
					d0*cp[1][j] + d1*cp[2][j];
	}
	VCOPY(ro, p[0]);		/* assign ray origin and direction */
	VSUB(rd, p[1], p[0]);
	return(normalize(rd));		/* return maximum inside distance */
}


double
hdinter(gc, r, ed, hp, ro, rd)	/* compute ray intersection with section */
register GCOORD	gc[2];	/* returned */
BYTE	r[2][2];	/* returned (optional) */
double	*ed;		/* returned (optional) */
register HOLO	*hp;
FVECT	ro, rd;		/* normalization of rd affects distances */
{
	FVECT	p[2], vt;
	double	d, t0, t1, d0, d1;
	register RREAL	*v;
	register int	i;
					/* first, intersect walls */
	gc[0].w = gc[1].w = -1;
	t0 = -FHUGE; t1 = FHUGE;
	VSUB(vt, ro, hp->orig);
	for (i = 0; i < 3; i++) {		/* for each wall pair */
		d = -DOT(rd, hp->wg[i]);	/* plane distance */
		if (d <= FTINY && d >= -FTINY)	/* check for parallel */
			continue;
		d1 = DOT(vt, hp->wg[i]);	/* ray distances */
		d0 = d1 / d;
		d1 = (d1 - hp->grid[i]) / d;
		if (d < 0) {			/* check against best */
			if (d0 > t0) {
				t0 = d0;
				gc[0].w = i<<1;
			}
			if (d1 < t1) {
				t1 = d1;
				gc[1].w = i<<1 | 1;
			}
		} else {
			if (d1 > t0) {
				t0 = d1;
				gc[0].w = i<<1 | 1;
			}
			if (d0 < t1) {
				t1 = d0;
				gc[1].w = i<<1;
			}
		}
	}
	if (gc[0].w < 0 | gc[1].w < 0)		/* paranoid check */
		return(FHUGE);
						/* compute intersections */
	VSUM(p[0], ro, rd, t0);
	VSUM(p[1], ro, rd, t1);
					/* now, compute grid coordinates */
	for (i = 0; i < 2; i++) {
		VSUB(vt, p[i], hp->orig);
		v = hp->wg[hdwg0[gc[i].w]];
		d = DOT(vt, v);
		if (d < 0 || d >= hp->grid[hdwg0[gc[i].w]])
			return(FHUGE);		/* outside wall */
		gc[i].i[0] = d;
		if (r != NULL)
			r[i][0] = 256. * (d - gc[i].i[0]);
		v = hp->wg[hdwg1[gc[i].w]];
		d = DOT(vt, v);
		if (d < 0 || d >= hp->grid[hdwg1[gc[i].w]])
			return(FHUGE);		/* outside wall */
		gc[i].i[1] = d;
		if (r != NULL)
			r[i][1] = 256. * (d - gc[i].i[1]);
	}
	if (ed != NULL)			/* assign distance to exit point */
		*ed = t1;
	return(t0);			/* return distance to entry point */
}
