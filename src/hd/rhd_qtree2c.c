#ifndef lint
static const char	RCSid[] = "$Id: rhd_qtree2c.c,v 3.3 2004/01/01 11:21:55 schorsch Exp $";
#endif
/*
 * Quadtree display support routines for cone output.
 */

#include "standard.h"
#include "rhd_qtree.h"

static void redraw(RTREE *tp, int x0, int y0, int x1, int y1, int l[2][2]);
static void cpaint(BYTE rgb[3], float *p, int x0, int y0, int x1, int y1);
static void update(BYTE ca[3], RTREE *tp, int x0, int y0, int x1, int y1);



static void
redraw(	/* mark portion of a tree for redraw */
	register RTREE	*tp,
	int	x0,
	int	y0,
	int	x1,
	int	y1,
	int	l[2][2]
)
{
	int	quads = CH_ANY;
	int	mx, my;
	register int	i;
					/* compute midpoint */
	mx = (x0 + x1) >> 1;
	my = (y0 + y1) >> 1;
					/* see what to do */
	if (l[0][0] > mx)
		quads &= ~(CHF(UL)|CHF(DL));
	else if (l[0][1] < mx)
		quads &= ~(CHF(UR)|CHF(DR));
	if (l[1][0] > my)
		quads &= ~(CHF(DR)|CHF(DL));
	else if (l[1][1] < my)
		quads &= ~(CHF(UR)|CHF(UL));
	tp->flgs |= quads;		/* mark quadrants for update */
					/* climb the branches */
	for (i = 0; i < 4; i++)
		if (tp->flgs & BRF(i) && quads & CHF(i))
			redraw(tp->k[i].b, i&01 ? mx : x0, i&02 ? my : y0,
					i&01 ? x1 : mx, i&02 ? y1 : my, l);
}


static void
cpaint(	/* paint a cone within a rectangle */
	BYTE	rgb[3],
	register float	*p,
	int	x0,
	int	y0,
	int	x1,
	int	y1
)
{
	static FVECT	ip, wp;
	double	rad;
						/* compute base radius */
	rad = (double)(y1 - y0 + x1 - x0)/(odev.hres + odev.vres);
						/* approximate apex pos? */
	if (p == NULL || y1-y0 <= qtMinNodesiz || x1-x0 <= qtMinNodesiz) {
		ip[0] = (double)(x0 + x1)/(odev.hres<<1);
		ip[1] = (double)(y0 + y1)/(odev.vres<<1);
		if (p != NULL) {
			wp[0] = p[0] - odev.v.vp[0];
			wp[1] = p[1] - odev.v.vp[1];
			wp[2] = p[2] - odev.v.vp[2];
			if (odev.v.type == VT_PER)
				ip[2] = DOT(wp,odev.v.vdir);
			else
				ip[2] = VLEN(wp);
		} else
			ip[2] = FHUGE;
	} else if (odev.v.type == VT_PER) {	/* special case (faster) */
		wp[0] = p[0] - odev.v.vp[0];
		wp[1] = p[1] - odev.v.vp[1];
		wp[2] = p[2] - odev.v.vp[2];
		ip[2] = DOT(wp,odev.v.vdir);
		ip[0] = DOT(wp,odev.v.hvec)/(ip[2]*odev.v.hn2) +
				0.5 - odev.v.hoff;
		ip[1] = DOT(wp,odev.v.vvec)/(ip[2]*odev.v.vn2) +
				0.5 - odev.v.voff;
	} else {				/* general case */
		VCOPY(wp, p);
		viewloc(ip, &odev.v, wp);
	}
	dev_cone(rgb, ip, rad);
}


static void
update(	/* update tree display as needed */
	BYTE	ca[3],		/* returned average color */
	register RTREE	*tp,
	int	x0,
	int	y0,
	int	x1,
	int	y1
)
{
	int	csm[3], nc;
	register BYTE	*cp;
	BYTE	rgb[3];
	int	gaps = 0;
	int	mx, my;
	register int	i;
					/* compute midpoint */
	mx = (x0 + x1) >> 1;
	my = (y0 + y1) >> 1;
					/* draw leaves */
	csm[0] = csm[1] = csm[2] = nc = 0;
	for (i = 0; i < 4; i++) {
		if (tp->flgs & LFF(i)) {
			cp = qtL.rgb[tp->k[i].li];
			csm[0] += cp[0]; csm[1] += cp[1]; csm[2] += cp[2];
			nc++;
			if (tp->flgs & CHF(i))
				cpaint(cp, qtL.wp[tp->k[i].li],
					i&01 ? mx : x0, i&02 ? my : y0,
					i&01 ? x1 : mx, i&02 ? y1 : my);
		} else if ((tp->flgs & CHBRF(i)) == CHF(i))
			gaps |= 1<<i;	/* empty stem */
	}
					/* do branches */
	for (i = 0; i < 4; i++)
		if ((tp->flgs & CHBRF(i)) == CHBRF(i)) {
			update(rgb, tp->k[i].b, i&01 ? mx : x0, i&02 ? my : y0,
					i&01 ? x1 : mx, i&02 ? y1 : my);
			csm[0] += rgb[0]; csm[1] += rgb[1]; csm[2] += rgb[2];
			nc++;
		}
	if (nc > 1) {
		ca[0] = csm[0]/nc; ca[1] = csm[1]/nc; ca[2] = csm[2]/nc;
	} else {
		ca[0] = csm[0]; ca[1] = csm[1]; ca[2] = csm[2];
	}
					/* fill in gaps with average */
	for (i = 0; gaps && i < 4; gaps >>= 1, i++)
		if (gaps & 01)
			cpaint(ca, (float *)NULL,
					i&01 ? mx : x0, i&02 ? my : y0,
					i&01 ? x1 : mx, i&02 ? y1 : my);
	tp->flgs &= ~CH_ANY;		/* all done */
}


extern void
qtRedraw(	/* redraw part or all of our screen */
	int	x0,
	int	y0,
	int	x1,
	int	y1
)
{
	int	lim[2][2];

	if (is_stump(&qtrunk))
		return;
	if (!(qtMapLeaves((lim[0][0]=x0) <= 0) & ((lim[1][0]=y0) <= 0) &
		((lim[0][1]=x1) >= odev.hres-1) & ((lim[1][1]=y1) >= odev.vres-1)))
		return;
	redraw(&qtrunk, 0, 0, odev.hres, odev.vres, lim);
}


extern void
qtUpdate(void)			/* update our tree display */
{
	BYTE	ca[3];

	if (is_stump(&qtrunk))
		return;
	if (!qtMapLeaves(0))
		return;
	update(ca, &qtrunk, 0, 0, odev.hres, odev.vres);
}
