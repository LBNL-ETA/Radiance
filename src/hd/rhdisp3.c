#ifndef lint
static const char	RCSid[] = "$Id: rhdisp3.c,v 3.13 2003/04/23 00:52:33 greg Exp $";
#endif
/*
 * Holodeck beam support for display process
 */

#include "rholo.h"
#include "rhdisp.h"
#include "view.h"

struct cellist {
	GCOORD	*cl;
	int	n;
};


int
npixels(vp, hr, vr, hp, bi)	/* compute appropriate nrays to evaluate */
register VIEW	*vp;
int	hr, vr;
HOLO	*hp;
int	bi;
{
	VIEW	vrev;
	GCOORD	gc[2];
	FVECT	cp[4], ip[4], pf, pb;
	double	af, ab, sf2, sb2, dfb2, df2, db2, penalty;
	register int	i;
					/* special case */
	if (hr <= 0 | vr <= 0)
		return(0);
					/* compute cell corners in image */
	if (!hdbcoord(gc, hp, bi))
		error(CONSISTENCY, "bad beam index in npixels");
	hdcell(cp, hp, gc+1);		/* find cell on front image */
	for (i = 3; i--; )		/* compute front center */
		pf[i] = 0.5*(cp[0][i] + cp[2][i]);
	sf2 = 0.25*dist2(cp[0], cp[2]);	/* compute half diagonal length */
	for (i = 0; i < 4; i++) {	/* compute visible quad */
		viewloc(ip[i], vp, cp[i]);
		if (ip[i][2] < 0.) {
			af = 0;
			goto getback;
		}
		ip[i][0] *= (double)hr;	/* scale by resolution */
		ip[i][1] *= (double)vr;
	}
					/* compute front area */
	af = (ip[1][0]-ip[0][0])*(ip[2][1]-ip[0][1]) -
		(ip[2][0]-ip[0][0])*(ip[1][1]-ip[0][1]);
	af += (ip[2][0]-ip[3][0])*(ip[1][1]-ip[3][1]) -
		(ip[1][0]-ip[3][0])*(ip[2][1]-ip[3][1]);
	af *= af >= 0 ? 0.5 : -0.5;
getback:
	copystruct(&vrev, vp);		/* compute reverse view */
	for (i = 0; i < 3; i++) {
		vrev.vdir[i] = -vp->vdir[i];
		vrev.vup[i] = -vp->vup[i];
		vrev.hvec[i] = -vp->hvec[i];
		vrev.vvec[i] = -vp->vvec[i];
	}
	hdcell(cp, hp, gc);		/* find cell on back image */
	for (i = 3; i--; )		/* compute rear center */
		pb[i] = 0.5*(cp[0][i] + cp[2][i]);
	sb2 = 0.25*dist2(cp[0], cp[2]);	/* compute half diagonal length */
	for (i = 0; i < 4; i++) {	/* compute visible quad */
		viewloc(ip[i], &vrev, cp[i]);
		if (ip[i][2] < 0.) {
			ab = 0;
			goto finish;
		}
		ip[i][0] *= (double)hr;	/* scale by resolution */
		ip[i][1] *= (double)vr;
	}
					/* compute back area */
	ab = (ip[1][0]-ip[0][0])*(ip[2][1]-ip[0][1]) -
		(ip[2][0]-ip[0][0])*(ip[1][1]-ip[0][1]);
	ab += (ip[2][0]-ip[3][0])*(ip[1][1]-ip[3][1]) -
		(ip[1][0]-ip[3][0])*(ip[2][1]-ip[3][1]);
	ab *= ab >= 0 ? 0.5 : -0.5;
finish:		/* compute penalty based on dist. sightline - viewpoint */
	df2 = dist2(vp->vp, pf);
	db2 = dist2(vp->vp, pb);
	dfb2 = dist2(pf, pb);
	penalty = dfb2 + df2 - db2;
	penalty = df2 - 0.25*penalty*penalty/dfb2;
	if (df2 > db2)	penalty /= df2 <= dfb2 ? sb2 : sb2*df2/dfb2;
	else		penalty /= db2 <= dfb2 ? sf2 : sf2*db2/dfb2;
	if (penalty < 1.) penalty = 1.;
					/* round off smaller non-zero area */
	if (ab <= FTINY || (af > FTINY && af <= ab))
		return((int)(af/penalty + 0.5));
	return((int)(ab/penalty + 0.5));
}


/*
 * The ray directions that define the pyramid in visit_cells() needn't
 * be normalized, but they must be given in clockwise order as seen
 * from the pyramid's apex (origin).
 * If no cell centers fall within the domain, the closest cell is visited.
 */
int
visit_cells(orig, pyrd, hp, vf, dp)	/* visit cells within a pyramid */
FVECT	orig, pyrd[4];		/* pyramid ray directions in clockwise order */
register HOLO	*hp;
int	(*vf)();
char	*dp;
{
	int	ncalls = 0, n = 0;
	int	inflags = 0;
	FVECT	gp, pn[4], lo, ld;
	double	po[4], lbeg, lend, d, t;
	GCOORD	gc, gc2[2];
	register int	i;
					/* figure out whose side we're on */
	hdgrid(gp, hp, orig);
	for (i = 0; i < 3; i++) {
		inflags |= (gp[i] > FTINY) << (i<<1);
		inflags |= (gp[i] < hp->grid[i]-FTINY) << (i<<1 | 1);
	}
					/* compute pyramid planes */
	for (i = 0; i < 4; i++) {
		fcross(pn[i], pyrd[i], pyrd[(i+1)&03]);
		po[i] = DOT(pn[i], orig);
	}
					/* traverse each wall */
	for (gc.w = 0; gc.w < 6; gc.w++) {
		if (!(inflags & 1<<gc.w))	/* origin on wrong side */
			continue;
					/* scanline algorithm */
		for (gc.i[1] = hp->grid[hdwg1[gc.w]]; gc.i[1]--; ) {
						/* compute scanline */
			gp[gc.w>>1] = gc.w&1 ? hp->grid[gc.w>>1] : 0;
			gp[hdwg0[gc.w]] = 0;
			gp[hdwg1[gc.w]] = gc.i[1] + 0.5;
			hdworld(lo, hp, gp);
			gp[hdwg0[gc.w]] = 1;
			hdworld(ld, hp, gp);
			ld[0] -= lo[0]; ld[1] -= lo[1]; ld[2] -= lo[2];
						/* find scanline limits */
			lbeg = 0; lend = hp->grid[hdwg0[gc.w]];
			for (i = 0; i < 4; i++) {
				t = DOT(pn[i], lo) - po[i];
				d = -DOT(pn[i], ld);
				if (d > FTINY) {		/* <- plane */
					if ((t /= d) < lend)
						lend = t;
				} else if (d < -FTINY) {	/* plane -> */
					if ((t /= d) > lbeg)
						lbeg = t;
				} else if (t < 0) {		/* outside */
					lend = -1;
					break;
				}
			}
			if (lbeg >= lend)
				continue;
			i = lend + .5;		/* visit cells on this scan */
			for (gc.i[0] = lbeg + .5; gc.i[0] < i; gc.i[0]++) {
				n += (*vf)(&gc, dp);
				ncalls++;
			}
		}
	}
	if (ncalls)			/* got one at least */
		return(n);
					/* else find closest cell */
	VSUM(ld, pyrd[0], pyrd[1], 1.);
	VSUM(ld, ld, pyrd[2], 1.);
	VSUM(ld, ld, pyrd[3], 1.);
#if 0
	if (normalize(ld) == 0.0)	/* technically not necessary */
		return(0);
#endif
	d = hdinter(gc2, NULL, &t, hp, orig, ld);
	if (d >= FHUGE || t <= 0.)
		return(0);
	return((*vf)(gc2+1, dp));	/* visit it */
}


sect_behind(hp, vp)		/* check if section is "behind" viewpoint */
register HOLO	*hp;
register VIEW	*vp;
{
	FVECT	hcent;
					/* compute holodeck section center */
	VSUM(hcent, hp->orig, hp->xv[0], 0.5);
	VSUM(hcent, hcent, hp->xv[1], 0.5);
	VSUM(hcent, hcent, hp->xv[2], 0.5);
					/* behind if center is behind */
	return(DOT(vp->vdir,hcent) < DOT(vp->vdir,vp->vp));
}


viewpyramid(org, dir, hp, vp)	/* compute view pyramid */
FVECT	org, dir[4];
HOLO	*hp;
VIEW	*vp;
{
	register int	i;
					/* check view type */
	if (vp->type == VT_PAR)
		return(0);
					/* in front or behind? */
	if (!sect_behind(hp, vp)) {
		if (viewray(org, dir[0], vp, 0., 0.) < -FTINY)
			return(0);
		if (viewray(org, dir[1], vp, 0., 1.) < -FTINY)
			return(0);
		if (viewray(org, dir[2], vp, 1., 1.) < -FTINY)
			return(0);
		if (viewray(org, dir[3], vp, 1., 0.) < -FTINY)
			return(0);
		return(1);
	}				/* reverse pyramid */
	if (viewray(org, dir[3], vp, 0., 0.) < -FTINY)
		return(0);
	if (viewray(org, dir[2], vp, 0., 1.) < -FTINY)
		return(0);
	if (viewray(org, dir[1], vp, 1., 1.) < -FTINY)
		return(0);
	if (viewray(org, dir[0], vp, 1., 0.) < -FTINY)
		return(0);
	for (i = 0; i < 3; i++) {
		dir[0][i] = -dir[0][i];
		dir[1][i] = -dir[1][i];
		dir[2][i] = -dir[2][i];
		dir[3][i] = -dir[3][i];
	}
	return(-1);
}


int
addcell(gcp, cl)		/* add a cell to a list */
GCOORD	*gcp;
register struct cellist	*cl;
{
	copystruct(cl->cl+cl->n, gcp);
	cl->n++;
	return(1);
}


int
cellcmp(gcp1, gcp2)		/* visit_cells() cell ordering */
register GCOORD	*gcp1, *gcp2;
{
	register int	c;

	if ((c = gcp1->w - gcp2->w))
		return(c);
	if ((c = gcp2->i[1] - gcp1->i[1]))	/* wg1 is reverse-ordered */
		return(c);
	return(gcp1->i[0] - gcp2->i[0]);
}


GCOORD *
getviewcells(np, hp, vp)	/* get ordered cell list for section view */
int	*np;		/* returned number of cells (negative if reversed) */
register HOLO	*hp;
VIEW	*vp;
{
	FVECT	org, dir[4];
	int	orient;
	struct cellist	cl;
					/* compute view pyramid */
	*np = 0;
	orient = viewpyramid(org, dir, hp, vp);
	if (!orient)
		return(NULL);
					/* allocate enough list space */
	cl.n = 2*(	hp->grid[0]*hp->grid[1] +
			hp->grid[0]*hp->grid[2] +
			hp->grid[1]*hp->grid[2]	);
	cl.cl = (GCOORD *)malloc(cl.n*sizeof(GCOORD));
	if (cl.cl == NULL)
		goto memerr;
	cl.n = 0;			/* add cells within pyramid */
	visit_cells(org, dir, hp, addcell, (char *)&cl);
	if (!cl.n) {
		free((void *)cl.cl);
		return(NULL);
	}
	*np = cl.n * orient;
#if 0
	/* We're just going to free this memory in a moment, and list is
	 * sorted automatically by visit_cells(), so we don't need this.
	 */
					/* optimize memory use */
	cl.cl = (GCOORD *)realloc((void *)cl.cl, cl.n*sizeof(GCOORD));
	if (cl.cl == NULL)
		goto memerr;
					/* sort the list */
	qsort((char *)cl.cl, cl.n, sizeof(GCOORD), cellcmp);
#endif
	return(cl.cl);
memerr:
	error(SYSTEM, "out of memory in getviewcells");
}


gridlines(f)			/* run through holodeck section grid lines */
int	(*f)();
{
	register int	hd, w, i;
	int	g0, g1;
	FVECT	wp[2], mov;
	double	d;
					/* do each wall on each section */
	for (hd = 0; hdlist[hd] != NULL; hd++)
		for (w = 0; w < 6; w++) {
			g0 = hdwg0[w];
			g1 = hdwg1[w];
			d = 1.0/hdlist[hd]->grid[g0];
			mov[0] = d * hdlist[hd]->xv[g0][0];
			mov[1] = d * hdlist[hd]->xv[g0][1];
			mov[2] = d * hdlist[hd]->xv[g0][2];
			if (w & 1) {
				VSUM(wp[0], hdlist[hd]->orig,
						hdlist[hd]->xv[w>>1], 1.);
				VSUM(wp[0], wp[0], mov, 1.);
			} else
				VCOPY(wp[0], hdlist[hd]->orig);
			VSUM(wp[1], wp[0], hdlist[hd]->xv[g1], 1.);
			for (i = hdlist[hd]->grid[g0]; ; ) {	/* g0 lines */
				(*f)(wp);
				if (!--i) break;
				wp[0][0] += mov[0]; wp[0][1] += mov[1];
				wp[0][2] += mov[2]; wp[1][0] += mov[0];
				wp[1][1] += mov[1]; wp[1][2] += mov[2];
			}
			d = 1.0/hdlist[hd]->grid[g1];
			mov[0] = d * hdlist[hd]->xv[g1][0];
			mov[1] = d * hdlist[hd]->xv[g1][1];
			mov[2] = d * hdlist[hd]->xv[g1][2];
			if (w & 1)
				VSUM(wp[0], hdlist[hd]->orig,
						hdlist[hd]->xv[w>>1], 1.);
			else
				VSUM(wp[0], hdlist[hd]->orig, mov, 1.);
			VSUM(wp[1], wp[0], hdlist[hd]->xv[g0], 1.);
			for (i = hdlist[hd]->grid[g1]; ; ) {	/* g1 lines */
				(*f)(wp);
				if (!--i) break;
				wp[0][0] += mov[0]; wp[0][1] += mov[1];
				wp[0][2] += mov[2]; wp[1][0] += mov[0];
				wp[1][1] += mov[1]; wp[1][2] += mov[2];
			}
		}
}
