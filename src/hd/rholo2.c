#ifndef lint
static const char	RCSid[] = "$Id: rholo2.c,v 3.27 2003/07/27 22:12:02 schorsch Exp $";
#endif
/*
 * Rtrace support routines for holodeck rendering
 */

#include "rholo.h"
#include "paths.h"
#include "random.h"


VIEWPOINT	myeye;		/* target view position */

struct gclim {
	HOLO	*hp;			/* holodeck pointer */
	GCOORD	gc;			/* grid cell */
	FVECT	egp;			/* eye grid point */
	double	erg2;			/* mean square eye grid range */
	double	gmin[2], gmax[2];	/* grid coordinate limits */
};				/* a grid coordinate range */


static
initeyelim(gcl, hp, gc)		/* initialize grid coordinate limits */
register struct gclim	*gcl;
register HOLO	*hp;
GCOORD	*gc;
{
	register RREAL	*v;
	register int	i;

	if (hp != NULL) {
		hdgrid(gcl->egp, gcl->hp = hp, myeye.vpt);
		gcl->erg2 = 0;
		for (i = 0, v = hp->wg[0]; i < 3; i++, v += 3)
			gcl->erg2 += DOT(v,v);
		gcl->erg2 *= (1./3.) * myeye.rng*myeye.rng;
	}
	if (gc != NULL)
		gcl->gc = *gc;
	gcl->gmin[0] = gcl->gmin[1] = FHUGE;
	gcl->gmax[0] = gcl->gmax[1] = -FHUGE;
}


static
groweyelim(gcl, gc, r0, r1, tight)	/* grow grid limits about eye point */
register struct gclim	*gcl;
GCOORD	*gc;
double	r0, r1;
int	tight;
{
	FVECT	gp, ab;
	double	ab2, od, cfact;
	double	sqcoef[3], ctcoef[3], licoef[3], cnst;
	int	gw, gi[2];
	double	wallpos, a, b, c, d, e, f;
	double	root[2], yex;
	int	n, i, j, nex;
						/* point/view cone */
	i = gc->w>>1;
	gp[i] = gc->w&1 ? gcl->hp->grid[i] : 0;
	gp[hdwg0[gc->w]] = gc->i[0] + r0;
	gp[hdwg1[gc->w]] = gc->i[1] + r1;
	VSUB(ab, gcl->egp, gp);
	ab2 = DOT(ab, ab);
	gw = gcl->gc.w>>1;
	if ((i==gw ? ab[gw]*ab[gw] : ab2)  <= gcl->erg2 + FTINY) {
		gcl->gmin[0] = gcl->gmin[1] = -FHUGE;
		gcl->gmax[0] = gcl->gmax[1] = FHUGE;
		return;			/* too close (to wall) */
	}
	ab2 = 1./ab2;				/* 1/norm2(ab) */
	od = DOT(gp, ab);			/* origin dot direction */
	cfact = 1./(1. - ab2*gcl->erg2);	/* tan^2 + 1 of cone angle */
	for (i = 0; i < 3; i++) {		/* compute cone equation */
		sqcoef[i] = ab[i]*ab[i]*cfact*ab2 - 1.;
		ctcoef[i] = 2.*ab[i]*ab[(i+1)%3]*cfact*ab2;
		licoef[i] = 2.*(gp[i] - ab[i]*cfact*od*ab2);
	}
	cnst = cfact*od*od*ab2 - DOT(gp,gp);
	/*
	 * CONE:	sqcoef[0]*x*x + sqcoef[1]*y*y + sqcoef[2]*z*z
	 *		+ ctcoef[0]*x*y + ctcoef[1]*y*z + ctcoef[2]*z*x
	 *		+ licoef[0]*x + licoef[1]*y + licoef[2]*z + cnst == 0
	 */
				/* equation for conic section in plane */
	gi[0] = hdwg0[gcl->gc.w];
	gi[1] = hdwg1[gcl->gc.w];
	wallpos = gcl->gc.w&1 ? gcl->hp->grid[gw] : 0;
	a = sqcoef[gi[0]];					/* x2 */
	b = ctcoef[gi[0]];					/* xy */
	c = sqcoef[gi[1]];					/* y2 */
	d = ctcoef[gw]*wallpos + licoef[gi[0]];			/* x */
	e = ctcoef[gi[1]]*wallpos + licoef[gi[1]];		/* y */
	f = wallpos*(wallpos*sqcoef[gw] + licoef[gw]) + cnst;
	for (i = 0; i < 2; i++) {
		if (i) {		/* swap x and y coefficients */
			register double	t;
			t = a; a = c; c = t;
			t = d; d = e; e = t;
		}
		nex = 0;		/* check global extrema */
		n = quadratic(root, a*(4.*a*c-b*b), 2.*a*(2.*c*d-b*e),
				d*(c*d-b*e) + f*b*b);
		while (n-- > 0) {
			if (gc->w>>1 == gi[i] &&
					(gc->w&1) ^ (root[n] < gp[gc->w>>1])) {
				if (gc->w&1)
					gcl->gmin[i] = -FHUGE;
				else
					gcl->gmax[i] = FHUGE;
				nex++;
				continue;		/* hyperbolic */
			}
			if (tight) {
				yex = (-2.*a*root[n] - d)/b;
				if (yex < gcl->gc.i[1-i] ||
						yex > gcl->gc.i[1-i]+1)
					continue;	/* outside cell */
			}
			if (root[n] < gcl->gmin[i])
				gcl->gmin[i] = root[n];
			if (root[n] > gcl->gmax[i])
				gcl->gmax[i] = root[n];
			nex++;
		}
					/* check local extrema */
		for (j = nex < 2 ? 2 : 0; j--; ) {
			yex = gcl->gc.i[1-i] + j;
			n = quadratic(root, a, b*yex+d, yex*(yex*c+e)+f);
			while (n-- > 0) {
				if (gc->w>>1 == gi[i] &&
					(gc->w&1) ^ (root[n] < gp[gc->w>>1]))
					continue;
				if (root[n] < gcl->gmin[i])
					gcl->gmin[i] = root[n];
				if (root[n] > gcl->gmax[i])
					gcl->gmax[i] = root[n];
			}
		}
	}
}


static int
clipeyelim(rrng, gcl)		/* clip eye limits to grid cell */
register short	rrng[2][2];
register struct gclim	*gcl;
{
	int	incell = 1;
	register int	i;

	for (i = 0; i < 2; i++) {
		if (gcl->gmin[i] < gcl->gc.i[i])
			gcl->gmin[i] = gcl->gc.i[i];
		if (gcl->gmax[i] > gcl->gc.i[i]+1)
			gcl->gmax[i] = gcl->gc.i[i]+1;
		if (gcl->gmax[i] > gcl->gmin[i]) {
			rrng[i][0] = 256.*(gcl->gmin[i] - gcl->gc.i[i]) +
					(1.-FTINY);
			rrng[i][1] = 256.*(gcl->gmax[i] - gcl->gc.i[i]) +
					(1.-FTINY) - rrng[i][0];
		} else
			rrng[i][0] = rrng[i][1] = 0;
		incell &= rrng[i][1] > 0;
	}
	return(incell);
}


packrays(rod, p)		/* pack ray origins and directions */
register float	*rod;
register PACKET	*p;
{
#if 0
	double	dist2sum = 0.;
	FVECT	vt;
#endif
	int	nretries = p->nr + 2;
	struct gclim	eyelim;
	short	rrng0[2][2], rrng1[2][2];
	int	useyelim;
	GCOORD	gc[2];
	FVECT	ro, rd;
	double	d;
	register int	i;

	if (!hdbcoord(gc, hdlist[p->hd], p->bi))
		error(CONSISTENCY, "bad beam index in packrays");
	if ((useyelim = myeye.rng > FTINY)) {
		initeyelim(&eyelim, hdlist[p->hd], gc);
		groweyelim(&eyelim, gc+1, 0., 0., 0);
		groweyelim(&eyelim, gc+1, 1., 1., 0);
		useyelim = clipeyelim(rrng0, &eyelim);
#ifdef DEBUG
		if (!useyelim)
			error(WARNING, "no eye overlap in packrays");
#endif
	}
	for (i = 0; i < p->nr; i++) {
	retry:
		if (useyelim) {
			initeyelim(&eyelim, NULL, gc+1);
			p->ra[i].r[0][0] = (int)(frandom()*rrng0[0][1])
						+ rrng0[0][0];
			p->ra[i].r[0][1] = (int)(frandom()*rrng0[1][1])
						+ rrng0[1][0];
			groweyelim(&eyelim, gc,
					(1./256.)*(p->ra[i].r[0][0]+.5),
					(1./256.)*(p->ra[i].r[0][1]+.5), 1);
			if (!clipeyelim(rrng1, &eyelim)) {
				useyelim = nretries-- > 0;
#ifdef DEBUG
				if (!useyelim)
					error(WARNING,
					"exceeded retry limit in packrays");
#endif
				goto retry;
			}
			p->ra[i].r[1][0] = (int)(frandom()*rrng1[0][1])
						+ rrng1[0][0];
			p->ra[i].r[1][1] = (int)(frandom()*rrng1[1][1])
						+ rrng1[1][0];
		} else {
			p->ra[i].r[0][0] = frandom() * 256.;
			p->ra[i].r[0][1] = frandom() * 256.;
			p->ra[i].r[1][0] = frandom() * 256.;
			p->ra[i].r[1][1] = frandom() * 256.;
		}
		d = hdray(ro, rd, hdlist[p->hd], gc, p->ra[i].r);
#if 0
		VSUM(vt, ro, rd, d);
		dist2sum += dist2line(myeye.vpt, ro, vt);
#endif
		if (p->offset != NULL) {
			if (!vdef(OBSTRUCTIONS))
				d *= frandom();		/* random offset */
			VSUM(ro, ro, rd, d);		/* advance ray */
			p->offset[i] = d;
		}
		VCOPY(rod, ro);
		rod += 3;
		VCOPY(rod, rd);
		rod += 3;
	}
#if 0
	fprintf(stderr, "%f RMS (%d retries)\t", sqrt(dist2sum/p->nr),
			p->nr + 2 - nretries);
#endif
}


donerays(p, rvl)		/* encode finished ray computations */
register PACKET	*p;
register float	*rvl;
{
	double	d;
	register int	i;

	for (i = 0; i < p->nr; i++) {
		setcolr(p->ra[i].v, rvl[0], rvl[1], rvl[2]);
		d = rvl[3];
		if (p->offset != NULL)
			d += p->offset[i];
		p->ra[i].d = hdcode(hdlist[p->hd], d);
		rvl += 4;
	}
	p->nc += p->nr;
}


int
done_rtrace()			/* clean up and close rtrace calculation */
{
	int	status;
					/* already closed? */
	if (!nprocs)
		return(0);
					/* flush beam queue */
	done_packets(flush_queue());
					/* sync holodeck */
	hdsync(NULL, 1);
					/* close rtrace */
	if ((status = end_rtrace()))
		error(WARNING, "bad exit status from rtrace");
	if (vdef(REPORT)) {		/* report time */
		eputs("rtrace process closed\n");
		report(0);
	}
	return(status);			/* return status */
}


new_rtrace()			/* restart rtrace calculation */
{
	char	combuf[128];

	if (nprocs > 0)			/* already running? */
		return;
	starttime = time(NULL);		/* reset start time and counts */
	npacksdone = nraysdone = 0L;
	if (vdef(TIME))			/* reset end time */
		endtime = starttime + vflt(TIME)*3600. + .5;
	if (vdef(RIF)) {		/* rerun rad to update octree */
		sprintf(combuf, "rad -v 0 -s -w %s", vval(RIF));
		if (system(combuf))
			error(WARNING, "error running rad");
	}
	if (start_rtrace() < 1)		/* start rtrace */
		error(WARNING, "cannot restart rtrace");
	else if (vdef(REPORT)) {
		eputs("rtrace process restarted\n");
		report(0);
	}
}


getradfile()			/* run rad and get needed variables */
{
	static short	mvar[] = {OCTREE,EYESEP,-1};
	static char	tf1[] = TEMPLATE;
	char	tf2[64];
	char	combuf[256];
	char	*pippt;
	register int	i;
	register char	*cp;
					/* check if rad file specified */
	if (!vdef(RIF))
		return(0);
					/* create rad command */
	mktemp(tf1);
	sprintf(tf2, "%s.rif", tf1);
	sprintf(combuf,
		"rad -v 0 -s -e -w %s OPTFILE=%s | egrep '^[ \t]*(NOMATCH",
			vval(RIF), tf1);
	cp = combuf;
	while (*cp){
		if (*cp == '|') pippt = cp;
		cp++;
	}				/* match unset variables */
	for (i = 0; mvar[i] >= 0; i++)
		if (!vdef(mvar[i])) {
			*cp++ = '|';
			strcpy(cp, vnam(mvar[i]));
			while (*cp) cp++;
			pippt = NULL;
		}
	if (pippt != NULL)
		strcpy(pippt, "> " NULL_DEVICE);	/* nothing to match */
	else
		sprintf(cp, ")[ \t]*=' > %s", tf2);
#ifdef DEBUG
	wputs(combuf); wputs("\n");
#endif
	system(combuf);				/* ignore exit code */
	if (pippt == NULL) {
		loadvars(tf2);			/* load variables */
		unlink(tf2);
	}
	rtargc += wordfile(rtargv+rtargc, tf1);	/* get rtrace options */
	unlink(tf1);			/* clean up */
	return(1);
}


report(t)			/* report progress so far */
time_t	t;
{
	static time_t	seconds2go = 1000000;

	if (t == 0L)
		t = time(NULL);
	sprintf(errmsg, "%ld packets (%ld rays) done after %.2f hours\n",
			npacksdone, nraysdone, (t-starttime)/3600.);
	eputs(errmsg);
	if (seconds2go == 1000000)
		seconds2go = vdef(REPORT) ? (long)(vflt(REPORT)*60. + .5) : 0L;
	if (seconds2go)
		reporttime = t + seconds2go;
}
