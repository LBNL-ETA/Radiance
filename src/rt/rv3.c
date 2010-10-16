#ifndef lint
static const char	RCSid[] = "$Id: rv3.c,v 2.35 2010/10/16 15:08:14 greg Exp $";
#endif
/*
 *  rv3.c - miscellaneous routines for rview.
 *
 *  External symbols declared in rpaint.h
 */

#include "copyright.h"

#include <string.h>

#include  "ray.h"
#include  "rpaint.h"
#include  "random.h"

#ifndef WFLUSH
#define WFLUSH		64		/* flush after this many primary rays */
#endif
#ifndef WFLUSH1
#define WFLUSH1		512		/* or this many total rays */
#endif


#ifdef  SMLFLT
#define  sscanvec(s,v)	(sscanf(s,"%f %f %f",v,v+1,v+2)==3)
#else
#define  sscanvec(s,v)	(sscanf(s,"%lf %lf %lf",v,v+1,v+2)==3)
#endif

static RNUMBER  niflush;		/* flushes since newimage() */

int
getrect(				/* get a box */
	char  *s,
	RECT  *r
)
{
	int  x0, y0, x1, y1;

	if (*s && !strncmp(s, "all", strlen(s))) {
		r->l = r->d = 0;
		r->r = hresolu;
		r->u = vresolu;
		return(0);
	}
	if (sscanf(s, "%d %d %d %d", &x0, &y0, &x1, &y1) != 4) {
		if (dev->getcur == NULL)
			return(-1);
		(*dev->comout)("Pick first corner\n");
		if ((*dev->getcur)(&x0, &y0) == ABORT)
			return(-1);
		(*dev->comout)("Pick second corner\n");
		if ((*dev->getcur)(&x1, &y1) == ABORT)
			return(-1);
	}
	if (x0 < x1) {
		r->l = x0;
		r->r = x1;
	} else {
		r->l = x1;
		r->r = x0;
	}
	if (y0 < y1) {
		r->d = y0;
		r->u = y1;
	} else {
		r->d = y1;
		r->u = y0;
	}
	if (r->l < 0) r->l = 0;
	if (r->d < 0) r->d = 0;
	if (r->r > hresolu) r->r = hresolu;
	if (r->u > vresolu) r->u = vresolu;
	if (r->l > r->r) r->l = r->r;
	if (r->d > r->u) r->d = r->u;
	return(0);
}


int
getinterest(		/* get area of interest */
	char  *s,
	int  direc,
	FVECT  vec,
	double  *mp
)
{
	int  x, y;
	RAY  thisray;
	int  i;

	if (sscanf(s, "%lf", mp) != 1)
		*mp = 1.0;
	else if (*mp < -FTINY)		/* negative zoom is reduction */
		*mp = -1.0 / *mp;
	else if (*mp <= FTINY) {	/* too small */
		error(COMMAND, "illegal magnification");
		return(-1);
	}
	if (!sscanvec(sskip(s), vec)) {
		if (dev->getcur == NULL)
			return(-1);
		(*dev->comout)("Pick view center\n");
		if ((*dev->getcur)(&x, &y) == ABORT)
			return(-1);
		if ((thisray.rmax = viewray(thisray.rorg, thisray.rdir,
			&ourview, (x+.5)/hresolu, (y+.5)/vresolu)) < -FTINY) {
			error(COMMAND, "not on image");
			return(-1);
		}
		if (!direc || ourview.type == VT_PAR) {
			rayorigin(&thisray, PRIMARY, NULL, NULL);
			if (!localhit(&thisray, &thescene)) {
				error(COMMAND, "not a local object");
				return(-1);
			}
		}
		if (direc)
			if (ourview.type == VT_PAR)
				for (i = 0; i < 3; i++)
					vec[i] = thisray.rop[i] - ourview.vp[i];
			else
				VCOPY(vec, thisray.rdir);
		else
			VCOPY(vec, thisray.rop);
	} else if (direc) {
		for (i = 0; i < 3; i++)
			vec[i] -= ourview.vp[i];
		if (normalize(vec) == 0.0) {
			error(COMMAND, "point at view origin");
			return(-1);
		}
	}
	return(0);
}


float *		/* keep consistent with COLOR typedef */
greyof(				/* convert color to greyscale */
	COLOR  col
)
{
	static COLOR  gcol;
	double  b;

	b = bright(col);
	setcolor(gcol, b, b, b);
	return(gcol);
}

static void
recolor(					/* recolor the given node */
	PNODE *p
)
{
	while (p->kid != NULL) {		/* need to propogate down */
		int  mx = (p->xmin + p->xmax) >> 1;
		int  my = (p->ymin + p->ymax) >> 1;
		int  ki;
		if (p->x >= mx)
			ki = (p->y >= my) ? UR : DR;
		else
			ki = (p->y >= my) ? UL : DL;
		pcopy(p, p->kid+ki);
		p = p->kid + ki;
	}

	(*dev->paintr)(greyscale?greyof(p->v):p->v,
			p->xmin, p->ymin, p->xmax, p->ymax);
}

int
paint(			/* compute and paint a rectangle */
	PNODE  *p
)
{
	static RAY  thisray;
	double  h, v;

	if ((p->xmax <= p->xmin) | (p->ymax <= p->ymin)) {	/* empty */
		p->x = p->xmin;
		p->y = p->ymin;
		setcolor(p->v, 0.0, 0.0, 0.0);
		return(0);
	}
						/* jitter ray direction */
	p->x = h = p->xmin + (p->xmax-p->xmin)*frandom();
	p->y = v = p->ymin + (p->ymax-p->ymin)*frandom();
	
	if ((thisray.rmax = viewray(thisray.rorg, thisray.rdir, &ourview,
			h/hresolu, v/vresolu)) < -FTINY) {
		setcolor(thisray.rcol, 0.0, 0.0, 0.0);
	} else if (nproc == 1) {		/* immediate mode */
		ray_trace(&thisray);
	} else {				/* queuing mode */
		int	rval;
		rayorigin(&thisray, PRIMARY, NULL, NULL);
		thisray.rno = (RNUMBER)p;
		rval = ray_pqueue(&thisray);
		if (!rval)
			return(0);
		if (rval < 0)
			return(-1);
						/* get node for returned ray */
		p = (PNODE *)thisray.rno;
	}

	copycolor(p->v, thisray.rcol);
	scalecolor(p->v, exposure);

	recolor(p);				/* paint it */

	if (dev->flush != NULL) {		/* shall we check for input? */
		static RNUMBER	lastflush = 0;
		RNUMBER		counter = raynum;
		int		flushintvl;
		if (nproc == 1) {
			counter = nrays;
			flushintvl = WFLUSH1;
		} else if (ambounce == 0)
			flushintvl = nproc*WFLUSH;
		else if (niflush < WFLUSH)
			flushintvl = nproc*niflush/(ambounce+1);
		else
			flushintvl = nproc*WFLUSH/(ambounce+1);
		if (lastflush > counter)
			lastflush = 0;		/* counter wrapped */

		if (counter - lastflush >= flushintvl) {
			lastflush = counter;
			(*dev->flush)();
			niflush++;
		}
	}
	return(1);
}


int
waitrays(void)					/* finish up pending rays */
{
	int	nwaited = 0;
	int	rval;
	RAY	raydone;

	if (nproc <= 1)				/* immediate mode? */
		return(0);
	while ((rval = ray_presult(&raydone, 0)) > 0) {
		PNODE  *p = (PNODE *)raydone.rno;
		copycolor(p->v, raydone.rcol);
		scalecolor(p->v, exposure);
		recolor(p);
		nwaited++;
	}
	if (rval < 0)
		return(-1);
	return(nwaited);
}


void
newimage(					/* start a new image */
	char *s
)
{
	extern int	ray_pnprocs;
	int		newnp = 0;
						/* # rendering procs arg? */
	if (s != NULL)
		sscanf(s, "%d", &newnp);
						/* free old image */
	freepkids(&ptrunk);
						/* compute resolution */
	hresolu = dev->xsiz;
	vresolu = dev->ysiz;
	normaspect(viewaspect(&ourview), &dev->pixaspect, &hresolu, &vresolu);
	ptrunk.xmin = ptrunk.ymin = pframe.l = pframe.d = 0;
	ptrunk.xmax = pframe.r = hresolu;
	ptrunk.ymax = pframe.u = vresolu;
	pdepth = 0;
						/* clear device */
	(*dev->clear)(hresolu, vresolu);

	if (newparam) {				/* (re)start rendering procs */
		if (ray_pnprocs > 0)
			ray_pclose(0);
		if (newnp > 0)
			nproc = newnp;
		if (nproc > 1)
			ray_popen(nproc);
		newparam = 0;
	} else if ((newnp > 0) & (newnp != nproc)) {
		if (newnp == 1)			/* change # rendering procs */
			ray_pclose(0);
		else if (newnp < ray_pnprocs)
			ray_pclose(ray_pnprocs - newnp);
		else
			ray_popen(newnp - ray_pnprocs);
		nproc = newnp;
	}
	niflush = 0;				/* get first value */
	paint(&ptrunk);
}


void
redraw(void)				/* redraw the image */
{
	(*dev->clear)(hresolu, vresolu);
	(*dev->comout)("redrawing...\n");
	repaint(0, 0, hresolu, vresolu);
	(*dev->comout)("\n");
}


void
repaint(				/* repaint a region */
	int  xmin,
	int  ymin,
	int  xmax,
	int  ymax
)
{
	RECT  reg;

	reg.l = xmin; reg.r = xmax;
	reg.d = ymin; reg.u = ymax;

	paintrect(&ptrunk, &reg);
}


void
paintrect(				/* paint picture rectangle */
	PNODE  *p,
	RECT  *r
)
{
	int  mx, my;

	if (p->xmax - p->xmin <= 0 || p->ymax - p->ymin <= 0)
		return;

	if (p->kid == NULL) {
		(*dev->paintr)(greyscale?greyof(p->v):p->v,
			p->xmin, p->ymin, p->xmax, p->ymax);	/* do this */
		return;
	}
	mx = (p->xmin + p->xmax) >> 1;				/* do kids */
	my = (p->ymin + p->ymax) >> 1;
	if (mx > r->l) {
		if (my > r->d)
			paintrect(p->kid+DL, r);
		if (my < r->u)
			paintrect(p->kid+UL, r);
	}
	if (mx < r->r) {
		if (my > r->d)
			paintrect(p->kid+DR, r);
		if (my < r->u)
			paintrect(p->kid+UR, r);
	}
}


PNODE *
findrect(				/* find a rectangle */
	int  x,
	int  y,
	PNODE  *p,
	int  pd
)
{
	int  mx, my;

	while (p->kid != NULL && pd--) {

		mx = (p->xmin + p->xmax) >> 1;
		my = (p->ymin + p->ymax) >> 1;

		if (x < mx) {
			if (y < my) {
				p = p->kid+DL;
			} else {
				p = p->kid+UL;
			}
		} else {
			if (y < my) {
				p = p->kid+DR;
			} else {
				p = p->kid+UR;
			}
		}
	}
	return(p);
}


void
compavg(				/* recompute averages */
	PNODE	*p
)
{
	int	i, navg;
	
	if (p->kid == NULL)
		return;

	setcolor(p->v, .0, .0, .0);
	navg = 0;
	for (i = 0; i < 4; i++) {
		if (p->kid[i].xmin >= p->kid[i].xmax) continue;
		if (p->kid[i].ymin >= p->kid[i].ymax) continue;
		compavg(p->kid+i);
		addcolor(p->v, p->kid[i].v);
		navg++;
	}
	if (navg > 1)
		scalecolor(p->v, 1./navg);
}


void
scalepict(				/* scale picture values */
	PNODE  *p,
	double  sf
)
{
	scalecolor(p->v, sf);		/* do this node */

	if (p->kid == NULL)
		return;
					/* do children */
	scalepict(p->kid+DL, sf);
	scalepict(p->kid+DR, sf);
	scalepict(p->kid+UL, sf);
	scalepict(p->kid+UR, sf);
}


void
getpictcolrs(				/* get scanline from picture */
	int  yoff,
	COLR  *scan,
	PNODE  *p,
	int  xsiz,
	int  ysiz
)
{
	int  mx;
	int  my;

	if (p->kid == NULL) {			/* do this node */
		setcolr(scan[0], colval(p->v,RED),
				colval(p->v,GRN),
				colval(p->v,BLU));
		for (mx = 1; mx < xsiz; mx++)
			copycolr(scan[mx], scan[0]);
		return;
	}
						/* do kids */
	mx = xsiz >> 1;
	my = ysiz >> 1;
	if (yoff < my) {
		getpictcolrs(yoff, scan, p->kid+DL, mx, my);
		getpictcolrs(yoff, scan+mx, p->kid+DR, xsiz-mx, my);
	} else {
		getpictcolrs(yoff-my, scan, p->kid+UL, mx, ysiz-my);
		getpictcolrs(yoff-my, scan+mx, p->kid+UR, xsiz-mx, ysiz-my);
	}
}


void
freepkids(				/* free pnode's children */
	PNODE  *p
)
{
	if (p->kid == NULL)
		return;
	freepkids(p->kid+DL);
	freepkids(p->kid+DR);
	freepkids(p->kid+UL);
	freepkids(p->kid+UR);
	free((void *)p->kid);
	p->kid = NULL;
}


void
newview(					/* change viewing parameters */
	VIEW  *vp
)
{
	char  *err;

	if ((err = setview(vp)) != NULL) {
		sprintf(errmsg, "view not set - %s", err);
		error(COMMAND, errmsg);
	} else if (memcmp((char *)vp, (char *)&ourview, sizeof(VIEW))) {
		oldview = ourview;
		ourview = *vp;
		newimage(NULL);
	}
}


void
moveview(					/* move viewpoint */
	double  angle,
	double  elev,
	double  mag,
	FVECT  vc
)
{
	double  d;
	FVECT  v1;
	VIEW  nv = ourview;
	int  i;

	spinvector(nv.vdir, ourview.vdir, ourview.vup, angle*(PI/180.));
	if (elev != 0.0) {
		fcross(v1, ourview.vup, nv.vdir);
		normalize(v1);
		spinvector(nv.vdir, nv.vdir, v1, elev*(PI/180.));
	}
	if (nv.type == VT_PAR) {
		nv.horiz /= mag;
		nv.vert /= mag;
		d = 0.0;			/* don't move closer */
		for (i = 0; i < 3; i++)
			d += (vc[i] - ourview.vp[i])*ourview.vdir[i];
	} else {
		d = sqrt(dist2(ourview.vp, vc)) / mag;
		if (nv.vfore > FTINY) {
			nv.vfore += d - d*mag;
			if (nv.vfore < 0.0) nv.vfore = 0.0;
		}
		if (nv.vaft > FTINY) {
			nv.vaft += d - d*mag;
			if (nv.vaft <= nv.vfore) nv.vaft = 0.0;
		}
		nv.vdist /= mag;
	}
	for (i = 0; i < 3; i++)
		nv.vp[i] = vc[i] - d*nv.vdir[i];
	newview(&nv);
}


void
pcopy(					/* copy paint node p1 into p2 */
	PNODE  *p1,
	PNODE  *p2
)
{
	copycolor(p2->v, p1->v);
	p2->x = p1->x;
	p2->y = p1->y;
}


void
zoomview(				/* zoom in or out */
	VIEW  *vp,
	double  zf
)
{
	switch (vp->type) {
	case VT_PAR:				/* parallel view */
		vp->horiz /= zf;
		vp->vert /= zf;
		return;
	case VT_ANG:				/* angular fisheye */
		vp->horiz /= zf;
		if (vp->horiz > 360.)
			vp->horiz = 360.;
		vp->vert /= zf;
		if (vp->vert > 360.)
			vp->vert = 360.;
		return;
	case VT_PLS:				/* planisphere fisheye */
		vp->horiz = sin((PI/180./2.)*vp->horiz) /
				(1.0 + cos((PI/180./2.)*vp->horiz)) / zf;
		vp->horiz *= vp->horiz;
		vp->horiz = (2.*180./PI)*acos((1. - vp->horiz) /
						(1. + vp->horiz));
		vp->vert = sin((PI/180./2.)*vp->vert) /
				(1.0 + cos((PI/180./2.)*vp->vert)) / zf;
		vp->vert *= vp->vert;
		vp->vert = (2.*180./PI)*acos((1. - vp->vert) /
						(1. + vp->vert));
		return;
	case VT_CYL:				/* cylindrical panorama */
		vp->horiz /= zf;
		if (vp->horiz > 360.)
			vp->horiz = 360.;
		vp->vert = atan(tan(vp->vert*(PI/180./2.))/zf) / (PI/180./2.);
		return;
	case VT_PER:				/* perspective view */
		vp->horiz = atan(tan(vp->horiz*(PI/180./2.))/zf) /
				(PI/180./2.);
		vp->vert = atan(tan(vp->vert*(PI/180./2.))/zf) /
				(PI/180./2.);
		return;
	case VT_HEM:				/* hemispherical fisheye */
		vp->horiz = sin(vp->horiz*(PI/180./2.))/zf;
		if (vp->horiz >= 1.0-FTINY)
			vp->horiz = 180.;
		else
			vp->horiz = asin(vp->horiz) / (PI/180./2.);
		vp->vert = sin(vp->vert*(PI/180./2.))/zf;
		if (vp->vert >= 1.0-FTINY)
			vp->vert = 180.;
		else
			vp->vert = asin(vp->vert) / (PI/180./2.);
		return;
	}
}
