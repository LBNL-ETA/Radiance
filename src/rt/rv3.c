#ifndef lint
static const char	RCSid[] = "$Id: rv3.c,v 2.13 2003/07/21 22:30:19 schorsch Exp $";
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
#ifdef SPEED
#define WFLUSH		(5*SPEED)
#else
#define WFLUSH		100		/* flush after this many rays */
#endif
#endif

#ifdef  SMLFLT
#define  sscanvec(s,v)	(sscanf(s,"%f %f %f",v,v+1,v+2)==3)
#else
#define  sscanvec(s,v)	(sscanf(s,"%lf %lf %lf",v,v+1,v+2)==3)
#endif


int
getrect(s, r)				/* get a box */
char  *s;
register RECT  *r;
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
getinterest(s, direc, vec, mp)		/* get area of interest */
char  *s;
int  direc;
FVECT  vec;
double  *mp;
{
	int  x, y;
	RAY  thisray;
	register int  i;

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
			rayorigin(&thisray, NULL, PRIMARY, 1.0);
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
	} else if (direc)
			for (i = 0; i < 3; i++)
				vec[i] -= ourview.vp[i];
	return(0);
}


float *		/* keep consistent with COLOR typedef */
greyof(col)				/* convert color to greyscale */
register COLOR  col;
{
	static COLOR  gcol;
	double  b;

	b = bright(col);
	setcolor(gcol, b, b, b);
	return(gcol);
}


void
paint(p, xmin, ymin, xmax, ymax)	/* compute and paint a rectangle */
register PNODE  *p;
int  xmin, ymin, xmax, ymax;
{
	static unsigned long  lastflush = 0;
	static RAY  thisray;
	double  h, v;

	if (xmax - xmin <= 0 || ymax - ymin <= 0) {	/* empty */
		p->x = xmin;
		p->y = ymin;
		setcolor(p->v, 0.0, 0.0, 0.0);
		return;
	}
						/* jitter ray direction */
	h = xmin + (xmax-xmin)*frandom();
	v = ymin + (ymax-ymin)*frandom();
	
	if ((thisray.rmax = viewray(thisray.rorg, thisray.rdir, &ourview,
			h/hresolu, v/vresolu)) < -FTINY) {
		setcolor(thisray.rcol, 0.0, 0.0, 0.0);
	} else {
		rayorigin(&thisray, NULL, PRIMARY, 1.0);
		samplendx++;
		rayvalue(&thisray);
	}

	p->x = h;
	p->y = v;
	copycolor(p->v, thisray.rcol);
	scalecolor(p->v, exposure);

	(*dev->paintr)(greyscale?greyof(p->v):p->v, xmin, ymin, xmax, ymax);

	if (dev->flush != NULL && nrays - lastflush >= WFLUSH) {
		lastflush = nrays;
		(*dev->flush)();
	}
}


void
newimage()				/* start a new image */
{
						/* free old image */
	freepkids(&ptrunk);
						/* save reserve memory */
	fillreserves();
						/* compute resolution */
	hresolu = dev->xsiz;
	vresolu = dev->ysiz;
	normaspect(viewaspect(&ourview), &dev->pixaspect, &hresolu, &vresolu);
	pframe.l = pframe.d = 0;
	pframe.r = hresolu; pframe.u = vresolu;
	pdepth = 0;
						/* clear device */
	(*dev->clear)(hresolu, vresolu);
						/* get first value */
	paint(&ptrunk, 0, 0, hresolu, vresolu);
}


void
redraw()				/* redraw the image */
{
	(*dev->clear)(hresolu, vresolu);
	(*dev->comout)("redrawing...\n");
	repaint(0, 0, hresolu, vresolu);
	(*dev->comout)("\n");
}


void
repaint(xmin, ymin, xmax, ymax)			/* repaint a region */
int  xmin, ymin, xmax, ymax;
{
	RECT  reg;

	reg.l = xmin; reg.r = xmax;
	reg.d = ymin; reg.u = ymax;

	paintrect(&ptrunk, 0, 0, hresolu, vresolu, &reg);
}


void
paintrect(p, xmin, ymin, xmax, ymax, r)		/* paint picture rectangle */
register PNODE  *p;
int  xmin, ymin, xmax, ymax;
register RECT  *r;
{
	int  mx, my;

	if (xmax - xmin <= 0 || ymax - ymin <= 0)
		return;

	if (p->kid == NULL) {
		(*dev->paintr)(greyscale?greyof(p->v):p->v,
				xmin, ymin, xmax, ymax);	/* do this */
		return;
	}
	mx = (xmin + xmax) >> 1;				/* do kids */
	my = (ymin + ymax) >> 1;
	if (mx > r->l) {
		if (my > r->d)
			paintrect(p->kid+DL, xmin, ymin, mx, my, r);
		if (my < r->u)
			paintrect(p->kid+UL, xmin, my, mx, ymax, r);
	}
	if (mx < r->r) {
		if (my > r->d)
			paintrect(p->kid+DR, mx, ymin, xmax, my, r);
		if (my < r->u)
			paintrect(p->kid+UR, mx, my, xmax, ymax, r);
	}
}


PNODE *
findrect(x, y, p, r, pd)		/* find a rectangle */
int  x, y;
register PNODE  *p;
register RECT  *r;
int  pd;
{
	int  mx, my;

	while (p->kid != NULL && pd--) {

		mx = (r->l + r->r) >> 1;
		my = (r->d + r->u) >> 1;

		if (x < mx) {
			r->r = mx;
			if (y < my) {
				r->u = my;
				p = p->kid+DL;
			} else {
				r->d = my;
				p = p->kid+UL;
			}
		} else {
			r->l = mx;
			if (y < my) {
				r->u = my;
				p = p->kid+DR;
			} else {
				r->d = my;
				p = p->kid+UR;
			}
		}
	}
	return(p);
}


void
scalepict(p, sf)			/* scale picture values */
register PNODE  *p;
double  sf;
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
getpictcolrs(yoff, scan, p, xsiz, ysiz)	/* get scanline from picture */
int  yoff;
register COLR  *scan;
register PNODE  *p;
int  xsiz, ysiz;
{
	register int  mx;
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
freepkids(p)				/* free pnode's children */
register PNODE  *p;
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
newview(vp)				/* change viewing parameters */
register VIEW  *vp;
{
	char  *err;

	if ((err = setview(vp)) != NULL) {
		sprintf(errmsg, "view not set - %s", err);
		error(COMMAND, errmsg);
	} else if (memcmp((char *)vp, (char *)&ourview, sizeof(VIEW))) {
		oldview = ourview;
		ourview = *vp;
		newimage();
	}
}


void
moveview(angle, elev, mag, vc)			/* move viewpoint */
double  angle, elev, mag;
FVECT  vc;
{
	double  d;
	FVECT  v1;
	VIEW  nv;
	register int  i;

	VCOPY(nv.vup, ourview.vup);
	nv.hoff = ourview.hoff; nv.voff = ourview.voff;
	spinvector(nv.vdir, ourview.vdir, ourview.vup, angle*(PI/180.));
	if (elev != 0.0) {
		fcross(v1, ourview.vup, nv.vdir);
		normalize(v1);
		spinvector(nv.vdir, nv.vdir, v1, elev*(PI/180.));
	}
	if ((nv.type = ourview.type) == VT_PAR) {
		nv.horiz = ourview.horiz / mag;
		nv.vert = ourview.vert / mag;
		d = 0.0;			/* don't move closer */
		for (i = 0; i < 3; i++)
			d += (vc[i] - ourview.vp[i])*ourview.vdir[i];
		nv.vfore = ourview.vfore;
		nv.vaft = ourview.vaft;
	} else {
		nv.horiz = ourview.horiz;
		nv.vert = ourview.vert;
		d = sqrt(dist2(ourview.vp, vc)) / mag;
		if ((nv.vfore = ourview.vfore) > FTINY) {
			nv.vfore += d - d*mag;
			if (nv.vfore < 0.0) nv.vfore = 0.0;
		}
		if ((nv.vaft = ourview.vaft) > FTINY) {
			nv.vaft += d - d*mag;
			if (nv.vaft <= nv.vfore) nv.vaft = 0.0;
		}
	}
	for (i = 0; i < 3; i++)
		nv.vp[i] = vc[i] - d*nv.vdir[i];
	newview(&nv);
}


void
pcopy(p1, p2)				/* copy paint node p1 into p2 */
register PNODE  *p1, *p2;
{
	copycolor(p2->v, p1->v);
	p2->x = p1->x;
	p2->y = p1->y;
}


void
zoomview(vp, zf)			/* zoom in or out */
register VIEW  *vp;
double  zf;
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
