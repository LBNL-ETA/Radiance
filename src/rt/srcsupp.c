#ifndef lint
static const char	RCSid[] = "$Id: srcsupp.c,v 2.11 2003/04/23 00:52:34 greg Exp $";
#endif
/*
 *  Support routines for source objects and materials
 *
 *  External symbols declared in source.h
 */

#include "copyright.h"

#include  "ray.h"

#include  "otypes.h"

#include  "source.h"

#include  "cone.h"

#include  "face.h"

#define SRCINC		4		/* realloc increment for array */

SRCREC  *source = NULL;			/* our list of sources */
int  nsources = 0;			/* the number of sources */

SRCFUNC  sfun[NUMOTYPE];		/* source dispatch table */


void
initstypes()			/* initialize source dispatch table */
{
	extern VSMATERIAL  mirror_vs, direct1_vs, direct2_vs;
	static SOBJECT  fsobj = {fsetsrc, flatpart, fgetplaneq, fgetmaxdisk};
	static SOBJECT  ssobj = {ssetsrc, nopart};
	static SOBJECT  sphsobj = {sphsetsrc, nopart};
	static SOBJECT  cylsobj = {cylsetsrc, cylpart};
	static SOBJECT  rsobj = {rsetsrc, flatpart, rgetplaneq, rgetmaxdisk};

	sfun[MAT_MIRROR].mf = &mirror_vs;
	sfun[MAT_DIRECT1].mf = &direct1_vs;
	sfun[MAT_DIRECT2].mf = &direct2_vs;
	sfun[OBJ_FACE].of = &fsobj;
	sfun[OBJ_SOURCE].of = &ssobj;
	sfun[OBJ_SPHERE].of = &sphsobj;
	sfun[OBJ_CYLINDER].of = &cylsobj;
	sfun[OBJ_RING].of = &rsobj;
}


int
newsource()			/* allocate new source in our array */
{
	if (nsources == 0)
		source = (SRCREC *)malloc(SRCINC*sizeof(SRCREC));
	else if (nsources%SRCINC == 0)
		source = (SRCREC *)realloc((void *)source,
				(unsigned)(nsources+SRCINC)*sizeof(SRCREC));
	if (source == NULL)
		return(-1);
	source[nsources].sflags = 0;
	source[nsources].nhits = 1;
	source[nsources].ntests = 2;	/* initial hit probability = 1/2 */
	return(nsources++);
}


void
setflatss(src)				/* set sampling for a flat source */
register SRCREC  *src;
{
	double  mult;
	register int  i;

	src->ss[SV][0] = src->ss[SV][1] = src->ss[SV][2] = 0.0;
	for (i = 0; i < 3; i++)
		if (src->snorm[i] < 0.6 && src->snorm[i] > -0.6)
			break;
	src->ss[SV][i] = 1.0;
	fcross(src->ss[SU], src->ss[SV], src->snorm);
	mult = .5 * sqrt( src->ss2 / DOT(src->ss[SU],src->ss[SU]) );
	for (i = 0; i < 3; i++)
		src->ss[SU][i] *= mult;
	fcross(src->ss[SV], src->snorm, src->ss[SU]);
}


void
fsetsrc(src, so)			/* set a face as a source */
register SRCREC  *src;
OBJREC  *so;
{
	register FACE  *f;
	register int  i, j;
	double  d;
	
	src->sa.success = 2*AIMREQT-1;		/* bitch on second failure */
	src->so = so;
						/* get the face */
	f = getface(so);
						/* find the center */
	for (j = 0; j < 3; j++) {
		src->sloc[j] = 0.0;
		for (i = 0; i < f->nv; i++)
			src->sloc[j] += VERTEX(f,i)[j];
		src->sloc[j] /= (double)f->nv;
	}
	if (!inface(src->sloc, f))
		objerror(so, USER, "cannot hit center");
	src->sflags |= SFLAT;
	VCOPY(src->snorm, f->norm);
	src->ss2 = f->area;
						/* find maximum radius */
	src->srad = 0.;
	for (i = 0; i < f->nv; i++) {
		d = dist2(VERTEX(f,i), src->sloc);
		if (d > src->srad)
			src->srad = d;
	}
	src->srad = sqrt(src->srad);
						/* compute size vectors */
	if (f->nv == 4)				/* parallelogram case */
		for (j = 0; j < 3; j++) {
			src->ss[SU][j] = .5*(VERTEX(f,1)[j]-VERTEX(f,0)[j]);
			src->ss[SV][j] = .5*(VERTEX(f,3)[j]-VERTEX(f,0)[j]);
		}
	else
		setflatss(src);
}


void
ssetsrc(src, so)			/* set a source as a source */
register SRCREC  *src;
register OBJREC  *so;
{
	double  theta;
	
	src->sa.success = 2*AIMREQT-1;		/* bitch on second failure */
	src->so = so;
	if (so->oargs.nfargs != 4)
		objerror(so, USER, "bad arguments");
	src->sflags |= SDISTANT;
	VCOPY(src->sloc, so->oargs.farg);
	if (normalize(src->sloc) == 0.0)
		objerror(so, USER, "zero direction");
	theta = PI/180.0/2.0 * so->oargs.farg[3];
	if (theta <= FTINY)
		objerror(so, USER, "zero size");
	src->ss2 = 2.0*PI * (1.0 - cos(theta));
					/* the following is approximate */
	src->srad = sqrt(src->ss2/PI);
	VCOPY(src->snorm, src->sloc);
	setflatss(src);			/* hey, whatever works */
	src->ss[SW][0] = src->ss[SW][1] = src->ss[SW][2] = 0.0;
}


void
sphsetsrc(src, so)			/* set a sphere as a source */
register SRCREC  *src;
register OBJREC  *so;
{
	register int  i;

	src->sa.success = 2*AIMREQT-1;		/* bitch on second failure */
	src->so = so;
	if (so->oargs.nfargs != 4)
		objerror(so, USER, "bad # arguments");
	if (so->oargs.farg[3] <= FTINY)
		objerror(so, USER, "illegal radius");
	VCOPY(src->sloc, so->oargs.farg);
	src->srad = so->oargs.farg[3];
	src->ss2 = PI * src->srad * src->srad;
	for (i = 0; i < 3; i++)
		src->ss[SU][i] = src->ss[SV][i] = src->ss[SW][i] = 0.0;
	for (i = 0; i < 3; i++)
		src->ss[i][i] = .7236 * so->oargs.farg[3];
}


void
rsetsrc(src, so)			/* set a ring (disk) as a source */
register SRCREC  *src;
OBJREC  *so;
{
	register CONE  *co;
	
	src->sa.success = 2*AIMREQT-1;		/* bitch on second failure */
	src->so = so;
						/* get the ring */
	co = getcone(so, 0);
	VCOPY(src->sloc, CO_P0(co));
	if (CO_R0(co) > 0.0)
		objerror(so, USER, "cannot hit center");
	src->sflags |= SFLAT;
	VCOPY(src->snorm, co->ad);
	src->srad = CO_R1(co);
	src->ss2 = PI * src->srad * src->srad;
	setflatss(src);
}


void
cylsetsrc(src, so)			/* set a cylinder as a source */
register SRCREC  *src;
OBJREC  *so;
{
	register CONE  *co;
	register int  i;
	
	src->sa.success = 4*AIMREQT-1;		/* bitch on fourth failure */
	src->so = so;
						/* get the cylinder */
	co = getcone(so, 0);
	if (CO_R0(co) > .2*co->al)		/* heuristic constraint */
		objerror(so, WARNING, "source aspect too small");
	src->sflags |= SCYL;
	for (i = 0; i < 3; i++)
		src->sloc[i] = .5 * (CO_P1(co)[i] + CO_P0(co)[i]);
	src->srad = .5*co->al;
	src->ss2 = 2.*CO_R0(co)*co->al;
						/* set sampling vectors */
	for (i = 0; i < 3; i++)
		src->ss[SU][i] = .5 * co->al * co->ad[i];
	src->ss[SV][0] = src->ss[SV][1] = src->ss[SV][2] = 0.0;
	for (i = 0; i < 3; i++)
		if (co->ad[i] < 0.6 && co->ad[i] > -0.6)
			break;
	src->ss[SV][i] = 1.0;
	fcross(src->ss[SW], src->ss[SV], co->ad);
	normalize(src->ss[SW]);
	for (i = 0; i < 3; i++)
		src->ss[SW][i] *= .8559 * CO_R0(co);
	fcross(src->ss[SV], src->ss[SW], co->ad);
}


SPOT *
makespot(m)			/* make a spotlight */
register OBJREC  *m;
{
	register SPOT  *ns;

	if ((ns = (SPOT *)m->os) != NULL)
		return(ns);
	if ((ns = (SPOT *)malloc(sizeof(SPOT))) == NULL)
		return(NULL);
	ns->siz = 2.0*PI * (1.0 - cos(PI/180.0/2.0 * m->oargs.farg[3]));
	VCOPY(ns->aim, m->oargs.farg+4);
	if ((ns->flen = normalize(ns->aim)) == 0.0)
		objerror(m, USER, "zero focus vector");
	m->os = (char *)ns;
	return(ns);
}


int
spotout(r, s)			/* check if we're outside spot region */
register RAY  *r;
register SPOT  *s;
{
	double  d;
	FVECT  vd;
	
	if (s == NULL)
		return(0);
	if (s->flen < -FTINY) {		/* distant source */
		vd[0] = s->aim[0] - r->rorg[0];
		vd[1] = s->aim[1] - r->rorg[1];
		vd[2] = s->aim[2] - r->rorg[2];
		d = DOT(r->rdir,vd);
		/*			wrong side?
		if (d <= FTINY)
			return(1);	*/
		d = DOT(vd,vd) - d*d;
		if (PI*d > s->siz)
			return(1);	/* out */
		return(0);	/* OK */
	}
					/* local source */
	if (s->siz < 2.0*PI * (1.0 + DOT(s->aim,r->rdir)))
		return(1);	/* out */
	return(0);	/* OK */
}


double
fgetmaxdisk(ocent, op)		/* get center and squared radius of face */
FVECT  ocent;
OBJREC  *op;
{
	double  maxrad2;
	double  d;
	register int  i, j;
	register FACE  *f;
	
	f = getface(op);
	if (f->area == 0.)
		return(0.);
	for (i = 0; i < 3; i++) {
		ocent[i] = 0.;
		for (j = 0; j < f->nv; j++)
			ocent[i] += VERTEX(f,j)[i];
		ocent[i] /= (double)f->nv;
	}
	d = DOT(ocent,f->norm);
	for (i = 0; i < 3; i++)
		ocent[i] += (f->offset - d)*f->norm[i];
	maxrad2 = 0.;
	for (j = 0; j < f->nv; j++) {
		d = dist2(VERTEX(f,j), ocent);
		if (d > maxrad2)
			maxrad2 = d;
	}
	return(maxrad2);
}


double
rgetmaxdisk(ocent, op)		/* get center and squared radius of ring */
FVECT  ocent;
OBJREC  *op;
{
	register CONE  *co;
	
	co = getcone(op, 0);
	VCOPY(ocent, CO_P0(co));
	return(CO_R1(co)*CO_R1(co));
}


double
fgetplaneq(nvec, op)			/* get plane equation for face */
FVECT  nvec;
OBJREC  *op;
{
	register FACE  *fo;

	fo = getface(op);
	VCOPY(nvec, fo->norm);
	return(fo->offset);
}


double
rgetplaneq(nvec, op)			/* get plane equation for ring */
FVECT  nvec;
OBJREC  *op;
{
	register CONE  *co;

	co = getcone(op, 0);
	VCOPY(nvec, co->ad);
	return(DOT(nvec, CO_P0(co)));
}


int
commonspot(sp1, sp2, org)	/* set sp1 to intersection of sp1 and sp2 */
register SPOT  *sp1, *sp2;
FVECT  org;
{
	FVECT  cent;
	double  rad2, cos1, cos2;

	cos1 = 1. - sp1->siz/(2.*PI);
	cos2 = 1. - sp2->siz/(2.*PI);
	if (sp2->siz >= 2.*PI-FTINY)		/* BIG, just check overlap */
		return(DOT(sp1->aim,sp2->aim) >= cos1*cos2 -
					sqrt((1.-cos1*cos1)*(1.-cos2*cos2)));
				/* compute and check disks */
	rad2 = intercircle(cent, sp1->aim, sp2->aim,
			1./(cos1*cos1) - 1.,  1./(cos2*cos2) - 1.);
	if (rad2 <= FTINY || normalize(cent) == 0.)
		return(0);
	VCOPY(sp1->aim, cent);
	sp1->siz = 2.*PI*(1. - 1./sqrt(1.+rad2));
	return(1);
}


int
commonbeam(sp1, sp2, dir)	/* set sp1 to intersection of sp1 and sp2 */
register SPOT  *sp1, *sp2;
FVECT  dir;
{
	FVECT  cent, c1, c2;
	double  rad2, d;
	register int  i;
					/* move centers to common plane */
	d = DOT(sp1->aim, dir);
	for (i = 0; i < 3; i++)
		c1[i] = sp1->aim[i] - d*dir[i];
	d = DOT(sp2->aim, dir);
	for (i = 0; i < 3; i++)
		c2[i] = sp2->aim[i] - d*dir[i];
					/* compute overlap */
	rad2 = intercircle(cent, c1, c2, sp1->siz/PI, sp2->siz/PI);
	if (rad2 <= FTINY)
		return(0);
	VCOPY(sp1->aim, cent);
	sp1->siz = PI*rad2;
	return(1);
}


int
checkspot(sp, nrm)		/* check spotlight for behind source */
register SPOT  *sp;	/* spotlight */
FVECT  nrm;		/* source surface normal */
{
	double  d, d1;

	d = DOT(sp->aim, nrm);
	if (d > FTINY)			/* center in front? */
		return(1);
					/* else check horizon */
	d1 = 1. - sp->siz/(2.*PI);
	return(1.-FTINY-d*d < d1*d1);
}


double
spotdisk(oc, op, sp, pos)	/* intersect spot with object op */
FVECT  oc;
OBJREC  *op;
register SPOT  *sp;
FVECT  pos;
{
	FVECT  onorm;
	double  offs, d, dist;
	register int  i;

	offs = getplaneq(onorm, op);
	d = -DOT(onorm, sp->aim);
	if (d >= -FTINY && d <= FTINY)
		return(0.);
	dist = (DOT(pos, onorm) - offs)/d;
	if (dist < 0.)
		return(0.);
	for (i = 0; i < 3; i++)
		oc[i] = pos[i] + dist*sp->aim[i];
	return(sp->siz*dist*dist/PI/(d*d));
}


double
beamdisk(oc, op, sp, dir)	/* intersect beam with object op */
FVECT  oc;
OBJREC  *op;
register SPOT  *sp;
FVECT  dir;
{
	FVECT  onorm;
	double  offs, d, dist;
	register int  i;

	offs = getplaneq(onorm, op);
	d = -DOT(onorm, dir);
	if (d >= -FTINY && d <= FTINY)
		return(0.);
	dist = (DOT(sp->aim, onorm) - offs)/d;
	for (i = 0; i < 3; i++)
		oc[i] = sp->aim[i] + dist*dir[i];
	return(sp->siz/PI/(d*d));
}


double
intercircle(cc, c1, c2, r1s, r2s)	/* intersect two circles */
FVECT  cc;			/* midpoint (return value) */
FVECT  c1, c2;			/* circle centers */
double  r1s, r2s;		/* radii squared */
{
	double  a2, d2, l;
	FVECT  disp;
	register int  i;

	for (i = 0; i < 3; i++)
		disp[i] = c2[i] - c1[i];
	d2 = DOT(disp,disp);
					/* circle within overlap? */
	if (r1s < r2s) {
		if (r2s >= r1s + d2) {
			VCOPY(cc, c1);
			return(r1s);
		}
	} else {
		if (r1s >= r2s + d2) {
			VCOPY(cc, c2);
			return(r2s);
		}
	}
	a2 = .25*(2.*(r1s+r2s) - d2 - (r2s-r1s)*(r2s-r1s)/d2);
					/* no overlap? */
	if (a2 <= 0.)
		return(0.);
					/* overlap, compute center */
	l = sqrt((r1s - a2)/d2);
	for (i = 0; i < 3; i++)
		cc[i] = c1[i] + l*disp[i];
	return(a2);
}
