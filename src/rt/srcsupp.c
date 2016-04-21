#ifndef lint
static const char	RCSid[] = "$Id: srcsupp.c,v 2.23 2016/04/21 00:40:35 greg Exp $";
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

#define SRCINC		32		/* realloc increment for array */

SRCREC  *source = NULL;			/* our list of sources */
int  nsources = 0;			/* the number of sources */

SRCFUNC  sfun[NUMOTYPE];		/* source dispatch table */


void
initstypes(void)			/* initialize source dispatch table */
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
newsource(void)			/* allocate new source in our array */
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
	source[nsources].ntests = 2;	/* initial hit probability = 50% */
#if SHADCACHE
	source[nsources].obscache = NULL;
#endif
	return(nsources++);
}


void
setflatss(				/* set sampling for a flat source */
	SRCREC  *src
)
{
	double  mult;
	int  i;

	getperpendicular(src->ss[SU], src->snorm, rand_samp);
	mult = .5 * sqrt( src->ss2 );
	for (i = 0; i < 3; i++)
		src->ss[SU][i] *= mult;
	fcross(src->ss[SV], src->snorm, src->ss[SU]);
}


void
fsetsrc(			/* set a face as a source */
	SRCREC  *src,
	OBJREC  *so
)
{
	FACE  *f;
	int  i, j;
	double  d;
	
	src->sa.success = 2*AIMREQT-1;		/* bitch on second failure */
	src->so = so;
						/* get the face */
	f = getface(so);
	if (f->area == 0.)
		objerror(so, USER, "zero source area");
						/* find the center */
	for (j = 0; j < 3; j++) {
		src->sloc[j] = 0.0;
		for (i = 0; i < f->nv; i++)
			src->sloc[j] += VERTEX(f,i)[j];
		src->sloc[j] /= (double)f->nv;
	}
	if (!inface(src->sloc, f))
		objerror(so, USER, "cannot hit source center");
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
ssetsrc(			/* set a source as a source */
	SRCREC  *src,
	OBJREC  *so
)
{
	double  theta;
	
	src->sa.success = 2*AIMREQT-1;		/* bitch on second failure */
	src->so = so;
	if (so->oargs.nfargs != 4)
		objerror(so, USER, "bad arguments");
	src->sflags |= (SDISTANT|SCIR);
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
sphsetsrc(			/* set a sphere as a source */
	SRCREC  *src,
	OBJREC  *so
)
{
	int  i;

	src->sa.success = 2*AIMREQT-1;		/* bitch on second failure */
	src->so = so;
	if (so->oargs.nfargs != 4)
		objerror(so, USER, "bad # arguments");
	if (so->oargs.farg[3] <= FTINY)
		objerror(so, USER, "illegal source radius");
	src->sflags |= SCIR;
	VCOPY(src->sloc, so->oargs.farg);
	src->srad = so->oargs.farg[3];
	src->ss2 = PI * src->srad * src->srad;
	for (i = 0; i < 3; i++)
		src->ss[SU][i] = src->ss[SV][i] = src->ss[SW][i] = 0.0;
	for (i = 0; i < 3; i++)
		src->ss[i][i] = 0.7236 * so->oargs.farg[3];
}


void
rsetsrc(			/* set a ring (disk) as a source */
	SRCREC  *src,
	OBJREC  *so
)
{
	CONE  *co;
	
	src->sa.success = 2*AIMREQT-1;		/* bitch on second failure */
	src->so = so;
						/* get the ring */
	co = getcone(so, 0);
	if (co == NULL)
		objerror(so, USER, "illegal source");
	if (CO_R1(co) <= FTINY)
		objerror(so, USER, "illegal source radius");
	VCOPY(src->sloc, CO_P0(co));
	if (CO_R0(co) > 0.0)
		objerror(so, USER, "cannot hit source center");
	src->sflags |= (SFLAT|SCIR);
	VCOPY(src->snorm, co->ad);
	src->srad = CO_R1(co);
	src->ss2 = PI * src->srad * src->srad;
	setflatss(src);
}


void
cylsetsrc(			/* set a cylinder as a source */
	SRCREC  *src,
	OBJREC  *so
)
{
	CONE  *co;
	int  i;
	
	src->sa.success = 4*AIMREQT-1;		/* bitch on fourth failure */
	src->so = so;
						/* get the cylinder */
	co = getcone(so, 0);
	if (co == NULL)
		objerror(so, USER, "illegal source");
	if (CO_R0(co) <= FTINY)
		objerror(so, USER, "illegal source radius");
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
	getperpendicular(src->ss[SW], co->ad, rand_samp);
	for (i = 0; i < 3; i++)
		src->ss[SW][i] *= .8559 * CO_R0(co);
	fcross(src->ss[SV], src->ss[SW], co->ad);
}


SPOT *
makespot(			/* make a spotlight */
	OBJREC  *m
)
{
	SPOT  *ns;

	if ((ns = (SPOT *)m->os) != NULL)
		return(ns);
	if ((ns = (SPOT *)malloc(sizeof(SPOT))) == NULL)
		return(NULL);
	if (m->oargs.farg[3] <= FTINY)
		objerror(m, USER, "zero angle");
	ns->siz = 2.0*PI * (1.0 - cos(PI/180.0/2.0 * m->oargs.farg[3]));
	VCOPY(ns->aim, m->oargs.farg+4);
	if ((ns->flen = normalize(ns->aim)) == 0.0)
		objerror(m, USER, "zero focus vector");
	m->os = (char *)ns;
	return(ns);
}


int
spotout(			/* check if we're outside spot region */
	RAY  *r,
	SPOT  *s
)
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
fgetmaxdisk(		/* get center and squared radius of face */
	FVECT  ocent,
	OBJREC  *op
)
{
	double  maxrad2;
	double  d;
	int  i, j;
	FACE  *f;
	
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
rgetmaxdisk(		/* get center and squared radius of ring */
	FVECT  ocent,
	OBJREC  *op
)
{
	CONE  *co;
	
	co = getcone(op, 0);
	if (co == NULL)
		return(0.);
	VCOPY(ocent, CO_P0(co));
	return(CO_R1(co)*CO_R1(co));
}


double
fgetplaneq(			/* get plane equation for face */
	FVECT  nvec,
	OBJREC  *op
)
{
	FACE  *fo;

	fo = getface(op);
	VCOPY(nvec, fo->norm);
	return(fo->offset);
}


double
rgetplaneq(			/* get plane equation for ring */
	FVECT  nvec,
	OBJREC  *op
)
{
	CONE  *co;

	co = getcone(op, 0);
	if (co == NULL) {
		memset(nvec, 0, sizeof(FVECT));
		return(0.);
	}
	VCOPY(nvec, co->ad);
	return(DOT(nvec, CO_P0(co)));
}


int
commonspot(		/* set sp1 to intersection of sp1 and sp2 */
	SPOT  *sp1,
	SPOT  *sp2,
	FVECT  org
)
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
commonbeam(		/* set sp1 to intersection of sp1 and sp2 */
	SPOT  *sp1,
	SPOT  *sp2,
	FVECT  dir
)
{
	FVECT  cent, c1, c2;
	double  rad2, d;
					/* move centers to common plane */
	d = DOT(sp1->aim, dir);
	VSUM(c1, sp1->aim, dir, -d);
	d = DOT(sp2->aim, dir);
	VSUM(c2, sp2->aim, dir, -d);
					/* compute overlap */
	rad2 = intercircle(cent, c1, c2, sp1->siz/PI, sp2->siz/PI);
	if (rad2 <= FTINY)
		return(0);
	VCOPY(sp1->aim, cent);
	sp1->siz = PI*rad2;
	return(1);
}


int
checkspot(			/* check spotlight for behind source */
	SPOT  *sp,	/* spotlight */
	FVECT  nrm	/* source surface normal */
)
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
spotdisk(		/* intersect spot with object op */
	FVECT  oc,
	OBJREC  *op,
	SPOT  *sp,
	FVECT  pos
)
{
	FVECT  onorm;
	double  offs, d, dist;

	offs = getplaneq(onorm, op);
	d = -DOT(onorm, sp->aim);
	if (d >= -FTINY && d <= FTINY)
		return(0.);
	dist = (DOT(pos, onorm) - offs)/d;
	if (dist < 0.)
		return(0.);
	VSUM(oc, pos, sp->aim, dist);
	return(sp->siz*dist*dist/PI/(d*d));
}


double
beamdisk(		/* intersect beam with object op */
	FVECT  oc,
	OBJREC  *op,
	SPOT  *sp,
	FVECT  dir
)
{
	FVECT  onorm;
	double  offs, d, dist;

	offs = getplaneq(onorm, op);
	d = -DOT(onorm, dir);
	if (d >= -FTINY && d <= FTINY)
		return(0.);
	dist = (DOT(sp->aim, onorm) - offs)/d;
	VSUM(oc, sp->aim, dir, dist);
	return(sp->siz/PI/(d*d));
}


double
intercircle(				/* intersect two circles */
	FVECT  cc,		/* midpoint (return value) */
	FVECT  c1,		/* circle centers */
	FVECT  c2,
	double  r1s,		/* radii squared */
	double	r2s
)
{
	double  a2, d2, l;
	FVECT  disp;

	VSUB(disp, c2, c1);
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
	VSUM(cc, c1, disp, l);
	return(a2);
}
