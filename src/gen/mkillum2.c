#ifndef lint
static const char	RCSid[] = "$Id: mkillum2.c,v 2.17 2004/09/19 08:42:22 greg Exp $";
#endif
/*
 * Routines to do the actual calculation for mkillum
 */

#include <string.h>

#include  "mkillum.h"
#include  "face.h"
#include  "cone.h"
#include  "random.h"
#include  "selcall.h"


static void mkaxes(FVECT u, FVECT v, FVECT n);
static void rounddir(FVECT dv, double alt, double azi);
static void flatdir(FVECT dv, double alt, double azi);


static void
rayclean(			/* finish all pending rays */
	struct rtproc *rt0
)
{
	rayflush(rt0, 1);
	while (raywait(rt0) != NULL)
		;
}


int /* XXX type conflict with otypes.h */
o_default(	/* default illum action */
	OBJREC  *ob,
	struct illum_args  *il,
	struct rtproc  *rt0,
	char  *nm
)
{
	sprintf(errmsg, "(%s): cannot make illum for %s \"%s\"",
			nm, ofun[ob->otype].funame, ob->oname);
	error(WARNING, errmsg);
	printobj(il->altmat, ob);
	return(1);
}


int
o_face(		/* make an illum face */
	OBJREC  *ob,
	struct illum_args  *il,
	struct rtproc  *rt0,
	char  *nm
)
{
#define MAXMISS		(5*n*il->nsamps)
	int  dim[3];
	int  n, nalt, nazi, h;
	float  *distarr;
	double  sp[2], r1, r2;
	FVECT  dn, org, dir;
	FVECT  u, v;
	double  ur[2], vr[2];
	int  nmisses;
	register FACE  *fa;
	register int  i, j;
				/* get/check arguments */
	fa = getface(ob);
	if (fa->area == 0.0) {
		freeface(ob);
		return(o_default(ob, il, rt0, nm));
	}
				/* set up sampling */
	if (il->sampdens <= 0)
		nalt = nazi = 1;
	else {
		n = PI * il->sampdens;
		nalt = sqrt(n/PI) + .5;
		nazi = PI*nalt + .5;
	}
	n = nalt*nazi;
	distarr = (float *)calloc(n, 3*sizeof(float));
	if (distarr == NULL)
		error(SYSTEM, "out of memory in o_face");
				/* take first edge longer than sqrt(area) */
	for (j = fa->nv-1, i = 0; i < fa->nv; j = i++) {
		u[0] = VERTEX(fa,i)[0] - VERTEX(fa,j)[0];
		u[1] = VERTEX(fa,i)[1] - VERTEX(fa,j)[1];
		u[2] = VERTEX(fa,i)[2] - VERTEX(fa,j)[2];
		if ((r1 = DOT(u,u)) >= fa->area-FTINY)
			break;
	}
	if (i < fa->nv) {	/* got one! -- let's align our axes */
		r2 = 1.0/sqrt(r1);
		u[0] *= r2; u[1] *= r2; u[2] *= r2;
		fcross(v, fa->norm, u);
	} else			/* oh well, we'll just have to wing it */
		mkaxes(u, v, fa->norm);
				/* now, find limits in (u,v) coordinates */
	ur[0] = vr[0] = FHUGE;
	ur[1] = vr[1] = -FHUGE;
	for (i = 0; i < fa->nv; i++) {
		r1 = DOT(VERTEX(fa,i),u);
		if (r1 < ur[0]) ur[0] = r1;
		if (r1 > ur[1]) ur[1] = r1;
		r2 = DOT(VERTEX(fa,i),v);
		if (r2 < vr[0]) vr[0] = r2;
		if (r2 > vr[1]) vr[1] = r2;
	}
	dim[0] = random();
				/* sample polygon */
	nmisses = 0;
	for (dim[1] = 0; dim[1] < nalt; dim[1]++)
	    for (dim[2] = 0; dim[2] < nazi; dim[2]++)
		for (i = 0; i < il->nsamps; i++) {
					/* random direction */
		    h = ilhash(dim, 3) + i;
		    multisamp(sp, 2, urand(h));
		    r1 = (dim[1] + sp[0])/nalt;
		    r2 = (dim[2] + sp[1] - .5)/nazi;
		    flatdir(dn, r1, r2);
		    for (j = 0; j < 3; j++)
			dir[j] = -dn[0]*u[j] - dn[1]*v[j] - dn[2]*fa->norm[j];
					/* random location */
		    do {
			multisamp(sp, 2, urand(h+4862+nmisses));
			r1 = ur[0] + (ur[1]-ur[0]) * sp[0];
			r2 = vr[0] + (vr[1]-vr[0]) * sp[1];
			for (j = 0; j < 3; j++)
			    org[j] = r1*u[j] + r2*v[j]
					+ fa->offset*fa->norm[j];
		    } while (!inface(org, fa) && nmisses++ < MAXMISS);
		    if (nmisses > MAXMISS) {
			objerror(ob, WARNING, "bad aspect");
			rayclean(rt0);
			freeface(ob);
			free((void *)distarr);
			return(o_default(ob, il, rt0, nm));
		    }
		    for (j = 0; j < 3; j++)
			org[j] += .001*fa->norm[j];
					/* send sample */
		    raysamp(distarr+3*(dim[1]*nazi+dim[2]), org, dir, rt0);
		}
	rayclean(rt0);
				/* write out the face and its distribution */
	if (average(il, distarr, nalt*nazi)) {
		if (il->sampdens > 0)
			flatout(il, distarr, nalt, nazi, u, v, fa->norm);
		illumout(il, ob);
	} else
		printobj(il->altmat, ob);
				/* clean up */
	freeface(ob);
	free((void *)distarr);
	return(0);
#undef MAXMISS
}


int
o_sphere(	/* make an illum sphere */
	register OBJREC  *ob,
	struct illum_args  *il,
	struct rtproc  *rt0,
	char  *nm
)
{
	int  dim[3];
	int  n, nalt, nazi;
	float  *distarr;
	double  sp[4], r1, r2, r3;
	FVECT  org, dir;
	FVECT  u, v;
	register int  i, j;
				/* check arguments */
	if (ob->oargs.nfargs != 4)
		objerror(ob, USER, "bad # of arguments");
				/* set up sampling */
	if (il->sampdens <= 0)
		nalt = nazi = 1;
	else {
		n = 4.*PI * il->sampdens;
		nalt = sqrt(2./PI*n) + .5;
		nazi = PI/2.*nalt + .5;
	}
	n = nalt*nazi;
	distarr = (float *)calloc(n, 3*sizeof(float));
	if (distarr == NULL)
		error(SYSTEM, "out of memory in o_sphere");
	dim[0] = random();
				/* sample sphere */
	for (dim[1] = 0; dim[1] < nalt; dim[1]++)
	    for (dim[2] = 0; dim[2] < nazi; dim[2]++)
		for (i = 0; i < il->nsamps; i++) {
					/* next sample point */
		    multisamp(sp, 4, urand(ilhash(dim,3)+i));
					/* random direction */
		    r1 = (dim[1] + sp[0])/nalt;
		    r2 = (dim[2] + sp[1] - .5)/nazi;
		    rounddir(dir, r1, r2);
					/* random location */
		    mkaxes(u, v, dir);		/* yuck! */
		    r3 = sqrt(sp[2]);
		    r2 = 2.*PI*sp[3];
		    r1 = r3*ob->oargs.farg[3]*cos(r2);
		    r2 = r3*ob->oargs.farg[3]*sin(r2);
		    r3 = ob->oargs.farg[3]*sqrt(1.01-r3*r3);
		    for (j = 0; j < 3; j++) {
			org[j] = ob->oargs.farg[j] + r1*u[j] + r2*v[j] +
					r3*dir[j];
			dir[j] = -dir[j];
		    }
					/* send sample */
		    raysamp(distarr+3*(dim[1]*nazi+dim[2]), org, dir, rt0);
		}
	rayclean(rt0);
				/* write out the sphere and its distribution */
	if (average(il, distarr, nalt*nazi)) {
		if (il->sampdens > 0)
			roundout(il, distarr, nalt, nazi);
		else
			objerror(ob, WARNING, "diffuse distribution");
		illumout(il, ob);
	} else
		printobj(il->altmat, ob);
				/* clean up */
	free((void *)distarr);
	return(1);
}


int
o_ring(		/* make an illum ring */
	OBJREC  *ob,
	struct illum_args  *il,
	struct rtproc  *rt0,
	char  *nm
)
{
	int  dim[3];
	int  n, nalt, nazi;
	float  *distarr;
	double  sp[4], r1, r2, r3;
	FVECT  dn, org, dir;
	FVECT  u, v;
	register CONE  *co;
	register int  i, j;
				/* get/check arguments */
	co = getcone(ob, 0);
				/* set up sampling */
	if (il->sampdens <= 0)
		nalt = nazi = 1;
	else {
		n = PI * il->sampdens;
		nalt = sqrt(n/PI) + .5;
		nazi = PI*nalt + .5;
	}
	n = nalt*nazi;
	distarr = (float *)calloc(n, 3*sizeof(float));
	if (distarr == NULL)
		error(SYSTEM, "out of memory in o_ring");
	mkaxes(u, v, co->ad);
	dim[0] = random();
				/* sample disk */
	for (dim[1] = 0; dim[1] < nalt; dim[1]++)
	    for (dim[2] = 0; dim[2] < nazi; dim[2]++)
		for (i = 0; i < il->nsamps; i++) {
					/* next sample point */
		    multisamp(sp, 4, urand(ilhash(dim,3)+i));
					/* random direction */
		    r1 = (dim[1] + sp[0])/nalt;
		    r2 = (dim[2] + sp[1] - .5)/nazi;
		    flatdir(dn, r1, r2);
		    for (j = 0; j < 3; j++)
			dir[j] = -dn[0]*u[j] - dn[1]*v[j] - dn[2]*co->ad[j];
					/* random location */
		    r3 = sqrt(CO_R0(co)*CO_R0(co) +
			    sp[2]*(CO_R1(co)*CO_R1(co) - CO_R0(co)*CO_R0(co)));
		    r2 = 2.*PI*sp[3];
		    r1 = r3*cos(r2);
		    r2 = r3*sin(r2);
		    for (j = 0; j < 3; j++)
			org[j] = CO_P0(co)[j] + r1*u[j] + r2*v[j] +
					.001*co->ad[j];

					/* send sample */
		    raysamp(distarr+3*(dim[1]*nazi+dim[2]), org, dir, rt0);
		}
	rayclean(rt0);
				/* write out the ring and its distribution */
	if (average(il, distarr, nalt*nazi)) {
		if (il->sampdens > 0)
			flatout(il, distarr, nalt, nazi, u, v, co->ad);
		illumout(il, ob);
	} else
		printobj(il->altmat, ob);
				/* clean up */
	freecone(ob);
	free((void *)distarr);
	return(1);
}


void
raysamp(	/* queue a ray sample */
	float  res[3],
	FVECT  org,
	FVECT  dir,
	struct rtproc *rt0
)
{
	register struct rtproc  *rt;
	register float  *fp;

	for (rt = rt0; rt != NULL; rt = rt->next)
		if (rt->nrays < rt->bsiz && rt->dest[rt->nrays] == NULL)
			break;
	if (rt == NULL)		/* need to free up buffer? */
		rt = raywait(rt0);
	if (rt == NULL)
		error(SYSTEM, "raywait() returned NULL in raysamp()");
	fp = rt->buf + 6*rt->nrays;
	*fp++ = org[0]; *fp++ = org[1]; *fp++ = org[2];
	*fp++ = dir[0]; *fp++ = dir[1]; *fp = dir[2];
	rt->dest[rt->nrays++] = res;
	if (rt->nrays == rt->bsiz)
		rayflush(rt, 0);
}


void
rayflush(			/* flush queued rays to rtrace */
	register struct rtproc  *rt,
	int  doall
)
{
	int	nw;

	do {
		if (rt->nrays <= 0)
			continue;
		memset(rt->buf+6*rt->nrays, 0, 6*sizeof(float));
		nw = 6*sizeof(float)*(rt->nrays+1);
		errno = 0;
		if (writebuf(rt->pd.w, (char *)rt->buf, nw) < nw)
			error(SYSTEM, "error writing to rtrace process");
		rt->nrays = 0;		/* flag buffer as flushed */
	} while (doall && (rt = rt->next) != NULL);
}


struct rtproc *
raywait(			/* retrieve rtrace results */
	struct rtproc *rt0
)
{
	fd_set	readset, errset;
	int	nr;
	struct rtproc	*rtfree;
	register struct rtproc	*rt;
	register int	n;
				/* prepare select call */
	FD_ZERO(&readset); FD_ZERO(&errset); n = 0;
	nr = 0;
	for (rt = rt0; rt != NULL; rt = rt->next) {
		if (rt->nrays == 0 && rt->dest[0] != NULL) {
			FD_SET(rt->pd.r, &readset);
			++nr;
		}
		FD_SET(rt->pd.r, &errset);
		if (rt->pd.r >= n)
			n = rt->pd.r + 1;
	}
	if (!nr)		/* no rays pending */
		return(NULL);
	if (nr > 1) {		/* call select for multiple processes */
		errno = 0;
		if (select(n, &readset, NULL, &errset, NULL) < 0)
			error(SYSTEM, "select call error in raywait()");
	} else
		FD_ZERO(&errset);
	rtfree = NULL;		/* read from ready process(es) */
	for (rt = rt0; rt != NULL; rt = rt->next) {
		if (!FD_ISSET(rt->pd.r, &readset) &&
				!FD_ISSET(rt->pd.r, &errset))
			continue;
		for (n = 0; n < rt->bsiz && rt->dest[n] != NULL; n++)
			;
		nr = 3*sizeof(float)*(n+1);
		if (readbuf(rt->pd.r, (char *)rt->buf, nr) < nr)
			error(USER, "rtrace process died");
		while (n-- > 0) {
			rt->dest[n][0] += rt->buf[3*n];
			rt->dest[n][1] += rt->buf[3*n+1];
			rt->dest[n][2] += rt->buf[3*n+2];
			rt->dest[n] = NULL;
		}
		rtfree = rt;
	}
	return(rtfree);
}


static void
mkaxes(			/* compute u and v to go with n */
	FVECT  u,
	FVECT  v,
	FVECT  n
)
{
	register int  i;

	v[0] = v[1] = v[2] = 0.0;
	for (i = 0; i < 3; i++)
		if (n[i] < 0.6 && n[i] > -0.6)
			break;
	v[i] = 1.0;
	fcross(u, v, n);
	normalize(u);
	fcross(v, n, u);
}


static void
rounddir(		/* compute uniform spherical direction */
	register FVECT  dv,
	double  alt,
	double  azi
)
{
	double  d1, d2;

	dv[2] = 1. - 2.*alt;
	d1 = sqrt(1. - dv[2]*dv[2]);
	d2 = 2.*PI * azi;
	dv[0] = d1*cos(d2);
	dv[1] = d1*sin(d2);
}


static void
flatdir(		/* compute uniform hemispherical direction */
	register FVECT  dv,
	double  alt,
	double  azi
)
{
	double  d1, d2;

	d1 = sqrt(alt);
	d2 = 2.*PI * azi;
	dv[0] = d1*cos(d2);
	dv[1] = d1*sin(d2);
	dv[2] = sqrt(1. - alt);
}
