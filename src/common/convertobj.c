#ifndef lint
static const char RCSid[] = "$Id: convertobj.c,v 2.6 2021/04/17 14:51:33 greg Exp $";
#endif
/*
 *  convertobj.c
 *
 *  Convert .OBJ scene to Radiance
 *
 *  Created by Greg Ward on Thu Apr 22 2004.
 */

#include "paths.h"
#include "fvect.h"
#include "tmesh.h"
#include "rterror.h"
#include "objutil.h"

#define TEXNAME		"T-nor"

static int      fcnt = 0;			/* face output counter */

/* Callback for face smoothing detection */
static int
checksmooth(Scene *sc, Face *f, void *ptr)
{
	int	nrev = 0;
	Normal	fnrm;
	int	i;

	f->flags &= ~FACE_RESERVED;		/* using reserved flag */

	for (i = f->nv; i--; )
		if (f->v[i].nid < 0)
			return(0);		/* missing normal */

	if (faceArea(sc, f, fnrm) == 0)		/* degenerate?? */
		return(0);
						/* check each normal */
	for (i = f->nv; i--; ) {
		float	*tnrm = sc->norm[f->v[i].nid];
		double	dprod = DOT(tnrm, fnrm);
		if (dprod >= COSTOL)		/* this one agrees? */
			continue;
		if (dprod > -COSTOL)		/* active smoothing? */
			break;
		++nrev;				/* count opposite face normal */
	}
	if ((i < 0) & !nrev)			/* all normals agree w/ face? */
		return(0);
	if (nrev == f->nv) {			/* all reversed? */
		for (i = f->nv/2; i--; ) {	/* swap vertices around */
			int	j = f->nv-1 - i;
			VertEnt	tve = f->v[i];
			f->v[i] = f->v[j];
			f->v[j] = tve;
		}
		return(0);
	}
	f->flags |= FACE_RESERVED;		/* else we got one to smooth */
	return(1);
}

/* Callback to write out smoothed Radiance triangle */
static int
trismooth(Scene *sc, Face *f, void *ptr)
{
	FILE		*fp = (FILE *)ptr;
	BARYCCM		bcm;
	FVECT		coor[3];
	int		i;

	if (f->nv != 3)
		return(0);			/* should never happen */
#ifdef  SMLFLT
	for (i = 3; i--; ) {
		double	*v = sc->vert[f->v[i].vid].p;
		VCOPY(coor[i], v);
	}
	if (comp_baryc(&bcm, coor[0], coor[1], coor[2]) < 0)
		return(0);			/* degenerate?? */
#else
	if (comp_baryc(&bcm, sc->vert[f->v[0].vid].p,
			sc->vert[f->v[1].vid].p, sc->vert[f->v[2].vid].p) < 0)
		return(0);			/* degenerate?? */
#endif
	for (i = 3; i--; ) {			/* assign BC normals */
		float	*tnrm = sc->norm[f->v[i].nid];
		coor[0][i] = tnrm[0];
		coor[1][i] = tnrm[1];
		coor[2][i] = tnrm[2];
	}					/* print texture */
	fprintf(fp, "\n%s texfunc %s\n4 dx dy dz %s\n0\n",
			sc->matname[f->mat], TEXNAME, TCALNAME);
	fput_baryc(&bcm, coor, 3, fp);		/* with BC normals */
	fprintf(fp, "\n%s polygon %s.%d\n0\n0\n9\n",
			TEXNAME, sc->grpname[f->grp], ++fcnt);
	for (i = 0; i < 3; i++) {		/* then triangle */
		double	*v = sc->vert[f->v[i].vid].p;
		fprintf(fp, "\t%18.12g %18.12g %18.12g\n", v[0], v[1], v[2]);
	}
	return(1);
}

/* Callback to convert face to Radiance */
static int
radface(Scene *sc, Face *f, void *ptr)
{
	FILE		*fp = (FILE *)ptr;
	int		i;
	
	fprintf(fp, "\n%s polygon %s.%d\n0\n0\n%d\n", sc->matname[f->mat],
				sc->grpname[f->grp], ++fcnt, 3*f->nv);
	for (i = 0; i < f->nv; i++) {
		Vertex  *vp = sc->vert + f->v[i].vid;
		fprintf(fp, "\t%18.12g %18.12g %18.12g\n",
				vp->p[0], vp->p[1], vp->p[2]);
	}
	return(1);
}

/* Convert indicated faces to Radiance, return # written or -1 on error */
int
toRadiance(Scene *sc, FILE *fp, int flreq, int flexc)
{
	int     n;

	if (sc == NULL || sc->nfaces <= 0 || fp == NULL)
		return(0);
						/* not passing empties */
	flexc |= FACE_DEGENERATE;
						/* reset counter if not stdout */
	fcnt *= (fp == stdout);
						/* write comments */
	for (n = 0; n < sc->ndescr; n++)
		fprintf(fp, "# %s\n", sc->descr[n]);
						/* flag faces to smooth */
	n = foreachFace(sc, checksmooth, flreq, flexc, NULL);
	if (n > 0) {				/* write smooth faces */
		Scene	*smoothies = dupScene(sc, flreq|FACE_RESERVED, flexc);
		if (!smoothies)
			return(-1);
		n = triangulateScene(smoothies);
		if (n >= 0)
			n = foreachFace(smoothies, trismooth, 0, 0, fp);
		freeScene(smoothies);
	}
	if (n < 0)
		return(-1);
						/* write flat faces */
	n += foreachFace(sc, radface, flreq, flexc|FACE_RESERVED, fp);

	if (fflush(fp) < 0) {
		error(SYSTEM, "Error writing Radiance scene data");
		return(-1);
	}
	return(n);
}

/* Convert faces to Radiance file, return # written or -1 on error */
int
writeRadiance(Scene *sc, const char *fspec, int flreq, int flexc)
{
	extern char     *progname;
	FILE		*fp;
	int		n;

	if (sc == NULL || sc->nfaces <= 0 || fspec == NULL || !*fspec)
		return(0);
#if POPEN_SUPPORT
	if (fspec[0] == '!') {
		if ((fp = popen(fspec+1, "w")) == NULL) {
			sprintf(errmsg, "%s: cannot execute", fspec);
			error(SYSTEM, errmsg);
			return(-1);
		}
	} else
#endif
	if ((fp = fopen(fspec, "w")) == NULL) {
		sprintf(errmsg, "%s: cannot open for writing", fspec);
		error(SYSTEM, errmsg);
		return(-1);
	}
						/* start off header */
	fprintf(fp, "# Radiance scene file converted from .OBJ by %s\n#\n",
						progname);
	n = toRadiance(sc, fp, flreq, flexc);   /* write file */
#if POPEN_SUPPORT
	if (fspec[0] == '!') {
		if (pclose(fp)) {
			sprintf(errmsg, "%s: error writing to command\n", fspec);
			error(SYSTEM, errmsg);
			return(-1);
		}
	} else
#endif
	fclose(fp);				/* close it (already flushed) */
	return(n);
}
