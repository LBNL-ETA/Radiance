#ifndef lint
static const char RCSid[] = "$Id: convertobj.c,v 2.3 2020/04/23 22:35:27 greg Exp $";
#endif
/*
 *  convertobj.c
 *
 *  Convert .OBJ scene to Radiance
 *
 *  Created by Greg Ward on Thu Apr 22 2004.
 */

#include "paths.h"
#include "rterror.h"
#include "objutil.h"

/* Callback to convert face to Radiance */
static int
radface(Scene *sc, Face *f, void *ptr)
{
	static int      fcnt = 0;
	FILE		*fp = (FILE *)ptr;
	int		i;
	
	if (f->flags & FACE_DEGENERATE)
		return(0);
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
						/* write comments */
	for (n = 0; n < sc->ndescr; n++)
		fprintf(fp, "# %s\n", sc->descr[n]);
						/* write faces */
	n = foreachFace(sc, radface, flreq, flexc, (void *)fp);
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
