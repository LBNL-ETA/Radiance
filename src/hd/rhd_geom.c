#ifndef lint
static const char	RCSid[] = "$Id: rhd_geom.c,v 3.14 2004/01/01 11:21:55 schorsch Exp $";
#endif
/*
 * Geometry drawing operations for OpenGL driver.
 */

#include <string.h>

#include "rtio.h"
#include "radogl.h"
#include "rhdriver.h"
#include "rhd_odraw.h"

#ifndef MAXGEO
#define	MAXGEO		8		/* maximum geometry list length */
#endif
#ifndef MAXPORT
#define MAXPORT		(MAXGEO*4)	/* maximum number of portal files */
#endif

int	gmPortals = 0;			/* current portal GL list id */
static int	Nlists = 0;		/* number of lists allocated */
static char	*curportlist[MAXPORT];	/* current portal list */
static char	*newportlist[MAXPORT];	/* new portal file list */

static struct gmEntry {
	char	*gfile;			/* geometry file name */
	FVECT	cent;			/* centroid */
	RREAL	rad;			/* radius */
	int	listid;			/* display list identifier */
	int	nlists;			/* number of lists allocated */
} gmCurrent[MAXGEO], gmNext[MAXGEO];	/* current and next list */

#define FORALLGEOM(ot,i)	for (i=0;i<MAXGEO&&ot[i].gfile!=NULL;i++)

#define FORALLPORT(pl,i)		for (i=0;i<MAXPORT&&pl[i]!=NULL;i++)




extern void
gmNewGeom(			/* add new geometry to next list */
	char	*file
)
{
	register int	i, j;
					/* check if already in next list */
	FORALLGEOM(gmNext, i)
		if (!strcmp(file, gmNext[i].gfile))
			return;
	if (i >= MAXGEO) {
		error(WARNING, "too many section octrees -- ignoring extra");
		return;
	}
					/* check if copy in current list */
	FORALLGEOM(gmCurrent, j)
		if (!strcmp(file, gmCurrent[j].gfile)) {
			gmNext[i] = gmCurrent[j];
			return;
		}
					/* else load new octree */
	gmNext[i].gfile = savestr(file);
	dolights = 0;
	domats = 1;
	gmNext[i].listid = rgl_octlist(file, gmNext[i].cent, &gmNext[i].rad,
					&gmNext[i].nlists);
	gmNext[i].rad *= 1.732;		/* go to corners */
#ifdef DEBUG
	fprintf(stderr, "Loaded octree \"%s\" into listID %d with radius %f\n",
			file, gmNext[i].listid, gmNext[i].rad);
#endif
}


extern void
gmEndGeom(void)			/* make next list current */
{
	register int	i, j;

	FORALLGEOM(gmCurrent, i) {
		FORALLGEOM(gmNext, j)
			if (gmNext[j].listid == gmCurrent[i].listid)
				break;
		if (j >= MAXGEO || gmNext[j].gfile == NULL) {
			glDeleteLists(gmCurrent[i].listid,	/* not found */
					gmCurrent[i].nlists);
			freestr(gmCurrent[i].gfile);
		}
	}
	memcpy((void *)gmCurrent, (void *)gmNext, sizeof(gmNext));
	memset((void *)gmNext, '\0', sizeof(gmNext));
}


extern int
gmDrawGeom(void)			/* draw current list of octrees */
{
	register int	n;

	FORALLGEOM(gmCurrent, n)
		glCallList(gmCurrent[n].listid);
	return(n);
}


extern void
gmDrawPortals(	/* draw portals with specific RGBA value */
	int	r,
	int	g,
	int	b,
	int	a
)
{
	if (!gmPortals || (r<0) & (g<0) & (b<0) & (a<0))
		return;
	glPushAttrib(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT|
			GL_POLYGON_BIT|GL_LIGHTING_BIT);
	glDisable(GL_LIGHTING);
	glDisable(GL_DITHER);
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
					/* don't actually write depth */
	glDepthMask(GL_FALSE);
					/* draw only selected channels */
	glColorMask(r>=0, g>=0, b>=0, a>=0);
	glColor4ub(r&0xff, g&0xff, b&0xff, a&0xff);
	glCallList(gmPortals);		/* draw them portals */
	glPopAttrib();
}


extern void
gmDepthLimit(	/* compute approximate depth limits for view */
	double	dl[2],
	FVECT	vorg,
	FVECT	vdir
)
{
	FVECT	v;
	double	dcent;
	register int	i;

	dl[0] = FHUGE; dl[1] = 0.;
	FORALLGEOM(gmCurrent, i) {
		VSUB(v, gmCurrent[i].cent, vorg);
		dcent = DOT(v, vdir);
		if (dl[0] > dcent-gmCurrent[i].rad)
			dl[0] = dcent-gmCurrent[i].rad;
		if (dl[1] < dcent+gmCurrent[i].rad)
			dl[1] = dcent+gmCurrent[i].rad;
	}
	if (dl[0] < 0.)
		dl[0] = 0.;
}


extern void
gmNewPortal(		/* add portal file(s) to our new list */
	char	*pflist
)
{
	register int	i;
	char	newfile[128];

	if (pflist == NULL)
		return;
	while ((pflist = nextword(newfile, sizeof(newfile), pflist)) != NULL) {
		FORALLPORT(newportlist,i)
			if (!strcmp(newportlist[i], newfile))
				goto endloop;	/* in list already */
		if (i >= MAXPORT) {
			error(WARNING, "too many portals -- ignoring extra");
			return;
		}
		newportlist[i] = savestr(newfile);
	endloop:;
	}
}


static int
sstrcmp(ss0, ss1)
char	**ss0, **ss1;
{
	return(strcmp(*ss0, *ss1));
}


extern int
gmEndPortal(void)			/* close portal list and return GL list */
{
	register int	n;

	FORALLPORT(newportlist, n);
	if (!n) {			/* free old GL list */
		if (gmPortals)
			glDeleteLists(gmPortals, Nlists);
		gmPortals = 0;
	} else
		qsort(newportlist, n, sizeof(char *), sstrcmp);
	FORALLPORT(newportlist, n)		/* compare sorted lists */
		if (curportlist[n] == NULL ||
				strcmp(curportlist[n],newportlist[n])) {
						/* load new list */
			if (gmPortals)
				glDeleteLists(gmPortals, 1);
			FORALLPORT(newportlist, n);
			dolights = 0;
			domats = 0;
			gmPortals = rgl_filelist(n, newportlist, &Nlists);
			break;
		}
	FORALLPORT(curportlist, n)		/* free old file list */
		freestr(curportlist[n]);
	memcpy((void *)curportlist, (void *)newportlist, sizeof(newportlist));
	memset((void *)newportlist, '\0', sizeof(newportlist));
	return(gmPortals);			/* return GL list id */
}
