#ifndef lint
static const char	RCSid[] = "$Id: viewbeams.c,v 3.6 2003/07/21 22:30:18 schorsch Exp $";
#endif
/*
 * Convert view to beam list.
 */

#include <string.h>

#include "rholo.h"
#include "view.h"
#include "random.h"

#ifndef MINRES
#define MINRES		12		/* section sample resolution */
#endif
#ifndef NVSAMPS
#define NVSAMPS		16384		/* target number of view samples */
#endif

#define OJITTER		0.25		/* amount to jitter sample origin */

#define BALLOCBLK	128		/* beam allocation block size */


static BEAMLIST	blist;


static
add2blist(hd, bi, nr)			/* add to beam sample list */
int	hd, bi, nr;
{
	register int	i;

	for (i = blist.nb; i--; )
		if (blist.bl[i].bi == bi && blist.bl[i].hd == hd) {
			blist.bl[i].nr += nr;	/* found it */
			return;
		}
	i = blist.nb++;				/* else add beam to list */
	if (i % BALLOCBLK == 0) {
		blist.bl = (PACKHEAD *)realloc((void *)blist.bl,
				(i+BALLOCBLK)*sizeof(PACKHEAD));
		CHECK(blist.bl==NULL, SYSTEM, "out of memory in add2blist");
	}
	blist.bl[i].hd = hd; blist.bl[i].bi = bi;
	blist.bl[i].nr = nr; blist.bl[i].nc = 0;
}


int16 *
viewbeams(vp, hr, vr, blp)		/* convert view into sections/beams */
VIEW	*vp;
int	hr, vr;
BEAMLIST	*blp;
{
	static int16	sectlist[HDMAX+1];
	int16	sectarr[MINRES+1][MINRES+1];
	double	d0, d1, mindist;
	GCOORD	gc[2];
	FVECT	rorg, rdir;
	int	shr, svr, sampquant;
	int	v;
	register int	h, hd;
						/* clear section flags */
	memset((char *)sectlist, '\0', sizeof(sectlist));
						/* identify view sections */
	for (v = 0; v <= MINRES; v++)
		for (h = 0; h <= MINRES; h++) {
			sectarr[v][h] = -1;
			mindist = 0.99*FHUGE;
			if (viewray(rorg, rdir, vp, (double)h/MINRES,
					(double)v/MINRES) < -FTINY)
				continue;
			for (hd = 0; hdlist[hd] != NULL; hd++) {
				d0 = hdinter(gc, NULL, &d1,
						hdlist[hd], rorg, rdir);
				if (d0 >= 0.99*FHUGE)
					continue;		/* missed */
				if (d0 <= 0. && d1 >= 0.) {
					sectarr[v][h] = hd;
					break;			/* inside */
				}
				if (d0 > 0. && d0 < mindist) {
					sectarr[v][h] = hd;
					mindist = d0;
				} else if (d1 < 0. && -d1 < mindist) {
					sectarr[v][h] = hd;
					mindist = -d1;
				}
			}
			if ((hd = sectarr[v][h]) >= 0)
				sectlist[hd]++;		/* flag section */
		}
						/* convert flags to list */
	for (h = hd = 0; hdlist[hd] != NULL; h++, hd++) {
		while (!sectlist[hd])
			if (hdlist[++hd] == NULL)
				goto loopexit;
		sectlist[h] = hd;
	}
loopexit:
	sectlist[h] = -1;			/* list terminator */
	if (blp == NULL)			/* if no beam list... */
		return(sectlist);			/* return early */
						/* else set up sampling */
	if (hr|vr && hr*vr <= NVSAMPS) {
		shr = hr; svr = vr;
		sampquant = 1;
	} else {
		shr = sqrt(NVSAMPS/viewaspect(vp)) + .5;
		svr = (NVSAMPS + shr/2)/shr;
		sampquant = (hr*vr + shr*svr/2)/(shr*svr);
	}
	blist.bl = NULL; blist.nb = 0;		/* sample view rays */
	for (v = svr; v--; )
		for (h = shr; h--; ) {
			hd = random()>>6 & 03;		/* get section */
			hd = sectarr[v*MINRES/svr + (hd&01)]
					[h*MINRES/shr + (hd>>1)];
			if (hd < 0)
				continue;
							/* get sample ray */
			if (viewray(rorg, rdir, vp, (h+frandom())/shr,
					(v+frandom())/svr) < -FTINY)
				continue;
#ifdef OJITTER
							/* jitter origin */
			d0 = OJITTER*(frandom() - .5) / hdlist[hd]->grid[0];
			VSUM(rorg, rorg, hdlist[hd]->xv[0], d0);
			d0 = OJITTER*(frandom() - .5) / hdlist[hd]->grid[1];
			VSUM(rorg, rorg, hdlist[hd]->xv[1], d0);
			d0 = OJITTER*(frandom() - .5) / hdlist[hd]->grid[2];
			VSUM(rorg, rorg, hdlist[hd]->xv[2], d0);
#endif
							/* intersect section */
			if (hdinter(gc, NULL, NULL, hdlist[hd], rorg, rdir)
					>= 0.99*FHUGE)
				continue;
							/* add samples */
			add2blist(hd, hdbindex(hdlist[hd],gc), sampquant);
		}
	*blp = blist;		/* transfer beam list */
	return(sectlist);			/* all done! */
}


int
nextview(vp, fp)			/* get next view from fp */
VIEW	*vp;
FILE	*fp;
{
	char	linebuf[256];

	while (fgets(linebuf, sizeof(linebuf), fp) != NULL)
		if (isview(linebuf) && sscanview(vp, linebuf) > 0)
			return(0);
	return(EOF);
}	
