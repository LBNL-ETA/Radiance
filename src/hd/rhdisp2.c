#ifndef lint
static const char	RCSid[] = "$Id: rhdisp2.c,v 3.36 2003/07/27 22:12:02 schorsch Exp $";
#endif
/*
 * Holodeck beam tracking for display process
 */

#include "rholo.h"
#include "rhdisp.h"
#include "rhdriver.h"

#ifndef MEYERNG
#define MEYERNG		0.2	/* target mean eye range (rel. to grid) */
#endif

#define CBEAMBLK	1024	/* cbeam allocation block size */

static PACKHEAD	*cbeam = NULL;	/* current beam list */
static int	ncbeams = 0;	/* number of sorted beams in cbeam */
static int	xcbeams = 0;	/* extra (unregistered) beams past ncbeams */
static int	maxcbeam = 0;	/* size of cbeam array */

VIEWPOINT	cureye;		/* current eye position */


int
newcbeam()		/* allocate new entry at end of cbeam array */
{
	int	i;

	if ((i = ncbeams + xcbeams++) >= maxcbeam) {	/* grow array */
		maxcbeam += CBEAMBLK;
		if (cbeam == NULL)
			cbeam = (PACKHEAD *)malloc(
					maxcbeam*sizeof(PACKHEAD) );
		else
			cbeam = (PACKHEAD *)realloc((void *)cbeam,
					maxcbeam*sizeof(PACKHEAD) );
		if (cbeam == NULL)
			error(SYSTEM, "out of memory in newcbeam");
	}
	return(i);
}


int
cbeamcmp(cb1, cb2)	/* compare two cbeam entries for sort: keep orphans */
register PACKHEAD	*cb1, *cb2;
{
	register int	c;

	if ((c = cb1->bi - cb2->bi))	/* sort on beam index first */
		return(c);
	return(cb1->hd - cb2->hd);	/* use hd to resolve matches */
}


int
cbeamcmp2(cb1, cb2)	/* compare two cbeam entries for sort: no orphans */
register PACKHEAD	*cb1, *cb2;
{
	register int	c;

	if (!cb1->nr)			/* put orphans at the end, unsorted */
		return(cb2->nr);
	if (!cb2->nr)
		return(-1);
	if ((c = cb1->bi - cb2->bi))	/* sort on beam index first */
		return(c);
	return(cb1->hd - cb2->hd);	/* use hd to resolve matches */
}


int
findcbeam(hd, bi)	/* find the specified beam in our sorted list */
int	hd, bi;
{
	PACKHEAD	cb;
	register PACKHEAD	*p;

	if (ncbeams <= 0)
		return(-1);
	cb.hd = hd; cb.bi = bi; cb.nr = cb.nc = 0;
	p = (PACKHEAD *)bsearch((char *)&cb, (char *)cbeam, ncbeams,
			sizeof(PACKHEAD), cbeamcmp);
	if (p == NULL)
		return(-1);
	return(p - cbeam);
}


int
getcbeam(hd, bi)	/* get the specified beam, allocating as necessary */
register int	hd;
int	bi;
{
	register int	n;
				/* first, look in sorted list */
	if ((n = findcbeam(hd, bi)) >= 0)
		return(n);
				/* linear search through xcbeams to be sure */
	for (n = ncbeams+xcbeams; n-- > ncbeams; )
		if (cbeam[n].bi == bi && cbeam[n].hd == hd)
			return(n);
				/* check legality */
	if ((hd < 0) | (hd >= HDMAX) || hdlist[hd] == NULL)
		error(INTERNAL, "illegal holodeck number in getcbeam");
	if ((bi < 1) | (bi > nbeams(hdlist[hd])))
		error(INTERNAL, "illegal beam index in getcbeam");
	n = newcbeam();		/* allocate and assign */
	cbeam[n].hd = hd; cbeam[n].bi = bi; cbeam[n].nr = cbeam[n].nc = 0;
	return(n);
}


cbeamsort(adopt)	/* sort our beam list, possibly turning out orphans */
int	adopt;
{
	register int	i;

	if (!(ncbeams += xcbeams))
		return;
	xcbeams = 0;
	qsort((char *)cbeam, ncbeams, sizeof(PACKHEAD),
			adopt ? cbeamcmp : cbeamcmp2);
	if (adopt)
		return;
	for (i = ncbeams; i--; )	/* identify orphans */
		if (cbeam[i].nr)
			break;
	xcbeams = ncbeams - ++i;	/* put orphans after ncbeams */
	ncbeams = i;
}


beam_init(fresh)		/* clear beam list for new view(s) */
int	fresh;
{
	register int	i;

	if (fresh)			/* discard old beams? */
		ncbeams = xcbeams = 0;
	else				/* else clear sample requests */
		for (i = ncbeams+xcbeams; i--; )
			cbeam[i].nr = 0;
	cureye.rng = 0.;
}


int16 *
beam_view(vn, hr, vr)		/* add beam view (if advisable) */
VIEW	*vn;
int	hr, vr;
{
	int16	*slist;
	BEAMLIST	blist;
	double	eravg, d;
	register HOLO	*hp;
	register int	i, n;
					/* compute beams for view */
	slist = viewbeams(vn, hr, vr, &blist);
	if (*slist < 0) {
		error(COMMAND, "no sections visible from this view");
		return(NULL);
	}
					/* sort current beam list */
	cbeamsort(1);
					/* add new beams to list */
	for (i = blist.nb; i--; ) {
		n = getcbeam(blist.bl[i].hd, blist.bl[i].bi);
		if (blist.bl[i].nr > cbeam[n].nr)
			cbeam[n].nr = blist.bl[i].nr;
	}
	free((void *)blist.bl);		/* free list */
	if (MEYERNG <= FTINY)
		return(slist);
					/* compute average eye range */
	eravg = 0.;
	for (n = 0; slist[n] >= 0; n++) {
		hp = hdlist[slist[n]];
		eravg += 	MEYERNG/3. / VLEN(hp->wg[0]) +
				MEYERNG/3. / VLEN(hp->wg[1]) +
				MEYERNG/3. / VLEN(hp->wg[2]) ;
	}
	eravg /= (double)n;
					/* add to current eye position */
	if (cureye.rng <= FTINY) {
		VCOPY(cureye.vpt, vn->vp);
		cureye.rng = eravg;
	} else if ((d = sqrt(dist2(vn->vp,cureye.vpt))) + eravg > cureye.rng) {
		for (i = 3; i--; )
			cureye.vpt[i] = 0.5*(cureye.vpt[i] + vn->vp[i]);
		cureye.rng = 0.5*(cureye.rng + eravg + d);
	}
	return(slist);
}


int
beam_sync(all)			/* update beam list on server */
int	all;
{
					/* set new eye position */
	serv_request(DR_VIEWPOINT, sizeof(VIEWPOINT), (char *)&cureye);
					/* sort list (put orphans at end) */
	cbeamsort(all < 0);
					/* send beam request */
	if (all) {
		if (ncbeams > 0)
			serv_request(DR_NEWSET,
					ncbeams*sizeof(PACKHEAD), cbeam);
	} else {
		if (ncbeams+xcbeams > 0)
			serv_request(DR_ADJSET,
				(ncbeams+xcbeams)*sizeof(PACKHEAD), cbeam);
	}
	xcbeams = 0;			/* truncate our list */
	return(ncbeams);
}


gridlines(f)			/* run through holodeck section grid lines */
int	(*f)();
{
	register int	hd, w, i;
	int	g0, g1;
	FVECT	wp[2], mov;
	double	d;
					/* do each wall on each section */
	for (hd = 0; hdlist[hd] != NULL; hd++)
		for (w = 0; w < 6; w++) {
			g0 = hdwg0[w];
			g1 = hdwg1[w];
			d = 1.0/hdlist[hd]->grid[g0];
			mov[0] = d * hdlist[hd]->xv[g0][0];
			mov[1] = d * hdlist[hd]->xv[g0][1];
			mov[2] = d * hdlist[hd]->xv[g0][2];
			if (w & 1) {
				VSUM(wp[0], hdlist[hd]->orig,
						hdlist[hd]->xv[w>>1], 1.);
				VSUM(wp[0], wp[0], mov, 1.);
			} else
				VCOPY(wp[0], hdlist[hd]->orig);
			VSUM(wp[1], wp[0], hdlist[hd]->xv[g1], 1.);
			for (i = hdlist[hd]->grid[g0]; ; ) {	/* g0 lines */
				(*f)(wp);
				if (!--i) break;
				wp[0][0] += mov[0]; wp[0][1] += mov[1];
				wp[0][2] += mov[2]; wp[1][0] += mov[0];
				wp[1][1] += mov[1]; wp[1][2] += mov[2];
			}
			d = 1.0/hdlist[hd]->grid[g1];
			mov[0] = d * hdlist[hd]->xv[g1][0];
			mov[1] = d * hdlist[hd]->xv[g1][1];
			mov[2] = d * hdlist[hd]->xv[g1][2];
			if (w & 1)
				VSUM(wp[0], hdlist[hd]->orig,
						hdlist[hd]->xv[w>>1], 1.);
			else
				VSUM(wp[0], hdlist[hd]->orig, mov, 1.);
			VSUM(wp[1], wp[0], hdlist[hd]->xv[g0], 1.);
			for (i = hdlist[hd]->grid[g1]; ; ) {	/* g1 lines */
				(*f)(wp);
				if (!--i) break;
				wp[0][0] += mov[0]; wp[0][1] += mov[1];
				wp[0][2] += mov[2]; wp[1][0] += mov[0];
				wp[1][1] += mov[1]; wp[1][2] += mov[2];
			}
		}
}
