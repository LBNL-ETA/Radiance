#ifndef lint
static const char	RCSid[] = "$Id: clumpbeams.c,v 3.7 2003/11/17 02:21:53 greg Exp $";
#endif
/*
 * Bundle holodeck beams together into clumps.
 */

#include "holo.h"

#define flgop(p,i,op)		((p)[(i)>>5] op (1L<<((i)&0x1f)))
#define issetfl(p,i)		flgop(p,i,&)
#define setfl(p,i)		flgop(p,i,|=)
#define clrfl(p,i)		flgop(p,i,&=~)

static int	bneighlist[9*9-1];
static int	bneighrem;

#define nextneigh()	(bneighrem<=0 ? 0 : bneighlist[--bneighrem])


gcshifti(gc, ia, di, hp)	/* shift cell row or column */
register GCOORD	*gc;
int	ia, di;
register HOLO	*hp;
{
	int	nw;

	if (di > 0) {
		if (++gc->i[ia] >= hp->grid[((gc->w>>1)+1+ia)%3]) {
			nw = ((gc->w&~1) + (ia<<1) + 3) % 6;
			gc->i[ia] = gc->i[1-ia];
			gc->i[1-ia] = gc->w&1 ? hp->grid[((nw>>1)+2-ia)%3]-1 : 0;
			gc->w = nw;
		}
	} else if (di < 0) {
		if (--gc->i[ia] < 0) {
			nw = ((gc->w&~1) + (ia<<1) + 2) % 6;
			gc->i[ia] = gc->i[1-ia];
			gc->i[1-ia] = gc->w&1 ? hp->grid[((nw>>1)+2-ia)%3]-1 : 0;
			gc->w = nw;
		}
	}
}


mkneighgrid(ng, hp, gc)		/* compute neighborhood for grid cell */
GCOORD	ng[3*3];
HOLO	*hp;
GCOORD	*gc;
{
	GCOORD	gci0;
	register int	i, j;

	for (i = 3; i--; ) {
		gci0 = *gc;
		gcshifti(&gci0, 0, i-1, hp);
		for (j = 3; j--; ) {
			*(ng+(3*i+j)) = gci0;
			gcshifti(ng+(3*i+j), gci0.w==gc->w, j-1, hp);
		}
	}
}


static int
firstneigh(hp, b)		/* initialize neighbor list and return first */
HOLO	*hp;
int	b;
{
	GCOORD	wg0[9], wg1[9], bgc[2];
	int	i, j;

	hdbcoord(bgc, hp, b);
	mkneighgrid(wg0, hp, bgc);
	mkneighgrid(wg1, hp, bgc+1);
	bneighrem = 0;
	for (i = 9; i--; )
		for (j = 9; j--; ) {
			if ((i == 4) & (j == 4))	/* don't copy starting beam */
				continue;
			if (wg0[i].w == wg1[j].w)
				continue;
			*bgc = *(wg0+i);
			*(bgc+1) = *(wg1+j);
			bneighlist[bneighrem++] = hdbindex(hp, bgc);
#ifdef DEBUG
			if (bneighlist[bneighrem-1] <= 0)
				error(CONSISTENCY, "bad beam in firstneigh");
#endif
		}
	return(nextneigh());
}


clumpbeams(hp, maxcnt, maxsiz, cf)	/* clump beams from hinp */
register HOLO	*hp;
int	maxcnt, maxsiz;
int	(*cf)();
{
	static short	primes[] = {9431,6803,4177,2659,1609,887,587,251,47,1};
	uint32	*bflags;
	int	*bqueue;
	int	bqlen;
	int32	bqtotal;
	int	bc, bci, bqc, myprime;
	register int	i;
					/* get clump size */
	if (maxcnt <= 1)
		maxcnt = nbeams(hp);
	maxsiz /= sizeof(RAYVAL);
					/* allocate beam queue */
	bqueue = (int *)malloc(maxcnt*sizeof(int));
	bflags = (uint32 *)calloc((nbeams(hp)>>5)+1,
			sizeof(uint32));
	if ((bqueue == NULL) | (bflags == NULL))
		error(SYSTEM, "out of memory in clumpbeams");
					/* mark empty beams as done */
	for (i = nbeams(hp); i > 0; i--)
		if (!bnrays(hp, i))
			setfl(bflags, i);
					/* pick a good prime step size */
	for (i = 0; primes[i]<<5 >= nbeams(hp); i++)
		;
	while ((myprime = primes[i++]) > 1)
		if (nbeams(hp) % myprime)
			break;
					/* add each input beam and neighbors */
	for (bc = bci = nbeams(hp); bc > 0; bc--,
			bci += bci>myprime ? -myprime : nbeams(hp)-myprime) {
		if (issetfl(bflags, bci))
			continue;
		bqueue[0] = bci;		/* initialize queue */
		bqlen = 1;
		bqtotal = bnrays(hp, i);
		setfl(bflags, bci);
						/* run through growing queue */
		for (bqc = 0; bqc < bqlen; bqc++) {
						/* add neighbors until full */
			for (i = firstneigh(hp,bqueue[bqc]); i > 0;
					i = nextneigh()) {
				if (issetfl(bflags, i))	/* done already? */
					continue;
				bqueue[bqlen++] = i;	/* add it */
				bqtotal += bnrays(hp, i);
				setfl(bflags, i);
				if (bqlen >= maxcnt ||
						(maxsiz && bqtotal >= maxsiz))
					break;		/* queue full */
			}
			if (i > 0)
				break;
		}
		(*cf)(hp, bqueue, bqlen);	/* transfer clump */
	}
					/* all done; clean up */
	free((void *)bqueue);
	free((void *)bflags);
}
