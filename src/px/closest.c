#ifndef lint
static const char	RCSid[] = "$Id: closest.c,v 2.3 2003/07/27 22:12:03 schorsch Exp $";
#endif
/*
CLOSEST - nearest-color lookup by locally ordered search
we use distance in rgb space as color metric

fillbucket subroutine makes trimmed, sorted lists of
possible closest colors to speed the closest searching.
About 6 times faster than exhaustive.

Paul Heckbert
*/

#include "ciq.h"

#define inf 3*256*256
#define key(r,g,b) (((r)&0xe0)<<1|((g)&0xe0)>>2|((b)&0xe0)>>5)   
    /* 9-bit hash key: rrrgggbbb */

static struct thing {int mindist,pv;} space[512*257],*next,*bucket[512];
					/* sorted lists of colors */

static int nfill,tests,calls;
static char filled[512];		/* has bucket[key] been filled yet? */
static int sq[511];

initializeclosest() {			/* reset the buckets */
    int k;
    nfill = tests = calls = 0;
    for (k=0; k<512; k++) filled[k] = 0;
    next = space;
}

closest(r,g,b)		/* find pv of colormap color closest to (r,g,b) */
int r,g,b;
{
    register struct thing *p;
    register int *rsq,*gsq,*bsq,dist,min;
    int best,k;
    struct thing *p0;

    k = key(r,g,b);
    if (!filled[k]) fillbucket(k);
    min = inf;
    rsq = sq+255-r;
    gsq = sq+255-g;
    bsq = sq+255-b;

    /* stop looking when best is closer than next could be */
    for (p0=p=bucket[k]; min>p->mindist; p++) {
	dist = rsq[color[0][p->pv]] + gsq[color[1][p->pv]] +
	    bsq[color[2][p->pv]];
	if (dist<min) {
	    best = p->pv;
	    min = dist;
	}
    }
    tests += p-p0;
    calls++;
    return best;
}

#define H 16	/* half-width of a bucket */

compare(a,b)
register struct thing *a,*b;
{
    return a->mindist-b->mindist;
}

fillbucket(k)	/* make list of colormap colors which could be closest */
int k;		/* matches for colors in bucket #k */
{
    register int r,g,b,j,dist,*rsq,*gsq,*bsq,min,best;
    struct thing *p,*q;

    if (!sq[0]) for (j= -255; j<256; j++) sq[j+255] = j*j;

    r = (k>>1&0xe0)|H;			/* center of 32x32x32 cubical bucket */
    g = (k<<2&0xe0)|H;
    b = (k<<5&0xe0)|H;
    rsq = sq+255-r;
    gsq = sq+255-g;
    bsq = sq+255-b;
    bucket[k] = p = next;
    min = inf;

    for (j=0; j<n; j++, p++) {
	p->pv = j;
	dist = 0;			/* compute distance from cube surface */
	     if (color[0][j]< r-H) dist += rsq[color[0][j]+H];
	else if (color[0][j]>=r+H) dist += rsq[color[0][j]-H+1];
	     if (color[1][j]< g-H) dist += gsq[color[1][j]+H];
	else if (color[1][j]>=g+H) dist += gsq[color[1][j]-H+1];
	     if (color[2][j]< b-H) dist += bsq[color[2][j]+H];
	else if (color[2][j]>=b+H) dist += bsq[color[2][j]-H+1];
	p->mindist = dist;		/* for termination test in closest */
	dist = rsq[color[0][j]]+gsq[color[1][j]]+bsq[color[2][j]];
	if (dist<min) {			/* find color closest to cube center */
	    best = j;
	    min = dist;
	}
    }

    dist = rsq[color[0][best]+(color[0][best]<r?1-H:H)] +
	   gsq[color[1][best]+(color[1][best]<g?1-H:H)] +
	   bsq[color[2][best]+(color[2][best]<b?1-H:H)];
    /* dist is an upper bound on mindist: maximum possible distance */
    /* between a color in the bucket and its nearest neighbor */

    /* eliminate colors which couldn't be closest */
    for (p=q=next, j=0; j<n; j++, p++)
	if (p->mindist<dist) {
	    q->mindist = p->mindist;
	    q->pv = p->pv;
	    q++;
	}
    if (q<=next) fprintf(stderr,"ERROR: empty bucket!\n");

    /* and sort the list by mindist */
    qsort(next,q-next,sizeof(struct thing),compare);
    q->mindist = inf;				/* list terminator */
    next = q+1;
    filled[k] = 1;
    nfill++;
}

/*endclosest() {
    printf("filled %d buckets, averaging %d colors per list, %d of which were tested\n",
	nfill,(next-space-nfill)/nfill,tests/calls);
}*/
