#ifndef lint
static const char	RCSid[] = "$Id: clrtab.c,v 2.12 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 * Simple median-cut color quantization based on colortab.c
 */

#include "standard.h"

#include "color.h"
				/* histogram resolution */
#define NRED		36
#define NGRN		48
#define NBLU		24
#define HMAX		NGRN
				/* minimum box count for adaptive partition */
#define MINSAMP		7
				/* color partition */
#define set_branch(p,c) ((c)<<2|(p))
#define part(cn)	((cn)>>2)
#define prim(cn)	((cn)&3)
				/* our color table (global) */
extern BYTE	clrtab[256][3];
				/* histogram of colors / color assignments */
static unsigned	histo[NRED][NGRN][NBLU];
#define cndx(c)		histo[((c)[RED]*NRED)>>8][((c)[GRN]*NGRN)>>8][((c)[BLU]*NBLU)>>8]
				/* initial color cube boundary */
static int	CLRCUBE[3][2] = {0,NRED,0,NGRN,0,NBLU};
				/* maximum propagated error during dithering */
#define MAXERR		20
				/* define CLOSEST to get closest colors */
#ifndef CLOSEST
#ifdef SPEED
#if  SPEED > 8
#define CLOSEST		1	/* this step takes a little longer */
#endif
#endif
#endif

static	cut(), mktabent(), closest(), addneigh(), setclosest();
static int	split();
static unsigned	dist();


new_histo(n)		/* clear our histogram */
int	n;
{
	bzero((char *)histo, sizeof(histo));
	return(0);
}


cnt_pixel(col)		/* add pixel to our histogram */
register BYTE	col[];
{
	cndx(col)++;
}


cnt_colrs(cs, n)		/* add a scanline to our histogram */
register COLR	*cs;
register int	n;
{
	while (n-- > 0) {
		cndx(cs[0])++;
		cs++;
	}
}


new_clrtab(ncolors)		/* make new color table using ncolors */
int	ncolors;
{
	if (ncolors < 1)
		return(0);
	if (ncolors > 256)
		ncolors = 256;
				/* partition color space */
	cut(CLRCUBE, 0, ncolors);
#ifdef CLOSEST
	closest(ncolors);	/* ensure colors picked are closest */
#endif
				/* reset dithering function */
	dith_colrs((BYTE *)NULL, (COLR *)NULL, 0);
				/* return new color table size */
	return(ncolors);
}


int
map_pixel(col)			/* get pixel for color */
register BYTE	col[];
{
    return(cndx(col));
}


map_colrs(bs, cs, n)		/* convert a scanline to color index values */
register BYTE	*bs;
register COLR	*cs;
register int	n;
{
	while (n-- > 0) {
		*bs++ = cndx(cs[0]);
		cs++;
	}
}


dith_colrs(bs, cs, n)		/* convert scanline to dithered index values */
register BYTE	*bs;
register COLR	*cs;
int	n;
{
	static short	(*cerr)[3] = NULL;
	static int	N = 0;
	int	err[3], errp[3];
	register int	x, i;

	if (n != N) {		/* get error propogation array */
		if (N) {
			free((void *)cerr);
			cerr = NULL;
		}
		if (n)
			cerr = (short (*)[3])malloc(3*n*sizeof(short));
		if (cerr == NULL) {
			N = 0;
			map_colrs(bs, cs, n);
			return;
		}
		N = n;
		bzero((char *)cerr, 3*N*sizeof(short));
	}
	err[0] = err[1] = err[2] = 0;
	for (x = 0; x < n; x++) {
		for (i = 0; i < 3; i++) {	/* dither value */
			errp[i] = err[i];
			err[i] += cerr[x][i];
#ifdef MAXERR
			if (err[i] > MAXERR) err[i] = MAXERR;
			else if (err[i] < -MAXERR) err[i] = -MAXERR;
#endif
			err[i] += cs[x][i];
			if (err[i] < 0) err[i] = 0;
			else if (err[i] > 255) err[i] = 255;
		}
		bs[x] = cndx(err);
		for (i = 0; i < 3; i++) {	/* propagate error */
			err[i] -= clrtab[bs[x]][i];
			err[i] /= 3;
			cerr[x][i] = err[i] + errp[i];
		}
	}
}


static
cut(box, c0, c1)			/* partition color space */
register int	box[3][2];
int	c0, c1;
{
	register int	branch;
	int	kb[3][2];
	
	if (c1-c0 <= 1) {		/* assign pixel */
		mktabent(c0, box);
		return;
	}
					/* split box */
	branch = split(box);
	bcopy((char *)box, (char *)kb, sizeof(kb));
						/* do left (lesser) branch */
	kb[prim(branch)][1] = part(branch);
	cut(kb, c0, (c0+c1)>>1);
						/* do right branch */
	kb[prim(branch)][0] = part(branch);
	kb[prim(branch)][1] = box[prim(branch)][1];
	cut(kb, (c0+c1)>>1, c1);
}


static int
split(box)				/* find median cut for box */
register int	box[3][2];
{
#define c0	r
	register int	r, g, b;
	int	pri;
	long	t[HMAX], med;
					/* find dominant axis */
	pri = RED;
	if (box[GRN][1]-box[GRN][0] > box[pri][1]-box[pri][0])
		pri = GRN;
	if (box[BLU][1]-box[BLU][0] > box[pri][1]-box[pri][0])
		pri = BLU;
					/* sum histogram over box */
	med = 0;
	switch (pri) {
	case RED:
		for (r = box[RED][0]; r < box[RED][1]; r++) {
			t[r] = 0;
			for (g = box[GRN][0]; g < box[GRN][1]; g++)
				for (b = box[BLU][0]; b < box[BLU][1]; b++)
					t[r] += histo[r][g][b];
			med += t[r];
		}
		break;
	case GRN:
		for (g = box[GRN][0]; g < box[GRN][1]; g++) {
			t[g] = 0;
			for (b = box[BLU][0]; b < box[BLU][1]; b++)
				for (r = box[RED][0]; r < box[RED][1]; r++)
					t[g] += histo[r][g][b];
			med += t[g];
		}
		break;
	case BLU:
		for (b = box[BLU][0]; b < box[BLU][1]; b++) {
			t[b] = 0;
			for (r = box[RED][0]; r < box[RED][1]; r++)
				for (g = box[GRN][0]; g < box[GRN][1]; g++)
					t[b] += histo[r][g][b];
			med += t[b];
		}
		break;
	}
	if (med < MINSAMP)		/* if too sparse, split at midpoint */
		return(set_branch(pri,(box[pri][0]+box[pri][1])>>1));
					/* find median position */
	med >>= 1;
	for (c0 = box[pri][0]; med > 0; c0++)
		med -= t[c0];
	if (c0 > (box[pri][0]+box[pri][1])>>1)	/* if past the midpoint */
		c0--;				/* part left of median */
	return(set_branch(pri,c0));
#undef c0
}


static
mktabent(p, box)	/* compute average color for box and assign */
int	p;
register int	box[3][2];
{
	unsigned long	sum[3];
	unsigned	r, g;
	unsigned long	n;
	register unsigned	b, c;
						/* sum pixels in box */
	n = 0;
	sum[RED] = sum[GRN] = sum[BLU] = 0;
	for (r = box[RED][0]; r < box[RED][1]; r++)
	    for (g = box[GRN][0]; g < box[GRN][1]; g++)
	    	for (b = box[BLU][0]; b < box[BLU][1]; b++) {
	    	    if (c = histo[r][g][b]) {
			n += c;
			sum[RED] += (long)c*r;
			sum[GRN] += (long)c*g;
			sum[BLU] += (long)c*b;
		    }
		    histo[r][g][b] = p;		/* assign pixel */
		}
	if (n >= (1L<<23)/HMAX) {		/* avoid overflow */
		sum[RED] /= n;
		sum[GRN] /= n;
		sum[BLU] /= n;
		n = 1;
	}
	if (n) {				/* compute average */
		clrtab[p][RED] = sum[RED]*256/NRED/n;
		clrtab[p][GRN] = sum[GRN]*256/NGRN/n;
		clrtab[p][BLU] = sum[BLU]*256/NBLU/n;
	} else {				/* empty box -- use midpoint */
		clrtab[p][RED] = (box[RED][0]+box[RED][1])*256/NRED/2;
		clrtab[p][GRN] = (box[GRN][0]+box[GRN][1])*256/NGRN/2;
		clrtab[p][BLU] = (box[BLU][0]+box[BLU][1])*256/NBLU/2;
	}
}


#ifdef CLOSEST
#define NBSIZ		32
static
closest(n)			/* make sure we have the closest colors */
int	n;
{
	BYTE	*neigh[256];
	register int	r, g, b;
#define i r
					/* get space for neighbor lists */
	for (i = 0; i < n; i++) {
		if ((neigh[i] = (BYTE *)malloc(NBSIZ)) == NULL) {
			while (i--)
				free(neigh[i]);
			return;			/* ENOMEM -- abandon effort */
		}
		neigh[i][0] = i;		/* identity is terminator */
	}
					/* make neighbor lists */
	for (r = 0; r < NRED; r++)
	    for (g = 0; g < NGRN; g++)
		for (b = 0; b < NBLU; b++) {
		    if (r < NRED-1 && histo[r][g][b] != histo[r+1][g][b])
			addneigh(neigh, histo[r][g][b], histo[r+1][g][b]);
		    if (g < NGRN-1 && histo[r][g][b] != histo[r][g+1][b])
			addneigh(neigh, histo[r][g][b], histo[r][g+1][b]);
		    if (b < NBLU-1 && histo[r][g][b] != histo[r][g][b+1])
			addneigh(neigh, histo[r][g][b], histo[r][g][b+1]);
		}
					/* assign closest values */
	for (r = 0; r < NRED; r++)
	    for (g = 0; g < NGRN; g++)
		for (b = 0; b < NBLU; b++)
		    setclosest(neigh, r, g, b);
					/* free neighbor lists */
	for (i = 0; i < n; i++)
		free(neigh[i]);
#undef i
}


static
addneigh(nl, i, j)		/* i and j are neighbors; add them to list */
register BYTE	*nl[];
register int	i;
int	j;
{
	int	nc;
	char	*nnl;
	register int	t;
	
	for (nc = 0; nc < 2; nc++) {		/* do both neighbors */
		for (t = 0; nl[i][t] != i; t++)
			if (nl[i][t] == j)
				break;		/* in list already */
		if (nl[i][t] == i) {		/* add to list */
			nl[i][t++] = j;
			if (t % NBSIZ == 0) {	/* enlarge list */
				if ((nnl = realloc(nl[i], t+NBSIZ)) == NULL)
					t--;
				else
					nl[i] = (BYTE *)nnl;
			}
			nl[i][t] = i;		/* terminator */
		}
		t = i; i = j; j = t;		/* swap and do it again */
	}
}


static unsigned
dist(col, r, g, b)		/* find distance from clrtab entry to r,g,b */
register BYTE	col[3];
int	r, g, b;
{
	register int	tmp;
	register unsigned	sum;
	
	tmp = col[RED]*NRED/256 - r;
	sum = tmp*tmp;
	tmp = col[GRN]*NGRN/256 - g;
	sum += tmp*tmp;
	tmp = col[BLU]*NBLU/256 - b;
	sum += tmp*tmp;
	return(sum);
}


static
setclosest(nl, r, g, b)		/* find index closest to color and assign */
BYTE	*nl[];
int	r, g, b;
{
	int	ident;
	unsigned	min;
	register unsigned	d;
	register BYTE	*p;
					/* get starting value */
	min = dist(clrtab[ident=histo[r][g][b]], r, g, b);
					/* find minimum */
	for (p = nl[ident]; *p != ident; p++)
		if ((d = dist(clrtab[*p], r, g, b)) < min) {
			min = d;
			histo[r][g][b] = *p;
		}
}
#endif
