#ifndef lint
static const char	RCSid[] = "$Id: colortab.c,v 2.6 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 * colortab.c - allocate and control dynamic color table.
 *
 *	We start off with a uniform partition of color space.
 *	As pixels are sent to the frame buffer, a histogram is built.
 *	When a new color table is requested, the histogram is used
 *	to make a pseudo-optimal partition, after which the
 *	histogram is cleared.  This algorithm
 *	performs only as well as the next drawing's color
 *	distribution is correlated to the last.
 *
 *  External symbols declared in drvier.h
 */

#include "copyright.h"

#include "standard.h"

#include "color.h"
				/* histogram resolution */
#define NRED		24
#define NGRN		32
#define NBLU		16
#define HMAX		NGRN
				/* minimum box count for adaptive partition */
#define MINSAMP		7
				/* maximum distance^2 before color reassign */
#define MAXDST2		12
				/* map a color */
#define map_col(c,p)	clrmap[p][ colval(c,p)<1. ? \
				(int)(colval(c,p)*256.) : 255 ]
				/* color partition tree */
#define CNODE		short
#define set_branch(p,c) ((c)<<2|(p))
#define set_pval(pv)	((pv)<<2|3)
#define is_branch(cn)	(((cn)&3)!=3)
#define is_pval(cn)	(((cn)&3)==3)
#define part(cn)	((cn)>>2)
#define prim(cn)	((cn)&3)
#define pval(cn)	((cn)>>2)
				/* our color table */
static struct tabent {
	long	sum[3];		/* sum of colors using this entry */
	int	n;		/* number of colors */
	BYTE	ent[3];		/* current table value */
}	*clrtab = NULL;
				/* color cube partition */
static CNODE	*ctree = NULL;
				/* our color correction map */
static BYTE	clrmap[3][256];
				/* histogram of colors used */
static unsigned short	histo[NRED][NGRN][NBLU];
				/* initial color cube boundary */
static int	CLRCUBE[3][2] = {{0,NRED},{0,NGRN},{0,NBLU}};

static int	split();
static void	cut();


int
new_ctab(ncolors)		/* start new color table with max ncolors */
int	ncolors;
{
	int	treesize;

	if (ncolors < 1)
		return(0);
				/* free old tables */
	if (clrtab != NULL)
		free((void *)clrtab);
	if (ctree != NULL)
		free((void *)ctree);
				/* get new tables */
	for (treesize = 1; treesize < ncolors; treesize <<= 1)
		;
	treesize <<= 1;
	clrtab = (struct tabent *)calloc(ncolors, sizeof(struct tabent));
	ctree = (CNODE *)malloc(treesize*sizeof(CNODE));
	if (clrtab == NULL || ctree == NULL)
		return(0);
				/* partition color space */
	cut(ctree, 0, CLRCUBE, 0, ncolors);
				/* clear histogram */
	bzero((char *)histo, sizeof(histo));
				/* return number of colors used */
	return(ncolors);
}


int
get_pixel(col, set_pixel)	/* get pixel for color */
COLOR	col;
void	(*set_pixel)();
{
	int	r, g, b;
	int	cv[3];
	register CNODE	*tp;
	register int	h;
						/* map color */
	r = map_col(col,RED);
	g = map_col(col,GRN);
	b = map_col(col,BLU);
						/* reduce resolution */
	cv[RED] = (r*NRED)>>8;
	cv[GRN] = (g*NGRN)>>8;
	cv[BLU] = (b*NBLU)>>8;
						/* add to histogram */
	histo[cv[RED]][cv[GRN]][cv[BLU]]++;
						/* find pixel in tree */
	for (tp = ctree, h = 0; is_branch(*tp); h++)
		if (cv[prim(*tp)] < part(*tp))
			tp += 1<<h;		/* left branch */
		else
			tp += 1<<(h+1);		/* right branch */
	h = pval(*tp);
						/* add to color table */
	clrtab[h].sum[RED] += r;
	clrtab[h].sum[GRN] += g;
	clrtab[h].sum[BLU] += b;
	clrtab[h].n++;
					/* recompute average */
	r = clrtab[h].sum[RED] / clrtab[h].n;
	g = clrtab[h].sum[GRN] / clrtab[h].n;
	b = clrtab[h].sum[BLU] / clrtab[h].n;
					/* check for movement */
	if (clrtab[h].n == 1 ||
			(r-clrtab[h].ent[RED])*(r-clrtab[h].ent[RED]) +
			(g-clrtab[h].ent[GRN])*(g-clrtab[h].ent[GRN]) +
			(b-clrtab[h].ent[BLU])*(b-clrtab[h].ent[BLU]) > MAXDST2) {
		clrtab[h].ent[RED] = r;
		clrtab[h].ent[GRN] = g; /* reassign pixel */
		clrtab[h].ent[BLU] = b;
#ifdef DEBUG
		sprintf(errmsg, "pixel %d = (%d,%d,%d) (%d refs)\n",
				h, r, g, b, clrtab[h].n);
		eputs(errmsg);
#endif
		(*set_pixel)(h, r, g, b);
	}
	return(h);				/* return pixel value */
}


void
make_gmap(gam)			/* make gamma correction map */
double	gam;
{
	register int	i;
	
	for (i = 0; i < 256; i++)
		clrmap[RED][i] =
		clrmap[GRN][i] =
		clrmap[BLU][i] = 256.0 * pow((i+0.5)/256.0, 1.0/gam);
}


void
set_cmap(rmap, gmap, bmap)	/* set custom color correction map */
BYTE	*rmap, *gmap, *bmap;
{
	bcopy((char *)rmap, (char *)clrmap[RED], 256);
	bcopy((char *)gmap, (char *)clrmap[GRN], 256);
	bcopy((char *)bmap, (char *)clrmap[BLU], 256);
}


void
map_color(rgb, col)		/* map a color to a byte triplet */
BYTE	rgb[3];
COLOR	col;
{
	rgb[RED] = map_col(col,RED);
	rgb[GRN] = map_col(col,GRN);
	rgb[BLU] = map_col(col,BLU);
}


static void
cut(tree, level, box, c0, c1)		/* partition color space */
register CNODE	*tree;
int	level;
register int	box[3][2];
int	c0, c1;
{
	int	kb[3][2];
	
	if (c1-c0 <= 1) {		/* assign pixel */
		*tree = set_pval(c0);
		return;
	}
					/* split box */
	*tree = split(box);
	bcopy((char *)box, (char *)kb, sizeof(kb));
						/* do left (lesser) branch */
	kb[prim(*tree)][1] = part(*tree);
	cut(tree+(1<<level), level+1, kb, c0, (c0+c1)>>1);
						/* do right branch */
	kb[prim(*tree)][0] = part(*tree);
	kb[prim(*tree)][1] = box[prim(*tree)][1];
	cut(tree+(1<<(level+1)), level+1, kb, (c0+c1)>>1, c1);
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
