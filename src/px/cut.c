#ifndef lint
static const char	RCSid[] = "$Id: cut.c,v 2.3 2003/07/27 22:12:03 schorsch Exp $";
#endif
/*
CUT - Median bisection (k-d tree) algorithm for color image quantization

flaw: doesn't always create as many colors as requested
because the binary subdivision is always done in a balanced tree and when a box
contains only one color it's unsplittable

Paul Heckbert
*/

#include "ciq.h"

struct index {short *o; int freq;};
/* offset and pixel count for a box */

struct pail {struct plum *first; int freak;};
/* each pail is a bucket containing plums */

struct plum {struct plum *next; int code;};
/* a color in a pail's linked list */

static short off[len];
/* color codes: offsets of colors in hist array */


makecm(nw,na)			/* subdivide colorspace to generate a colormap*/
int nw,*na;			/* number of colors wanted, actual */
{
    int i,m,n,freq,weight,sr,sg,sb;
    short *o;
    struct index index[257];	/* pointers into off */
    /* index[i].o is pointer to code for first color in box i */

    for (o=off, *na=i=0; i<len; i++)
	if (hist[i]) {
	    (*na)++;
	    *o++ = i;			/* initialize offsets */
	}
    index[0].o = off; /* start with one big box which contains all *na colors */
    index[0].freq = xmax*ymax;		/* and (sum of hist[i]) pixels */
    index[1].o = off+*na;

    /* breadth-first binary subdivision: try to make nw boxes */
    for (m=1; m<nw; m<<=1) {
	/*printf("%d ",m);*/
	index[m<<1].o = index[m].o;
	for (i=m-1; i>=0; --i) splitbox(index+i,index+(i<<1));
    }
    for (n=i=0; i<m && n<nw; i++) if (index[i+1].o>index[i].o) {
	sr = sg = sb = (freq = index[i].freq)>>1;
	for (o=index[i].o; o<index[i+1].o; o++) {
	    weight = hist[*o];
	    /* find weighted average of colors in box i */
	    sr += red(*o)*weight;
	    sg += gre(*o)*weight;
	    sb += blu(*o)*weight;
	}
	color[0][n] = sr/freq;
	color[1][n] = sg/freq;
	color[2][n] = sb/freq;
	/* printf("color[%d]=(%3d,%3d,%3d) freq=%d\n",
	    n,color[0][n],color[1][n],color[2][n],freq); */
	n++;
    }
    /*printf("[%d empty boxes] from %d to %d colors\n",nw-n,*na,n);*/
    return n;
}

splitbox(ii,io)		/* split a box in two: half of the pixels from */
struct index *ii,*io;	/* box ii will go into io, the other half into io+1 */
{
    register short *o,*o1,*o2;
    register int shift,count,k,freq,r,g,b,r1,g1,b1,r2,g2,b2;
    struct pail pail[32],*p,*q;			/* buckets for sorting colors */
    struct plum plum[len],*pl=plum,*u;		/* colors */

    o1 = ii[0].o;
    o2 = ii[1].o;
    freq = ii[0].freq;		/* number of pixels with color in this box */
    r1 = g1 = b1 = 256; r2 = g2 = b2 = -1;
    for (o=o1; o<o2; o++) {
	r = red(*o);
	g = gre(*o);
	b = blu(*o);
	if (r<r1) r1 = r; if (r>r2) r2 = r;	/* find min&max r, g, and b */
	if (g<g1) g1 = g; if (g>g2) g2 = g;
	if (b<b1) b1 = b; if (b>b2) b2 = b;
    }

    /* adaptive partitioning: find dominant dimension */
    shift = r2-r1>g2-g1?r2-r1>=b2-b1?10:0:g2-g1>=b2-b1?5:0;
    /* printf("splitbox %3d-%3d %3dx%3dx%3d %c ",
	o1-off,o2-off,r2-r1,g2-g1,b2-b1,shift==10?'R':shift==5?'G':'B'); */

    for (p=pail, k=0; k<32; k++, p++) {		/* start with empty pails */
	p->first = 0;
	p->freak = 0;
    }
    for (o=o1; o<o2; o++) {
	p = pail+(*o>>shift&31);		/* find appropriate pail */
	pl->code = *o;
	pl->next = p->first;
	p->first = pl++;			/* add new plum to pail */
	p->freak += hist[*o];
    }

    /* find median along dominant dimension */
    for (count=freq>>1, p=pail; count>0; p++) count -= p->freak;

    if (p>pail && count+(p-1)->freak<-count)	/* back up one? */
	count += (--p)->freak;

    io[1].o = o1;			/* in case of empty box */
    for (o=o1, q=pail; o<o2;) {	/* dump the pails to reorder colors in box */
	for (u=q->first; u; u=u->next) *o++ = u->code;
	if (++q==p) io[1].o = o;	/* place partition at offset o */
    }
    io[0].o = o1;
    io[0].freq = (freq>>1)-count;
    io[1].freq = ((freq+1)>>1)+count;
    /* printf(" at %3d %d=%d+%d\n",io[1].o-off,freq,io[0].freq,io[1].freq); */
}
