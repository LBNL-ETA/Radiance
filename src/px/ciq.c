#ifndef lint
static const char	RCSid[] = "$Id: ciq.c,v 2.4 2003/06/30 14:59:12 schorsch Exp $";
#endif
/*
CIQ - main program for color image quantization
options for Floyd-Steinberg dithering by minimization of accumulated error
or undithered quantization

Paul Heckbert	16 April 82, cleaned up 8 June 86
Greg Ward	1 March 88, modified for arbitrary picture sizes
*/

#include <string.h>

#include "standard.h"
#include "ciq.h"

#define table(m,r,g,b) hist[m[0][r]|m[1][g]|m[2][b]]  /* histogram/pvtable */

int hist[len];		/* list of frequencies or pixelvalues for coded color */

colormap color;		/* quantization colormap */
int n;			/* number of colors in it */

#define ERRMAX 20	/* maximum allowed propagated error,
			   if >=255 then no clamping */

/*----------------------------------------------------------------------*/

ciq(dith,nw,synth,cm)
int dith;		/* is dithering desired? 0=no, 1=yes */
int nw;			/* number of colors wanted in output image */
int synth;		/* synthesize colormap? 0=no, 1=yes */
colormap cm;		/* quantization colormap */
			/* read if synth=0; always written */
{
    int na,i;
    colormap ocm;

    picreadcm(ocm);	/* read original picture's colormap (usually identity)*/
    sample(ocm);	/* find histogram */

    if (synth)
	n = makecm(nw,&na);	/* analyze histogram and synthesize colormap */
    else {
	memcpy((void *)color,(void *)cm,sizeof color);
	n = nw;
	na = 0;
	for (i=0; i<len; i++) if (hist[i]) na++;
	/*printf("from %d to %d colors\n",na,n);*/
    }
    picwritecm(color);

    initializeclosest();
    if (dith)
	draw_dith(ocm);
    else {
	for (i=0; i<len; i++)
	    if (hist[i]) hist[i] = closest(red(i),gre(i),blu(i));
	draw_nodith(ocm);
    }

    memcpy((void *)cm,(void *)color,sizeof color);
    /*endclosest();*/
}

/*----------------------------------------------------------------------*/

sample(ocm)		/* make color frequency table (histogram) */
colormap ocm;
{
    register int x,y;
    register rgbpixel *lp;
    rgbpixel *line;
    colormap map;

    for (x = 0; x < len; x++)	/* clear histogram */
	hist[x] = 0;

    line = line3alloc(xmax);
    convertmap(ocm,map);

    for (y=0; y<ymax; y++) {
	picreadline3(y,line);
	for (lp=line, x=0; x<xmax; x++, lp++)
	    table(map,lp->r,lp->g,lp->b)++;
    }
    free((void *)line);
}

convertmap(in,out)		/* convert to shifted color map */
register colormap in,out;
{
    register int j;

    for (j=0; j<256; j++) {
	out[0][j] = (in[0][j]&0xf8)<<7;
	out[1][j] = (in[1][j]&0xf8)<<2;
	out[2][j] = (in[2][j]&0xf8)>>3;
    }
}

draw_nodith(ocm)	/* quantize without dithering */
colormap ocm;			/* colormap for orig */
{
    register int x,y;
    register rgbpixel *lp;
    rgbpixel *iline;
    pixel *oline;
    colormap map;

    iline = line3alloc(xmax);
    oline = linealloc(xmax);
    convertmap(ocm,map);

    for (y=0; y<ymax; y++) {
	picreadline3(y,iline);
	for (lp=iline, x=0; x<xmax; x++, lp++)
	    oline[x] = table(map,lp->r,lp->g,lp->b);
	picwriteline(y,oline);
    }
    free((void *)iline);
    free((void *)oline);
}

draw_dith(ocm)			/* quantize with dithering */
colormap ocm;			/* colormap for orig */
{
    register int *p,*q,rr,gg,bb,v,x,y,*P,*Q,*R;
    int *buf;
    rgbpixel *iline;		/* input scanline */
    pixel *oline;		/* output scanline */

    buf = (int *)ecalloc(2,3*sizeof(int)*(xmax+1));
    iline = line3alloc(xmax);
    oline = linealloc(xmax);
    /* reuse hist array as quantization table */
    for (v=0; v<len; v++) hist[v] = -1;

    P = &buf[0];
    Q = &buf[3*(xmax+1)];
    for (y=0; y<ymax; y++, R=P, P=Q, Q=R) {
	Q[0] = Q[1] = Q[2] = 0;
	picreadline3(y,iline);
	for (p=P, q=Q, x=0; x<xmax; x++, p+=3, q+=3) {
	    rr = p[0];
	    if (rr<-ERRMAX) rr = -ERRMAX; else if (rr>ERRMAX) rr = ERRMAX;
	    gg = p[1];
	    if (gg<-ERRMAX) gg = -ERRMAX; else if (gg>ERRMAX) gg = ERRMAX;
	    bb = p[2];
	    if (bb<-ERRMAX) bb = -ERRMAX; else if (bb>ERRMAX) bb = ERRMAX;
	    /* now rr,gg,bb is propagated color */
	    rr += ocm[0][iline[x].r];
	    gg += ocm[1][iline[x].g];	/* ideal */
	    bb += ocm[2][iline[x].b];
	    if (rr<0) rr = 0; else if (rr>255) rr = 255;
	    if (gg<0) gg = 0; else if (gg>255) gg = 255;
	    if (bb<0) bb = 0; else if (bb>255) bb = 255;
	    v = (rr&0xf8)<<7 | (gg&0xf8)<<2 | (bb&0xf8)>>3;
	    if (hist[v]<0)
		hist[v] = closest(rr,gg,bb);	/* nearest to ideal */
	    oline[x] = v = hist[v];
	    rr -= color[0][v];
	    gg -= color[1][v];			/* error */
	    bb -= color[2][v];
	    v = rr*3>>3; p[3] += v; q[0] += v; q[3] = rr-(v<<1);
	    v = gg*3>>3; p[4] += v; q[1] += v; q[4] = gg-(v<<1);
	    v = bb*3>>3; p[5] += v; q[2] += v; q[5] = bb-(v<<1);
	}
	picwriteline(y,oline);
    }
    free((void *)iline);
    free((void *)oline);
    free((void *)buf);
}
