#ifndef lint
static const char	RCSid[] = "$Id: biq.c,v 2.3 2003/05/13 17:58:33 greg Exp $";
#endif
/*
 *  biq.c - simple greyscale quantization.
 *
 *	9/19/88
 */

#include "standard.h"
#include "ciq.h"


biq(dith,nw,synth,cm)
int dith;		/* is dithering desired? 0=no, 1=yes */
int nw;			/* number of colors wanted in output image */
int synth;		/* synthesize colormap? 0=no, 1=yes */
colormap cm;		/* quantization colormap */
			/* read if synth=0; always written */
{
    colormap ocm;

    picreadcm(ocm);	/* read original picture's colormap (usually identity)*/

    for (n = 0; n < nw; n++)			/* also sets n to nw */
	color[0][n] = color[1][n] = color[2][n] = (n*256L+128)/nw;

    picwritecm(color);

    draw_grey(ocm);

    bcopy((void *)color,(void *)cm,sizeof color);
}

/*----------------------------------------------------------------------*/

draw_grey(ocm)
colormap ocm;
{
    register rgbpixel *linin;
    register pixel *linout;
    rgbpixel intmp;
    int outtmp;
    int y;
    register int x;

    linin = line3alloc(xmax);
    linout = linealloc(xmax);

    for (y = 0; y < ymax; y++) {
	picreadline3(y, linin);
	for (x = 0; x < xmax; x++) {
		intmp.r = ocm[0][linin[x].r];
		intmp.g = ocm[1][linin[x].g];
		intmp.b = ocm[2][linin[x].b];
		outtmp = rgb_bright(&intmp);
		linout[x] = (outtmp*n+n/2)/256;
	}
	picwriteline(y, linout);
    }
    free((void *)linin);
    free((void *)linout);
}
