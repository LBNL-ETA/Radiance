#ifndef lint
static const char	RCSid[] = "$Id: implot.c,v 1.2 2003/11/15 02:13:37 schorsch Exp $";
#endif
/*
 *   Impress plotting functions
 *
 *      1/2/86
 */


#include  "meta.h"

#include  "plot.h"

#include  "imPfuncs.h"

#include  "imPcodes.h"

#include  <ctype.h>


#define  mapx(x)  CONV(x,dxsize)
#define  mapy(y)  CONV(y,dysize)

#define  DPI  296		/* dots per inch for imagen */

#define  TexFamily  90		/* family for texture glyphs */

				/* operation codes for graphics */
#define  whiteOp  0
#define  shadeOp  3
#define  orOp  7
#define  blackOp  15

				/* imagen resident fonts */
static struct {
	char  *name;		/* font name */
	int  iwspace;		/* interword spacing */
	short  loaded;		/* boolean true if font loaded */
} font[] = {
	{"cour07", 21, FALSE},
	{"cour08", 24, FALSE},
	{"cour10", 30, FALSE},
	{"cour12", 36, FALSE},
	{"cour14", 42, FALSE},
	{"zurm20", 60, FALSE}
};

				/* map from our matrix string types */
static int  fontmap[16] = {3,4,4,5,2,4,4,5,1,2,2,3,0,2,2,3};

int  dxsize, dysize;		/* page size */

static void settex(int  pp);
static void setpen(register int  ps);



void
imInit(void)		/* initialize imagen */
{

    imout = stdout;

    fputs("@document(language imPress)", imout);
    fputs("@document(pagereversal off)", imout);

    imSetPum(0);			/* new path deletes old one */
    dxsize = dysize = 8 * DPI;		/* square on 8.5 X 11 in. paper */
    imSetAbsH(DPI/4);			/* .25 in. on left and right */
    imSetAbsV(dysize);
    imSetHVSystem(3, 3, 4);		/* conventional orientation */

}


void
printstr(		/* output a string */
	register PRIMITIVE  *p
)

{
    static int  curfont = -1;
    int  fn;
    register char  *cp;

    fn = fontmap[(p->arg0 >> 2) & 017];		/* get font number */
						/* set font */
    if (fn != curfont) {
	imSetFamily(fn);
	if (!font[fn].loaded) {
	    imCreateFamilyTable(fn, 1, 0, font[fn].name);
	    font[fn].loaded = TRUE;
	}
	imSetSp(font[fn].iwspace);
	curfont = fn;
    }

    imSetAbsH(mapx(p->xy[XMN]));	/* set page position */
    imSetAbsV(mapy(p->xy[YMN]));

    for (cp = p->args; *cp; cp++)	/* write out the string */
	if (isspace(*cp))
	    imSp();
	else
	    im_putbyte(*cp);

}




void
plotlseg(		/* plot a line segment */
	register PRIMITIVE  *p
)

{
    int  x1, x2, y1, y2;
    short  pp;

    pp = (p->arg0 >> 4) & 03;
    if (p->arg0 & 0100 && pp != 0)
	pp += 03;
    settex(pp);

    setpen((p->arg0 >> 2) & 03);

    x1 = mapx(p->xy[XMN]);
    x2 = mapx(p->xy[XMX]);

    if (p->arg0 & 0100) {
	y1 = mapy(p->xy[YMX]);
	y2 = mapy(p->xy[YMN]);
    } else {
	y1 = mapy(p->xy[YMN]);
	y2 = mapy(p->xy[YMX]);
    }

    imCreatePath(2, x1, y1, x2, y2);
    imDrawPath(orOp);

}



void
setfill(			/* set filling mode */
	register int  a0
)

{

    settex(pati[(a0 >> 2) & 03]);

}


void
setpen(			/* set pen size to ps */
	register int  ps
)

{
    static int  curpsiz = -1;		/* current pen size */

    if (ps == curpsiz)
	return;

    curpsiz = ps;
    ps = WIDTH(ps);
    ps = CONV(ps, dxsize);
    if (ps < 2)
	ps = 2;
    else if (ps > 20)
	ps = 20;
    imSetPen(ps);

}


void
settex(			/* set texture pattern to pp */
	int  pp
)

{
    static int  curtex = -1;		/* current texture */
    static short  ploaded[NPATS];	/* booleans for loaded patterns */
    char  span[4];
    int  i, k;
    register int  j;

    if (pp == curtex)			/* already set */
	return;

    if (pp == 0) {			/* all black */
	imSetTexture(0, 0);		/* special case */
	curtex = 0;
	return;
    }

    if (pp >= NPATS || !ploaded[pp]) {	/* download texture glyph */
	im_77(imP_BGLY, TexFamily, pp);	/* family and member */
	im_putword(32);			/* advance-width */
	im_putword(32);			/* width */
	im_putword(0);			/* left-offset */
	im_putword(32);			/* height */
	im_putword(0);			/* top-offset */
	for (i = PATSIZE-1; i >= 0; i--) {
	    for (j = 0; j < 4; j++)
		span[j] = 0;
	    for (j = 0; j < 32; j++)
		span[j>>3] |= ((pattern[pp][i>>3][j/(32/PATSIZE)]
				    >> (i&07)) & 01) << (7-(j&07));
	    for (k = 0; k < 32/PATSIZE; k++)
		for (j = 0; j < 4; j++)
		    im_putbyte(span[j]);
	}
	if (pp < NPATS)
	    ploaded[pp]++;
    }
    imSetTexture(TexFamily, pp);
    curtex = pp;

}


void
fillrect(			/* fill a rectangle */
	register PRIMITIVE  *p
)

{
    int  left, right, top, bottom;

    left = mapx(p->xy[XMN]);
    right = mapx(p->xy[XMX]);
    top = mapy(p->xy[YMX]);
    bottom = mapy(p->xy[YMN]);

    setfill(p->arg0);
    imCreatePath(4, left, bottom, right, bottom, right, top, left, top);
    imFillPath(orOp);

}


void
filltri(			/* fill a triangle */
	register PRIMITIVE  *p
)

{
    int  left, right, top, bottom;

    left = mapx(p->xy[XMN]);
    right = mapx(p->xy[XMX]);
    top = mapy(p->xy[YMX]);
    bottom = mapy(p->xy[YMN]);

    setfill(p->arg0);

    switch (p->arg0 & 060) {
    case 0:			/* right (& down) */
	imCreatePath(3, left, bottom, right, bottom, right, top);
	break;
    case 020:			/* up */
	imCreatePath(3, right, bottom, right, top, left, top);
	break;
    case 040:			/* left */
	imCreatePath(3, right, top, left, top, left, bottom);
	break;
    case 060:			/* down */
	imCreatePath(3, left, top, left, bottom, right, bottom);
	break;
    }
    
    imFillPath(orOp);

}



void
xform(		/* transform a point according to p */
	register int  *xp,
	register int  *yp,
	register PRIMITIVE  *p
)

{
    int  x = 0, y = 0;

    switch (p->arg0 & 060) {
	case 0:			/* right */
	    x = *xp;
	    y = *yp;
	    break;
	case 020:		/* up */
	    x = (XYSIZE-1) - *yp;
	    y = *xp;
	    break;
	case 040:		/* left */
	    x = (XYSIZE-1) - *xp;
	    y = (XYSIZE-1) - *yp;
	    break;
	case 060:		/* down */
	    x = *yp;
	    y = (XYSIZE-1) - *xp;
	    break;
    }
    
    *xp = CONV(x, p->xy[XMX] - p->xy[XMN]) + p->xy[XMN];
    *yp = CONV(y, p->xy[YMX] - p->xy[YMN]) + p->xy[YMN];

}


void
fillpoly(			/* fill a polygon */
	register PRIMITIVE  *p
)

{
    int  points[128];
    register int  *pp;
    register char  *s;

    pp = points;

    if ((s = nextscan(nextscan(p->args,"%d",(char*)pp),"%d",(char*)pp+1))==NULL)
        error(USER, "illegal polygon spec in fillpoly");
    
    xform(pp, pp+1, p);
    pp[0] = mapx(pp[0]); pp[1] = mapy(pp[1]);
    pp += 2;
    
    while ((s = nextscan(nextscan(s,"%d",(char*)pp),"%d",(char*)pp+1))!=NULL) {
        xform(pp, pp+1, p);
	pp[0] = mapx(pp[0]); pp[1] = mapy(pp[1]);
	pp += 2;
    }
    pp[0] = points[0]; pp[1] = points[1];
    pp += 2;
    
    settex(p->arg0 & 077);
    imCreatePathV((pp-points)/2, points);
    imFillPath(orOp);

    if (p->arg0 & 0100) {		/* draw border */
	settex(0);
	setpen(0);
	imDrawPath(orOp);
    }
    
}
