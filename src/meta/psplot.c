#ifndef lint
static const char	RCSid[] = "$Id: psplot.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  psplot.c - routines for generating PostScript output.
 *
 *	9/23/88
 */

#include "meta.h"

#include <ctype.h>

#ifdef SMALLAREA
#define HMARGIN		(0.2*72)		/* horizontal margin */
#define VMARGIN		(0.2*72)		/* vertical margin */
#define PWIDTH		(6*72-2*HMARGIN)	/* width of device */
#define PHEIGHT		(6*72-2*VMARGIN)	/* height of device */
#else
#define HMARGIN		(.6*72)			/* horizontal margin */
#define VMARGIN		(.6*72)			/* vertical margin */
#define PWIDTH		(8.5*72-2*HMARGIN)	/* width of device */
#define PHEIGHT		(11*72-2*VMARGIN)	/* height of device */
#endif

#define PSQUARE		PWIDTH			/* maximal square */

#define MAXPATHL	700			/* maximum path length */

#define PREFIX		"seg_"		/* name prefix */
#define ESCCHAR		'#'		/* to introduce hex in name */

#define checkline()	if (inaline) endline()


static int	savelvl = 0;

static int	inaline = 0;

static char	hexdigit[] = "0123456789ABCDEF";

static char	buf[512];


static char *
convertname(s)				/* convert to segment name */
char	*s;
{
	register char	*cp, *cp0;

	for (cp = buf, cp0 = PREFIX; *cp0; )
		*cp++ = *cp0++;
	for (cp0 = s; *cp0; cp0++) {
		if (isalnum(*cp0))
			*cp++ = *cp0;
		else {
			*cp++ = ESCCHAR;
			*cp++ = hexdigit[*cp0 >> 4];
			*cp++ = hexdigit[*cp0 & 0xf];
		}
	}
	*cp = '\0';
	return(buf);
}


static char *
convertstring(s)		/* convert to acceptable string */
register char	*s;
{
	register char	*cp;
	
	for (cp = buf; *s; s++) {
		if (*s == '(' || *s == ')')
			*cp++ = '\\';
		if (isprint(*s))
			*cp++ = *s;
		else {
			*cp++ = '\\';
			*cp++ = (*s>>6) + '0';
			*cp++ = (*s>>3 & 07) + '0';
			*cp++ = (*s & 07) + '0';
		}
	}
	*cp = '\0';
	return(buf);
}


init(id)			/* initialize */
char	*id;
{
	printf("%%!PS-Adobe-2.0 EPSF-2.0\n");
	printf("%%%%BoundingBox: %.0f %.0f %.0f %.0f\n",
			HMARGIN+(PWIDTH-PSQUARE)/2.,
			VMARGIN+(PHEIGHT-PSQUARE)/2.,
			HMARGIN+(PWIDTH-PSQUARE)/2.+PSQUARE,
			VMARGIN+(PHEIGHT-PSQUARE)/2.+PSQUARE);
	printf("%%%%DocumentFonts: Helvetica Courier\n");
	if (id != NULL)
		printf("%%%%Creator: %s\n", id);
	printf("gsave\n");
	printf("save\n");
	printf("256 dict begin\n");
	printf("%f %f translate\n", HMARGIN+(PWIDTH-PSQUARE)/2.,
			VMARGIN+(PHEIGHT-PSQUARE)/2.);
	printf("%f %f scale\n", (double)PSQUARE/XYSIZE,
			(double)PSQUARE/XYSIZE);
	printf("/mapping matrix currentmatrix def\n");
	printf("/segment {\t\t%% name lth pat col dir xmn ymn xmx ymx\n");
	printf("\t/startstate save def\n");
	printf("\t128 dict begin\t\t%% get parameters\n");
	printf("\t/ymax exch def /xmax exch def\n");
	printf("\t/ymin exch def /xmin exch def\n");
	printf("\t/orient exch def\n");
	printf("\t/colndx exch def /patndx exch def /lthick exch def\n");
	printf("\t%% segment name left on stack\n");
	printf("\txmax xmin add 2 div ymax ymin add 2 div translate\n");
	printf("\txmax xmin sub %d div ymax ymin sub %d div scale\n",
			XYSIZE, XYSIZE);
	printf("\torient 90 mul rotate\n");
	printf("\t%d %d translate\n", -XYSIZE/2, -XYSIZE/2);
	printf("\tcvx exec\t\t\t%% execute segment\n");
	printf("\tend\n");
	printf("\tstartstate restore\n");
	printf("} def\n");
	printf("/vstr {\t\t\t%% string\n");
	printf("\t0 setcolor\n");
	printf("\t/Helvetica findfont setfont\n");
	printf("\tdup stringwidth pop %d exch div %d scale\n", XYSIZE, XYSIZE);
	printf("\tcurrentfont charorigin translate\n");
	printf("\t0 0 moveto show\n");
	printf("} def\n");
	printf("/charorigin {\t\t\t%% font charorigin xorg yorg\n");
	printf("\tdup /FontBBox get aload pop pop pop\n");
	printf("\t3 2 roll /FontMatrix get transform\n");
	printf("\texch neg exch neg\n");
	printf("} def\n");
	printf("/newline {\t\t\t%% lpat lthick col x0 y0\n");
	printf("\tnewpath moveto\n");
	printf("\tsetcolor setlthick setlpat\n");
	printf("} def\n");
	printf("/rfilldict 4 dict def\n");
	printf("/rfill {\t\t\t%% pat col xmn ymn xmx ymx\n");
	printf("\trfilldict begin\n");
	printf("\t/ymax exch def /xmax exch def\n");
	printf("\t/ymin exch def /xmin exch def\n");
	printf("\tsetcolor setpattern\n");
	printf("\tnewpath\n");
	printf("\txmin ymin moveto\n");
	printf("\txmax ymin lineto\n");
	printf("\txmax ymax lineto\n");
	printf("\txmin ymax lineto\n");
	printf("\tclosepath patfill\n");
	printf("\tend\n");
	printf("} def\n");
	printf("/fillpoly {\t\t\t%% pat col mark x1 y1 .. xn yn\n");
	printf("\tnewpath moveto\n");
	printf("\tcounttomark 2 idiv { lineto } repeat cleartomark\n");
	printf("\tclosepath\n");
	printf("\tsetcolor setpattern\n");
	printf("\t{ gsave 0 setlinewidth stroke grestore } if\n");
	printf("\tpatfill\n");
	printf("} def\n");
	printf("/mstr {\t\t\t%% fnt col x y\n");
	printf("\t100 add moveto setcolor setmfont show\n");
	printf("} def\n");
	printf("/patfill {\n");
	printf("\tfill\t\t\t\t%% unimplemented\n");
	printf("} def\n");
	printf("/setmfont {\t\t\t%% fontndx\n");
	printf("\tmfontab exch get setfont\n");
	printf("} def\n");
	printf("/setcolor {\t\t\t%% colndx\n");
	printf("\tcolndx 0 ne { pop colndx } if\n");
	printf("\trgbcoltab exch get aload pop setrgbcolor\n");
	printf("} def\n");
	printf("/setlthick {\t\t\t%% lthick\n");
	printf("\tlthick 0 ne { pop lthick } if\n");
	printf("\tsetlinewidth\n");
	printf("} def\n");
	printf("/setlpat {\t\t\t%% lpatndx\n");
	printf("\tdashtab exch get 0 setdash\n");
	printf("} def\n");
	printf("/setpattern {\t\t\t%% patndx\n");
	printf("\tpop\t\t\t\t%% unimplemented\n");
	printf("} def\n");
	printf("/canonfont\t\t\t%% canonical matrix string font\n");
	printf("\t/Courier findfont\n");
	printf("\tdup charorigin matrix translate\n");
	printf("\tmakefont\n");
	printf("def\n");
	printf("/mfontab [\t\t\t%% hardware font table\n");
	printf("\t[\n");
	printf("\t\t[ 340 0 0 340 0 -340 ]\n");
	printf("\t\t[ 681 0 0 340 0 -340 ]\n");
	printf("\t\t[ 340 0 0 681 0 -681 ]\n");
	printf("\t\t[ 681 0 0 681 0 -681 ]\n");
	printf("\t\t[ 282 0 0 282 0 -282 ]\n");
	printf("\t\t[ 564 0 0 282 0 -282 ]\n");
	printf("\t\t[ 282 0 0 564 0 -564 ]\n");
	printf("\t\t[ 564 0 0 564 0 -564 ]\n");
	printf("\t\t[ 199 0 0 199 0 -199 ]\n");
	printf("\t\t[ 398 0 0 199 0 -199 ]\n");
	printf("\t\t[ 199 0 0 398 0 -398 ]\n");
	printf("\t\t[ 398 0 0 398 0 -398 ]\n");
	printf("\t\t[ 169 0 0 169 0 -169 ]\n");
	printf("\t\t[ 339 0 0 169 0 -169 ]\n");
	printf("\t\t[ 169 0 0 339 0 -339 ]\n");
	printf("\t\t[ 339 0 0 339 0 -339 ]\n");
	printf("\t]\n");
	printf("\t{ canonfont exch makefont }\n");
	printf("\tforall\n");
	printf("] def\n");
	printf("/dashtab [ [ ] [ 200 80 ] [ 80 80 ] [ 200 80 80 80 ] ] def\n");
	printf("/rgbcoltab [ [ 0 0 0 ] [ 1 0 0 ] [ 0 1 0 ] [ 0 0 1 ] ] def\n");
	printf("/colndx 0 def /patndx 0 def /lthick 0 def\n");
	printf("%%%%EndProlog\n");
}


done()				/* done with graphics */
{
	printf("end\nrestore\ngrestore\n");
}


endpage()			/* done with page */
{
	checkline();

	printf("showpage\nmapping setmatrix\n");
}


segopen(s)			/* open a segment */
char	*s;
{
	checkline();

	printf("/%s {\n", convertname(s));
}


segclose()			/* close a segment */
{
	checkline();

	printf("} def\n");
}


doseg(p)			/* instance of a segment */
register PRIMITIVE	*p;
{
	checkline();

	printf("/%s %d %d %d ", convertname(p->args),
			WIDTH(p->arg0>>2 & 03),
			p->arg0>>2 & 03, p->arg0 & 03);
	printf("%d %d %d %d %d segment\n", p->arg0>>4,
			p->xy[XMN], p->xy[YMN], p->xy[XMX], p->xy[YMX]);
}


printstr(p)			/* print a string */
register PRIMITIVE	*p;
{
	checkline();
	printf("(%s) %d %d %d %d mstr\n", convertstring(p->args),
			p->arg0>>2 & 017, p->arg0 & 03,
			p->xy[XMN], p->xy[YMX]);
}


plotvstr(p)			/* print a vector string */
register PRIMITIVE	*p;
{
	checkline();

	printf("(%s) /vstr %d 0 %d ", convertstring(p->args),
			WIDTH(p->arg0>>2 & 03), p->arg0 & 03);
	printf("%d %d %d %d %d segment\n", p->arg0>>4,
			p->xy[XMN], p->xy[YMN], p->xy[XMX], p->xy[YMX]);
}


plotlseg(p)			/* plot a line segment */
register PRIMITIVE	*p;
{
	static int	right = FALSE;
	static int	curx, cury;
	static int	curlpat, curlthick, curlcol;
	int	y1, y2;
	int	lpat, lthick, lcol;

	if (p->arg0 & 0100) {
		y1 = p->xy[YMX];
		y2 = p->xy[YMN];
	} else {
		y1 = p->xy[YMN];
		y2 = p->xy[YMX];
	}
	lpat = p->arg0>>4 & 03;
	lthick = WIDTH(p->arg0>>2 & 03);
	lcol = p->arg0 & 03;
	if (!inaline || lpat != curlpat ||
			lthick != curlthick || lcol != curlcol) {
		checkline();
		printf("%d %d %d %d %d newline\n",
				curlpat = lpat,
				curlthick = lthick,
				curlcol = lcol,
				curx = p->xy[XMN],
				cury = y1);
	}
	if (curx == p->xy[XMN] && cury == y1) {
		printf("%d %d lineto\n", curx = p->xy[XMX], cury = y2);
	} else if (curx == p->xy[XMX] && cury == y2) {
		printf("%d %d lineto\n", curx = p->xy[XMN], cury = y1);
	} else if (right = !right) {
		printf("%d %d moveto ", p->xy[XMN], y1);
		printf("%d %d lineto\n", curx = p->xy[XMX], cury = y2);
	} else {
		printf("%d %d moveto ", p->xy[XMX], y2);
		printf("%d %d lineto\n", curx = p->xy[XMN], cury = y1);
	}
	if (++inaline >= MAXPATHL)
		endline();
}


endline()			/* close current line */
{
	printf("stroke\n");
	inaline = 0;
}


fillrect(p)			/* fill a rectangle */
register PRIMITIVE	*p;
{
	checkline();

	printf("%d %d %d %d %d %d rfill\n", p->arg0>>2 & 03, p->arg0 & 03,
			p->xy[XMN], p->xy[YMN], p->xy[XMX], p->xy[YMX]);
}


filltri(p)			/* fill a triangle */
register PRIMITIVE	*p;
{
	static short	corn[4][2] = {XMN,YMX,XMN,YMN,XMX,YMN,XMX,YMX};
	int		orient;
	register int	i;
	
	checkline();

	printf("false %d %d mark\n", p->arg0>>2 & 03, p->arg0 & 03);
	orient = p->arg0>>4 & 03;
	for (i = 0; i < 4; i++)
		if (i != orient)
			printf("%d %d ", p->xy[corn[i][0]],
					p->xy[corn[i][1]]);
	printf("fillpoly\n");
}


xform(xp, yp, p)		/* transform a point according to p */
register int  *xp, *yp;
register PRIMITIVE  *p;
{
    int  x, y;

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


fillpoly(p)			/* fill a polygon */
register PRIMITIVE	*p;
{
	extern char	*nextscan();
	register char	*s;
	int	curx, cury;

	checkline();

	printf("%s %d %d mark\n", p->arg0 & 0100 ? "false" : "true",
			p->arg0>>2 & 03, p->arg0 & 03);
	s = p->args;
	while ((s = nextscan(nextscan(s, "%d", &curx), "%d", &cury)) != NULL) {
		xform(&curx, &cury, p);
		printf("%d %d ", curx, cury);
	}
	printf("fillpoly\n");
}


set(attrib, val)		/* set an attribute or context */
int	attrib;
char	*val;
{
	checkline();

	switch (attrib) {
	case SALL:
		printf("save\n");
		savelvl++;
		break;
	case SPAT0:
	case SPAT1:
	case SPAT2:
	case SPAT3:
		break;
	default:
		error(WARNING, "illegal set command");
		break;
	}
}


unset(attrib)			/* unset an attribute or context */
int	attrib;
{
	checkline();

	switch (attrib) {
	case SALL:
		if (savelvl > 0) {
			printf("restore\n");
			savelvl--;
		}
		break;
	case SPAT0:
	case SPAT1:
	case SPAT2:
	case SPAT3:
		break;
	default:
		error(WARNING, "illegal unset command");
		break;
	}
}


reset(attrib)			/* reset an attribute or context */
int	attrib;
{
	checkline();

	switch (attrib) {
	case SALL:
		while (savelvl > 0) {
			printf("restore\n");
			savelvl--;
		}
		printf("restore save\n");
		break;
	case SPAT0:
	case SPAT1:
	case SPAT2:
	case SPAT3:
		break;
	default:
		error(WARNING, "illegal reset command");
		break;
	}
}
