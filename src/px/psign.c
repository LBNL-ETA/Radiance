#ifndef lint
static const char	RCSid[] = "$Id: psign.c,v 2.24 2003/10/22 02:06:35 greg Exp $";
#endif
/*
 *  psign.c - produce picture from text.
 *
 *	7/1/87
 */

#include  "copyright.h"

#include  <string.h>

#include  "platform.h"
#include  "standard.h"
#include  "color.h"
#include  "font.h"

#ifndef	 SSS
#define	 SSS			3	/* super-sample size */
#endif

#define	 MAXLINE		512	/* longest allowable line */

char  *fontfile = "helvet.fnt";		/* our font file */

COLOR  bgcolor = WHTCOLOR;		/* background color */
COLOR  fgcolor = BLKCOLOR;		/* foreground color */

int  direct = 'r';			/* direction (right, up, left, down) */

int  cheight = 32*SSS;			/* character height */
double	aspect = 1.67;			/* height/width character aspect */
double	spacing = 0.0;			/* character spacing */
int  cwidth;				/* computed character width */

unsigned char  *ourbitmap;		/* our output bitmap */
int  xsiz, ysiz;			/* bitmap dimensions */
int  xdim;				/* size of horizontal scan (bytes) */

				/* conflicting def's in param.h */
#undef	tstbit
#undef	setbit
#undef	clrbit
#undef	tglbit

#define	 bitop(x,y,op)		(ourbitmap[(y)*xdim+((x)>>3)] op (1<<((x)&7)))
#define	 tstbit(x,y)		bitop(x,y,&)
#define	 setbit(x,y)		bitop(x,y,|=)
#define	 clrbit(x,y)		bitop(x,y,&=~)
#define	 tglbit(x,y)		bitop(x,y,^=)

FONT  *ourfont;				/* our font */

typedef struct line {
	char  *s;		/* line w/o LF */
	short  *sp;		/* character spacing */
	struct line  *next;	/* next line up */
} LINE;

LINE  *ourtext;				/* our text */
int  nlines, maxline;			/* text dimensions */
int  maxwidth;				/* maximum line width (dvi) */


main(argc, argv)
int  argc;
char  *argv[];
{
	int  an;
	SET_FILE_BINARY(stdout);
	for (an = 1; an < argc && argv[an][0] == '-'; an++)
		switch (argv[an][1]) {
		case 'c':			/* color */
			switch (argv[an][2]) {
			case 'f':			/* foreground */
				setcolor(fgcolor, atof(argv[an+1]),
						atof(argv[an+2]),
						atof(argv[an+3]));
				an += 3;
				break;
			case 'b':			/* background */
				setcolor(bgcolor, atof(argv[an+1]),
						atof(argv[an+2]),
						atof(argv[an+3]));
				an += 3;
				break;
			default:
				goto unkopt;
			}
			break;
		case 'f':			/* font */
			fontfile = argv[++an];
			break;
		case 'd':			/* direction */
			switch (argv[an][2]) {
			case 'r':			/* right */
			case 'u':			/* up */
			case 'l':			/* left */
			case 'd':			/* down */
				direct = argv[an][2];
				break;
			default:
				goto unkopt;
			}
			break;
		case 'x':			/* x resolution */
			xsiz = atoi(argv[++an])*SSS;
			break;
		case 'y':
			ysiz = atoi(argv[++an])*SSS;
			break;
		case 'h':			/* height of characters */
			cheight = atoi(argv[++an])*SSS;
			break;
		case 'a':			/* aspect ratio */
			aspect = atof(argv[++an]);
			break;
		case 's':			/* spacing */
			spacing = atof(argv[++an]);
			break;
		default:;
unkopt:
			fprintf(stderr, "%s: unknown option: %s\n",
					argv[0], argv[an]);
			exit(1);
		}
					/* load font file */
	ourfont = getfont(fontfile);
					/* get text */
	if (an == argc)
		gettext(stdin);
	else
		arg_text(argc-an, argv+an);

					/* create bit map */
	makemap();
					/* convert text to bitmap */
	maptext();
					/* print header */
	newheader("RADIANCE", stdout);
	printargs(argc, argv, stdout);
	fputformat(COLRFMT, stdout);
	putchar('\n');
					/* write out bitmap */
	writemap(stdout);

	exit(0);
}


makemap()			/* create the bit map */
{
	double	pictaspect;
	
	if (direct == 'r' || direct == 'l') {
		if (xsiz <= 0) {
			cwidth = cheight/aspect + 0.5;
			xsiz = (long)maxwidth*cwidth >> 8;
			ysiz = nlines*cheight;
		} else if (aspect > FTINY) {
			if (ysiz <= 0)
				ysiz = cheight*nlines;
			pictaspect = 256*nlines*aspect/maxwidth;
			if (pictaspect*xsiz < ysiz)
				ysiz = pictaspect*xsiz + 0.5;
			else
				xsiz = ysiz/pictaspect + 0.5;
			cheight = ysiz/nlines;
			cwidth = cheight/aspect;
		} else {
			if (ysiz <= 0)
				ysiz = cheight*nlines;
			pictaspect = (double)ysiz/xsiz;
			aspect = pictaspect*maxwidth/(256*nlines);
			cheight = ysiz/nlines;
			cwidth = cheight/aspect;
		}
	} else {			/* reverse orientation */
		if (ysiz <= 0) {
			cwidth = cheight/aspect + 0.5;
			xsiz = nlines*cheight;
			ysiz = (long)maxwidth*cwidth >> 8;
		} else if (aspect > FTINY) {
			if (xsiz <= 0)
				xsiz = cheight*nlines;
			pictaspect = maxwidth/(256*nlines*aspect);
			if (pictaspect*xsiz < ysiz)
				ysiz = pictaspect*xsiz + 0.5;
			else
				xsiz = ysiz/pictaspect + 0.5;
			cheight = xsiz/nlines;
			cwidth = cheight/aspect;
		} else {
			if (xsiz <= 0)
				xsiz = cheight*nlines;
			pictaspect = (double)ysiz/xsiz;
			aspect = maxwidth/(256*nlines*pictaspect);
			cheight = xsiz/nlines;
			cwidth = cheight/aspect;
		}
	}
	if (xsiz % SSS)
		xsiz += SSS - xsiz%SSS;
	if (ysiz % SSS)
		ysiz += SSS - ysiz%SSS;
	xdim = (xsiz+7)/8;
	ourbitmap = (BYTE *)bmalloc(ysiz*xdim);
	if (ourbitmap == NULL)
		error(SYSTEM, "Out of memory in makemap");
	memset((char *)ourbitmap, '\0', ysiz*xdim);
}


gettext(fp)			/* get text from a file */
FILE  *fp;
{
	char  buf[MAXLINE];
	register LINE  *curl;
	int  len;

	maxline = 0;
	maxwidth = 0;
	nlines = 0;
	while (fgets(buf, MAXLINE, fp) != NULL) {
		curl = (LINE *)malloc(sizeof(LINE));
		if (curl == NULL)
			goto memerr;
		len = strlen(buf);
		curl->s = (char *)malloc(len);
		curl->sp = (short *)malloc(sizeof(short)*len--);
		if ((curl->s == NULL) | (curl->sp == NULL))
			goto memerr;
		if (len > maxline)
			maxline = len;
		strncpy(curl->s, buf, len);
		curl->s[len] = '\0';
		if (spacing < -1./256.)
			len = squeeztext(curl->sp, curl->s, ourfont,
					(int)(spacing*-256.));
		else if (spacing > 1./256.)
			len = proptext(curl->sp, curl->s, ourfont,
					(int)(spacing*256.), 3);
		else
			len = uniftext(curl->sp, curl->s, ourfont);
		if (len > maxwidth)
			maxwidth = len;
		curl->next = ourtext;
		ourtext = curl;
		nlines++;
	}
	return;
memerr:
	error(SYSTEM, "Out of memory in gettext");
}


arg_text(ac, av)			/* get text from arguments */
int  ac;
char  *av[];
{
	register char  *cp;

	ourtext = (LINE *)malloc(sizeof(LINE));
	if (ourtext == NULL)
		goto memerr;
	ourtext->s = (char *)malloc(MAXLINE);
	if (ourtext->s == NULL)
		goto memerr;
	for (cp = ourtext->s; ac-- > 0; av++) {
		strcpy(cp, *av);
		cp += strlen(*av);
		*cp++ = ' ';
	}
	*--cp = '\0';
	ourtext->next = NULL;
	maxline = strlen(ourtext->s);
	ourtext->sp = (short *)malloc(sizeof(short)*(maxline+1));
	if (ourtext->sp == NULL)
		goto memerr;
	if (spacing < -1./256.)
		maxwidth = squeeztext(ourtext->sp, ourtext->s, ourfont,
				(int)(spacing*-256.));
	else if (spacing > 1./256.)
		maxwidth = proptext(ourtext->sp, ourtext->s, ourfont,
				(int)(spacing*256.), 3);
	else
		maxwidth = uniftext(ourtext->sp, ourtext->s, ourfont);
	nlines = 1;
	return;
memerr:
	error(SYSTEM, "Out of memory in arg_text");
}


maptext()			/* map our text */
{
	register LINE  *curl;
	int  l, len;
	register int  i, c;

	for (l = 0, curl = ourtext; curl != NULL; l += 256, curl = curl->next) {
		len = strlen(curl->s); c = 0;
		for (i = 0; i < len; i++) {
			c += curl->sp[i];
			mapglyph(ourfont->fg[curl->s[i]&0xff], c, l);
		}
	}
}


mapglyph(gl, tx0, ty0)		/* convert a glyph */
GLYPH  *gl;
int  tx0, ty0;
{
	int  n;
	register GORD  *gp;
	int  p0[2], p1[2];

	if (gl == NULL)
		return;

	n = gl->nverts;
	gp = gvlist(gl);
	mapcoord(p0, gp[2*n-2]+tx0, gp[2*n-1]+ty0);
	while (n--) {
		mapcoord(p1, gp[0]+tx0, gp[1]+ty0);
		mapedge(p0[0], p0[1], p1[0]-p0[0], p1[1]-p0[1]);
		p0[0] = p1[0]; p0[1] = p1[1];
		gp += 2;
	}
}


mapcoord(p, tx, ty)		/* map text to picture coordinates */
int  p[2], tx, ty;
{
	tx = (long)tx*cwidth >> 8;
	ty = (long)ty*cheight >> 8;

	switch (direct) {
	case 'r':			/* right */
		p[0] = tx;
		p[1] = ty;
		return;
	case 'u':			/* up */
		p[0] = xsiz-1-ty;
		p[1] = tx;
		return;
	case 'l':			/* left */
		p[0] = xsiz-1-tx;
		p[1] = ysiz-1-ty;
		return;
	case 'd':			/* down */
		p[0] = ty;
		p[1] = ysiz-1-tx;
		return;
	}
}


mapedge(x, y, run, rise)		/* map an edge */
register int  x, y;
int  run, rise;
{
	int  xstep;
	int  rise2, run2;
	int  n;

	if (rise == 0)
		return;
						/* always draw up */
	if (rise < 0) {
		x += run;
		y += rise;
		rise = -rise;
		run = -run;
	}
	if (run < 0) {
		xstep = -1;
		run = -run;
	} else
		xstep = 1;
	n = rise;
	run2 = rise2 = 0;
	while (n)
		if (rise2 >= run2) {
			tglbit(x, y);
			n--;
			y++;
			run2 += run;
		} else {
			x += xstep;
			rise2 += rise;
		}
}


writemap(fp)			/* write out bitmap */
FILE  *fp;
{
	COLR  pixval[SSS*SSS+1];	/* possible pixel values */
	COLOR  ctmp0, ctmp1;
	double	d;
	COLR  *scanout;
	int  x, y;
	register int  i, j;
	int  cnt;
	register int  inglyph;

	fprintf(fp, "-Y %d +X %d\n", ysiz/SSS, xsiz/SSS);

	scanout = (COLR *)malloc(xsiz/SSS*sizeof(COLR));
	if (scanout == NULL)
		error(SYSTEM, "Out of memory in writemap");
	for (i = 0; i <= SSS*SSS; i++) {	/* compute possible values */
		copycolor(ctmp0, fgcolor);
		d = (double)i/(SSS*SSS);
		scalecolor(ctmp0, d);
		copycolor(ctmp1, bgcolor);
		d = 1.0 - d;
		scalecolor(ctmp1, d);
		addcolor(ctmp0, ctmp1);
		setcolr(pixval[i], colval(ctmp0,RED),
				colval(ctmp0,GRN), colval(ctmp0,BLU));
	}
	for (y = ysiz/SSS-1; y >= 0; y--) {
		inglyph = 0;
		for (x = 0; x < xsiz/SSS; x++) {
			cnt = 0;
			for (j = 0; j < SSS; j++)
				for (i = 0; i < SSS; i++) {
					if (tstbit(x*SSS+i, y*SSS+j))
						inglyph ^= 1<<j;
					if (inglyph & 1<<j)
						cnt++;
				}
			copycolr(scanout[x], pixval[cnt]);
		}
		if (fwritecolrs(scanout, xsiz/SSS, fp) < 0) {
			fprintf(stderr, "write error in writemap\n");
			exit(1);
		}
	}
	free((void *)scanout);
}
