#ifndef lint
static const char	RCSid[] = "$Id: metacalls.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  metacalls.c - functional interface to metafile.
 *
 *     2/24/86
 */

#include  "meta.h"


#define  RIGHT		0
#define  UP		1
#define  LEFT		2
#define  DOWN		3

#define  pflush()	if (inpoly) closepoly()

#define  putdec(v)	if (v) decival(v); else *cap++ = '0';


static int  curx = 0;
static int  cury = 0;
static int  cura0 = 0;
static int  inpoly = FALSE;
static char  curargs[MAXARGS] = "";
static char  *cap;


mendpage()			/* end of page */
{
    pflush();
    pglob(PEOP, 0200, NULL);
}


mdone()				/* end of graphics metafile */
{
    pflush();
    pglob(PEOF, 0200, NULL);
}


minclude(fname)			/* include a file */
char  *fname;
{
    pflush();
    pglob(PINCL, 1, fname);
}


msetpat(pn, pat)		/* set a pattern */
int  pn;
char  *pat;
{
    pflush();
    pglob(PSET, pn+4, pat);
}


mopenseg(sname)			/* open a segment */
char  *sname;
{
    pflush();
    pglob(POPEN, 0, sname);
}


mcloseseg()			/* close current segment */
{
    pflush();
    pglob(PCLOSE, 0200, NULL);
}


mline(x, y, type, thick, color)		/* start a line */
int  x, y;
int  type, thick, color;
{
    pflush();
    cura0 = (type<<4 & 060) | (thick<<2 & 014) | (color & 03);
    curx = x;
    cury = y;
}


mrectangle(xmin, ymin, xmax, ymax, pat, color)	/* fill a rectangle */
int  xmin, ymin, xmax, ymax;
int  pat, color;
{
    pflush();
    cura0 = (pat<<2 & 014) | (color & 03);
    pprim(PRFILL, cura0, xmin, ymin, xmax, ymax, NULL);
}


mtriangle(xmin, ymin, xmax, ymax, d, pat, color)	/* fill a triangle */
int  xmin, ymin, xmax, ymax;
int  d, pat, color;
{
    pflush();
    cura0 = (let_dir(d)<<4 & 060) | (pat<<2 & 014) | (color & 03);
    pprim(PTFILL, cura0, xmin, ymin, xmax, ymax, NULL);
}


mpoly(x, y, border, pat, color)		/* start a polygon */
int  x, y;
int  border, pat, color;
{
    pflush();
    cura0 = (border<<6 & 0100) | (pat<<2 & 014) | (color & 03);
    cap = curargs;
    inpoly = TRUE;
    putdec(x);
    *cap++ = ' ';
    putdec(y);
}


mtext(x, y, s, cpi, color)		/* matrix string */
int  x, y;
char  *s;
int  cpi;
int  color;
{
    pflush();
    cura0 = (color & 03);
    if (cpi < 10) {
	cura0 += 04;
	cpi *= 2;
    }
    if (cpi > 11)
	cura0 += 020;
    if (cpi > 14)
	cura0 += 020;
    if (cpi > 18)
	cura0 += 020;
    pprim(PMSTR, cura0, x, y, x, y, s);
}


mvstr(xmin, ymin, xmax, ymax, s, d, thick, color)	/* vector string */
int  xmin, ymin, xmax, ymax;
char  *s;
int  d, thick, color;
{
    pflush();
    cura0 = (let_dir(d)<<4 & 060) | (thick<<2 & 014) | (color & 03);
    pprim(PVSTR, cura0, xmin, ymin, xmax, ymax, s);
}


msegment(xmin, ymin, xmax, ymax, sname, d, thick, color)	/* segment */
int  xmin, ymin, xmax, ymax;
char  *sname;
int  d, thick, color;
{
    pflush();
    cura0 = (let_dir(d)<<4 & 060) | (thick<<2 & 014) | (color & 03);
    pprim(PSEG, cura0, xmin, ymin, xmax, ymax, sname);
}


mdraw(x, y)				/* draw to next point */
int  x, y;
{
    if (inpoly) {
	polyval(x, y);
    } else if (x != curx || y != cury) {
	plseg(cura0, curx, cury, x, y);
	curx = x;
	cury = y;
    }
}


static
decival(v)				/* add value to polygon */
register int  v;
{
    if (!v)
	return;
    decival(v/10);
    *cap++ = v%10 + '0';
}


static
polyval(x, y)				/* add vertex to a polygon */
register int  x, y;
{
    *cap++ = ' ';
    putdec(x);
    *cap++ = ' ';
    putdec(y);
}


static
closepoly()				/* close current polygon */
{
    *cap = '\0';
    pprim(PPFILL, cura0, 0, 0, XYSIZE-1, XYSIZE-1, curargs);
    inpoly = FALSE;
}


static int
let_dir(c)		/* convert letter to corresponding direction */
register int  c;
{
    switch (c) {
    case 'R':
    case 'r':
	return(RIGHT);
    case 'U':
    case 'u':
	return(UP);
    case 'L':
    case 'l':
	return(LEFT);
    case 'D':
    case 'd':
	return(DOWN);
    }
    return(0);
}
