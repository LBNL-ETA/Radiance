#ifndef lint
static const char	RCSid[] = "$Id: metacalls.c,v 1.5 2004/03/22 02:24:23 greg Exp $";
#endif
/*
 *  metacalls.c - functional interface to metafile.
 *
 *     2/24/86
 */

#include  "meta.h"
#include  "plot.h"


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


static int let_dir( register int  c);
static void closepoly(void);
static void polyval( register int  x, register int y);
static void decival( register int  v);


void
mendpage(void)			/* end of page */
{
    pflush();
    pglob(PEOP, 0200, NULL);
}

void
mdone(void)				/* end of graphics metafile */
{
    pflush();
    pglob(PEOF, 0200, NULL);
}

void
minclude(			/* include a file */
char  *fname
)
{
    pflush();
    pglob(PINCL, 1, fname);
}

void
msetpat(		/* set a pattern */
int  pn,
char  *pat
)
{
    pflush();
    pglob(PSET, pn+4, pat);
}


void
mopenseg(			/* open a segment */
char  *sname
)
{
    pflush();
    pglob(POPEN, 0, sname);
}


void
mcloseseg(void)			/* close current segment */
{
    pflush();
    pglob(PCLOSE, 0200, NULL);
}


void
mline(		/* start a line */
int  x, int y,
int  type, int thick, int color
)
{
    pflush();
    cura0 = (type<<4 & 060) | (thick<<2 & 014) | (color & 03);
    curx = x;
    cury = y;
}


void
mrectangle(	/* fill a rectangle */
int  xmin, int ymin, int xmax, int ymax,
int  pat, int color
)
{
    pflush();
    cura0 = (pat<<2 & 014) | (color & 03);
    pprim(PRFILL, cura0, xmin, ymin, xmax, ymax, NULL);
}


void
mtriangle(	/* fill a triangle */
int  xmin, int ymin, int xmax, int ymax,
int  d, int pat, int color
)
{
    pflush();
    cura0 = (let_dir(d)<<4 & 060) | (pat<<2 & 014) | (color & 03);
    pprim(PTFILL, cura0, xmin, ymin, xmax, ymax, NULL);
}


void
mpoly(		/* start a polygon */
int  x, int y,
int  border, int pat, int color
)
{
    pflush();
    cura0 = (border<<6 & 0100) | (pat<<2 & 014) | (color & 03);
    cap = curargs;
    inpoly = TRUE;
    putdec(x);
    *cap++ = ' ';
    putdec(y);
}


void
mtext(		/* matrix string */
int  x, int y,
char  *s,
int  cpi,
int  color
)
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


void
mvstr(	/* vector string */
int  xmin, int ymin, int xmax, int ymax,
char  *s,
int  d, int thick, int color
)
{
    pflush();
    cura0 = (let_dir(d)<<4 & 060) | (thick<<2 & 014) | (color & 03);
    pprim(PVSTR, cura0, xmin, ymin, xmax, ymax, s);
}


void
msegment(	/* segment */
int  xmin, int ymin, int xmax, int ymax,
char  *sname,
int  d, int thick, int color
)
{
    pflush();
    cura0 = (let_dir(d)<<4 & 060) | (thick<<2 & 014) | (color & 03);
    pprim(PSEG, cura0, xmin, ymin, xmax, ymax, sname);
}


void
mdraw(				/* draw to next point */
int  x, int y
)
{
    if (inpoly) {
	polyval(x, y);
    } else if (x != curx || y != cury) {
	plseg(cura0, curx, cury, x, y);
	curx = x;
	cury = y;
    }
}


static void
decival(				/* add value to polygon */
register int  v
)
{
    if (!v)
	return;
    decival(v/10);
    *cap++ = v%10 + '0';
}


static void
polyval(				/* add vertex to a polygon */
register int  x,
register int y
)
{
    *cap++ = ' ';
    putdec(x);
    *cap++ = ' ';
    putdec(y);
}


static void
closepoly(void)				/* close current polygon */
{
    *cap = '\0';
    pprim(PPFILL, cura0, 0, 0, XYSIZE-1, XYSIZE-1, curargs);
    inpoly = FALSE;
}


static int
let_dir(		/* convert letter to corresponding direction */
register int  c
)
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
