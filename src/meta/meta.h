/* RCSid: $Id: meta.h,v 1.7 2003/08/01 14:14:24 schorsch Exp $ */
/*
 *   Standard meta-file definitions and limits
 */
#ifndef _RAD_META_H_
#define _RAD_META_H_

#include "copyright.h"

#include  <stdio.h>
#include  <stdlib.h>
#include  <ctype.h>

#include  "rterror.h"

#ifdef __cplusplus
extern "C" {
#endif

#define  TRUE  1
#define  FALSE  0

#define  PEOF  'F'		/* end of file global */
#define  PEOP  'E'		/* end of page global */
#define  PPAUS  'P'		/* pause global */
#define  PDRAW  'D'		/* draw global */
#define  POPEN  'O'		/* open segment */
#define  PCLOSE  'C'		/* close segment */
#define  PSET  'S'		/* set global */
#define  PUNSET  'U'		/* unset global */
#define  PRESET  'R'		/* reset global to default */
#define  PINCL  'I'		/* include file */

#define  PLSEG  'l'		/* line segment command */
#define  PRFILL  'r'		/* rectangle fill command */
#define  PTFILL  't'		/* triangle fill command */
#define  PMSTR  'm'		/* matrix string command */
#define  PVSTR  'v'		/* vector string command */
#define  PSEG  's'		/* print segment command */
#define  PPFILL  'p'		/* polygon fill command */

#define  NCOMMANDS  17		/* number of commands */

#define  COML  "lrtmsvpOCESURPDIF"	/* command letters */

#define  ADELIM  '`'		/* additional argument delimiter */
#define  CDELIM  '#'		/* comment delimiter */

#define  MAXARGS  512		/* maximum # of arguments for primitive */

#define  SALL  0		/* set all */
#define  SPAT0  04		/* set pattern 0 */
#define  SPAT1  05		/* set pattern 1 */
#define  SPAT2  06		/* set pattern 2 */
#define  SPAT3  07		/* set pattern 3 */


#ifdef  _WIN32  /* XXX */
#define MDIR "c\\tmp\\" /* XXX we just need something to compile for now */
#define TTY "CON:"   /* XXX this probably doesn't work */
#else  /* XXX */

#define  TDIR  "/tmp/"		/* directory for temporary files */
#ifndef  MDIR
#define  MDIR  "/usr/local/lib/meta/"	/* directory for metafiles */
#endif
#define  TTY  "/dev/tty"	/* console name */
#endif

#define  MAXFNAME  64		/* maximum file name length */

#define  XYSIZE  (1<<14)	/* metafile coordinate size */

#ifndef max
#define  max(x, y)  ((x) > (y) ? (x) : (y))
#endif
#ifndef min
#define  min(x, y)  ((x) < (y) ? (x) : (y))
#endif

#define  abs(x)  ((x) < 0 ? -(x) : (x))

#define  iscom(c)  (comndx(c) != -1)
#define  isglob(c)  isupper(c)
#define  isprim(c)  islower(c)

#define  WIDTH(wspec)  ((wspec)==0 ? 0 : 12*(1<<(wspec)))
#define  CONV(coord, size)  ((int)(((long)(coord)*(size))>>14))
#define  ICONV(dcoord, size)  ((int)(((long)(dcoord)<<14)/(size)))

#define  XMN  0			/* index in xy array for xmin */
#define  YMN  1			/* index in xy array for ymin */
#define  XMX  2			/* index in xy array for xmax */
#define  YMX  3			/* index in xy array for ymax */


/*
 *    Structure definitions for primitives
 */

struct primitive  {			/* output primitive */
                   short  com,		/* command (0 - 127) */
			  arg0;		/* first argument (1 byte) */
		   int  xy[4];		/* extent=(xmin,ymin,xmax,ymax) */
		   char  *args;		/* additional arguments */
		   struct primitive  *pnext;	/* next primitive */
		   };

typedef struct primitive  PRIMITIVE;

struct plist  {				/* list of primitives */
	       PRIMITIVE  *ptop, *pbot;
	       };

typedef struct plist  PLIST;


/*
 *   External declarations
 */

char  *savestr();

PRIMITIVE *pop();

FILE  *efopen(), *mfopen();

extern char  coms[];
extern char  errmsg[];
extern char  *progname;

	/* expand.c */
extern void expand(FILE *infp, short *exlist);
	/* palloc.c */
extern PRIMITIVE *palloc(void);
extern void pfree(register PRIMITIVE *p);
extern void plfree(register PLIST *pl);
	/* sort.c */
extern void sort(FILE *infp, int (*pcmp)());
extern void pmergesort(FILE *fi[], int nf, PLIST *pl, int (*pcmp)(), FILE *ofp);
	/* metacalls.c */
extern void mdraw(int x, int y);
extern void msegment(int xmin, int ymin, int xmax, int ymax, char *sname,
		int d, int thick, int color);
extern void mvstr(int xmin, int ymin, int xmax, int ymax, char *s,
		int d, int thick, int color);
extern void mtext(int x, int y, char *s, int cpi, int color);
extern void mpoly(int x, int y, int border, int pat, int color);
extern void mtriangle(int xmin, int ymin, int xmax, int ymax,
		int d, int pat, int color);
extern void mrectangle(int xmin, int ymin, int xmax, int ymax,
		int pat, int color);
extern void mline(int x, int y, int type, int thick, int color);
extern void mcloseseg(void);
extern void mopenseg(char *sname);
extern void msetpat(int pn, char *pat);
extern void minclude(char *fname);
extern void mdone(void);
extern void mendpage(void);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_META_H_ */

