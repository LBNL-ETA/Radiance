/* RCSid: $Id: meta.h,v 1.1 2003/02/22 02:07:26 greg Exp $ */
/*
 *   Standard meta-file definitions and limits
 */


#include  <stdio.h>

#include  <ctype.h>


#define  TRUE  1
#define  FALSE  0

#define  PEOF  'F'		/* end of file global */
#define  PEOP  'E'		/* end of page global */
#define  PPAUSE  'P'		/* pause global */
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

#define  SYSTEM	 0		/* system error, internal, fatal */
#define  USER  1		/* user error, fatal */
#define  WARNING  2		/* user error, not fatal */

#ifdef  UNIX
#define  TDIR  "/tmp/"		/* directory for temporary files */
#ifndef  MDIR
#define  MDIR  "/usr/local/lib/meta/"	/* directory for metafiles */
#endif
#define  TTY  "/dev/tty"	/* console name */
#endif

#ifdef  MAC
#define  TDIR  "5:tmp/"		/* directory for temporary files */
#define  MDIR  "/meta/"		/* directory for metafiles */
#define  TTY  ".con"		/* console name */
#endif

#ifdef  CPM
#define  TDIR  "0/"		/* directory for temporary files */
#define  MDIR  "0/"		/* directory for metafiles */
#define  TTY  "CON:"		/* console name */
#endif

#define  MAXFNAME  64		/* maximum file name length */

#define  XYSIZE  (1<<14)	/* metafile coordinate size */

#define  max(x, y)  ((x) > (y) ? (x) : (y))
#define  min(x, y)  ((x) < (y) ? (x) : (y))

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

char  *malloc(), *savestr();

PRIMITIVE  *palloc(), *pop();

FILE  *efopen(), *mfopen();

extern char  coms[];

extern char  errmsg[];

extern char  *progname;
