/* RCSid: $Id: mgvars.h,v 1.4 2003/11/15 02:13:37 schorsch Exp $ */
/*
 *  mgvars.h - header file for graphing routines using variables.
 *
 *     6/23/86
 *
 *     Greg Ward Larson
 */
#ifndef _RAD_MGVARS_H_
#define _RAD_MGVARS_H_

#include <errno.h>

#include "calcomp.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 *  Data arrays are used to store point data.
 */

typedef struct {
	int  size;		/* the size of the array */
	float  *data;		/* pointer to the array */
}  DARRAY;

/*
 *  Intermediate variables are not referenced
 *  directly by the program, but may be used
 *  for defining other variables.
 */

typedef struct ivar {
	char  *name;		/* the variable name */
	char  *dfn;		/* its definition */
	struct ivar  *next;	/* next in list */
}  IVAR;		/* an intermediate variable */

/*
 *  Variables are used for all graph parameters.  The four variable
 *  types are:  REAL, FUNCTION, STRING and DATA.
 *  Of these, only STRING and DATA must be interpreted by us; an
 *  expression parser will handle the rest.
 */
				/* standard variables */
#define  FTHICK		0		/* frame thickness */
#define  GRID		1		/* grid on? */
#define  LEGEND		2		/* legend title */
#define  OTHICK		3		/* origin thickness */
#define  PERIOD		4		/* period of polar plot */
#define  SUBTITLE	5		/* subtitle */
#define  SYMFILE	6		/* symbol file */
#define  TSTYLE		7		/* tick mark style */
#define  TITLE		8		/* title */
#define  XLABEL		9		/* x axis label */
#define  XMAP		10		/* x axis mapping */
#define  XMAX		11		/* x axis maximum */
#define  XMIN		12		/* x axis minimum */
#define  XSTEP		13		/* x axis step */
#define  YLABEL		14		/* y axis label */
#define  YMAP		15		/* y axis mapping */
#define  YMAX		16		/* y axis maximum */
#define  YMIN		17		/* y axis minimum */
#define  YSTEP		18		/* y axis step */

#define  NVARS		19	/* number of standard variables */

				/* curve variables */
#define  C		0		/* the curve function */
#define  CCOLOR		1		/* the curve color */
#define  CDATA		2		/* the curve data */
#define  CLABEL		3		/* the curve label */
#define  CLINTYPE	4		/* the curve line type */
#define  CNPOINTS	5		/* number of curve points to plot */
#define  CSYMSIZE	6		/* the curve symbol size */
#define  CSYMTYPE	7		/* the curve symbol type */
#define  CTHICK		8		/* the curve line thickness */

#define  NCVARS		9	/* number of curve variables */

#define  MAXCUR		8	/* maximum number of curves */

				/* types */
#define  REAL		1		/* floating point */
#define  STRING		2		/* character array */
#define  FUNCTION	3		/* function definition */
#define  DATA		4		/* a set of real points */

				/* flags */
#define  DEFINED	01		/*  variable is defined */

typedef struct {
	char  *name;		/* the variable name */
	short  type;		/* the variable type */
	char  *descrip;		/* brief description */
	unsigned short  flags;	/* DEFINED, etc. */
	union {
		char  *s;		/* STRING value */
		DARRAY  d;		/* DATA value */
		char  *dfn;		/* variable definition */
	} v;			/* value */
}  VARIABLE;		/* a variable */

extern IVAR  *ivhead;				/* the intermediate list */

extern VARIABLE  gparam[NVARS];			/* the graph variables */
extern VARIABLE  cparam[MAXCUR][NCVARS];	/* the curve variables */

extern VARIABLE  *vlookup();

#define  mgclear(vname)		undefine(vlookup(vname))

extern void mgclearall(void);
extern void mgload(char *file);
extern void mgsave(char *file);
extern void setmgvar(char *fname, FILE *fp, char *string);
extern int mgcurve(int c, void (*f)(int c, double x, double y));
extern void mgtoa(register char *s, VARIABLE *vp);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_MGVARS_H_ */

