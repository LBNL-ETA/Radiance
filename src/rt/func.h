/* RCSid $Id: func.h,v 2.3 2003/02/25 02:47:22 greg Exp $ */
/*
 * Header file for modifiers using function files.
 *
 * Include after ray.h
 */

#include "copyright.h"

#include  "calcomp.h"

#define  MAXEXPR	9	/* maximum expressions in modifier */

typedef struct {
	EPNODE  *ep[MAXEXPR+1];		/* NULL-terminated expression list */
	char  *ctx;			/* context (from file name) */
	XF  *f, *b;			/* forward and backward transforms */
} MFUNC;			/* material function */

extern XF  unitxf;		/* identity transform */
extern XF  funcxf;		/* current transform */

#ifdef NOPROTO

extern MFUNC	*getfunc();
extern void	freefunc();
extern int	setfunc();
extern void	loadfunc();

#else

extern MFUNC	*getfunc(OBJREC *m, int ff, unsigned int ef, int dofwd);
extern void	freefunc(OBJREC *m);
extern int	setfunc(OBJREC *m, RAY *r);
extern void	loadfunc(char *fname);

#endif
