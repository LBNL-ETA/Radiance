/* RCSid $Id: func.h,v 2.8 2011/02/18 00:40:25 greg Exp $ */
/*
 * Header file for modifiers using function files.
 *
 * Include after ray.h
 */
#ifndef _RAD_FUNC_H_
#define _RAD_FUNC_H_

#include "calcomp.h"

#ifdef __cplusplus
extern "C" {
#endif

#define  MAXEXPR	9	/* maximum expressions in modifier */

typedef struct {
	EPNODE  *ep[MAXEXPR+1];		/* NULL-terminated expression list */
	char  *ctx;			/* context (from file name) */
	XF  *f, *b;			/* forward and backward transforms */
} MFUNC;			/* material function */

extern XF  unitxf;		/* identity transform */
extern XF  funcxf;		/* current transform */


extern MFUNC	*getfunc(OBJREC *m, int ff, unsigned int ef, int dofwd);
extern void	freefunc(OBJREC *m);
extern int	setfunc(OBJREC *m, RAY *r);
extern void	loadfunc(char *fname);

	/* defined in noise3.c */
extern void	setnoisefuncs(void);

	/* defined in fprism.c */
extern void	setprismfuncs(void);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_FUNC_H_ */

