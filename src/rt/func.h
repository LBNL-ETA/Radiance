/* RCSid $Id: func.h,v 2.5 2003/06/27 06:53:22 greg Exp $ */
/*
 * Header file for modifiers using function files.
 *
 * Include after ray.h
 */
#ifndef _RAD_FUNC_H_
#define _RAD_FUNC_H_
#ifdef __cplusplus
extern "C" {
#endif

#include  "calcomp.h"

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


#ifdef __cplusplus
}
#endif
#endif /* _RAD_FUNC_H_ */

