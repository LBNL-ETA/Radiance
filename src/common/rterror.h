/* RCSid $Id: rterror.h,v 3.2 2003/07/14 22:23:59 schorsch Exp $ */
/*
 * Header for Radiance error-handling routines
 */

#ifndef _RAD_RTERROR_H_
#define _RAD_RTERROR_H_

#include <errno.h>

#ifdef __cplusplus
extern "C" {
#endif

				/* error codes */
#define	 WARNING	0		/* non-fatal error */
#define	 USER		1		/* fatal user-caused error */
#define	 SYSTEM		2		/* fatal system-related error */
#define	 INTERNAL	3		/* fatal program-related error */
#define	 CONSISTENCY	4		/* bad consistency check, abort */
#define	 COMMAND	5		/* interactive error */
#define  NERRS		6
				/* error struct */
extern struct erract {
	char	pre[16];		/* prefix message */
	void	(*pf)();		/* put function (resettable) */
	int	ec;			/* exit code (0 means non-fatal) */
} erract[NERRS];	/* list of error actions */

#define  ERRACT_INIT	{	{"warning - ", wputs, 0}, \
				{"fatal - ", eputs, 1}, \
				{"system - ", eputs, 2}, \
				{"internal - ", eputs, 3}, \
				{"consistency - ", eputs, -1}, \
				{"", NULL, 0}	}

extern char  errmsg[];			/* global buffer for error messages */

					/* custom version of assert(3) */
#define  CHECK(be,et,em)	if (be) error(et,em); else
#ifdef  DEBUG
#define  DCHECK			CHECK
#else
#define  DCHECK(be,et,em)	(void)0
#endif
					/* defined in error.c */
extern void	error(int etype, char *emsg);
					/* error & warning output & exit */
extern void	eputs(char *s);
extern void	wputs(char *s);
extern void	quit(int code);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RTERROR_H_ */

