/* RCSid $Id: rtmisc.h,v 3.4 2003/07/17 09:21:29 schorsch Exp $ */
/*
 *	Miscellaneous Radiance definitions
 */
#ifndef _RAD_RTMISC_H_
#define _RAD_RTMISC_H_

#include  <stdlib.h>
					/* memory operations */
#ifdef	NOSTRUCTASS
#include  <string.h>
#define	 copystruct(d,s)	memcpy((void *)(d),(void *)(s),sizeof(*(d)))
#else
#define	 copystruct(d,s)	(*(d) = *(s))
#endif

#ifdef __cplusplus
extern "C" {
#endif

					/* defined in bmalloc.c */
extern char	*bmalloc(unsigned int n);
extern void	bfree(char *p, unsigned int n);
					/*  defined in myhostname.c */
extern char	*myhostname(void);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RTMISC_H_ */

