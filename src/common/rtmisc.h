/* RCSid $Id: rtmisc.h,v 3.3 2003/07/14 22:23:59 schorsch Exp $ */
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

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RTMISC_H_ */

