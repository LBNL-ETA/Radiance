/* RCSid $Id: rtmisc.h,v 3.2 2003/06/30 14:59:11 schorsch Exp $ */
/*
 *	Miscellaneous Radiance definitions
 */
#ifndef _RAD_RTMISC_H_
#define _RAD_RTMISC_H_
#ifdef __cplusplus
extern "C" {
#endif

#include  <stdlib.h>
					/* memory operations */
#ifdef	NOSTRUCTASS
#include  <string.h>
#define	 copystruct(d,s)	memcpy((void *)(d),(void *)(s),sizeof(*(d)))
#else
#define	 copystruct(d,s)	(*(d) = *(s))
#endif

					/* defined in bmalloc.c */
extern char	*bmalloc(unsigned int n);
extern void	bfree(char *p, unsigned int n);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RTMISC_H_ */

