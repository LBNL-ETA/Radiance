/* RCSid $Id: rtmisc.h,v 3.1 2003/06/27 06:53:21 greg Exp $ */
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
#define	 copystruct(d,s)	bcopy((void *)(s),(void *)(d),sizeof(*(d)))
#else
#define	 copystruct(d,s)	(*(d) = *(s))
#endif

#ifndef BSD
#define	 bcopy(s,d,n)		(void)memcpy(d,s,n)
#define	 bzero(d,n)		(void)memset(d,0,n)
#define	 bcmp(b1,b2,n)		memcmp(b1,b2,n)
#endif

#ifdef _WIN32
#define NIX 1
#endif
#ifdef AMIGA
#define NIX 1
#endif
					/* defined in bmalloc.c */
extern char	*bmalloc(unsigned int n);
extern void	bfree(char *p, unsigned int n);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RTMISC_H_ */

