/* RCSid $Id: rtmisc.h,v 3.5 2003/07/21 22:30:17 schorsch Exp $ */
/*
 *	Miscellaneous Radiance definitions
 */
#ifndef _RAD_RTMISC_H_
#define _RAD_RTMISC_H_

#include  <stdlib.h>

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

