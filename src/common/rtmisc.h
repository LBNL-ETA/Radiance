/* RCSid $Id: rtmisc.h,v 3.6 2004/03/28 20:33:12 schorsch Exp $ */
/*
 *	Miscellaneous Radiance definitions
 */
#ifndef _RAD_RTMISC_H_
#define _RAD_RTMISC_H_

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

					/* defined in bmalloc.c */
extern void	*bmalloc(size_t size);
extern void	bfree(void *ptr, size_t size);

					/* defined in ealloc.c */
extern void *emalloc(size_t size);
extern void *ecalloc(size_t nmemb, size_t size);
extern void *erealloc(void *ptr, size_t size);
extern void	efree(void *ptr);

					/* defined in myhostname.c */
extern char	*myhostname(void);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RTMISC_H_ */

