/* RCSid $Id: random.h,v 2.23 2016/03/06 01:13:17 schorsch Exp $ */
/*
 *  random.h - header file for random(3) and urand() function.
 */
#ifndef _RAD_RANDOM_H_
#define _RAD_RANDOM_H_

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) || defined(_WIN64)
#if (RAND_MAX <= 65536)
#define random()	((long)rand()<<16^(long)rand()<<6^(long)rand()>>4)
#else
#define random()	rand()
#endif
#define srandom(s)	srand((unsigned)(s))

#define frandom()	(rand()*(1./(RAND_MAX+.5)))

#else

#define	 random()	lrand48()
#define  srandom(s)	srand48((long)(s))
#define	 frandom()	drand48()

#endif

extern unsigned short	*urperm;
extern int	urmask;

#define	 urand(i)	(urmask ? (urperm[(i)&urmask]+frandom())/(urmask+1.) \
				: frandom())

extern int	initurand(int size);

				/* defined in urand.c */
extern int	ilhash(int *d, int n);
				/* defined in urind.c */
extern int	urind(int s, int i);
				/* defined in multisamp.c */
extern void	multisamp(double t[], int n, double r);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_RANDOM_H_ */

