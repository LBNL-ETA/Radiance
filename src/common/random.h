/* RCSid $Id$ */
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
				/* defined in random.c */
extern void	srandom(unsigned long s);
extern long	random(void);
extern double	frandom(void);

#else

#define	 random()	lrand48()
#define  srandom(s)	srand48((long)(s))
#define	 frandom()	drand48()

#endif

extern unsigned short	*urperm;
extern int	urmask;

#define	 urand(i)	(urmask ? (urperm[(i)&urmask]+frandom())/(urmask+1.) \
				: frandom())

				/* defined in urand.c */
extern long	irandom(long modulus);

extern int	initurand(int size);

extern int	ilhash(int *d, int n);

extern int	urind(int s, int i);
				/* defined in multisamp.c */
extern void	multisamp(double t[], int n, double r);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_RANDOM_H_ */
