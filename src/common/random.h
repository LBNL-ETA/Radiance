/* RCSid $Id: random.h,v 2.15 2005/06/13 20:07:56 greg Exp $ */
/*
 *  random.h - header file for random(3) and urand() function.
 */
#ifndef _RAD_RANDOM_H_
#define _RAD_RANDOM_H_
#ifdef __cplusplus
extern "C" {
#endif

#ifdef	_WIN32

#if (RAND_MAX <= 65536)
#define random()	((long)rand()<<16^(long)rand()<<6^(long)rand()>>4)
#else
#define random()	rand()
#endif
#define srandom(s)	srand((unsigned)(s))

#define frandom()	(rand()*(1./RAND_MAX))

#else
#ifdef	BSD

extern long  random();

#define	 frandom()	(random()*(1./2147483648.))

#else

extern long  lrand48();
extern double  drand48();

#define	 random()	lrand48()
#define  srandom(s)	srand48((long)(s))
#define	 frandom()	drand48()

#endif
#endif

extern unsigned short	*urperm;
extern int	urmask;

#define	 urand(i)	((urperm[(i)&urmask]+frandom())/(urmask+1))

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

