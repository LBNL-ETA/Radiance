/* RCSid $Id: random.h,v 2.9 2003/02/25 02:47:21 greg Exp $ */
/*
 *  random.h - header file for random(3) and urand() function.
 */

#include "copyright.h"

#ifdef  NORANDOM

#undef  random
#define  random()	1073741820
#undef  srandom
#define  srandom(s)	(s)
#define  frandom()	0.5
#define  urand(i)	0.5
#define  initurand(n)	(n)

#else

#ifdef	MSDOS

#define random()	((long)rand()<<16^(long)rand()<<6^(long)rand()>>4)
#define srandom(s)	srand((unsigned)(s))

#define frandom()	(rand()*(1./32768.))

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

#ifdef  MC

#define  urand(i)	frandom()
#define  initurand(n)	(n)

#else

extern unsigned short	*urperm;
extern int	urmask;

#define	 urand(i)	((urperm[(i)&urmask]+frandom())/(urmask+1))

#endif

#endif

#ifdef NOPROTO

extern int	initurand();
extern int	ilhash();
extern int	urind();
extern void	multisamp();

#else
				/* defined in urand.c */
extern int	initurand(int size);
extern int	ilhash(int *d, int n);
				/* defined in urind.c */
extern int	urind(int s, int i);
				/* defined in multisamp.c */
extern void	multisamp(double t[], int n, double r);

#endif
