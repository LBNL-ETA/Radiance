#ifndef lint
static const char	RCSid[] = "$Id: portio.c,v 2.13 2005/05/03 05:10:13 greg Exp $";
#endif
/*
 * Portable i/o for binary files
 *
 * External symbols declared in standard.h
 */

#include "copyright.h"

#include <stdio.h>

#include "rtio.h"

#include <math.h>

#ifdef getc_unlocked		/* avoid horrendous overhead of flockfile */
#undef getc
#undef putc
#define getc    getc_unlocked
#define putc    putc_unlocked
#endif


void
putstr(s, fp)			/* write null-terminated string to fp */
register char  *s;
register FILE  *fp;
{
	do
		putc(*s, fp);
	while (*s++);
}


void
putint(i, siz, fp)		/* write a siz-byte integer to fp */
long  i;
register int  siz;
register FILE  *fp;
{
	while (siz--)
		putc((int)(i>>(siz<<3) & 0xff), fp);
}


void
putflt(f, fp)			/* put out floating point number */
double	f;
FILE  *fp;
{
	long  m;
	int  e;

	m = frexp(f, &e) * 0x7fffffff;
	if (e > 127) {			/* overflow */
		m = m > 0 ? (long)0x7fffffff : -(long)0x7fffffff;
		e = 127;
	} else if (e < -128) {		/* underflow */
		m = 0;
		e = 0;
	}
	putint(m, 4, fp);
	putint((long)e, 1, fp);
}


char *
getstr(s, fp)			/* get null-terminated string */
char  *s;
register FILE  *fp;
{
	register char  *cp;
	register int  c;

	cp = s;
	while ((c = getc(fp)) != EOF)
		if ((*cp++ = c) == '\0')
			return(s);
		
	return(NULL);
}


long
getint(siz, fp)			/* get a siz-byte integer */
int  siz;
register FILE  *fp;
{
	register int  c;
	register long  r;

	if ((c = getc(fp)) == EOF)
		return(EOF);
	r = 0x80&c ? -1<<8|c : c;		/* sign extend */
	while (--siz > 0) {
		if ((c = getc(fp)) == EOF)
			return(EOF);
		r <<= 8;
		r |= c;
	}
	return(r);
}


double
getflt(fp)			/* get a floating point number */
FILE  *fp;
{
	long	l;
	double	d;

	l = getint(4, fp);
	if (l == 0) {
		getc(fp);		/* exactly zero -- ignore exponent */
		return(0.0);
	}
	d = (l + (l > 0 ? .5 : -.5)) * (1./0x7fffffff);
	return(ldexp(d, (int)getint(1, fp)));
}
