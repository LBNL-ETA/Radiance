#ifndef lint
static const char	RCSid[] = "$Id: portio.c,v 2.21 2019/07/05 03:04:22 greg Exp $";
#endif
/*
 * Portable i/o for binary files
 *
 * External symbols declared in rtio.h
 */

#include "copyright.h"

#include "rtio.h"

#include <math.h>


void
putstr(				/* write null-terminated string to fp */
	char  *s,
	FILE  *fp
)
{
	do
		putc(*s, fp);
	while (*s++);
}


void
putint(				/* write a siz-byte integer to fp */
	long  i,
	int  siz,
	FILE  *fp
)
{
	siz <<= 3;
	while ((siz -= 8) >= 0)
		putc((int)(i>>siz & 0xff), fp);
}


void
putflt(				/* put out floating point number */
	double	f,
	FILE  *fp
)
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


int
putbinary(			/* fwrite() replacement for small objects */
	const void *p,
	int elsiz,
	int nel,
	FILE *fp)
{
	const char	*s = (const char *)p;
	int		nbytes = elsiz*nel;

	if (nbytes > 256)
		return(fwrite(p, elsiz, nel, fp));
	
	while (nbytes-- > 0)
		if (putc(*s++, fp) == EOF)
			return((elsiz*nel - nbytes)/elsiz);

	return(nel);
}


char *
getstr(				/* get null-terminated string */
	char  *s,
	FILE  *fp
)
{
	char  *cp;
	int  c;

	cp = s;
	while ((c = getc(fp)) != EOF)
		if ((*cp++ = c) == '\0')
			return(s);
		
	return(NULL);
}


long
getint(				/* get a siz-byte integer */
	int  siz,
	FILE  *fp
)
{
	int  c;
	long  r;

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
getflt(				/* get a floating point number */
	FILE  *fp
)
{
	long	l;
	double	d;

	l = getint(4, fp);
	if (l == EOF && feof(fp))	/* EOF? */
		return((double)EOF);
	if (l == 0) {
		getc(fp);		/* exactly zero -- ignore exponent */
		return(0.0);
	}
	d = (l + (l > 0 ? .5 : -.5)) * (1./0x7fffffff);
	return(ldexp(d, (int)getint(1, fp)));
}


int
getbinary(			/* fread() replacement for small objects */
	void *p,
	int elsiz,
	int nel,
	FILE *fp)
{
	char	*s = (char *)p;
	int	nbytes = elsiz*nel;
	int	c;

	if (nbytes > 256)
		return(fread(p, elsiz, nel, fp));
	
	while (nbytes-- > 0) {
		if ((c = getc(fp)) == EOF)
			return((elsiz*nel - nbytes)/elsiz);
		*s++ = c;
	}
	return(nel);
}
