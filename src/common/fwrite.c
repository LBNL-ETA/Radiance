#ifndef lint
static const char	RCSid[] = "$Id: fwrite.c,v 2.3 2003/02/25 02:47:21 greg Exp $";
#endif
#include "copyright.h"

#include  <stdio.h>


int
fwrite(ptr, size, count, fp)		/* system version is broken (*sigh*) */
register char  *ptr;
int  size, count;
register FILE  *fp;
{
    register long  n = size*count;

    while (n--)
	putc(*ptr++, fp); 

    return(ferror(fp) ? 0 : count);
}
