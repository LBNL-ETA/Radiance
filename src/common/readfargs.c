#ifndef lint
static const char	RCSid[] = "$Id: readfargs.c,v 2.9 2003/06/26 00:58:09 schorsch Exp $";
#endif
/*
 * Allocate, read and free object arguments
 */

#include "copyright.h"

#include "standard.h"

#include "object.h"



int
readfargs(fa, fp)		/* read function arguments from stream */
register FUNARGS  *fa;
FILE  *fp;
{
#define getstr(s)	(fgetword(s,sizeof(s),fp)!=NULL)
#define getint(s)	(getstr(s) && isint(s))
#define getflt(s)	(getstr(s) && isflt(s))
	char  sbuf[MAXSTR];
	register int  n, i;

	if (!getint(sbuf) || (n = atoi(sbuf)) < 0)
		return(0);
	if (fa->nsargs = n) {
		fa->sarg = (char **)malloc(n*sizeof(char *));
		if (fa->sarg == NULL)
			return(-1);
		for (i = 0; i < fa->nsargs; i++) {
			if (!getstr(sbuf))
				return(0);
			fa->sarg[i] = savestr(sbuf);
		}
	} else
		fa->sarg = NULL;
	if (!getint(sbuf) || (n = atoi(sbuf)) < 0)
		return(0);
#ifdef  IARGS
	if (fa->niargs = n) {
		fa->iarg = (long *)malloc(n*sizeof(long));
		if (fa->iarg == NULL)
			return(-1);
		for (i = 0; i < n; i++) {
			if (!getint(sbuf))
				return(0);
			fa->iarg[i] = atol(sbuf);
		}
	} else
		fa->iarg = NULL;
#else
	if (n != 0)
		return(0);
#endif
	if (!getint(sbuf) || (n = atoi(sbuf)) < 0)
		return(0);
	if (fa->nfargs = n) {
		fa->farg = (RREAL *)malloc(n*sizeof(RREAL));
		if (fa->farg == NULL)
			return(-1);
		for (i = 0; i < n; i++) {
			if (!getflt(sbuf))
				return(0);
			fa->farg[i] = atof(sbuf);
		}
	} else
		fa->farg = NULL;
	return(1);
#undef getflt
#undef getint
#undef getstr
}


void
freefargs(fa)			/* free object arguments */
register FUNARGS  *fa;
{
	register int  i;

	if (fa->nsargs) {
		for (i = 0; i < fa->nsargs; i++)
			freestr(fa->sarg[i]);
		free((void *)fa->sarg);
		fa->sarg = NULL;
		fa->nsargs = 0;
	}
#ifdef  IARGS
	if (fa->niargs) {
		free((void *)fa->iarg);
		fa->iarg = NULL;
		fa->niargs = 0;
	}
#endif
	if (fa->nfargs) {
		free((void *)fa->farg);
		fa->farg = NULL;
		fa->nfargs = 0;
	}
}
