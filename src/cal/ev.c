#ifndef lint
static const char	RCSid[] = "$Id: ev.c,v 1.6 2003/11/14 17:31:24 schorsch Exp $";
#endif
/*
 *  ev.c - program to evaluate expression arguments
 *
 *     1/29/87
 */

#include  <stdlib.h>
#include  <stdio.h>
#include  <errno.h>

#include  "calcomp.h"
#include  "rterror.h"


int
main(argc, argv)
int  argc;
char  *argv[];
{
	int  i;

	esupport |= E_FUNCTION;
	esupport &= ~(E_VARIABLE|E_INCHAN|E_OUTCHAN|E_RCONST);

#ifdef  BIGGERLIB
	biggerlib();
#endif

	errno = 0;
	for (i = 1; i < argc; i++)
		printf("%.9g\n", eval(argv[i]));

	quit(errno ? 2 : 0);
	return (errno ? 2 : 0); /* pro forma return */
}


void
eputs(msg)
char  *msg;
{
	fputs(msg, stderr);
}


void
wputs(msg)
char  *msg;
{
	eputs(msg);
}


void
quit(code)
int  code;
{
	exit(code);
}
