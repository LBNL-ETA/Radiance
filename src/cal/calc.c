#ifndef lint
static const char	RCSid[] = "$Id: calc.c,v 1.5 2003/11/14 17:31:24 schorsch Exp $";
#endif
/*
 *  calc.c - simple algebraic desk calculator program.
 *
 *     4/1/86
 */

#include  <stdlib.h>
#include  <stdio.h>
#include  <string.h>
#include  <setjmp.h>
#include  <ctype.h>

#include  "rterror.h"
#include  "calcomp.h"

#define  MAXRES		100

double  result[MAXRES];
int	nres = 0;

jmp_buf  env;
int  recover = 0;


int
main(argc, argv)
int  argc;
char  *argv[];
{
	char  expr[512];
	FILE  *fp;
	int  i;
	register char  *cp;

	esupport |= E_VARIABLE|E_INCHAN|E_FUNCTION;
	esupport &= ~(E_REDEFW|E_RCONST|E_OUTCHAN);
#ifdef  BIGGERLIB
	biggerlib();
#endif
	varset("PI", ':', 3.14159265358979323846);

	for (i = 1; i < argc; i++)
		fcompile(argv[i]);

	setjmp(env);
	recover = 1;
	eclock++;

	while (fgets(expr, sizeof(expr), stdin) != NULL) {
		for (cp = expr; *cp && *cp != '\n'; cp++)
			;
		*cp = '\0';
		switch (expr[0]) {
		case '\0':
			continue;
		case '?':
			for (cp = expr+1; isspace(*cp); cp++)
				;
			if (*cp)
				dprint(cp, stdout);
			else
				dprint(NULL, stdout);
			continue;
		case '>':
			for (cp = expr+1; isspace(*cp); cp++)
				;
			if (!*cp) {
				eputs("file name required\n");
				continue;
			}
			if ((fp = fopen(cp, "w")) == NULL) {
				eputs(cp);
				eputs(": cannot open\n");
				continue;
			}
			dprint(NULL, fp);
			fclose(fp);
			continue;
		case '<':
			for (cp = expr+1; isspace(*cp); cp++)
				;
			if (!*cp) {
				eputs("file name required\n");
				continue;
			}
			fcompile(cp);
			eclock++;
			continue;
		}
		if ((cp = strchr(expr, '=')) != NULL ||
				(cp = strchr(expr, ':')) != NULL) {
			if (cp[1])
				scompile(expr, NULL, 0);
			else if (*cp == '=') {
				*cp = '\0';
				if (!strcmp(expr, "*"))
					dcleanup(1);
				else
					dclear(expr);
			} else {
				*cp = '\0';
				if (!strcmp(expr, "*"))
					dcleanup(2);
				else
					dremove(expr);
			}
			eclock++;
		} else {
			printf("$%d=%.9g\n", nres+1,
					result[nres%MAXRES] = eval(expr));
			nres++;
		}
	}

	recover = 0;
	quit(0);
	return 0; /* pro forma exit */
}


double
chanvalue(n)			/* return channel value */
int  n;
{
	if (n == 0)
		n = nres;
	else if (n > nres || nres-n >= MAXRES) {
		fprintf(stderr, "$%d: illegal result\n", n);
		return(0.0);
	}
	return(result[(n-1)%MAXRES]);
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
	if (recover)			/* a cavalier approach */
		longjmp(env, 1);
	exit(code);
}
