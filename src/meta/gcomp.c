#ifndef lint
static const char	RCSid[] = "$Id: gcomp.c,v 1.3 2003/11/15 02:13:37 schorsch Exp $";
#endif
/*
 *  gcomp.c - program to calculate things from graph files.
 *
 *     7/7/86
 *
 *     Greg Ward Larson
 */

#include  <stdlib.h>
#include  <stdio.h>

#include "rterror.h"
#include "meta.h"
#include "mgvars.h"

#define  istyp(s)	(s[0] == '-')

#define  isvar(s)	(s[0] == '+')

char  *progname;

char  *libpath[4];

static void dofile(int  optc, char  *optv[], char  *file);

int
main(
	int  argc,
	char  *argv[]
)
{
	char  *getenv();
	int  i, file0;

	progname = argv[0];
	libpath[0] = "./";
	if ((libpath[i=1] = getenv("MDIR")) != NULL)
		i++;
	libpath[i++] = MDIR;
	libpath[i] = NULL;

	for (file0 = 1; file0 < argc; )
		if (istyp(argv[file0]))
			file0++;
		else if (isvar(argv[file0]) && file0 < argc-1)
			file0 += 2;
		else
			break;

	if (file0 >= argc)
		dofile(argc-1, argv+1, NULL);
	else
		for (i = file0; i < argc; i++)
			dofile(file0-1, argv+1, argv[i]);

	quit(0);
	return 0; /* pro forma return */
}


void
dofile(		/* plot a file */
	int  optc,
	char  *optv[],
	char  *file
)
{
	char  types[16], stmp[256], *strcat();
	int  i;

	mgclearall();

	mgload(file);

	types[0] = '\0';
	for (i = 0; i < optc; i++)
		if (istyp(optv[i]))
			strcat(types, optv[i]+1);
		else {
			sprintf(stmp, "%s=%s", optv[i]+1, optv[i+1]);
			setmgvar("command line", stdin, stmp);
			i++;
		}

	gcalc(types);
}


void
eputs(				/* print error message */
	char  *msg
)
{
	fputs(msg, stderr);
}


void
quit(				/* quit program */
	int  code
)
{
	exit(code);
}
