#ifndef lint
static const char	RCSid[] = "$Id: gcomp.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  gcomp.c - program to calculate things from graph files.
 *
 *     7/7/86
 *
 *     Greg Ward Larson
 */

#include  <stdio.h>

#define  istyp(s)	(s[0] == '-')

#define  isvar(s)	(s[0] == '+')

char  *progname;

char  *libpath[4];


main(argc, argv)
int  argc;
char  *argv[];
{
#if  UNIX || MAC
	char  *getenv();
#endif
	int  i, file0;

#ifdef  CPM
	fixargs("gcomp", &argc, &argv);
	progname = argv[0];
	libpath[0] = "";
	libpath[1] = "0/";
	libpath[2] = NULL;
#endif
#ifdef  MAC
	progname = argv[0];
	libpath[0] = "./";
	if ((libpath[i=1] = getenv("MDIR")) != NULL)
		i++;
	libpath[i++] = "/meta/";
	libpath[i] = NULL;
#endif
#ifdef  UNIX
	progname = argv[0];
	libpath[0] = "./";
	if ((libpath[i=1] = getenv("MDIR")) != NULL)
		i++;
	libpath[i++] = MDIR;
	libpath[i] = NULL;
#endif

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
}


dofile(optc, optv, file)		/* plot a file */
int  optc;
char  *optv[];
char  *file;
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


eputs(msg)				/* print error message */
char  *msg;
{
	fputs(msg, stderr);
}


quit(code)				/* quit program */
int  code;
{
	exit(code);
}
