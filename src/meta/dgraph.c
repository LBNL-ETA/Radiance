#ifndef lint
static const char	RCSid[] = "$Id: dgraph.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  dgraph.c - program to send plots to dumb terminal.
 *
 *     7/7/86
 *
 *     Greg Ward Larson
 */

#include  <stdio.h>


#define  isopt(s)	(s[0] == '+' || s[0] == '-')

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
	fixargs("dgraph", &argc, &argv);
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
	libpath[i++] = "/usr/local/lib/meta/";
	libpath[i] = NULL;
#endif

	for (file0 = 1; file0 < argc-1; file0 += 2)
		if (!isopt(argv[file0]))
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
	int  width = 79;
	int  length = 21;
	char  stmp[256];
	int  i;
						/* start fresh */
	mgclearall();
						/* load file */
	mgload(file);
						/* do options */
	for (i = 0; i < optc; i += 2)
		if (optv[i][0] == '+') {
			sprintf(stmp, "%s=%s", optv[i]+1, optv[i+1]);
			setmgvar("command line", stdin, stmp);
		} else
			switch (optv[i][1]) {
			case 'w':
				width = atoi(optv[i+1]);
				break;
			case 'l':
				length = atoi(optv[i+1]);
				break;
			default:
				fprintf(stderr, "%s: unknown option: %s\n",
						progname, optv[i]);
				quit(1);
			}

						/* graph it */
	cgraph(width, length);
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
