#ifndef lint
static const char	RCSid[] = "$Id: igraph.c,v 1.2 2003/10/27 10:28:59 schorsch Exp $";
#endif
/*
 *  igraph.c - interactive graphing program.
 *
 *     6/30/86
 *
 *     Greg Ward Larson
 */

#include  <stdio.h>

#include  <ctype.h>

#include  <setjmp.h>

#include  "rtprocess.h"
#include  "mgvars.h"

typedef struct {
	char  *descrip;		/* description */
	char  *value;		/* value */
}  SMAP;

#define  NDEV		10	/* number of devices */

SMAP  dev[NDEV] = {
	{ "Tektronix 4014", "t4014" },
	{ "AED 512", "output aed5" },
	{ "X10 Window System", "xmeta" },
	{ "X11 Window System", "x11meta" },
	{ "Epson Printer", "output mx80" },
	{ "Mannesman-Tally Printer", "output mt160" },
	{ "Mannesman-Tally High Density", "output mt160l" },
	{ "Apple Imagewriter", "output imagew" },
	{ "Imagen Printer", "impress | ipr" },
	{ "Postscript Printer", "psmeta | lpr" },
};

#define  NTYP		7	/* number of plot types */

SMAP  typ[NTYP] = {
	{ "Cartesian Plot", "cartesian.plt" },
	{ "Polar Plot (degrees)", "polar.plt" },
	{ "Curve Plot (symbols & lines)", "curve.plt" },
	{ "Scatter Plot (symbols only)", "scatter.plt" },
	{ "Line Plot (lines only)", "line.plt" },
	{ "Function Plot (lines)", "function.plt" },
	{ "Box and Whisker Plot", "boxw.plt" },
};

#define  NCAL		3	/* number of calculation types */

SMAP  cal[NCAL] = {
	{ "Extrema, Average", "nma" },
	{ "Integral", "ni" },
	{ "Linear Regression", "nl" },
};

char  *progname;

char  *libpath[4];

static int  recover = 0;	/* try to recover? */

static jmp_buf  env;		/* setjmp buffer */

static char  curfile[64] = "temp.plt";	/* current file name */

char  *getfile(), *gets(), *strcpy();


main(argc, argv)
int  argc;
char  *argv[];
{
	char  *getenv();
	int  i;

	progname = argv[0];
	libpath[0] = "./";
	if ((libpath[i=1] = getenv("MDIR")) != NULL)
		i++;
	libpath[i++] = MDIR;
	libpath[i] = NULL;

	for (i = 1; i < argc; i++)
		mgload(strcpy(curfile, argv[i]));

	mainmenu();

	quit(0);
}


eputs(msg)				/* print error message */
char  *msg;
{
	fputs(msg, stderr);
}


quit(code)				/* recover or quit */
int  code;
{
	if (code && recover--)
		longjmp(env, 1);
	exit(code);
}


mainmenu()			/* the main menu loop */
{
	char  sbuf[128];

	for ( ; ; ) {
		printf("\nGRAPH MENU\n");
		printf("\t1	New Plot\n");
		printf("\t2	Set Plot Type\n");
		printf("\t3	Load Plot File\n");
		printf("\t4	Save Plot File\n");
		printf("\t5	Change Plot Variable\n");
		printf("\t6	Change Curve Variable\n");
		printf("\t7	Output Plot\n");
		printf("\t8	Escape Command\n");
		printf("\t9	Computation\n");
		printf("\t10	Quit\n");

		printf("\nChoice: ");
		clearerr(stdin); gets(sbuf);
		if (feof(stdin))
			return;
		switch (atoi(sbuf)) {
		case 1:				/* new plot */
			mgclearall();
			break;
		case 2:				/* set plot type */
			settype();
			break;
		case 3:				/* load plot file */
			if (setjmp(env) == 0) {
				recover = 1;
				mgload(getfile());
				recover = 0;
			} else
				uwait();
			break;
		case 4:				/* save plot file */
			if (setjmp(env) == 0) {
			    recover = 1;
			    mgsave(getfile());
			    recover = 0;
			} else
				uwait();
			break;
		case 5:				/* set plot variable */
			setvars(NVARS, gparam);
			break;
		case 6:				/* set curve variable */
			printf("\nCurve (A-%c): ", MAXCUR-1+'A');
			clearerr(stdin); gets(sbuf);
			if (islower(sbuf[0]))
				sbuf[0] -= 'a';
			else if (isupper(sbuf[0]))
				sbuf[0] -= 'A';
			else
				break;
			if (sbuf[0] < MAXCUR)
				setvars(NCVARS, cparam[sbuf[0]]);
			break;
		case 7:				/* output plot */
			plotout();
			break;
		case 8:				/* escape command */
			printf("\nCommand: ");
			clearerr(stdin); gets(sbuf);
			if (sbuf[0]) {
				system(sbuf);
				uwait();
			}
			break;
		case 9:				/* computation */
			gcompute();
			break;
		case 10:			/* quit */
			return;
		}
	}
}


char *
getfile()			/* get file name from user */
{
	char  sbuf[128];

	printf("\nFile (%s): ", curfile);
	clearerr(stdin); gets(sbuf);
	if (sbuf[0] == '-')
		return(NULL);
	if (sbuf[0])
		strcpy(curfile, sbuf);
	return(curfile);
}


settype()			/* set plot type */
{
	char  sbuf[128];
	int  i;

	printf("\nPLOT TYPE\n");

	for (i = 0; i < NTYP; i++)
		printf("\t%d	%s\n", i+1, typ[i].descrip);
	
	printf("\nChoice (0): ");
	clearerr(stdin); gets(sbuf);
	i = atoi(sbuf) - 1;
	if (i < 0 || i >= NTYP)
		return;
	sprintf(sbuf, "include=%s", typ[i].value);
	setmgvar("plot type", stdin, sbuf);
}


plotout()			/* output our graph */
{
	extern FILE  *pout;
	char  sbuf[128];
	char  *command;
	int  i;

	printf("\nOUTPUT PLOT\n");

	for (i = 0; i < NDEV; i++)
		printf("\t%d	%s\n", i+1, dev[i].descrip);
	printf("\t%d	Dumb Terminal\n", NDEV+1);
	printf("\t%d	Other\n", NDEV+2);
	
	printf("\nChoice (0): ");
	clearerr(stdin); gets(sbuf);
	i = atoi(sbuf) - 1;
	if (i < 0 || i > NDEV+1)
		return;
	if (i < NDEV)
		command = dev[i].value;
	else if (i == NDEV) {
		if (setjmp(env) == 0) {
			recover = 1;
			cgraph(79, 22);
			recover = 0;
		}
		uwait();
		return;
	} else if (i == NDEV+1) {
		printf("\nDriver command: ");
		clearerr(stdin); gets(sbuf);
		if (!sbuf[0])
			return;
		command = sbuf;
	}
	if ((pout = popen(command, "w")) == NULL) {
		fprintf(stderr, "%s: cannot execute device driver: %s\n",
				progname, dev[i].value);
		quit(1);
	}
	if (setjmp(env) == 0) {
		recover = 1;
		mgraph();
		recover = 0;
	} else
		uwait();
	mdone();
	pclose(pout);
}


setvars(nvar, vars)		/* set variables */
int  nvar;
VARIABLE  vars[];
{
	char  sbuf[128];
	int  i, j;

	for ( ; ; ) {
		printf("\nVARIABLES\n");

		for (i = 0; i < nvar; i++)
			printf("\t%d	%s\n", i+1, vars[i].descrip == NULL ?
					vars[i].name : vars[i].descrip);

		printf("\nChoice (0): ");
		clearerr(stdin); gets(sbuf);
		i = atoi(sbuf) - 1;
		if (i < 0 || i >= nvar)
			return;
		printf("\n%s", vars[i].name);
		if (vars[i].type == FUNCTION)
			printf("(x)");
		if (vars[i].flags & DEFINED) {
			mgtoa(sbuf, &vars[i]);
			printf(" (%s)", sbuf);
		}
		printf(": ");
		if (vars[i].type == FUNCTION)
			sprintf(sbuf, "%s(x)=", vars[i].name);
		else
			sprintf(sbuf, "%s=", vars[i].name);
		j = strlen(sbuf);
		clearerr(stdin); gets(sbuf+j);
		if (sbuf[j]) {
			if (vars[i].type == DATA && sbuf[j] == '-')
				sbuf[j] = '\0';

			if (setjmp(env) == 0) {
				recover = 1;
				setmgvar("variable set", stdin, sbuf);
				recover = 0;
			} else
				uwait();
		}
	}
}


gcompute()			/* do a computation */
{
	char  sbuf[64];
	int  i;

	printf("\nCOMPUTE\n");

	for (i = 0; i < NCAL; i++)
		printf("\t%d	%s\n", i+1, cal[i].descrip);
	
	printf("\nChoice (0): ");
	clearerr(stdin); gets(sbuf);
	i = atoi(sbuf) - 1;
	if (i < 0 || i >= NCAL)
		return;
	printf("\n");
	if (setjmp(env) == 0) {
		recover = 1;
		gcalc(cal[i].value);
		recover = 0;
	}
	printf("\n");
	uwait();
}


uwait()				/* wait for user after command, error */
{
	printf("Hit return to continue: ");
	while (getchar() != '\n')
		;
}
