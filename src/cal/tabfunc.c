#ifndef lint
static const char	RCSid[] = "$Id: tabfunc.c,v 1.2 2003/06/08 12:03:09 schorsch Exp $";
#endif
/*
 * Put tabular data into functions suitable for cal programs.
 *
 *	2/2/95	Greg Ward
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#ifdef _WIN32
 #include <process.h> /* getpid() */
#else
 #include <sys/types.h>
 #include <unistd.h>
#endif

#include "standard.h"
#include "platform.h"

#define  isdelim(c)	(isspace(c) || (c)==',')

#define  MAXTAB		1024		/* maximum number of data rows */
#define  MAXLINE	4096		/* maximum line width (characters) */
/*#define  FLOAT		float	*/	/* real type (precision) */
#define  OUTFMT		"%.7g"		/* output format conversion string */

int	interpolate = 0;
char	*progname;
char	**func;
int	nfuncs;

FLOAT	abscissa[MAXTAB];		/* independent values (first column) */
FLOAT	(*ordinate)[MAXTAB];		/* dependent values (other columns) */
int	tabsize = 0;			/* final table size (number of rows) */
char	locID[16];			/* local identifier (for uniqueness) */

/*extern char	*fgets(), *fskip(), *absc_exp();*/

static void load_data(FILE *fp);
static void print_funcs(char *xe);
static void putlist(register FLOAT *av, int al, register int pos);
static char * absc_exp(void);

int
main(
int	argc,
char	**argv
)
{
	progname = argv[0];
	argv++;
	argc--;
	if (argc && !strcmp(argv[0], "-i")) {
		interpolate++;
		puts("interp_arr2`(i,x,f):(i+1-x)*f(i)+(x-i)*f(i+1);");
		puts("interp_arr`(x,f):if(x-1,if(f(0)-x,interp_arr2`(floor(x),x,f),f(f(0))),f(1));");
		argv++;
		argc--;
	}
	if (!argc || argv[0][0] == '-') {
		fprintf(stderr, "Usage: %s [-i] func1 [func2 ..]\n", progname);
		exit(1);
	}
	func = argv;
	nfuncs = argc;
	ordinate = (FLOAT (*)[MAXTAB])malloc(nfuncs*MAXTAB*sizeof(FLOAT));
	if (ordinate == NULL) {
		fprintf(stderr, "%s: not enough memory\n", progname);
		exit(1);
	}
	sprintf(locID, "p%d", getpid());
	load_data(stdin);
	print_funcs(absc_exp());
	exit(0);
}


static void
load_data(			/* load tabular data from fp */
FILE	*fp
)
{
	int	lineno;
	char	*err;
	char	inpbuf[MAXLINE];
	register char	*cp;
	register int	i;

	tabsize = lineno = 0;
	inpbuf[MAXLINE-2] = '\n';
	while (fgets(inpbuf, MAXLINE, fp) != NULL) {
		lineno++;
		if (inpbuf[MAXLINE-2] != '\n') {
			err = "line too long";
			goto fatal;
		}
		if (tabsize >= MAXTAB-1) {
			err = "too many rows";
			goto fatal;
		}
		if ((cp = fskip(inpbuf)) == NULL)
			continue;		/* skip non-data lines */
		abscissa[tabsize] = atof(inpbuf);
		for (i = 0; i < nfuncs; i++) {
			while (isdelim(*cp))
				cp++;
			if (!*cp) {
				err = "too few columns";
				goto fatal;
			}
			ordinate[i][tabsize] = atof(cp);
			if ((cp = fskip(cp)) == NULL) {
				err = "bad floating-point format";
				goto fatal;
			}
		}
		tabsize++;
	}
	return;
fatal:
	fprintf(stderr, "%s: input line %d: %s\n", progname, lineno, err);
	exit(1);
}


static char *
absc_exp(void)			/* produce expression for abscissa */
{
	static char	ourexp[64];
	double	step, eps;
	int	uniform, increasing;
	register int	i;

	if (tabsize < 2)
		return("1");
	step = abscissa[1] - abscissa[0];
	eps = ((increasing = (step > 0)) ? 1e-3 : -1e-3) * step;
	uniform = 1;
	for (i = 2; i < tabsize; i++) {
		if (uniform && fabs((abscissa[i]-abscissa[i-1]) - step) > eps)
			uniform = 0;
		if (!uniform && (abscissa[i-1] < abscissa[i]) != increasing) {
			fprintf(stderr, "%s: input not a function\n",
					progname);
			exit(1);
		}
	}
	if (uniform) {
		if (increasing && fabs(step - 1) < eps) {
			if (fabs(abscissa[0] - 1) < eps)
				strcpy(ourexp, "x");
			else
				sprintf(ourexp, "x-%g", abscissa[0]-1);
		} else
			sprintf(ourexp, "(x-%g)/%g+1", abscissa[0], step);
	} else {
		printf("X`%s(i):select(i,", locID);
		putlist(abscissa, tabsize, 20);
		puts(");");
		if (increasing) {
			printf("fx`%s(x):if(x-%g,if(%g-x,fx2`%s(x,%d),%d),1);\n",
					locID, abscissa[0], abscissa[tabsize-1],
					locID, tabsize, tabsize);
			printf("fx2`%s(x,i):if(x-X`%s(i),\n", locID, locID);
		} else {
			printf("fx`%s(x):if(%g-x,if(x-%g,fx2`%s(x,%d),%d),1);\n",
					locID, abscissa[0], abscissa[tabsize-1],
					locID, tabsize, tabsize);
			printf("fx2`%s(x,i):if(X`%s(i)-x,\n", locID, locID);
		}
		printf("\ti+(x-X`%s(i))/(X`%s(i+1)-X`%s(i)),\n",
				locID, locID, locID);
		printf("\tfx2`%s(x,i-1));\n", locID);
		sprintf(ourexp, "fx`%s(x)", locID);
	}
	return(ourexp);
}


static void
print_funcs(			/* print functions */
char	*xe
)
{
	int	xelen;
	register int	i;

	xelen = strlen(xe);
	for (i = 0; i < nfuncs; i++) {
		if (func[i][0] == '\0' | func[i][0] == '0')
			continue;
		if (interpolate) {
			printf("%s`%s(i):select(i,", func[i], locID);
			putlist(ordinate[i], tabsize,
					27+strlen(func[i]));
			puts(");");
			printf("%s(x):interp_arr`(%s,%s`%s);\n",
					func[i], xe, func[i], locID);
		} else {
			printf("%s(x):select(%s,", func[i], xe);
			putlist(ordinate[i], tabsize, strlen(func[i])+xelen+12);
			puts(");");
		}
	}
}


static void
putlist(		/* put out array of values */
register FLOAT	*av,
int	al,
register int	pos
)
{
	char	obuf[32];
	int	len;

	while (al--) {
		sprintf(obuf, OUTFMT, *av++);
		pos += (len = strlen(obuf)+1);
		if (pos >= 64) {
			putchar('\n');
			pos = len;
		}
		fputs(obuf, stdout);
		if (al)
			putchar(',');
	}
}
