#ifndef lint
static const char	RCSid[] = "$Id: total.c,v 1.2 2003/06/08 12:03:09 schorsch Exp $";
#endif
/*
 *  total.c - program to reduce columns of data.
 *
 *	5/18/88
 */

#include  <stdio.h>
#include  <stdlib.h>
#include  <ctype.h>
#include  <math.h>

#define  MAXCOL		256		/* maximum number of columns */

#define  ADD		0		/* add numbers */
#define  MULT		1		/* multiply numbers */
#define  MAX		2		/* maximum */
#define  MIN		3		/* minimum */

double  init_val[] = {0., 1., -1e12, 1e12};	/* initial values */

int  func = ADD;			/* default function */
double  power = 0.0;			/* power for sum */
int  mean = 0;				/* compute mean */
int  count = 0;				/* number to sum */
int  tabc = '\t';			/* default separator */
int  subtotal = 0;			/* produce subtotals? */

static int execute(char *fname);

int
main(
int  argc,
char  *argv[]
)
{
	int  status;
	int  a;

	for (a = 1; a < argc; a++)
		if (argv[a][0] == '-')
			if (argv[a][1] >= '0' && argv[a][1] <= '9')
				count = atoi(argv[a]+1);
			else
				switch (argv[a][1]) {
				case 'm':
					mean = !mean;
					break;
				case 'p':
					func = MULT;
					break;
				case 's':
					func = ADD;
					power = atof(argv[a]+2);
					break;
				case 'u':
					func = MAX;
					break;
				case 'l':
					func = MIN;
					break;
				case 't':
					tabc = argv[a][2];
					break;
				case 'r':
					subtotal = !subtotal;
					break;
				default:
					fprintf(stderr,
		"Usage: %s [-m][-sE|p|u|l][-tC][-count [-r]] [file..]\n",
							argv[0]);
					exit(1);
				}
		else
			break;
	if (mean) {
		if (func == MAX) {
			fprintf(stderr, "%s: average maximum?!\n", argv[0]);
			exit(1);
		}
		if (func == MIN) {
			fprintf(stderr, "%s: average minimum?!\n", argv[0]);
			exit(1);
		}
	}
	status = 0;
	if (a >= argc)
		status = execute(NULL) == -1 ? 1 : status;
	else
		for ( ; a < argc; a++)
			status = execute(argv[a]) == -1 ? 2 : status;
	exit(status);
}


static int
execute(			/* compute result */
char  *fname
)
{
	double  result[MAXCOL];
	char  buf[16*MAXCOL];
	register char  *cp, *num;
	register int  n;
	double  d;
	int  ncol;
	long  nlin, ltotal;
	FILE  *fp;
							/* open file */
	if (fname == NULL)
		fp = stdin;
	else if ((fp = fopen(fname, "r")) == NULL) {
		fprintf(stderr, "%s: cannot open\n", fname);
		return(-1);
	}
	ltotal = 0;
	while (!feof(fp)) {
		if (ltotal == 0)			/* initialize */
			if (func == MULT)	/* special case */
				for (n = 0; n < MAXCOL; n++)
					result[n] = 0.0;
			else
				for (n = 0; n < MAXCOL; n++)
					result[n] = init_val[func];
		ncol = 0;
		for (nlin = 0; (count <= 0 || nlin < count) &&
				(cp = fgets(buf, sizeof(buf), fp)) != NULL;
				nlin++) {
							/* compute */
			for (n = 0; n < MAXCOL; n++) {
				while (isspace(*cp))
					cp++;
				if (!*cp)
					break;
				num = cp;
				while (*cp && !(isspace(tabc)?isspace(*cp):*cp==tabc))
					cp++;
				if (*cp)
					*cp++ = '\0';
				if (!*num || isspace(*num))
					d = init_val[func];
				else
					d = atof(num);
				switch (func) {
				case ADD:
					if (d == 0.0)
						break;
					if (power != 0.0)
						result[n] += pow(fabs(d),power);
					else
						result[n] += d;
					break;
				case MULT:
					if (d == 0.0)
						break;
					result[n] += log(fabs(d));
					break;
				case MAX:
					if (d > result[n])
						result[n] = d;
					break;
				case MIN:
					if (d < result[n])
						result[n] = d;
					break;
				}
			}
			if (n > ncol)
				ncol = n;
		}
		if (nlin == 0)
			break;
						/* compute and print */
		ltotal += nlin;
		for (n = 0; n < ncol; n++) {
			d = result[n];
			if (mean) {
				d /= (double)ltotal;
				if (func == ADD && power != 0.0 && d != 0.0)
					d = pow(d, 1.0/power);
			}
			if (func == MULT)
				d = exp(d);
			printf("%.9g%c", d, tabc);
		}
		putchar('\n');
		if (!subtotal)
			ltotal = 0;
	}
							/* close file */
	return(fclose(fp));
}
