#ifndef lint
static const char	RCSid[] = "$Id: total.c,v 1.14 2019/12/28 18:05:13 greg Exp $";
#endif
/*
 *  total.c - program to reduce columns of data.
 *
 *	5/18/88
 */

#include  <stdlib.h>
#include  <ctype.h>
#include  <math.h>
#include  "platform.h"
#include  "rtio.h"

#define  MAXCOL		8192		/* maximum number of columns */

#define  ADD		0		/* add numbers */
#define  MULT		1		/* multiply numbers */
#define  MAX		2		/* maximum */
#define  MIN		3		/* minimum */

double  init_val[] = {0., 1., -1e12, 1e12};	/* initial values */

long  incnt = 0;			/* limit number of input records? */
long  outcnt = 0;			/* limit number of output records? */

int  func = ADD;			/* default function */
double  power = 0.0;			/* power for sum */
int  mean = 0;				/* compute mean */
int  count = 0;				/* number to sum */
int  nbicols = 0;			/* number of binary input columns */
int  bocols = 0;			/* produce binary output columns */
int  tabc = '\t';			/* default separator */
int  subtotal = 0;			/* produce subtotals? */
int  nrecsout = 0;			/* number of records produced */

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
				case 'i':
					switch (argv[a][2]) {
					case 'n':
						incnt = atol(argv[++a]);
						break;
					case 'a':
						nbicols = 0;
						break;
					case 'd':
						if (isdigit(argv[a][3]))
							nbicols = atoi(argv[a]+3);
						else
							nbicols = 1;
						if (nbicols > MAXCOL) {
							fprintf(stderr,
							"%s: too many input columns\n",
								argv[0]);
							exit(1);
						}
						break;
					case 'f':
						if (isdigit(argv[a][3]))
							nbicols = -atoi(argv[a]+3);
						else
							nbicols = -1;
						if (-nbicols > MAXCOL) {
							fprintf(stderr,
							"%s: too many input columns\n",
								argv[0]);
							exit(1);
						}
						break;
					default:
						goto userr;
					}
					break;
				case 'o':
					switch (argv[a][2]) {
					case 'n':
						outcnt = atol(argv[++a]);
						break;
					case 'a':
						bocols = 0;
						break;
					case 'd':
						bocols = 1;
						break;
					case 'f':
						bocols = -1;
						break;
					default:
						goto userr;
					}
					break;
				default:;
				userr:
					fprintf(stderr,
"Usage: %s [-m][-sE|p|u|l][-tC][-i{f|d}[N]][-o{f|d}][-count [-r]] [file..]\n",
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
	if (bocols)
		SET_FILE_BINARY(stdout);
	status = 0;
	if (a >= argc)
		status = execute(NULL) == -1 ? 1 : status;
	else
		for ( ; a < argc && (outcnt <= 0 || nrecsout < outcnt); a++)
			status = execute(argv[a]) == -1 ? 2 : status;
	exit(status);
}


static int
getrecord(			/* read next input record */
	double field[MAXCOL],
	FILE *fp
)
{
	char  buf[16*MAXCOL];
	char  *cp, *num;
	int   nf;
						/* reading binary input? */
	if (nbicols > 0)
		return(getbinary(field, sizeof(double), nbicols, fp));
	if (nbicols < 0) {
		float	*fbuf = (float *)buf;
		int	i;
		nf = getbinary(fbuf, sizeof(float), -nbicols, fp);
		for (i = nf; i-- > 0; )
			field[i] = fbuf[i];
		return(nf);
	}
						/* reading ASCII input */
	cp = fgets(buf, sizeof(buf), fp);
	if (cp == NULL)
		return(0);
	for (nf = 0; nf < MAXCOL; nf++) {
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
			field[nf] = init_val[func];
		else
			field[nf] = atof(num);
	}
	return(nf);
}


static void
putrecord(			/* write out results record */
	const double *field,
	int n,
	FILE *fp
)
{
						/* binary output? */
	if (bocols > 0) {
		putbinary(field, sizeof(double), n, fp);
		return;
	}
	if (bocols < 0) {
		float	fv;
		while (n-- > 0) {
			fv = *field++;
			putbinary(&fv, sizeof(float), 1, fp);
		}
		return;
	}
						/* ASCII output */
	while (n-- > 0) {
		fprintf(fp, "%.9g", *field++);
		if (n) fputc(tabc, fp);
	}
	fputc('\n', fp);
}


static int
execute(			/* compute result */
char  *fname
)
{
	double	inpval[MAXCOL];
	double	tally[MAXCOL];
	short	rsign[MAXCOL];
	double  result[MAXCOL];
	int  n;
	int  nread, ncol;
	long  nlin, ltotal;
	FILE  *fp;
							/* open file */
	if (fname == NULL)
		fp = stdin;
	else if ((fp = fopen(fname, "r")) == NULL) {
		fprintf(stderr, "%s: cannot open\n", fname);
		return(-1);
	}
	if (nbicols)
		SET_FILE_BINARY(fp);
#ifdef getc_unlocked				/* avoid lock/unlock overhead */
	flockfile(fp);
#endif

	ltotal = 0;
	while (!feof(fp)) {
		if (ltotal == 0) {			/* initialize */
			if (func == MULT)	/* special case */
				for (n = 0; n < MAXCOL; n++) {
					tally[n] = 0.0;
					rsign[n] = 1;
				}
			else
				for (n = 0; n < MAXCOL; n++)
					tally[n] = init_val[func];
		}
		ncol = 0;
		for (nlin = 0; (count <= 0 || nlin < count) &&
				(incnt <= 0 || nlin < incnt) &&
				(nread = getrecord(inpval, fp)) > 0;
				nlin++) {
							/* compute */
			for (n = 0; n < nread; n++)
				switch (func) {
				case ADD:
					if (inpval[n] == 0.0)
						break;
					if (power != 0.0)
						tally[n] += pow(fabs(inpval[n]),power);
					else
						tally[n] += inpval[n];
					break;
				case MULT:
					if (inpval[n] == 0.0)
						rsign[n] = 0;
					else if (inpval[n] < 0.0) {
						rsign[n] = -rsign[n];
						inpval[n] = -inpval[n];
					}
					if (rsign[n])
						tally[n] += log(inpval[n]);
					break;
				case MAX:
					if (inpval[n] > tally[n])
						tally[n] = inpval[n];
					break;
				case MIN:
					if (inpval[n] < tally[n])
						tally[n] = inpval[n];
					break;
				}
			if (nread > ncol)
				ncol = nread;
		}
		if (nlin == 0)
			break;
						/* compute and print */
		ltotal += nlin;
		for (n = 0; n < ncol; n++) {
			result[n] = tally[n];
			if (mean) {
				result[n] /= (double)ltotal;
				if (func == ADD && power != 0.0 && result[n] != 0.0)
					result[n] = pow(result[n], 1.0/power);
			}
			if (func == MULT)
				result[n] = rsign[n] * exp(result[n]);
		}
		putrecord(result, ncol, stdout);
		++nrecsout;
		if (outcnt > 0 && nrecsout >= outcnt)
			break;
		if (!subtotal)
			ltotal = 0;
		if (incnt > 0 && nlin >= incnt)
			break;
	}
							/* close input */
	return(fclose(fp));
}
