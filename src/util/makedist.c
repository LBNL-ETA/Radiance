#ifndef lint
static const char	RCSid[] = "$Id: makedist.c,v 2.5 2004/03/26 21:36:19 schorsch Exp $";
#endif
/*
 *  makedist.c - program to make a source distribution.
 *
 *     8/18/86
 *
 *	Program uses angular coordinates as follows:
 *
 *		alpha - degrees measured from x1 axis.
 *		beta - degress of projection into x2x3 plane,
 *			measured from x2 axis towards x3 axis.
 *
 *		Default axes are (x1,x2,x3) = (x,y,z).
 */

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>

#include  "rtmath.h"
#include  "random.h"
#include  "setscan.h"

#ifndef BSD
#define vfork		fork
#endif

#define  atorad(a)	((PI/180.0) * (a))

char  *rtargv[128] = {"rtrace", "-h"};
int  rtargc = 2;

#define  passarg(s)	(rtargv[rtargc++] = s)

ANGLE  alpha[181] = {10, 25, 40, 55, 70, 85, AEND};
ANGLE  beta[361] = {0,30,60,90,120,150,180,210,240,270,300,330,AEND};

					/* coordinate system */
double  axis[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

double  targetw = 1.0;			/* target width */
double  targeth = 1.0;			/* target height */
double  targetd = 1.0;			/* target distance */
double  targetp[3] = {0,0,0};		/* target center */
int  xres = 16;				/* x sample resolution */
int  yres = 16;				/* y sample resolution */

int  userformat = 1;			/* output 1=nice,0=plain,-1=raw */
int  header = 1;			/* print header? */

char  *progname;			/* program name */

static void printargs(int ac, char **av, FILE *fp);
static void fixaxes(void);
static int doscan(void);
static FILE * scanstart(void);
static int scanend(FILE *fp);
static void sendsamples(FILE *fp);
static void writesample(FILE *fp, ANGLE a, ANGLE b);
static double readsample(FILE *fp);


int
main(
	int  argc,
	char  *argv[]
)
{
	int  i;

	progname = argv[0];

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 't':			/* target */
			switch (argv[i][2]) {
			case 'w':			/* width */
				targetw = atof(argv[++i]);
				break;
			case 'h':			/* height */
				targeth = atof(argv[++i]);
				break;
			case 'd':			/* distance */
				targetd = atof(argv[++i]);
				break;
			case 'c':			/* center */
				targetp[0] = atof(argv[++i]);
				targetp[1] = atof(argv[++i]);
				targetp[2] = atof(argv[++i]);
				break;
			default:
				goto badopt;
			}
			break;
		case 'a':
			if (!strcmp(argv[i], "-alpha")) {
				if (setscan(alpha, argv[++i]) < 0) {
					fprintf(stderr,
						"%s: bad alpha spec: %s\n",
							progname, argv[i]);
					exit(1);
				}
			} else
				switch (argv[i][2]) {
				case 'd':		/* ambient divisions */
				case 's':		/* ambient samples */
				case 'r':		/* ambient resolution */
				case 'a':		/* ambient accuracy */
				case 'b':		/* ambient bounces */
				case 'i':		/* include */
				case 'e':		/* exclude */
				case 'f':		/* file */
					passarg(argv[i]);
					passarg(argv[++i]);
					break;
				case 'v':		/* ambient value */
					passarg(argv[i]);
					passarg(argv[++i]);
					passarg(argv[++i]);
					passarg(argv[++i]);
					break;
				default:
					goto badopt;
				}
			break;
		case 'b':
			if (!strcmp(argv[i], "-beta")) {
				if (setscan(beta, argv[++i]) < 0) {
					fprintf(stderr,
						"%s: bad beta spec: %s\n",
							progname, argv[i]);
					exit(1);
				}
			} else
				goto badopt;
			break;
		case 'x':
			switch (argv[i][2]) {
			case '\0':			/* x resolution */
				xres = atoi(argv[++i]);
				break;
			case '1':
			case '2':			/* axis spec */
			case '3':
				axis[argv[i][2]-'1'][0] = atof(argv[i+1]);
				axis[argv[i][2]-'1'][1] = atof(argv[i+2]);
				axis[argv[i][2]-'1'][2] = atof(argv[i+3]);
				i += 3;
				break;
			default:
				goto badopt;
			}
			break;
		case 'y':			/* y resolution */
			yres = atoi(argv[++i]);
			break;
		case 'l':
			switch (argv[i][2]) {
			case 'w':			/* limit weight */
			case 'r':			/* limit recursion */
				passarg(argv[i]);
				passarg(argv[++i]);
				break;
			default:
				goto badopt;
			}
			break;
		case 'u':			/* user friendly format */
			userformat = 1;
			break;
		case 'o':			/* custom format */
			passarg(argv[i]);
			userformat = -1;
			break;
		case 'd':
			switch (argv[i][2]) {
			case '\0':			/* data format */
				userformat = 0;
				break;
			case 'j':			/* direct jitter */
			case 't':			/* direct threshold */
			case 'c':			/* direct certainty */
				passarg(argv[i]);
				passarg(argv[++i]);
				break;
			default:
				goto badopt;
			}
			break;
		case 'h':			/* no header */
			header = 0;
			break;
		default:;
badopt:
			fprintf(stderr, "%s: bad option: %s\n",
					progname, argv[i]);
			exit(1);
		}

	if (userformat >= 0) {		/* we cook */
		passarg("-ov");
		passarg("-fff");
	} else				/* raw */
		passarg("-ffa");

	if (i != argc-1) {
		fprintf(stderr, "%s: single octree required\n", progname);
		exit(1);
	}

	passarg(argv[i]);
	passarg(NULL);

	fixaxes();

	if (header) {
		if (userformat > 0) {
			printf("Alpha\tBeta\tRadiance\n");
		} else {
			printargs(argc, argv, stdout);
			putchar('\n');
		}
	}

	if (doscan() == -1)
		exit(1);

	exit(0);
}


static void
printargs(		/* print arguments to a file */
	int  ac,
	char  **av,
	FILE  *fp
)
{
	while (ac-- > 0) {
		fputs(*av++, fp);
		putc(' ', fp);
	}
	putc('\n', fp);
}


static void
fixaxes(void)				/* make sure axes are OK */
{
	double  d;
	register int  i;

	for (i = 0; i < 3; i++)
		if (normalize(axis[i]) == 0.0)
			goto axerr;

	for (i = 0; i < 3; i++) {
		d = fdot(axis[i], axis[(i+1)%3]);
		if (d > FTINY || d < -FTINY)
			goto axerr;
	}
	return;
axerr:
	fprintf(stderr, "%s: bad axis specification\n", progname);
	exit(1);
}


static int
doscan(void)				/* do scan for target */
{
	FILE  *fopen(), *scanstart();
	double  readsample();
	FILE  *fp;
	ANGLE  *a, *b;

	fp = scanstart();		/* returns in child if not raw */
	for (a = alpha; *a != AEND; a++)		/* read data */
		for (b = beta; *b != AEND; b++)
			if (userformat > 0)
				printf("%d\t%d\t%e\n", *a, *b, readsample(fp));
			else
				printf("%e\n", readsample(fp));
	return(scanend(fp));
}


static FILE *
scanstart(void)			/* open scanner pipeline */
{
	int  status;
	char  *fgets();
	int  p1[2], p2[2];
	int  pid;
	FILE  *fp;

	fflush(stdout);				/* empty output buffer */
	if (pipe(p1) < 0)			/* get pipe1 */
		goto err;
	if ((pid = fork()) == 0) {		/* if child1 */
		close(p1[1]);
		if (p1[0] != 0) {		/* connect pipe1 to stdin */
			dup2(p1[0], 0);
			close(p1[0]);
		}
		if (userformat >= 0) {		/* if processing output */
			if (pipe(p2) < 0)		/* get pipe2 */
				goto err;
			if ((pid = vfork()) == 0) {	/* if child2 */
				close(p2[0]);
				if (p2[1] != 1) {	/* pipe2 to stdout */
					dup2(p2[1], 1);
					close(p2[1]);
				}
				execvp(rtargv[0], rtargv);	/* rtrace */
				goto err;
			}
			if (pid == -1)
				goto err;
			close(p2[1]);
			if ((fp = fdopen(p2[0], "r")) == NULL)
				goto err;
			return(fp);
		} else {				/* if running raw */
			execvp(rtargv[0], rtargv);		/* rtrace */
			goto err;
		}
	}
	if (pid == -1)
		goto err;
	close(p1[0]);
	if ((fp = fdopen(p1[1], "w")) == NULL)
		goto err;
	sendsamples(fp);		/* parent writes output to pipe1 */
	fclose(fp);
	if (wait(&status) != -1)	/* wait for child1 */
		_exit(status);
err:
	perror(progname);
	_exit(1);
}


static int
scanend(		/* done with scanner input */
	FILE  *fp
)
{
	int  status;

	fclose(fp);
	if (wait(&status) == -1) {	/* wait for child2 */
		perror(progname);
		exit(1);
	}
	return(status ? -1 : 0);
}


static void
sendsamples(			/* send our samples to fp */
	FILE  *fp
)
{
	ANGLE  *a, *b;

	for (a = alpha; *a != AEND; a++)
		for (b = beta; *b != AEND; b++)
			writesample(fp, *a, *b);
}


static void
writesample(		/* write out sample ray grid */
	FILE	*fp,
	ANGLE	a,
	ANGLE	b
)
{
	double  sin(), cos();
	float  sp[6];
	double  dv[3];
	double  xpos, ypos;
	int  i, j, k;
					/* get direction in their system */
	dv[0] = -cos(atorad(a));
	dv[1] = dv[2] = -sin(atorad(a));
	dv[1] *= cos(atorad(b));
	dv[2] *= sin(atorad(b));
					/* do sample grid in our system */
	for (j = 0; j < yres; j++) {
		for (i = 0; i < xres; i++) {
			xpos = ((i + frandom())/xres - 0.5) * targetw;
			ypos = ((j + frandom())/yres - 0.5) * targeth;
			for (k = 0; k < 3; k++) {
				sp[k] = targetp[k];
				sp[k] += xpos*axis[2][k];
				sp[k] += ypos*axis[1][k];
				sp[k+3] = dv[0]*axis[0][k] +
						dv[1]*axis[1][k] +
						dv[2]*axis[2][k];
				sp[k] -= sp[k+3]*targetd;
			}
			if (fwrite(sp, sizeof(*sp), 6, fp) != 6) {
				fprintf(stderr, "%s: data write error\n",
						progname);
				_exit(1);
			}
		}
	}
}


static double
readsample(			/* read in sample ray grid */
	FILE  *fp
)
{
	float  col[3];
	double  sum;
	int  i;

	sum = 0.0;
	i = xres*yres;
	while (i--) {
		if (fread(col, sizeof(*col), 3, fp) != 3) {
			fprintf(stderr, "%s: data read error\n", progname);
			exit(1);
		}
		sum += .3*col[0] + .59*col[1] + .11*col[2];
	}
	return(sum / (xres*yres));
}
