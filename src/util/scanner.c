#ifndef lint
static const char	RCSid[] = "$Id: scanner.c,v 1.5 2004/03/26 23:34:23 schorsch Exp $";
#endif
/*
 *  scanner.c - program to simulate bi-directional scanner.
 *
 *     6/10/86
 */

#include  <stdio.h>
#include  <string.h>
#include  <stdlib.h>
#include  <ctype.h>
#include  <signal.h>
#include  <math.h>
#include  <unistd.h> /* XXX pipe() fork() etc, use common/process.c instead */

#include  "random.h"


#define  PI		3.14159265358979324

#define  ANGLE		short
#define  AEND		-1

#define  atorad(a)	((PI/180.0) * (a))

char  otcom[128] = "oconv";
char  rtcom[256] = "rtrace -bf -ar .25 -ad 128";

char  *sourcetemp,			/* temp files */
      *octreetemp;

ANGLE  alpha[181] = {90, AEND};
ANGLE  beta[181] = {30, 60, 90, 120, 150, AEND};
ANGLE  sgamma[181] = {45, AEND}; /* gamma() is a function from gnu <math.h> */

char  *target;				/* target file name */
double  targetw = 3.0;			/* target width (inches) */
double  targeth = 3.0;			/* target height (inches) */
int  xres = 16;				/* x sample resolution */
int  yres = 16;				/* y sample resolution */

static void  quit(int code);
static void setscan(ANGLE  *ang, char  *arg);
static void doscan(void);
static void scanend(FILE  *fp);
static void makesource(ANGLE  g);
static void sendsamples(FILE  *fp);
static void writesample(FILE	*fp, ANGLE	a, ANGLE	b);
static FILE * scanstart(ANGLE  g);
static double readsample(FILE  *fp);

int
main(
	int  argc,
	char  *argv[]
)
{
	int  i;

#ifdef SIGHUP
	signal(SIGHUP, quit);
#endif
	signal(SIGINT, quit);
	signal(SIGTERM, quit);
#ifdef SIGXCPU
	signal(SIGXCPU, SIG_IGN);
	signal(SIGXFSZ, SIG_IGN);
#endif

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
			default:
				goto badopt;
			}
			break;
		case 'a':
			if (!strcmp(argv[i], "-alpha"))
				setscan(alpha, argv[++i]);
			else
				switch (argv[i][2]) {
				case 'd':		/* ambient divisions */
					strcat(rtcom, " -ad ");
					strcat(rtcom, argv[++i]);
					break;
				case 'r':		/* ambient radius */
					strcat(rtcom, " -ar ");
					strcat(rtcom, argv[++i]);
					break;
				default:
					goto badopt;
				}
			break;
		case 'b':
			if (!strcmp(argv[i], "-beta"))
				setscan(beta, argv[++i]);
			else
				goto badopt;
			break;
		case 'g':
			if (!strcmp(argv[i], "-gamma"))
				setscan(sgamma, argv[++i]);
			else
				goto badopt;
			break;
		case 'x':			/* x resolution */
			xres = atoi(argv[++i]);
			break;
		case 'y':			/* y resolution */
			yres = atoi(argv[++i]);
			break;
		case 'l':
			switch (argv[i][2]) {
			case 'w':			/* limit weight */
				strcat(rtcom, " -lw ");
				strcat(rtcom, argv[++i]);
				break;
			case 'r':			/* limit recursion */
				strcat(rtcom, " -lr ");
				strcat(rtcom, argv[++i]);
				break;
			default:
				goto badopt;
			}
			break;
		default:;
badopt:
			fprintf(stderr, "%s: bad option: %s\n",
					argv[0], argv[i]);
			quit(1);
		}

	if (i >= argc) {
		fprintf(stderr, "%s: missing target file\n", argv[0]);
		quit(1);
	} else if (i < argc-1) {
		fprintf(stderr, "%s: too many target files\n", argv[0]);
		quit(1);
	}

	target = argv[i];
	sourcetemp = mktemp("/usr/tmp/soXXXXXX");
	octreetemp = mktemp("/usr/tmp/ocXXXXXX");
	strcat(strcat(otcom, " "), target);
	strcat(strcat(otcom, " "), sourcetemp);
	strcat(strcat(otcom, " > "), octreetemp);
	strcat(strcat(rtcom, " "), octreetemp);

	doscan();

	quit(0);
	return 0; /* pro forma return */
}


void
quit(			/* unlink temp files and exit */
	int  code
)
{
	unlink(sourcetemp);
	unlink(octreetemp);
	exit(code);
}


static void
setscan(			/* set up scan according to arg */
	register ANGLE  *ang,
	register char  *arg
)
{
	int  start = 0, finish = -1, step = 1;

	for ( ; ; ) {
		switch (*arg) {
		case '\0':
		case ',':
			while (start <= finish) {
				*ang++ = start;
				start += step;
			}
			if (*arg) {
				arg++;
				continue;
			} else {
				*ang = AEND;
				return;
			}
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			start = atoi(arg);
			finish = start;
			step = 1;
			break;
		case '-':
			finish = atoi(++arg);
			break;
		case ':':
			step = atoi(++arg);
			break;
		default:
			fprintf(stderr, "Scan specification error\n");
			quit(1);
		}
		while (isdigit(*arg))
			arg++;
	}
}


static void
doscan(void)				/* do scan for target */
{
	FILE  *fp;
	ANGLE  *a, *b, *g;

	printf("Alpha\tBeta\tGamma\tDistribution for \"%s\"\n", target);
	for (g = sgamma; *g != AEND; g++) {
		fp = scanstart(*g);
		for (b = beta; *b != AEND; b++)		/* read data */
			for (a = alpha; *a != AEND; a++)
				printf("%d\t%d\t%d\t%f\n",
						*a, *b, *g, readsample(fp));
		scanend(fp);
	}
}


static FILE *
scanstart(		/* open scanner pipeline */
	ANGLE  g
)
{
	int  p1[2], p2[2];
	int  pid;
	FILE  *fp;
	char  buf[128];

	makesource(g);			/* make source */
	if (system(otcom))		/* generate new octree */
		goto err;

	if (pipe(p1) < 0)		/* get pipe1 */
		goto err;
	if ((pid = fork()) == 0) {		/* if child1 */
		close(p1[0]);
		if (p1[1] != 1) {		/* connect stdout to pipe1 */
			dup2(p1[1], 1);
			close(p1[1]);
		}
		if (pipe(p2) < 0)		/* get pipe2 */
			goto err;
		if ((pid = fork()) == 0) {	/* if child2 */
			close(p2[0]);
			if ((fp = fdopen(p2[1], "w")) == NULL)
				goto err;
			sendsamples(fp);	/* send our samples to pipe2 */
			fclose(fp);
			_exit(0);
		}
		if (pid == -1)
			goto err;
		close(p2[1]);
			
		if (p2[0] != 0) {		/* for child1 pipe2 to stdin */
			dup2(p2[0], 0);
			close(p2[0]);
		}
		execl("/bin/sh", "sh", "-c", rtcom, 0);	/* run ray tracer */
		goto err;
	}
	if (pid == -1)
		goto err;
	close(p1[1]);
	if ((fp = fdopen(p1[0], "r")) == NULL)	/* parent reads pipe1 */
		goto err;
	while (fgets(buf, sizeof(buf), fp) != NULL && buf[0] != '\n')
		/* discard header */;
	return(fp);
err:
	fprintf(stderr, "System error in scanstart\n");
	quit(1);
	return NULL; /* pro forma return */
}


static void
scanend(		/* done with scanner input */
	FILE  *fp
)
{
	fclose(fp);
}


static void
makesource(		/* make a source (output normalized) */
	ANGLE  g
)
{
	FILE  *fp;
	double  srcdir[3];

	srcdir[0] = sin(atorad(g));
	srcdir[1] = 0.0;
	srcdir[2] = -cos(atorad(g));
	if ((fp = fopen(sourcetemp, "w")) == NULL) {
		fprintf(stderr, "%s: cannot create\n", sourcetemp);
		quit(1);
	}
	fprintf(fp, "0 0 0 1\n");
	fprintf(fp, "\nvoid light scanner_light\n");
	fprintf(fp, "0\n0\n3 3283 3283 3283\n");
	fprintf(fp, "\nscanner_light source scanner_source\n");
	fprintf(fp, "0\n0\n4 %f %f %f 2\n", srcdir[0], srcdir[1], srcdir[2]);
	fclose(fp);
}


static void
sendsamples(			/* send our samples to fp */
	FILE  *fp
)
{
	ANGLE  *a, *b;

	for (b = beta; *b != AEND; b++)
		for (a = alpha; *a != AEND; a++)
			writesample(fp, *a, *b);
}


static void
writesample(		/* write out sample ray grid */
	FILE	*fp,
	ANGLE	a,
	ANGLE	b
)
{
	float  sp[6];
	int  i, j;

	sp[3] = -sin(atorad(b)-(PI/2.0)) * sin(atorad(a));
	sp[4] = -cos(atorad(a));
	sp[5] = -sin(atorad(b)) * sin(atorad(a));

	for (j = 0; j < yres; j++)
		for (i = 0; i < xres; i++) {
			sp[0] = -100.0*sp[3] +
					(i - xres/2.0 + frandom())*targetw/xres;
			sp[1] = -100.0*sp[4] +
					(j - yres/2.0 + frandom())*targeth/yres;
			sp[2] = -100.0*sp[5];
			fwrite(sp, sizeof(*sp), 6, fp);
			if (ferror(fp)) {
				fprintf(stderr, "Data write error\n");
				quit(1);
			}
		}
}


static double
readsample(			/* read in sample ray grid */
	FILE  *fp
)
{
	double  sum;
	float  col[3];
	int  i;

	sum = 0.0;
	i = xres*yres;
	while (i--) {
		if (fread(col, sizeof(*col), 3, fp) != 3) {
			fprintf(stderr, "Data read error\n");
			quit(1);
		}
		sum += col[1];
	}
	return(sum / (xres*yres));
}
