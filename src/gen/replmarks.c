#ifndef lint
static const char RCSid[] = "$Id";
#endif
/*
 * Replace markers in Radiance scene description with objects or instances.
 *
 *	Created:	17 Feb 1991	Greg Ward
 */

#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>

#include "platform.h"
#include "rtio.h"
#include "rtprocess.h"
#include "fvect.h"

#ifdef  M_PI
#define  PI		((double)M_PI)
#else
#define  PI		3.14159265358979323846
#endif

#define  FEQ(a,b)	((a)-(b) <= 1e-7 && (b)-(a) <= 1e-7)

#define  MAXVERT	6	/* maximum number of vertices for markers */
#define  MAXMARK	32	/* maximum number of markers */

typedef struct {
	short	beg, end;		/* beginning and ending vertex */
	float	len2;			/* length squared */
}  EDGE;			/* a marker edge */

struct mrkr {
	char	*modout;		/* output modifier */
	double	mscale;			/* scale by this to get unit */
	char	*modin;			/* input modifier indicating marker */
	char	*objname;		/* output object file or octree */
	int	doxform;		/* true if xform, false if instance */
}  marker[MAXMARK];		/* array of markers */
int	nmarkers = 0;		/* number of markers */

int	expand;			/* expand commands? */

char	*progname;

static void convert(char *name, FILE *fin);
static void cvcomm(char *fname, FILE *fin);
static void cvobject(char *fname, FILE *fin);
static void replace(char *fname, struct mrkr *m, char *mark, FILE *fin);
static int edgecmp(const void *e1, const void *e2);
static int buildxf(char *xf, double markscale, FILE *fin);
static int addrot(char *xf, FVECT xp, FVECT yp, FVECT zp);


int
main(
	int	argc,
	char	*argv[]
)
{
	FILE	*fp;
	int	i, j;

	progname = argv[0];
	i = 1;
	while (i < argc && argv[i][0] == '-') {
		do {
			switch (argv[i][1]) {
			case 'i':
				marker[nmarkers].doxform = 0;
				marker[nmarkers].objname = argv[++i];
				break;
			case 'x':
				marker[nmarkers].doxform = 1;
				marker[nmarkers].objname = argv[++i];
				break;
			case 'e':
				expand = 1;
				break;
			case 'm':
				marker[nmarkers].modout = argv[++i];
				break;
			case 's':
				marker[nmarkers].mscale = atof(argv[++i]);
				break;
			default:
				goto userr;
			}
			if (++i >= argc)
				goto userr;
		} while (argv[i][0] == '-');
		if (marker[nmarkers].objname == NULL)
			goto userr;
		marker[nmarkers++].modin = argv[i++];
		if (nmarkers >= MAXMARK)
			break;
		marker[nmarkers].mscale = marker[nmarkers-1].mscale;
	}
	if (nmarkers == 0)
		goto userr;
					/* simple header */
	putchar('#');
	for (j = 0; j < i; j++) {
		putchar(' ');
		fputs(argv[j], stdout);
	}
	putchar('\n');
	if (i == argc)
		convert("<stdin>", stdin);
	else
		for ( ; i < argc; i++) {
			if ((fp = fopen(argv[i], "r")) == NULL) {
				perror(argv[i]);
				exit(1);
			}
			convert(argv[i], fp);
			fclose(fp);
		}
	return 0;
userr:
	fprintf(stderr,
"Usage: %s [-e][-s size][-m modout] {-x objfile|-i octree} modname .. [file ..]\n",
		progname);
	return 1;
}


void
convert(		/* replace marks in a stream */
	char	*name,
	register FILE	*fin
)
{
	register int	c;

	while ((c = getc(fin)) != EOF) {
		if (isspace(c))				/* blank */
			continue;
		if (c == '#') {				/* comment */
			putchar(c);
			do {
				if ((c = getc(fin)) == EOF)
					return;
				putchar(c);
			} while (c != '\n');
		} else if (c == '!') {			/* command */
			ungetc(c, fin);
			cvcomm(name, fin);
		} else {				/* object */
			ungetc(c, fin);
			cvobject(name, fin);
		}
	}
}


void
cvcomm(		/* convert a command */
	char	*fname,
	FILE	*fin
)
{
	FILE	*pin;
	char	buf[512], *fgetline();

	fgetline(buf, sizeof(buf), fin);
	if (expand) {
		if ((pin = popen(buf+1, "r")) == NULL) {
			fprintf(stderr,
			"%s: (%s): cannot execute \"%s\"\n",
					progname, fname, buf);
			exit(1);
		}
		convert(buf, pin);
		pclose(pin);
	} else
		printf("\n%s\n", buf);
}


void
cvobject(		/* convert an object */
	char	*fname,
	FILE	*fin
)
{
	extern char	*fgetword();
	char	buf[128], typ[16], nam[128];
	int	i, n;
	register int	j;

	if (fscanf(fin, "%s %s %s", buf, typ, nam) != 3)
		goto readerr;
	if (!strcmp(typ, "polygon"))
		for (j = 0; j < nmarkers; j++)
			if (!strcmp(buf, marker[j].modin)) {
				replace(fname, &marker[j], nam, fin);
				return;
			}
	printf("\n%s %s %s\n", buf, typ, nam);
	if (!strcmp(typ, "alias")) {		/* alias special case */
		if (fscanf(fin, "%s", buf) != 1)
			goto readerr;
		printf("\t%s\n", buf);
		return;
	}
	for (i = 0; i < 3; i++) {		/* pass along arguments */
		if (fscanf(fin, "%d", &n) != 1)
			goto readerr;
		printf("%d", n);
		for (j = 0; j < n; j++) {
			if (fgetword(buf, sizeof(buf), fin) == NULL)
				goto readerr;
			if (j%3 == 0)
				putchar('\n');
			putchar('\t');
			fputword(buf, stdout);
		}
		putchar('\n');
	}
	return;
readerr:
	fprintf(stderr, "%s: (%s): read error for %s \"%s\"\n",
			progname, fname, typ, nam);
}


void
replace(		/* replace marker */
	char	*fname,
	register struct mrkr	*m,
	char	*mark,
	FILE	*fin
)
{
	int	n;
	char	buf[256];

	buf[0] = '\0';			/* bug fix thanks to schorsch */
	if (m->doxform) {
		sprintf(buf, "xform -e -n %s", mark);
		if (m->modout != NULL)
			sprintf(buf+strlen(buf), " -m %s", m->modout);
		if (buildxf(buf+strlen(buf), m->mscale, fin) < 0)
			goto badxf;
		sprintf(buf+strlen(buf), " %s", m->objname);
		if (expand) {
			fflush(stdout);
			system(buf);
		} else
			printf("\n!%s\n", buf);
	} else {
		if ((n = buildxf(buf, m->mscale, fin)) < 0)
			goto badxf;
		printf("\n%s instance %s\n",
				m->modout==NULL?"void":m->modout, mark);
		printf("%d %s%s\n0\n0\n", n+1, m->objname, buf);
	}
	return;
badxf:
	fprintf(stderr, "%s: (%s): bad arguments for marker \"%s\"\n",
			progname, fname, mark);
	exit(1);
}


int
edgecmp(			/* compare two edges, descending order */
	const void *e1,
	const void *e2
)
{
	if (((EDGE*)e1)->len2 > ((EDGE*)e2)->len2)
		return(-1);
	if (((EDGE*)e1)->len2 < ((EDGE*)e2)->len2)
		return(1);
	return(0);
}


int
buildxf(		/* build transform for marker */
	register char	*xf,
	double	markscale,
	FILE	*fin
)
{
	static FVECT	vlist[MAXVERT];
	static EDGE	elist[MAXVERT];
	FVECT	xvec, yvec, zvec;
	double	xlen;
	int	n;
	register int	i;
	/*
	 * Read and sort vectors:  longest is hypotenuse,
	 *	second longest is x' axis,
	 *	third longest is y' axis (approximate),
	 *	other vectors are ignored.
	 * It is an error if the x' and y' axes do
	 *	not share a common vertex (their origin).
	 */
	if (fscanf(fin, "%*d %*d %d", &n) != 1)
		return(-1);
	if (n%3 != 0)
		return(-1);
	n /= 3;
	if (n < 3 || n > MAXVERT)
		return(-1);
					/* sort edges in descending order */
	for (i = 0; i < n; i++) {
		if (fscanf(fin, "%lf %lf %lf", &vlist[i][0],
				&vlist[i][1], &vlist[i][2]) != 3)
			return(-1);
		if (i) {
			elist[i].beg = i-1;
			elist[i].end = i;
			elist[i].len2 = dist2(vlist[i-1],vlist[i]);
		}
	}
	elist[0].beg = n-1;
	elist[0].end = 0;
	elist[0].len2 = dist2(vlist[n-1],vlist[0]);
	qsort(elist, n, sizeof(EDGE), edgecmp);
					/* find x' and y' */
	if (elist[1].end == elist[2].beg || elist[1].end == elist[2].end) {
		i = elist[1].beg;
		elist[1].beg = elist[1].end;
		elist[1].end = i;
	}
	if (elist[2].end == elist[1].beg) {
		i = elist[2].beg;
		elist[2].beg = elist[2].end;
		elist[2].end = i;
	}
	if (elist[1].beg != elist[2].beg)
		return(-1);		/* x' and y' not connected! */
	for (i = 0; i < 3; i++) {
		xvec[i] = vlist[elist[1].end][i] - vlist[elist[1].beg][i];
		yvec[i] = vlist[elist[2].end][i] - vlist[elist[2].beg][i];
	}
	if ((xlen = normalize(xvec)) == 0.0)
		return(-1);
	fcross(zvec, xvec, yvec);
	if (normalize(zvec) == 0.0)
		return(-1);
	fcross(yvec, zvec, xvec);
	n = 0;				/* start transformation... */
	if (markscale > 0.0) {		/* add scale factor */
		sprintf(xf, " -s %f", xlen*markscale);
		n += 2;
		while (*xf) ++xf;
	}
					/* add rotation */
	n += addrot(xf, xvec, yvec, zvec);
	while (*xf) ++xf;
					/* add translation */
	n += 4;
	sprintf(xf, " -t %f %f %f", vlist[elist[1].beg][0],
			vlist[elist[1].beg][1], vlist[elist[1].beg][2]);
	return(n);			/* all done */
}


int
addrot(		/* compute rotation (x,y,z) => (xp,yp,zp) */
	register char	*xf,
	FVECT xp,
	FVECT yp,
	FVECT zp
)
{
	int	n;
	double	theta;

	n = 0;
	theta = atan2(yp[2], zp[2]);
	if (!FEQ(theta,0.0)) {
		sprintf(xf, " -rx %f", theta*(180./PI));
		while (*xf) ++xf;
		n += 2;
	}
	theta = asin(-xp[2]);
	if (!FEQ(theta,0.0)) {
		sprintf(xf, " -ry %f", theta*(180./PI));
		while (*xf) ++xf;
		n += 2;
	}
	theta = atan2(xp[1], xp[0]);
	if (!FEQ(theta,0.0)) {
		sprintf(xf, " -rz %f", theta*(180./PI));
		/* while (*xf) ++xf; */
		n += 2;
	}
	return(n);
}
