#ifndef lint
static const char	RCSid[] = "$Id: ra_rgbe.c,v 2.23 2024/09/10 20:24:42 greg Exp $";
#endif
/*
 *  program to convert from RADIANCE RLE to flat format
 */

#include  <math.h>

#include  "platform.h"
#include  "rtio.h"
#include  "paths.h"
#include  "color.h"
#include  "resolu.h"

#define dumpheader(fp)	putbinary(headlines, 1, headlen, fp)

int  bradj = 0;				/* brightness adjustment */
int  doflat = 1;			/* produce flat file */
int  force = 0;				/* force file overwrite? */
int  findframe = 0;			/* find a specific frame? */
int  frameno = 0;			/* current frame number */
int  fmterr = 0;			/* got input format error */
char  *headlines = NULL;		/* current header info. */
int  headlen1 = 0;			/* length of initial frame header */
int  headlen = 0;			/* current header length */
char  fmt[MAXFMTLEN];			/* input format */

char  *progname;

static gethfunc addhline;
static int transfer(char *ospec);
static int loadheader(FILE *fp);


int
main(int  argc, char  *argv[])
{
	char	*ospec;
	int  i;

	progname = argv[0];

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'r':
				doflat = !doflat;
				break;
			case 'e':
				if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
					goto userr;
				bradj = atoi(argv[++i]);
				break;
			case 'n':
				findframe = atoi(argv[++i]);
				break;
			case 'f':
				force++;
				break;
			case '\0':
				goto gotfile;
			default:
				goto userr;
			}
		else
			break;
gotfile:
	if (i < argc-2)
		goto userr;
	if (i <= argc-1 && strcmp(argv[i], "-") &&
			freopen(argv[i], "r", stdin) == NULL) {
		fprintf(stderr, "%s: can't open input \"%s\"\n",
				progname, argv[i]);
		exit(1);
	}
	SET_FILE_BINARY(stdin);
	ospec = i==argc-2 ? argv[i+1] : (char *)NULL;
	while (transfer(ospec))
		;
	exit(0);
userr:
	fprintf(stderr,
	"Usage: %s [-r][-e +/-stops][-f][-n frame] [input [outspec]]\n",
			progname);
	exit(1);
}


static int
transfer(			/* transfer a Radiance picture */
	char	*ospec
)
{
	char	oname[PATH_MAX];
	FILE	*fp;
	int	order;
	int	xmax, ymax;
	COLR	*scanin;
	int	y;
					/* get header info. */
	if (!(y = loadheader(stdin)))
		return(0);
	if (y < 0 || (order = fgetresolu(&xmax, &ymax, stdin)) < 0) {
		fprintf(stderr, "%s: bad input format\n", progname);
		exit(1);
	}
					/* did we pass the target frame? */
	if (findframe && findframe < frameno)
		return(0);
					/* allocate scanline */
	scanin = (COLR *)malloc(xmax*sizeof(COLR));
	if (scanin == NULL) {
		perror(progname);
		exit(1);
	}
					/* skip frame? */
	if (findframe > frameno) {
		if (NCSAMP > 3) {
			if (fseek(stdin, ymax*xmax*LSCOLR, SEEK_CUR) < 0) {
				perror(progname);
				exit(1);
			}
		} else
		    for (y = ymax; y--; )
			if (freadcolrs(scanin, xmax, stdin) < 0) {
				fprintf(stderr,
					"%s: error reading input picture\n",
						progname);
				exit(1);
			}
		free(scanin);
		return(1);
	}
					/* open output file/command */
	if (ospec == NULL) {
		strcpy(oname, "<stdout>");
		fp = stdout;
	} else {
		sprintf(oname, ospec, frameno);
		if (oname[0] == '!') {
			if ((fp = popen(oname+1, "w")) == NULL) {
				fprintf(stderr, "%s: cannot start \"%s\"\n",
						progname, oname);
				exit(1);
			}
		} else {
			if (!force && access(oname, 0) >= 0) {
				fprintf(stderr,
					"%s: output file \"%s\" exists\n",
						progname, oname);
				exit(1);
			}
			if ((fp = fopen(oname, "w")) == NULL) {
				fprintf(stderr, "%s: ", progname);
				perror(oname);
				exit(1);
			}
		}
	}
	SET_FILE_BINARY(fp);
	newheader("RADIANCE", fp);		/* put out header */
	dumpheader(fp);
	fputs(progname, fp);
	if (bradj)
		fprintf(fp, " -e %+d", bradj);
	if (!doflat)
		fputs(" -r", fp);
	fputc('\n', fp);
	if (bradj)
		fputexpos(pow(2.0, (double)bradj), fp);
	if (frameno)
		fprintf(fp, "FRAME=%d\n", frameno);
	if (fmt[0])
		fputformat(fmt, fp);
	fputc('\n', fp);
	fputresolu(order, xmax, ymax, fp);
					/* transfer picture */
	for (y = ymax; y--; ) {
		if (fread2colrs(scanin, xmax, stdin, NCSAMP, WLPART) < 0) {
			fprintf(stderr, "%s: error reading input picture\n",
					progname);
			exit(1);
		}
		if (bradj)
			shiftcolrs(scanin, xmax, bradj);
		if (doflat ? (putbinary(scanin, sizeof(COLR), xmax, fp) != xmax) :
				(fwritecolrs(scanin, xmax, fp) < 0))
			goto writerr;
	}
	free(scanin);			/* clean up */
	if (fflush(fp) == EOF)
		goto writerr;
	if (oname[0] == '!')
		pclose(fp);
	else if (ospec != NULL)
		fclose(fp);
	return(1);
writerr:
	fprintf(stderr, "%s: error writing output to \"%s\"\n",
			progname, oname);
	exit(1);
}


static int
addhline(			/* add a line to our info. header */
	char	*s,
	void	*p
)
{
	int	n;

	if (isheadid(s))
		return(0);
	if (!strncmp(s, "FRAME=", 6)) {
		frameno = atoi(s+6);
		return(0);
	}
	if (formatval(fmt, s)) {
		if (!strcmp(fmt, SPECFMT))
			strcpy(fmt, COLRFMT);
		else
			fmterr += !globmatch(PICFMT, fmt);
		return(0);
	}
	if (isncomp(s)) {
		NCSAMP = ncompval(s);
		return(NCSAMP - 3);
	}
	if (iswlsplit(s)) {
		wlsplitval(WLPART, s);
		return(0);
	}
	n = strlen(s);
	if (headlen)
		headlines = (char *)realloc((void *)headlines, headlen+n+1);
	else
		headlines = (char *)malloc(n+1);
	if (headlines == NULL) {
		perror(progname);
		exit(1);
	}
	strcpy(headlines+headlen, s);
	headlen += n;
	return(0);
}


static int
loadheader(			/* load an info. header into memory */
	FILE	*fp
)
{
	fmterr = 0; frameno = 0;
				/* revert to initial header length */
	if (!headlen1) headlen1 = headlen;
	else headlen = headlen1;

	if (getheader(fp, addhline, NULL) < 0)
		return(0);
	if (fmterr)
		return(-1);
	return(1);
}
