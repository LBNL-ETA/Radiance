#ifndef lint
static const char	RCSid[] = "$Id: ra_rgbe.c,v 2.13 2003/07/03 22:41:44 schorsch Exp $";
#endif
/*
 *  program to convert from RADIANCE RLE to flat format
 */

#include  <stdio.h>
#include  <math.h>
#include  <time.h>
#include  <string.h>

#include  "platform.h"
#include  "color.h"
#include  "resolu.h"

extern int	addhline();

#define dumpheader(fp)	fwrite(headlines, 1, headlen, fp)

int  bradj = 0;				/* brightness adjustment */

int  doflat = 1;			/* produce flat file */

int  force = 0;				/* force file overwrite? */

int  findframe = 0;			/* find a specific frame? */

int  frameno = 0;			/* current frame number */
int  fmterr = 0;			/* got input format error */
char  *headlines;			/* current header info. */
int  headlen;				/* current header length */

char  *progname;


main(argc, argv)
int  argc;
char  *argv[];
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


transfer(ospec)			/* transfer a Radiance picture */
char	*ospec;
{
	char	oname[128];
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
	scanin = (COLR *)tempbuffer(xmax*sizeof(COLR));
	if (scanin == NULL) {
		perror(progname);
		exit(1);
	}
					/* skip frame? */
	if (findframe > frameno) {
		for (y = ymax; y--; )
			if (freadcolrs(scanin, xmax, stdin) < 0) {
				fprintf(stderr,
					"%s: error reading input picture\n",
						progname);
				exit(1);
			}
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
	dumpheader(fp);			/* put out header */
	fputs(progname, fp);
	if (bradj)
		fprintf(fp, " -e %+d", bradj);
	if (!doflat)
		fputs(" -r", fp);
	fputc('\n', fp);
	if (bradj)
		fputexpos(pow(2.0, (double)bradj), fp);
	fputc('\n', fp);
	fputresolu(order, xmax, ymax, fp);
					/* transfer picture */
	for (y = ymax; y--; ) {
		if (freadcolrs(scanin, xmax, stdin) < 0) {
			fprintf(stderr, "%s: error reading input picture\n",
					progname);
			exit(1);
		}
		if (bradj)
			shiftcolrs(scanin, xmax, bradj);
		if (doflat)
			fwrite((char *)scanin, sizeof(COLR), xmax, fp);
		else
			fwritecolrs(scanin, xmax, fp);
		if (ferror(fp)) {
			fprintf(stderr, "%s: error writing output to \"%s\"\n",
					progname, oname);
			exit(1);
		}
	}
					/* clean up */
	if (oname[0] == '!')
		pclose(fp);
	else if (ospec != NULL)
		fclose(fp);
	return(1);
}


int
addhline(s)			/* add a line to our info. header */
char	*s;
{
	char	fmt[32];
	int	n;

	if (formatval(fmt, s))
		fmterr += !globmatch(PICFMT, fmt);
	else if (!strncmp(s, "FRAME=", 6))
		frameno = atoi(s+6);
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


loadheader(fp)			/* load an info. header into memory */
FILE	*fp;
{
	fmterr = 0; frameno = 0;
	if (headlen) {			/* free old header */
		free(headlines);
		headlen = 0;
	}
	if (getheader(fp, addhline, NULL) < 0)
		return(0);
	if (fmterr)
		return(-1);
	return(1);
}
