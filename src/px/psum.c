#ifndef lint
static const char	RCSid[] = "$Id: psum.c,v 2.5 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 *  psum.c - program to sum pictures.
 *
 *     10/14/86
 */

#include  <stdio.h>

#include  <math.h>

#include  "color.h"


#define  MAXFILE	8

int  xsiz, ysiz;

char  *progname;

char  *fname[MAXFILE];			/* the file names */
FILE  *fptr[MAXFILE];			/* the file pointers */
COLOR  scale[MAXFILE];			/* scaling factors */
int  nfile;				/* number of files */


int
tabputs(s)			/* print line preceded by a tab */
char  *s;
{
	putc('\t', stdout);
	return(fputs(s, stdout));
}


main(argc, argv)
int  argc;
char  *argv[];
{
	double  d;
	int  xres, yres;
	int  an;

	progname = argv[0];

	nfile = 0;
	setcolor(scale[0], 1.0, 1.0, 1.0);

	for (an = 1; an < argc; an++) {
		if (nfile >= MAXFILE) {
			fprintf(stderr, "%s: too many files\n", progname);
			quit(1);
		}
		if (argv[an][0] == '-')
			switch (argv[an][1]) {
			case 's':
				d = atof(argv[an+1]);
				switch (argv[an][2]) {
				case '\0':
					scalecolor(scale[nfile], d);
					break;
				case 'r':
					colval(scale[nfile],RED) *= d;
					break;
				case 'g':
					colval(scale[nfile],GRN) *= d;
					break;
				case 'b':
					colval(scale[nfile],BLU) *= d;
					break;
				default:
					goto unkopt;
				}
				an++;
				continue;
			case '\0':
				fptr[nfile] = stdin;
				fname[nfile] = "<stdin>";
				break;
			default:;
			unkopt:
				fprintf(stderr, "%s: unknown option: %s\n",
						progname, argv[an]);
				quit(1);
			}
		else if ((fptr[nfile] = fopen(argv[an], "r")) == NULL) {
			fprintf(stderr, "%s: can't open file: %s\n",
						progname, argv[an]);
			quit(1);
		} else
			fname[nfile] = argv[an];
						/* get header */
		fputs(fname[nfile], stdout);
		fputs(":\n", stdout);
		getheader(fptr[nfile], tabputs);
						/* get picture size */
		if (fgetresolu(&xres, &yres, fptr[nfile]) != (YMAJOR|YDECR)) {
			fprintf(stderr, "%s: bad picture size\n", progname);
			quit(1);
		} else if (nfile == 0) {
			xsiz = xres;
			ysiz = yres;
		} else if (xres != xsiz || yres != ysiz) {
			fprintf(stderr, "%s: pictures are different sizes\n",
					progname);
			quit(1);
		}
		nfile++;
		setcolor(scale[nfile], 1.0, 1.0, 1.0);
	}
					/* add new header info. */
	printargs(argc, argv, stdout);
	putchar('\n');
	fputresolu(YMAJOR|YDECR, xsiz, ysiz, stdout);

	psum();
	
	quit(0);
}


psum()				/* sum the files */
{
	COLOR  *scanin, *scanout;
	int  y, i;
	register int  x;

	scanin = (COLOR *)malloc(xsiz*sizeof(COLOR));
	scanout = (COLOR *)malloc(xsiz*sizeof(COLOR));
	if (scanin == NULL || scanout == NULL) {
		fprintf(stderr, "%s: out of memory\n", progname);
		quit(1);
	}
	for (y = ysiz-1; y >= 0; y--) {
		for (x = 0; x < xsiz; x++)
			setcolor(scanout[x], 0.0, 0.0, 0.0);
		for (i = 0; i < nfile; i++) {
			if (freadscan(scanin, xsiz, fptr[i]) < 0) {
				fprintf(stderr, "%s: read error on file: %s\n",
						progname, fname[i]);
				quit(1);
			}
			for (x = 0; x < xsiz; x++)
				multcolor(scanin[x], scale[i]);
			for (x = 0; x < xsiz; x++)
				addcolor(scanout[x], scanin[x]);
		}
		if (fwritescan(scanout, xsiz, stdout) < 0) {
			fprintf(stderr, "%s: write error\n", progname);
			quit(1);
		}
	}
	free((void *)scanin);
	free((void *)scanout);
}


void
quit(code)		/* exit gracefully */
int  code;
{
	exit(code);
}
