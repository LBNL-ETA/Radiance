#ifndef lint
static const char	RCSid[] = "$Id: protate.c,v 2.10 2004/04/30 20:29:24 greg Exp $";
#endif
/*
 * prot.c - program to rotate picture file 90 degrees clockwise.
 *
 *	2/26/88
 */

#include "standard.h"

#include "color.h"

#include <time.h>

#include "resolu.h"

int	order;				/* input scanline order */
int	xres, yres;			/* input resolution */

int	correctorder = 0;		/* order correction? */
int	ccw = 0;			/* rotate CCW? */

#ifdef SMLMEM
char	buf[1<<20];			/* output buffer */
#else
char	buf[1<<22];			/* output buffer */
#endif

int	nrows;				/* number of rows output at once */

#define scanbar		((COLR *)buf)

char	*progname;

short	ordertab[4][2] = {
	{0,0}, {XDECR,XDECR|YDECR}, {XDECR|YDECR,YDECR}, {YDECR,XDECR}
};


static int neworder(void);
static void rotatecw(FILE *fp);
static void rotateccw(FILE *fp);



static int
neworder(void)		/* return corrected order */
{
	static short    ordercw[8];
	register int	i;

	if (correctorder)
		return(order);
	if (!ordercw[0]) {
		ordercw[YMAJOR|YDECR] = 0;
		ordercw[0] = YMAJOR|XDECR;
		ordercw[YMAJOR|XDECR] = XDECR|YDECR;
		ordercw[XDECR|YDECR] = YMAJOR|YDECR;
		ordercw[YMAJOR|XDECR|YDECR] = XDECR;
		ordercw[XDECR] = YMAJOR;
		ordercw[YMAJOR] = YDECR;
		ordercw[YDECR] = YMAJOR|XDECR|YDECR;
	}
	if (!ccw)
		return(ordercw[order]);
	for (i = 8; i--; )
		if (ordercw[i] == order)
			return(i);
	fputs("Order botch!\n", stderr);
	exit(2);
}

int
main(
	int	argc,
	char	*argv[]
)
{
	static char	picfmt[LPICFMT+1] = PICFMT;
	int	rval;
	FILE	*fin;

	progname = argv[0];

	while (argc > 2 && argv[1][0] == '-') {
		switch (argv[1][1]) {
		case 'c':
			correctorder = 1;
			break;
		case 'r':
			ccw = 1;
			break;
		default:
			goto userr;
		}
		argc--; argv++;
	}
	if (argc != 2 && argc != 3)
		goto userr;
	if ((fin = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "%s: cannot open\n", argv[1]);
		exit(1);
	}
	if (argc == 3 && freopen(argv[2], "w", stdout) == NULL) {
		fprintf(stderr, "%s: cannot open\n", argv[2]);
		exit(1);
	}
					/* transfer header */
	if ((rval = checkheader(fin, picfmt, stdout)) < 0) {
		fprintf(stderr, "%s: not a Radiance picture\n", progname);
		exit(1);
	}
	if (rval)
		fputformat(picfmt, stdout);
					/* add new header info. */
	fputs(progname, stdout);
	if (ccw) fputs(" -r", stdout);
	if (correctorder) fputs(" -c", stdout);
	fputs("\n\n", stdout);
					/* get picture size */
	if ((order = fgetresolu(&xres, &yres, fin)) < 0) {
		fprintf(stderr, "%s: bad picture size\n", progname);
		exit(1);
	}
					/* write new picture size */
	fputresolu(neworder(), yres, xres, stdout);
					/* compute buffer capacity */
	nrows = sizeof(buf)/sizeof(COLR)/yres;
	if (ccw)			/* rotate the image */
		rotateccw(fin);
	else
		rotatecw(fin);
	exit(0);
userr:
	fprintf(stderr, "Usage: %s [-r][-c] infile [outfile]\n", progname);
	exit(1);
}


static void
rotatecw(			/* rotate picture clockwise */
	FILE	*fp
)
{
	register COLR	*inln;
	register int	xoff, inx, iny;
	long	start, ftell();

	if ((inln = (COLR *)malloc(xres*sizeof(COLR))) == NULL) {
		fprintf(stderr, "%s: out of memory\n", progname);
		exit(1);
	}
	start = ftell(fp);
	for (xoff = 0; xoff < xres; xoff += nrows) {
		if (fseek(fp, start, 0) < 0) {
			fprintf(stderr, "%s: seek error\n", progname);
			exit(1);
		}
		for (iny = yres-1; iny >= 0; iny--) {
			if (freadcolrs(inln, xres, fp) < 0) {
				fprintf(stderr, "%s: read error\n", progname);
				exit(1);
			}
			for (inx = 0; inx < nrows && xoff+inx < xres; inx++)
				copycolr(scanbar[inx*yres+iny],
						inln[xoff+inx]);
		}
		for (inx = 0; inx < nrows && xoff+inx < xres; inx++)
			if (fwritecolrs(scanbar+inx*yres, yres, stdout) < 0) {
				fprintf(stderr, "%s: write error\n", progname);
				exit(1);
			}
	}
	free((void *)inln);
}


static void
rotateccw(			/* rotate picture counter-clockwise */
	FILE	*fp
)
{
	register COLR	*inln;
	register int	xoff, inx, iny;
	long	start, ftell();

	if ((inln = (COLR *)malloc(xres*sizeof(COLR))) == NULL) {
		fprintf(stderr, "%s: out of memory\n", progname);
		exit(1);
	}
	start = ftell(fp);
	for (xoff = xres-1; xoff >= 0; xoff -= nrows) {
		if (fseek(fp, start, 0) < 0) {
			fprintf(stderr, "%s: seek error\n", progname);
			exit(1);
		}
		for (iny = 0; iny < yres; iny++) {
			if (freadcolrs(inln, xres, fp) < 0) {
				fprintf(stderr, "%s: read error\n", progname);
				exit(1);
			}
			for (inx = 0; inx < nrows && xoff-inx >= 0; inx++)
				copycolr(scanbar[inx*yres+iny],
						inln[xoff-inx]);
		}
		for (inx = 0; inx < nrows && xoff-inx >= 0; inx++)
			if (fwritecolrs(scanbar+inx*yres, yres, stdout) < 0) {
				fprintf(stderr, "%s: write error\n", progname);
				exit(1);
			}
	}
	free((void *)inln);
}
