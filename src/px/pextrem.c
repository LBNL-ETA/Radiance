#ifndef lint
static const char	RCSid[] = "$Id: pextrem.c,v 2.11 2005/06/02 18:13:47 greg Exp $";
#endif
/*
 * Find extrema points in a Radiance picture.
 */

#include  <stdio.h>
#include  <math.h>
#include  <string.h>

#include  "platform.h"
#include  "color.h"
#include  "resolu.h"


int  orig = 0;

int  wrongformat = 0;

COLOR  expos = WHTCOLOR;

static gethfunc headline;


static int
headline(			/* check header line */
	char  *s,
	void	*p
)
{
	char	fmt[32];
	double	d;
	COLOR	ctmp;

	if (isformat(s)) {			/* format */
		formatval(fmt, s);
		wrongformat = !globmatch(PICFMT, fmt);
	}
	if (!orig)
		return(0);
	if (isexpos(s)) {			/* exposure */
		d = exposval(s);
		scalecolor(expos, d);
	} else if (iscolcor(s)) {		/* color correction */
		colcorval(ctmp, s);
		multcolor(expos, ctmp);
	}
	return(0);
}


int
main(
	int  argc,
	char  *argv[]
)
{
	int  i;
	int  xres, yres;
	int  y;
	register int  x;
	COLR  *scan;
	COLR  cmin, cmax;
	int  xmin, ymin, xmax, ymax;
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	for (i = 1; i < argc; i++)	/* get options */
		if (!strcmp(argv[i], "-o"))
			orig++;
		else
			break;

	if (i == argc-1 && freopen(argv[i], "r", stdin) == NULL) {
		fprintf(stderr, "%s: can't open input \"%s\"\n",
				argv[0], argv[i]);
		exit(1);
	}
					/* get our header */
	if (getheader(stdin, headline, NULL) < 0 || wrongformat ||
			fgetresolu(&xres, &yres, stdin) < 0) {
		fprintf(stderr, "%s: bad picture format\n", argv[0]);
		exit(1);
	}
	if ((scan = (COLR *)malloc(xres*sizeof(COLR))) == NULL) {
		fprintf(stderr, "%s: out of memory\n", argv[0]);
		exit(1);
	}
	setcolr(cmin, 1e30, 1e30, 1e30);
	setcolr(cmax, 0., 0., 0.); xmax=ymax=0;
					/* find extrema */
	for (y = yres-1; y >= 0; y--) {
		if (freadcolrs(scan, xres, stdin) < 0) {
			fprintf(stderr, "%s: read error on input\n", argv[0]);
			exit(1);
		}
		for (x = xres; x-- > 0; ) {
			if (scan[x][EXP] > cmax[EXP] ||
					(scan[x][EXP] == cmax[EXP] &&
					 normbright(scan[x]) >
						normbright(cmax))) {
				copycolr(cmax, scan[x]);
				xmax = x; ymax = y;
			}
			if (scan[x][EXP] < cmin[EXP] ||
					(scan[x][EXP] == cmin[EXP] &&
					 normbright(scan[x]) <
						normbright(cmin))) {
				copycolr(cmin, scan[x]);
				xmin = x; ymin = y;
			}
		}
	}
	free((void *)scan);
	printf("%d %d\t%e %e %e\n", xmin, ymin,
			colrval(cmin,RED)/colval(expos,RED),
			colrval(cmin,GRN)/colval(expos,GRN),
			colrval(cmin,BLU)/colval(expos,BLU));
	printf("%d %d\t%e %e %e\n", xmax, ymax,
			colrval(cmax,RED)/colval(expos,RED),
			colrval(cmax,GRN)/colval(expos,GRN),
			colrval(cmax,BLU)/colval(expos,BLU));
	exit(0);
}
