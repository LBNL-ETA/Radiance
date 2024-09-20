#ifndef lint
static const char	RCSid[] = "$Id: pextrem.c,v 2.16 2024/09/20 17:39:12 greg Exp $";
#endif
/*
 * Find extrema points in a Radiance picture (RGBE, XYZE, or HyperSpectral)
 */

#include  <math.h>

#include  "rtio.h"
#include  "platform.h"
#include  "color.h"
#include  "resolu.h"


int  orig = 0;

COLOR  expos = WHTCOLOR;

char	fmt[MAXFMTLEN];

static gethfunc headline;


static int
headline(			/* check header line */
	char  *s,
	void	*p
)
{
	double	d;
	COLOR	ctmp;

	if (formatval(fmt, s))			/* format */
		return(0);
	if (iswlsplit(s)) {			/* wavelength splits */
		wlsplitval(WLPART, s);
		return(0);
	}
	if (isncomp(s)) {			/* # spectral components */
		NCSAMP = ncompval(s);
		return(0);
	}
	if (!orig)				/* don't undo exposure? */
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
	int  x;
	COLRV  *scan;
	SCOLR  cmin, cmax;
	COLOR  tcol;
	int  xmin, ymin, xmax, ymax;

	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	for (i = 1; i < argc; i++)	/* get options */
		if (!strcmp(argv[i], "-o"))
			orig = 1;
		else if (!strcmp(argv[i], "-O"))
			orig = -1;
		else
			break;

	if (i == argc-1 && freopen(argv[i], "r", stdin) == NULL) {
		fprintf(stderr, "%s: can't open input \"%s\"\n",
				argv[0], argv[i]);
		return(1);
	}
					/* get our header */
	if (getheader(stdin, headline, NULL) < 0 ||
			(!globmatch(PICFMT, fmt) && strcmp(fmt, SPECFMT)) ||
			fgetresolu(&xres, &yres, stdin) < 0) {
		fprintf(stderr, "%s: bad picture format\n", argv[0]);
		return(1);
	}
	if (setspectrsamp(CNDX, WLPART) < 0) {
		fprintf(stderr, "%s: bad wavelength split or component count",
				argv[0]);
		return(1);
	}
	if (orig < 0 && !strcmp(CIEFMT, fmt))
		scalecolor(expos, 1./WHTEFFICACY);
	if ((scan = (COLRV *)malloc(xres*sizeof(COLRV)*(NCSAMP+1))) == NULL) {
		fprintf(stderr, "%s: out of memory\n", argv[0]);
		return(1);
	}
	setscolr(cmin, 1e30, 1e30, 1e30); xmin=ymin=0;
	scolrblack(cmax); xmax=ymax=0;
					/* find extrema */
	for (y = yres-1; y >= 0; y--) {
		if (freadscolrs(scan, NCSAMP, xres, stdin) < 0) {
			fprintf(stderr, "%s: read error on input\n", argv[0]);
			return(1);
		}
		for (x = xres; x-- > 0; ) {
			const COLRV *	sclr = scan + x*(NCSAMP+1);
			if (sclr[CNDX[EXP]] > cmax[CNDX[EXP]] ||
					(sclr[CNDX[EXP]] == cmax[CNDX[EXP]] &&
					 normpbright(sclr) > normpbright(cmax))) {
				copyscolr(cmax, sclr);
				xmax = x; ymax = y;
			}
			if (sclr[CNDX[EXP]] < cmin[CNDX[EXP]] ||
				(sclr[CNDX[EXP]] == cmin[CNDX[EXP]] &&
					 normpbright(sclr) < normpbright(cmin))) {
				copyscolr(cmin, sclr);
				xmin = x; ymin = y;
			}
		}
	}
	free(scan);
	scolr_color(tcol, cmin);
	printf("%d %d\t%.2e %.2e %.2e\n", xmin, ymin,
			colval(tcol,RED)/colval(expos,RED),
			colval(tcol,GRN)/colval(expos,GRN),
			colval(tcol,BLU)/colval(expos,BLU));
	scolr_color(tcol, cmax);
	printf("%d %d\t%.2e %.2e %.2e\n", xmax, ymax,
			colval(tcol,RED)/colval(expos,RED),
			colval(tcol,GRN)/colval(expos,GRN),
			colval(tcol,BLU)/colval(expos,BLU));
	return(0);
}
