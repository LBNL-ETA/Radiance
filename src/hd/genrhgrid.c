#ifndef lint
static const char	RCSid[] = "$Id: genrhgrid.c,v 3.8 2004/01/01 11:21:55 schorsch Exp $";
#endif
/*
 * Generate renderable grids from a holodeck file
 */

#include <stdio.h>

#include "platform.h"
#include "resolu.h"
#include "holo.h"

char	*progname;		/* global argv[0] */

char	*mat, *name;		/* material and object id */
double	rad;			/* grid line radius */

static void gridsect(char *fname, int sect);
static void putgrid(HOLO *hp);
static void putline(FVECT wp[2]);


int
main(
	int	argc,
	char	*argv[]
)
{
	int	sect;

	progname = argv[0];
	if ((argc < 5) | (argc > 6))
		goto userr;
	mat = argv[1];
	name = argv[2];
	rad = atof(argv[3]);
	sect = argc==5 ? -1 : atoi(argv[5]);
	fputs("# ", stdout);
	printargs(argc, argv, stdout);
	gridsect(argv[4], sect);
	quit(0);
userr:
	fprintf(stderr, "Usage: %s mat name rad input.hdk [section]\n",
			progname);
	exit(1);
}


void
gridsect(		/* get specified section(s) and print grids */
	char	*fname,
	int	sect
)
{
	FILE	*fp;
	HOLO	hdsect;
	int	fd;
	int32	nextloc;
	int	n;
					/* open holodeck file */
	if ((fp = fopen(fname, "r")) == NULL) {
		sprintf(errmsg, "cannot open \"%s\"", fname);
		error(SYSTEM, errmsg);
	}
					/* check header and magic number */
	if (checkheader(fp, HOLOFMT, NULL) < 0 || getw(fp) != HOLOMAGIC) {
		sprintf(errmsg, "file \"%s\" not in holodeck format", fname);
		error(USER, errmsg);
	}
	fd = dup(fileno(fp));			/* dup file handle */
	nextloc = ftell(fp);			/* get stdio position */
	fclose(fp);				/* done with stdio */
	for (n = 0; nextloc > 0L; n++) {	/* get the section(s) */
		lseek(fd, (off_t)nextloc, SEEK_SET);
		read(fd, (char *)&nextloc, sizeof(nextloc));
		if ((sect < 0) | (n == sect)) {
			read(fd, (char *)&hdsect, sizeof(HDGRID));
			hdcompgrid(&hdsect);
			putgrid(&hdsect);	/* print grid */
		}
	}
}


void
putgrid(			/* run through holodeck section grid lines */
	register HOLO	*hp
)
{
	register int	w, i;
	int	g0, g1;
	FVECT	wp[2], mov;
	double	d;
					/* do each wall on this section */
	for (w = 0; w < 6; w++) {
		g0 = hdwg0[w];
		g1 = hdwg1[w];
		d = 1.0/hp->grid[g0];
		mov[0] = d * hp->xv[g0][0];
		mov[1] = d * hp->xv[g0][1];
		mov[2] = d * hp->xv[g0][2];
		if (w & 1) {
			VSUM(wp[0], hp->orig, hp->xv[w>>1], 1.);
			VSUM(wp[0], wp[0], mov, 1.);
		} else
			VCOPY(wp[0], hp->orig);
		VSUM(wp[1], wp[0], hp->xv[g1], 1.);
		for (i = hp->grid[g0]; ; ) {	/* g0 lines */
			putline(wp);
			if (!--i) break;
			wp[0][0] += mov[0]; wp[0][1] += mov[1];
			wp[0][2] += mov[2]; wp[1][0] += mov[0];
			wp[1][1] += mov[1]; wp[1][2] += mov[2];
		}
		d = 1.0/hp->grid[g1];
		mov[0] = d * hp->xv[g1][0];
		mov[1] = d * hp->xv[g1][1];
		mov[2] = d * hp->xv[g1][2];
		if (w & 1)
			VSUM(wp[0], hp->orig, hp->xv[w>>1], 1.);
		else
			VSUM(wp[0], hp->orig, mov, 1.);
		VSUM(wp[1], wp[0], hp->xv[g0], 1.);
		for (i = hp->grid[g1]; ; ) {	/* g1 lines */
			putline(wp);
			if (!--i) break;
			wp[0][0] += mov[0]; wp[0][1] += mov[1];
			wp[0][2] += mov[2]; wp[1][0] += mov[0];
			wp[1][1] += mov[1]; wp[1][2] += mov[2];
		}
	}
}


void
putline(		/* put out a line */
	register FVECT	wp[2]
)
{
	static int	cnt = 0;

	printf("\n%s cylinder %s.%d\n0\n0\n7\n", mat, name, ++cnt);
	printf("\t%.4e %.4e %.4e\n", wp[0][0], wp[0][1], wp[0][2]);
	printf("\t%.4e %.4e %.4e\n", wp[1][0], wp[1][1], wp[1][2]);
	printf("\t%.4e\n", rad);
}
