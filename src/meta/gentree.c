#ifndef lint
static const char	RCSid[] = "$Id: gentree.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  gentree.c - program to generate 2-D Christmas trees.
 *
 *     12/4/85
 *
 *  cc gentree.c libmeta.a
 */


#include  "meta.h"


#define  ratio  1/3


char  *progname;

int  nbranch;


main(argc, argv)
int  argc;
char  *argv[];
{
#ifdef  CPM
	fixargs("gentree", &argc, &argv);
#endif
	progname = argv[0];
	if (argc != 7) {
		sprintf(errmsg, "Usage: %s x0 y0 x1 y1 nbranch depth", progname);
		error(USER, errmsg);
	}
	nbranch = atoi(argv[5]);
	gentree(atoi(argv[1]), atoi(argv[2]),
			atoi(argv[3]), atoi(argv[4]), atoi(argv[6]));
	pglob(PEOP, 0200, NULL);
	writeof(stdout);
	exit(0);
}


gentree(x0, y0, x1, y1, depth)
int  x0, y0, x1, y1, depth;
{
	int  xstart, ystart;
	int  xend, yend;
	register int  i;
	
	plseg(0, x0, y0, x1, y1);
	if (depth <= 0)
		return;
	for (i = 1; i < nbranch+1; i++) {
		xstart = (x1-x0)*(long)i/(nbranch+2)+x0;
		ystart = (y1-y0)*(long)i/(nbranch+2)+y0;
		xend = xstart+(x1-x0)/(nbranch+2);
		yend = ystart+(y1-y0)/(nbranch+2);
		gentree(xstart, ystart,
			xend+(y1-yend)*ratio, yend+(xend-x1)*ratio,
			depth-1);
		gentree(xstart, ystart,
			xend+(yend-y1)*ratio, yend+(x1-xend)*ratio,
			depth-1);
	}
}
