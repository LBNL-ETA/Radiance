#ifndef lint
static const char	RCSid[] = "$Id: vwright.c,v 2.4 2003/02/22 02:07:30 greg Exp $";
#endif
/*
 * Move a viewpoint the given distance to the right
 */


#include "standard.h"

#include "view.h"

#include <ctype.h>

VIEW	vw = STDVIEW;

char	*progname;


main(argc, argv)
int	argc;
char	*argv[];
{
	char	linebuf[256];
	char	*err;
	int	gotview = 0;
	FVECT	hvn, vvn;
	double	dist;
	register int	i;

	progname = argv[0];
	if (argc != 2)
		goto userr;
	while (fgets(linebuf, sizeof(linebuf), stdin) != NULL) {
		if (linebuf[0] == '\n')
			break;
		if (isview(linebuf) && sscanview(&vw, linebuf) > 0)
			gotview++;
	}
	if (!gotview) {
		fprintf(stderr, "%s: no view on standard input\n", progname);
		exit(1);
	}
	if ((err= setview(&vw)) != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	VCOPY(hvn, vw.hvec);
	normalize(hvn);
	VCOPY(vvn, vw.vvec);
	normalize(vvn);
	if (isalpha(argv[1][0])) {	/* print out variables */
		switch (vw.type) {
			case VT_PER: i=1; break;
			case VT_PAR: i=2; break;
			case VT_ANG: i=3; break;
			case VT_HEM: i=4; break;
			case VT_CYL: i=5; break;
			default: i=0; break;
		}
		printf("%st:%d;", argv[1], i);
		printf("%spx:%g;%spy:%g;%spz:%g;", argv[1], vw.vp[0],
				argv[1], vw.vp[1], argv[1], vw.vp[2]);
		printf("%sdx:%g;%sdy:%g;%sdz:%g;", argv[1], vw.vdir[0],
				argv[1], vw.vdir[1], argv[1], vw.vdir[2]);
		printf("%sux:%g;%suy:%g;%suz:%g;", argv[1], vw.vup[0],
				argv[1], vw.vup[1], argv[1], vw.vup[2]);
		printf("%sh:%g;%sv:%g;", argv[1], vw.horiz,
				argv[1], vw.vert);
		printf("%ss:%g;%sl:%g;%so:%g;%sa:%g;",
				argv[1], vw.hoff, argv[1], vw.voff,
				argv[1], vw.vfore, argv[1], vw.vaft);
		printf("%shx:%g;%shy:%g;%shz:%g;%shn:%g;",
				argv[1], hvn[0], argv[1], hvn[1],
				argv[1], hvn[2], argv[1], sqrt(vw.hn2));
		printf("%svx:%g;%svy:%g;%svz:%g;%svn:%g;\n",
				argv[1], vvn[0], argv[1], vvn[1],
				argv[1], vvn[2], argv[1], sqrt(vw.vn2));
		exit(0);
	}
	if (!isflt(argv[1]))
		goto userr;
	dist = atof(argv[1]);
	for (i = 0; i < 3; i++)
		vw.vp[i] += dist*hvn[i];
	fprintview(&vw, stdout);
	putchar('\n');
	exit(0);
userr:
	fprintf(stderr, "Usage: %s {offset|name}\n", progname);
	exit(1);
}
