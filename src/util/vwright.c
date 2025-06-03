#ifndef lint
static const char	RCSid[] = "$Id: vwright.c,v 2.11 2025/06/03 21:31:51 greg Exp $";
#endif
/*
 * Move a viewpoint the given distance to the right
 */


#include "standard.h"

#include "view.h"

#include <ctype.h>

VIEW	vw = STDVIEW;


int
main(
	int	argc,
	char	*argv[]
)
{
	char	*err;
	int	gotview = 0;
	FVECT	hvn, vvn;
	double	dist;
	int	i;

	fixargv0(argv[0]);		/* sets global progname */
	++argv; --argc;
	while (argc && argv[0][0] == '-' && argv[0][1] == 'v') {
		int	rv;
		if (argc > 2 && argv[0][2] == 'f') {
			rv = viewfile(argv[1], &vw, NULL);
			if (rv <= 0) {
				fprintf(stderr, "%s: %s file '%s'\n",
						progname,
						rv < 0 ? "cannot open view" :
						"no view in", argv[1]);
				exit(1);
			}
			++gotview;
			++argv; --argc;
		} else {
			rv = getviewopt(&vw, argc, argv);
			if (rv < 0) {
				fprintf(stderr, "%s: bad view option at '%s'\n",
						progname, argv[0]);
				exit(1);
			}
			++gotview;
			argv += rv; argc -= rv;
		}
		++argv; --argc;
	}
	if (argc != 1)
		goto userr;
	if (!gotview) {
		gotview = viewfile(NULL, &vw, NULL);
		if (gotview <= 0) {
			fprintf(stderr, "%s: no view on standard input\n",
						progname);
			exit(1);
		}
	}
	if ((err= setview(&vw)) != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	VCOPY(hvn, vw.hvec);
	normalize(hvn);
	VCOPY(vvn, vw.vvec);
	normalize(vvn);
	if (isalpha(argv[0][0])) {	/* print out variables */
		switch (vw.type) {
			case VT_PER: i=1; break;
			case VT_PAR: i=2; break;
			case VT_ANG: i=3; break;
			case VT_HEM: i=4; break;
			case VT_CYL: i=5; break;
			case VT_PLS: i=6; break;
			default: i=0; break;
		}
		printf("%st:%d;", argv[0], i);
		printf("%spx:%g;%spy:%g;%spz:%g;", argv[0], vw.vp[0],
				argv[0], vw.vp[1], argv[0], vw.vp[2]);
		printf("%sdx:%g;%sdy:%g;%sdz:%g;", argv[0], vw.vdir[0],
				argv[0], vw.vdir[1], argv[0], vw.vdir[2]);
		printf("%sd:%g;", argv[0], vw.vdist);
		printf("%sux:%g;%suy:%g;%suz:%g;", argv[0], vw.vup[0],
				argv[0], vw.vup[1], argv[0], vw.vup[2]);
		printf("%sh:%g;%sv:%g;", argv[0], vw.horiz,
				argv[0], vw.vert);
		printf("%ss:%g;%sl:%g;%so:%g;%sa:%g;",
				argv[0], vw.hoff, argv[0], vw.voff,
				argv[0], vw.vfore, argv[0], vw.vaft);
		printf("%shx:%g;%shy:%g;%shz:%g;%shn:%g;",
				argv[0], hvn[0], argv[0], hvn[1],
				argv[0], hvn[2], argv[0], sqrt(vw.hn2));
		printf("%svx:%g;%svy:%g;%svz:%g;%svn:%g;\n",
				argv[0], vvn[0], argv[0], vvn[1],
				argv[0], vvn[2], argv[0], sqrt(vw.vn2));
		exit(0);
	}
	if (!isflt(argv[0]))
		goto userr;
	dist = atof(argv[0]);
	VSUM(vw.vp, vw.vp, hvn, dist);
	fprintview(&vw, stdout);
	fputc('\n', stdout);
	exit(0);
userr:
	fprintf(stderr, "Usage: %s [view options] {offset|name}\n", progname);
	exit(1);
}
