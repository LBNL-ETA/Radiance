#ifndef lint
static const char	RCSid[] = "$Id: mgf2meta.c,v 2.12 2011/02/22 16:45:12 greg Exp $";
#endif
/*
 * Convert MGF (Materials and Geometry Format) to Metafile 2-d graphics
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "meta.h"
#include "random.h"
#include "mgf_parser.h"
#include "plocate.h" /* XXX shouldn't this rather be in rtmath.h? */

#define MSIZE	((1<<14)-1)
#define	MX(v)	(int)(MSIZE*(v)[(proj_axis+1)%3])
#define	MY(v)	(int)(MSIZE*(v)[(proj_axis+2)%3])

int	proj_axis;
double	limit[3][2];
int	layer;
long	rthresh = 1;

extern int	mg_nqcdivs;

static int r_face(int ac, char **av);
static void newlayer(void);
static int doline(int v1x, int v1y, int v2x, int v2y);


int
main(		/* convert files to stdout */
	int	argc,
	char	*argv[]
)
{
	int	i;
				/* initialize dispatch table */
	mg_ehand[MG_E_FACE] = r_face;
	mg_ehand[MG_E_POINT] = c_hvertex;
	mg_ehand[MG_E_VERTEX] = c_hvertex;
	mg_ehand[MG_E_XF] = xf_handler;
	mg_nqcdivs = 3;		/* reduce object subdivision */
	mg_init();		/* initialize the parser */
					/* get arguments */
	if (argc > 9 && !strcmp(argv[1], "-t")) {
		rthresh = atof(argv[2])*MSIZE + 0.5;
		rthresh *= rthresh;
		argv += 2;
		argc -= 2;
	}
	if (argc < 8 || (proj_axis = argv[1][0]-'x') < 0 || proj_axis > 2)
		goto userr;
	limit[0][0] = atof(argv[2]); limit[0][1] = atof(argv[3]);
	limit[1][0] = atof(argv[4]); limit[1][1] = atof(argv[5]);
	limit[2][0] = atof(argv[6]); limit[2][1] = atof(argv[7]);
	
	if (argc == 8) {		/* convert stdin */
		if (mg_load(NULL) != MG_OK)
			exit(1);
	} else				/* convert each file */
		for (i = 8; i < argc; i++) {
			if (mg_load(argv[i]) != MG_OK)
				exit(1);
			newlayer();
		}
	mendpage();			/* print page */
	mdone();			/* close output */
	exit(0);
userr:
	fputs("Usage: mgf2meta [-t thresh] {x|y|z} xmin xmax ymin ymax zmin zmax [file.mgf] ..\n",
			stderr);
	exit(1);
}


int
r_face(			/* convert a face */
	int	ac,
	char	**av
)
{
	static FVECT	bbmin = {0,0,0}, bbmax = {1,1,1};
	register int	i, j;
	register C_VERTEX	*cv;
	FVECT	v1, v2, vo;

	if (ac < 4)
		return(MG_EARGC);
						/* connect to last point */
	if ((cv = c_getvert(av[ac-1])) == NULL)
		return(MG_EUNDEF);
	xf_xfmpoint(vo, cv->p);
	for (j = 0; j < 3; j++)
		vo[j] = (vo[j] - limit[j][0])/(limit[j][1]-limit[j][0]);
	for (i = 1; i < ac; i++) {		/* go around face */
		if ((cv = c_getvert(av[i])) == NULL)
			return(MG_EUNDEF);
		xf_xfmpoint(v2, cv->p);
		for (j = 0; j < 3; j++)
			v2[j] = (v2[j] - limit[j][0])/(limit[j][1]-limit[j][0]);
		VCOPY(v1, vo);
		VCOPY(vo, v2);
		if (clip(v1, v2, bbmin, bbmax))
			doline(MX(v1), MY(v1), MX(v2), MY(v2));
	}
	return(MG_OK);
}


#define  HTBLSIZ	16381		/* prime hash table size */

short	hshtab[HTBLSIZ][4];		/* done line segments */

#define  hash(mx1,my1,mx2,my2)	((long)(mx1)<<15 ^ (long)(my1)<<10 ^ \
					(long)(mx2)<<5 ^ (long)(my2))


void
newlayer(void)				/* start a new layer */
{
	(void)memset((char *)hshtab, '\0', sizeof(hshtab));
	if (++layer >= 16) {
		mendpage();
		layer = 0;
	}
}


int
doline(		/* draw line conditionally */
	int	v1x,
	int v1y,
	int v2x,
	int v2y
)
{
	register int	h;

	if (v1x > v2x || (v1x == v2x && v1y > v2y)) {	/* sort endpoints */
		h=v1x; v1x=v2x; v2x=h;
		h=v1y; v1y=v2y; v2y=h;
	}
	h = hash(v1x, v1y, v2x, v2y) % HTBLSIZ;
	if (hshtab[h][0] == v1x && hshtab[h][1] == v1y &&
			hshtab[h][2] == v2x && hshtab[h][3] == v2y)
		return(0);
	hshtab[h][0] = v1x; hshtab[h][1] = v1y;
	hshtab[h][2] = v2x; hshtab[h][3] = v2y;
	if ((long)(v2x-v1x)*(v2x-v1x) + (long)(v2y-v1y)*(v2y-v1y)
			<= random() % rthresh)
		return(0);
	mline(v1x, v1y, layer/4, 0, layer%4);
	mdraw(v2x, v2y);
	return(1);
}
