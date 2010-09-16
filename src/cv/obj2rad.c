#ifndef lint
static const char	RCSid[] = "$Id: obj2rad.c,v 2.26 2010/09/16 03:28:21 greg Exp $";
#endif
/*
 * Convert a Wavefront .obj file to Radiance format.
 *
 * Currently, we support only polygonal geometry.  Non-planar
 * faces are broken rather haphazardly into triangles.
 * Also, texture map indices only work for triangles, though
 * I'm not sure they work correctly.  (Taken out -- see TEXMAPS defines.)
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "rtmath.h"
#include "rtio.h"
#include "resolu.h"
#include "trans.h"
#include "tmesh.h"


#define PATNAME		"M-pat"		/* mesh pattern name (reused) */
#define TEXNAME		"M-nor"		/* mesh texture name (reused) */
#define DEFOBJ		"unnamed"	/* default object name */
#define DEFMAT		"white"		/* default material name */

#define pvect(v)	printf("%18.12g %18.12g %18.12g\n",(v)[0],(v)[1],(v)[2])

FVECT	*vlist;			/* our vertex list */
int	nvs;			/* number of vertices in our list */
FVECT	*vnlist;		/* vertex normal list */
int	nvns;
RREAL	(*vtlist)[2];		/* map vertex list */
int	nvts;

int	ndegen = 0;		/* count of degenerate faces */
int	n0norm = 0;		/* count of zero normals */

typedef int	VNDX[3];	/* vertex index (point,map,normal) */

#define CHUNKSIZ	1024	/* vertex allocation chunk size */

#define MAXARG		512	/* maximum # arguments in a statement */

				/* qualifiers */
#define Q_MTL		0
#define Q_MAP		1
#define Q_GRP		2
#define Q_OBJ		3
#define Q_FAC		4
#define NQUALS		5

char	*qname[NQUALS] = {
	"Material",
	"Map",
	"Group",
	"Object",
	"Face",
};

QLIST	qlist = {NQUALS, qname};
				/* valid qualifier ids */
IDLIST	qual[NQUALS];
				/* mapping rules */
RULEHD	*ourmapping = NULL;

char	*defmat = DEFMAT;	/* default (starting) material name */
char	*defobj = DEFOBJ;	/* default (starting) object name */

int	flatten = 0;		/* discard surface normal information */

char	mapname[128];		/* current picture file */
char	matname[64];		/* current material name */
char	group[16][32];		/* current group names */
char	objname[128];		/* current object name */
char	*inpfile;		/* input file name */
int	lineno;			/* current line number */
int	faceno;			/* current face number */

static void getnames(FILE *fp);
static void convert(FILE *fp);
static int getstmt(char *av[MAXARG], FILE *fp);
static char * getmtl(void);
static char * getonm(void);
static int matchrule(RULEHD *rp);
static int cvtndx(VNDX vi, char *vs);
static int nonplanar(int ac, char **av);
static int putface(int ac, char **av);
static int puttri(char *v1, char *v2, char *v3);
static void freeverts(void);
static int newv(double x, double y, double z);
static int newvn(double x, double y, double z);
static int newvt(double x, double y);
static void syntax(char *er);


int
main(		/* read in .obj file and convert */
	int	argc,
	char	*argv[]
)
{
	int	donames = 0;
	int	i;

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'o':		/* object name */
			defobj = argv[++i];
			break;
		case 'n':		/* just produce name list */
			donames++;
			break;
		case 'm':		/* use custom mapfile */
			ourmapping = getmapping(argv[++i], &qlist);
			break;
		case 'f':		/* flatten surfaces */
			flatten++;
			break;
		default:
			goto userr;
		}
	if ((i > argc) | (i < argc-1))
		goto userr;
	if (i == argc)
		inpfile = "<stdin>";
	else if (freopen(inpfile=argv[i], "r", stdin) == NULL) {
		fprintf(stderr, "%s: cannot open\n", inpfile);
		exit(1);
	}
	if (donames) {				/* scan for ids */
		getnames(stdin);
		printf("filename \"%s\"\n", inpfile);
		printf("filetype \"Wavefront\"\n");
		write_quals(&qlist, qual, stdout);
		printf("qualifier %s begin\n", qlist.qual[Q_FAC]);
		printf("[%d:%d]\n", 1, faceno);
		printf("end\n");
	} else {				/* translate file */
		printf("# ");
		printargs(argc, argv, stdout);
		convert(stdin);
	}
	if (ndegen)
		printf("# %d degenerate faces\n", ndegen);
	if (n0norm)
		printf("# %d invalid (zero) normals\n", n0norm);
	exit(0);
userr:
	fprintf(stderr, "Usage: %s [-o obj][-m mapping][-n][-f] [file.obj]\n",
			argv[0]);
	exit(1);
}


void
getnames(			/* get valid qualifier names */
	FILE	*fp
)
{
	char	*argv[MAXARG];
	int	argc;
	ID	tmpid;
	register int	i;

	while ( (argc = getstmt(argv, fp)) )
		switch (argv[0][0]) {
		case 'f':				/* face */
			if (!argv[0][1])
				faceno++;
			break;
		case 'u':
			if (!strcmp(argv[0], "usemtl")) {	/* material */
				if (argc < 2)
					break;		/* not fatal */
				tmpid.number = 0;
				tmpid.name = argv[1];
				findid(&qual[Q_MTL], &tmpid, 1);
			} else if (!strcmp(argv[0], "usemap")) {/* map */
				if (argc < 2 || !strcmp(argv[1], "off"))
					break;		/* not fatal */
				tmpid.number = 0;
				tmpid.name = argv[1];
				findid(&qual[Q_MAP], &tmpid, 1);
			}
			break;
		case 'o':		/* object name */
			if (argv[0][1] || argc < 2)
				break;
			tmpid.number = 0;
			tmpid.name = argv[1];
			findid(&qual[Q_OBJ], &tmpid, 1);
			break;
		case 'g':		/* group name(s) */
			if (argv[0][1])
				break;
			tmpid.number = 0;
			for (i = 1; i < argc; i++) {
				tmpid.name = argv[i];
				findid(&qual[Q_GRP], &tmpid, 1);
			}
			break;
		}
}


void
convert(			/* convert a T-mesh */
	FILE	*fp
)
{
	char	*argv[MAXARG];
	int	argc;
	int	nstats, nunknown;
	register int	i;

	nstats = nunknown = 0;
					/* scan until EOF */
	while ( (argc = getstmt(argv, fp)) ) {
		switch (argv[0][0]) {
		case 'v':		/* vertex */
			switch (argv[0][1]) {
			case '\0':			/* point */
				if (badarg(argc-1,argv+1,"fff"))
					syntax("Bad vertex");
				newv(atof(argv[1]), atof(argv[2]),
						atof(argv[3]));
				break;
			case 'n':			/* normal */
				if (argv[0][2])
					goto unknown;
				if (badarg(argc-1,argv+1,"fff"))
					syntax("Bad normal");
				if (!newvn(atof(argv[1]), atof(argv[2]),
						atof(argv[3])))
					syntax("Zero normal");
				break;
			case 't':			/* texture map */
				if (argv[0][2])
					goto unknown;
				if (badarg(argc-1,argv+1,"ff"))
					goto unknown;
				newvt(atof(argv[1]), atof(argv[2]));
				break;
			default:
				goto unknown;
			}
			break;
		case 'f':				/* face */
			if (argv[0][1])
				goto unknown;
			faceno++;
			switch (argc-1) {
			case 0: case 1: case 2:
				syntax("Too few vertices");
				break;
			case 3:
				if (!puttri(argv[1], argv[2], argv[3]))
					syntax("Bad triangle");
				break;
			default:
				if (!putface(argc-1, argv+1))
					syntax("Bad face");
				break;
			}
			break;
		case 'u':
			if (!strcmp(argv[0], "usemtl")) {	/* material */
				if (argc < 2)
					break;		/* not fatal */
				strcpy(matname, argv[1]);
			} else if (!strcmp(argv[0], "usemap")) {/* map */
				if (argc < 2)
					break;		/* not fatal */
				if (!strcmp(argv[1], "off"))
					mapname[0] = '\0';
				else
					sprintf(mapname, "%s.hdr", argv[1]);
			} else
				goto unknown;
			break;
		case 'o':		/* object name */
			if (argv[0][1])
				goto unknown;
			if (argc < 2)
				break;		/* not fatal */
			strcpy(objname, argv[1]);
			break;
		case 'g':		/* group name(s) */
			if (argv[0][1])
				goto unknown;
			for (i = 1; i < argc; i++)
				strcpy(group[i-1], argv[i]);
			group[i-1][0] = '\0';
			break;
		case '#':		/* comment */
			printargs(argc, argv, stdout);
			break;
		default:;		/* something we don't deal with */
		unknown:
			nunknown++;
			break;
		}
		nstats++;
	}
	printf("\n# Done processing file: %s\n", inpfile);
	printf("# %d lines, %d statements, %d unrecognized\n",
			lineno, nstats, nunknown);
}


int
getstmt(				/* read the next statement from fp */
	register char	*av[MAXARG],
	FILE	*fp
)
{
	static char	sbuf[MAXARG*16];
	register char	*cp;
	register int	i;

	do {
		if (fgetline(cp=sbuf, sizeof(sbuf), fp) == NULL)
			return(0);
		i = 0;
		for ( ; ; ) {
			while (isspace(*cp) || *cp == '\\') {
				if (*cp == '\n')
					lineno++;
				*cp++ = '\0';
			}
			if (!*cp)
				break;
			if (i >= MAXARG-1) {
				fprintf(stderr,
			"warning: line %d: too many arguments (limit %d)\n",
					lineno+1, MAXARG-1);
				break;
			}
			av[i++] = cp;
			while (*++cp && !isspace(*cp))
				;
		}
		av[i] = NULL;
		lineno++;
	} while (!i);

	return(i);
}


char *
getmtl(void)				/* figure material for this face */
{
	register RULEHD	*rp = ourmapping;

	if (rp == NULL) {		/* no rule set */
		if (matname[0])
			return(matname);
		if (group[0][0])
			return(group[0]);
		return(defmat);
	}
					/* check for match */
	do {
		if (matchrule(rp)) {
			if (!strcmp(rp->mnam, VOIDID))
				return(NULL);	/* match is null */
			return(rp->mnam);
		}
		rp = rp->next;
	} while (rp != NULL);
					/* no match found */
	return(NULL);
}


char *
getonm(void)				/* invent a good name for object */
{
	static char	name[64];
	register char	*cp1, *cp2;
	register int	i;
					/* check for preset */
	if (objname[0])
		return(objname);
	if (!group[0][0])
		return(defobj);
	cp1 = name;			/* else make name out of groups */
	for (i = 0; group[i][0]; i++) {
		cp2 = group[i];
		if (cp1 > name)
			*cp1++ = '.';
		while ( (*cp1 = *cp2++) )
			if (++cp1 >= name+sizeof(name)-2) {
				*cp1 = '\0';
				return(name);
			}
	}
	return(name);
}


int
matchrule(				/* check for a match on this rule */
	register RULEHD	*rp
)
{
	ID	tmpid;
	int	gotmatch;
	register int	i;

	if (rp->qflg & FL(Q_MTL)) {
		if (!matname[0])
			return(0);
		tmpid.number = 0;
		tmpid.name = matname;
		if (!matchid(&tmpid, &idm(rp)[Q_MTL]))
			return(0);
	}
	if (rp->qflg & FL(Q_MAP)) {
		if (!mapname[0])
			return(0);
		tmpid.number = 0;
		tmpid.name = mapname;
		if (!matchid(&tmpid, &idm(rp)[Q_MAP]))
			return(0);
	}
	if (rp->qflg & FL(Q_GRP)) {
		tmpid.number = 0;
		gotmatch = 0;
		for (i = 0; group[i][0]; i++) {
			tmpid.name = group[i];
			gotmatch |= matchid(&tmpid, &idm(rp)[Q_GRP]);
		}
		if (!gotmatch)
			return(0);
	}
	if (rp->qflg & FL(Q_OBJ)) {
		if (!objname[0])
			return(0);
		tmpid.number = 0;
		tmpid.name = objname;
		if (!matchid(&tmpid, &idm(rp)[Q_OBJ]))
			return(0);
	}
	if (rp->qflg & FL(Q_FAC)) {
		tmpid.name = NULL;
		tmpid.number = faceno;
		if (!matchid(&tmpid, &idm(rp)[Q_FAC]))
			return(0);
	}
	return(1);
}


int
cvtndx(				/* convert vertex string to index */
	register VNDX	vi,
	register char	*vs
)
{
					/* get point */
	vi[0] = atoi(vs);
	if (vi[0] > 0) {
		if (vi[0]-- > nvs)
			return(0);
	} else if (vi[0] < 0) {
		vi[0] += nvs;
		if (vi[0] < 0)
			return(0);
	} else
		return(0);
					/* get map */
	while (*vs)
		if (*vs++ == '/')
			break;
	vi[1] = atoi(vs);
	if (vi[1] > 0) {
		if (vi[1]-- > nvts)
			return(0);
	} else if (vi[1] < 0) {
		vi[1] += nvts;
		if (vi[1] < 0)
			return(0);
	} else
		vi[1] = -1;
					/* get normal */
	while (*vs)
		if (*vs++ == '/')
			break;
	vi[2] = atoi(vs);
	if (vi[2] > 0) {
		if (vi[2]-- > nvns)
			return(0);
	} else if (vi[2] < 0) {
		vi[2] += nvns;
		if (vi[2] < 0)
			return(0);
	} else
		vi[2] = -1;
					/* zero normal is not normal */
	if (vi[2] >= 0 && DOT(vnlist[vi[2]],vnlist[vi[2]]) <= FTINY)
		vi[2] = -1;
	return(1);
}


int
nonplanar(			/* are vertices non-planar? */
	register int	ac,
	register char	**av
)
{
	VNDX	vi;
	RREAL	*p0, *p1;
	FVECT	v1, v2, nsum, newn;
	double	d;
	register int	i;

	if (!cvtndx(vi, av[0]))
		return(0);
	if (!flatten && vi[2] >= 0)
		return(1);		/* has interpolated normals */
	if (ac < 4)
		return(0);		/* it's a triangle! */
					/* set up */
	p0 = vlist[vi[0]];
	if (!cvtndx(vi, av[1]))
		return(0);		/* error gets caught later */
	nsum[0] = nsum[1] = nsum[2] = 0.;
	p1 = vlist[vi[0]];
	fvsum(v2, p1, p0, -1.0);
	for (i = 2; i < ac; i++) {
		VCOPY(v1, v2);
		if (!cvtndx(vi, av[i]))
			return(0);
		p1 = vlist[vi[0]];
		fvsum(v2, p1, p0, -1.0);
		fcross(newn, v1, v2);
		if (normalize(newn) == 0.0) {
			if (i < 3)
				return(1);	/* can't deal with this */
			fvsum(nsum, nsum, nsum, 1./(i-2));
			continue;
		}
		d = fdot(newn,nsum);
		if (d >= 0) {
			if (d < (1.0-FTINY)*(i-2))
				return(1);
			fvsum(nsum, nsum, newn, 1.0);
		} else {
			if (d > -(1.0-FTINY)*(i-2))
				return(1);
			fvsum(nsum, nsum, newn, -1.0);
		}
	}
	return(0);
}


int
putface(				/* put out an N-sided polygon */
	int	ac,
	register char	**av
)
{
	VNDX	vi;
	char	*cp;
	register int	i;

	if (nonplanar(ac, av)) {	/* break into triangles */
		while (ac > 2) {
			if (!puttri(av[0], av[1], av[2]))
				return(0);
			ac--;		/* remove vertex & rotate */
			cp = av[0];
			for (i = 0; i < ac-1; i++)
				av[i] = av[i+2];
			av[i] = cp;
		}
		return(1);
	}
	if ((cp = getmtl()) == NULL)
		return(0);
	printf("\n%s polygon %s.%d\n", cp, getonm(), faceno);
	printf("0\n0\n%d\n", 3*ac);
	for (i = 0; i < ac; i++) {
		if (!cvtndx(vi, av[i]))
			return(0);
		pvect(vlist[vi[0]]);
	}
	return(1);
}


int
puttri(			/* put out a triangle */
	char *v1,
	char *v2,
	char *v3
)
{
	char	*mod;
	VNDX	v1i, v2i, v3i;
	BARYCCM	bvecs;
	RREAL	bcoor[3][3];
	int	texOK = 0, patOK;
	int	flatness;
	register int	i;

	if ((mod = getmtl()) == NULL)
		return(0);

	if (!cvtndx(v1i, v1) || !cvtndx(v2i, v2) || !cvtndx(v3i, v3))
		return(0);
					/* compute barycentric coordinates */
	if (v1i[2]>=0 && v2i[2]>=0 && v3i[2]>=0)
		flatness = flat_tri(vlist[v1i[0]], vlist[v2i[0]], vlist[v3i[0]],
				vnlist[v1i[2]], vnlist[v2i[2]], vnlist[v3i[2]]);
	else
		flatness = ISFLAT;

	switch (flatness) {
	case DEGEN:			/* zero area */
		ndegen++;
		return(-1);
	case RVFLAT:			/* reversed normals, but flat */
	case ISFLAT:			/* smoothing unnecessary */
		texOK = 0;
		break;
	case RVBENT:			/* reversed normals with smoothing */
	case ISBENT:			/* proper smoothing */
		texOK = 1;
		break;
	}
	if (flatten)
		texOK = 0;
#ifdef TEXMAPS
	patOK = mapname[0] && (v1i[1]>=0 && v2i[1]>=0 && v3i[1]>=0);
#else
	patOK = 0;
#endif
	if (texOK | patOK)
		if (comp_baryc(&bvecs, vlist[v1i[0]], vlist[v2i[0]],
				vlist[v3i[0]]) < 0)
			texOK = patOK = 0;
					/* put out texture (if any) */
	if (texOK) {
		printf("\n%s texfunc %s\n", mod, TEXNAME);
		mod = TEXNAME;
		printf("4 dx dy dz %s\n", TCALNAME);
		printf("0\n");
		for (i = 0; i < 3; i++) {
			bcoor[i][0] = vnlist[v1i[2]][i];
			bcoor[i][1] = vnlist[v2i[2]][i];
			bcoor[i][2] = vnlist[v3i[2]][i];
		}
		put_baryc(&bvecs, bcoor, 3);
	}
#ifdef TEXMAPS
					/* put out pattern (if any) */
	if (patOK) {
		printf("\n%s colorpict %s\n", mod, PATNAME);
		mod = PATNAME;
		printf("7 noneg noneg noneg %s %s u v\n", mapname, TCALNAME);
		printf("0\n");
		for (i = 0; i < 2; i++) {
			bcoor[i][0] = vtlist[v1i[1]][i];
			bcoor[i][1] = vtlist[v2i[1]][i];
			bcoor[i][2] = vtlist[v3i[1]][i];
		}
		put_baryc(&bvecs, bcoor, 2);
	}
#endif
					/* put out (reversed) triangle */
	printf("\n%s polygon %s.%d\n", mod, getonm(), faceno);
	printf("0\n0\n9\n");
	if (flatness == RVFLAT || flatness == RVBENT) {
		pvect(vlist[v3i[0]]);
		pvect(vlist[v2i[0]]);
		pvect(vlist[v1i[0]]);
	} else {
		pvect(vlist[v1i[0]]);
		pvect(vlist[v2i[0]]);
		pvect(vlist[v3i[0]]);
	}
	return(1);
}


void
freeverts(void)			/* free all vertices */
{
	if (nvs) {
		free((void *)vlist);
		nvs = 0;
	}
	if (nvts) {
		free((void *)vtlist);
		nvts = 0;
	}
	if (nvns) {
		free((void *)vnlist);
		nvns = 0;
	}
}


int
newv(			/* create a new vertex */
	double	x,
	double	y,
	double	z
)
{
	if (!(nvs%CHUNKSIZ)) {		/* allocate next block */
		if (nvs == 0)
			vlist = (FVECT *)malloc(CHUNKSIZ*sizeof(FVECT));
		else
			vlist = (FVECT *)realloc((void *)vlist,
					(nvs+CHUNKSIZ)*sizeof(FVECT));
		if (vlist == NULL) {
			fprintf(stderr,
			"Out of memory while allocating vertex %d\n", nvs);
			exit(1);
		}
	}
					/* assign new vertex */
	vlist[nvs][0] = x;
	vlist[nvs][1] = y;
	vlist[nvs][2] = z;
	return(++nvs);
}


int
newvn(			/* create a new vertex normal */
	double	x,
	double	y,
	double	z
)
{
	if (!(nvns%CHUNKSIZ)) {		/* allocate next block */
		if (nvns == 0)
			vnlist = (FVECT *)malloc(CHUNKSIZ*sizeof(FVECT));
		else
			vnlist = (FVECT *)realloc((void *)vnlist,
					(nvns+CHUNKSIZ)*sizeof(FVECT));
		if (vnlist == NULL) {
			fprintf(stderr,
			"Out of memory while allocating normal %d\n", nvns);
			exit(1);
		}
	}
					/* assign new normal */
	vnlist[nvns][0] = x;
	vnlist[nvns][1] = y;
	vnlist[nvns][2] = z;
	n0norm += (normalize(vnlist[nvns]) == 0.0);
	return(++nvns);
}


int
newvt(			/* create a new texture map vertex */
	double	x,
	double	y
)
{
	if (!(nvts%CHUNKSIZ)) {		/* allocate next block */
		if (nvts == 0)
			vtlist = (RREAL (*)[2])malloc(CHUNKSIZ*2*sizeof(RREAL));
		else
			vtlist = (RREAL (*)[2])realloc((void *)vtlist,
					(nvts+CHUNKSIZ)*2*sizeof(RREAL));
		if (vtlist == NULL) {
			fprintf(stderr,
			"Out of memory while allocating texture vertex %d\n",
					nvts);
			exit(1);
		}
	}
					/* assign new vertex */
	vtlist[nvts][0] = x;
	vtlist[nvts][1] = y;
	return(++nvts);
}


void
syntax(			/* report syntax error and exit */
	char	*er
)
{
	fprintf(stderr, "%s: Wavefront syntax error near line %d: %s\n",
			inpfile, lineno, er);
	exit(1);
}
