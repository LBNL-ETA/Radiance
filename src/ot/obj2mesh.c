#ifndef lint
static const char RCSid[] = "$Id: obj2mesh.c,v 2.9 2004/04/22 17:33:48 greg Exp $";
#endif
/*
 *  Main program to compile a Wavefront .OBJ file into a Radiance mesh
 */

#include "copyright.h"
#include "platform.h"
#include "standard.h"
#include "resolu.h"
#include "cvmesh.h"
#include "otypes.h"

extern int	o_face(); /* XXX should go to a header file */

int	o_default() { return(O_MISS); }

FUN  ofun[NUMOTYPE] = INIT_OTYPE;	/* needed for link resolution */

char  *progname;			/* argv[0] */

int  nowarn = 0;			/* supress warnings? */

int  objlim = 15;			/* # of objects before split */

int  resolu = 16384;			/* octree resolution limit */

double	mincusize;			/* minimum cube size from resolu */

static void addface(CUBE  *cu, OBJECT	obj);
static void add2full(CUBE  *cu, OBJECT	obj);


int
main(		/* compile a .OBJ file into a mesh */
	int  argc,
	char  *argv[]
)
{
	int  nmatf = 0;
	char  *matinp[32];
	int  i, j;

	progname = argv[0];
	ofun[OBJ_FACE].funp = o_face;

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'n':				/* set limit */
			objlim = atoi(argv[++i]);
			break;
		case 'r':				/* resolution limit */
			resolu = atoi(argv[++i]);
			break;
		case 'a':				/* material file */
			matinp[nmatf++] = argv[++i];
			break;
		case 'w':				/* supress warnings */
			nowarn = 1;
			break;
		default:
			sprintf(errmsg, "unknown option: '%s'", argv[i]);
			error(USER, errmsg);
			break;
		}

	if (i < argc-2)
		error(USER, "too many file arguments");
					/* initialize mesh */
	cvinit(i==argc-2 ? argv[i+1] : "<stdout>");
					/* load material input */
	for (j = 0; j < nmatf; j++)
		readobj(matinp[j]);
					/* read .OBJ file into triangles */
	if (i == argc)
		wfreadobj(NULL);
	else
		wfreadobj(argv[i]);
	
	cvmeshbounds();			/* set octree boundaries */

	if (i == argc-2)		/* open output file */
		if (freopen(argv[i+1], "w", stdout) == NULL)
			error(SYSTEM, "cannot open output file");
	SET_FILE_BINARY(stdout);
	newheader("RADIANCE", stdout);	/* new binary file header */
	printargs(i<argc ? i+1 : argc, argv, stdout);
	fputformat(MESHFMT, stdout);
	fputc('\n', stdout);

	mincusize = ourmesh->mcube.cusize / resolu - FTINY;

	for (i = 0; i < nobjects; i++)	/* add triangles to octree */
		if (objptr(i)->otype == OBJ_FACE)
			addface(&ourmesh->mcube, i);

					/* optimize octree */
	ourmesh->mcube.cutree = combine(ourmesh->mcube.cutree);

	if (ourmesh->mcube.cutree == EMPTY)
		error(WARNING, "mesh is empty");
	
	cvmesh();			/* convert mesh and leaf nodes */

	writemesh(ourmesh, stdout);	/* write mesh to output */
	
	/* printmeshstats(ourmesh, stderr); */

	quit(0);
	return 0; /* pro forma return */
}


void
quit(				/* exit program */
	int  code
)
{
	exit(code);
}


void
cputs(void)					/* interactive error */
{
	/* referenced, but not used */
}


void
wputs(				/* warning message */
	char  *s
)
{
	if (!nowarn)
		eputs(s);
}


void
eputs(				/* put string to stderr */
	register char  *s
)
{
	static int  inln = 0;

	if (!inln++) {
		fputs(progname, stderr);
		fputs(": ", stderr);
	}
	fputs(s, stderr);
	if (*s && s[strlen(s)-1] == '\n')
		inln = 0;
}


static void
addface(			/* add a face to a cube */
	register CUBE  *cu,
	OBJECT	obj
)
{

	if (o_face(objptr(obj), cu) == O_MISS)
		return;

	if (istree(cu->cutree)) {
		CUBE  cukid;			/* do children */
		int  i, j;
		cukid.cusize = cu->cusize * 0.5;
		for (i = 0; i < 8; i++) {
			cukid.cutree = octkid(cu->cutree, i);
			for (j = 0; j < 3; j++) {
				cukid.cuorg[j] = cu->cuorg[j];
				if ((1<<j) & i)
					cukid.cuorg[j] += cukid.cusize;
			}
			addface(&cukid, obj);
			octkid(cu->cutree, i) = cukid.cutree;
		}
		return;
	}
	if (isempty(cu->cutree)) {
		OBJECT	oset[2];		/* singular set */
		oset[0] = 1; oset[1] = obj;
		cu->cutree = fullnode(oset);
		return;
	}
					/* add to full node */
	add2full(cu, obj);
}


static void
add2full(			/* add object to full node */
	register CUBE  *cu,
	OBJECT	obj
)
{
	OCTREE	ot;
	OBJECT	oset[MAXSET+1];
	CUBE  cukid;
	register int  i, j;

	objset(oset, cu->cutree);
	cukid.cusize = cu->cusize * 0.5;

	if (inc==O_IN || oset[0] < objlim || cukid.cusize <
			(oset[0] < MAXSET ? mincusize : mincusize/256.0)) {
						/* add to set */
		if (oset[0] >= MAXSET) {
			sprintf(errmsg, "set overflow in addobject (%s)",
					objptr(obj)->oname);
			error(INTERNAL, errmsg);
		}
		insertelem(oset, obj);
		cu->cutree = fullnode(oset);
		return;
	}
					/* subdivide cube */
	if ((ot = octalloc()) == EMPTY)
		error(SYSTEM, "out of octree space");
					/* assign subcubes */
	for (i = 0; i < 8; i++) {
		cukid.cutree = EMPTY;
		for (j = 0; j < 3; j++) {
			cukid.cuorg[j] = cu->cuorg[j];
			if ((1<<j) & i)
				cukid.cuorg[j] += cukid.cusize;
		}
		for (j = 1; j <= oset[0]; j++)
			addface(&cukid, oset[j]);
		addface(&cukid, obj);
					/* returned node */
		octkid(ot, i) = cukid.cutree;
	}
	cu->cutree = ot;
}
