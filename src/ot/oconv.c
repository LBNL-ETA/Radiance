#ifndef lint
static const char RCSid[] = "$Id: oconv.c,v 2.18 2003/10/22 02:06:34 greg Exp $";
#endif
/*
 *  oconv.c - main program for object to octree conversion.
 *
 *     7/29/85
 */

#include  "platform.h"
#include  "standard.h"
#include  "octree.h"
#include  "object.h"
#include  "otypes.h"
#include  "paths.h"

#define	 OMARGIN	(10*FTINY)	/* margin around global cube */

#define	 MAXOBJFIL	127		/* maximum number of scene files */

char  *progname;			/* argv[0] */

int  nowarn = 0;			/* supress warnings? */

int  objlim = 6;			/* # of objects before split */

int  resolu = 16384;			/* octree resolution limit */

CUBE  thescene = {EMPTY, {0.0, 0.0, 0.0}, 0.0};		/* our scene */

char  *ofname[MAXOBJFIL+1];		/* object file names */
int  nfiles = 0;			/* number of object files */

double	mincusize;			/* minimum cube size from resolu */

void  (*addobjnotify[])() = {NULL};	/* new object notifier functions */


main(argc, argv)		/* convert object files to an octree */
int  argc;
char  *argv[];
{
	FVECT  bbmin, bbmax;
	char  *infile = NULL;
	int  inpfrozen = 0;
	int  outflags = IO_ALL;
	OBJECT	startobj;
	int  i;

	progname = argv[0] = fixargv0(argv[0]);

	initotypes();

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case '\0':				/* scene from stdin */
			goto breakopt;
		case 'i':				/* input octree */
			infile = argv[++i];
			break;
		case 'b':				/* bounding cube */
			thescene.cuorg[0] = atof(argv[++i]) - OMARGIN;
			thescene.cuorg[1] = atof(argv[++i]) - OMARGIN;
			thescene.cuorg[2] = atof(argv[++i]) - OMARGIN;
			thescene.cusize = atof(argv[++i]) + 2*OMARGIN;
			break;
		case 'n':				/* set limit */
			objlim = atoi(argv[++i]);
			break;
		case 'r':				/* resolution limit */
			resolu = atoi(argv[++i]);
			break;
		case 'f':				/* freeze octree */
			outflags &= ~IO_FILES;
			break;
		case 'w':				/* supress warnings */
			nowarn = 1;
			break;
		default:
			sprintf(errmsg, "unknown option: '%s'", argv[i]);
			error(USER, errmsg);
			break;
		}
breakopt:
	SET_FILE_BINARY(stdout);
	if (infile != NULL) {		/* get old octree & objects */
		if (thescene.cusize > FTINY)
			error(USER, "only one of '-b' or '-i'");
		nfiles = readoct(infile, IO_ALL, &thescene, ofname);
		if (nfiles == 0)
			inpfrozen++;
	} else
		newheader("RADIANCE", stdout);	/* new binary file header */
	printargs(argc, argv, stdout);
	fputformat(OCTFMT, stdout);
	printf("\n");

	startobj = nobjects;		/* previous objects already converted */

	for ( ; i < argc; i++)		/* read new scene descriptions */
		if (!strcmp(argv[i], "-")) {	/* from stdin */
			readobj(NULL);
			outflags &= ~IO_FILES;
		} else {			/* from file */
			if (nfiles >= MAXOBJFIL)
				error(INTERNAL, "too many scene files");
			readobj(ofname[nfiles++] = argv[i]);
		}

	ofname[nfiles] = NULL;

	if (inpfrozen && outflags & IO_FILES) {
		error(WARNING, "frozen octree");
		outflags &= ~IO_FILES;
	}
						/* find bounding box */
	bbmin[0] = bbmin[1] = bbmin[2] = FHUGE;
	bbmax[0] = bbmax[1] = bbmax[2] = -FHUGE;
	for (i = startobj; i < nobjects; i++)
		add2bbox(objptr(i), bbmin, bbmax);
						/* set/check cube */
	if (thescene.cusize == 0.0) {
		if (bbmin[0] <= bbmax[0]) {
			for (i = 0; i < 3; i++) {
				bbmin[i] -= OMARGIN;
				bbmax[i] += OMARGIN;
			}
			for (i = 0; i < 3; i++)
				if (bbmax[i] - bbmin[i] > thescene.cusize)
					thescene.cusize = bbmax[i] - bbmin[i];
			for (i = 0; i < 3; i++)
				thescene.cuorg[i] =
					(bbmax[i]+bbmin[i]-thescene.cusize)*.5;
		}
	} else {
		for (i = 0; i < 3; i++)
			if (bbmin[i] < thescene.cuorg[i] ||
				bbmax[i] > thescene.cuorg[i] + thescene.cusize)
				error(USER, "boundary does not encompass scene");
	}

	mincusize = thescene.cusize / resolu - FTINY;

	for (i = startobj; i < nobjects; i++)		/* add new objects */
		addobject(&thescene, i);

	thescene.cutree = combine(thescene.cutree);	/* optimize */

	writeoct(outflags, &thescene, ofname);	/* write structures to stdout */

	quit(0);
}


void
quit(code)				/* exit program */
int  code;
{
	exit(code);
}


void
cputs()					/* interactive error */
{
	/* referenced, but not used */
}


void
wputs(s)				/* warning message */
char  *s;
{
	if (!nowarn)
		eputs(s);
}


void
eputs(s)				/* put string to stderr */
register char  *s;
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

				/* conflicting def's in param.h */
#undef  tstbit
#undef  setbit
#undef  clrbit
#undef  tglbit

#define	 bitop(f,i,op)		(f[((i)>>3)] op (1<<((i)&7)))
#define	 tstbit(f,i)		bitop(f,i,&)
#define	 setbit(f,i)		bitop(f,i,|=)
#define	 clrbit(f,i)		bitop(f,i,&=~)
#define	 tglbit(f,i)		bitop(f,i,^=)


addobject(cu, obj)			/* add an object to a cube */
register CUBE  *cu;
OBJECT	obj;
{
	int  inc;

	inc = (*ofun[objptr(obj)->otype].funp)(objptr(obj), cu);

	if (inc == O_MISS)
		return;				/* no intersection */

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
			addobject(&cukid, obj);
			octkid(cu->cutree, i) = cukid.cutree;
		}
		return;
	}
	if (isempty(cu->cutree)) {
		OBJECT  oset[2];		/* singular set */
		oset[0] = 1; oset[1] = obj;
		cu->cutree = fullnode(oset);
		return;
	}
					/* add to full node */
	add2full(cu, obj, inc);
}


add2full(cu, obj, inc)			/* add object to full node */
register CUBE  *cu;
OBJECT	obj;
int  inc;
{
	OCTREE	ot;
	OBJECT	oset[MAXSET+1];
	CUBE  cukid;
	unsigned char  inflg[(MAXSET+7)/8], volflg[(MAXSET+7)/8];
	register int  i, j;

	objset(oset, cu->cutree);
	cukid.cusize = cu->cusize * 0.5;

	if (inc==O_IN || oset[0] < objlim || cukid.cusize < mincusize) {
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
					/* mark volumes */
	j = (oset[0]+7)>>3;
	while (j--)
		volflg[j] = inflg[j] = 0;
	for (j = 1; j <= oset[0]; j++)
		if (isvolume(objptr(oset[j])->otype)) {
			setbit(volflg,j-1);
			if ((*ofun[objptr(oset[j])->otype].funp)
					(objptr(oset[j]), cu) == O_IN)
				setbit(inflg,j-1);
		}
					/* assign subcubes */
	for (i = 0; i < 8; i++) {
		cukid.cutree = EMPTY;
		for (j = 0; j < 3; j++) {
			cukid.cuorg[j] = cu->cuorg[j];
			if ((1<<j) & i)
				cukid.cuorg[j] += cukid.cusize;
		}
					/* surfaces first */
		for (j = 1; j <= oset[0]; j++)
			if (!tstbit(volflg,j-1))
				addobject(&cukid, oset[j]);
					/* then this object */
		addobject(&cukid, obj);
					/* then partial volumes */
		for (j = 1; j <= oset[0]; j++)
			if (tstbit(volflg,j-1) &&
					!tstbit(inflg,j-1))
				addobject(&cukid, oset[j]);
					/* full volumes last */
		for (j = 1; j <= oset[0]; j++)
			if (tstbit(inflg,j-1))
				addobject(&cukid, oset[j]);
					/* returned node */
		octkid(ot, i) = cukid.cutree;
	}
	cu->cutree = ot;
}
