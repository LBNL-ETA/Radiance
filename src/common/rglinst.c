#ifndef lint
static const char	RCSid[] = "$Id: rglinst.c,v 3.6 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 * Routines for reading instances and converting to OpenGL.
 */

#include "copyright.h"

#include "radogl.h"
#include "octree.h"

#define MAXLEVEL	16		/* maximum instance hierarchy level */

typedef struct {
	int	listid;				/* our list id */
	short	localmatl;			/* uses local material only */
	FVECT	cent;				/* center of octree cube */
	char	octfile[256];			/* octree file path */
} OCTINST;				/* octree to instantiate */

static double  ogetflt();
static long  ogetint();
static char  *ogetstr();
static int  loadobj();
static void  skiptree();
static void  octerror();
static OCTINST	*getoct();

static char  *infn;			/* input file name */
static FILE  *infp;			/* input file stream */
static int  objsize;			/* size of stored OBJECT's */
static short  otypmap[NUMOTYPE+8];	/* object type map */

static unsigned long	imhash(mod) char *mod; {return((unsigned long)mod);}
static LUTAB	imtab = {imhash,NULL,NULL,NULL,0,NULL,0};

static LUTAB	ottab = LU_SINIT(free,free);


int
o_instance(o)				/* convert instance to list call */
register OBJREC	*o;
{
	XF	xfs;
	register OCTINST	*ot;
					/* set up */
	if (o->oargs.nsargs < 1)
		objerror(o, USER, "missing octree");
	setmaterial(NULL, NULL, 0);
					/* put out transform (if any) */
	if (o->oargs.nsargs > 1) {
		if (xf(&xfs, o->oargs.nsargs-1, o->oargs.sarg+1) !=
				o->oargs.nsargs-1)
			objerror(o, USER, "bad transform");
		glPushAttrib(GL_TRANSFORM_BIT);
		if (xfs.sca < 1.-FTINY | xfs.sca > 1.+FTINY)
			glEnable(GL_NORMALIZE);
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
					/* matrix order works out to same */
#ifdef SMLFLT
		glMultMatrixf((GLfloat *)xfs.xfm);
#else
		glMultMatrixd((GLdouble *)xfs.xfm);
#endif
	}
	ot = getoct(o->oargs.sarg[0]);	/* get octree reference */
	if (ot->localmatl &= o->os != NULL)	/* set material */
		setmaterial((MATREC *)o->os, ot->cent, 0);
					/* call the assigned list */
	glCallList(ot->listid);

	if (o->oargs.nsargs > 1) {	/* end transform */
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();
		glPopAttrib();
	}
	rgl_checkerr("creating instance");
}


static int
buildoctlist(lp)			/* build octree list */
LUENT	*lp;
{
	int	old_dolights = dolights, old_domats = domats;
	register OCTINST	*op = (OCTINST *)lp->data;

	domats = !op->localmatl;	/* do materials only if needed */
	dolights = 0;			/* never do light sources */
	glNewList(op->listid, GL_COMPILE);
	loadoct(op->octfile);		/* load objects into display list */
	surfclean();			/* clean up */
	glEndList();
	dolights = old_dolights;	/* restore */
	domats = old_domats;
	return(1);			/* return success */
}


int
loadoctrees()				/* load octrees we've saved up */
{
	int	levelsleft = MAXLEVEL;
	int	nocts = 0;
	LUTAB	looptab;
				/* loop through new octree references */
	while (ottab.tsiz) {
		if (!levelsleft--)
			error(USER, "too many octree levels -- instance loop?");
		copystruct(&looptab, &ottab);
		ottab.tsiz = 0;
		nocts += lu_doall(&looptab, buildoctlist);
		lu_done(&looptab);
	}
	return(nocts);
}


static OCTINST *
getoct(name)				/* get/assign octree list id */
char	*name;
{
	extern char	*getpath(), *getlibpath();
	char	*path;
	register LUENT	*lp;
	register OCTINST	*op;

	if ((lp = lu_find(&ottab, name)) == NULL)
		goto memerr;
	if (lp->key == NULL) {
		lp->key = (char *)malloc(strlen(name)+1);
		if (lp->key == NULL)
			goto memerr;
		strcpy(lp->key, name);
	}
	if ((op = (OCTINST *)lp->data) == NULL) {
		path = getpath(name, getlibpath(), R_OK);
		if (path == NULL) {
			sprintf(errmsg, "cannot find octree \"%s\"", name);
			error(USER, errmsg);
		}
		op = (OCTINST *)(lp->data = (char *)malloc(sizeof(OCTINST)));
		strcpy(op->octfile, path);
		checkoct(op->octfile, op->cent);
		op->listid = newglist();
		op->localmatl = ~0;
	}
	return(op);
memerr:
	error(SYSTEM, "out of memory in getoct");
}


double
checkoct(fname, cent)			/* check octree file for validity */
char	*fname;
FVECT	cent;
{
	char  sbuf[64];
	FILE	*fp = infp;
	char	*fn = infn;
	double	siz = 0.;
	register int  i;
	
	if ((infp = fopen(infn=fname, "r")) == NULL) {
		sprintf(errmsg, "cannot open octree file \"%s\"", fname);
		error(SYSTEM, errmsg);
	}
#ifdef MSDOS
	setmode(fileno(infp), O_BINARY);
#endif
					/* get header */
	if (checkheader(infp, OCTFMT, NULL) < 0)
		octerror(USER, "not an octree");
					/* check format */
	if ((objsize = ogetint(2)-OCTMAGIC) <= 0 ||
			objsize > MAXOBJSIZ || objsize > sizeof(long))
		octerror("incompatible octree format");
	if (cent != NULL) {		/* get boundaries (compute center) */
		for (i = 0; i < 3; i++)
			cent[i] = atof(ogetstr(sbuf));
		siz = atof(ogetstr(sbuf))*.5;
		cent[0] += siz; cent[1] += siz; cent[2] += siz;
	} else {			/* get size (radius) only */
		for (i = 0; i < 3; i++)
			ogetstr(sbuf);
		siz = atof(ogetstr(sbuf))*.5;
	}
	fclose(infp);
	infp = fp;
	infn = fn;
	return(siz);
}


int
loadoct(fname)				/* read in objects from octree */
char  *fname;
{
	OBJECT  fnobjects;
	char  sbuf[256];
	int  nf;
	register int  i;
	long  m;
	
	infn = fname;
	infp = fopen(fname, "r");	/* assume already checked */
#ifdef MSDOS
	setmode(fileno(infp), O_BINARY);
#endif
					/* skip header */
	getheader(infp, NULL, NULL);
					/* get format */
	objsize = ogetint(2)-OCTMAGIC;
					/* skip boundaries */
	for (i = 0; i < 4; i++)
		ogetstr(sbuf);
	nf = 0;				/* load object files */
	while (*ogetstr(sbuf)) {
		rgl_load(sbuf);
		nf++;
	}
					/* get number of objects */
	fnobjects = m = ogetint(objsize);
	if (fnobjects != m)
		octerror(USER, "too many objects");

	if (nf == 0) {
		skiptree();
		for (i = 0; *ogetstr(sbuf); i++)
			if ((otypmap[i] = otype(sbuf)) < 0) {
				sprintf(errmsg, "unknown type \"%s\"", sbuf);
				octerror(WARNING, errmsg);
			}
		lu_init(&imtab, 1000); nobjects = 0;
		while (loadobj() != OVOID)
			;
		lu_done(&imtab);
		if (nobjects != fnobjects)
			octerror(USER, "inconsistent object count");
	}
	fclose(infp);
	return(nf);
}


static char *
ogetstr(s)			/* get null-terminated string */
char  *s;
{
	extern char  *getstr();

	if (getstr(s, infp) == NULL)
		octerror(USER, "truncated octree");
	return(s);
}


static long
ogetint(siz)			/* get a siz-byte integer */
int  siz;
{
	extern long  getint();
	register long  r;

	r = getint(siz, infp);
	if (feof(infp))
		octerror(USER, "truncated octree");
	return(r);
}


static double
ogetflt()			/* get a floating point number */
{
	extern double  getflt();
	double	r;

	r = getflt(infp);
	if (feof(infp))
		octerror(USER, "truncated octree");
	return(r);
}
	

static void
skiptree()				/* skip octree on input */
{
	register int  i;
	
	switch (getc(infp)) {
	case OT_EMPTY:
		return;
	case OT_FULL:
		for (i = ogetint(objsize)*objsize; i-- > 0; )
			if (getc(infp) == EOF)
				octerror(USER, "truncated octree");
		return;
	case OT_TREE:
		for (i = 0; i < 8; i++)
			skiptree();
		return;
	case EOF:
		octerror(USER, "truncated octree");
	default:
		octerror(USER, "damaged octree");
	}
}


static int
loadobj()				/* get next object */
{
	static OBJREC	 ob;
	char  idbuf[MAXSTR], sbuf[MAXSTR];
	register LUENT	*lep;
	register int  i;
	register long  m;
					/* get type */
	i = ogetint(1);
	if (i == -1)
		return(OVOID);		/* terminator */
	if ((ob.otype = otypmap[i]) < 0)
		octerror(USER, "reference to unknown type");
					/* get modifier */
	if ((m = ogetint(objsize)) != OVOID && (OBJECT)m != m)
		octerror(USER, "too many objects");
	if ((ob.omod = m) != OVOID && domats) {
		if ((lep = lu_find(&imtab, (char *)m)) == NULL)
			goto memerr;
		ob.os = lep->data;
	} else
		ob.os = NULL;
					/* get name id */
	ob.oname = ogetstr(idbuf);
					/* get string arguments */
	if (ob.oargs.nsargs = ogetint(2)) {
		ob.oargs.sarg = (char **)malloc
				(ob.oargs.nsargs*sizeof(char *));
		if (ob.oargs.sarg == NULL)
			goto memerr;
		for (i = 0; i < ob.oargs.nsargs; i++)
			ob.oargs.sarg[i] = savestr(ogetstr(sbuf));
	} else
		ob.oargs.sarg = NULL;
					/* get integer arguments */
#ifdef	IARGS
	if (ob.oargs.niargs = ogetint(2)) {
		ob.oargs.iarg = (long *)malloc
				(ob.oargs.niargs*sizeof(long));
		if (ob.oargs.iarg == NULL)
			goto memerr;
		for (i = 0; i < ob.oargs.niargs; i++)
			ob.oargs.iarg[i] = ogetint(4);
	} else
		ob.oargs.iarg = NULL;
#endif
					/* get real arguments */
	if (ob.oargs.nfargs = ogetint(2)) {
		ob.oargs.farg = (FLOAT *)malloc
				(ob.oargs.nfargs*sizeof(FLOAT));
		if (ob.oargs.farg == NULL)
			goto memerr;
		for (i = 0; i < ob.oargs.nfargs; i++)
			ob.oargs.farg[i] = ogetflt();
	} else
		ob.oargs.farg = NULL;
					/* process object */
	(*ofun[ob.otype].funp)(&ob);
					/* record material if modifier */
	if (ismodifier(ob.otype)) {
		if ((lep = lu_find(&imtab, (char *)nobjects)) == NULL)
			goto memerr;
		lep->key = (char *)nobjects;
		lep->data = (char *)getmatp(ob.oname);
	}
	freefargs(&ob.oargs);		/* free arguments */
	return(nobjects++);		/* return object id */
memerr:
	error(SYSTEM, "out of memory in loadobj");
}


static void
octerror(etyp, msg)			/* octree error */
int  etyp;
char  *msg;
{
	char  msgbuf[128];

	sprintf(msgbuf, "(%s): %s", infn, msg);
	error(etyp, msgbuf);
}
