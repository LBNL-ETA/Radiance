#ifndef lint
static const char	RCSid[] = "$Id: readoct.c,v 2.21 2003/07/21 22:30:17 schorsch Exp $";
#endif
/*
 *  readoct.c - routines to read octree information.
 */

#include "copyright.h"

#include  <stdio.h>
#include  <time.h>

#include  "standard.h"
#include  "platform.h"
#include  "octree.h"
#include  "object.h"
#include  "otypes.h"
#include  "resolu.h"

static double  ogetflt(void);
static long  ogetint(int);
static char  *ogetstr(char *);
static int  nonsurfinset(OBJECT *);
static void  octerror(int  etyp, char  *msg);
static void  skiptree(void);
static OCTREE  getfullnode(void), gettree(void);

static char  *infn;			/* input file specification */
static FILE  *infp;			/* input file stream */
static int  objsize;			/* size of stored OBJECT's */
static OBJECT  objorig;			/* zeroeth object */
static OBJECT  fnobjects;		/* number of objects in this file */


int
readoct(inpspec, load, scene, ofn)	/* read in octree file or stream */
char  *inpspec;
int  load;
CUBE  *scene;
char  *ofn[];
{
	char  sbuf[512];
	int  nf;
	int  i;
	long  m;
	
	if (inpspec == NULL) {
		infn = "standard input";
		infp = stdin;
	} else if (inpspec[0] == '!') {
		infn = inpspec;
		if ((infp = popen(inpspec+1, "r")) == NULL) {
			sprintf(errmsg, "cannot execute \"%s\"", inpspec);
			error(SYSTEM, errmsg);
		}
	} else {
		infn = inpspec;
		if ((infp = fopen(inpspec, "r")) == NULL) {
			sprintf(errmsg, "cannot open octree file \"%s\"",
					inpspec);
			error(SYSTEM, errmsg);
		}
	}
	SET_FILE_BINARY(infp);
					/* get header */
	if (checkheader(infp, OCTFMT, load&IO_INFO ? stdout : (FILE *)NULL) < 0)
		octerror(USER, "not an octree");
					/* check format */
	if ((objsize = ogetint(2)-OCTMAGIC) <= 0 ||
			objsize > MAXOBJSIZ || objsize > sizeof(long))
		octerror(USER, "incompatible octree format");
					/* get boundaries */
	if (load & IO_BOUNDS) {
		for (i = 0; i < 3; i++)
			scene->cuorg[i] = atof(ogetstr(sbuf));
		scene->cusize = atof(ogetstr(sbuf));
	} else {
		for (i = 0; i < 4; i++)
			ogetstr(sbuf);
	}
	objorig = nobjects;		/* set object offset */
	nf = 0;				/* get object files */
	while (*ogetstr(sbuf)) {
		if (load & IO_SCENE)
			readobj(sbuf);
		if (load & IO_FILES)
			ofn[nf] = savqstr(sbuf);
		nf++;
	}
	if (load & IO_FILES)
		ofn[nf] = NULL;
					/* get number of objects */
	fnobjects = m = ogetint(objsize);
	if (fnobjects != m)
		octerror(USER, "too many objects");

	if (load & IO_TREE)		/* get the octree */
		scene->cutree = gettree();
	else if (load & IO_SCENE && nf == 0)
		skiptree();
		
	if (load & IO_SCENE) {		/* get the scene */
	    if (nf == 0) {
					/* load binary scene data */
		readscene(infp, objsize);

	    } else {			/* consistency checks */
				/* check object count */
		if (nobjects != objorig+fnobjects)
			octerror(USER, "bad object count; octree stale?");
				/* check for non-surfaces */
		if (dosets(nonsurfinset))
			octerror(USER, "modifier in tree; octree stale?");
	    }
	}
				/* close the input */
	if (infn[0] == '!')
		pclose(infp);
	else
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


static OCTREE
getfullnode()			/* get a set, return fullnode */
{
	OBJECT	set[MAXSET+1];
	register int  i;
	register long  m;

	if ((set[0] = ogetint(objsize)) > MAXSET)
		octerror(USER, "bad set in getfullnode");
	for (i = 1; i <= set[0]; i++) {
		m = ogetint(objsize) + objorig;
		if ((set[i] = m) != m)
			octerror(USER, "too many objects");
	}
	return(fullnode(set));
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
	

static OCTREE
gettree()			/* get a pre-ordered octree */
{
	register OCTREE	 ot;
	register int  i;
	
	switch (getc(infp)) {
	case OT_EMPTY:
		return(EMPTY);
	case OT_FULL:
		return(getfullnode());
	case OT_TREE:
		if ((ot = octalloc()) == EMPTY)
			octerror(SYSTEM, "out of tree space in gettree");
		for (i = 0; i < 8; i++)
			octkid(ot, i) = gettree();
		return(ot);
	case EOF:
		octerror(USER, "truncated octree");
	default:
		octerror(USER, "damaged octree");
	}
	return NULL; /* pro forma return */
}


static int
nonsurfinset(os)			/* check set for modifier */
register OBJECT  *os;
{
	register OBJECT  s;
	register int  i;

	for (i = *os; i-- > 0; )
		if ((s = *++os) >= objorig && s < objorig+fnobjects &&
				ismodifier(objptr(s)->otype))
			return(1);
	return(0);
}


static void
skiptree(void)				/* skip octree on input */
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


static void
octerror(etyp, msg)			/* octree error */
int  etyp;
char  *msg;
{
	char  msgbuf[128];

	sprintf(msgbuf, "(%s): %s", infn, msg);
	error(etyp, msgbuf);
}
