#ifndef lint
static const char	RCSid[] = "$Id: func.c,v 2.20 2003/08/04 22:37:53 greg Exp $";
#endif
/*
 *  func.c - interface to calcomp functions.
 */

#include "copyright.h"

#include  "ray.h"

#include  "otypes.h"

#include  "func.h"


#define  INITFILE	"rayinit.cal"
#define  CALSUF		".cal"
#define  LCALSUF	4
char  REFVNAME[] = "`FILE_REFCNT";

XF  unitxf = {			/* identity transform */
	{{1.0, 0.0, 0.0, 0.0},
	{0.0, 1.0, 0.0, 0.0},
	{0.0, 0.0, 1.0, 0.0},
	{0.0, 0.0, 0.0, 1.0}},
	1.0
};

XF  funcxf;			/* current transformation */
static OBJREC  *fobj = NULL;	/* current function object */
static RAY  *fray = NULL;	/* current function ray */

static double  l_erf(char *), l_erfc(char *), l_arg(char *);


MFUNC *
getfunc(m, ff, ef, dofwd)	/* get function for this modifier */
OBJREC  *m;
int  ff;
unsigned int  ef;
int  dofwd;
{
	static char  initfile[] = INITFILE;
	char  sbuf[MAXSTR];
	register char  **arg;
	register MFUNC  *f;
	int  ne, na;
	register int  i;
					/* check to see if done already */
	if ((f = (MFUNC *)m->os) != NULL)
		return(f);
	fobj = NULL; fray = NULL;
	if (initfile[0]) {		/* initialize on first call */
		esupport |= E_VARIABLE|E_FUNCTION|E_INCHAN|E_RCONST|E_REDEFW;
		esupport &= ~(E_OUTCHAN);
		setcontext("");
		scompile("Dx=$1;Dy=$2;Dz=$3;", NULL, 0);
		scompile("Nx=$4;Ny=$5;Nz=$6;", NULL, 0);
		scompile("Px=$7;Py=$8;Pz=$9;", NULL, 0);
		scompile("T=$10;Ts=$25;Rdot=$11;", NULL, 0);
		scompile("S=$12;Tx=$13;Ty=$14;Tz=$15;", NULL, 0);
		scompile("Ix=$16;Iy=$17;Iz=$18;", NULL, 0);
		scompile("Jx=$19;Jy=$20;Jz=$21;", NULL, 0);
		scompile("Kx=$22;Ky=$23;Kz=$24;", NULL, 0);
		scompile("Lu=$26;Lv=$27;", NULL, 0);
		funset("arg", 1, '=', l_arg);
		funset("erf", 1, ':', l_erf);
		funset("erfc", 1, ':', l_erfc);
		setnoisefuncs();
		setprismfuncs();
		loadfunc(initfile);
		initfile[0] = '\0';
	}
	if ((na = m->oargs.nsargs) <= ff)
		goto toofew;
	arg = m->oargs.sarg;
	if ((f = (MFUNC *)calloc(1, sizeof(MFUNC))) == NULL)
		goto memerr;
	i = strlen(arg[ff]);			/* set up context */
	if (i == 1 && arg[ff][0] == '.')
		setcontext(f->ctx = "");	/* "." means no file */
	else {
		strcpy(sbuf,arg[ff]);		/* file name is context */
		if (i > LCALSUF && !strcmp(sbuf+i-LCALSUF, CALSUF))
			sbuf[i-LCALSUF] = '\0';	/* remove suffix */
		setcontext(f->ctx = savestr(sbuf));
		if (!vardefined(REFVNAME)) {	/* file loaded? */
			loadfunc(arg[ff]);
			varset(REFVNAME, '=', 1.0);
		} else				/* reference_count++ */
			varset(REFVNAME, '=', varvalue(REFVNAME)+1.0);
	}
	curfunc = NULL;			/* parse expressions */
	sprintf(sbuf, "%s \"%s\"", ofun[m->otype].funame, m->oname);
	for (i=0, ne=0; ef && i < na; i++, ef>>=1)
		if (ef & 1) {			/* flagged as an expression? */
			if (ne >= MAXEXPR)
				objerror(m, INTERNAL, "too many expressions");
			initstr(arg[i], sbuf, 0);
			f->ep[ne++] = getE1();
			if (nextc != EOF)
				syntax("unexpected character");
		}
	if (ef)
		goto toofew;
	if (i <= ff)			/* find transform args */
		i = ff+1;
	while (i < na && arg[i][0] != '-')
		i++;
	if (i == na)			/* no transform */
		f->f = f->b = &unitxf;
	else {				/* get transform */
		if ((f->b = (XF *)malloc(sizeof(XF))) == NULL)
			goto memerr;
		if (invxf(f->b, na-i, arg+i) != na-i)
			objerror(m, USER, "bad transform");
		if (f->b->sca < 0.0)
			f->b->sca = -f->b->sca;
		if (dofwd) {			/* do both transforms */
			if ((f->f = (XF *)malloc(sizeof(XF))) == NULL)
				goto memerr;
			xf(f->f, na-i, arg+i);
			if (f->f->sca < 0.0)
				f->f->sca = -f->f->sca;
		}
	}
	m->os = (char *)f;
	return(f);
toofew:
	objerror(m, USER, "too few string arguments");
memerr:
	error(SYSTEM, "out of memory in getfunc");
}


void
freefunc(m)			/* free memory associated with modifier */
OBJREC  *m;
{
	register MFUNC  *f;
	register int  i;

	if ((f = (MFUNC *)m->os) == NULL)
		return;
	for (i = 0; f->ep[i] != NULL; i++)
		epfree(f->ep[i]);
	if (f->ctx[0]) {			/* done with definitions */
		setcontext(f->ctx);
		i = varvalue(REFVNAME)-.5;	/* reference_count-- */
		if (i > 0)
			varset(REFVNAME, '=', (double)i);
		else
			dcleanup(2);		/* remove definitions */
		freestr(f->ctx);
	}
	if (f->b != &unitxf)
		free((void *)f->b);
	if (f->f != NULL && f->f != &unitxf)
		free((void *)f->f);
	free((void *)f);
	m->os = NULL;
}


int
setfunc(m, r)			/* set channels for function call */
OBJREC  *m;
register RAY  *r;
{
	static unsigned long  lastrno = ~0;
	register MFUNC  *f;
					/* get function */
	if ((f = (MFUNC *)m->os) == NULL)
		objerror(m, CONSISTENCY, "setfunc called before getfunc");
					/* set evaluator context */
	setcontext(f->ctx);
					/* check to see if matrix set */
	if (m == fobj && r->rno == lastrno)
		return(0);
	fobj = m;
	fray = r;
	if (r->rox != NULL)
		if (f->b != &unitxf) {
			funcxf.sca = r->rox->b.sca * f->b->sca;
			multmat4(funcxf.xfm, r->rox->b.xfm, f->b->xfm);
		} else
			funcxf = r->rox->b;
	else
		funcxf = *(f->b);
	lastrno = r->rno;
	eclock++;		/* notify expression evaluator */
	return(1);
}


void
loadfunc(fname)			/* load definition file */
char  *fname;
{
	char  *ffname;

	if ((ffname = getpath(fname, getrlibpath(), R_OK)) == NULL) {
		sprintf(errmsg, "cannot find function file \"%s\"", fname);
		error(USER, errmsg);
	}
	fcompile(ffname);
}


static double
l_arg(char *nm)			/* return nth real argument */
{
	register int  n;

	if (fobj == NULL)
		syntax("arg(n) used in constant expression");

	n = argument(1) + .5;		/* round to integer */

	if (n < 1)
		return(fobj->oargs.nfargs);

	if (n > fobj->oargs.nfargs) {
		sprintf(errmsg, "missing real argument %d", n);
		objerror(fobj, USER, errmsg);
	}
	return(fobj->oargs.farg[n-1]);
}


static double
l_erf(char *nm)			/* error function */
{
	extern double  erf();

	return(erf(argument(1)));
}


static double
l_erfc(char *nm)		/* cumulative error function */
{
	extern double  erfc();

	return(erfc(argument(1)));
}


double
chanvalue(n)			/* return channel n to calcomp */
register int  n;
{
	if (fray == NULL)
		syntax("ray parameter used in constant expression");

	if (--n < 0)
		goto badchan;

	if (n <= 2)			/* ray direction */

		return( (	fray->rdir[0]*funcxf.xfm[0][n] +
				fray->rdir[1]*funcxf.xfm[1][n] +
				fray->rdir[2]*funcxf.xfm[2][n]	)
			 / funcxf.sca );

	if (n <= 5)			/* surface normal */

		return( (	fray->ron[0]*funcxf.xfm[0][n-3] +
				fray->ron[1]*funcxf.xfm[1][n-3] +
				fray->ron[2]*funcxf.xfm[2][n-3]	)
			 / funcxf.sca );

	if (n <= 8)			/* intersection */

		return( fray->rop[0]*funcxf.xfm[0][n-6] +
				fray->rop[1]*funcxf.xfm[1][n-6] +
				fray->rop[2]*funcxf.xfm[2][n-6] +
					     funcxf.xfm[3][n-6] );

	if (n == 9)			/* total distance */
		return(raydist(fray,PRIMARY) * funcxf.sca);

	if (n == 10)			/* dot product (range [-1,1]) */
		return(	fray->rod <= -1.0 ? -1.0 :
			fray->rod >= 1.0 ? 1.0 :
			fray->rod );

	if (n == 11)			/* scale */
		return(funcxf.sca);

	if (n <= 14)			/* origin */
		return(funcxf.xfm[3][n-12]);

	if (n <= 17)			/* i unit vector */
		return(funcxf.xfm[0][n-15] / funcxf.sca);

	if (n <= 20)			/* j unit vector */
		return(funcxf.xfm[1][n-18] / funcxf.sca);

	if (n <= 23)			/* k unit vector */
		return(funcxf.xfm[2][n-21] / funcxf.sca);

	if (n == 24)			/* single ray (shadow) distance */
		return((fray->rot+raydist(fray->parent,SHADOW)) * funcxf.sca);

	if (n <= 26)			/* local (u,v) coordinates */
		return(fray->uv[n-25]);
badchan:
	error(USER, "illegal channel number");
	return(0.0);
}
