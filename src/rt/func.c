#ifndef lint
static const char	RCSid[] = "$Id$";
#endif
/*
 *  func.c - interface to calcomp functions.
 */

#include "copyright.h"

#include  "ray.h"
#include  "paths.h"
#include  "otypes.h"
#include  "func.h"
#include <ctype.h>


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

static char  rayinitcal[] = INITFILE;

static double  l_erf(char *), l_erfc(char *), l_arg(char *),
		l_source_corr(char *), l_source_angle(char *);


void
initfunc(void)	/* initialize function evaluation */
{
	if (!rayinitcal[0])	/* already done? */
		return;
	esupport |= E_VARIABLE|E_FUNCTION|E_INCHAN|E_RCONST|E_REDEFW;
	esupport &= ~(E_OUTCHAN);
	calcontext("");
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
	funset("source_corr", 2, '=', l_source_corr);
	funset("source_theta", 1, '=', l_source_angle);
	funset("source_phi", 1, '=', l_source_angle);
	setnoisefuncs();
	setprismfuncs();
	loadfunc(rayinitcal);
	rayinitcal[0] = '\0';
	doptimize(1);		/* optimize definitions */
}


/* Set parameters for current evaluation */
void
set_eparams(const char *prms)
{
	static const char	*last_params = NULL;
	char			vname[RMAXWORD];
	double			value;
	char			*cpd;
					/* check if already set */
	if (prms == NULL || !*prms)
		return;
	if (prms == last_params || (last_params != NULL &&
					!strcmp(prms, last_params)))
		return;
	last_params = prms;		/* XXX assumes static string */
					/* assign each variable */
	while (*prms) {
		if (isspace(*prms)) {
			++prms; continue;
		}
		if (!isalpha(*prms))
			goto bad_params;
		cpd = vname;
		while (*prms && (*prms != '=') & !isspace(*prms)) {
			if (!isid(*prms) | (cpd-vname >= RMAXWORD-1))
				goto bad_params;
			*cpd++ = *prms++;
		}
		*cpd = '\0';
		while (isspace(*prms)) prms++;
		if (*prms++ != '=')
			goto bad_params;
		value = atof(prms);
		if ((prms = fskip((char *)prms)) == NULL)
			goto bad_params;
		while (isspace(*prms)) prms++;
		prms += (*prms == ',') | (*prms == ';') | (*prms == ':');
		varset(vname, '=', value);
	}
	eclock++;		/* notify expression evaluator */
	return;
bad_params:
	sprintf(errmsg, "bad parameter list '%s'", last_params);
	error(USER, errmsg);
}


MFUNC *
getfunc(	/* get function for this modifier */
	OBJREC  *m,
	int  ff,
	unsigned int  ef,
	int  dofwd
)
{
	char  sbuf[MAXSTR];
	char  **arg;
	MFUNC  *f;
	int  ne, na;
	int  i;
					/* check to see if done already */
	if ((f = (MFUNC *)m->os) != NULL)
		return(f);
	fobj = NULL; fray = NULL;
	if ((na = m->oargs.nsargs) <= ff)
		goto toofew;
	arg = m->oargs.sarg;
	if ((f = (MFUNC *)calloc(1, sizeof(MFUNC))) == NULL)
		goto memerr;
	i = strlen(arg[ff]);			/* set up context */
	if (i == 1 && arg[ff][0] == '.') {
		calcontext(f->ctx = "");	/* "." means no file */
	} else {
		strcpy(sbuf,arg[ff]);		/* file name is context */
		if (i > LCALSUF && !strcasecmp(sbuf+i-LCALSUF, CALSUF))
			sbuf[i-LCALSUF] = '\0';	/* remove suffix */
		calcontext(f->ctx = savestr(sbuf));
		if (!vardefined(REFVNAME)) {	/* file loaded? */
			loadfunc(arg[ff]);
			varset(REFVNAME, '=', 1.0);
		} else				/* reference_count++ */
			varset(REFVNAME, '=', varvalue(REFVNAME)+1.0);
	}
	ecurfunc = NULL;			/* parse expressions */
	sprintf(sbuf, "%s \"%s\"", ofun[m->otype].funame, m->oname);
	for (i=0, ne=0; ef && i < na; i++, ef>>=1)
		if (ef & 1) {			/* flagged as an expression? */
			if (ne >= MAXEXPR)
				objerror(m, INTERNAL, "too many expressions");
			initstr(arg[i], sbuf, 0);
			f->ep[ne++] = getE1();
			if (nextc != EOF)
				esyntax("unexpected character");
		}
	if (ef)
		goto toofew;
	if (i <= ff)			/* find transform args */
		i = ff+1;
	while (i < na && !isxfopt(arg[i]))
		i++;
	if (i == na) {			/* no transform */
		f->fxp = f->bxp = &unitxf;
	} else {			/* get transform */
		if ((f->bxp = (XF *)malloc(sizeof(XF))) == NULL)
			goto memerr;
		if (invxf(f->bxp, na-i, arg+i) != na-i)
			objerror(m, USER, "bad transform");
		if (f->bxp->sca < 0.0)
			f->bxp->sca = -f->bxp->sca;
		if (dofwd) {			/* do both transforms */
			if ((f->fxp = (XF *)malloc(sizeof(XF))) == NULL)
				goto memerr;
			xf(f->fxp, na-i, arg+i);
			if (f->fxp->sca < 0.0)
				f->fxp->sca = -f->fxp->sca;
		}
	}
	m->os = (char *)f;
	return(f);
toofew:
	objerror(m, USER, "too few string arguments");
memerr:
	error(SYSTEM, "out of memory in getfunc");
	return NULL; /* pro forma return */
}


void
freefunc(			/* free memory associated with modifier */
	OBJREC  *m
)
{
	MFUNC  *f;
	int  i;

	if ((f = (MFUNC *)m->os) == NULL)
		return;
	for (i = 0; f->ep[i] != NULL; i++)
		epfree(f->ep[i],1);
	if (f->ctx[0]) {			/* done with definitions */
		calcontext(f->ctx);
		i = varvalue(REFVNAME)-.5;	/* reference_count-- */
		if (i > 0)
			varset(REFVNAME, '=', (double)i);
		else
			dcleanup(2);		/* remove definitions */
		freestr(f->ctx);
	}
	if (f->bxp != &unitxf)
		free((void *)f->bxp);
	if ((f->fxp != NULL) & (f->fxp != &unitxf))
		free((void *)f->fxp);
	free((void *)f);
	m->os = NULL;
}


int
setfunc(			/* set channels for function call */
	OBJREC	*m,
	RAY	*r
)
{
	static RNUMBER	lastrno = ~0;
	MFUNC		*f;
					/* get function if any */
	if ((f = (MFUNC *)m->os) == NULL)
		objerror(m, CONSISTENCY, "setfunc called before getfunc");
		
	calcontext(f->ctx);		/* set evaluator context */
					/* check to see if matrix set */
	if ((m == fobj) & (r->rno == lastrno))
		return(0);
	fobj = m;
	fray = r;
	if (r->rox != NULL) {
		if (f->bxp != &unitxf) {
			funcxf.sca = r->rox->b.sca * f->bxp->sca;
			multmat4(funcxf.xfm, r->rox->b.xfm, f->bxp->xfm);
		} else
			funcxf = r->rox->b;
	} else
		funcxf = *f->bxp;
	lastrno = r->rno;
	eclock++;		/* notify expression evaluator */
	return(1);
}


int
worldfunc(			/* special function context sans object */
	const char	*ctx,
	RAY	*r
)
{
	static RNUMBER	lastrno = ~0;
					/* set evaluator context */
	calcontext((char *)ctx);
					/* check if ray already set */
	if ((fobj == NULL) & (r->rno == lastrno))
		return(0);
	fobj = NULL;
	fray = r;
	funcxf = unitxf;
	lastrno = r->rno;
	eclock++;		/* notify expression evaluator */
	return(1);
}


void
loadfunc(			/* load definition file */
	char  *fname
)
{
	char  *ffname;

	if ((ffname = getpath(fname, getrlibpath(), R_OK)) == NULL) {
		sprintf(errmsg, "cannot find function file \"%s\"", fname);
		error(SYSTEM, errmsg);
	}
	fcompile(ffname);
}


static double
l_arg(char *nm)			/* return nth real argument */
{
	int  n;

	if (fobj == NULL)
		error(USER, "bad call to arg(n) - illegal constant in .cal file?");

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
	return(erf(argument(1)));
}


static double
l_erfc(char *nm)		/* cumulative error function */
{
	return(erfc(argument(1)));
}


static double
l_source_corr(char *nm)		/* photometry correction */
{
#define rarg	(fobj->oargs.farg)
#define nrargs	(fobj->oargs.nfargs)
	double	d, d1, rv;

	if ((fobj == NULL) | (fray == NULL))
		error(USER,
		    "bad call to source_corr() -- illegal constant in .cal file?");

	switch((int)(argument(2)+.5)) {
	case 0:			/* basic correction */
		return(nrargs > 0 ? rarg[0]*argument(1) : argument(1));
	case 1:			/* flat emitter correction */
		return( (nrargs > 0 ? rarg[0]*argument(1) : argument(1))
				/ fray->rod );
	case 2:			/* box source correction */
		if (nrargs < 4)
			goto notenough;
		return( argument(1) * rarg[0] /
			(fabs(chanvalue(1))*rarg[2]*rarg[3] +
			 fabs(chanvalue(2))*rarg[1]*rarg[3] +
			 fabs(chanvalue(3))*rarg[1]*rarg[2]) );
	case 3:			/* local box source correction */
		if (nrargs < 4)
			goto notenough;
		d = chanvalue(25);	/* Ts (shadow distance) */
		d1 = fabs(chanvalue(7) - chanvalue(1)*d) - rarg[1]*.5;
		rv = d1*rarg[2]*rarg[3]*(d1 > 0);
		d1 = fabs(chanvalue(8) - chanvalue(2)*d) - rarg[2]*.5;
		rv += d1*rarg[1]*rarg[3]*(d1 > 0);
		d1 = fabs(chanvalue(9) - chanvalue(3)*d) - rarg[3]*.5;
		rv += d1*rarg[1]*rarg[2]*(d1 > 0);
		return( argument(1) * rarg[0] / rv );
	case 4:			/* cylindrical source correction */
		if (nrargs < 3)
			goto notenough;
		d = chanvalue(3);	/* Dz */
		if (d >= 1.) d = 1.-FTINY;
		rv = rarg[1]*rarg[2]*sqrt(1. - d*d) +
				.25*PI*rarg[1]*rarg[1]*fabs(d);
		return( argument(1) * rarg[0] / rv );
	}
	objerror(fobj, USER, "illegal first argument in source_corr()");
notenough:
	objerror(fobj, USER, "missing real argument(s) for source_corr()");
	return(0.0);	/* pro forma */
#undef rarg
#undef nrargs
}


static double
l_source_angle(char *nm)	/* photometry angle */
{
	double	angle;

	if ((fray == NULL) | (fobj == NULL)) {
		sprintf(errmsg,
		    "bad call to %s() -- illegal constant in .cal file?", nm);
		error(USER, errmsg);
	}
	if (!strcmp(nm, "source_theta")) {
		double	Dz = chanvalue(3);
		if (Dz >= 1.) return(0.);
		if (Dz <= -1.) return(180.);
		return(180./PI * acos(Dz));
	}
	angle = 180./PI * atan2(-chanvalue(2), -chanvalue(1));
	angle += 360.*(angle < 0);

	switch ((int)(argument(1) + .5)) {
	case 90:		/* quad symmetry */
		angle -= 180.*(angle >= 180.);
		return(angle < 90. ? angle : 180.-angle);
	case 180:		/* bilateral symmetry */
		return(angle < 180. ? angle : 360.-angle);
	case 360:		/* no symmetry */
		return(angle);
	}
	objerror(fobj, USER, "bad call to source_phi()");
	return(0.0);
}


double
chanvalue(			/* return channel n to calcomp */
	int  n
)
{
	if (fray == NULL)
		esyntax("ray parameter used in constant expression");

	if (--n < 0)
		goto badchan;

	if (n < 3)			/* ray direction */

		return( (	fray->rdir[0]*funcxf.xfm[0][n] +
				fray->rdir[1]*funcxf.xfm[1][n] +
				fray->rdir[2]*funcxf.xfm[2][n]	)
			 / funcxf.sca );

	if (n < 6)			/* surface normal */

		return( (	fray->ron[0]*funcxf.xfm[0][n-3] +
				fray->ron[1]*funcxf.xfm[1][n-3] +
				fray->ron[2]*funcxf.xfm[2][n-3]	)
			 / funcxf.sca );

	if (n < 9)			/* intersection point */

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

	if (n < 15)			/* origin */
		return(funcxf.xfm[3][n-12]);

	if (n < 18)			/* i unit vector */
		return(funcxf.xfm[0][n-15] / funcxf.sca);

	if (n < 21)			/* j unit vector */
		return(funcxf.xfm[1][n-18] / funcxf.sca);

	if (n < 24)			/* k unit vector */
		return(funcxf.xfm[2][n-21] / funcxf.sca);

	if (n == 24)			/* single ray (shadow) distance */
		return((fray->rot+raydist(fray->parent,SHADOW)) * funcxf.sca);

	if (n < 27)			/* local (u,v) coordinates */
		return(fray->uv[n-25]);
badchan:
	error(USER, "illegal channel number");
	return(0.0);
}
