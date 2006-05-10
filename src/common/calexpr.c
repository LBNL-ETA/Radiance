#ifndef lint
static const char	RCSid[] = "$Id: calexpr.c,v 2.31 2006/05/10 15:21:20 greg Exp $";
#endif
/*
 *  Compute data values using expression parser
 *
 *  7/1/85  Greg Ward
 *
 *  11/11/85  Made channel input conditional with (INCHAN) compiles.
 *
 *  4/2/86  Added conditional compiles for function definitions (FUNCTION).
 *
 *  1/29/87  Made variables conditional (VARIABLE)
 *
 *  5/19/88  Added constant subexpression elimination (RCONST)
 *
 *  2/19/03	Eliminated conditional compiles in favor of esupport extern.
 */

#include "copyright.h"

#include  <stdio.h>
#include  <string.h>
#include  <ctype.h>
#include  <errno.h>
#include  <math.h>
#include  <stdlib.h>

#include  "rtmisc.h"
#include  "rtio.h"
#include  "rterror.h"
#include  "calcomp.h"

#define	 MAXLINE	256		/* maximum line length */

#define	 newnode()	(EPNODE *)ecalloc(1, sizeof(EPNODE))

#define	 isdecimal(c)	(isdigit(c) || (c) == '.')

static double  euminus(EPNODE *), eargument(EPNODE *), enumber(EPNODE *);
static double  echannel(EPNODE *);
static double  eadd(EPNODE *), esubtr(EPNODE *),
               emult(EPNODE *), edivi(EPNODE *),
               epow(EPNODE *);
static double  ebotch(EPNODE *);

unsigned int  esupport =		/* what to support */
		E_VARIABLE | E_FUNCTION ;

int  eofc = 0;				/* optional end-of-file character */
int  nextc;				/* lookahead character */

double	(*eoper[])(EPNODE *) = {	/* expression operations */
	ebotch,
	evariable,
	enumber,
	euminus,
	echannel,
	efunc,
	eargument,
	ebotch,
	ebotch,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	emult,
	eadd,
	0,
	esubtr,
	0,
	edivi,
	0,0,0,0,0,0,0,0,0,0,
	ebotch,
	0,0,
	ebotch,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	epow,
};

static FILE  *infp;			/* input file pointer */
static char  *linbuf;			/* line buffer */
static char  *infile;			/* input file name */
static int  lineno;			/* input line number */
static int  linepos;			/* position in buffer */


EPNODE *
eparse(			/* parse an expression string */
    char  *expr
)
{
    EPNODE  *ep;

    initstr(expr, NULL, 0);
    curfunc = NULL;
    ep = getE1();
    if (nextc != EOF)
	syntax("unexpected character");
    return(ep);
}


double
eval(			/* evaluate an expression string */
    char  *expr
)
{
    register EPNODE  *ep;
    double  rval;

    ep = eparse(expr);
    rval = evalue(ep);
    epfree(ep);
    return(rval);
}


int
epcmp(			/* compare two expressions for equivalence */
    register EPNODE  *ep1,
    register EPNODE  *ep2
)
{
	double  d;

	if (ep1->type != ep2->type)
		return(1);

	switch (ep1->type) {

	case VAR:
		return(ep1->v.ln != ep2->v.ln);

	case NUM:
		if (ep2->v.num == 0)
			return(ep1->v.num != 0);
		d = ep1->v.num / ep2->v.num;
		return((d > 1.000000000001) | (d < 0.999999999999));

	case CHAN:
	case ARG:
		return(ep1->v.chan != ep2->v.chan);

	case '=':
	case ':':
		return(epcmp(ep1->v.kid->sibling, ep2->v.kid->sibling));

	case TICK:
	case SYM:			/* should never get this one */
		return(0);

	default:
		ep1 = ep1->v.kid;
		ep2 = ep2->v.kid;
		while (ep1 != NULL) {
			if (ep2 == NULL)
				return(1);
			if (epcmp(ep1, ep2))
				return(1);
			ep1 = ep1->sibling;
			ep2 = ep2->sibling;
		}
		return(ep2 != NULL);
	}
}


void
epfree(			/* free a parse tree */
    register EPNODE	 *epar
)
{
    register EPNODE  *ep;

    switch (epar->type) {

	case VAR:
	    varfree(epar->v.ln);
	    break;
	    
	case SYM:
	    freestr(epar->v.name);
	    break;

	case NUM:
	case CHAN:
	case ARG:
	case TICK:
	    break;

	default:
	    while ((ep = epar->v.kid) != NULL) {
		epar->v.kid = ep->sibling;
		epfree(ep);
	    }
	    break;

    }

    efree((char *)epar);
}

				/* the following used to be a switch */
static double
eargument(
    EPNODE	*ep
)
{
    return(argument(ep->v.chan));
}

static double
enumber(
    EPNODE	*ep
)
{
    return(ep->v.num);
}

static double
euminus(
    EPNODE	*ep
)
{
    register EPNODE  *ep1 = ep->v.kid;

    return(-evalue(ep1));
}

static double
echannel(
    EPNODE	*ep
)
{
    return(chanvalue(ep->v.chan));
}

static double
eadd(
    EPNODE	*ep
)
{
    register EPNODE  *ep1 = ep->v.kid;

    return(evalue(ep1) + evalue(ep1->sibling));
}

static double
esubtr(
    EPNODE	*ep
)
{
    register EPNODE  *ep1 = ep->v.kid;

    return(evalue(ep1) - evalue(ep1->sibling));
}

static double
emult(
    EPNODE	*ep
)
{
    register EPNODE  *ep1 = ep->v.kid;

    return(evalue(ep1) * evalue(ep1->sibling));
}

static double
edivi(
    EPNODE	*ep
)
{
    register EPNODE  *ep1 = ep->v.kid;
    double  d;

    d = evalue(ep1->sibling);
    if (d == 0.0) {
	wputs("Division by zero\n");
	errno = ERANGE;
	return(0.0);
    }
    return(evalue(ep1) / d);
}

static double
epow(
    EPNODE	*ep
)
{
    register EPNODE  *ep1 = ep->v.kid;
    double  d;
    int	 lasterrno;

    lasterrno = errno;
    errno = 0;
    d = pow(evalue(ep1), evalue(ep1->sibling));
#ifdef  isnan
    if (errno == 0)
	if (isnan(d))
	    errno = EDOM;
	else if (isinf(d))
	    errno = ERANGE;
#endif
    if (errno == EDOM || errno == ERANGE) {
	wputs("Illegal power\n");
	return(0.0);
    }
    errno = lasterrno;
    return(d);
}

static double
ebotch(
    EPNODE	*ep
)
{
    eputs("Bad expression!\n");
    quit(1);
	return 0.0; /* pro forma return */
}


EPNODE *
ekid(			/* return pointer to a node's nth kid */
    register EPNODE	 *ep,
    register int  n
)
{

    for (ep = ep->v.kid; ep != NULL; ep = ep->sibling)
	if (--n < 0)
	    break;

    return(ep);
}


int
nekids(			/* return # of kids for node ep */
    register EPNODE	 *ep
)
{
    register int  n = 0;

    for (ep = ep->v.kid; ep != NULL; ep = ep->sibling)
	n++;

    return(n);
}


void
initfile(		/* prepare input file */
    FILE  *fp,
    char  *fn,
    int  ln
)
{
    static char	 inpbuf[MAXLINE];

    infp = fp;
    linbuf = inpbuf;
    infile = fn;
    lineno = ln;
    linepos = 0;
    inpbuf[0] = '\0';
    scan();
}


void
initstr(		/* prepare input string */
    char  *s,
    char  *fn,
    int  ln
)
{
    infp = NULL;
    infile = fn;
    lineno = ln;
    linbuf = s;
    linepos = 0;
    scan();
}


void
getscanpos(	/* return current scan position */
    char  **fnp,
    int  *lnp,
    char  **spp,
    FILE  **fpp
)
{
    if (fnp != NULL) *fnp = infile;
    if (lnp != NULL) *lnp = lineno;
    if (spp != NULL) *spp = linbuf+linepos;
    if (fpp != NULL) *fpp = infp;
}


int
scan(void)		/* scan next character, return literal next */
{
    register int  lnext = 0;

    do {
	if (linbuf[linepos] == '\0')
	    if (infp == NULL || fgets(linbuf, MAXLINE, infp) == NULL)
		nextc = EOF;
	    else {
		nextc = linbuf[0];
		lineno++;
		linepos = 1;
	    }
	else
	    nextc = linbuf[linepos++];
	if (!lnext)
		lnext = nextc;
	if (nextc == eofc) {
		nextc = EOF;
		break;
	}
	if (nextc == '{') {
	    scan();
	    while (nextc != '}')
		if (nextc == EOF)
		    syntax("'}' expected");
		else
		    scan();
	    scan();
	}
    } while (isspace(nextc));
    return(lnext);
}


char *
long2ascii(			      /* convert long to ascii */
    long  l
)
{
    static char	 buf[16];
    register char  *cp;
    int	 neg = 0;

    if (l == 0)
	return("0");
    if (l < 0) {
	l = -l;
	neg++;
    }
    cp = buf + sizeof(buf);
    *--cp = '\0';
    while (l) {
	*--cp = l % 10 + '0';
	l /= 10;
    }
    if (neg)
	*--cp = '-';
    return(cp);
}


void
syntax(			/* report syntax error and quit */
    char  *err
)
{
    register int  i;

    if (infile != NULL || lineno != 0) {
	if (infile != NULL) eputs(infile);
	if (lineno != 0) {
	    eputs(infile != NULL ? ", line " : "line ");
	    eputs(long2ascii((long)lineno));
	}
	eputs(":\n");
    }
    eputs(linbuf);
    if (linbuf[strlen(linbuf)-1] != '\n')
	eputs("\n");
    for (i = 0; i < linepos-1; i++)
	eputs(linbuf[i] == '\t' ? "\t" : " ");
    eputs("^ ");
    eputs(err);
    eputs("\n");
    quit(1);
}


void
addekid(			/* add a child to ep */
    register EPNODE	 *ep,
    EPNODE	*ekid
)
{
    if (ep->v.kid == NULL)
	ep->v.kid = ekid;
    else {
	for (ep = ep->v.kid; ep->sibling != NULL; ep = ep->sibling)
	    ;
	ep->sibling = ekid;
    }
    ekid->sibling = NULL;
}


char *
getname(void)			/* scan an identifier */
{
    static char	 str[RMAXWORD+1];
    register int  i, lnext;

    lnext = nextc;
    for (i = 0; i < RMAXWORD && isid(lnext); i++, lnext = scan())
	str[i] = lnext;
    str[i] = '\0';
    while (isid(lnext))		/* skip rest of name */
	lnext = scan();

    return(str);
}


int
getinum(void)			/* scan a positive integer */
{
    register int  n, lnext;

    n = 0;
    lnext = nextc;
    while (isdigit(lnext)) {
	n = n * 10 + lnext - '0';
	lnext = scan();
    }
    return(n);
}


double
getnum(void)			/* scan a positive float */
{
    register int  i, lnext;
    char  str[RMAXWORD+1];

    i = 0;
    lnext = nextc;
    while (isdigit(lnext) && i < RMAXWORD) {
	str[i++] = lnext;
	lnext = scan();
    }
    if (lnext == '.' && i < RMAXWORD) {
	str[i++] = lnext;
	lnext = scan();
	if (i == 1 && !isdigit(lnext))
	    syntax("badly formed number");
	while (isdigit(lnext) && i < RMAXWORD) {
	    str[i++] = lnext;
	    lnext = scan();
	}
    }
    if ((lnext == 'e') | (lnext == 'E') && i < RMAXWORD) {
	str[i++] = lnext;
	lnext = scan();
	if ((lnext == '-') | (lnext == '+') && i < RMAXWORD) {
	    str[i++] = lnext;
	    lnext = scan();
	}
	if (!isdigit(lnext))
	    syntax("missing exponent");
	while (isdigit(lnext) && i < RMAXWORD) {
	    str[i++] = lnext;
	    lnext = scan();
	}
    }
    str[i] = '\0';

    return(atof(str));
}


EPNODE *
getE1(void)				/* E1 -> E1 ADDOP E2 */
				/*	 E2 */
{
    register EPNODE  *ep1, *ep2;

    ep1 = getE2();
    while (nextc == '+' || nextc == '-') {
	ep2 = newnode();
	ep2->type = nextc;
	scan();
	addekid(ep2, ep1);
	addekid(ep2, getE2());
	if (esupport&E_RCONST &&
			ep1->type == NUM && ep1->sibling->type == NUM)
		ep2 = rconst(ep2);
	ep1 = ep2;
    }
    return(ep1);
}


EPNODE *
getE2(void)				/* E2 -> E2 MULOP E3 */
				/*	 E3 */
{
    register EPNODE  *ep1, *ep2;

    ep1 = getE3();
    while (nextc == '*' || nextc == '/') {
	ep2 = newnode();
	ep2->type = nextc;
	scan();
	addekid(ep2, ep1);
	addekid(ep2, getE3());
	if (esupport&E_RCONST &&
			ep1->type == NUM && ep1->sibling->type == NUM)
		ep2 = rconst(ep2);
	ep1 = ep2;
    }
    return(ep1);
}


EPNODE *
getE3(void)				/* E3 -> E4 ^ E3 */
				/*	 E4 */
{
    register EPNODE  *ep1, *ep2;

    ep1 = getE4();
    if (nextc == '^') {
	ep2 = newnode();
	ep2->type = nextc;
	scan();
	addekid(ep2, ep1);
	addekid(ep2, getE3());
	if (esupport&E_RCONST &&
			ep1->type == NUM && ep1->sibling->type == NUM)
		ep2 = rconst(ep2);
	return(ep2);
    }
    return(ep1);
}


EPNODE *
getE4(void)				/* E4 -> ADDOP E5 */
				/*	 E5 */
{
    register EPNODE  *ep1, *ep2;

    if (nextc == '-') {
	scan();
	ep2 = getE5();
	if (ep2->type == NUM) {
		ep2->v.num = -ep2->v.num;
		return(ep2);
	}
	if (ep2->type == UMINUS) {	/* don't generate -(-E5) */
	    efree((char *)ep2);
	    return(ep2->v.kid);
	}
	ep1 = newnode();
	ep1->type = UMINUS;
	addekid(ep1, ep2);
	return(ep1);
    }
    if (nextc == '+')
	scan();
    return(getE5());
}


EPNODE *
getE5(void)			/* E5 -> (E1) */
				/*	 VAR */
				/*	 NUM */
				/*	 $N */
				/*	 FUNC(E1,..) */
				/*	 ARG */
{
	int	 i;
	char  *nam;
	register EPNODE  *ep1, *ep2;

	if (nextc == '(') {
		scan();
		ep1 = getE1();
		if (nextc != ')')
			syntax("')' expected");
		scan();
		return(ep1);
	}

	if (esupport&E_INCHAN && nextc == '$') {
		scan();
		ep1 = newnode();
		ep1->type = CHAN;
		ep1->v.chan = getinum();
		return(ep1);
	}

	if (esupport&(E_VARIABLE|E_FUNCTION) &&
			(isalpha(nextc) || nextc == CNTXMARK)) {
		nam = getname();
		ep1 = NULL;
		if ((esupport&(E_VARIABLE|E_FUNCTION)) == (E_VARIABLE|E_FUNCTION)
				&& curfunc != NULL)
			for (i = 1, ep2 = curfunc->v.kid->sibling;
					ep2 != NULL; i++, ep2 = ep2->sibling)
				if (!strcmp(ep2->v.name, nam)) {
					ep1 = newnode();
					ep1->type = ARG;
					ep1->v.chan = i;
					break;
				}
		if (ep1 == NULL) {
			ep1 = newnode();
			ep1->type = VAR;
			ep1->v.ln = varinsert(nam);
		}
		if (esupport&E_FUNCTION && nextc == '(') {
			ep2 = newnode();
			ep2->type = FUNC;
			addekid(ep2, ep1);
			ep1 = ep2;
			do {
				scan();
				addekid(ep1, getE1());
			} while (nextc == ',');
			if (nextc != ')')
				syntax("')' expected");
			scan();
		} else if (!(esupport&E_VARIABLE))
			syntax("'(' expected");
		if (esupport&E_RCONST && isconstvar(ep1))
			ep1 = rconst(ep1);
		return(ep1);
	}

	if (isdecimal(nextc)) {
		ep1 = newnode();
		ep1->type = NUM;
		ep1->v.num = getnum();
		return(ep1);
	}
	syntax("unexpected character");
	return NULL; /* pro forma return */
}


EPNODE *
rconst(			/* reduce a constant expression */
    register EPNODE	 *epar
)
{
    register EPNODE  *ep;

    ep = newnode();
    ep->type = NUM;
    errno = 0;
    ep->v.num = evalue(epar);
    if (errno == EDOM || errno == ERANGE)
	syntax("bad constant expression");
    epfree(epar);
 
    return(ep);
}


int
isconstvar(			/* is ep linked to a constant expression? */
    register EPNODE	 *ep
)
{
    register EPNODE  *ep1;

    if (esupport&E_FUNCTION && ep->type == FUNC) {
	if (!isconstfun(ep->v.kid))
		return(0);
	for (ep1 = ep->v.kid->sibling; ep1 != NULL; ep1 = ep1->sibling)
	    if (ep1->type != NUM && !isconstfun(ep1))
		return(0);
	return(1);
    }
    if (ep->type != VAR)
	return(0);
    ep1 = ep->v.ln->def;
    if (ep1 == NULL || ep1->type != ':')
	return(0);
    if (esupport&E_FUNCTION && ep1->v.kid->type != SYM)
	return(0);
    return(1);
}


int
isconstfun(			/* is ep linked to a constant function? */
    register EPNODE	 *ep
)
{
    register EPNODE  *dp;
    register LIBR  *lp;

    if (ep->type != VAR)
	return(0);
    if ((dp = ep->v.ln->def) != NULL) {
	if (dp->v.kid->type == FUNC)
	    return(dp->type == ':');
	else
	    return(0);		/* don't identify masked library functions */
    }
    if ((lp = ep->v.ln->lib) != NULL)
	return(lp->atyp == ':');
    return(0);
}
