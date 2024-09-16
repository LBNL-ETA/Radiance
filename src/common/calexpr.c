#ifndef lint
static const char	RCSid[] = "$Id: calexpr.c,v 2.51 2024/09/16 17:31:14 greg Exp $";
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

#define	 isdecimal(c)	(isdigit(c) | ((c) == '.'))

#define  envalue(ep)	((ep)->type==NUM ? (ep)->v.num : evalue(ep))

static double  euminus(EPNODE *), enumber(EPNODE *);
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
    ecurfunc = NULL;
    ep = getE1();
    if (nextc != EOF)
	esyntax("unexpected character");
    return(ep);
}


double
eval(			/* evaluate an expression string */
    char  *expr
)
{
    int  prev_support = esupport;
    EPNODE  *ep;
    double  rval;

    esupport &= ~E_RCONST;	/* don't bother reducing constant expr */
    ep = eparse(expr);
    esupport = prev_support;	/* as you were */
    rval = evalue(ep);
    epfree(ep,1);
    return(rval);
}


int
epcmp(			/* compare two expressions for equivalence */
    EPNODE  *ep1,
    EPNODE  *ep2
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

	case CLKT:
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
    EPNODE	 *epar,
    int		frep
)
{
    EPNODE	*ep;

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
	case CLKT:
	    break;

	default:
	    if (epar->nkids < 0) {
	    	ep = epar->v.kid - epar->nkids;
	    	while (ep > epar->v.kid)
	    		epfree(--ep, 0);
		efree(ep);	/* free array space */
	    } else
	    	while ((ep = epar->v.kid) != NULL) {
		    epar->v.kid = ep->sibling;
		    epfree(ep, 1);
		}
	    break;

    }
    if (frep)
    	efree(epar);
    else
    	memset(epar, 0, sizeof(EPNODE));
}


static void
epflatten(			/* flatten hierarchies for '+', '*' */
	EPNODE *epar
)
{
    EPNODE	*ep;

    if (epar->nkids < 0)	/* can't handle array allocations */
    	return;

    for (ep = epar->v.kid; ep != NULL; ep = ep->sibling)
    	while (ep->type == epar->type && ep->nkids > 0) {
	    EPNODE	*ep1 = ep->v.kid;
	    while (ep1->sibling != NULL)
	    	ep1 = ep1->sibling;
	    ep1->sibling = ep->sibling;
	    epar->nkids += ep->nkids - 1;
	    ep1 = ep->v.kid;
	    *ep = *ep1;
	    efree(ep1);		/* not epfree()! */
	}
}


void
epoptimize(			/* flatten operations, lists -> arrays */
	EPNODE	*epar
)
{
    EPNODE	*ep;

    if ((epar->type == '+') | (epar->type == '*'))
    	epflatten(epar);	/* flatten associative operations */

    if (epar->nkids)		/* do children if any */
    	for (ep = epar->v.kid; ep != NULL; ep = ep->sibling)
	    epoptimize(ep);

    if (epar->nkids > 4) {	/* make list into array if > 4 kids */
        int	n = 1;
    	epar->v.kid = (EPNODE *)erealloc(epar->v.kid,
    					sizeof(EPNODE)*epar->nkids);
    	while (n < epar->nkids) {
	    ep = epar->v.kid[n-1].sibling;
	    epar->v.kid[n] = *ep;
	    efree(ep);		/* not epfree()! */
	    epar->v.kid[n-1].sibling = epar->v.kid + n;
	    n++;
	}
	epar->nkids = -n;
    }
}

				/* the following used to be a switch */
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
    EPNODE  *ep1 = ep->v.kid;

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
    double  sum = 0;
    EPNODE  *ep1 = ep->v.kid;

    do
    	sum += envalue(ep1);
    while ((ep1 = ep1->sibling) != NULL);

    return(sum);
}

static double
esubtr(
    EPNODE	*ep
)
{
    EPNODE  *ep1 = ep->v.kid;
    EPNODE  *ep2 = ep1->sibling;

    return(envalue(ep1) - envalue(ep2));
}

static double
emult(
    EPNODE	*ep
)
{
    double  prod = 1;
    EPNODE  *ep1 = ep->v.kid;

    do
    	prod *= envalue(ep1);
    while ((ep1 = ep1->sibling) != NULL);

    return(prod);
}

static double
edivi(
    EPNODE	*ep
)
{
    EPNODE  *ep1 = ep->v.kid;
    double  den = evalue(ep1->sibling);

    if (den == 0.0) {
	wputs("Division by zero\n");
	errno = ERANGE;
	return(0.0);
    }
    return(envalue(ep1) / den);
}

static double
epow(
    EPNODE	*ep
)
{
    EPNODE  *ep1 = ep->v.kid;
    double  d;
    int	 lasterrno;

    lasterrno = errno;
    errno = 0;
    d = pow(evalue(ep1), evalue(ep1->sibling));
#ifdef  isnan
    if (errno == 0) {
	if (isnan(d))
	    errno = EDOM;
	else if (isinf(d))
	    errno = ERANGE;
    }
#endif
    if ((errno == EDOM) | (errno == ERANGE)) {
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
    EPNODE	 *ep,
    int  n
)
{
    if (ep->nkids < 0) {	/* allocated array? */
    	if (n >= -ep->nkids)
	    return(NULL);
    	return(ep->v.kid + n);
    }
    ep = ep->v.kid;		/* else get from list */
    while (n-- > 0)
    	if ((ep = ep->sibling) == NULL)
		break;
    return(ep);
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
    escan();
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
    escan();
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
escan(void)		/* scan next character, return literal next */
{
    int  lnext = 0;

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
	    escan();
	    while (nextc != '}')
		if (nextc == EOF)
		    esyntax("'}' expected");
		else
		    escan();
	    escan();
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
    char  *cp;
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
esyntax(			/* report syntax error and quit */
    char  *err
)
{
    int  i;

    if ((infile != NULL) | (lineno != 0)) {
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
    EPNODE	 *ep,
    EPNODE	*ek
)
{
    if (ep->nkids < 0) {
    	eputs("Cannot add kid to EPNODE array\n");
    	quit(1);
    }
    ep->nkids++;
    if (ep->v.kid == NULL)
	ep->v.kid = ek;
    else {
	for (ep = ep->v.kid; ep->sibling != NULL; ep = ep->sibling)
	    ;
	ep->sibling = ek;
    }
    ek->sibling = NULL;		/* shouldn't be necessary */
}


char *
getname(void)			/* scan an identifier */
{
    static char	 str[RMAXWORD+1];
    int  i, lnext;

    lnext = nextc;
    for (i = 0; i < RMAXWORD && isid(lnext); i++, lnext = escan())
	str[i] = lnext;
    str[i] = '\0';
    while (isid(lnext))		/* skip rest of name */
	lnext = escan();

    return(str);
}


int
getinum(void)			/* scan a positive integer */
{
    int  n, lnext;

    n = 0;
    lnext = nextc;
    while (isdigit(lnext)) {
	n = n * 10 + lnext - '0';
	lnext = escan();
    }
    return(n);
}


double
getnum(void)			/* scan a positive float */
{
    int  i, lnext;
    char  str[RMAXWORD+1];

    i = 0;
    lnext = nextc;
    while (isdigit(lnext) && i < RMAXWORD) {
	str[i++] = lnext;
	lnext = escan();
    }
    if ((lnext == '.') & (i < RMAXWORD)) {
	str[i++] = lnext;
	lnext = escan();
	if (i == 1 && !isdigit(lnext))
	    esyntax("badly formed number");
	while (isdigit(lnext) && i < RMAXWORD) {
	    str[i++] = lnext;
	    lnext = escan();
	}
    }
    if ((lnext == 'e') | (lnext == 'E') && i < RMAXWORD) {
	str[i++] = lnext;
	lnext = escan();
	if ((lnext == '-') | (lnext == '+') && i < RMAXWORD) {
	    str[i++] = lnext;
	    lnext = escan();
	}
	if (!isdigit(lnext))
	    esyntax("missing exponent");
	while (isdigit(lnext) && i < RMAXWORD) {
	    str[i++] = lnext;
	    lnext = escan();
	}
    }
    str[i] = '\0';

    return(atof(str));
}


EPNODE *
getE1(void)			/* E1 -> E1 ADDOP E2 */
				/*	 E2 */
{
    EPNODE  *ep1, *ep2;

    ep1 = getE2();
    while ((nextc == '+') | (nextc == '-')) {
	ep2 = newnode();
	ep2->type = nextc;
	escan();
	addekid(ep2, ep1);
	addekid(ep2, getE2());
	if (esupport&E_RCONST &&
			(ep1->type == NUM) & (ep1->sibling->type == NUM))
		ep2 = rconst(ep2);
	ep1 = ep2;
    }
    return(ep1);
}


EPNODE *
getE2(void)			/* E2 -> E2 MULOP E3 */
				/*	 E3 */
{
    EPNODE  *ep1, *ep2;

    ep1 = getE3();
    while ((nextc == '*') | (nextc == '/')) {
	ep2 = newnode();
	ep2->type = nextc;
	escan();
	addekid(ep2, ep1);
	addekid(ep2, getE3());
	if (esupport&E_RCONST) {
		EPNODE	*ep3 = ep1->sibling;
		if ((ep1->type == NUM) & (ep3->type == NUM)) {
			ep2 = rconst(ep2);
		} else if (ep3->type == NUM) {
			if (ep2->type == '/') {
				if (ep3->v.num == 0)
					esyntax("divide by zero constant");
				ep2->type = '*';	/* for speed */
				ep3->v.num = 1./ep3->v.num;
			} else if (ep3->v.num == 0) {
				ep1->sibling = NULL;	/* (E2 * 0) */
				epfree(ep2,1);
				ep2 = ep3;
			}
		} else if (ep1->type == NUM && ep1->v.num == 0) {
			epfree(ep3,1);		/* (0 * E3) or (0 / E3) */
			ep1->sibling = NULL;
			efree(ep2);
			ep2 = ep1;
		}
	}
	ep1 = ep2;
    }
    return(ep1);
}


EPNODE *
getE3(void)			/* E3 -> E4 ^ E3 */
				/*	 E4 */
{
	EPNODE  *ep1, *ep2;

	ep1 = getE4();
	if (nextc != '^')
		return(ep1);
	ep2 = newnode();
	ep2->type = nextc;
	escan();
	addekid(ep2, ep1);
	addekid(ep2, getE3());
	if (esupport&E_RCONST) {
		EPNODE	*ep3 = ep1->sibling;
		if ((ep1->type == NUM) & (ep3->type == NUM)) {
			ep2 = rconst(ep2);
		} else if (ep1->type == NUM && ep1->v.num == 0) {
			epfree(ep3,1);		/* (0 ^ E3) */
			ep1->sibling = NULL;
			efree(ep2);
			ep2 = ep1;
		} else if ((ep3->type == NUM && ep3->v.num == 0) |
				(ep1->type == NUM && ep1->v.num == 1)) {
			epfree(ep2,0);		/* (E4 ^ 0) or (1 ^ E3) */
			ep2->type = NUM;
			ep2->v.num = 1;
		} else if (ep3->type == NUM && ep3->v.num == 1) {
			efree(ep3);	/* (E4 ^ 1) */
			ep1->sibling = NULL;
			efree(ep2);
			ep2 = ep1;
		}
	}
	return(ep2);
}


EPNODE *
getE4(void)			/* E4 -> ADDOP E5 */
				/*	 E5 */
{
    EPNODE  *ep1, *ep2;

    if (nextc == '-') {
	escan();
	ep2 = getE5();
	if (ep2->type == NUM) {
		ep2->v.num = -ep2->v.num;
		return(ep2);
	}
	if (ep2->type == UMINUS) {	/* don't generate -(-E5) */
	    ep1 = ep2->v.kid;
	    efree(ep2);
	    return(ep1);
	}
	ep1 = newnode();
	ep1->type = UMINUS;
	addekid(ep1, ep2);
	return(ep1);
    }
    if (nextc == '+')
	escan();
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
	EPNODE  *ep1, *ep2;

	if (nextc == '(') {
		escan();
		ep1 = getE1();
		if (nextc != ')')
			esyntax("')' expected");
		escan();
		return(ep1);
	}
	if (esupport&E_INCHAN && nextc == '$') {
		escan();
		ep1 = newnode();
		ep1->type = CHAN;
		ep1->v.chan = getinum();
		return(ep1);
	}
	if (esupport&(E_VARIABLE|E_FUNCTION) &&
			(isalpha(nextc) | (nextc == CNTXMARK))) {
		nam = getname();
		ep1 = NULL;
		if ((esupport&(E_VARIABLE|E_FUNCTION)) == (E_VARIABLE|E_FUNCTION)
				&& ecurfunc != NULL)
			for (i = 1, ep2 = ecurfunc->v.kid->sibling;
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
				escan();
				addekid(ep1, getE1());
			} while (nextc == ',');
			if (nextc != ')')
				esyntax("')' expected");
			escan();
		} else if (!(esupport&E_VARIABLE))
			esyntax("'(' expected");
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
	esyntax("unexpected character");
	return NULL; /* pro forma return */
}


EPNODE *
rconst(			/* reduce a constant expression */
    EPNODE	 *epar
)
{
    EPNODE  *ep;

    ep = newnode();
    ep->type = NUM;
    errno = 0;
    ep->v.num = evalue(epar);
    if ((errno == EDOM) | (errno == ERANGE))
	esyntax("bad constant expression");
    epfree(epar,1);
 
    return(ep);
}


int
isconstvar(			/* is ep linked to a constant expression? */
    EPNODE	 *ep
)
{
    EPNODE  *ep1;

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
    EPNODE	 *ep
)
{
    EPNODE  *dp;
    ELIBR  *lp;

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
