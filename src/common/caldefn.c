#ifndef lint
static const char	RCSid[] = "$Id: caldefn.c,v 2.23 2008/05/01 16:42:05 greg Exp $";
#endif
/*
 *  Store variable definitions.
 *
 *  7/1/85  Greg Ward
 *
 *  11/11/85  Added conditional compiles (OUTCHAN) for control output.
 *
 *  4/2/86  Added conditional compiles for function definitions (FUNCTION).
 *
 *  1/15/88  Added clock for caching of variable values.
 *
 *  11/16/88  Added VARDEF structure for hard linking.
 *
 *  5/31/90  Added conditional compile (REDEFW) for redefinition warning.
 *
 *  4/23/91  Added ':' assignment for constant expressions
 *
 *  8/7/91  Added optional context path to append to variable names
 *
 *  5/17/2001  Fixed clock counter wrapping behavior
 *
 *  2/19/03	Eliminated conditional compiles in favor of esupport extern.
 */

#include "copyright.h"

#include  <stdio.h>
#include  <string.h>
#include  <ctype.h>

#include  "rterror.h"
#include  "rtio.h"
#include  "rtmisc.h"
#include  "calcomp.h"

#ifndef	 NHASH
#define	 NHASH		521		/* hash size (a prime!) */
#endif

#define  hash(s)	(shash(s)%NHASH)

#define	 newnode()	(EPNODE *)ecalloc(1, sizeof(EPNODE))

static double  dvalue(char  *name, EPNODE *d);

#define  MAXCLOCK	(1L<<31)	/* clock wrap value */

unsigned long  eclock = 0;		/* value storage timer */

#define  MAXCNTX	1023		/* maximum context length */

static char  context[MAXCNTX+1];	/* current context path */

static VARDEF  *hashtbl[NHASH];		/* definition list */
static int  htndx;			/* index for */		
static VARDEF  *htpos;			/* ...dfirst() and */
static EPNODE  *ochpos;			/* ...dnext */
static EPNODE  *outchan;

EPNODE	*curfunc = NULL;
#define	 dname(ep)	((ep)->v.kid->type == SYM ? \
			(ep)->v.kid->v.name : \
			(ep)->v.kid->v.kid->v.name)


void
fcompile(			/* get definitions from a file */
	char  *fname
)
{
    FILE  *fp;

    if (fname == NULL)
	fp = stdin;
    else if ((fp = fopen(fname, "r")) == NULL) {
	eputs(fname);
	eputs(": cannot open\n");
	quit(1);
    }
    initfile(fp, fname, 0);
    while (nextc != EOF)
	getstatement();
    if (fname != NULL)
	fclose(fp);
}


void
scompile(		/* get definitions from a string */
	char  *str,
	char  *fn,
	int  ln
)
{
    initstr(str, fn, ln);
    while (nextc != EOF)
	getstatement();
}


double
varvalue(			/* return a variable's value */
	char  *vname
)
{
    return(dvalue(vname, dlookup(vname)));
}


double
evariable(			/* evaluate a variable */
	EPNODE	*ep
)
{
    register VARDEF  *dp = ep->v.ln;

    return(dvalue(dp->name, dp->def));
}


void
varset(		/* set a variable's value */
	char  *vname,
	int  assign,
	double	val
)
{
    char  *qname;
    register EPNODE  *ep1, *ep2;
					/* get qualified name */
    qname = qualname(vname, 0);
					/* check for quick set */
    if ((ep1 = dlookup(qname)) != NULL && ep1->v.kid->type == SYM) {
	ep2 = ep1->v.kid->sibling;
	if (ep2->type == NUM) {
	    ep2->v.num = val;
	    ep1->type = assign;
	    return;
	}
    }
					/* hand build definition */
    ep1 = newnode();
    ep1->type = assign;
    ep2 = newnode();
    ep2->type = SYM;
    ep2->v.name = savestr(vname);
    addekid(ep1, ep2);
    ep2 = newnode();
    ep2->type = NUM;
    ep2->v.num = val;
    addekid(ep1, ep2);
    dremove(qname);
    dpush(qname, ep1);
}


void
dclear(			/* delete variable definitions of name */
	char  *name
)
{
    register EPNODE  *ep;

    while ((ep = dpop(name)) != NULL) {
	if (ep->type == ':') {
	    dpush(name, ep);		/* don't clear constants */
	    return;
	}
	epfree(ep);
    }
}


void
dremove(			/* delete all definitions of name */
	char  *name
)
{
    register EPNODE  *ep;

    while ((ep = dpop(name)) != NULL)
	epfree(ep);
}


int
vardefined(	/* return non-zero if variable defined */
	char  *name
)
{
    register EPNODE  *dp;

    return((dp = dlookup(name)) != NULL && dp->v.kid->type == SYM);
}


char *
setcontext(			/* set a new context path */
	register char  *ctx
)
{
    register char  *cpp;

    if (ctx == NULL)
	return(context);		/* just asking */
    while (*ctx == CNTXMARK)
	ctx++;				/* skip past marks */
    if (!*ctx) {
	context[0] = '\0';		/* empty means clear context */
	return(context);
    }
    cpp = context;			/* start context with mark */
    *cpp++ = CNTXMARK;
    do {				/* carefully copy new context */
	if (cpp >= context+MAXCNTX)
	    break;			/* just copy what we can */
	if (isid(*ctx))
	    *cpp++ = *ctx++;
	else {
	    *cpp++ = '_'; ctx++;
	}
    } while (*ctx);
    while (cpp[-1] == CNTXMARK)		/* cannot end in context mark */
	cpp--;
    *cpp = '\0';
    return(context);
}


char *
pushcontext(		/* push on another context */
	char  *ctx
)
{
    char  oldcontext[MAXCNTX+1];
    register int  n;

    strcpy(oldcontext, context);	/* save old context */
    setcontext(ctx);			/* set new context */
    n = strlen(context);		/* tack on old */
    if (n+strlen(oldcontext) > MAXCNTX) {
	strncpy(context+n, oldcontext, MAXCNTX-n);
	context[MAXCNTX] = '\0';
    } else
	strcpy(context+n, oldcontext);
    return(context);
}


char *
popcontext(void)			/* pop off top context */
{
    register char  *cp1, *cp2;

    if (!context[0])			/* nothing left to pop */
	return(context);
    cp2 = context;			/* find mark */
    while (*++cp2 && *cp2 != CNTXMARK)
	;
    cp1 = context;			/* copy tail to front */
    while ( (*cp1++ = *cp2++) )
	;
    return(context);
}


char *
qualname(		/* get qualified name */
	register char  *nam,
	int  lvl
)
{
    static char	 nambuf[RMAXWORD+1];
    register char  *cp = nambuf, *cpp;
				/* check for explicit local */
    if (*nam == CNTXMARK)
	if (lvl > 0)		/* only action is to refuse search */
	    return(NULL);
	else
	    nam++;
    else if (nam == nambuf)	/* check for repeat call */
	return(lvl > 0 ? NULL : nam);
				/* copy name to static buffer */
    while (*nam) {
	if (cp >= nambuf+RMAXWORD)
		goto toolong;
	*cp++ = *nam++;
    }
				/* check for explicit global */
    if (cp > nambuf && cp[-1] == CNTXMARK) {
	if (lvl > 0)
	    return(NULL);
	*--cp = '\0';
	return(nambuf);		/* already qualified */
    }
    cpp = context;		/* else skip the requested levels */
    while (lvl-- > 0) {
	if (!*cpp)
	    return(NULL);	/* return NULL if past global level */
	while (*++cpp && *cpp != CNTXMARK)
	    ;
    }
    while (*cpp) {		/* copy context to static buffer */
	if (cp >= nambuf+RMAXWORD)
	    goto toolong;
	*cp++ = *cpp++;
    }
toolong:
    *cp = '\0';
    return(nambuf);		/* return qualified name */
}


int
incontext(			/* is qualified name in current context? */
	register char  *qn
)
{
    if (!context[0])			/* global context accepts all */
	return(1);
    while (*qn && *qn != CNTXMARK)	/* find context mark */
	qn++;
    return(!strcmp(qn, context));
}


void
chanout(			/* set output channels */
	void  (*cs)(int n, double v)
)
{
    register EPNODE  *ep;

    for (ep = outchan; ep != NULL; ep = ep->sibling)
	(*cs)(ep->v.kid->v.chan, evalue(ep->v.kid->sibling));

}


void
dcleanup(		/* clear definitions (0->vars,1->output,2->consts) */
	int  lvl
)
{
    register int  i;
    register VARDEF  *vp;
    register EPNODE  *ep;
				/* if context is global, clear all */
    for (i = 0; i < NHASH; i++)
	for (vp = hashtbl[i]; vp != NULL; vp = vp->next)
	    if (incontext(vp->name)) {
		if (lvl >= 2)
		    dremove(vp->name);
		else
		    dclear(vp->name);
	    }
    if (lvl >= 1) {
	for (ep = outchan; ep != NULL; ep = ep->sibling)
	    epfree(ep);
	outchan = NULL;
    }
}


EPNODE *
dlookup(			/* look up a definition */
	char  *name
)
{
    register VARDEF  *vp;
    
    if ((vp = varlookup(name)) == NULL)
	return(NULL);
    return(vp->def);
}


VARDEF *
varlookup(			/* look up a variable */
	char  *name
)
{
    int	 lvl = 0;
    register char  *qname;
    register VARDEF  *vp;
				/* find most qualified match */
    while ((qname = qualname(name, lvl++)) != NULL)
	for (vp = hashtbl[hash(qname)]; vp != NULL; vp = vp->next)
	    if (!strcmp(vp->name, qname))
		return(vp);
    return(NULL);
}


VARDEF *
varinsert(			/* get a link to a variable */
	char  *name
)
{
    register VARDEF  *vp;
    int	 hv;
    
    if ((vp = varlookup(name)) != NULL) {
	vp->nlinks++;
	return(vp);
    }
    vp = (VARDEF *)emalloc(sizeof(VARDEF));
    vp->lib = liblookup(name);
    if (vp->lib == NULL)		/* if name not in library */
	name = qualname(name, 0);	/* use fully qualified version */
    hv = hash(name);
    vp->name = savestr(name);
    vp->nlinks = 1;
    vp->def = NULL;
    vp->next = hashtbl[hv];
    hashtbl[hv] = vp;
    return(vp);
}


void
libupdate(			/* update library links */
	char  *fn
)
{
    register int  i;
    register VARDEF  *vp;
					/* if fn is NULL then relink all */
    for (i = 0; i < NHASH; i++)
	for (vp = hashtbl[i]; vp != NULL; vp = vp->next)
	    if (vp->lib != NULL || fn == NULL || !strcmp(fn, vp->name))
		vp->lib = liblookup(vp->name);
}


void
varfree(				/* release link to variable */
	register VARDEF	 *ln
)
{
    register VARDEF  *vp;
    int	 hv;

    if (--ln->nlinks > 0)
	return;				/* still active */

    hv = hash(ln->name);
    vp = hashtbl[hv];
    if (vp == ln)
	hashtbl[hv] = vp->next;
    else {
	while (vp->next != ln)		/* must be in list */
		vp = vp->next;
	vp->next = ln->next;
    }
    freestr(ln->name);
    efree((char *)ln);
}


EPNODE *
dfirst(void)			/* return pointer to first definition */
{
    htndx = 0;
    htpos = NULL;
    ochpos = outchan;
    return(dnext());
}


EPNODE *
dnext(void)				/* return pointer to next definition */
{
    register EPNODE  *ep;
    register char  *nm;

    while (htndx < NHASH) {
	if (htpos == NULL)
		htpos = hashtbl[htndx++];
	while (htpos != NULL) {
	    ep = htpos->def;
	    nm = htpos->name;
	    htpos = htpos->next;
	    if (ep != NULL && incontext(nm))
		return(ep);
	}
    }
    if ((ep = ochpos) != NULL)
	ochpos = ep->sibling;
    return(ep);
}


EPNODE *
dpop(			/* pop a definition */
	char  *name
)
{
    register VARDEF  *vp;
    register EPNODE  *dp;
    
    if ((vp = varlookup(name)) == NULL || vp->def == NULL)
	return(NULL);
    dp = vp->def;
    vp->def = dp->sibling;
    varfree(vp);
    return(dp);
}


void
dpush(			/* push on a definition */
	char  *nm,
	register EPNODE	 *ep
)
{
    register VARDEF  *vp;

    vp = varinsert(nm);
    ep->sibling = vp->def;
    vp->def = ep;
}


void
addchan(			/* add an output channel assignment */
	EPNODE	*sp
)
{
    int	 ch = sp->v.kid->v.chan;
    register EPNODE  *ep, *epl;

    for (epl = NULL, ep = outchan; ep != NULL; epl = ep, ep = ep->sibling)
	if (ep->v.kid->v.chan >= ch) {
	    if (epl != NULL)
		epl->sibling = sp;
	    else
		outchan = sp;
	    if (ep->v.kid->v.chan > ch)
		sp->sibling = ep;
	    else {
		sp->sibling = ep->sibling;
		epfree(ep);
	    }
	    return;
	}
    if (epl != NULL)
	epl->sibling = sp;
    else
	outchan = sp;
    sp->sibling = NULL;

}


void
getstatement(void)			/* get next statement */
{
    register EPNODE  *ep;
    char  *qname;
    register VARDEF  *vdef;

    if (nextc == ';') {		/* empty statement */
	scan();
	return;
    }
    if (esupport&E_OUTCHAN &&
		nextc == '$') {		/* channel assignment */
	ep = getchan();
	addchan(ep);
    } else {				/* ordinary definition */
	ep = getdefn();
	qname = qualname(dname(ep), 0);
	if (esupport&E_REDEFW && (vdef = varlookup(qname)) != NULL) {
	    if (vdef->def != NULL && epcmp(ep, vdef->def)) {
		wputs(qname);
		if (vdef->def->type == ':')
		    wputs(": redefined constant expression\n");
		else
		    wputs(": redefined\n");
	    } else if (ep->v.kid->type == FUNC && vdef->lib != NULL) {
		wputs(qname);
		wputs(": definition hides library function\n");
	    }
	}
	if (ep->type == ':')
	    dremove(qname);
	else
	    dclear(qname);
	dpush(qname, ep);
    }
    if (nextc != EOF) {
	if (nextc != ';')
	    syntax("';' expected");
	scan();
    }
}


EPNODE *
getdefn(void)
	/* A -> SYM = E1 */
	/*	SYM : E1 */
	/*	FUNC(SYM,..) = E1 */
	/*	FUNC(SYM,..) : E1 */
{
    register EPNODE  *ep1, *ep2;

    if (!isalpha(nextc) && nextc != CNTXMARK)
	syntax("illegal variable name");

    ep1 = newnode();
    ep1->type = SYM;
    ep1->v.name = savestr(getname());

    if (esupport&E_FUNCTION && nextc == '(') {
	ep2 = newnode();
	ep2->type = FUNC;
	addekid(ep2, ep1);
	ep1 = ep2;
	do {
	    scan();
	    if (!isalpha(nextc))
		syntax("illegal variable name");
	    ep2 = newnode();
	    ep2->type = SYM;
	    ep2->v.name = savestr(getname());
	    addekid(ep1, ep2);
	} while (nextc == ',');
	if (nextc != ')')
	    syntax("')' expected");
	scan();
	curfunc = ep1;
    }

    if (nextc != '=' && nextc != ':')
	syntax("'=' or ':' expected");

    ep2 = newnode();
    ep2->type = nextc;
    scan();
    addekid(ep2, ep1);
    addekid(ep2, getE1());

    if (ep1->type == SYM && ep1->sibling->type != NUM) {
	ep1 = newnode();
	ep1->type = CLKT;
	ep1->v.tick = 0;
	addekid(ep2, ep1);
	ep1 = newnode();
	ep1->type = NUM;
	addekid(ep2, ep1);
    }
    curfunc = NULL;

    return(ep2);
}


EPNODE *
getchan(void)			/* A -> $N = E1 */
{
    register EPNODE  *ep1, *ep2;

    if (nextc != '$')
	syntax("missing '$'");
    scan();

    ep1 = newnode();
    ep1->type = CHAN;
    ep1->v.chan = getinum();

    if (nextc != '=')
	syntax("'=' expected");
    scan();

    ep2 = newnode();
    ep2->type = '=';
    addekid(ep2, ep1);
    addekid(ep2, getE1());

    return(ep2);
}



/*
 *  The following routines are for internal use only:
 */


static double			/* evaluate a variable */
dvalue(char  *name, EPNODE	*d)
{
    register EPNODE  *ep1, *ep2;
    
    if (d == NULL || d->v.kid->type != SYM) {
	eputs(name);
	eputs(": undefined variable\n");
	quit(1);
    }
    ep1 = d->v.kid->sibling;			/* get expression */
    if (ep1->type == NUM)
	return(ep1->v.num);			/* return if number */
    ep2 = ep1->sibling;				/* check time */
    if (eclock >= MAXCLOCK)
	eclock = 1;				/* wrap clock counter */
    if (ep2->v.tick < MAXCLOCK &&
		(ep2->v.tick == 0) | (ep2->v.tick != eclock)) {
	ep2->v.tick = d->type == ':' ? MAXCLOCK : eclock;
	ep2 = ep2->sibling;
	ep2->v.num = evalue(ep1);		/* needs new value */
    } else
	ep2 = ep2->sibling;			/* else reuse old value */

    return(ep2->v.num);
}
