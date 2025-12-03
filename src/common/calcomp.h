/* RCSid $Id$ */
/*
 *  calcomp.h - header file for expression parser.
 */
#ifndef _RAD_CALCOMP_H_
#define _RAD_CALCOMP_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define	 VAR		1
#define	 NUM		2
#define	 UMINUS		3
#define	 CHAN		4
#define	 FUNC		5
#define	 ARG		6
#define	 CLKT		7
#define	 SYM		8
				/* also: '+', '-', '*', '/', '^', '=', ':' */

typedef struct {
    char  *fname;		/* function name */
    short  nargs;		/* # of required arguments */
    short  atyp;		/* assignment type (':' or '=') */
    double  (*f)(char *);	/* pointer to function */
}  ELIBR;		/* a library function */

typedef struct vardef {
	    char  *name;		/* variable name */
	    int	 nlinks;		/* number of references */
	    struct epnode  *def;	/* definition */
	    ELIBR  *lib;		/* library definition */
	    struct vardef  *next;	/* next in hash list */
}  VARDEF;		/* a variable definition */

typedef struct epnode {
    union {
	struct epnode  *kid;	/* first child */
	double	num;		/* number */
	char  *name;		/* symbol name */
	int  chan;		/* channel number */
	unsigned long  tick;	/* timestamp */
	VARDEF  *ln;		/* variable definition link */
    } v;		/* value */
    struct epnode  *sibling;	/* next child this level */
    short  type;		/* node type */
    short  nkids;		/* child count (neg if array) */
}  EPNODE;	/* an expression node */

#define  nekids(ep)	abs((ep)->nkids)

#define	 RMAXWORD	127		/* maximum word/id length */
#define	 CNTXMARK	'`'		/* context mark */

#define	 isid(c)	(isalnum(c) || (c) == '_' || \
			(c) == '.' || (c) == CNTXMARK)

#define	 evalue(ep)	(*eoper[(ep)->type])(ep)

#define	 dfn_name(ep)	((ep)->v.kid->type == SYM ? \
			(ep)->v.kid->v.name : \
			(ep)->v.kid->v.kid->v.name)

					/* flags to set in esupport */
#define  E_VARIABLE	001
#define  E_FUNCTION	002
#define  E_INCHAN	004
#define  E_OUTCHAN	010
#define  E_RCONST	020
#define  E_REDEFW	040

extern double  (*eoper[])(EPNODE *);
extern unsigned long  eclock;
extern unsigned int  esupport;
extern EPNODE	*ecurfunc;
extern int  nextc;
extern int  eofc;
					/* defined in biggerlib.c */
extern void biggerlib(void);
					/* defined in caldefn.c */
extern void	fcompile(char *fname);
extern void	scompile(char *str, char *fname, int ln);
extern double	varvalue(char *vname);
extern double	evariable(EPNODE *ep);
extern void	varset(char *vname, int assign, double val);
extern void	dclear(char *name);
extern void	dremove(char *name);
extern int	vardefined(char *name);
extern char	*calcontext(char *ctx);
extern char	*pushcontext(char *ctx);
extern char	*popcontext(void);
extern char	*qualname(char *nam, int lvl);
extern int	incontext(char *qn);
extern void	chanout(void (*cs)(int n, double v));
extern void	doptimize(int activate);
extern void	dcleanup(int lvl);
extern EPNODE	*dlookup(char *name);
extern VARDEF	*varlookup(char *name);
extern VARDEF	*varinsert(char *name);
extern void	varfree(VARDEF *ln);
extern EPNODE	*dfirst(void);
extern EPNODE	*dnext(void);
extern EPNODE	*dpop(char *name);
extern void	dpush(char *nm, EPNODE *ep);
extern void	eaddchan(EPNODE *sp);
extern void	egetstatement(void);
extern EPNODE	*egetdefn(void);
extern EPNODE	*egetchan(void);
					/* defined in calexpr.c */
extern EPNODE	*eparse(char *expr);
extern double	eval(char *expr);
extern int	epcmp(EPNODE *ep1, EPNODE *ep2);
extern void	epfree(EPNODE *epar, int frep);
extern void	epoptimize(EPNODE *epar);
extern EPNODE	*ekid(EPNODE *ep, int n);
extern void	initfile(FILE *fp, char *fn, int ln);
extern void	initstr(char *s, char *fn, int ln);
extern void	getscanpos(char **fnp, int *lnp, char **spp, FILE **fpp);
extern int	escan(void);
extern char	*long2ascii(long l);
extern void	esyntax(char *err);
extern void	addekid(EPNODE *ep, EPNODE *ek);
extern char	*getname(void);
extern int	getinum(void);
extern double	getnum(void);
extern EPNODE	*getE1(void);
extern EPNODE	*getE2(void);
extern EPNODE	*getE3(void);
extern EPNODE	*getE4(void);
extern EPNODE	*getE5(void);
extern EPNODE	*rconst(EPNODE *epar);
extern int	isconstvar(EPNODE *ep);
extern int	isconstfun(EPNODE *ep);
					/* defined in calfunc.c */
extern int	fundefined(char *fname);
extern double	funvalue(char *fname, int n, double *a);
extern void	funset(char *fname, int nargs, int assign,
				double (*fptr)(char *));
extern int	nargum(void);
extern double	argument(int n);
extern VARDEF	*eargf(int n);
extern char	*eargfun(int n);
extern double	efunc(EPNODE *ep);
extern double	eargument(EPNODE *ep);
extern ELIBR	*eliblookup(char *fname);
extern void	elibupdate(char *fn);
					/* defined in calprnt.c */
extern void	eprint(EPNODE *ep, FILE *fp);
extern void	dprint(char *name, FILE *fp);
					/* defined by client */
extern double	chanvalue(int n);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_CALCOMP_H_ */

