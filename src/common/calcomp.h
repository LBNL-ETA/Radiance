/* RCSid $Id: calcomp.h,v 2.15 2003/08/04 19:20:26 greg Exp $ */
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
#define	 TICK		7
#define	 SYM		8
				/* also: '+', '-', '*', '/', '^', '=', ':' */

typedef struct {
    char  *fname;		/* function name */
    short  nargs;		/* # of required arguments */
    short  atyp;		/* assignment type (':' or '=') */
    double  (*f)();		/* pointer to function */
}  LIBR;		/* a library function */

typedef struct epnode {
    int	 type;			/* node type */
    struct epnode  *sibling;	/* next child this level */
    union {
	struct epnode  *kid;	/* first child */
	double	num;		/* number */
	char  *name;		/* symbol name */
	int  chan;		/* channel number */
	unsigned long  tick;	/* timestamp */
	struct vardef {
	    char  *name;		/* variable name */
	    int	 nlinks;		/* number of references */
	    struct epnode  *def;	/* definition */
	    LIBR  *lib;			/* library definition */
	    struct vardef  *next;	/* next in hash list */
	}  *ln;			/* link */
    } v;		/* value */
}  EPNODE;	/* an expression node */

typedef struct vardef  VARDEF;	/* a variable definition */

#define	 RMAXWORD	127		/* maximum word/id length */
#define	 CNTXMARK	'`'		/* context mark */

#define	 isid(c)	(isalnum(c) || (c) == '_' || \
			(c) == '.' || (c) == CNTXMARK)

#define	 evalue(ep)	(*eoper[(ep)->type])(ep)

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
extern EPNODE	*curfunc;
extern int  nextc;

					/* defined in biggerlib.c */
extern void biggerlib(void);

					/* defined in caldefn.c */
extern void	fcompile(char *fname);
extern void	scompile(char *str, char *fname, int ln);
extern double	varvalue(char *vname);
extern double	evariable(EPNODE *ep);
extern void	varset(char *fname, int assign, double val);
extern void	dclear(char *name);
extern void	dremove(char *name);
extern int	vardefined(char *name);
extern char	*setcontext(char *ctx);
extern char	*pushcontext(char *ctx);
extern char	*popcontext(void);
extern char	*qualname(char *nam, int lvl);
extern int	incontext(char *qn);
extern void	chanout(void (*cs)(int n, double v));
extern void	dcleanup(int lvl);
extern EPNODE	*dlookup(char *name);
extern VARDEF	*varlookup(char *name);
extern VARDEF	*varinsert(char *name);
extern void	varfree(VARDEF *ln);
extern EPNODE	*dfirst(void);
extern EPNODE	*dnext(void);
extern EPNODE	*dpop(char *name);
extern void	dpush(char *nm, EPNODE *ep);
extern void	addchan(EPNODE *sp);
extern void	getstatement(void);
extern EPNODE	*getdefn(void);
extern EPNODE	*getchan(void);
					/* defined in calexpr.c */
extern EPNODE	*eparse(char *expr);
extern double	eval(char *expr);
extern int	epcmp(EPNODE *ep1, EPNODE *ep2);
extern void	epfree(EPNODE *epar);
extern EPNODE	*ekid(EPNODE *ep, int n);
extern int	nekids(EPNODE *ep);
extern void	initfile(FILE *fp, char *fn, int ln);
extern void	initstr(char *s, char *fn, int ln);
extern void	getscanpos(char **fnp, int *lnp, char **spp, FILE **fpp);
extern int	scan(void);
extern char	*long2ascii(long l);
extern void	syntax(char *err);
extern void	addekid(EPNODE *ep, EPNODE *ekid);
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
extern void	funset(char *fname, int nargs, int assign, double (*fptr)());
extern int	nargum(void);
extern double	argument(int n);
extern VARDEF	*argf(int n);
extern char	*argfun(int n);
extern double	efunc(EPNODE *ep);
extern LIBR	*liblookup(char *fname);
extern void	libupdate(char *fn);
					/* defined in calprnt.c */
extern void	eprint(EPNODE *ep, FILE *fp);
extern void	dprint(char *name, FILE *fp);
					/* defined in savestr.c */
extern char	*savestr(char *str);
extern void	freestr(char *s);
extern int	shash(char *s);
					/* defined in ealloc.c */
extern char	*emalloc(unsigned int n);
extern char	*ecalloc(unsigned int ne, unsigned int es);
extern char	*erealloc(char *cp, unsigned int n);
extern void	efree(char *cp);
					/* defined by client */
extern double	chanvalue(int n);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_CALCOMP_H_ */

