/* RCSid: $Id: calcomp.h,v 2.6 2003/02/22 02:07:21 greg Exp $ */
/*
 *  calcomp.h - header file for expression parser.
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */
				/* EPNODE types */
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

#define	 MAXWORD	127		/* maximum word/id length */
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

extern double  (*eoper[])();
extern unsigned long  eclock;
extern unsigned int  esupport;
extern EPNODE	*curfunc;
extern int  nextc;

#ifdef NOPROTO

extern void	fcompile();
extern void	scompile();
extern double	varvalue();
extern double	evariable();
extern void	varset();
extern void	dclear();
extern void	dremove();
extern int	vardefined();
extern char	*setcontext();
extern char	*pushcontext();
extern char	*popcontext();
extern char	*qualname();
extern int	incontext();
extern void	chanout();
extern void	dcleanup();
extern EPNODE	*dlookup();
extern VARDEF	*varlookup();
extern VARDEF	*varinsert();
extern void	varfree();
extern EPNODE	*dfirst();
extern EPNODE	*dnext();
extern EPNODE	*dpop();
extern void	dpush();
extern void	addchan();
extern void	getstatement();
extern EPNODE	*getdefn();
extern EPNODE	*getchan();
extern EPNODE	*eparse();
extern double	eval();
extern int	epcmp();
extern void	epfree();
extern EPNODE	*ekid();
extern int	nekids();
extern void	initfile();
extern void	initstr();
extern void	getscanpos();
extern int	scan();
extern char	*long2ascii();
extern void	syntax();
extern void	addekid();
extern char	*getname();
extern int	getinum();
extern double	getnum();
extern EPNODE	*getE1();
extern EPNODE	*getE2();
extern EPNODE	*getE3();
extern EPNODE	*getE4();
extern EPNODE	*getE5();
extern EPNODE	*rconst();
extern int	isconstvar();
extern int	isconstfun();
extern int	fundefined();
extern double	funvalue();
extern void	funset();
extern int	nargum();
extern double	argument();
extern VARDEF	*argf();
extern char	*argfun();
extern double	efunc();
extern LIBR	*liblookup();
extern void	libupdate();
extern void	eprint();
extern void	dprint();
extern char	*savestr();
extern void	freestr();
extern int	shash();
extern char	*emalloc();
extern char	*ecalloc();
extern char	*erealloc();
extern void	efree();
extern void	eputs();
extern void	wputs();
extern void	quit();

extern double	chanvalue();

#else
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
extern void	chanout(int (*cs)());
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
extern void	getstatement();
extern EPNODE	*getdefn();
extern EPNODE	*getchan();
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
					/* miscellaneous */
extern void	eputs(char *s);
extern void	wputs(char *s);
extern void	quit(int code);
					/* defined by client */
extern double	chanvalue(int n);

#endif
