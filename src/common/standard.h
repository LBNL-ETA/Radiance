/* RCSid $Id: standard.h,v 2.22 2003/05/13 17:58:32 greg Exp $ */
/*
 *	Miscellaneous definitions required by many routines.
 */

#include "copyright.h"

#include  <stdio.h>

#include  <sys/types.h>

#include  <fcntl.h>

#include  <math.h>

#include  <errno.h>

#include  <stdlib.h>

#include  <string.h>

#include  "mat4.h"
				/* regular transformation */
typedef struct {
	MAT4  xfm;				/* transform matrix */
	FLOAT  sca;				/* scalefactor */
}  XF;
				/* complemetary tranformation */
typedef struct {
	XF  f;					/* forward */
	XF  b;					/* backward */
}  FULLXF;

#ifndef  PI
#ifdef	M_PI
#define	 PI		((double)M_PI)
#else
#define	 PI		3.14159265358979323846
#endif
#endif

#ifndef	 F_OK			/* mode bits for access(2) call */
#define	 R_OK		4		/* readable */
#define	 W_OK		2		/* writable */
#define	 X_OK		1		/* executable */
#define	 F_OK		0		/* exists */
#endif

#ifndef  int2
#define  int2		short		/* two-byte integer */
#endif
#ifndef  int4
#define  int4		int		/* four-byte integer */
#endif

				/* error codes */
#define	 WARNING	0		/* non-fatal error */
#define	 USER		1		/* fatal user-caused error */
#define	 SYSTEM		2		/* fatal system-related error */
#define	 INTERNAL	3		/* fatal program-related error */
#define	 CONSISTENCY	4		/* bad consistency check, abort */
#define	 COMMAND	5		/* interactive error */
#define  NERRS		6
				/* error struct */
extern struct erract {
	char	pre[16];		/* prefix message */
	void	(*pf)();		/* put function (resettable) */
	int	ec;			/* exit code (0 means non-fatal) */
} erract[NERRS];	/* list of error actions */

#define  ERRACT_INIT	{	{"warning - ", wputs, 0}, \
				{"fatal - ", eputs, 1}, \
				{"system - ", eputs, 2}, \
				{"internal - ", eputs, 3}, \
				{"consistency - ", eputs, -1}, \
				{"", NULL, 0}	}

extern char  errmsg[];			/* global buffer for error messages */

#ifdef  FASTMATH
#define  tcos			cos
#define  tsin			sin
#define  ttan			tan
#else
extern double	tcos();			/* table-based cosine approximation */
#define  tsin(x)		tcos((x)-(PI/2.))
#define  ttan(x)		(tsin(x)/tcos(x))
#endif
					/* custom version of assert(3) */
#define  CHECK(be,et,em)	if (be) error(et,em); else
#ifdef  DEBUG
#define  DCHECK			CHECK
#else
#define  DCHECK(be,et,em)	(void)0
#endif
					/* memory operations */
#ifdef	NOSTRUCTASS
#define	 copystruct(d,s)	bcopy((void *)(s),(void *)(d),sizeof(*(d)))
#else
#define	 copystruct(d,s)	(*(d) = *(s))
#endif

#ifndef BSD
#define	 bcopy(s,d,n)		(void)memcpy(d,s,n)
#define	 bzero(d,n)		(void)memset(d,0,n)
#define	 bcmp(b1,b2,n)		memcmp(b1,b2,n)
#define	 index			strchr
#define	 rindex			strrchr
#endif
extern off_t  lseek();

#ifdef MSDOS
#define NIX 1
#endif
#ifdef AMIGA
#define NIX 1
#endif

#ifdef NOPROTO

extern int	badarg();
extern char	*bmalloc();
extern void	bfree();
extern void	error();
extern int	expandarg();
extern time_t	fdate();
extern int	setfdate();
extern char	*fgetline();
extern int	fgetval();
extern char	*fgetword();
extern void	fputword();
extern char	*fixargv0();
extern FILE	*frlibopen();
extern char	*getrlibpath();
extern char	*getpath();
extern void	putstr();
extern void	putint();
extern void	putflt();
extern char	*getstr();
extern long	getint();
extern double	getflt();
extern int	open_process();
extern int	process();
extern int	close_process();
extern int	readbuf();
extern int	writebuf();
extern int	ecompile();
extern char	*expsave();
extern void	expset();
extern char	*eindex();
extern char	*savestr();
extern void	freestr();
extern int	shash();
extern char	*savqstr();
extern void	freeqstr();
extern double	tcos();
extern int	wordfile();
extern int	wordstring();
extern char	*atos();
extern char	*nextword();
extern char	*sskip();
extern char	*sskip2();
extern char	*iskip();
extern char	*fskip();
extern int	isint();
extern int	isintd();
extern int	isflt();
extern int	isfltd();
extern int	xf();
extern int	invxf();
extern int	fullxf();
extern int	quadtratic();
extern int4	encodedir();
extern void	decodedir();
extern double	dir2diff();
extern double	fdir2diff();
extern void	eputs();
extern void	wputs();
extern void	quit();

#else
					/* defined in badarg.c */
extern int	badarg(int ac, char **av, char *fl);
					/* defined in bmalloc.c */
extern char	*bmalloc(unsigned int n);
extern void	bfree(char *p, unsigned int n);
					/* defined in error.c */
extern void	error(int etype, char *emsg);
					/* defined in expandarg.c */
extern int	expandarg(int *acp, char ***avp, int n);
					/* defined in fdate.c */
extern time_t	fdate(char *fname);
extern int	setfdate(char *fname, long ftim);
					/* defined in fgetline.c */
extern char	*fgetline(char *s, int n, FILE *fp);
					/* defined in fgetval.c */
extern int	fgetval(FILE *fp, int ty, char *vp);
					/* defined in fgetword.c */
extern char	*fgetword(char *s, int n, FILE *fp);
					/* defined in fputword.c */
extern void	fputword(char *s, FILE *fp);
					/* defined in fixargv0.c */
extern char	*fixargv0(char *av0);
					/* defined in fropen.c */
extern FILE	*frlibopen(char *fname);
					/* defined in getlibpath.c */
extern char	*getrlibpath(void);
					/* defined in getpath.c */
extern char	*getpath(char *fname, char *searchpath, int mode);
					/* defined in portio.c */
extern void	putstr(char *s, FILE *fp);
extern void	putint(long i, int siz, FILE *fp);
extern void	putflt(double f, FILE *fp);
extern char	*getstr(char *s, FILE *fp);
extern long	getint(int siz, FILE *fp);
extern double	getflt(FILE *fp);
					/* defined in process.c */
extern int	open_process(int pd[3], char *av[]);
extern int	process(int pd[3], char *recvbuf, char *sendbuf,
				int nbr, int nbs);
extern int	close_process(int pd[3]);
extern int	readbuf(int fd, char *bpos, int siz);
extern int	writebuf(int fd, char *bpos, int siz);
					/* defined in rexpr.c */
extern int	ecompile(char *sp, int iflg, int wflag);
extern char	*expsave(void);
extern void	expset(char *ep);
extern char	*eindex(char *sp);
					/* defined in savestr.c */
extern char	*savestr(char *str);
extern void	freestr(char *s);
extern int	shash(char *s);
					/* defined in savqstr.c */
extern char	*savqstr(char *s);
extern void	freeqstr(char *s);
					/* defined in tcos.c */
extern double	tcos(double x);
					/* defined in wordfile.c */
extern int	wordfile(char **words, char *fname);
extern int	wordstring(char **avl, char *str);
					/* defined in words.c */
extern char	*atos(char *rs, int nb, char *s);
extern char	*nextword(char *cp, int nb, char *s);
extern char	*sskip(char *s);
extern char	*sskip2(char *s, int n);
extern char	*iskip(char *s);
extern char	*fskip(char *s);
extern int	isint(char *s);
extern int	isintd(char *s, char *ds);
extern int	isflt(char *s);
extern int	isfltd(char *s, char *ds);
					/* defined in xf.c */
extern int	xf(XF *ret, int ac, char *av[]);
extern int	invxf(XF *ret, int ac, char *av[]);
extern int	fullxf(FULLXF *fx, int ac, char *av[]);
					/* defined in zeroes.c */
extern int	quadtratic(double *r, double a, double b, double c);
					/* defined in dircode.c */
extern int4	encodedir(FVECT dv);
extern void	decodedir(FVECT dv, int4 dc);
extern double	dir2diff(int4 dc1, int4 dc2);
extern double	fdir2diff(int4 dc1, FVECT v2);
					/* miscellaneous */
extern void	eputs(char *s);
extern void	wputs(char *s);
extern void	quit(int code);

#endif
