/* RCSid $Id: rtio.h,v 3.25 2021/02/19 16:15:23 greg Exp $ */
/*
 *	Radiance i/o and string routines
 */

#ifndef _RAD_RTIO_H_
#define _RAD_RTIO_H_

#include  <stdio.h>
#include  <sys/types.h>
#include  <fcntl.h>
#include  <string.h>
#include  <time.h>

#ifdef getc_unlocked		/* avoid horrendous overhead of flockfile */
#undef getc
#undef getchar
#undef putc
#undef putchar
#undef feof
#undef ferror
#define getc    getc_unlocked
#define getchar	getchar_unlocked
#define putc    putc_unlocked
#define putchar	putchar_unlocked
#ifndef __cplusplus
#define feof	feof_unlocked
#define ferror	ferror_unlocked
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif
					/* identify header lines */
#define  MAXFMTLEN	64
#define  isheadid(s)	headidval(NULL,s)
#define  isformat(s)	formatval(NULL,s)
#define  isdate(s)	dateval(NULL,s)
#define  isgmt(s)	gmtval(NULL,s)

#define  LATLONSTR	"LATLONG="
#define  LLATLONSTR	8
#define  islatlon(hl)		(!strncmp(hl,LATLONSTR,LLATLONSTR))
#define  latlonval(ll,hl)	sscanf((hl)+LLATLONSTR, "%f %f", \
						&(ll)[0],&(ll)[1])
#define  fputlatlon(lat,lon,fp)	fprintf(fp,"%s %.6f %.6f\n",LATLONSTR,lat,lon)
					/* defined in header.c */
extern void	newheader(const char *t, FILE *fp);
extern int	headidval(char *r, const char *s);
extern int	dateval(time_t *t, const char *s);
extern int	gmtval(time_t *t, const char *s);
extern void	fputdate(time_t t, FILE *fp);
extern void	fputnow(FILE *fp);
extern void	printargs(int ac, char **av, FILE *fp);
extern int	formatval(char fmt[MAXFMTLEN], const char *s);
extern void	fputformat(const char *s, FILE *fp);
extern int	nativebigendian(void);
extern int	isbigendian(const char *s);
extern void	fputendian(FILE *fp);
typedef int gethfunc(char *s, void *p); /* callback to process header lines */
extern int	getheader(FILE *fp, gethfunc *f, void *p);
extern int	globmatch(const char *pat, const char *str);
extern int	checkheader(FILE *fin, char fmt[MAXFMTLEN], FILE *fout);
					/* defined in fltdepth.c */
extern int	open_float_depth(const char *fname, long expected_length);
					/* defined in badarg.c */
extern int	badarg(int ac, char **av, char *fl);
					/* defined in expandarg.c */
extern int	envexpchr, filexpchr;
extern int	expandarg(int *acp, char ***avp, int n);
					/* defined in fdate.c */
extern time_t	fdate(char *fname);
extern int	setfdate(char *fname, long ftim);
					/* defined in fgetline.c */
extern char	*fgetline(char *s, int n, FILE *fp);
					/* defined in fgetval.c */
extern int	fgetval(FILE *fp, int ty, void *vp);
					/* defined in fgetword.c */
extern char	*fgetword(char *s, int n, FILE *fp);
					/* defined in fputword.c */
extern void	fputword(char *s, FILE *fp);
					/* defined in fropen.c */
extern FILE	*frlibopen(char *fname);
					/* defined in getlibpath.c */
extern char	*getrlibpath(void);
					/* defined in gethomedir.c */
extern char	*gethomedir(char *uname, char *path, int plen);
					/* defined in getpath.c */
extern char	*getpath(char *fname, char *searchpath, int mode);
					/* defined in byteswap.c */
extern void	swap16(char *wp, size_t n);
extern void	swap32(char *wp, size_t n);
extern void	swap64(char *wp, size_t n);
					/* defined in portio.c */
extern int	putstr(char *s, FILE *fp);
extern int	putint(long i, int siz, FILE *fp);
extern int	putflt(double f, FILE *fp);
extern size_t	putbinary(const void *s, size_t elsiz, size_t nel, FILE *fp);
extern char	*getstr(char *s, FILE *fp);
extern long	getint(int siz, FILE *fp);
extern double	getflt(FILE *fp);
extern size_t	getbinary(void *s, size_t elsiz, size_t nel, FILE *fp);
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
					/* defined in wordfile.c */
extern int	wordfile(char **words, int nargs, char *fname);
extern int	wordstring(char **avl, int nargs, char *str);
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
					/* defined in lamp.c */
extern float *	matchlamp(char *s);
extern int	loadlamps(char *file);
extern void	freelamps(void);

#ifndef strlcpy				/* defined in option strlcpy.c */
extern size_t	strlcpy(char *dst, const char *src, size_t siz);
extern size_t	strlcat(char *dst, const char *src, size_t siz);
#endif

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RTIO_H_ */

