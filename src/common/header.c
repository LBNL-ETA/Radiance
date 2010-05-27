#ifndef lint
static const char	RCSid[] = "$Id: header.c,v 2.26 2010/05/27 19:32:12 greg Exp $";
#endif
/*
 *  header.c - routines for reading and writing information headers.
 *
 *  Externals declared in resolu.h
 *
 *  newheader(t,fp)	start new information header identified by string t
 *  headidval(r,s)	copy header identifier value in s to r
 *  dateval(t,s)	get capture date value as UTC
 *  gmtval(t,s)		get GMT as UTC
 *  fputdate(t,fp)	put out the given UTC
 *  fputnow(fp)		put out the current date and time
 *  printargs(ac,av,fp) print an argument list to fp, followed by '\n'
 *  formatval(r,s)	copy the format value in s to r
 *  fputformat(s,fp)	write "FORMAT=%s" to fp
 *  getheader(fp,f,p)	read header from fp, calling f(s,p) on each line
 *  globmatch(pat, str)	check for glob match of str against pat
 *  checkheader(i,p,o)	check header format from i against p and copy to o
 *
 *  To copy header from input to output, use getheader(fin, fputs, fout)
 */

#include "copyright.h"

#include  <ctype.h>

#include  "rtio.h"
#include  "resolu.h"

#define	 MAXLINE	2048

extern time_t		timegm(struct tm *tm);

const char  HDRSTR[] = "#?";		/* information header magic number */

const char  FMTSTR[] = "FORMAT=";	/* format identifier */

const char  TMSTR[] = "CAPDATE=";	/* capture date identifier */
const char  GMTSTR[] = "GMT=";		/* GMT identifier */

static gethfunc mycheck;


extern void
newheader(		/* identifying line of information header */
	char  *s,
	FILE  *fp
)
{
	fputs(HDRSTR, fp);
	fputs(s, fp);
	putc('\n', fp);
}


extern int
headidval(			/* get header id (return true if is id) */
	char  *r,
	char	*s
)
{
	const char  *cp = HDRSTR;

	while (*cp) if (*cp++ != *s++) return(0);
	if (r == NULL) return(1);
	while (*s && !isspace(*s)) *r++ = *s++;
	*r = '\0';
	return(1);
}


extern int
dateval(		/* convert capture date line to UTC */
	time_t	*tloc,
	char	*s
)
{
	struct tm	tms;
	const char	*cp = TMSTR;

	while (*cp) if (*cp++ != *s++) return(0);
	while (isspace(*s)) s++;
	if (!*s) return(0);
	if (sscanf(s, "%d:%d:%d %d:%d:%d",
			&tms.tm_year, &tms.tm_mon, &tms.tm_mday,
			&tms.tm_hour, &tms.tm_min, &tms.tm_sec) != 6)
		return(0);
	if (tloc == NULL)
		return(1);
	tms.tm_mon--;
	tms.tm_year -= 1900;
	tms.tm_isdst = -1;	/* ask mktime() to figure out DST */
	*tloc = mktime(&tms);
	return(1);
}


extern int
gmtval(			/* convert GMT date line to UTC */
	time_t	*tloc,
	char	*s
)
{
	struct tm	tms;
	const char	*cp = GMTSTR;

	while (*cp) if (*cp++ != *s++) return(0);
	while (isspace(*s)) s++;
	if (!*s) return(0);
	if (sscanf(s, "%d:%d:%d %d:%d:%d",
			&tms.tm_year, &tms.tm_mon, &tms.tm_mday,
			&tms.tm_hour, &tms.tm_min, &tms.tm_sec) != 6)
		return(0);
	if (tloc == NULL)
		return(1);
	tms.tm_mon--;
	tms.tm_year -= 1900;
	*tloc = timegm(&tms);
	return(1);
}


extern void
fputdate(		/* write out the given time value (local & GMT) */
	time_t	tv,
	FILE	*fp
)
{
	struct tm	*tms;

	tms = localtime(&tv);
	if (tms != NULL)
		fprintf(fp, "%s %04d:%02d:%02d %02d:%02d:%02d\n", TMSTR,
				tms->tm_year+1900, tms->tm_mon+1, tms->tm_mday,
				tms->tm_hour, tms->tm_min, tms->tm_sec);
	tms = gmtime(&tv);
	if (tms != NULL)
		fprintf(fp, "%s %04d:%02d:%02d %02d:%02d:%02d\n", GMTSTR,
				tms->tm_year+1900, tms->tm_mon+1, tms->tm_mday,
				tms->tm_hour, tms->tm_min, tms->tm_sec);
}


extern void
fputnow(			/* write out the current time */
	FILE	*fp
)
{
	time_t	tv;
	time(&tv);
	fputdate(tv, fp);
}


extern void
printargs(		/* print arguments to a file */
	int  ac,
	char  **av,
	FILE  *fp
)
{
	while (ac-- > 0) {
		fputword(*av++, fp);
		fputc(ac ? ' ' : '\n', fp);
	}
}


extern int
formatval(			/* get format value (return true if format) */
	char  *r,
	char  *s
)
{
	const char  *cp = FMTSTR;

	while (*cp) if (*cp++ != *s++) return(0);
	while (isspace(*s)) s++;
	if (!*s) return(0);
	if (r == NULL) return(1);
	do
		*r++ = *s++;
	while(*s && !isspace(*s));
	*r = '\0';
	return(1);
}


extern void
fputformat(		/* put out a format value */
	char  *s,
	FILE  *fp
)
{
	fputs(FMTSTR, fp);
	fputs(s, fp);
	putc('\n', fp);
}


extern int
getheader(		/* get header from file */
	FILE  *fp,
	gethfunc *f,
	void  *p
)
{
	char  buf[MAXLINE];

	for ( ; ; ) {
		buf[MAXLINE-2] = '\n';
		if (fgets(buf, MAXLINE, fp) == NULL)
			return(-1);
		if (buf[0] == '\n')
			return(0);
#ifdef MSDOS
		if (buf[0] == '\r' && buf[1] == '\n')
			return(0);
#endif
		if (buf[MAXLINE-2] != '\n') {
			ungetc(buf[MAXLINE-2], fp);	/* prevent false end */
			buf[MAXLINE-2] = '\0';
		}
		if (f != NULL && (*f)(buf, p) < 0)
			return(-1);
	}
}


struct check {
	FILE	*fp;
	char	fs[64];
};


static int
mycheck(			/* check a header line for format info. */
	char  *s,
	void  *cp
)
{
	if (!formatval(((struct check*)cp)->fs, s)
			&& ((struct check*)cp)->fp != NULL) {
		fputs(s, ((struct check*)cp)->fp);
	}
	return(0);
}


extern int
globmatch(			/* check for match of s against pattern p */
	char	*p,
	char	*s
)
{
	int	setmatch;

	do {
		switch (*p) {
		case '?':			/* match any character */
			if (!*s++)
				return(0);
			break;
		case '*':			/* match any string */
			while (p[1] == '*') p++;
			do
				if ( (p[1]=='?' || p[1]==*s) &&
						globmatch(p+1,s) )
					return(1);
			while (*s++);
			return(0);
		case '[':			/* character set */
			setmatch = *s == *++p;
			if (!*p)
				return(0);
			while (*++p != ']') {
				if (!*p)
					return(0);
				if (*p == '-') {
					setmatch += p[-1] <= *s && *s <= p[1];
					if (!*++p)
						break;
				} else
					setmatch += *p == *s;
			}
			if (!setmatch)
				return(0);
			s++;
			break;
		case '\\':			/* literal next */
			p++;
		/* fall through */
		default:			/* normal character */
			if (*p != *s)
				return(0);
			s++;
			break;
		}
	} while (*p++);
	return(1);
}


/*
 * Checkheader(fin,fmt,fout) returns a value of 1 if the input format
 * matches the specification in fmt, 0 if no input format was found,
 * and -1 if the input format does not match or there is an
 * error reading the header.  If fmt is empty, then -1 is returned
 * if any input format is found (or there is an error), and 0 otherwise.
 * If fmt contains any '*' or '?' characters, then checkheader
 * does wildcard expansion and copies a matching result into fmt.
 * Be sure that fmt is big enough to hold the match in such cases,
 * and that it is not a static, read-only string!
 * The input header (minus any format lines) is copied to fout
 * if fout is not NULL.
 */

extern int
checkheader(
	FILE  *fin,
	char  *fmt,
	FILE  *fout
)
{
	struct check	cdat;
	char	*cp;

	cdat.fp = fout;
	cdat.fs[0] = '\0';
	if (getheader(fin, mycheck, &cdat) < 0)
		return(-1);
	if (!cdat.fs[0])
		return(0);
	for (cp = fmt; *cp; cp++)		/* check for globbing */
		if ((*cp == '?') | (*cp == '*')) {
			if (globmatch(fmt, cdat.fs)) {
				strcpy(fmt, cdat.fs);
				return(1);
			} else
				return(-1);
		}
	return(strcmp(fmt, cdat.fs) ? -1 : 1);	/* literal match */
}
