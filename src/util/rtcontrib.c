#ifndef lint
static const char RCSid[] = "$Id: rtcontrib.c,v 1.6 2005/05/28 22:27:54 greg Exp $";
#endif
/*
 * Gather rtrace output to compute contributions from particular sources
 */

#include  "standard.h"
#include  <ctype.h>
#include  "platform.h"
#include  "rtprocess.h"
#include  "selcall.h"
#include  "color.h"
#include  "resolu.h"
#include  "lookup.h"
#include  "calcomp.h"

#define	MAXMODLIST	1024		/* maximum modifiers we'll track */

int	treebufsiz = BUFSIZ;		/* current tree buffer size */

typedef double	DCOLOR[3];		/* double-precision color */

/*
 * The modcont structure is used to accumulate ray contributions
 * for a particular modifier, which may be subdivided into bins
 * if binv is non-NULL.  If outspec contains a %s in it, this will
 * be replaced with the modifier name.  If outspec contains a %d in it,
 * this will be used to create one output file per bin, otherwise all bins
 * will be written to the same file, in order.  If the global outfmt
 * is 'c', then a 4-byte RGBE pixel will be output for each bin value
 * and the file will conform to a RADIANCE image if xres & yres are set.
 */
typedef struct {
	const char	*outspec;	/* output file specification */
	const char	*modname;	/* modifier name */
	EPNODE		*binv;		/* bin value expression */
	int		nbins;		/* number of accumulation bins */
	DCOLOR		cbin[1];	/* contribution bins (extends struct) */
} MODCONT;			/* modifier contribution */

static void mcfree(void *p) { epfree((*(MODCONT *)p).binv); free(p); }

LUTAB	modconttab = LU_SINIT(NULL,mcfree);	/* modifier lookup table */

/* close output stream */
static void closefile(void *p) { fclose((FILE *)p); }

LUTAB	ofiletab = LU_SINIT(free,closefile);	/* output file table */

FILE *getofile(const char *ospec, const char *mname, int bn);

/*
 * The rcont structure is used to manage i/o with a particular
 * rtrace child process.  Input is passed unchanged from stdin,
 * and output is processed in input order and accumulated according
 * to the corresponding modifier and bin number.
 */
struct rtproc {
	struct rtproc	*next;		/* next in list of processes */
	SUBPROC		pd;		/* rtrace pipe descriptors */
	unsigned long	raynum;		/* ray number for this tree */
	int		bsiz;		/* ray tree buffer length */
	char		*buf;		/* ray tree buffer */
	int		nbr;		/* number of bytes from rtrace */
};				/* rtrace process buffer */

					/* rtrace command and defaults */
char		*rtargv[256] = { "rtrace", "-dj", ".5", "-dr", "3",
				"-ab", "1", "-ad", "128", };
int  rtargc = 9;
					/* overriding rtrace options */
char		*myrtopts[] = { "-o~~TmWdp", "-h-", "-x", "1", "-y", "0",
				"-dt", "0", "-as", "0", "-aa", "0", NULL };

struct rtproc	rt0;			/* head of rtrace process list */

struct rtproc	*rt_unproc = NULL;	/* unprocessed ray trees */

char	persistfn[] = "pfXXXXXX";	/* persist file name */

int		gargc;			/* global argc */
char		**gargv;		/* global argv */
#define  progname	gargv[0]

char		*octree;		/* global octree argument */

int		inpfmt = 'a';		/* input format */
int		outfmt = 'a';		/* output format */

int		header = 1;		/* output header? */
int		xres = 0;		/* horiz. output resolution */
int		yres = 0;		/* vert. output resolution */

long		raysleft;		/* number of rays left to trace */
long		waitflush;		/* how long until next flush */

unsigned long	lastray = 0;		/* last ray number sent */
unsigned long	lastdone = 0;		/* last ray processed */

int		using_stdout = 0;	/* are we using stdout? */

const char	*modname[MAXMODLIST];	/* ordered modifier name list */
int		nmods = 0;		/* number of modifiers */

MODCONT *addmodifier(char *modn, char *outf, char *binv);

void init(int np);
int done_rprocs(struct rtproc *rtp);
void trace_contribs(FILE *fp);
struct rtproc *wait_rproc(void);
struct rtproc *get_rproc(void);
void queue_raytree(struct rtproc *rtp);
void process_queue(void);

void putcontrib(const DCOLOR cnt, FILE *fout);
void add_contrib(const char *modn);
void done_contrib(void);

/* set input/output format */
static void
setformat(const char *fmt)
{
	switch (fmt[0]) {
	case 'a':
	case 'f':
	case 'd':
		inpfmt = fmt[0];
		break;
	default:
		goto fmterr;
	}
	switch (fmt[1]) {
	case '\0':
		outfmt = inpfmt;
		return;
	case 'a':
	case 'f':
	case 'd':
	case 'c':
		outfmt = fmt[1];
		break;
	default:
		goto fmterr;
	}
	if (!fmt[2])
		return;
fmterr:
	sprintf(errmsg, "Illegal i/o format: -f%s", fmt);
	error(USER, errmsg);
}

/* gather rays from rtrace and output contributions */
int
main(int argc, char *argv[])
{
	int	nprocs = 1;
	char	*curout = NULL;
	char	*binval = NULL;
	char	fmt[8];
	int	i, j;
				/* global program name */
	gargv = argv;
				/* set up calcomp mode */
	esupport |= E_VARIABLE|E_FUNCTION|E_INCHAN|E_RCONST|E_REDEFW;
	esupport &= ~(E_OUTCHAN);
				/* get our options */
	for (i = 1; i < argc-1; i++) {
						/* expand arguments */
		while ((j = expandarg(&argc, &argv, i)) > 0)
			;
		if (j < 0) {
			fprintf(stderr, "%s: cannot expand '%s'",
					argv[0], argv[i]);
			exit(1);
		}
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'n':		/* number of processes */
				if (argv[i][2] || i >= argc-1) break;
				nprocs = atoi(argv[++i]);
				if (nprocs <= 0)
					error(USER, "illegal number of processes");
				continue;
			case 'h':		/* output header? */
				switch (argv[i][2]) {
				case '\0':
					header = !header;
					continue;
				case '+': case '1':
				case 'T': case 't':
				case 'Y': case 'y':
					header = 1;
					continue;
				case '-': case '0':
				case 'F': case 'f':
				case 'N': case 'n':
					header = 0;
					continue;
				}
				break;
			case 'f':		/* file or i/o format */
				if (!argv[i][2]) {
					if (i >= argc-1) break;
					fcompile(argv[++i]);
					continue;
				}
				setformat(argv[i]+2);
				continue;
			case 'e':		/* expression */
				if (argv[i][2] || i >= argc-1) break;
				scompile(argv[++i], NULL, 0);
				continue;
			case 'o':		/* output file spec. */
				if (argv[i][2] || i >= argc-1) break;
				curout = argv[++i];
				continue;
			case 'x':		/* horiz. output resolution */
				if (argv[i][2] || i >= argc-1) break;
				xres = atoi(argv[++i]);
				continue;
			case 'y':		/* vert. output resolution */
				if (argv[i][2] || i >= argc-1) break;
				yres = atoi(argv[++i]);
				continue;
			case 'b':		/* bin expression */
				if (argv[i][2] || i >= argc-1) break;
				binval = argv[++i];
				continue;
			case 'm':		/* modifier name */
				if (argv[i][2] || i >= argc-1) break;
				rtargv[rtargc++] = "-ti";
				rtargv[rtargc++] = argv[++i];
				addmodifier(argv[i], curout, binval);
				continue;
			}
		rtargv[rtargc++] = argv[i];	/* assume rtrace option */
	}
				/* set global argument list */
	gargc = argc; gargv = argv;
				/* add "mandatory" rtrace settings */
	for (j = 0; myrtopts[j] != NULL; j++)
		rtargv[rtargc++] = myrtopts[j];
				/* just asking for defaults? */
	if (!strcmp(argv[i], "-defaults")) {
		char	sxres[16], syres[16];
		char	*rtpath;
		printf("-n  %-2d\t\t\t\t# number of processes\n", nprocs);
		fflush(stdout);			/* report OUR options */
		rtargv[rtargc++] = header ? "-h+" : "-h-";
		sprintf(fmt, "-f%c%c", inpfmt, outfmt);
		rtargv[rtargc++] = fmt;
		rtargv[rtargc++] = "-x";
		sprintf(sxres, "%d", xres);
		rtargv[rtargc++] = sxres;
		rtargv[rtargc++] = "-y";
		sprintf(syres, "%d", yres);
		rtargv[rtargc++] = syres;
		rtargv[rtargc++] = "-oTW";
		rtargv[rtargc++] = "-defaults";
		rtargv[rtargc] = NULL;
		rtpath = getpath(rtargv[0], getenv("PATH"), X_OK);
		if (rtpath == NULL) {
			eputs(rtargv[0]);
			eputs(": command not found\n");
			exit(1);
		}
		execv(rtpath, rtargv);
		perror(rtpath);	/* execv() should not return */
		exit(1);
	}
	if (nprocs > 1) {	/* add persist file if parallel */
		rtargv[rtargc++] = "-PP";
		rtargv[rtargc++] = mktemp(persistfn);
	} 
				/* add format string */
	sprintf(fmt, "-f%cf", inpfmt);
	rtargv[rtargc++] = fmt;
				/* octree argument is last */
	if (i <= 0 || i != argc-1 || argv[i][0] == '-')
		error(USER, "missing octree argument");
	rtargv[rtargc++] = octree = argv[i];
	rtargv[rtargc] = NULL;
				/* start rtrace & compute contributions */
	init(nprocs);
	trace_contribs(stdin);
	quit(0);
}

/* kill persistent rtrace process */
static void
killpersist(void)
{
	FILE	*fp = fopen(persistfn, "r");
	int	pid;

	if (fp == NULL)
		return;
	if (fscanf(fp, "%*s %d", &pid) != 1 || kill(pid, SIGALRM) < 0)
		unlink(persistfn);
	fclose(fp);
}

/* close rtrace processes and clean up */
int
done_rprocs(struct rtproc *rtp)
{
	int	st0, st1 = 0;

	if (rtp->next != NULL) {	/* close last opened first! */
		st1 = done_rprocs(rtp->next);
		free((void *)rtp->next);
		rtp->next = NULL;
	}
	st0 = close_process(&rtp->pd);
	if (st0 < 0)
		error(WARNING, "unknown return status from rtrace process");
	else if (st0 > 0)
		return(st0);
	return(st1);
}

/* exit with status */
void
quit(int status)
{
	int	rtstat;

	if (rt0.next != NULL)		/* terminate persistent rtrace */
		killpersist();
					/* clean up rtrace process(es) */
	rtstat = done_rprocs(&rt0);
	if (status == 0)
		status = rtstat;
	exit(status);			/* flushes all output streams */
}

/* start rtrace processes and initialize */
void
init(int np)
{
	struct rtproc	*rtp;
	int	i;
	int	maxbytes;
					/* make sure we have something to do */
	if (!nmods)
		error(USER, "No modifiers specified");
					/* assign ray variables */
	scompile("Dx=$1;Dy=$2;Dz=$3;", NULL, 0);
	scompile("Px=$4;Py=$5;Pz=$6;", NULL, 0);
					/* set up signal handling */
	signal(SIGINT, quit);
#ifdef SIGHUP
	signal(SIGHUP, quit);
#endif
#ifdef SIGTERM
	signal(SIGTERM, quit);
#endif
#ifdef SIGPIPE
	signal(SIGPIPE, quit);
#endif
	rtp = &rt0;			/* start rtrace process(es) */
	for (i = 0; i++ < np; ) {
		errno = 0;
		maxbytes = open_process(&rtp->pd, rtargv);
		if (maxbytes == 0) {
			eputs(rtargv[0]);
			eputs(": command not found\n");
			exit(1);
		}
		if (maxbytes < 0)
			error(SYSTEM, "cannot start rtrace process");
		if (maxbytes > treebufsiz)
			treebufsiz = maxbytes;
		rtp->raynum = 0;
		rtp->bsiz = 0;
		rtp->buf = NULL;
		rtp->nbr = 0;
		if (i == np)		/* last process? */
			break;
		if (i == 1)
			sleep(2);	/* wait for persist file */
		rtp->next = (struct rtproc *)malloc(sizeof(struct rtproc));
		if (rtp->next == NULL)
			error(SYSTEM, "out of memory in init");
		rtp = rtp->next;
	}
	rtp->next = NULL;		/* terminate list */
	if (yres > 0) {
		if (xres > 0)
			raysleft = xres*yres;
		else
			raysleft = yres;
	} else
		raysleft = 0;
	waitflush = xres;
}

/* add modifier to our list to track */
MODCONT *
addmodifier(char *modn, char *outf, char *binv)
{
	LUENT	*lep = lu_find(&modconttab, modn);
	MODCONT	*mp;
	
	if (lep->data != NULL) {
		sprintf(errmsg, "duplicate modifier '%s'", modn);
		error(USER, errmsg);
	}
	if (nmods >= MAXMODLIST)
		error(USER, "too many modifiers");
	modname[nmods++] = modn;	/* XXX assumes static string */
	lep->key = modn;		/* XXX assumes static string */
	mp = (MODCONT *)malloc(sizeof(MODCONT));
	if (mp == NULL)
		error(SYSTEM, "out of memory in addmodifier");
	lep->data = (char *)mp;
	mp->outspec = outf;		/* XXX assumes static string */
	mp->modname = modn;		/* XXX assumes static string */
	if (binv != NULL)
		mp->binv = eparse(binv);
	else
		mp->binv = eparse("0");
	mp->nbins = 1;
	setcolor(mp->cbin[0], 0., 0., 0.);
	return mp;
}

/* put string to stderr */
void
eputs(char  *s)
{
	static int  midline = 0;

	if (!*s) return;
	if (!midline) {
		fputs(progname, stderr);
		fputs(": ", stderr);
	}
	fputs(s, stderr);
	midline = s[strlen(s)-1] != '\n';
}

/* write header to the given output stream */
void
printheader(FILE *fout)
{
	extern char	VersionID[];
	FILE		*fin = fopen(octree, "r");
	
	if (fin == NULL)
		quit(1);
	checkheader(fin, "ignore", fout);	/* copy octree header */
	fclose(fin);
	printargs(gargc-1, gargv, fout);	/* add our command */
	fprintf(fout, "SOFTWARE= %s\n", VersionID);
	fputnow(fout);
	switch (outfmt) {			/* add output format */
	case 'a':
		fputformat("ascii", fout);
		break;
	case 'f':
		fputformat("float", fout);
		break;
	case 'd':
		fputformat("double", fout);
		break;
	case 'c':
		fputformat(COLRFMT, fout);
		break;
	}
	fputc('\n', fout);
	if (xres > 0) {
		if (yres > 0)			/* resolution string */
			fprtresolu(xres, yres, fout);
		fflush(fout);
	}
}

/* Get output file pointer (open and write header if new) */
FILE *
getofile(const char *ospec, const char *mname, int bn)
{
	const char	*mnp = NULL;
	const char	*bnp = NULL;
	const char	*cp;
	char		ofname[1024];
	LUENT		*lep;
	
	if (ospec == NULL) {			/* use stdout? */
		if (!using_stdout && header)
			printheader(stdout);
		using_stdout = 1;
		return stdout;
	}
	for (cp = ospec; *cp; cp++)		/* check format position(s) */
		if (*cp == '%') {
			do
				++cp;
			while (isdigit(*cp));
			switch (*cp) {
			case '%':
				break;
			case 's':
				if (mnp != NULL)
					goto badspec;
				mnp = cp;
				break;
			case 'd':
				if (bnp != NULL)
					goto badspec;
				bnp = cp;
				break;
			default:
				goto badspec;
			}
		}
	if (mnp != NULL) {			/* create file name */
		if (bnp != NULL) {
			if (bnp > mnp)
				sprintf(ofname, ospec, mname, bn);
			else
				sprintf(ofname, ospec, bn, mname);
		} else
			sprintf(ofname, ospec, mname);
	} else if (bnp != NULL)
		sprintf(ofname, ospec, bn);
	else
		strcpy(ofname, ospec);
	lep = lu_find(&ofiletab, ofname);	/* look it up */
	if (lep->key == NULL)			/* new entry */
		lep->key = strcpy((char *)malloc(strlen(ofname)+1), ofname);
	if (lep->data == NULL) {		/* open output file */
		FILE		*fp = fopen(ofname, "w");
		int		i;
		if (fp == NULL) {
			sprintf(errmsg, "cannot open '%s' for writing", ofname);
			error(SYSTEM, errmsg);
		}
		if (header)
			printheader(fp);
						/* play catch-up */
		for (i = 0; i < lastdone; i++) {
			static const DCOLOR	nocontrib = BLKCOLOR;
			putcontrib(nocontrib, fp);
			if (outfmt == 'a')
				putc('\n', fp);
		}
		if (xres > 0)
			fflush(fp);
		lep->data = (char *)fp;
	}
	return (FILE *)lep->data;		/* return open file pointer */
badspec:
	sprintf(errmsg, "bad output format '%s'", ospec);
	error(USER, errmsg);
	return NULL;		/* pro forma return */
}

/* read input ray into buffer */
int
getinp(char *buf, FILE *fp)
{
	char	*cp;
	int	i;

	switch (inpfmt) {
	case 'a':
		cp = buf;		/* make sure we get 6 floats */
		for (i = 0; i < 6; i++) {
			if (fgetword(cp, buf+127-cp, fp) == NULL)
				return 0;
			if ((cp = fskip(cp)) == NULL || *cp)
				return 0;
			*cp++ = ' ';
		}
		getc(fp);		/* get/put eol */
		*cp-- = '\0'; *cp = '\n';
		return strlen(buf);
	case 'f':
		if (fread(buf, sizeof(float), 6, fp) < 6)
			return 0;
		return sizeof(float)*6;
	case 'd':
		if (fread(buf, sizeof(double), 6, fp) < 6)
			return 0;
		return sizeof(double)*6;
	}
	error(INTERNAL, "botched input format");
	return 0;	/* pro forma return */
}

static float	rparams[9];		/* traced ray parameters */

/* return channel (ray) value */
double
chanvalue(int n)
{
	if (--n < 0 || n >= 6)
		error(USER, "illegal channel number ($N)");
	return rparams[n+3];
}

/* add current ray contribution to the appropriate modifier bin */
void
add_contrib(const char *modn)
{
	LUENT	*le = lu_find(&modconttab, modn);
	MODCONT	*mp = (MODCONT *)le->data;
	int	bn;

	if (mp == NULL) {
		sprintf(errmsg, "unexpected modifier '%s' from rtrace", modn);
		error(USER, errmsg);
	}
	eclock++;			/* get bin number */
	bn = (int)(evalue(mp->binv) + .5);
	if (bn <= 0)
		bn = 0;
	else if (bn > mp->nbins) {	/* new bin */
		mp = (MODCONT *)realloc(mp, sizeof(MODCONT) +
						bn*sizeof(DCOLOR));
		if (mp == NULL)
			error(SYSTEM, "out of memory in add_contrib");
		memset(mp->cbin+mp->nbins, 0, sizeof(DCOLOR)*(bn+1-mp->nbins));
		mp->nbins = bn+1;
		le->data = (char *)mp;
	}
	addcolor(mp->cbin[bn], rparams);
}

/* output newline to ASCII file and/or flush as requested */
static int
puteol(const LUENT *e, void *p)
{
	FILE	*fp = (FILE *)e->data;

	if (outfmt == 'a')
		putc('\n', fp);
	if (!waitflush)
		fflush(fp);
	if (ferror(fp)) {
		sprintf(errmsg, "write error on file '%s'", e->key);
		error(SYSTEM, errmsg);
	}
	return 0;
}

/* put out ray contribution to file */
void
putcontrib(const DCOLOR cnt, FILE *fout)
{
	float	fv[3];
	COLR	cv;

	switch (outfmt) {
	case 'a':
		fprintf(fout, "%.6e\t%.6e\t%.6e\t", cnt[0], cnt[1], cnt[2]);
		break;
	case 'f':
		fv[0] = cnt[0];
		fv[1] = cnt[1];
		fv[2] = cnt[2];
		fwrite(fv, sizeof(float), 3, fout);
		break;
	case 'd':
		fwrite(cnt, sizeof(double), 3, fout);
		break;
	case 'c':
		setcolr(cv, cnt[0], cnt[1], cnt[2]);
		fwrite(cv, sizeof(cv), 1, fout);
		break;
	default:
		error(INTERNAL, "botched output format");
	}
}

/* output ray tallies and clear for next primary */
void
done_contrib(void)
{
	int	i, j;
	MODCONT	*mp;
						/* output modifiers in order */
	for (i = 0; i < nmods; i++) {
		mp = (MODCONT *)lu_find(&modconttab,modname[i])->data;
		for (j = 0; j < mp->nbins; j++)
			putcontrib(mp->cbin[j],
					getofile(mp->outspec, mp->modname, j));
						/* clear for next ray tree */
		memset(mp->cbin, 0, sizeof(DCOLOR)*mp->nbins);
	}
	--waitflush;				/* terminate records */
	lu_doall(&ofiletab, puteol, NULL);
	if (using_stdout & (outfmt == 'a'))
		putc('\n', stdout);
	if (!waitflush) {
		waitflush = xres;
		if (using_stdout)
			fflush(stdout);
	}
}

/* queue completed ray tree produced by rtrace process */
void
queue_raytree(struct rtproc *rtp)
{
	struct rtproc	*rtu, *rtl = NULL;
					/* insert following ray order */
	for (rtu = rt_unproc; rtu != NULL; rtu = (rtl=rtu)->next)
		if (rtp->raynum < rtu->raynum)
			break;
	rtu = (struct rtproc *)malloc(sizeof(struct rtproc));
	if (rtu == NULL)
		error(SYSTEM, "out of memory in queue_raytree");
	*rtu = *rtp;
	if (rtl == NULL) {
		rtu->next = rt_unproc;
		rt_unproc = rtu;
	} else {
		rtu->next = rtl->next;
		rtl->next = rtu;
	}
	rtp->raynum = 0;		/* clear path for next ray tree */
	rtp->bsiz = 0;
	rtp->buf = NULL;
	rtp->nbr = 0;
}

/* process completed ray trees from our queue */
void
process_queue(void)
{
	char	modname[128];
					/* ray-ordered queue */
	while (rt_unproc != NULL && rt_unproc->raynum == lastdone+1) {
		struct rtproc	*rtp = rt_unproc;
		int		n = rtp->nbr;
		const char	*cp = rtp->buf;
		while (n > 0) {		/* process rays */
			register char	*mnp = modname;
					/* skip leading tabs */
			while (n > 0 && *cp == '\t') {
				cp++; n--;
			}
			if (!n || !(isalpha(*cp) | (*cp == '_')))
				error(USER, "bad modifier name from rtrace");
					/* get modifier name */
			while (n > 0 && *cp != '\t') {
				*mnp++ = *cp++; n--;
			}
			*mnp = '\0';
			cp++; n--;	/* eat following tab */
			if (n < (int)(sizeof(float)*9))
				error(USER, "incomplete ray value from rtrace");
					/* add ray contribution */
			memcpy(rparams, cp, sizeof(float)*9);
			cp += sizeof(float)*9; n -= sizeof(float)*9;
			add_contrib(modname);
		}
		done_contrib();		/* sum up contributions & output */
		lastdone = rtp->raynum;
		free(rtp->buf);		/* free up buffer space */
		rt_unproc = rtp->next;
		free(rtp);		/* done with this ray tree */
	}
}

/* wait for rtrace process to finish with ray tree */
struct rtproc *
wait_rproc(void)
{
	struct rtproc	*rtfree = NULL;
	fd_set		readset, errset;
	int		nr;
	struct rtproc	*rt;
	int		n;
	
	do {
		nr = 0;				/* prepare select call */
		FD_ZERO(&readset); FD_ZERO(&errset); n = 0;
		for (rt = &rt0; rt != NULL; rt = rt->next) {
			if (rt->raynum) {
				FD_SET(rt->pd.r, &readset);
				++nr;
			}
			FD_SET(rt->pd.r, &errset);
			if (rt->pd.r >= n)
				n = rt->pd.r + 1;
		}
		if (!nr)			/* no rays pending */
			break;
		if (nr > 1) {			/* call select for multiple proc's */
			errno = 0;
			if (select(n, &readset, NULL, &errset, NULL) < 0)
				error(SYSTEM, "select call error in wait_rproc()");
		} else
			FD_ZERO(&errset);
		nr = 0;
		for (rt = &rt0; rt != NULL; rt = rt->next) {
			if (!FD_ISSET(rt->pd.r, &readset) &&
					!FD_ISSET(rt->pd.r, &errset))
				continue;
			if (rt->buf == NULL) {
				rt->bsiz = treebufsiz;
				rt->buf = (char *)malloc(treebufsiz);
			} else if (rt->nbr + BUFSIZ > rt->bsiz) {
				if (rt->bsiz + BUFSIZ <= treebufsiz)
					rt->bsiz = treebufsiz;
				else
					rt->bsiz = treebufsiz += BUFSIZ;
				rt->buf = (char *)realloc(rt->buf, rt->bsiz);
			}
			if (rt->buf == NULL)
				error(SYSTEM, "out of memory in wait_rproc");
			nr = read(rt->pd.r, rt->buf+rt->nbr, rt->bsiz-rt->nbr);
			if (nr <= 0)
				error(USER, "rtrace process died");
			rt->nbr += nr;		/* advance & check */
			if (rt->nbr >= 4 && !memcmp(rt->buf+rt->nbr-4,
							"~\t~\t", 4)) {
				rt->nbr -= 4;	/* elide terminator */
				queue_raytree(rt);
				rtfree = rt;	/* ready for next ray */
			}
		}
	} while ((rtfree == NULL) & (nr > 0));	/* repeat until ready or out */
	return rtfree;
}

/* return next available rtrace process */
struct rtproc *
get_rproc(void)
{
	struct rtproc	*rtp;
						/* check for idle rtrace */
	for (rtp = &rt0; rtp != NULL; rtp = rtp->next)
		if (!rtp->raynum)
			return rtp;
	return wait_rproc();			/* need to wait for one */
}

/* trace ray contributions (main loop) */
void
trace_contribs(FILE *fin)
{
	char		inpbuf[128];
	int		iblen;
	struct rtproc	*rtp;
						/* loop over input */
	while ((iblen = getinp(inpbuf, fin)) > 0) {
		if (lastray+1 < lastray) {	/* counter rollover? */
			while (wait_rproc() != NULL)
				process_queue();
			lastdone = lastray = 0;
		}
		rtp = get_rproc();		/* get avail. rtrace process */
		rtp->raynum = ++lastray;	/* assign ray to it */
		writebuf(rtp->pd.w, inpbuf, iblen);
		if (!--raysleft)
			break;
		process_queue();		/* catch up with results */
	}
	while (wait_rproc() != NULL)		/* process outstanding rays */
		process_queue();
	if (raysleft > 0)
		error(USER, "unexpected EOF on input");
}
