#ifndef lint
static const char RCSid[] = "$Id: rtcontrib.c,v 1.58 2010/10/05 18:05:22 greg Exp $";
#endif
/*
 * Gather rtrace output to compute contributions from particular sources
 */

#include  "standard.h"
#include  <ctype.h>
#include  <signal.h>
#include  "platform.h"
#include  "rtprocess.h"
#include  "selcall.h"
#include  "color.h"
#include  "resolu.h"
#include  "lookup.h"
#include  "random.h"
#include  "calcomp.h"

#ifndef	MAXMODLIST
#define	MAXMODLIST	1024		/* maximum modifiers we'll track */
#endif

size_t	treebufsiz = BUFSIZ;		/* current tree buffer size */

typedef double	DCOLOR[3];		/* double-precision color */

/*
 * The MODCONT structure is used to accumulate ray contributions
 * for a particular modifier, which may be subdivided into bins
 * if binv evaluates > 0.  If outspec contains a %s in it, this will
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

/*
 * The STREAMOUT structure holds an open FILE pointer and a count of
 * the number of RGB triplets per record, or 0 if unknown.
 */
typedef struct {
	FILE		*ofp;		/* output file pointer */
	int		outpipe;	/* output is to a pipe */
	int		reclen;		/* triplets/record */
	int		xr, yr;		/* output resolution for picture */
} STREAMOUT;

/* close output stream and free record */
static void
closestream(void *p)
{
	STREAMOUT	*sop = (STREAMOUT *)p;
	int		status;
	if (sop->outpipe)
		status = pclose(sop->ofp);
	else
		status = fclose(sop->ofp);
	if (status)
		error(SYSTEM, "error closing output stream");
	free(p);
}

LUTAB	ofiletab = LU_SINIT(free,closestream);	/* output file table */

#define OF_MODIFIER	01
#define OF_BIN		02

STREAMOUT *getostream(const char *ospec, const char *mname, int bn, int noopen);
int ofname(char *oname, const char *ospec, const char *mname, int bn);
void printheader(FILE *fout, const char *info);
void printresolu(FILE *fout, int xr, int yr);

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
	size_t		bsiz;		/* ray tree buffer length */
	char		*buf;		/* ray tree buffer */
	size_t		nbr;		/* number of bytes from rtrace */
};				/* rtrace process buffer */

					/* rtrace command and defaults */
char		*rtargv[256+2*MAXMODLIST] = { "rtrace",
				"-dj", ".9", "-dr", "3",
				"-ab", "1", "-ad", "350", };

int  rtargc = 9;
					/* overriding rtrace options */
char		*myrtopts[] = { "-h-", "-x", "1", "-y", "0",
				"-dt", "0", "-as", "0", "-aa", "0", NULL };

#define	RTCOEFF		"-o~~TmWdp"	/* compute coefficients only */
#define RTCONTRIB	"-o~~TmVdp"	/* compute ray contributions */

struct rtproc	rt0;			/* head of rtrace process list */

struct rtproc	*rt_unproc = NULL;	/* unprocessed ray trees */

#define PERSIST_NONE	0		/* no persist file */
#define PERSIST_SINGLE	1		/* user set -P persist */
#define PERSIST_PARALL	2		/* user set -PP persist */
#define PERSIST_OURS	3		/* -PP persist belongs to us */
int	persist_state = PERSIST_NONE;	/* persist file state */
char	persistfn[] = "pfXXXXXX";	/* our persist file name, if set */

int		gargc;			/* global argc */
char		**gargv;		/* global argv */
#define  progname	gargv[0]

char		*octree;		/* global octree argument */

int		inpfmt = 'a';		/* input format */
int		outfmt = 'a';		/* output format */

int		header = 1;		/* output header? */
int		force_open = 0;		/* truncate existing output? */
int		recover = 0;		/* recover previous output? */
int		accumulate = 1;		/* input rays per output record */
int		xres = 0;		/* horiz. output resolution */
int		yres = 0;		/* vert. output resolution */

int		account;		/* current accumulation count */
unsigned long	raysleft;		/* number of rays left to trace */
long		waitflush;		/* how long until next flush */

unsigned long	lastray = 0;		/* last ray number sent */
unsigned long	lastdone = 0;		/* last ray processed */

int		using_stdout = 0;	/* are we using stdout? */

const char	*modname[MAXMODLIST];	/* ordered modifier name list */
int		nmods = 0;		/* number of modifiers */

#define queue_length()		(lastray - lastdone)

MODCONT *growmodifier(MODCONT *mp, int nb);
MODCONT *addmodifier(char *modn, char *outf, char *binv, int bincnt);
void addmodfile(char *fname, char *outf, char *binv, int bincnt);

void init(int np);
int done_rprocs(struct rtproc *rtp);
void reload_output(void);
void recover_output(FILE *fin);
void trace_contribs(FILE *fin);
struct rtproc *wait_rproc(void);
struct rtproc *get_rproc(void);
void queue_raytree(struct rtproc *rtp);
void process_queue(void);

void put_contrib(const DCOLOR cnt, FILE *fout);
void add_contrib(const char *modn);
void done_contrib(int navg);

/* return number of open rtrace processes */
static int
nrtprocs(void)
{
	int	nrtp = 0;
	struct rtproc	*rtp;

	for (rtp = &rt0; rtp != NULL; rtp = rtp->next)
		nrtp += rtp->pd.running;
	return(nrtp);
}

/* set input/output format */
static void
setformat(const char *fmt)
{
	switch (fmt[0]) {
	case 'f':
	case 'd':
		SET_FILE_BINARY(stdin);
		/* fall through */
	case 'a':
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
	int	contrib = 0;
	int	nprocs = 1;
	char	*curout = NULL;
	char	*binval = NULL;
	int	bincnt = 0;
	char	fmt[8];
	int	i, j;
				/* need at least one argument */
	if (argc < 2) {
		fprintf(stderr,
"Usage: %s [-n nprocs][-V][-r][-e expr][-f source][-o ospec][-b binv] {-m mod | -M file} [rtrace options] octree\n",
			argv[0]);
		exit(1);
	}
				/* global program name */
	gargv = argv;
				/* initialize calcomp routines */
	esupport |= E_VARIABLE|E_FUNCTION|E_INCHAN|E_RCONST|E_REDEFW;
	esupport &= ~(E_OUTCHAN);
	varset("PI", ':', PI);
				/* get our options */
	for (i = 1; i < argc-1; i++) {
						/* expand arguments */
		while ((j = expandarg(&argc, &argv, i)) > 0)
			;
		if (j < 0) {
			fprintf(stderr, "%s: cannot expand '%s'\n",
					argv[0], argv[i]);
			exit(1);
		}
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'n':		/* number of processes */
				if (argv[i][2] || i >= argc-2) break;
				nprocs = atoi(argv[++i]);
				if (nprocs <= 0)
					error(USER, "illegal number of processes");
				continue;
			case 'V':		/* output contributions */
				switch (argv[i][2]) {
				case '\0':
					contrib = !contrib;
					continue;
				case '+': case '1':
				case 'T': case 't':
				case 'Y': case 'y':
					contrib = 1;
					continue;
				case '-': case '0':
				case 'F': case 'f':
				case 'N': case 'n':
					contrib = 0;
					continue;
				}
				break;
			case 'c':		/* input rays per output */
				if (argv[i][2] || i >= argc-2) break;
				accumulate = atoi(argv[++i]);
				continue;
			case 'r':		/* recover output */
				if (argv[i][2]) break;
				recover = 1;
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
			case 'f':		/* file or force or format */
				if (!argv[i][2]) {
					char	*fpath;
					if (i >= argc-2) break;
					fpath = getpath(argv[++i],
							getrlibpath(), R_OK);
					if (fpath == NULL) {
						sprintf(errmsg,
							"cannot find file '%s'",
								argv[i]);
						error(USER, errmsg);
					}
					fcompile(fpath);
					continue;
				}
				if (argv[i][2] == 'o') {
					force_open = 1;
					continue;
				}
				setformat(argv[i]+2);
				continue;
			case 'e':		/* expression */
				if (argv[i][2] || i >= argc-2) break;
				scompile(argv[++i], NULL, 0);
				continue;
			case 'o':		/* output file spec. */
				if (argv[i][2] || i >= argc-2) break;
				curout = argv[++i];
				continue;
			case 'x':		/* horiz. output resolution */
				if (argv[i][2] || i >= argc-2) break;
				xres = atoi(argv[++i]);
				continue;
			case 'y':		/* vert. output resolution */
				if (argv[i][2] || i >= argc-2) break;
				yres = atoi(argv[++i]);
				continue;
			case 'b':		/* bin expression/count */
				if (i >= argc-2) break;
				if (argv[i][2] == 'n') {
					bincnt = (int)(eval(argv[++i]) + .5);
					continue;
				}
				if (argv[i][2]) break;
				binval = argv[++i];
				continue;
			case 'm':		/* modifier name */
				if (argv[i][2] || i >= argc-2) break;
				rtargv[rtargc++] = "-ti";
				rtargv[rtargc++] = argv[++i];
				addmodifier(argv[i], curout, binval, bincnt);
				continue;
			case 'M':		/* modifier file */
				if (argv[i][2] || i >= argc-2) break;
				rtargv[rtargc++] = "-tI";
				rtargv[rtargc++] = argv[++i];
				addmodfile(argv[i], curout, binval, bincnt);
				continue;
			case 'P':		/* persist file */
				if (i >= argc-2) break;
				persist_state = (argv[i][2] == 'P') ?
						PERSIST_PARALL : PERSIST_SINGLE;
				rtargv[rtargc++] = argv[i];
				rtargv[rtargc++] = argv[++i];
				continue;
			}
		rtargv[rtargc++] = argv[i];	/* assume rtrace option */
	}
	if (accumulate <= 0)	/* no output flushing for single record */
		xres = yres = 0;
				/* set global argument list */
	gargc = argc; gargv = argv;
				/* add "mandatory" rtrace settings */
	for (j = 0; myrtopts[j] != NULL; j++)
		rtargv[rtargc++] = myrtopts[j];
	rtargv[rtargc++] = contrib ? RTCONTRIB : RTCOEFF;
				/* just asking for defaults? */
	if (!strcmp(argv[i], "-defaults")) {
		char	nps[8], sxres[16], syres[16];
		char	*rtpath;
		printf("-c %-5d\t\t\t# accumulated rays per record\n",
				accumulate);
		printf("-V%c\t\t\t\t# output %s\n", contrib ? '+' : '-',
				contrib ? "contributions" : "coefficients");
		fflush(stdout);			/* report OUR options */
		rtargv[rtargc++] = "-n";
		sprintf(nps, "%d", nprocs);
		rtargv[rtargc++] = nps;
		rtargv[rtargc++] = header ? "-h+" : "-h-";
		sprintf(fmt, "-f%c%c", inpfmt, outfmt);
		rtargv[rtargc++] = fmt;
		rtargv[rtargc++] = "-x";
		sprintf(sxres, "%d", xres);
		rtargv[rtargc++] = sxres;
		rtargv[rtargc++] = "-y";
		sprintf(syres, "%d", yres);
		rtargv[rtargc++] = syres;
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
		if (persist_state == PERSIST_SINGLE)
			error(USER, "use -PP option for multiple processes");
		if (persist_state == PERSIST_NONE) {
			rtargv[rtargc++] = "-PP";
			rtargv[rtargc++] = mktemp(persistfn);
			persist_state = PERSIST_OURS;
		}
	} 
				/* add format string */
	sprintf(fmt, "-f%cf", inpfmt);
	rtargv[rtargc++] = fmt;
				/* octree argument is last */
	if (i <= 0 || i != argc-1 || argv[i][0] == '-')
		error(USER, "missing octree argument");
	rtargv[rtargc++] = octree = argv[i];
	rtargv[rtargc] = NULL;
				/* start rtrace & recover if requested */
	init(nprocs);
				/* compute contributions */
	trace_contribs(stdin);
				/* clean up */
	quit(0);
}

#ifndef SIGALRM
#define SIGALRM SIGTERM
#endif
/* kill persistent rtrace process */
static void
killpersist(void)
{
	FILE	*fp = fopen(persistfn, "r");
	RT_PID	pid;

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

	if (persist_state == PERSIST_OURS)  /* terminate waiting rtrace */
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
			raysleft = (unsigned long)xres*yres;
		else
			raysleft = yres;
	} else
		raysleft = 0;
	if ((account = accumulate) > 0)
		raysleft *= accumulate;
	waitflush = (yres > 0) & (xres > 1) ? 0 : xres;
	if (!recover)
		return;
					/* recover previous values */
	if (accumulate <= 0)
		reload_output();
	else
		recover_output(stdin);
}

/* grow modifier to accommodate more bins */
MODCONT *
growmodifier(MODCONT *mp, int nb)
{
	if (nb <= mp->nbins)
		return mp;
	mp = (MODCONT *)realloc(mp, sizeof(MODCONT) + sizeof(DCOLOR)*(nb-1));
	if (mp == NULL)
		error(SYSTEM, "out of memory in growmodifier");
	memset(mp->cbin+mp->nbins, 0, sizeof(DCOLOR)*(nb-mp->nbins));
	mp->nbins = nb;
	return mp;
}

/* add modifier to our list to track */
MODCONT *
addmodifier(char *modn, char *outf, char *binv, int bincnt)
{
	LUENT	*lep = lu_find(&modconttab, modn);
	MODCONT	*mp;
	int	i;
	
	if (lep->data != NULL) {
		sprintf(errmsg, "duplicate modifier '%s'", modn);
		error(USER, errmsg);
	}
	if (nmods >= MAXMODLIST)
		error(INTERNAL, "too many modifiers");
	modname[nmods++] = modn;	/* XXX assumes static string */
	lep->key = modn;		/* XXX assumes static string */
	mp = (MODCONT *)malloc(sizeof(MODCONT));
	if (mp == NULL)
		error(SYSTEM, "out of memory in addmodifier");
	mp->outspec = outf;		/* XXX assumes static string */
	mp->modname = modn;		/* XXX assumes static string */
	if (binv == NULL)
		binv = "0";		/* use single bin if unspecified */
	mp->binv = eparse(binv);
	if (mp->binv->type == NUM) {	/* check value if constant */
		bincnt = (int)(evalue(mp->binv) + 1.5);
		if (bincnt != 1) {
			sprintf(errmsg, "illegal non-zero constant for bin (%s)",
					binv);
			error(USER, errmsg);
		}
	}
	mp->nbins = 1;			/* initialize results holder */
	setcolor(mp->cbin[0], 0., 0., 0.);
	if (bincnt > 1)
		mp = growmodifier(mp, bincnt);
	lep->data = (char *)mp;
					/* allocate output streams */
	for (i = bincnt; i-- > 0; )
		getostream(mp->outspec, mp->modname, i, 1);
	return mp;
}

/* add modifiers from a file list */
void
addmodfile(char *fname, char *outf, char *binv, int bincnt)
{
	char	*mname[MAXMODLIST];
	int	i;
					/* find the file & store strings */
	if (wordfile(mname, getpath(fname, getrlibpath(), R_OK)) < 0) {
		sprintf(errmsg, "cannot find modifier file '%s'", fname);
		error(SYSTEM, errmsg);
	}
	for (i = 0; mname[i]; i++)	/* add each one */
		addmodifier(mname[i], outf, binv, bincnt);
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

/* construct output file name and return flags whether modifier/bin present */
int
ofname(char *oname, const char *ospec, const char *mname, int bn)
{
	const char	*mnp = NULL;
	const char	*bnp = NULL;
	const char	*cp;
	
	if (ospec == NULL)
		return -1;
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
					return -1;
				mnp = cp;
				break;
			case 'd':
			case 'i':
			case 'o':
			case 'x':
			case 'X':
				if (bnp != NULL)
					return -1;
				bnp = cp;
				break;
			default:
				return -1;
			}
		}
	if (mnp != NULL) {			/* create file name */
		if (bnp != NULL) {
			if (bnp > mnp)
				sprintf(oname, ospec, mname, bn);
			else
				sprintf(oname, ospec, bn, mname);
			return OF_MODIFIER|OF_BIN;
		}
		sprintf(oname, ospec, mname);
		return OF_MODIFIER;
	}
	if (bnp != NULL) {
		sprintf(oname, ospec, bn);
		return OF_BIN;
	}
	strcpy(oname, ospec);
	return 0;
}

/* write header to the given output stream */
void
printheader(FILE *fout, const char *info)
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
	if (info != NULL)			/* add extra info if given */
		fputs(info, fout);
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
}

/* write resolution string to given output stream */
void
printresolu(FILE *fout, int xr, int yr)
{
	if ((xr > 0) & (yr > 0))	/* resolution string */
		fprtresolu(xr, yr, fout);
}

/* Get output stream pointer (open and write header if new and noopen==0) */
STREAMOUT *
getostream(const char *ospec, const char *mname, int bn, int noopen)
{
	static const DCOLOR	nocontrib = BLKCOLOR;
	static STREAMOUT	stdos;
	int			ofl;
	char			oname[1024];
	LUENT			*lep;
	STREAMOUT		*sop;
	
	if (ospec == NULL) {			/* use stdout? */
		if (!noopen && !using_stdout) {
			if (outfmt != 'a')
				SET_FILE_BINARY(stdout);
			if (header)
				printheader(stdout, NULL);
			printresolu(stdout, xres, yres);
			if (waitflush > 0)
				fflush(stdout);
			stdos.xr = xres; stdos.yr = yres;
			using_stdout = 1;
		}
		stdos.ofp = stdout;
		stdos.reclen += noopen;
		return &stdos;
	}
	ofl = ofname(oname, ospec, mname, bn);	/* get output name */
	if (ofl < 0) {
		sprintf(errmsg, "bad output format '%s'", ospec);
		error(USER, errmsg);
	}
	lep = lu_find(&ofiletab, oname);	/* look it up */
	if (lep->key == NULL)			/* new entry */
		lep->key = strcpy((char *)malloc(strlen(oname)+1), oname);
	sop = (STREAMOUT *)lep->data;
	if (sop == NULL) {			/* allocate stream */
		sop = (STREAMOUT *)malloc(sizeof(STREAMOUT));
		if (sop == NULL)
			error(SYSTEM, "out of memory in getostream");
		sop->outpipe = oname[0] == '!';
		sop->reclen = 0;
		sop->ofp = NULL;		/* open iff noopen==0 */
		sop->xr = xres; sop->yr = yres;
		lep->data = (char *)sop;
		if (!sop->outpipe & !force_open & !recover &&
				access(oname, F_OK) == 0) {
			errno = EEXIST;		/* file exists */
			goto openerr;
		}
	}
	if (!noopen && sop->ofp == NULL) {	/* open output stream */
		long		i;
		if (oname[0] == '!')		/* output to command */
			sop->ofp = popen(oname+1, "w");
		else				/* else open file */
			sop->ofp = fopen(oname, "w");
		if (sop->ofp == NULL)
			goto openerr;
		if (outfmt != 'a')
			SET_FILE_BINARY(sop->ofp);
		if (header) {
			char	info[512];
			char	*cp = info;
			if (ofl & OF_MODIFIER || sop->reclen == 1) {
				sprintf(cp, "MODIFIER=%s\n", mname);
				while (*cp) ++cp;
			}
			if (ofl & OF_BIN) {
				sprintf(cp, "BIN=%d\n", bn);
				while (*cp) ++cp;
			}
			*cp = '\0';
			printheader(sop->ofp, info);
		}
		if (accumulate > 0) {		/* global resolution */
			sop->xr = xres; sop->yr = yres;
		}
		printresolu(sop->ofp, sop->xr, sop->yr);
						/* play catch-up */
		for (i = accumulate > 0 ? lastdone/accumulate : 0; i--; ) {
			int	j = sop->reclen;
			if (j <= 0) j = 1;
			while (j--)
				put_contrib(nocontrib, sop->ofp);
			if (outfmt == 'a')
				putc('\n', sop->ofp);
		}
		if (waitflush > 0)
			fflush(sop->ofp);
	}
	sop->reclen += noopen;			/* add to length if noopen */
	return sop;				/* return output stream */
openerr:
	sprintf(errmsg, "cannot open '%s' for writing", oname);
	error(SYSTEM, errmsg);
	return NULL;	/* pro forma return */
}

/* read input ray into buffer */
int
getinp(char *buf, FILE *fp)
{
	double	dv[3], *dvp;
	float	*fvp;
	char	*cp;
	int	i;

	switch (inpfmt) {
	case 'a':
		cp = buf;		/* make sure we get 6 floats */
		for (i = 0; i < 6; i++) {
			if (fgetword(cp, buf+127-cp, fp) == NULL)
				return -1;
			if (i >= 3)
				dv[i-3] = atof(cp);
			if ((cp = fskip(cp)) == NULL || *cp)
				return -1;
			*cp++ = ' ';
		}
		getc(fp);		/* get/put eol */
		*cp-- = '\0'; *cp = '\n';
		if (DOT(dv,dv) <= FTINY*FTINY)
			return 0;	/* dummy ray */
		return strlen(buf);
	case 'f':
		if (fread(buf, sizeof(float), 6, fp) < 6)
			return -1;
		fvp = (float *)buf + 3;
		if (DOT(fvp,fvp) <= FTINY*FTINY)
			return 0;	/* dummy ray */
		return sizeof(float)*6;
	case 'd':
		if (fread(buf, sizeof(double), 6, fp) < 6)
			return -1;
		dvp = (double *)buf + 3;
		if (DOT(dvp,dvp) <= FTINY*FTINY)
			return 0;	/* dummy ray */
		return sizeof(double)*6;
	}
	error(INTERNAL, "botched input format");
	return -1;	/* pro forma return */
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
	else if (bn >= mp->nbins)	/* new bin */
		le->data = (char *)(mp = growmodifier(mp, bn+1));
	addcolor(mp->cbin[bn], rparams);
}

/* output newline to ASCII file and/or flush as requested */
static int
puteol(const LUENT *e, void *p)
{
	STREAMOUT	*sop = (STREAMOUT *)e->data;

	if (outfmt == 'a')
		putc('\n', sop->ofp);
	if (!waitflush)
		fflush(sop->ofp);
	if (ferror(sop->ofp)) {
		sprintf(errmsg, "write error on file '%s'", e->key);
		error(SYSTEM, errmsg);
	}
	return 0;
}

/* put out ray contribution to file */
void
put_contrib(const DCOLOR cnt, FILE *fout)
{
	float	fv[3];
	COLR	cv;

	switch (outfmt) {
	case 'a':
		fprintf(fout, "%.6e\t%.6e\t%.6e\t", cnt[0], cnt[1], cnt[2]);
		break;
	case 'f':
		copycolor(fv, cnt);
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
	if (waitflush < 0 && frandom() < 0.001)
		fflush(fout);			/* staggers writes */
}

/* output ray tallies and clear for next accumulation */
void
done_contrib(int navg)
{
	double		sf = 1.;
	int		i, j;
	MODCONT		*mp;
	STREAMOUT	*sop;
						/* set average scaling */
	if (navg > 1)
		sf = 1. / (double)navg;
						/* output modifiers in order */
	for (i = 0; i < nmods; i++) {
		mp = (MODCONT *)lu_find(&modconttab,modname[i])->data;
		if (navg > 1)			/* average scaling */
			for (j = mp->nbins; j--; )
				scalecolor(mp->cbin[j], sf);
		sop = getostream(mp->outspec, mp->modname, 0,0);
		put_contrib(mp->cbin[0], sop->ofp);
		if (mp->nbins > 3 &&		/* minor optimization */
				sop == getostream(mp->outspec, mp->modname, 1,0))
			for (j = 1; j < mp->nbins; j++)
				put_contrib(mp->cbin[j], sop->ofp);
		else
			for (j = 1; j < mp->nbins; j++)
				put_contrib(mp->cbin[j],
				    getostream(mp->outspec,mp->modname,j,0)->ofp);
						/* clear for next time */
		memset(mp->cbin, 0, sizeof(DCOLOR)*mp->nbins);
	}
	--waitflush;				/* terminate records */
	lu_doall(&ofiletab, puteol, NULL);
	if (using_stdout & (outfmt == 'a'))
		putc('\n', stdout);
	if (!waitflush) {
		waitflush = (yres > 0) & (xres > 1) ? 0 : xres;
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
			while (n > 1 && *cp != '\t') {
				if (mnp - modname >= sizeof(modname)-2)
					error(INTERNAL, "modifier name too long");
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
					/* time to produce record? */
		if (account > 0 && !--account)
			done_contrib(account = accumulate);
		lastdone = rtp->raynum;
		if (rtp->buf != NULL)	/* free up buffer space */
			free(rtp->buf);
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
	ssize_t		nr;
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
					treebufsiz = rt->bsiz += BUFSIZ;
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
	static int	ignore_warning_given = 0;
	char		inpbuf[128];
	int		iblen;
	struct rtproc	*rtp;
						/* loop over input */
	while ((iblen = getinp(inpbuf, fin)) >= 0) {
		if (!iblen && accumulate != 1) {
			if (!ignore_warning_given++)
				error(WARNING,
				"dummy ray(s) ignored during accumulation\n");
			continue;
		}
		if (!iblen ||			/* need flush/reset? */
				queue_length() > 10*nrtprocs() ||
				lastray+1 < lastray) {
			while (wait_rproc() != NULL)
				process_queue();
			if (lastray+1 < lastray)
				lastdone = lastray = 0;
		}
		rtp = get_rproc();		/* get avail. rtrace process */
		rtp->raynum = ++lastray;	/* assign ray */
		if (iblen) {			/* trace ray if valid */
			writebuf(rtp->pd.w, inpbuf, iblen);
		} else {			/* else bypass dummy ray */
			queue_raytree(rtp);	/* queue empty ray/record */
			if ((yres <= 0) | (xres <= 0))
				waitflush = 1;	/* flush right after */
		}
		process_queue();		/* catch up with results */
		if (raysleft && !--raysleft)
			break;			/* preemptive EOI */
	}
	while (wait_rproc() != NULL)		/* process outstanding rays */
		process_queue();
	if (accumulate <= 0)
		done_contrib(0);		/* output tallies */
	else if (account < accumulate) {
		error(WARNING, "partial accumulation in final record");
		done_contrib(accumulate - account);
	}
	if (raysleft)
		error(USER, "unexpected EOF on input");
	lu_done(&ofiletab);			/* close output files */
}

/* get ray contribution from previous file */
static int
get_contrib(DCOLOR cnt, FILE *finp)
{
	COLOR	fv;
	COLR	cv;

	switch (outfmt) {
	case 'a':
		return(fscanf(finp,"%lf %lf %lf",&cnt[0],&cnt[1],&cnt[2]) == 3);
	case 'f':
		if (fread(fv, sizeof(fv[0]), 3, finp) != 3)
			return(0);
		copycolor(cnt, fv);
		return(1);
	case 'd':
		return(fread(cnt, sizeof(cnt[0]), 3, finp) == 3);
	case 'c':
		if (fread(cv, sizeof(cv), 1, finp) != 1)
			return(0);
		colr_color(fv, cv);
		copycolor(cnt, fv);
		return(1);
	default:
		error(INTERNAL, "botched output format");
	}
	return(0);	/* pro forma return */
}

/* close output file opened for input */
static int
myclose(const LUENT *e, void *p)
{
	STREAMOUT	*sop = (STREAMOUT *)e->data;
	
	if (sop->ofp == NULL)
		return(0);
	fclose(sop->ofp);
	sop->ofp = NULL;
	return(0);
}

/* load previously accumulated values */
void
reload_output(void)
{
	int		i, j;
	MODCONT		*mp;
	int		ofl;
	char		oname[1024];
	char		*fmode = "rb";
	char		*outvfmt;
	LUENT		*ment, *oent;
	int		xr, yr;
	STREAMOUT	sout;
	DCOLOR		rgbv;

	switch (outfmt) {
	case 'a':
		outvfmt = "ascii";
		fmode = "r";
		break;
	case 'f':
		outvfmt = "float";
		break;
	case 'd':
		outvfmt = "double";
		break;
	case 'c':
		outvfmt = COLRFMT;
		break;
	default:
		error(INTERNAL, "botched output format");
		return;
	}
						/* reload modifier values */
	for (i = 0; i < nmods; i++) {
		ment = lu_find(&modconttab,modname[i]);
		mp = (MODCONT *)ment->data;
		if (mp->outspec == NULL)
			error(USER, "cannot reload from stdout");
		if (mp->outspec[0] == '!')
			error(USER, "cannot reload from command");
		for (j = 0; ; j++) {		/* load each modifier bin */
			ofl = ofname(oname, mp->outspec, mp->modname, j);
			if (ofl < 0)
				error(USER, "bad output file specification");
			oent = lu_find(&ofiletab, oname);
			if (oent->data != NULL) {
				sout = *(STREAMOUT *)oent->data;
			} else {
				sout.reclen = 0;
				sout.outpipe = 0;
				sout.xr = xres; sout.yr = yres;
				sout.ofp = NULL;
			}
			if (sout.ofp == NULL) {	/* open output as input */
				sout.ofp = fopen(oname, fmode);
				if (sout.ofp == NULL) {
					if (j)
						break;	/* assume end of modifier */
					sprintf(errmsg, "missing reload file '%s'",
							oname);
					error(WARNING, errmsg);
					break;
				}
				if (header && checkheader(sout.ofp, outvfmt, NULL) != 1) {
					sprintf(errmsg, "format mismatch for '%s'",
							oname);
					error(USER, errmsg);
				}
				if ((sout.xr > 0) & (sout.yr > 0) &&
						(!fscnresolu(&xr, &yr, sout.ofp) ||
							(xr != sout.xr) |
							(yr != sout.yr))) {
					sprintf(errmsg, "resolution mismatch for '%s'",
							oname);
					error(USER, errmsg);
				}
			}
							/* read in RGB value */
			if (!get_contrib(rgbv, sout.ofp)) {
				if (!j) {
					fclose(sout.ofp);
					break;		/* ignore empty file */
				}
				if (j < mp->nbins) {
					sprintf(errmsg, "missing data in '%s'",
							oname);
					error(USER, errmsg);
				}
				break;
			}
			if (j >= mp->nbins)		/* grow modifier size */
				ment->data = (char *)(mp = growmodifier(mp, j+1));
			copycolor(mp->cbin[j], rgbv);
			if (oent->key == NULL)		/* new file entry */
				oent->key = strcpy((char *)
						malloc(strlen(oname)+1), oname);
			if (oent->data == NULL)
				oent->data = (char *)malloc(sizeof(STREAMOUT));
			*(STREAMOUT *)oent->data = sout;
		}
	}
	lu_doall(&ofiletab, myclose, NULL);	/* close all files */
}

/* seek on the given output file */
static int
myseeko(const LUENT *e, void *p)
{
	STREAMOUT	*sop = (STREAMOUT *)e->data;
	off_t		nbytes = *(off_t *)p;
	
	if (sop->reclen > 1)
		nbytes = nbytes * sop->reclen;
	if (fseeko(sop->ofp, nbytes, SEEK_CUR) < 0) {
		sprintf(errmsg, "seek error on file '%s'", e->key);
		error(SYSTEM, errmsg);
	}
	return 0;
}

/* recover output if possible */
void
recover_output(FILE *fin)
{
	off_t		lastout = -1;
	int		outvsiz, recsiz;
	char		*outvfmt;
	int		i, j;
	MODCONT		*mp;
	int		ofl;
	char		oname[1024];
	LUENT		*ment, *oent;
	STREAMOUT	sout;
	off_t		nvals;
	int		xr, yr;

	switch (outfmt) {
	case 'a':
		error(USER, "cannot recover ASCII output");
		return;
	case 'f':
		outvsiz = sizeof(float)*3;
		outvfmt = "float";
		break;
	case 'd':
		outvsiz = sizeof(double)*3;
		outvfmt = "double";
		break;
	case 'c':
		outvsiz = sizeof(COLR);
		outvfmt = COLRFMT;
		break;
	default:
		error(INTERNAL, "botched output format");
		return;
	}
						/* check modifier outputs */
	for (i = 0; i < nmods; i++) {
		ment = lu_find(&modconttab,modname[i]);
		mp = (MODCONT *)ment->data;
		if (mp->outspec == NULL)
			error(USER, "cannot recover from stdout");
		if (mp->outspec[0] == '!')
			error(USER, "cannot recover from command");
		for (j = 0; ; j++) {		/* check each bin's file */
			ofl = ofname(oname, mp->outspec, mp->modname, j);
			if (ofl < 0)
				error(USER, "bad output file specification");
			oent = lu_find(&ofiletab, oname);
			if (oent->data != NULL) {
				sout = *(STREAMOUT *)oent->data;
			} else {
				sout.reclen = 0;
				sout.outpipe = 0;
				sout.ofp = NULL;
			}
			if (sout.ofp != NULL) {	/* already open? */
				if (ofl & OF_BIN)
					continue;
				break;
			}
						/* open output */
			sout.ofp = fopen(oname, "rb+");
			if (sout.ofp == NULL) {
				if (j)
					break;	/* assume end of modifier */
				sprintf(errmsg, "missing recover file '%s'",
						oname);
				error(WARNING, errmsg);
				break;
			}
			nvals = lseek(fileno(sout.ofp), 0, SEEK_END);
			if (nvals <= 0) {
				lastout = 0;	/* empty output, quit here */
				fclose(sout.ofp);
				break;
			}
			if (!sout.reclen) {
				if (!(ofl & OF_BIN)) {
					sprintf(errmsg,
						"need -bn to recover file '%s'",
							oname);
					error(USER, errmsg);
				}
				recsiz = outvsiz;
			} else
				recsiz = outvsiz * sout.reclen;

			lseek(fileno(sout.ofp), 0, SEEK_SET);
			if (header && checkheader(sout.ofp, outvfmt, NULL) != 1) {
				sprintf(errmsg, "format mismatch for '%s'",
						oname);
				error(USER, errmsg);
			}
			sout.xr = xres; sout.yr = yres;
			if ((sout.xr > 0) & (sout.yr > 0) &&
					(!fscnresolu(&xr, &yr, sout.ofp) ||
						(xr != sout.xr) |
						(yr != sout.yr))) {
				sprintf(errmsg, "resolution mismatch for '%s'",
						oname);
				error(USER, errmsg);
			}
			nvals = (nvals - (off_t)ftell(sout.ofp)) / recsiz;
			if ((lastout < 0) | (nvals < lastout))
				lastout = nvals;
			if (oent->key == NULL)	/* new entry */
				oent->key = strcpy((char *)
						malloc(strlen(oname)+1), oname);
			if (oent->data == NULL)
				oent->data = (char *)malloc(sizeof(STREAMOUT));
			*(STREAMOUT *)oent->data = sout;
			if (!(ofl & OF_BIN))
				break;		/* no bin separation */
		}
		if (!lastout) {			/* empty output */
			error(WARNING, "no previous data to recover");
			lu_done(&ofiletab);	/* reclose all outputs */
			return;
		}
		if (j > mp->nbins)		/* reallocate modifier bins */
			ment->data = (char *)(mp = growmodifier(mp, j));
	}
	if (lastout < 0) {
		error(WARNING, "no output files to recover");
		return;
	}
	if (raysleft && lastout >= raysleft/accumulate) {
		error(WARNING, "output appears to be complete");
		/* XXX should read & discard input? */
		quit(0);
	}
						/* seek on all files */
	nvals = lastout * outvsiz;
	lu_doall(&ofiletab, myseeko, &nvals);
						/* skip repeated input */
	for (nvals = 0; nvals < lastout; nvals++)
		if (getinp(oname, fin) < 0)
			error(USER, "unexpected EOF on input");
	lastray = lastdone = (unsigned long)lastout * accumulate;
	if (raysleft)
		raysleft -= lastray;
}
