#ifndef lint
static const char	RCSid[] = "$Id: rcmain.c,v 2.8 2012/06/27 15:32:58 greg Exp $";
#endif
/*
 *  rcmain.c - main for rtcontrib ray contribution tracer
 */

#include "copyright.h"

#include <signal.h>
#include "rcontrib.h"
#include "random.h"
#include "source.h"
#include "ambient.h"

int	gargc;				/* global argc */
char	**gargv;			/* global argv */
char	*octname;			/* global octree name */
char	*progname;			/* global argv[0] */

char	*sigerr[NSIG];			/* signal error messages */

int	nproc = 1;			/* number of processes requested */
int	nchild = 0;			/* number of children (-1 in child) */

int	inpfmt = 'a';			/* input format */
int	outfmt = 'a';			/* output format */

int	header = 1;			/* output header? */
int	force_open = 0;			/* truncate existing output? */
int	recover = 0;			/* recover previous output? */
int	accumulate = 1;			/* input rays per output record */
int	contrib = 0;			/* computing contributions? */

int	xres = 0;			/* horizontal (scan) size */
int	yres = 0;			/* vertical resolution */

int	using_stdout = 0;		/* are we using stdout? */

int	imm_irrad = 0;			/* compute immediate irradiance? */
int	lim_dist = 0;			/* limit distance? */

const char	*modname[MAXMODLIST];	/* ordered modifier name list */
int		nmods = 0;		/* number of modifiers */

void	(*addobjnotify[8])() = {ambnotify, NULL};

char	RCCONTEXT[] = "RCONTRIB";	/* our special evaluation context */


static void
printdefaults(void)			/* print default values to stdout */
{
	char  *cp;

	printf("-c %-5d\t\t\t# accumulated rays per record\n", accumulate);
	printf("-V%c\t\t\t\t# output %s\n", contrib ? '+' : '-',
			contrib ? "contributions" : "coefficients");
	if (imm_irrad)
		printf("-I+\t\t\t\t# immediate irradiance on\n");
	printf("-n %-2d\t\t\t\t# number of rendering processes\n", nproc);
	printf("-x %-9d\t\t\t# %s\n", xres,
			yres && xres ? "x resolution" : "flush interval");
	printf("-y %-9d\t\t\t# y resolution\n", yres);
	printf(lim_dist ? "-ld+\t\t\t\t# limit distance on\n" :
			"-ld-\t\t\t\t# limit distance off\n");
	printf("-h%c\t\t\t\t# %s header\n", header ? '+' : '-',
			header ? "output" : "no");
	printf("-f%c%c\t\t\t\t# format input/output = %s/%s\n",
			inpfmt, outfmt, formstr(inpfmt), formstr(outfmt));
	printf(erract[WARNING].pf != NULL ?
			"-w+\t\t\t\t# warning messages on\n" :
			"-w-\t\t\t\t# warning messages off\n");
	print_rdefaults();
}


static void
onsig(				/* fatal signal */
	int  signo
)
{
	static int  gotsig = 0;

	if (gotsig++)			/* two signals and we're gone! */
		_exit(signo);

#ifdef SIGALRM
	alarm(15);			/* allow 15 seconds to clean up */
	signal(SIGALRM, SIG_DFL);	/* make certain we do die */
#endif
	eputs("signal - ");
	eputs(sigerr[signo]);
	eputs("\n");
	quit(3);
}


static void
sigdie(			/* set fatal signal */
	int  signo,
	char  *msg
)
{
	if (signal(signo, onsig) == SIG_IGN)
		signal(signo, SIG_IGN);
	sigerr[signo] = msg;
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


/* Set overriding options */
static void
override_options(void)
{
	shadthresh = 0;
	ambssamp = 0;
	ambacc = 0;
	if (accumulate <= 0)	/* no output flushing for single record */
		xres = yres = 0;
}


int
main(int argc, char *argv[])
{
#define	 check(ol,al)		if (argv[i][ol] || \
				badarg(argc-i-1,argv+i+1,al)) \
				goto badopt
#define	 bool(olen,var)		switch (argv[i][olen]) { \
				case '\0': var = !var; break; \
				case 'y': case 'Y': case 't': case 'T': \
				case '+': case '1': var = 1; break; \
				case 'n': case 'N': case 'f': case 'F': \
				case '-': case '0': var = 0; break; \
				default: goto badopt; }
	char	*curout = NULL;
	char	*binval = NULL;
	int	bincnt = 0;
	int	rval;
	int	i;
					/* global program name */
	progname = argv[0] = fixargv0(argv[0]);
	gargv = argv;
	gargc = argc;
					/* initialize calcomp routines early */
	initfunc();
	setcontext(RCCONTEXT);
					/* option city */
	for (i = 1; i < argc; i++) {
						/* expand arguments */
		while ((rval = expandarg(&argc, &argv, i)) > 0)
			;
		if (rval < 0) {
			sprintf(errmsg, "cannot expand '%s'", argv[i]);
			error(SYSTEM, errmsg);
		}
		if (argv[i] == NULL || argv[i][0] != '-')
			break;			/* break from options */
		if (!strcmp(argv[i], "-version")) {
			puts(VersionID);
			quit(0);
		}
		if (!strcmp(argv[i], "-defaults") ||
				!strcmp(argv[i], "-help")) {
			override_options();
			printdefaults();
			quit(0);
		}
		rval = getrenderopt(argc-i, argv+i);
		if (rval >= 0) {
			i += rval;
			continue;
		}
		switch (argv[i][1]) {
		case 'n':			/* number of cores */
			check(2,"i");
			nproc = atoi(argv[++i]);
			if (nproc <= 0)
				error(USER, "bad number of processes");
			break;
		case 'V':			/* output contributions */
			bool(2,contrib);
			break;
		case 'x':			/* x resolution */
			check(2,"i");
			xres = atoi(argv[++i]);
			break;
		case 'y':			/* y resolution */
			check(2,"i");
			yres = atoi(argv[++i]);
			break;
		case 'w':			/* warnings */
			rval = (erract[WARNING].pf != NULL);
			bool(2,rval);
			if (rval) erract[WARNING].pf = wputs;
			else erract[WARNING].pf = NULL;
			break;
		case 'e':			/* expression */
			check(2,"s");
			scompile(argv[++i], NULL, 0);
			break;
		case 'l':			/* limit distance */
			if (argv[i][2] != 'd')
				goto badopt;
			bool(3,lim_dist);
			break;
		case 'I':			/* immed. irradiance */
			bool(2,imm_irrad);
			break;
		case 'f':			/* file or force or format */
			if (!argv[i][2]) {
				check(2,"s");
				loadfunc(argv[++i]);
				break;
			}
			if (argv[i][2] == 'o') {
				bool(3,force_open);
				break;
			}
			setformat(argv[i]+2);
			break;
		case 'o':			/* output */
			check(2,"s");
			curout = argv[++i];
			break;
		case 'c':			/* input rays per output */
			check(2,"i");
			accumulate = atoi(argv[++i]);
			break;
		case 'r':			/* recover output */
			bool(2,recover);
			break;
		case 'h':			/* header output */
			bool(2,header);
			break;
		case 'b':			/* bin expression/count */
			if (argv[i][2] == 'n') {
				check(3,"s");
				bincnt = (int)(eval(argv[++i]) + .5);
				break;
			}
			check(2,"s");
			binval = argv[++i];
			break;
		case 'm':			/* modifier name */
			check(2,"s");
			addmodifier(argv[++i], curout, binval, bincnt);
			break;
		case 'M':			/* modifier file */
			check(2,"s");
			addmodfile(argv[++i], curout, binval, bincnt);
			break;
		default:
			goto badopt;
		}
	}
					/* override some option settings */
	override_options();
					/* initialize object types */
	initotypes();
					/* initialize urand */
	if (rand_samp) {
		srandom((long)time(0));
		initurand(0);
	} else {
		srandom(0L);
		initurand(2048);
	}
					/* set up signal handling */
	sigdie(SIGINT, "Interrupt");
#ifdef SIGHUP
	sigdie(SIGHUP, "Hangup");
#endif
	sigdie(SIGTERM, "Terminate");
#ifdef SIGPIPE
	sigdie(SIGPIPE, "Broken pipe");
#endif
#ifdef SIGALRM
	sigdie(SIGALRM, "Alarm clock");
#endif
#ifdef	SIGXCPU
	sigdie(SIGXCPU, "CPU limit exceeded");
	sigdie(SIGXFSZ, "File size exceeded");
#endif
#ifdef	NICE
	nice(NICE);			/* lower priority */
#endif
					/* get octree */
	if (i == argc)
		octname = NULL;
	else if (i == argc-1)
		octname = argv[i];
	else
		goto badopt;
	if (octname == NULL)
		error(USER, "missing octree argument");

	readoct(octname, ~(IO_FILES|IO_INFO), &thescene, NULL);
	nsceneobjs = nobjects;

	marksources();			/* find and mark sources */

	setambient();			/* initialize ambient calculation */

	rcontrib();			/* trace ray contributions (loop) */

	ambsync();			/* flush ambient file */

	quit(0);	/* exit clean */

badopt:
	fprintf(stderr,
"Usage: %s [-n nprocs][-V][-r][-e expr][-f source][-o ospec][-b binv][-bn N] {-m mod | -M file} [rtrace options] octree\n",
			progname);
	sprintf(errmsg, "command line error at '%s'", argv[i]);
	error(USER, errmsg);
	return(1);	/* pro forma return */

#undef	check
#undef	bool
}


void
wputs(				/* warning output function */
	char	*s
)
{
	int  lasterrno = errno;
	eputs(s);
	errno = lasterrno;
}


void
eputs(				/* put string to stderr */
	char  *s
)
{
	static int  midline = 0;

	if (!*s)
		return;
	if (!midline++) {
		fputs(progname, stderr);
		fputs(": ", stderr);
	}
	fputs(s, stderr);
	if (s[strlen(s)-1] == '\n') {
		fflush(stderr);
		midline = 0;
	}
}
