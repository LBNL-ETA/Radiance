#ifndef lint
static const char	RCSid[] = "$Id: rxcmain.cpp,v 2.13 2024/12/24 16:58:13 greg Exp $";
#endif
/*
 *  rxcmain.c - main for rxcontrib ray contribution tracer
 */

#include "copyright.h"

#include <signal.h>
#include <time.h>
#include <ctype.h>
#include "RcontribSimulManager.h"
#include "platform.h"
#include "func.h"

const char	*sigerr[NSIG];		/* signal error messages */

int	nproc = 1;			/* number of processes requested */

int	inpfmt = 'a';			/* input format */
int	outfmt = 'f';			/* output format */

int	report_intvl = 0;		/* reporting interval (seconds) */

extern char *	progname;		// global argv[0]

RcontribSimulManager	myRCmanager;	// global rcontrib simulation manager

#define RCONTRIB_FEATURES	"Multiprocessing\n" \
				"Accumulation\nRecovery\n" \
				"ImmediateIrradiance\n" \
				"ProgressReporting\nDistanceLimiting\n" \
				"InputFormats=a,f,d\nOutputFormats=f,d,c\n" \
				"Outputs=V,W\n" \
				"OutputCS=RGB,spec\n"

static void	rxcontrib(const int rstart = 0);

static void
printdefaults(void)			/* print default values to stdout */
{
	printf("-c %-5d\t\t\t# accumulated rays per record\n", myRCmanager.accum);
	printf(myRCmanager.HasFlag(RCcontrib) ?
			"-V+\t\t\t\t# output contributions\n" :
			"-V-\t\t\t\t# output coefficients\n");
	if (myRCmanager.HasFlag(RTimmIrrad))
		printf("-I+\t\t\t\t# immediate irradiance on\n");
	printf("-n %-2d\t\t\t\t# number of rendering processes\n", nproc);
	if (myRCmanager.xres > 0)
		printf("-x %-9d\t\t\t# x resolution\n", myRCmanager.xres);
	printf("-y %-9d\t\t\t# y resolution\n", myRCmanager.yres);
	printf(myRCmanager.HasFlag(RTlimDist) ?
			"-ld+\t\t\t\t# limit distance on\n" :
			"-ld-\t\t\t\t# limit distance off\n");
	printf("-f%c%c\t\t\t\t# format input/output = %s/%s\n",
			inpfmt, outfmt, formstr(inpfmt), formstr(outfmt));
	if (report_intvl > 0)
		printf("-t %-9d\t\t\t# time between reports\n", report_intvl);
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
	alarm(180);			/* allow 3 minutes to clean up */
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
	const char  *msg
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
		if (inpfmt == 'a')
			goto fmterr;
		outfmt = inpfmt;
		return;
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
	sprintf(errmsg, "Unsupported i/o format: -f%s", fmt);
	error(USER, errmsg);
}

/* Set default options */
static void
default_options(void)
{
	rand_samp = 1;
	dstrsrc = 0.9;
	directrelay = 3;
	vspretest = 512;
	srcsizerat = .2;
	specthresh = .02;
	specjitter = 1.;
	maxdepth = -10;
	minweight = 2e-3;
	ambres = 256;
	ambdiv = 350;
	ambounce = 1;
}

/* Set overriding options */
static void
override_options(void)
{
	shadthresh = 0;
	ambssamp = 0;
	ambacc = 0;
}

int
main(int argc, char *argv[])
{
#define	 check(ol,al)		if (argv[i][ol] || \
				badarg(argc-i-1,argv+i+1,al)) \
				goto badopt
#define	 check_bool(olen,var)		switch (argv[i][olen]) { \
				case '\0': var = !var; break; \
				case 'y': case 'Y': case 't': case 'T': \
				case '+': case '1': var = 1; break; \
				case 'n': case 'N': case 'f': case 'F': \
				case '-': case '0': var = 0; break; \
				default: goto badopt; }
	int	force_open = 0;
	int	recover = 0;
	char	*curout = NULL;
	char	*prms = NULL;
	char	*binval = NULL;
	int	bincnt = 0;
	int	rval;
	int	i;
					/* global program name */
	progname = argv[0];
					/* feature check only? */
	strcat(RFeatureList, RCONTRIB_FEATURES);
	if (argc > 1 && !strcmp(argv[1], "-features"))
		return feature_status(argc-2, argv+2);
					/* initialize calcomp routines early */
	initfunc();
	calcontext(RCCONTEXT);
					/* set rcontrib defaults */
	default_options();
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
		case 'n':			/* number of processes */
			check(2,"i");
			nproc = atoi(argv[++i]);
			if (nproc < 0 && (nproc += RadSimulManager::GetNCores()) <= 0)
				nproc = 1;
			break;
		case 'V':			/* output contributions? */
			rval = myRCmanager.HasFlag(RCcontrib);
			check_bool(2,rval);
			myRCmanager.SetFlag(RCcontrib, rval);
			break;
		case 'x':			/* x resolution */
			check(2,"i");
			myRCmanager.xres = atoi(argv[++i]);
			break;
		case 'y':			/* y resolution */
			check(2,"i");
			myRCmanager.yres = atoi(argv[++i]);
			break;
		case 'w':			/* warnings on/off */
			rval = (erract[WARNING].pf != NULL);
			check_bool(2,rval);
			if (rval) erract[WARNING].pf = wputs;
			else erract[WARNING].pf = NULL;
			break;
		case 'e':			/* .cal expression */
			check(2,"s");
			scompile(argv[++i], NULL, 0);
			break;
		case 'l':			/* limit distance */
			if (argv[i][2] != 'd')
				goto badopt;
			rval = myRCmanager.HasFlag(RTlimDist);
			check_bool(3,rval);
			myRCmanager.SetFlag(RTlimDist, rval);
			break;
		case 'I':			/* immed. irradiance */
			rval = myRCmanager.HasFlag(RTimmIrrad);
			check_bool(2,rval);
			myRCmanager.SetFlag(RTimmIrrad, rval);
			break;
		case 'f':			/* .cal file or force or format */
			if (!argv[i][2]) {
				check(2,"s");
				loadfunc(argv[++i]);
				break;
			}
			if (argv[i][2] == 'o') {
				check_bool(3,force_open);
				break;
			}
			setformat(argv[i]+2);
			myRCmanager.SetDataFormat(outfmt);
			break;
		case 'o':			/* output file */
			check(2,"s");
			curout = argv[++i];
			break;
		case 'r':			/* recover output */
			check_bool(2,recover);
			break;
		case 'p':			/* parameter setting(s) */
			check(2,"s");
			set_eparams(prms = argv[++i]);
			break;
		case 'c':			/* sample count */
			check(2,"i");
			myRCmanager.accum = atoi(argv[++i]);
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
			myRCmanager.AddModifier(argv[++i], curout, prms, binval, bincnt);
			break;
		case 'M':			/* file of modifier names */
			check(2,"s");
			myRCmanager.AddModFile(argv[++i], curout, prms, binval, bincnt);
			break;
		case 't':			/* reporting interval */
			check(2,"i");
			report_intvl = atoi(argv[++i]);
			break;
		default:
			goto badopt;
		}
	}
	if (i != argc-1)
		error(USER, "expected single octree argument");

	override_options();		/* override some option settings */

	if (!myRCmanager.GetOutput())	// check that we have work to do
		error(USER, "missing required modifier argument");
					// get ready to rock...
	if (setspectrsamp(CNDX, WLPART) < 0)
		error(USER, "unsupported spectral sampling");
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
					// load octree
	myRCmanager.LoadOctree(argv[argc-1]);
					// add to header
	myRCmanager.AddHeader(argc-1, argv);
					// prepare output files
	if (recover)
		myRCmanager.outOp = RCOrecover;
	else if (force_open)
		myRCmanager.outOp = RCOforce;
	else
		myRCmanager.outOp = RCOnew;
					// rval = # rows recovered
	rval = myRCmanager.PrepOutput();
					// check if recovered everything
	if (rval >= myRCmanager.GetRowMax()) {
		error(WARNING, "nothing left to compute");
		quit(0);
	}
	rxcontrib(rval);		/* trace ray contributions (loop) */

	quit(0);	/* exit clean */

badopt:
	fprintf(stderr,
"Usage: %s [-V][-c count][-r][-e expr][-f source][-o ospec][-p p1=V1,p2=V2][-b binv][-bn N] {-m mod | -M file} [rtrace options] octree\n",
			progname);
	sprintf(errmsg, "command line error at '%s'", argv[i]);
	error(USER, errmsg);
	return(1);	/* pro forma return */

#undef	check
#undef	check_bool
}


// skip specified number of bytes, return false if EOF
static bool
skipBytes(int n2skip)
{
	while (n2skip-- > 0)
		if (getchar() == EOF)
			return false;
	return true;
}


// skip specified number of whitespace-separated words, return false if EOF
static bool
skipWords(int n2skip)
{
	int	c;

	while (n2skip-- > 0) {
		do {
			c = getchar();
		} while (isspace(c));
		do {
			if (c == EOF) return false;
			c = getchar();
		} while (!isspace(c));
	}
	return true;
}


// read a bundle of myRCmanager.accum ray origins and directions
static bool
getRayBundle(FVECT *orig_dir = NULL)
{
	int	n2go = myRCmanager.accum;

	switch (inpfmt) {
	case 'a':			// ASCII input
		if (!orig_dir)
			return skipWords(6*n2go);
		while (n2go-- > 0) {
			if (scanf(FVFORMAT, &orig_dir[0][0],
					&orig_dir[0][1], &orig_dir[0][2]) != 3)
				return false;
			if (scanf(FVFORMAT, &orig_dir[1][0],
					&orig_dir[1][1], &orig_dir[1][2]) != 3)
				return false;
			orig_dir += 2;
		}
		break;
	case 'f':			// float input
		if (!orig_dir)
			return skipBytes(6*sizeof(float)*n2go);
#ifdef SMLFLT
		if (getbinary(orig_dir, sizeof(FVECT), 2*n2go, stdin) != 2*n2go)
			return false;
		orig_dir += 2*n2go;
#else
		while (n2go-- > 0) {
			float	fvecs[6];
			if (getbinary(fvecs, sizeof(fvecs), 1, stdin) != 1)
				return false;
			for (int i = 6; i--; )
				orig_dir[0][i] = fvecs[i];
			orig_dir += 2;
		}
#endif
		break;
	case 'd':			// double input
		if (!orig_dir)
			return skipBytes(6*sizeof(double)*n2go);
#ifndef SMLFLT
		if (getbinary(orig_dir, sizeof(FVECT), 2*n2go, stdin) != 2*n2go)
			return false;
		orig_dir += 2*n2go;
#else
		while (n2go-- > 0) {
			double	dvecs[6];
			if (getbinary(dvecs, sizeof(dvecs), 1, stdin) != 1)
				return false;
			for (int i = 6; i--; )
				orig_dir[0][i] = dvecs[i];
			orig_dir += 2;
		}
#endif
		break;
	default:
		error(INTERNAL, "unsupported format in getRayBundle()");
		return false;
	}
	n2go = myRCmanager.accum;	// normalize directions
	while (n2go-- > 0) {
		orig_dir -= 2;
		normalize(orig_dir[1]);
	}
	return true;
}


// Run loop to load data, report progress, etc.
void
rxcontrib(const int rstart)
{
	const int	totRows = myRCmanager.GetRowMax();
	FVECT *		odarr = (FVECT *)emalloc(sizeof(FVECT)*2*myRCmanager.accum);
	time_t		tstart, last_report;
	int		r = 0;

	while (r < rstart) {		// skip input rays already done
		if (!getRayBundle())
			goto readerr;
		r++;
	}
	if (report_intvl > 0) {		// set up reporting
		if (r > 0) {
			sprintf(errmsg, "recovered %.2f%% of total\n",
					100.*r/totRows);
			eputs(errmsg);
		}
		last_report = tstart = time(0);
	}
					// start children as requested
	myRCmanager.SetThreadCount(nproc);

	while (r < totRows) {		// loop until done
		time_t	tnow;
		if (!getRayBundle(odarr))
			goto readerr;
		if (myRCmanager.ComputeRecord(odarr) <= 0)
			return;		// error reported, hopefully...
		r++;
		if (report_intvl <= 0)
			continue;
		if (r == totRows)	// need to finish up?
			myRCmanager.SetThreadCount(1);
		tnow = time(0);
		if ((r < totRows) & (tnow < last_report+report_intvl))
			continue;
		sprintf(errmsg, "%.2f%% done after %.3f hours\n",
				100.*myRCmanager.GetRowFinished()/totRows,
				(1./3600.)*(tnow - tstart));
		eputs(errmsg);
		last_report = tnow;
	}
	efree(odarr);
	return;
readerr:
	sprintf(errmsg, "unexpected EOF on standard input (record %d of %d)",
			r, totRows);
	error(USER, errmsg);
}


void
wputs(				/* warning output function */
	const char	*s
)
{
	if (!erract[WARNING].pf) return;
	int  lasterrno = errno;
	eputs(s);
	errno = lasterrno;
}


void
eputs(				/* put string to stderr */
	const char  *s
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


/* Exit program */
void
quit(
	int  code
)
{
	if (!code && myRCmanager.Ready())	// clean up on normal exit
		code = myRCmanager.Cleanup();

	exit(code);
}
