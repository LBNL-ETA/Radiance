#ifndef lint
static const char	RCSid[] = "$Id: rxpmain.cpp,v 2.7 2025/01/02 20:19:48 greg Exp $";
#endif
/*
 *  rxpmain.cpp - main for rxpict batch rendering program
 */

#include "copyright.h"

#include  <time.h>
#include  <signal.h>

#include  "rtprocess.h" /* getpid() */
#include  "platform.h"
#include  "RpictSimulManager.h"

extern char  *progname;			/* argv[0] */
const char  *sigerr[NSIG];		/* signal error messages */
char  *errfile = NULL;			/* error output file */

VIEW  ourview = STDVIEW;		/* view parameters */
int  hresolu = 512;			/* horizontal resolution */
int  vresolu = 512;			/* vertical resolution */
double	pixaspect = 1.0;		/* pixel aspect ratio */
int  hres, vres;			/* current image resolution for srcdraw.c */

int  psample = 4;			/* pixel sample size */
double	maxdiff = .05;			/* max. difference for interpolation */
double	dstrpix = 0.67;			/* square pixel distribution */

double  mblur = 0.;			/* motion blur parameter */

double  dblur = 0.;			/* depth-of-field blur parameter */

int  nproc = 1;				/* number of processes to run */

int  ralrm = 0;				/* seconds between reports */

double	pctdone = 0.0;			/* percentage done */
time_t  tlastrept = 0L;			/* time at last report */
time_t  tstart;				/* starting time */

RenderDataType	dtype = RDTrgbe;	// output data flags

RpictSimulManager	myRPmanager;	// global simulation manager

static void onsig(int signo);
static void sigdie(int  signo, const char  *msg);
static void printdefaults(void);

					/* rxpict additional features */
#define RXPICT_FEATURES	"Recovery\nIrradianceCalc\nViewTypes=v,l,a,h,s,c\n" \
		"ParticipatingMedia=Mist\n" \
		"HessianAmbientCache\nAmbientAveraging\nAmbientValueSharing\n" \
		"PixelJitter\nPixelSampling\nPixelMotion\nPixelDepthOfField\n" \
		"SmallSourceDrawing\nViewSequence\nProgressReporting\n" \
		"AdaptiveShadowTesting\nOutputs=v,l\n" \
		"OutputCS=RGB,XYZ,prims,spec\n"

void
quit(int code)			/* quit program */
{
	exit(code);		// don't bother about freeing anything
}

/* Set default options */
static void
default_options(void)
{
	shadthresh = .05;
	shadcert = .5;
	srcsizerat = .25;
	directrelay = 1;
	ambacc = 0.2;
	ambres = 64;
	ambdiv = 512;
	ambssamp = 128;
	maxdepth = 7;
}

int
main(int  argc, char  *argv[])
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
	RGBPRIMS  our_prims;		/* private output color primitives */
	int  seqstart = 0;
	char  *recover = NULL;
	char  *outfile = NULL;
	char  *zfile = NULL;
	int	outfmt = 'c';
	int  rval;
	int  i;
					/* record start time */
	tstart = time(NULL);
					/* global program name */
	progname = argv[0];
					/* feature check only? */
	strcat(RFeatureList, RXPICT_FEATURES);
	if (argc > 1 && !strcmp(argv[1], "-features"))
		return feature_status(argc-2, argv+2);
					/* set defaults */
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
			printdefaults();
			quit(0);
		}
		rval = getrenderopt(argc-i, argv+i);
		if (rval >= 0) {
			i += rval;
			continue;
		}
		rval = getviewopt(&ourview, argc-i, argv+i);
		if (rval >= 0) {
			i += rval;
			continue;
		}
						/* rxpict options */
		switch (argv[i][1]) {
		case 'v':				/* view file */
			if (argv[i][2] != 'f')
				goto badopt;
			check(3,"s");
			rval = viewfile(argv[++i], &ourview, NULL);
			if (rval < 0) {
				sprintf(errmsg,
				"cannot open view file \"%s\"",
						argv[i]);
				error(SYSTEM, errmsg);
			} else if (rval == 0) {
				sprintf(errmsg,
					"bad view file \"%s\"",
						argv[i]);
				error(USER, errmsg);
			}
			break;
		case 'n':				/* number of processes */
			check(2,"i");
			nproc = atoi(argv[++i]);
			if (nproc < 0 && (nproc += RadSimulManager::GetNCores()) <= 0)
				nproc = 1;
			break;
		case 'f':				/* output format */
			if ((argv[i][2] != 'c') & (argv[i][2] != 'f')
					|| argv[i][3])
				goto badopt;
			outfmt = argv[i][2];
			break;
		case 'p':				/* pixel */
			switch (argv[i][2]) {
			case 's':				/* sample */
				check(3,"i");
				psample = atoi(argv[++i]);
				if (psample < 1) psample = 1;
				break;
			case 't':				/* threshold */
				check(3,"f");
				maxdiff = atof(argv[++i]);
				break;
			case 'j':				/* jitter */
				check(3,"f");
				dstrpix = atof(argv[++i]);
				break;
			case 'a':				/* aspect */
				check(3,"f");
				pixaspect = atof(argv[++i]);
				break;
			case 'm':				/* motion */
				check(3,"f");
				mblur = atof(argv[++i]);
				mblur *= (mblur > 0);
				break;
			case 'd':				/* aperture */
				check(3,"f");
				dblur = atof(argv[++i]);
				dblur *= (dblur > 0);
				break;
			case 'R':				/* standard RGB output */
				if (strcmp(argv[i]+2, "RGB"))
					goto badopt;
				myRPmanager.prims = stdprims;
				dtype = RDTnewCT(dtype, RDTrgbe);
				break;
			case 'X':				/* XYZ output */
				if (strcmp(argv[i]+2, "XYZ"))
					goto badopt;
				myRPmanager.prims = xyzprims;
				dtype = RDTnewCT(dtype, RDTxyze);
				break;
			case 'c':				/* chromaticities */
				check(3,"ffffffff");
				rval = 0;
				for (int j = 0; j < 8; j++) {
					our_prims[0][j] = atof(argv[++i]);
					rval |= fabs(our_prims[0][j]-stdprims[0][j]) > .001;
				}
				if (rval) {
					if (!colorprimsOK(our_prims))
						error(USER, "illegal primary chromaticities");
					myRPmanager.prims = our_prims;
				} else
					myRPmanager.prims = stdprims;
				dtype = RDTnewCT(dtype, RDTrgbe);
				break;
			default:
				goto badopt;
			}
			break;
		case 'd':				/* reference depth */
			if (argv[i][2] || !myRPmanager.SetReferenceDepth(argv[++i]))
				goto badopt;
			dtype = RDTnewDT(dtype, RDTdshort);
			break;
		case 'x':				/* x resolution */
			check(2,"i");
			hresolu = atoi(argv[++i]);
			break;
		case 'y':				/* y resolution */
			check(2,"i");
			vresolu = atoi(argv[++i]);
			break;
		case 'S':				/* start index */
			check(2,"i");
			seqstart = atoi(argv[++i]);
			seqstart *= (seqstart > 0);
			break;
		case 'o':				/* output file */
			check(2,"s");
			outfile = argv[++i];
			break;
		case 'z':				/* z file */
			check(2,"s");
			zfile = argv[++i];
			break;
		case 'r':				/* recover file */
			if (argv[i][2] == 'o') {		/* +output */
				check(3,"s");
				outfile = argv[i+1];
			} else
				check(2,"s");
			recover = argv[++i];
			break;
#if MAXCSAMP>3
		case 'c':				/* output spectral results */
			if (argv[i][2] != 'o')
				goto badopt;
			rval = (myRPmanager.prims == NULL);
			check_bool(3,rval);
			if (rval)
				myRPmanager.prims = NULL;
			else if (myRPmanager.prims == NULL)
				myRPmanager.prims = stdprims;
			dtype = RDTnewCT(dtype, rval ? RDTscolr : RDTrgbe);
			break;
#endif
		case 't':				/* timer */
			check(2,"i");
			ralrm = atoi(argv[++i]);
			break;
		case 'w':				/* warnings */
			rval = erract[WARNING].pf != NULL;
			check_bool(2,rval);
			if (rval) erract[WARNING].pf = wputs;
			else erract[WARNING].pf = NULL;
			break;
		case 'e':				/* error file */
			check(2,"s");
			errfile = argv[++i];
			break;
		default:
			goto badopt;
		}
	}
	if (maxdiff <= FTINY)		/* check for useless sampling */
		psample = 1;
	if (zfile == NULL)		/* set up depth output */
		dtype = RDTnewDT(dtype, RDTnone);
	else if (!RDTdepthT(dtype))
		dtype = RDTnewDT(dtype, RDTdfloat);
					/* check pixel output type */
	if ((myRPmanager.prims == NULL) & (NCSAMP == 3)) {
		myRPmanager.prims = stdprims;
		dtype = RDTnewCT(dtype, RDTrgbe);
	}
	if (outfmt == 'f')
		switch (RDTcolorT(dtype)) {
		case RDTrgbe:
			dtype = RDTnewCT(dtype, RDTrgb);
			break;
		case RDTxyze:
			dtype = RDTnewCT(dtype, RDTxyz);
			break;
		case RDTscolr:
			dtype = RDTnewCT(dtype, RDTscolor);
			break;
		case RDTrgb:
		case RDTxyz:
		case RDTscolor:
			break;
		default:
			error(INTERNAL, "botched color output type");
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
					/* open error file */
	if (errfile != NULL) {
		if (freopen(errfile, "a", stderr) == NULL)
			quit(2);
		fprintf(stderr, "**************\n*** PID %5d: ",
				getpid());
		printargs(argc, argv, stderr);
		putc('\n', stderr);
		fflush(stderr);
	}
#ifdef	NICE
	nice(NICE);			/* lower priority */
#endif
	if (i < argc-1)
		goto badopt;
					// load octree
	if (!myRPmanager.LoadOctree(argv[i]))
		error(USER, "missing octree argument");
					// add new header info
	myRPmanager.AddHeader(i, argv);
	{
		char	buf[128] = "SOFTWARE= ";
		strcpy(buf+10, VersionID);
		myRPmanager.AddHeader(buf);
	}
					// start our engines
	nproc = myRPmanager.SetThreadCount(nproc);
					// batch render picture(s)
	rpict(seqstart, outfile, zfile, recover);

	quit(0);			// clean up and exit

badopt:
	sprintf(errmsg, "command line error at '%s'", argv[i]);
	error(USER, errmsg);
	return 1; /* pro forma return */

#undef	check
#undef	check_bool
}


void
wputs(				/* warning output function */
	const char	*s
)
{
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


static void
onsig(				/* fatal signal */
	int  signo
)
{
	static int  gotsig = 0;

	if (gotsig++)			/* two signals and we're gone! */
		_exit(signo);

#ifdef SIGALRM /* XXX how critical is this? */
	alarm(30);			/* allow 30 seconds to clean up */
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


static void
printdefaults(void)			/* print default values to stdout */
{
	printf("-n %-2d\t\t\t\t# number of rendering processes\n", nproc);
	printf("-vt%c\t\t\t\t# view type %s\n", ourview.type,
			ourview.type==VT_PER ? "perspective" :
			ourview.type==VT_PAR ? "parallel" :
			ourview.type==VT_HEM ? "hemispherical" :
			ourview.type==VT_ANG ? "angular" :
			ourview.type==VT_CYL ? "cylindrical" :
			ourview.type==VT_PLS ? "planisphere" :
			"unknown");
	printf("-vp %f %f %f\t# view point\n",
			ourview.vp[0], ourview.vp[1], ourview.vp[2]);
	printf("-vd %f %f %f\t# view direction\n",
			ourview.vdir[0], ourview.vdir[1], ourview.vdir[2]);
	printf("-vu %f %f %f\t# view up\n",
			ourview.vup[0], ourview.vup[1], ourview.vup[2]);
	printf("-vh %f\t\t\t# view horizontal size\n", ourview.horiz);
	printf("-vv %f\t\t\t# view vertical size\n", ourview.vert);
	printf("-vo %f\t\t\t# view fore clipping plane\n", ourview.vfore);
	printf("-va %f\t\t\t# view aft clipping plane\n", ourview.vaft);
	printf("-vs %f\t\t\t# view shift\n", ourview.hoff);
	printf("-vl %f\t\t\t# view lift\n", ourview.voff);
	printf("-x  %-9d\t\t\t# x resolution\n", hresolu);
	printf("-y  %-9d\t\t\t# y resolution\n", vresolu);
	if (myRPmanager.prims == stdprims)
		printf("-pRGB\t\t\t\t# standard RGB color output\n");
	else if (myRPmanager.prims == xyzprims)
		printf("-pXYZ\t\t\t\t# CIE XYZ color output\n");
	else if (myRPmanager.prims != NULL)
		printf("-pc %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\t# output color primaries and white point\n",
				myRPmanager.prims[RED][0], myRPmanager.prims[RED][1],
				myRPmanager.prims[GRN][0], myRPmanager.prims[GRN][1],
				myRPmanager.prims[BLU][0], myRPmanager.prims[BLU][1],
				myRPmanager.prims[WHT][0], myRPmanager.prims[WHT][1]);
	if (NCSAMP > 3)
		printf(myRPmanager.prims != NULL ? "-co-\t\t\t\t# output tristimulus colors\n" :
				"-co+\t\t\t\t# output spectral values\n");
	printf("-pa %f\t\t\t# pixel aspect ratio\n", pixaspect);
	printf("-pj %f\t\t\t# pixel jitter\n", dstrpix);
	printf("-pm %f\t\t\t# pixel motion\n", mblur);
	printf("-pd %f\t\t\t# pixel depth-of-field\n", dblur);
	printf("-ps %-9d\t\t\t# pixel sample\n", psample);
	printf("-pt %f\t\t\t# pixel threshold\n", maxdiff);
	printf("-t  %-9d\t\t\t# time between reports\n", ralrm);
	printf(erract[WARNING].pf != NULL ?
			"-w+\t\t\t\t# warning messages on\n" :
			"-w-\t\t\t\t# warning messages off\n");
	print_rdefaults();
}


// my progress report call-back
static void
progReporter(double pct)
{
	static time_t	lastReportTime = 0;
	time_t		tnow = time(NULL);

	if (pct < 100.-FTINY && tnow - lastReportTime < ralrm)
		return;			// too soon, my Precious...

	sprintf(errmsg, "%7.3f%% done after %7.3f hours\n", pct, (tnow-tstart)/3600.);
	eputs(errmsg);
					// reset timer at 100%
	lastReportTime = tnow * (pct < 100.-FTINY);
}


// Principal function for rpict execution (prototype in ray.h for some reason)
void
rpict(int seq, char *pout, char *zout, char *prvr)
/*
 * If seq is greater than zero, we will render a sequence of
 * images based on view parameter strings read from the standard input.
 * If pout is NULL, then all images will be sent to the standard ouput.
 * If seq is greater than zero and prvr is an integer, then it is the
 * frame number at which rendering should begin.  Preceeding view parameter
 * strings will be skipped in the input.
 * Note that pout and zout should contain %d format specifications for
 * sequenced file naming.
 */
{
	char  fbuf[256], dbuf[256];

	if (!zout ^ !RDTdepthT(dtype))
		error(INTERNAL, "output depth type requires spec and vice versa");

	if (prvr && isint(prvr)) {	// sequence recovery?
		if (!seq)
			error(USER, "sequence recovery requires starting frame");
		if (!pout)
			error(USER, "need output spec for sequence recovery");
		int	rno = atoi(prvr);
		if (rno < seq)
			error(USER, "recovery frame before starting frame");
		while (seq <= rno) {
			if (!fgets(fbuf, sizeof(fbuf), stdin))
				error(USER, "recovery frame past end of sequence");
			seq += (isview(fbuf) && sscanview(&ourview, fbuf) > 0);
		}
		sprintf(prvr=fbuf, pout, rno);
	}
	if (ralrm > 0)
		myRPmanager.prCB = progReporter;
	else
		myRPmanager.prCB = NULL;
	if (prvr) {			// recovering partial render?
		if ((seq > 0) & (prvr != fbuf))
			error(USER, "recover spec must be number in sequence");
		if (zout) sprintf(dbuf, zout, seq);
		if (ralrm > 0) {
			sprintf(errmsg, "resuming partial rendering '%s'\n", prvr);
			eputs(errmsg);
		}
		dtype = myRPmanager.ResumeFrame(prvr, zout ? dbuf : zout);
		if (!dtype)
			error(USER, "ResumeFrame() failed");
		if (!seq)
			return;		// all done if not a sequence
	}
	do {
		if (prvr)		// have view from sequence recovery?
			prvr = NULL;
		else if (seq) {		// else read next view in sequence
			while (fgets(fbuf, sizeof(fbuf), stdin))
				if (isview(fbuf) && sscanview(&ourview, fbuf) > 0)
					break;
			if (feof(stdin))
				return;	// reached end of view input
		}
					// get (indexed) output file name(s)
		if (pout) sprintf(fbuf, pout, seq);
		if (zout) sprintf(dbuf, zout, seq);

		if (ralrm > 0) {	// start report if requested
			if (myRPmanager.NThreads() > 1)
				sprintf(errmsg, "%2d processes ", myRPmanager.NThreads());
			else
				sprintf(errmsg, "PID %6d ", getpid());
			eputs(errmsg);
			if (seq) {
				sprintf(errmsg, "rendering frame %5d ", seq);
				eputs(errmsg);
			} else
				eputs("rendering picture ");
			if (pout) {
				sprintf(errmsg, "to '%s'\n", fbuf);
				eputs(errmsg);
			} else
				eputs("to stdout\n");
			if (zout) {
				sprintf(errmsg, "\twith %s depth map in '%s'\n",
						RDTdepthT(dtype)==RDTdshort ?
							"encoded 16-bit" : "raw float",
						dbuf);
				eputs(errmsg);
			}
		}
					// set up view and size
		int	xydim[2] = {hresolu, vresolu};
		if (!myRPmanager.NewFrame(ourview, xydim, &pixaspect))
			error(USER, "NewFrame() failed");

		myRPmanager.frameNo = seq;

		errno = 0;		// render frame, skipping if it exists
		if (!myRPmanager.RenderFrame(pout ? fbuf : pout, dtype,
							zout ? dbuf : zout)
					&& !seq | (errno != EEXIST))
			error(USER, "RenderFrame() failed");
	} while (seq++);		// all done if not a sequence
}
