#ifndef lint
static const char	RCSid[] = "$Id: rxtmain.cpp,v 2.9 2024/09/16 19:20:09 greg Exp $";
#endif
/*
 *  rxtmain.cpp - main for per-ray calculation program
 */

#include "copyright.h"

#include  <signal.h>

#include  "rtprocess.h" /* getpid() */
#include  "platform.h"
#include  "RtraceSimulManager.h"

extern char	*progname;		/* global argv[0] */

static const char  *sigerr[NSIG];	/* signal error messages */
char  *errfile = NULL;			/* error output file */

extern const char  *formstr(int f);	/* string from format */
extern int  setrtoutput(const char *outvals);	/* set output values */

int  inform = 'a';			/* input format */
int  outform = 'a';			/* output format */

int  hresolu = 0;			/* horizontal (scan) size */
int  vresolu = 0;			/* vertical resolution */

RtraceSimulManager	myRTmanager;	// global simulation manager

static const char  *outvals = "v";	/* output specification */
static 	int  nproc = 1;			/* number of requested threads */
static int  doheader = 1;		/* include information header? */

#ifndef	MAXMODLIST
#define	MAXMODLIST	1024		/* maximum modifiers we'll track */
#endif

extern void  (*addobjnotify[])(OBJECT);	/* object notification calls */
extern void  tranotify(OBJECT obj);

char  *tralist[MAXMODLIST];		/* list of modifers to trace (or no) */
int  traincl = -1;			/* include == 1, exclude == 0 */

double  (*sens_curve)(const SCOLOR scol) = NULL;	/* spectral conversion for 1-channel */
double  out_scalefactor = 1;		/* output calibration scale factor */
RGBPRIMP  out_prims = stdprims;		/* output color primitives (NULL if spectral) */
static RGBPRIMS  our_prims;		/* private output color primitives */

static void onsig(int  signo);
static void sigdie(int  signo, const char  *msg);
static void printdefaults(void);

#define RXTRACE_FEATURES	"IrradianceCalc\nIrradianceCalc\nDistanceLimiting\n" \
				"HessianAmbientCache\nAmbientAveraging\n" \
				"AmbientValueSharing\nAdaptiveShadowTesting\n" \
				"Outputs=o,d,v,V,w,W,l,L,c,p,n,N,s,m,M,r,x,R,X,~\n" \
				"OutputCS=RGB,XYZ,Y,S,M,prims,spec\n"

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
	char  **tralp = NULL;
	int  rval;
	int  i;
					/* global program name */
	progname = argv[0];
					/* feature check only? */
	strcat(RFeatureList, RXTRACE_FEATURES);
	if (argc > 1 && !strcmp(argv[1], "-features"))
		return feature_status(argc-2, argv+2);
					/* add trace notify function */
	for (i = 0; addobjnotify[i] != NULL; i++)
		;
	addobjnotify[i] = tranotify;
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
		switch (argv[i][1]) {
		case 'n':				/* number of cores */
			check(2,"i");
			nproc = atoi(argv[++i]);
			break;
		case 'x':				/* x resolution */
			check(2,"i");
			hresolu = atoi(argv[++i]);
			break;
		case 'y':				/* y resolution */
			check(2,"i");
			vresolu = atoi(argv[++i]);
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
		case 'l':				/* limit distance */
			if (argv[i][2] != 'd')
				goto badopt;
			rval = myRTmanager.rtFlags & RTlimDist;
			check_bool(3,rval);
			if (rval) myRTmanager.rtFlags |= RTlimDist;
			else myRTmanager.rtFlags &= ~RTlimDist;
			break;
		case 'I':				/* immed. irradiance */
			rval = myRTmanager.rtFlags & RTimmIrrad;
			check_bool(3,rval);
			if (rval) myRTmanager.rtFlags |= RTimmIrrad;
			else myRTmanager.rtFlags &= ~RTimmIrrad;
			break;
		case 'f':				/* format i/o */
			switch (argv[i][2]) {
			case 'a':				/* ascii */
			case 'f':				/* float */
			case 'd':				/* double */
				inform = argv[i][2];
				break;
			default:
				goto badopt;
			}
			switch (argv[i][3]) {
			case '\0':
				outform = inform;
				break;
			case 'a':				/* ascii */
			case 'f':				/* float */
			case 'd':				/* double */
			case 'c':				/* color */
				check(4,"");
				outform = argv[i][3];
				break;
			default:
				goto badopt;
			}
			break;
		case 'o':				/* output */
			outvals = argv[i]+2;
			break;
		case 'h':				/* header output */
			check_bool(2,doheader);
			break;
		case 't':				/* trace */
			switch (argv[i][2]) {
			case 'i':				/* include */
			case 'I':
				check(3,"s");
				if (traincl != 1) {
					traincl = 1;
					tralp = tralist;
				}
				if (argv[i][2] == 'I') {	/* file */
					rval = wordfile(tralp, MAXMODLIST-(tralp-tralist),
					getpath(argv[++i],getrlibpath(),R_OK));
					if (rval < 0) {
						sprintf(errmsg,
				"cannot open trace include file \"%s\"",
								argv[i]);
						error(SYSTEM, errmsg);
					}
					tralp += rval;
				} else {
					*tralp++ = argv[++i];
					*tralp = NULL;
				}
				break;
			case 'e':				/* exclude */
			case 'E':
				check(3,"s");
				if (traincl != 0) {
					traincl = 0;
					tralp = tralist;
				}
				if (argv[i][2] == 'E') {	/* file */
					rval = wordfile(tralp, MAXMODLIST-(tralp-tralist),
					getpath(argv[++i],getrlibpath(),R_OK));
					if (rval < 0) {
						sprintf(errmsg,
				"cannot open trace exclude file \"%s\"",
								argv[i]);
						error(SYSTEM, errmsg);
					}
					tralp += rval;
				} else {
					*tralp++ = argv[++i];
					*tralp = NULL;
				}
				break;
			default:
				goto badopt;
			}
			break;
		case 'p':				/* value output */
			switch (argv[i][2]) {
			case 'R':			/* standard RGB output */
				if (strcmp(argv[i]+2, "RGB"))
					goto badopt;
				out_prims = stdprims;
				out_scalefactor = 1;
				sens_curve = NULL;
				break;
			case 'X':			/* XYZ output */
				if (strcmp(argv[i]+2, "XYZ"))
					goto badopt;
				out_prims = xyzprims;
				out_scalefactor = WHTEFFICACY;
				sens_curve = NULL;
				break;
			case 'c': {
				int	j;
				check(3,"ffffffff");
				rval = 0;
				for (j = 0; j < 8; j++) {
					our_prims[0][j] = atof(argv[++i]);
					rval |= fabs(our_prims[0][j]-stdprims[0][j]) > .001;
				}
				if (rval) {
					if (!colorprimsOK(our_prims))
						error(USER, "illegal primary chromaticities");
					out_prims = our_prims;
				} else
					out_prims = stdprims;
				out_scalefactor = 1;
				sens_curve = NULL;
				} break;
			case 'Y':			/* photopic response */
				if (argv[i][3])
					goto badopt;
				sens_curve = scolor_photopic;
				out_scalefactor = WHTEFFICACY;
				break;
			case 'S':			/* scotopic response */
				if (argv[i][3])
					goto badopt;
				sens_curve = scolor_scotopic;
				out_scalefactor = WHTSCOTOPIC;
				break;
			case 'M':			/* melanopic response */
				if (argv[i][3])
					goto badopt;
				sens_curve = scolor_melanopic;
				out_scalefactor = WHTMELANOPIC;
				break;
			default:
				goto badopt;
			}
			break;
#if MAXCSAMP>3
		case 'c':				/* output spectral results */
			if (argv[i][2] != 'o')
				goto badopt;
			rval = (out_prims == NULL) & (sens_curve == NULL);
			check_bool(3,rval);
			if (rval) {
				out_prims = NULL;
				sens_curve = NULL;
			} else if (out_prims == NULL)
				out_prims = stdprims;
			break;
#endif
		default:
			goto badopt;
		}
	}
					/* set/check spectral sampling */
	rval = setspectrsamp(CNDX, WLPART);
	if (rval < 0)
		error(USER, "unsupported spectral sampling");
	if (out_prims != NULL) {
		if (!rval)
			error(WARNING, "spectral range incompatible with color output");
	} else if (NCSAMP == 3)
		out_prims = stdprims;	/* 3 samples do not a spectrum make */
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
					/* get octree name */
	if (i == argc)
		error(USER, "missing octree argument");
	if (i != argc-1)
		goto badopt;
					/* set output options */
	rval = setrtoutput(outvals);
					/* load octree */
	if (!myRTmanager.LoadOctree(argv[i]))
		quit(1);
					/* set up output */
	if (outform != 'a')
		SET_FILE_BINARY(stdout);
	if (doheader) {			/* print header? */
		newheader("RADIANCE", stdout);
		fputs(myRTmanager.GetHeadStr(), stdout);
		printargs(i, argv, stdout);
		printf("SOFTWARE= %s\n", VersionID);
		fputnow(stdout);
		if (rval > 0)		/* saved from setrtoutput() call */
			printf("NCOMP=%d\n", rval);
		if ((outform == 'f') | (outform == 'd'))
			fputendian(stdout);
		fputformat(formstr(outform), stdout);
		fputc('\n', stdout);	/* end of header */
	}
	rtrace(NULL, nproc);		/* trace rays */
	quit(0);			/* clean up & exit */

badopt:
	sprintf(errmsg, "command line error at '%s'", argv[i]);
	error(USER, errmsg);
	return 1; /* pro forma return */

#undef	check
#undef	check_bool
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
	const char  *cp;

	if (myRTmanager.rtFlags & RTimmIrrad)
		printf("-I+\t\t\t\t# immediate irradiance on\n");
	printf("-n %-2d\t\t\t\t# number of rendering processes\n", nproc);
	printf("-x %-9d\t\t\t# %s\n", hresolu,
			vresolu && hresolu ? "x resolution" : "flush interval");
	printf("-y %-9d\t\t\t# y resolution\n", vresolu);
	printf(myRTmanager.rtFlags&RTlimDist ? "-ld+\t\t\t\t# limit distance on\n" :
			"-ld-\t\t\t\t# limit distance off\n");
	printf("-h%c\t\t\t\t# %s header\n", doheader ? '+' : '-',
			doheader ? "output" : "no");
	printf("-f%c%c\t\t\t\t# format input/output = %s/%s\n",
			inform, outform, formstr(inform), formstr(outform));
	printf("-o%-9s\t\t\t# output", outvals);
	for (cp = outvals; *cp; cp++)
		switch (*cp) {
		case 't': case 'T': printf(" trace"); break;
		case 'o': printf(" origin"); break;
		case 'd': printf(" direction"); break;
		case 'r': printf(" reflect_contrib"); break;
		case 'R': printf(" reflect_length"); break;
		case 'x': printf(" unreflect_contrib"); break;
		case 'X': printf(" unreflect_length"); break;
		case 'v': printf(" value"); break;
		case 'V': printf(" contribution"); break;
		case 'l': printf(" length"); break;
		case 'L': printf(" first_length"); break;
		case 'p': printf(" point"); break;
		case 'n': printf(" normal"); break;
		case 'N': printf(" unperturbed_normal"); break;
		case 's': printf(" surface"); break;
		case 'w': printf(" weight"); break;
		case 'W': printf(" coefficient"); break;
		case 'm': printf(" modifier"); break;
		case 'M': printf(" material"); break;
		case '~': printf(" tilde"); break;
		}
	putchar('\n');
	if (sens_curve == scolor_photopic)
		printf("-pY\t\t\t\t# photopic output\n");
	else if (sens_curve == scolor_scotopic)
		printf("-pS\t\t\t\t# scotopic output\n");
	else if (sens_curve == scolor_melanopic)
		printf("-pM\t\t\t\t# melanopic output\n");
	else if (out_prims == stdprims)
		printf("-pRGB\t\t\t\t# standard RGB color output\n");
	else if (out_prims == xyzprims)
		printf("-pXYZ\t\t\t\t# CIE XYZ color output\n");
	else if (out_prims != NULL)
		printf("-pc %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\t# output color primaries and white point\n",
				out_prims[RED][0], out_prims[RED][1],
				out_prims[GRN][0], out_prims[GRN][1],
				out_prims[BLU][0], out_prims[BLU][1],
				out_prims[WHT][0], out_prims[WHT][1]);
	if ((sens_curve == NULL) & (NCSAMP > 3))
		printf(out_prims != NULL ? "-co-\t\t\t\t# output tristimulus colors\n" :
				"-co+\t\t\t\t# output spectral values\n");
	printf(erract[WARNING].pf != NULL ?
			"-w+\t\t\t\t# warning messages on\n" :
			"-w-\t\t\t\t# warning messages off\n");
	print_rdefaults();
}
