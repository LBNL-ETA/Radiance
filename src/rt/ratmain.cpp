#ifndef lint
static const char	RCSid[] = "$Id: ratmain.cpp,v 2.1 2023/02/08 17:41:48 greg Exp $";
#endif
/*
 *  rtmain.c - main for rtrace per-ray calculation program
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

static void onsig(int  signo);
static void sigdie(int  signo, const char  *msg);
static void printdefaults(void);

#define RATRACE_FEATURES	"IrradianceCalc\nIrradianceCalc\nDistanceLimiting\n" \
			"HessianAmbientCache\nAmbientAveraging\n" \
			"AmbientValueSharing\nAdaptiveShadowTesting\n" \
			"Outputs=o,d,v,V,w,W,l,L,c,p,n,N,s,m,M,r,x,R,X,~\n"

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
	progname = argv[0] = fixargv0(argv[0]);
					/* feature check only? */
	strcat(RFeatureList, RATRACE_FEATURES);
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
			if (nproc < 0)
				error(USER, "bad number of processes");
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
		default:
			goto badopt;
		}
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
		static char	fmt[] = OCTFMT;
		FILE *		octfp = fopen(argv[i], "rb");
		if (checkheader(octfp, fmt, stdout) < 0)
			error(USER, "bad octree header");
		fclose(octfp);
		printargs(i, argv, stdout);
		printf("SOFTWARE= %s\n", VersionID);
		fputnow(stdout);
		if (rval > 0)		/* saved from setrtoutput() call */
			printf("NCOMP=%d\n", rval);
		if ((outform == 'f') | (outform == 'd'))
			fputendian(stdout);
		fputformat(formstr(outform), stdout);
		putchar('\n');
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
	printf(erract[WARNING].pf != NULL ?
			"-w+\t\t\t\t# warning messages on\n" :
			"-w-\t\t\t\t# warning messages off\n");
	print_rdefaults();
}
