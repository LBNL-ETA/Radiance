#ifndef lint
static const char	RCSid[] = "$Id: rtmain.c,v 2.18 2009/12/12 19:01:00 greg Exp $";
#endif
/*
 *  rtmain.c - main for rtrace per-ray calculation program
 */

#include "copyright.h"

#include  <signal.h>

#include  "platform.h"
#include  "rtprocess.h" /* getpid() */
#include  "resolu.h"
#include  "ray.h"
#include  "source.h"
#include  "ambient.h"
#include  "random.h"
#include  "paths.h"

extern char	*progname;		/* global argv[0] */

extern char	*shm_boundary;		/* boundary of shared memory */

					/* persistent processes define */
#ifdef  F_SETLKW
#define  PERSIST	1		/* normal persist */
#define  PARALLEL	2		/* parallel persist */
#define  PCHILD		3		/* child of normal persist */
#endif

char  *sigerr[NSIG];			/* signal error messages */
char  *errfile = NULL;			/* error output file */

int  nproc = 1;				/* number of processes */

extern char  *formstr();		/* string from format */
int  inform = 'a';			/* input format */
int  outform = 'a';			/* output format */
char  *outvals = "v";			/* output specification */

int  hresolu = 0;			/* horizontal (scan) size */
int  vresolu = 0;			/* vertical resolution */

int  imm_irrad = 0;			/* compute immediate irradiance? */
int  lim_dist = 0;			/* limit distance? */

#ifndef	MAXMODLIST
#define	MAXMODLIST	1024		/* maximum modifiers we'll track */
#endif

extern void  (*addobjnotify[])();	/* object notification calls */
extern void  tranotify(OBJECT obj);

char  *tralist[MAXMODLIST];		/* list of modifers to trace (or no) */
int  traincl = -1;			/* include == 1, exclude == 0 */

static int  loadflags = ~IO_FILES;	/* what to load from octree */

static void onsig(int  signo);
static void sigdie(int  signo, char  *msg);
static void printdefaults(void);


int
main(int  argc, char  *argv[])
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
	int  persist = 0;
	char  *octnm = NULL;
	char  **tralp;
	int  duped1;
	int  rval;
	int  i;
					/* global program name */
	progname = argv[0] = fixargv0(argv[0]);
					/* add trace notify function */
	for (i = 0; addobjnotify[i] != NULL; i++)
		;
	addobjnotify[i] = tranotify;
					/* set our defaults */
	maxdepth = -10;
	minweight = 2e-3;
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
			if (nproc <= 0)
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
			bool(2,rval);
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
			bool(3,lim_dist);
			break;
		case 'I':				/* immed. irradiance */
			bool(2,imm_irrad);
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
			rval = loadflags & IO_INFO;
			bool(2,rval);
			loadflags = rval ? loadflags | IO_INFO :
					loadflags & ~IO_INFO;
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
					rval = wordfile(tralp,
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
					rval = wordfile(tralp,
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
#ifdef  PERSIST
		case 'P':				/* persist file */
			if (argv[i][2] == 'P') {
				check(3,"s");
				persist = PARALLEL;
			} else {
				check(2,"s");
				persist = PERSIST;
			}
			persistfile(argv[++i]);
			break;
#endif
		default:
			goto badopt;
		}
	}
	if (nproc > 1) {
		if (persist)
			error(USER, "multiprocessing incompatible with persist file");
		if (imm_irrad)
			error(USER, "multiprocessing incompatible with immediate irradiance");
		if (hresolu > 0 && hresolu < nproc)
			error(WARNING, "number of cores should not exceed horizontal resolution");
		if (trace != NULL)
			error(WARNING, "multiprocessing does not work properly with trace mode");
	}
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
					/* get octree */
	if (i == argc)
		octnm = NULL;
	else if (i == argc-1)
		octnm = argv[i];
	else
		goto badopt;
	if (octnm == NULL)
		error(USER, "missing octree argument");
					/* set up output */
#ifdef  PERSIST
	if (persist) {
		duped1 = dup(fileno(stdout));	/* don't lose our output */
		openheader();
	}
#endif
	if (outform != 'a')
		SET_FILE_BINARY(stdout);
	readoct(octnm, loadflags, &thescene, NULL);
	nsceneobjs = nobjects;

	if (loadflags & IO_INFO) {	/* print header */
		printargs(i, argv, stdout);
		printf("SOFTWARE= %s\n", VersionID);
		fputnow(stdout);
		fputformat(formstr(outform), stdout);
		putchar('\n');
	}

	marksources();			/* find and mark sources */

	setambient();			/* initialize ambient calculation */

#ifdef  PERSIST
	if (persist) {
		fflush(stdout);
						/* reconnect stdout */
		dup2(duped1, fileno(stdout));
		close(duped1);
		if (persist == PARALLEL) {	/* multiprocessing */
			preload_objs();		/* preload scene */
			shm_boundary = (char *)malloc(16);
			strcpy(shm_boundary, "SHM_BOUNDARY");
			while ((rval=fork()) == 0) {	/* keep on forkin' */
				pflock(1);
				pfhold();
				ambsync();		/* load new values */
			}
			if (rval < 0)
				error(SYSTEM, "cannot fork child for persist function");
			pfdetach();		/* parent will run then exit */
		}
	}
runagain:
	if (persist)
		dupheader();			/* send header to stdout */
#endif
	if (nproc > 1)			/* start multiprocessing */
		ray_popen(nproc);
					/* trace rays */
	rtrace(NULL);
					/* flush ambient file */
	ambsync();
#ifdef  PERSIST
	if (persist == PERSIST) {	/* first run-through */
		if ((rval=fork()) == 0) {	/* child loops until killed */
			pflock(1);
			persist = PCHILD;
		} else {			/* original process exits */
			if (rval < 0)
				error(SYSTEM, "cannot fork child for persist function");
			pfdetach();		/* parent exits */
		}
	}
	if (persist == PCHILD) {	/* wait for a signal then go again */
		pfhold();
		raynum = nrays = 0;		/* reinitialize */
		goto runagain;
	}
#endif
	quit(0);

badopt:
	sprintf(errmsg, "command line error at '%s'", argv[i]);
	error(USER, errmsg);
	return 1; /* pro forma return */

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
	register char  *s
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
	char  *msg
)
{
	if (signal(signo, onsig) == SIG_IGN)
		signal(signo, SIG_IGN);
	sigerr[signo] = msg;
}


static void
printdefaults(void)			/* print default values to stdout */
{
	register char  *cp;

	if (imm_irrad)
		printf("-I+\t\t\t\t# immediate irradiance on\n");
	printf("-n %-2d\t\t\t\t# number of rendering processes\n", nproc);
	printf("-x %-9d\t\t\t# x resolution (flush interval)\n", hresolu);
	printf("-y %-9d\t\t\t# y resolution\n", vresolu);
	printf(lim_dist ? "-ld+\t\t\t\t# limit distance on\n" :
			"-ld-\t\t\t\t# limit distance off\n");
	printf("-h%c\t\t\t\t# %s header\n", loadflags & IO_INFO ? '+' : '-',
			loadflags & IO_INFO ? "output" : "no");
	printf("-f%c%c\t\t\t\t# format input/output = %s/%s\n",
			inform, outform, formstr(inform), formstr(outform));
	printf("-o%-9s\t\t\t# output", outvals);
	for (cp = outvals; *cp; cp++)
		switch (*cp) {
		case 't': case 'T': printf(" trace"); break;
		case 'o': printf(" origin"); break;
		case 'd': printf(" direction"); break;
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
		case '-': printf(" stroke"); break;
		}
	putchar('\n');
	printf(erract[WARNING].pf != NULL ?
			"-w+\t\t\t\t# warning messages on\n" :
			"-w-\t\t\t\t# warning messages off\n");
	print_rdefaults();
}
