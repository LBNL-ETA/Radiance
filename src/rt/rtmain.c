#ifndef lint
static const char	RCSid[] = "$Id: rtmain.c,v 2.6 2003/06/30 14:59:12 schorsch Exp $";
#endif
/*
 *  rtmain.c - main for rtrace per-ray calculation program
 */

#include "copyright.h"

#include  <sys/types.h>
#include  <signal.h>
#ifdef _WIN32
 #include <process.h> /* getpid() */
#endif

#include  "platform.h"
#include  "ray.h"
#include  "source.h"
#include  "ambient.h"
#include  "random.h"
#include  "paths.h"

					/* persistent processes define */
#ifdef  F_SETLKW
#define  PERSIST	1		/* normal persist */
#define  PARALLEL	2		/* parallel persist */
#define  PCHILD		3		/* child of normal persist */
#endif

char  *progname;			/* argv[0] */

char  *octname;				/* octree name */

char  *sigerr[NSIG];			/* signal error messages */

char  *shm_boundary = NULL;		/* boundary of shared memory */

char  *errfile = NULL;			/* error output file */

extern char  *formstr();		/* string from format */
extern int  inform;			/* input format */
extern int  outform;			/* output format */
extern char  *outvals;			/* output values */

extern int  hresolu;			/* horizontal resolution */
extern int  vresolu;			/* vertical resolution */

extern int  imm_irrad;			/* compute immediate irradiance? */
extern int  lim_dist;			/* limit distance? */

extern char  *tralist[];		/* list of modifers to trace (or no) */
extern int  traincl;			/* include == 1, exclude == 0 */

void	onsig();
void	sigdie();
void	printdefaults();


int
main(argc, argv)
int  argc;
char  *argv[];
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
	int  loadflags = ~IO_FILES;
	int  persist = 0;
	char  **tralp;
	int  duped1;
	int  rval;
	int  i;
					/* global program name */
	progname = argv[0] = fixargv0(argv[0]);
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
					/* initialize object types */
	initotypes();
					/* initialize urand */
	initurand(2048);
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
		octname = NULL;
	else if (i == argc-1)
		octname = argv[i];
	else
		goto badopt;
	if (octname == NULL)
		error(USER, "missing octree argument");
					/* set up output */
#ifdef  PERSIST
	if (persist) {
		duped1 = dup(fileno(stdout));	/* don't lose our output */
		openheader();
	}
#endif
#ifdef	_WIN32
	if (outform != 'a')
		SET_FILE_BINARY(stdout);
	if (octname == NULL)
		SET_FILE_BINARY(stdin);
#endif
	readoct(octname, loadflags, &thescene, NULL);
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
			}
			if (rval < 0)
				error(SYSTEM, "cannot fork child for persist function");
			pfdetach();		/* parent exits */
		}
	}
runagain:
	if (persist)
		dupheader();			/* send header to stdout */
#endif
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
		close(duped1);			/* release output handle */
		pfhold();
		raynum = nrays = 0;		/* reinitialize */
		goto runagain;
	}
#endif
	quit(0);

badopt:
	sprintf(errmsg, "command line error at '%s'", argv[i]);
	error(USER, errmsg);

#undef	check
#undef	bool
}


void
wputs(s)				/* warning output function */
char	*s;
{
	int  lasterrno = errno;
	eputs(s);
	errno = lasterrno;
}


void
eputs(s)				/* put string to stderr */
register char  *s;
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


void
onsig(signo)				/* fatal signal */
int  signo;
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


void
sigdie(signo, msg)			/* set fatal signal */
int  signo;
char  *msg;
{
	if (signal(signo, onsig) == SIG_IGN)
		signal(signo, SIG_IGN);
	sigerr[signo] = msg;
}


void
printdefaults()			/* print default values to stdout */
{
	register char  *cp;

	if (imm_irrad)
		printf("-I+\t\t\t\t# immediate irradiance on\n");
	printf("-x  %-9d\t\t\t# x resolution\n", hresolu);
	printf("-y  %-9d\t\t\t# y resolution\n", vresolu);
	printf(lim_dist ? "-ld+\t\t\t\t# limit distance on\n" :
			"-ld-\t\t\t\t# limit distance off\n");
	printf("-f%c%c\t\t\t\t# format input/output = %s/%s\n",
			inform, outform, formstr(inform), formstr(outform));
	printf("-o%s\t\t\t\t# output", outvals);
	for (cp = outvals; *cp; cp++)
		switch (*cp) {
		case 't': printf(" trace"); break;
		case 'o': printf(" origin"); break;
		case 'd': printf(" direction"); break;
		case 'v': printf(" value"); break;
		case 'l': printf(" length"); break;
		case 'L': printf(" first_length"); break;
		case 'p': printf(" point"); break;
		case 'n': printf(" normal"); break;
		case 'N': printf(" unperturbed_normal"); break;
		case 's': printf(" surface"); break;
		case 'w': printf(" weight"); break;
		case 'm': printf(" modifier"); break;
		}
	putchar('\n');
	printf(erract[WARNING].pf != NULL ?
			"-w+\t\t\t\t# warning messages on\n" :
			"-w-\t\t\t\t# warning messages off\n");
	print_rdefaults();
}
