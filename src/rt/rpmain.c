#ifndef lint
static const char	RCSid[] = "$Id: rpmain.c,v 2.5 2003/07/03 18:03:58 greg Exp $";
#endif
/*
 *  rpmain.c - main for rpict batch rendering program
 */

#include "copyright.h"

#include  <sys/types.h>
#include  <signal.h>
#ifdef _WIN32
  #include <process.h> /* getpid */
#endif

#include  "platform.h"
#include  "ray.h"
#include  "source.h"
#include  "ambient.h"
#include  "random.h"
#include  "paths.h"
#include  "view.h"

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

extern time_t  time();
extern time_t  tstart;			/* start time */

extern int  ralrm;			/* seconds between reports */

extern VIEW  ourview;			/* viewing parameters */

extern int  hresolu;			/* horizontal resolution */
extern int  vresolu;			/* vertical resolution */
extern double  pixaspect;		/* pixel aspect ratio */

extern int  psample;			/* pixel sample size */
extern double  maxdiff;			/* max. sample difference */
extern double  dstrpix;			/* square pixel distribution */

extern double  mblur;			/* motion blur parameter */

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
	char  *err;
	char  *recover = NULL;
	char  *outfile = NULL;
	char  *zfile = NULL;
	int  loadflags = ~IO_FILES;
	int  seqstart = 0;
	int  persist = 0;
	int  duped1;
	int  rval;
	int  i;
					/* record start time */
	tstart = time((time_t *)NULL);
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
		rval = getviewopt(&ourview, argc-i, argv+i);
		if (rval >= 0) {
			i += rval;
			continue;
		}
						/* rpict options */
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
		case 'p':				/* pixel */
			switch (argv[i][2]) {
			case 's':				/* sample */
				check(3,"i");
				psample = atoi(argv[++i]);
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
				break;
			default:
				goto badopt;
			}
			break;
		case 'x':				/* x resolution */
			check(2,"i");
			hresolu = atoi(argv[++i]);
			break;
		case 'y':				/* y resolution */
			check(2,"i");
			vresolu = atoi(argv[++i]);
			break;
		case 'S':				/* slave index */
			check(2,"i");
			seqstart = atoi(argv[++i]);
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
		case 't':				/* timer */
			check(2,"i");
			ralrm = atoi(argv[++i]);
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
		default:
			goto badopt;
		}
	}
	err = setview(&ourview);	/* set viewing parameters */
	if (err != NULL)
		error(USER, err);
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
	if (seqstart > 0 && octname == NULL)
		error(USER, "missing octree argument");
					/* set up output */
#ifdef  PERSIST
	if (persist) {
		if (recover != NULL)
			error(USER, "persist option used with recover file");
		if (seqstart <= 0)
			error(USER, "persist option only for sequences");
		if (outfile == NULL)
		duped1 = dup(fileno(stdout));	/* don't lose our output */
		openheader();
	} else
#endif
	if (outfile != NULL)
		openheader();
#ifdef	_WIN32
	SET_FILE_BINARY(stdout);
	if (octname == NULL)
		SET_FILE_BINARY(stdin);
#endif
	readoct(octname, loadflags, &thescene, NULL);
	nsceneobjs = nobjects;

	if (loadflags & IO_INFO) {	/* print header */
		printargs(i, argv, stdout);
		printf("SOFTWARE= %s\n", VersionID);
	}

	marksources();			/* find and mark sources */

	setambient();			/* initialize ambient calculation */

#ifdef  PERSIST
	if (persist) {
		fflush(stdout);
		if (outfile == NULL) {		/* reconnect stdout */
			dup2(duped1, fileno(stdout));
			close(duped1);
		}
		if (persist == PARALLEL) {	/* multiprocessing */
			preload_objs();		/* preload scene */
			shm_boundary = (char *)malloc(16);
			strcpy(shm_boundary, "SHM_BOUNDARY");
			while ((rval=fork()) == 0) {	/* keep on forkin' */
				pflock(1);
				pfhold();
				tstart = time((time_t *)NULL);
			}
			if (rval < 0)
				error(SYSTEM, "cannot fork child for persist function");
			pfdetach();		/* parent exits */
		}
	}
runagain:
	if (persist)
		if (outfile == NULL)			/* if out to stdout */
			dupheader();			/* send header */
		else					/* if out to file */
			duped1 = dup(fileno(stdout));	/* hang onto pipe */
#endif
					/* batch render picture(s) */
	rpict(seqstart, outfile, zfile, recover);
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
		if (outfile != NULL)
			close(duped1);		/* release output handle */
		pfhold();
		tstart = time((time_t *)NULL);	/* reinitialize */
		raynum = nrays = 0;
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

#ifdef SIGALRM /* XXX how critical is this? */
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

	printf("-vt%c\t\t\t\t# view type %s\n", ourview.type,
			ourview.type==VT_PER ? "perspective" :
			ourview.type==VT_PAR ? "parallel" :
			ourview.type==VT_HEM ? "hemispherical" :
			ourview.type==VT_ANG ? "angular" :
			ourview.type==VT_CYL ? "cylindrical" :
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
	printf("-pa %f\t\t\t# pixel aspect ratio\n", pixaspect);
	printf("-pj %f\t\t\t# pixel jitter\n", dstrpix);
	printf("-pm %f\t\t\t# pixel motion\n", mblur);
	printf("-ps %-9d\t\t\t# pixel sample\n", psample);
	printf("-pt %f\t\t\t# pixel threshold\n", maxdiff);
	printf("-t  %-9d\t\t\t# time between reports\n", ralrm);
	printf(erract[WARNING].pf != NULL ?
			"-w+\t\t\t\t# warning messages on\n" :
			"-w-\t\t\t\t# warning messages off\n");
	print_rdefaults();
}
