#ifndef lint
static const char	RCSid[] = "$Id: rpiece.c,v 2.51 2010/05/17 17:03:53 greg Exp $";
#endif
/*
 * Generate sections of a picture.
 */
 

#include <stdio.h>
#include <signal.h>
#include <sys/types.h>

#include "platform.h"
#ifndef NON_POSIX /* XXX need abstraction for process management */
 #include <sys/wait.h>
#endif

#include "standard.h"
#include "color.h"
#include "view.h"
#include "rtprocess.h"

#ifndef F_SETLKW

int
main(
	int argc,
	char *argv[]
)
{
	fprintf(stderr, "%s: no NFS lock manager on this machine\n", argv[0]);
	exit(1);
}

#else

#ifndef NFS
#define  NFS			1
#endif
				/* set the following to 0 to forgo forking */
#ifndef MAXFORK
#if NFS
#define  MAXFORK		3	/* allotment of duped processes */
#else
#define  MAXFORK		0
#endif
#endif
					/* protection from SYSV signals(!) */
#if defined(sgi)
#define guard_io()	sighold(SIGALRM)
#define unguard()	sigrelse(SIGALRM)
#endif
#ifndef guard_io
#define guard_io()	
#define unguard()	
#endif

extern char  *strerror();

				/* rpict command */
char  *rpargv[128] = {"rpict", "-S", "1"};
int  rpargc = 3;
FILE  *torp, *fromrp;
COLR  *pbuf;
				/* our view parameters */
VIEW  ourview = STDVIEW;
double  pixaspect = 1.0;
int  hres = 1024, vres = 1024, hmult = 4, vmult = 4;
				/* output file */
char  *outfile = NULL;
int  outfd;
long  scanorig;
FILE  *syncfp = NULL;		/* synchronization file pointer */
int  synclst = F_UNLCK;		/* synchronization file lock status */
int  nforked = 0;

#define  sflock(t)	if ((t)!=synclst) dolock(fileno(syncfp),synclst=t)

char  *progname;
int  verbose = 0;
unsigned  timelim = 0;
int  rvrlim = -1;

int  gotalrm = 0;
void  onalrm(int i) { gotalrm++; }

static void dolock(int  fd, int  ltyp);
static void init(int  ac, char  **av);
static int nextpiece(int	*xp, int	*yp);
static int rvrpiece(int	*xp, int	*yp);
static int cleanup(int  rstat);
static void rpiece(void);
static int putpiece(int	xpos, int	ypos);
static void filerr(char  *t);


int
main(
	int  argc,
	char  *argv[]
)
{
	register int  i, rval;
	
	progname = argv[0];
	for (i = 1; i < argc; i++) {
						/* expand arguments */
		while ((rval = expandarg(&argc, &argv, i)) > 0)
			;
		if (rval < 0) {
			fprintf(stderr, "%s: cannot expand '%s'\n",
					argv[0], argv[i]);
			exit(1);
		}
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'v':
				switch (argv[i][2]) {
				case '\0':	/* verbose option */
					verbose = !verbose;
					continue;
				case 'f':	/* view file */
					if (viewfile(argv[++i], &ourview, NULL) <= 0) {
						fprintf(stderr,
					"%s: not a view file\n", argv[i]);
						exit(1);
					}
					continue;
				default:	/* view option? */
					rval = getviewopt(&ourview, argc-i, argv+i);
					if (rval >= 0) {
						i += rval;
						continue;
					}
					break;
				}
				break;
			case 'p':		/* pixel aspect ratio? */
				if (argv[i][2] == 'm') {
					fprintf(stderr, "%s: -pm unsupported\n",
						argv[0]);
					++i;
					continue;
				}
				if (argv[i][2] != 'a' || argv[i][3])
					break;
				pixaspect = atof(argv[++i]);
				continue;
			case 'T':		/* time limit (hours) */
				if (argv[i][2])
					break;
				timelim = atof(argv[++i])*3600. + .5;
				continue;
			case 'x':		/* overall x resolution */
				if (argv[i][2])
					break;
				hres = atoi(argv[++i]);
				continue;
			case 'y':		/* overall y resolution */
				if (argv[i][2])
					break;
				vres = atoi(argv[++i]);
				continue;
			case 'X':		/* horizontal multiplier */
				if (argv[i][2])
					break;
				hmult = atoi(argv[++i]);
				continue;
			case 'Y':		/* vertical multiplier */
				if (argv[i][2])
					break;
				vmult = atoi(argv[++i]);
				continue;
			case 'R':		/* recover */
				if (argv[i][2])
					break;
				rvrlim = 0;
			/* fall through */
			case 'F':		/* syncronization file */
				if (argv[i][2])
					break;
				if ((syncfp =
		fdopen(open(argv[++i],O_RDWR|O_CREAT,0666),"r+")) == NULL) {
					fprintf(stderr, "%s: cannot open\n",
							argv[i]);
					exit(1);
				}
				continue;
			case 'z':		/* z-file ist verbotten */
				fprintf(stderr, "%s: -z option not allowed\n",
						argv[0]);
				exit(1);
			case 'o':		/* output file */
				if (argv[i][2])
					break;
				outfile = argv[++i];
				continue;
			} else if (i >= argc-1)
				break;
		rpargv[rpargc++] = argv[i];
	}
	if (i >= argc) {
		fprintf(stderr, "%s: missing octree argument\n", argv[0]);
		exit(1);
	}
	if (outfile == NULL) {
		fprintf(stderr, "%s: missing output file\n", argv[0]);
		exit(1);
	}
	init(argc, argv);
	rpiece();
	exit(cleanup(0));
}


static void
dolock(		/* lock or unlock a file */
	int  fd,
	int  ltyp
)
{
	static struct flock  fls;	/* static so initialized to zeroes */

	fls.l_type = ltyp;
	if (fcntl(fd, F_SETLKW, &fls) < 0) {
		fprintf(stderr, "%s: cannot lock/unlock file: %s\n",
				progname, strerror(errno));
		exit(1);
	}
}


static void
init(			/* set up output file and start rpict */
	int  ac,
	char  **av
)
{
	static char  hrbuf[16], vrbuf[16];
	extern char  VersionID[];
	char  *err;
	FILE  *fp;
	int  hr, vr;
	SUBPROC  rpd; /* since we don't close_process(), this can be local */
					/* set up view */
	if ((err = setview(&ourview)) != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	if (syncfp != NULL) {
		sflock(F_RDLCK);
		fscanf(syncfp, "%d %d", &hmult, &vmult);
		sflock(F_UNLCK);
	}
					/* compute piece size */
	hres /= hmult;
	vres /= vmult;
	if (hres <= 0 || vres <= 0) {
		fprintf(stderr, "%s: illegal resolution/subdivision\n", progname);
		exit(1);
	}
	normaspect(viewaspect(&ourview)*hmult/vmult, &pixaspect, &hres, &vres);
	sprintf(hrbuf, "%d", hres);
	rpargv[rpargc++] = "-x"; rpargv[rpargc++] = hrbuf;
	sprintf(vrbuf, "%d", vres);
	rpargv[rpargc++] = "-y"; rpargv[rpargc++] = vrbuf;
	rpargv[rpargc++] = "-pa"; rpargv[rpargc++] = "0";
	rpargv[rpargc++] = av[ac-1];
	rpargv[rpargc] = NULL;
					/* open output file */
	if ((outfd = open(outfile, O_WRONLY|O_CREAT|O_EXCL, 0666)) >= 0) {
		dolock(outfd, F_WRLCK);
		if ((fp = fdopen(dup(outfd), "w")) == NULL)
			goto filerr;
		newheader("RADIANCE", fp);	/* create header */
		printargs(ac, av, fp);
		fprintf(fp, "SOFTWARE= %s\n", VersionID);
		fputs(VIEWSTR, fp);
		fprintview(&ourview, fp);
		putc('\n', fp);
		if (pixaspect < .99 || pixaspect > 1.01)
			fputaspect(pixaspect, fp);
		fputformat(COLRFMT, fp);
		putc('\n', fp);
		fprtresolu(hres*hmult, vres*vmult, fp);
	} else if ((outfd = open(outfile, O_RDWR)) >= 0) {
		dolock(outfd, F_RDLCK);
		if ((fp = fdopen(dup(outfd), "r+")) == NULL)
			goto filerr;
		getheader(fp, NULL, NULL);	/* skip header */
		if (!fscnresolu(&hr, &vr, fp) ||	/* check resolution */
				hr != hres*hmult || vr != vres*vmult) {
			fprintf(stderr, "%s: resolution mismatch on file \"%s\"\n",
					progname, outfile);
			exit(1);
		}
	} else {
		fprintf(stderr, "%s: cannot open file \"%s\"\n",
				progname, outfile);
		exit(1);
	}
	scanorig = ftell(fp);		/* record position of first scanline */
	if (fclose(fp) == -1)		/* done with stream i/o */
		goto filerr;
	dolock(outfd, F_UNLCK);
					/* start rpict process */
	if (open_process(&rpd, rpargv) <= 0) {
		fprintf(stderr, "%s: cannot start %s\n", progname, rpargv[0]);
		exit(1);
	}
	if ((fromrp = fdopen(rpd.r, "r")) == NULL ||
			(torp = fdopen(rpd.w, "w")) == NULL) {
		fprintf(stderr, "%s: cannot open stream to %s\n",
				progname, rpargv[0]);
		exit(1);
	}
	if ((pbuf = (COLR *)bmalloc(hres*vres*sizeof(COLR))) == NULL) {
		fprintf(stderr, "%s: out of memory\n", progname);
		exit(1);
	}
	signal(SIGALRM, onalrm);
	if (timelim)
		alarm(timelim);
	return;
filerr:
	fprintf(stderr, "%s: i/o error on file \"%s\"\n", progname, outfile);
	exit(1);
}


static int
nextpiece(		/* get next piece assignment */
	int	*xp,
	int	*yp
)
{
	if (gotalrm)			/* someone wants us to quit */
		return(0);
	if (syncfp != NULL) {		/* use sync file */
		/*
		 * So we don't necessarily have to lock and unlock the file
		 * multiple times (very slow), we establish an exclusive
		 * lock at the beginning on our synchronization file and
		 * maintain it in the subroutine rvrpiece().
		 */
		sflock(F_WRLCK);
		fseek(syncfp, 0L, 0);		/* read position */
		if (fscanf(syncfp, "%*d %*d %d %d", xp, yp) < 2) {
			*xp = hmult-1;
			*yp = vmult;
		}
		if (rvrlim == 0)		/* initialize recovery limit */
			rvrlim = *xp*vmult + *yp;
		if (rvrpiece(xp, yp)) {		/* do stragglers first */
			sflock(F_UNLCK);
			return(1);
		}
		if (--(*yp) < 0) {		/* decrement position */
			*yp = vmult-1;
			if (--(*xp) < 0) {	/* all done */
				sflock(F_UNLCK);
				return(0);
			}
		}
		fseek(syncfp, 0L, 0);		/* write new position */
		fprintf(syncfp, "%4d %4d\n%4d %4d\n\n", hmult, vmult, *xp, *yp);
		fflush(syncfp);
		sflock(F_UNLCK);		/* release sync file */
		return(1);
	}
	return(scanf("%d %d", xp, yp) == 2);	/* use stdin */
}


static int
rvrpiece(		/* check for recoverable pieces */
	register int	*xp,
	register int	*yp
)
{
	static char  *pdone = NULL;	/* which pieces are done */
	static long  readpos = -1;	/* how far we've read */
	register int  i;
	/*
	 * This routine is called by nextpiece() with an
	 * exclusive lock on syncfp and the file pointer at the
	 * appropriate position to read in the finished pieces.
	 */
	if (rvrlim < 0)
		return(0);		/* only check if asked */
	if (pdone == NULL)		/* first call */
		pdone = calloc(hmult*vmult, sizeof(char));
	if (pdone == NULL) {
		fprintf(stderr, "%s: out of memory\n", progname);
		exit(1);
	}
	if (readpos != -1)		/* mark what's been done */
		fseek(syncfp, readpos, 0);
	while (fscanf(syncfp, "%d %d", xp, yp) == 2)
		pdone[*xp*vmult+*yp] = 1;
	if (!feof(syncfp)) {
		fprintf(stderr, "%s: format error in sync file\n", progname);
		exit(1);
	}
	readpos = ftell(syncfp);
	i = hmult*vmult;		/* find an unaccounted for piece */
	while (i-- > rvrlim)
		if (!pdone[i]) {
			*xp = i / vmult;
			*yp = i % vmult;
			pdone[i] = 1;	/* consider it done */
			return(1);
		}
	rvrlim = -1;			/* nothing left to recover */
	free(pdone);
	pdone = NULL;
	return(0);
}


static int
cleanup(			/* close rpict process and clean up */
	int  rstat
)
{
	int  status;

	bfree((char *)pbuf, hres*vres*sizeof(COLR));
	fclose(torp);
	fclose(fromrp);
	while (wait(&status) != -1)
		if (rstat == 0)
			rstat = status>>8 & 0xff;
	return(rstat);
}


static void
rpiece(void)			/* render picture piece by piece */
{
	VIEW  pview;
	int  xorg, yorg;
					/* compute view parameters */
	pview = ourview;
	switch (ourview.type) {
	case VT_PER:
		pview.horiz = (2.*180./PI)*atan(
				tan((PI/180./2.)*ourview.horiz)/hmult );
		pview.vert = (2.*180./PI)*atan(
				tan((PI/180./2.)*ourview.vert)/vmult );
		break;
	case VT_PAR:
	case VT_ANG:
		pview.horiz = ourview.horiz / hmult;
		pview.vert = ourview.vert / vmult;
		break;
	case VT_CYL:
		pview.horiz = ourview.horiz / hmult;
		pview.vert = (2.*180./PI)*atan(
				tan((PI/180./2.)*ourview.vert)/vmult );
		break;
	case VT_HEM:
		pview.horiz = (2.*180./PI)*asin(
				sin((PI/180./2.)*ourview.horiz)/hmult );
		pview.vert = (2.*180./PI)*asin(
				sin((PI/180./2.)*ourview.vert)/vmult );
		break;
	case VT_PLS:
		pview.horiz = sin((PI/180./2.)*ourview.horiz) /
				(1.0 + cos((PI/180./2.)*ourview.horiz)) / hmult;
		pview.horiz *= pview.horiz;
		pview.horiz = (2.*180./PI)*acos((1. - pview.horiz) /
						(1. + pview.horiz));
		pview.vert = sin((PI/180./2.)*ourview.vert) /
				(1.0 + cos((PI/180./2.)*ourview.vert)) / vmult;
		pview.vert *= pview.vert;
		pview.vert = (2.*180./PI)*acos((1. - pview.vert) /
						(1. + pview.vert));
		break;
	default:
		fprintf(stderr, "%s: unknown view type '-vt%c'\n",
				progname, ourview.type);
		exit(cleanup(1));
	}
					/* render each piece */
	while (nextpiece(&xorg, &yorg)) {
		pview.hoff = ourview.hoff*hmult + xorg - 0.5*(hmult-1);
		pview.voff = ourview.voff*vmult + yorg - 0.5*(vmult-1);
		fputs(VIEWSTR, torp);
		fprintview(&pview, torp);
		putc('\n', torp);
		fflush(torp);			/* assigns piece to rpict */
		putpiece(xorg, yorg);		/* place piece in output */
	}
}


static int
putpiece(		/* get next piece from rpict */
int	xpos,
int	ypos
)
{
	struct flock  fls;
	int  pid, status;
	int  hr, vr;
	register int  y;
				/* check bounds */
	if ((xpos < 0) | (ypos < 0) | (xpos >= hmult) | (ypos >= vmult)) {
		fprintf(stderr, "%s: requested piece (%d,%d) out of range\n",
				progname, xpos, ypos);
		exit(cleanup(1));
	}
				/* check header from rpict */
	guard_io();
	getheader(fromrp, NULL, NULL);
	if (!fscnresolu(&hr, &vr, fromrp) || (hr != hres) | (vr != vres)) {
		fprintf(stderr, "%s: resolution mismatch from %s\n",
				progname, rpargv[0]);
		exit(cleanup(1));
	}
	if (verbose) {				/* notify caller */
		printf("%d %d begun\n", xpos, ypos);
		fflush(stdout);
	}
	unguard();
				/* load new piece into buffer */
	for (y = 0; y < vr; y++) {
		guard_io();
		if (freadcolrs(pbuf+y*hr, hr, fromrp) < 0) {
			fprintf(stderr, "%s: read error from %s\n",
					progname, rpargv[0]);
			exit(cleanup(1));
		}
		unguard();
	}
#if MAXFORK
				/* fork so we don't slow rpict down */
	if ((pid = fork()) > 0) {
		if (++nforked >= MAXFORK) {
			wait(&status);		/* reap a child */
			if (status)
				exit(cleanup(status>>8 & 0xff));
			nforked--;
		}
		return(pid);
	}
#else
	pid = -1;		/* no forking */
#endif
	fls.l_start = scanorig +
		((long)(vmult-1-ypos)*vres*hmult+xpos)*hres*sizeof(COLR);
#if NFS
	fls.l_len = ((long)(vres-1)*hmult+1)*hres*sizeof(COLR);
				/* lock file section so NFS doesn't mess up */
	fls.l_whence = 0;
	fls.l_type = F_WRLCK;
	if (fcntl(outfd, F_SETLKW, &fls) < 0)
		filerr("lock");
#endif
				/* write new piece to file */
	if (lseek(outfd, (off_t)fls.l_start, SEEK_SET) < 0)
		filerr("seek");
	if (hmult == 1) {
		if (writebuf(outfd, (char *)pbuf,
				vr*hr*sizeof(COLR)) != vr*hr*sizeof(COLR))
			filerr("write");
	} else
		for (y = 0; y < vr; y++) {
			if (writebuf(outfd, (char *)(pbuf+y*hr),
					hr*sizeof(COLR)) != hr*sizeof(COLR))
				filerr("write");
			if (y < vr-1 && lseek(outfd,
					(off_t)(hmult-1)*hr*sizeof(COLR),
					SEEK_CUR) < 0)
				filerr("seek");
		}
#if NFS
	fls.l_type = F_UNLCK;		/* release lock */
	if (fcntl(outfd, F_SETLKW, &fls) < 0)
		filerr("lock");
#endif
	if (verbose) {				/* notify caller */
		printf("%d %d done\n", xpos, ypos);
		fflush(stdout);
	}
	if (syncfp != NULL) {			/* record what's been done */
		sflock(F_WRLCK);
		fseek(syncfp, 0L, 2);		/* append index */
		fprintf(syncfp, "%4d %4d\n", xpos, ypos);
		fflush(syncfp);
				/*** Unlock not necessary, since
		sflock(F_UNLCK);	_exit() or nextpiece() is next ***/
	}
	if (pid == -1)		/* didn't fork or fork failed */
		return(0);
	_exit(0);		/* else exit child process (releasing locks) */
}


static void
filerr(			/* report file error and exit */
	char  *t
)
{
	fprintf(stderr, "%s: %s error on file \"%s\": %s\n",
			progname, t, outfile, strerror(errno));
	_exit(1);
}

#endif
