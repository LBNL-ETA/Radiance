#ifndef lint
static const char	RCSid[] = "$Id: rholo.c,v 3.71 2005/12/16 21:43:07 greg Exp $";
#endif
/*
 * Radiance holodeck generation controller
 */

#include <stdio.h>
#include <time.h>
#include <signal.h>
#include <sys/stat.h>
#include <string.h>

#include "platform.h"
#include "rterror.h"
#include "resolu.h"
#include "rholo.h"
#include "random.h"

#ifndef FRAGWARN
#define FRAGWARN	20		/* fragmentation for warning (%) */
#endif
#ifndef MAXQTIME
#define MAXQTIME	5		/* target maximum seconds in queue */
#endif
					/* manual cache flushing frequency */
#ifndef RTFLUSH
#if MAXQTIME
#define RTFLUSH		(300/MAXQTIME*totqlen)	/* <= 5 minutes */
#else
#define RTFLUSH		(50*totqlen)		/* just guess */
#endif
#endif
			/* the following must be consistent with rholo.h */
int	NVARS = NRHVARS;		/* total number of variables */

VARIABLE	vv[] = RHVINIT;		/* variable-value pairs */

char	*progname;		/* our program name */
char	*hdkfile;		/* holodeck file name */
char	froot[256];		/* root file name */

int	ncprocs = 0;		/* desired number of compute processes */

char	*outdev = NULL;		/* output device name */

int	readinp = 0;		/* read commands from stdin */

int	force = 0;		/* allow overwrite of holodeck (-1 == read-only) */

time_t	starttime;		/* time we got started */
time_t	endtime;		/* time we should end by */
time_t	reporttime;		/* time for next report */

off_t	maxdisk;		/* maximum file space (bytes) */

int	rtargc = 1;		/* rtrace command */
char	*rtargv[128] = {"rtrace", NULL};

int	orig_mode = -1;		/* original file mode (-1 if unchanged) */

long	nraysdone = 0L;		/* number of rays done */
long	npacksdone = 0L;	/* number of packets done */

PACKET	*freepacks;		/* available packets */
int	totqlen;		/* maximum queue length (number of packets) */

char  *sigerr[NSIG];		/* signal error messages */

extern int	nowarn;		/* turn warnings off? */

static void onsig(int  signo);
static void sigdie(int  signo, char  *msg);
static int resfmode(int	fd, int	mod);
static void initrholo(void);
static int rholo(void);
static void setdefaults(HDGRID	*gp);
static void creatholo(HDGRID	*gp);
static gethfunc headline;
static void loadholo(void);
static void rootname(char	*rn, char	*fn);
static void badvalue(int	vc);


int
main(
	int	argc,
	char	*argv[]
)
{
	int	i;

	progname = argv[0];			/* get arguments */
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'w':			/* turn off warnings */
			nowarn++;
			break;
		case 'f':			/* force overwrite */
			force = 1;
			break;
		case 'r':			/* read-only mode */
			force = -1;
			break;
		case 'i':			/* read input from stdin */
			readinp++;
			break;
		case 'n':			/* compute processes */
			if (i >= argc-2)
				goto userr;
			ncprocs = atoi(argv[++i]);
			break;
		case 'o':			/* output display */
			if (i >= argc-2)
				goto userr;
			outdev = argv[++i];
			break;
		default:
			goto userr;
		}
						/* get root file name */
	if (i >= argc)
		goto userr;
	rootname(froot, hdkfile=argv[i++]);
						/* load variables? */
	if (i < argc)
		if (argv[i][0] != '-' && argv[i][0] != '+')
			loadvars(argv[i]);	/* load variables from file */

	if (i >= argc || argv[i][0] == '+')
		loadholo();			/* load existing holodeck */

	while (++i < argc)			/* get command line settings */
		if (setvariable(argv[i], matchvar) < 0) {
			sprintf(errmsg, "unknown variable: %s", argv[i]);
			error(USER, errmsg);
		}
						/* check settings */
	checkvalues();
						/* load rad input file */
	getradfile();

	if (hdlist[0] == NULL) {		/* create new holodeck */
		HDGRID	hdg[HDMAX];
							/* set defaults */
		setdefaults(hdg);
							/* check read-only */
		if (force < 0)
			error(USER, "cannot create read-only holodeck");
							/* holodeck exists? */
		if (!force && access(hdkfile, R_OK|W_OK) == 0)
			error(USER,
				"holodeck file exists -- use -f to overwrite");
							/* create holodeck */
		creatholo(hdg);
	} else					/* else just set defaults */
		setdefaults(NULL);
						/* initialize */
	initrholo();
						/* main loop */
	while (rholo())
		;
						/* done */
	quit(0);
userr:
	fprintf(stderr,
"Usage: %s [-n nprocs][-o disp][-w][-r|-f] output.hdk [control.hif|+|- [VAR=val ..]]\n",
			progname);
	quit(1);
	return 1; /* pro forma return */
}


static void
onsig(				/* fatal signal */
	int  signo
)
{
	static int  gotsig = 0;

	if (gotsig > 1)			/* we're going as fast as we can! */
		return;
	if (gotsig++) {			/* two signals and we split */
		hdsync(NULL, 0);	/* don't leave w/o saying goodbye */
		_exit(signo);
	}
	alarm(300);			/* allow 5 minutes to clean up */
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


static int
resfmode(		/* restrict open file access mode */
	int	fd,
	int	mod
)
{
	struct stat	stbuf;
					/* get original mode */
	if (fstat(fd, &stbuf) < 0)
		error(SYSTEM, "cannot stat open holodeck file");
	mod &= stbuf.st_mode;		/* always more restrictive */
	if (mod == (stbuf.st_mode & 0777))
		return(-1);		/* already set */
					/* else change it */
	if (fchmod(fd, mod) < 0) {
		error(WARNING, "cannot change holodeck file access mode");
		return(-1);
	}
	return(stbuf.st_mode);		/* return original mode */
}


static void
initrholo(void)			/* get our holodeck running */
{
	extern int	global_packet();
	register int	i;
						/* close holodeck on exec() */
	fcntl(hdlist[0]->fd, F_SETFD, FD_CLOEXEC);

	if (outdev != NULL)			/* open output device */
		disp_open(outdev);
	else if (ncprocs > 0)			/* else use global ray feed */
		init_global();
						/* record disk space limit */
	if (!vdef(DISKSPACE))
		maxdisk = ((off_t)1<<(sizeof(off_t)*8-2)) - 1024;
	else
		maxdisk = 1024.*1024.*vflt(DISKSPACE);
						/* set up memory cache */
	if (outdev == NULL)
		hdcachesize = 0;		/* manual flushing */
	else if (vdef(CACHE))
		hdcachesize = 1024.*1024.*vflt(CACHE);
						/* open report file */
	if (vdef(REPORT)) {
		register char	*s = sskip2(vval(REPORT), 1);
		if (*s && freopen(s, "a", stderr) == NULL)
			quit(2);
	}
						/* mark the starting time */
	starttime = time(NULL);
						/* compute end time */
	if (!vdef(TIME) || vflt(TIME) <= FTINY)
		endtime = 0;
	else
		endtime = starttime + vflt(TIME)*3600. + .5;
						/* start rtrace */
	if (ncprocs > 0) {
		totqlen = i = start_rtrace();
		if (i < 1)
			error(USER, "cannot start rtrace process(es)");
		if (vdef(REPORT)) {		/* make first report */
			printargs(rtargc, rtargv, stderr);
			report(0);
		}
						/* allocate packets */
		freepacks = (PACKET *)bmalloc(i*sizeof(PACKET));
		if (freepacks == NULL)
			goto memerr;
		freepacks[--i].nr = 0;
		freepacks[i].next = NULL;
		if (!vdef(OBSTRUCTIONS) || !vbool(OBSTRUCTIONS)) {
			freepacks[i].offset = (float *)bmalloc(
					RPACKSIZ*sizeof(float)*(i+1) );
			if (freepacks[i].offset == NULL)
				goto memerr;
		} else
			freepacks[i].offset = NULL;
		while (i--) {
			freepacks[i].nr = 0;
			freepacks[i].offset = freepacks[i+1].offset == NULL ?
					NULL : freepacks[i+1].offset+RPACKSIZ ;
			freepacks[i].next = &freepacks[i+1];
		}
	}
					/* set up signal handling */
	sigdie(SIGINT, "Interrupt");
	sigdie(SIGTERM, "Terminate");
#ifdef SIGHUP
	sigdie(SIGHUP, "Hangup");
#endif
#ifdef SIGPIPE
	sigdie(SIGPIPE, "Broken pipe");
#endif
#ifdef SIGALRM
	sigdie(SIGALRM, "Alarm clock");
#endif
#ifdef	SIGXCPU
	sigdie(SIGXCPU, "CPU limit exceeded");
#endif
#ifdef	SIGXFSZ
	sigdie(SIGXFSZ, "File size exceeded");
#endif
					/* protect holodeck file */
	orig_mode = resfmode(hdlist[0]->fd, (ncprocs>0) & (force>=0) ? 0 : 0444);
	return;
memerr:
	error(SYSTEM, "out of memory in initrholo");
}


static int
rholo(void)				/* holodeck main loop */
{
	static off_t	nextfragwarn = 100L<<20;
	static int	idle = 0;
	PACKET	*pl = NULL, *plend;
	off_t	fsiz;
	int	pksiz;
	register PACKET	*p;
	time_t	t;
					/* check display */
	if (nprocs <= 0)
		idle = 1;
	if (outdev != NULL) {
		if (!disp_check(idle))
			return(0);	/* quit request */
		if (nprocs <= 0)
			return(1);
	} else if (idle)
		return(0);		/* all done */
	fsiz = hdfilen(hdlist[0]->fd);	/* check file size */
	if (maxdisk > 0 && fsiz >= maxdisk) {
		error(USER, "file limit exceeded");
		done_rtrace();
		return(1);	/* comes back */
	}
#if FRAGWARN
	if (fsiz >= nextfragwarn) {
		double	pctfrag = 100.*(fsiz-hdfiluse(hdlist[0]->fd))/fsiz;
	       	if (pctfrag >= (double)FRAGWARN) {
			sprintf(errmsg, "holodeck file fragmentation is %.0f%%",
					pctfrag);
			error(WARNING, errmsg);
			nextfragwarn = fsiz + (fsiz>>2);
		} else
			nextfragwarn = fsiz + (10L<<20);
	}
#endif
	t = time(NULL);			/* check time */
	if (endtime > 0 && t >= endtime) {
		error(USER, "time limit exceeded");
		done_rtrace();
		return(1);	/* comes back */
	}
	if (reporttime > 0 && t >= reporttime)
		report(t);
					/* figure out good packet size */
	pksiz = RPACKSIZ;
#if MAXQTIME
	if (!chunkycmp) {
		pksiz = nraysdone*MAXQTIME/(totqlen*(t - starttime + 1L));
		if (pksiz < 1) pksiz = 1;
		else if (pksiz > RPACKSIZ) pksiz = RPACKSIZ;
	}
#endif
	idle = 0;			/* get packets to process */
	while (freepacks != NULL) {
		p = freepacks; freepacks = p->next; p->next = NULL;
		if (!next_packet(p, pksiz)) {
			p->next = freepacks; freepacks = p;
			idle = 1;
			break;
		}
		if (pl == NULL) pl = p;
		else plend->next = p;
		plend = p;
	}
					/* process packets */
	done_packets(do_packets(pl));
	return(1);			/* and continue */
}


static void
setdefaults(			/* set default values */
	register HDGRID	*gp
)
{
	extern char	*atos();
	register int	i;
	int	n;
	double	len[3], d;

	if (!vdef(SECTION)) {
		sprintf(errmsg, "%s must be defined", vnam(SECTION));
		error(USER, errmsg);
	}
	if (!vdef(OCTREE)) {
		if ((vval(OCTREE) = bmalloc(strlen(froot)+5)) == NULL)
			error(SYSTEM, "out of memory");
		sprintf(vval(OCTREE), "%s.oct", froot);
		vdef(OCTREE)++;
	}
	if (!vdef(VDIST)) {
		vval(VDIST) = "F";
		vdef(VDIST)++;
	}
				/* append rendering options */
	if (vdef(RENDER))
		rtargc += wordstring(rtargv+rtargc, vval(RENDER));
	
	if (gp == NULL)		/* already initialized? */
		return;
				/* set grid parameters */
	for (n = 0; n < vdef(SECTION); n++, gp++) {
		gp->grid[0] = gp->grid[1] = gp->grid[2] = 0;
		if (sscanf(nvalue(SECTION, n),
		"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %hd %hd %hd",
				&gp->orig[0], &gp->orig[1], &gp->orig[2],
				&gp->xv[0][0], &gp->xv[0][1], &gp->xv[0][2],
				&gp->xv[1][0], &gp->xv[1][1], &gp->xv[1][2],
				&gp->xv[2][0], &gp->xv[2][1], &gp->xv[2][2],
				&gp->grid[0], &gp->grid[1], &gp->grid[2]) < 12)
			badvalue(SECTION);
		for (i = 0; i < 3; i++)
			len[i] = VLEN(gp->xv[i]);
		if (!vdef(GRID)) {
			d = 2/2e5*( len[0]*len[0]*(len[1]*len[1] +
					len[2]*len[2] + 4*len[1]*len[2])
				+ len[1]*len[1]*len[2]*(len[2] + 4*len[0])
				+ 4*len[0]*len[1]*len[2]*len[2] );
			d = sqrt(sqrt(d));
		} else if ((d = vflt(GRID)) <= FTINY)
			badvalue(GRID);
		for (i = 0; i < 3; i++)
			if (gp->grid[i] <= 0)
				gp->grid[i] = len[i]/d + (1.-FTINY);
	}
}


static void
creatholo(			/* create a holodeck output file */
	HDGRID	*gp
)
{
	extern char	VersionID[];
	int32	lastloc, nextloc;
	int	n;
	int	fd;
	FILE	*fp;
					/* open & truncate file */
	if ((fp = fopen(hdkfile, "w+")) == NULL) {
		sprintf(errmsg, "cannot open \"%s\" for writing", hdkfile);
		error(SYSTEM, errmsg);
	}
					/* write information header */
	newheader("RADIANCE", fp);
	fprintf(fp, "SOFTWARE= %s\n", VersionID);
	printvars(fp);
	fputformat(HOLOFMT, fp);
	fputc('\n', fp);
	putw(HOLOMAGIC, fp);		/* put magic number */
	fd = dup(fileno(fp));
	fclose(fp);			/* flush and close stdio stream */
	lastloc = lseek(fd, (off_t)0, SEEK_END);
	for (n = vdef(SECTION); n--; gp++) {	/* initialize each section */
		nextloc = 0L;
		write(fd, (char *)&nextloc, sizeof(nextloc));
		hdinit(fd, gp);			/* writes beam index */
		if (!n)
			break;
		nextloc = hdfilen(fd);		/* write section pointer */
		if (lseek(fd, (off_t)lastloc, SEEK_SET) < 0)
			error(SYSTEM,
				"cannot seek on holodeck file in creatholo");
		write(fd, (char *)&nextloc, sizeof(nextloc));
		lseek(fd, (off_t)(lastloc=nextloc), SEEK_SET);
	}
}


static int
headline(			/* process information header line */
	char	*s,
	void	*p
)
{
	extern char	FMTSTR[];
	register char	*cp;
	char	fmt[32];

	if (formatval(fmt, s)) {
		if (strcmp(fmt, HOLOFMT)) {
			sprintf(errmsg, "%s file \"%s\" has %s%s",
					HOLOFMT, hdkfile, FMTSTR, fmt);
			error(USER, errmsg);
		}
		return(0);
	}
	for (cp = s; *cp; cp++)		/* take off any comments */
		if (*cp == '#') {
			*cp = '\0';
			break;
		}
	setvariable(s, matchvar);	/* don't flag errors */
	return(0);
}


static void
loadholo(void)			/* start loading a holodeck from fname */
{
	FILE	*fp;
	int	fd;
	int	n;
	int32	nextloc;
	
	if ((ncprocs > 0) & (force >= 0))
		fp = fopen(hdkfile, "r+");
	else
		fp = NULL;
	if (fp == NULL) {
		if ((fp = fopen(hdkfile, "r")) == NULL) {
			sprintf(errmsg, "cannot open \"%s\"", hdkfile);
			error(SYSTEM, errmsg);
		}
		if (ncprocs > 0) {
			sprintf(errmsg,
			"\"%s\" opened read-only; new rays will be discarded",
					hdkfile);
			error(WARNING, errmsg);
			force = -1;
		}
	}
					/* load variables from header */
	getheader(fp, headline, NULL);
					/* check magic number */
	if (getw(fp) != HOLOMAGIC) {
		sprintf(errmsg, "bad magic number in holodeck file \"%s\"",
				hdkfile);
		error(USER, errmsg);
	}
	nextloc = ftell(fp);			/* get stdio position */
	fd = dup(fileno(fp));
	fclose(fp);				/* done with stdio */
	for (n = 0; nextloc > 0L; n++) {	/* initialize each section */
		lseek(fd, (off_t)nextloc, SEEK_SET);
		read(fd, (char *)&nextloc, sizeof(nextloc));
		hdinit(fd, NULL);
	}
	if (n != vdef(SECTION)) {
		sprintf(errmsg, "number of sections does not match %s setting",
				vnam(SECTION));
		error(WARNING, errmsg);
	}
}


extern void
done_packets(		/* handle finished packets */
	PACKET	*pl
)
{
	static int	n2flush = 0;
	register PACKET	*p;

	while (pl != NULL) {
		p = pl; pl = p->next; p->next = NULL;
		if (p->nr > 0) {		/* add to holodeck */
			memcpy( (void *)hdnewrays(hdlist[p->hd],p->bi,p->nr),
				(void *)p->ra,
				p->nr*sizeof(RAYVAL));
			if (outdev != NULL)	/* display it */
				disp_packet((PACKHEAD *)p);
			if (hdcachesize <= 0)
				n2flush++;
			nraysdone += p->nr;
			npacksdone++;
			p->nr = 0;
		}
		p->next = freepacks;		/* push onto free list */
		freepacks = p;
	}
	if (n2flush >= RTFLUSH) {
		if (outdev != NULL)
			hdsync(NULL, 1);
		else
			hdflush(NULL);
		n2flush = 0;
	}
}


static void
rootname(		/* remove tail from end of fn */
	register char	*rn,
	register char	*fn
)
{
	char	*tp, *dp;

	for (tp = NULL, dp = rn; (*rn = *fn++); rn++) {
		if (*rn == '/')
			dp = rn;
		else if (*rn == '.')
			tp = rn;
	}
	if (tp != NULL && tp > dp)
		*tp = '\0';
}


static void
badvalue(			/* report bad variable value and exit */
	int	vc
)
{
	sprintf(errmsg, "bad value for variable '%s'", vnam(vc));
	error(USER, errmsg);
}


void
eputs(s)			/* put error message to stderr */
register char  *s;
{
	static int  midline = 0;

	if (!*s)
		return;
	if (!midline++) {	/* prepend line with program name */
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
quit(ec)			/* exit program gracefully */
int	ec;
{
	int	status = 0;

	if (hdlist[0] != NULL) {	/* close holodeck */
		if (nprocs > 0)
			status = done_rtrace();		/* calls hdsync() */
		if ((ncprocs > 0) & (force >= 0) && vdef(REPORT)) {
			off_t	fsiz, fuse;
			fsiz = hdfilen(hdlist[0]->fd);
			fuse = hdfiluse(hdlist[0]->fd);
			fprintf(stderr,
			"%s: %.1f Mbyte holodeck file, %.1f%% fragmentation\n",
					hdkfile, fsiz/(1024.*1024.),
					100.*(fsiz-fuse)/fsiz);
		}
	}
	if (orig_mode >= 0)		/* reset holodeck access mode */
		fchmod(hdlist[0]->fd, orig_mode);
	if (outdev != NULL)		/* close display */
		disp_close();
	exit(ec ? ec : status);		/* exit */
}
