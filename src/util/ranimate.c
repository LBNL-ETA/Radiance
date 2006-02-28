#ifndef lint
static const char RCSid[] = "$Id: ranimate.c,v 2.51 2006/02/28 03:44:54 greg Exp $";
#endif
/*
 * Radiance animation control program
 *
 * The main difference between this program and ranimove is that
 * we have many optimizations here for camera motion in static
 * environments, calling rpict and pinterp on multiple processors,
 * where ranimove puts its emphasis on object motion, and does
 * not use any external programs for image generation.
 *
 * See the ranimate(1) man page for further details.
 */

#include "copyright.h"

#include <stdio.h>
#include <ctype.h>
#include <sys/stat.h>
#include <time.h>
#include <signal.h>

#include "platform.h"
#include "rtprocess.h"
#include "paths.h"
#include "standard.h"
#include "view.h"
#include "vars.h"
#include "netproc.h"
				/* default blur samples */
#ifndef DEF_NBLUR
#define DEF_NBLUR	5
#endif
				/* default remote shell */
#ifndef REMSH
#define REMSH		"ssh"
#endif
				/* input variables (alphabetical by name) */
#define ANIMATE		0		/* animation command */
#define ARCHIVE		1		/* archiving command */
#define BASENAME	2		/* output image base name */
#define DBLUR		3		/* depth of field blur */
#define DIRECTORY	4		/* working (sub)directory */
#define DISKSPACE	5		/* how much disk space to use */
#define END		6		/* ending frame number */
#define EXPOSURE	7		/* how to compute exposure */
#define HOST		8		/* rendering host machine */
#define INTERP		9		/* # frames to interpolate */
#define MBLUR		10		/* motion blur parameters */
#define NEXTANIM	11		/* next animation file */
#define OCTREE		12		/* octree file name */
#define OVERSAMP	13		/* # times to oversample image */
#define PFILT		14		/* pfilt options */
#define PINTERP		15		/* pinterp options */
#define RENDER		16		/* rendering options */
#define RESOLUTION	17		/* desired final resolution */
#define RIF		18		/* rad input file */
#define RSH		19		/* remote shell script or program */
#define RTRACE		20		/* use rtrace with pinterp? */
#define START		21		/* starting frame number */
#define TRANSFER	22		/* frame transfer command */
#define VIEWFILE	23		/* animation frame views */

int	NVARS = 24;		/* total number of variables */

VARIABLE	vv[] = {		/* variable-value pairs */
	{"ANIMATE",	2,	0,	NULL,	onevalue},
	{"ARCHIVE",	2,	0,	NULL,	onevalue},
	{"BASENAME",	3,	0,	NULL,	onevalue},
	{"DBLUR",	2,	0,	NULL,	onevalue},
	{"DIRECTORY",	3,	0,	NULL,	onevalue},
	{"DISKSPACE",	3,	0,	NULL,	fltvalue},
	{"END",		2,	0,	NULL,	intvalue},
	{"EXPOSURE",	3,	0,	NULL,	onevalue},
	{"host",	4,	0,	NULL,	NULL},
	{"INTERPOLATE",	3,	0,	NULL,	intvalue},
	{"MBLUR",	2,	0,	NULL,	onevalue},
	{"NEXTANIM",	3,	0,	NULL,	onevalue},
	{"OCTREE",	3,	0,	NULL,	onevalue},
	{"OVERSAMPLE",	2,	0,	NULL,	fltvalue},
	{"pfilt",	2,	0,	NULL,	catvalues},
	{"pinterp",	2,	0,	NULL,	catvalues},
	{"render",	3,	0,	NULL,	catvalues},
	{"RESOLUTION",	3,	0,	NULL,	onevalue},
	{"RIF",		3,	0,	NULL,	onevalue},
	{"RSH",		3,	0,	NULL,	onevalue},
	{"RTRACE",	2,	0,	NULL,	boolvalue},
	{"START",	2,	0,	NULL,	intvalue},
	{"TRANSFER",	2,	0,	NULL,	onevalue},
	{"VIEWFILE",	2,	0,	NULL,	onevalue},
};

#define SFNAME	"STATUS"		/* status file name */

struct {
	char	host[64];		/* control host name */
	RT_PID	pid;			/* control process id */
	char	cfname[128];		/* control file name */
	int	rnext;			/* next frame to render */
	int	fnext;			/* next frame to filter */
	int	tnext;			/* next frame to transfer */
}	astat;			/* animation status */

char	*progname;		/* our program name */
char	*cfname;		/* our control file name */

int	nowarn = 0;		/* turn warnings off? */
int	silent = 0;		/* silent mode? */
int	noaction = 0;		/* take no action? */

char	*remsh;			/* remote shell program/script */
char	rendopt[2048];		/* rendering options */
char	rresopt[32];		/* rendering resolution options */
char	fresopt[32];		/* filter resolution options */
int	pfiltalways;		/* always use pfilt? */

char	arcargs[10240];		/* files to archive */
char	*arcfirst, *arcnext;	/* pointers to first and next argument */

struct pslot {
	RT_PID	pid;			/* process ID (0 if empty) */
	int	fout;			/* output frame number */
	int	(*rcvf)();		/* recover function */
}	*pslot;			/* process slots */
int	npslots;		/* number of process slots */

#define phostname(ps)	((ps)->hostname[0] ? (ps)->hostname : astat.host)
PSERVER	*lastpserver;		/* last process server with error */

static struct pslot * findpslot(RT_PID pid);
static void checkdir(void);
static VIEW * getview(int n);

static char * dirfile(char *df, register char *path);
static char * getexp(int n);
static int getblur(double *mbf, double *dbf);
static int getastat(void);
static void getradfile(char *rfargs);
static void badvalue(int vc);
static int rmfile(char *fn);
static int runcom(char *cs);
static int pruncom(char *com, char *ppins, int maxcopies);
static void bwait(int ncoms);
static RT_PID bruncom(char *com, int fout, int (*rf)());
static int serverdown(void);
static pscompfunc donecom;
static int countviews(void);
static int dofilt(int frame, int rvr);
static void archive(void);
static int frecover(int frame);
static int recover(int frame);
static void sethosts(void);
static void walkwait(int first, int last, char *vfn);
static void animrend(int frame, VIEW *vp);
static void transferframes(void);
static void filterframes(void);
static void renderframes(int nframes);
static void animate(void);
static void setdefaults(void);
static void putastat(void);


int
main(
	int	argc,
	char	*argv[]
)
{
	int	explicate = 0;
	int	i;

	progname = argv[0];			/* get arguments */
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'e':			/* print variables */
			explicate++;
			break;
		case 'w':			/* turn off warnings */
			nowarn++;
			break;
		case 's':			/* silent mode */
			silent++;
			break;
		case 'n':			/* take no action */
			noaction++;
			break;
		default:
			goto userr;
		}
	if (i != argc-1)
		goto userr;
	cfname = argv[i];
						/* load variables */
	loadvars(cfname);
						/* check variables */
	checkvalues();
						/* did we get DIRECTORY? */
	checkdir();
						/* check status */
	if (getastat() < 0) {
		fprintf(stderr, "%s: exiting\n", progname);
		quit(1);
	}
						/* pfilt always if options given */
	pfiltalways = vdef(PFILT);
						/* load RIF if any */
	if (vdef(RIF))
		getradfile(vval(RIF));
						/* set defaults */
	setdefaults();
						/* print variables */
	if (explicate)
		printvars(stdout);
						/* set up process servers */
	sethosts();
						/* run animation */
	animate();
						/* all done */
	if (vdef(NEXTANIM)) {
		char	*fullp;
		argv[i] = vval(NEXTANIM);	/* just change input file */
		if (!silent)
			printargs(argc, argv, stdout);
		fflush(stdout);
		if ((fullp = getpath(argv[0],getenv("PATH"),X_OK)) == NULL)
			fprintf(stderr, "%s: command not found\n", argv[0]);
		else
			execv(fullp, argv);
		quit(1);
	}
	quit(0);
	return(0); /* pro forma return */
userr:
	fprintf(stderr, "Usage: %s [-s][-n][-w][-e] anim_file\n", progname);
	quit(1);
	return 1; /* pro forma return */
}


static int
getastat(void)			/* check/set animation status */
{
	char	sfname[256];
	FILE	*fp;

	sprintf(sfname, "%s/%s", vval(DIRECTORY), SFNAME);
	if ((fp = fopen(sfname, "r")) == NULL) {
		if (errno != ENOENT) {
			perror(sfname);
			return(-1);
		}
		astat.rnext = astat.fnext = astat.tnext = 0;
		goto setours;
	}
	if (fscanf(fp, "Control host: %s\n", astat.host) != 1)
		goto fmterr;
	if (fscanf(fp, "Control PID: %d\n", &astat.pid) != 1)
		goto fmterr;
	if (fscanf(fp, "Control file: %s\n", astat.cfname) != 1)
		goto fmterr;
	if (fscanf(fp, "Next render: %d\n", &astat.rnext) != 1)
		goto fmterr;
	if (fscanf(fp, "Next filter: %d\n", &astat.fnext) != 1)
		goto fmterr;
	if (fscanf(fp, "Next transfer: %d\n", &astat.tnext) != 1)
		goto fmterr;
	fclose(fp);
	if (astat.pid != 0) {			/* thinks it's still running */
		if (strcmp(myhostname(), astat.host)) {
			fprintf(stderr,
			"%s: process %d may still be running on host %s\n",
					progname, astat.pid, astat.host);
			return(-1);
		}
		if (kill(astat.pid, 0) != -1 || errno != ESRCH) {
			fprintf(stderr, "%s: process %d is still running\n",
					progname, astat.pid);
			return(-1);
		}
		/* assume it is dead */
	}
	if (strcmp(cfname, astat.cfname) && astat.pid != 0) {	/* other's */
		fprintf(stderr, "%s: unfinished job \"%s\"\n",
				progname, astat.cfname);
		return(-1);
	}
						/* check control file mods. */
	if (!nowarn && fdate(cfname) > fdate(sfname))
		fprintf(stderr,
			"%s: warning - control file modified since last run\n",
				progname);
setours:					/* set our values */
	strcpy(astat.host, myhostname());
	astat.pid = getpid();
	strcpy(astat.cfname, cfname);
	return(0);
fmterr:
	fprintf(stderr, "%s: format error in status file \"%s\"\n",
			progname, sfname);
	fclose(fp);
	return(-1);
}


static void
putastat(void)			/* put out current status */
{
	char	buf[256];
	FILE	*fp;

	if (noaction)
		return;
	sprintf(buf, "%s/%s", vval(DIRECTORY), SFNAME);
	if ((fp = fopen(buf, "w")) == NULL) {
		perror(buf);
		quit(1);
	}
	fprintf(fp, "Control host: %s\n", astat.host);
	fprintf(fp, "Control PID: %d\n", astat.pid);
	fprintf(fp, "Control file: %s\n", astat.cfname);
	fprintf(fp, "Next render: %d\n", astat.rnext);
	fprintf(fp, "Next filter: %d\n", astat.fnext);
	fprintf(fp, "Next transfer: %d\n", astat.tnext);
	fclose(fp);
}


static void
checkdir(void)			/* make sure we have our directory */
{
	struct stat	stb;

	if (!vdef(DIRECTORY)) {
		fprintf(stderr, "%s: %s undefined\n",
				progname, vnam(DIRECTORY));
		quit(1);
	}
	if (stat(vval(DIRECTORY), &stb) == -1) {
		if (errno == ENOENT && mkdir(vval(DIRECTORY), 0777) == 0)
			return;
		perror(vval(DIRECTORY));
		quit(1);
	}
	if (!(stb.st_mode & S_IFDIR)) {
		fprintf(stderr, "%s: not a directory\n", vval(DIRECTORY));
		quit(1);
	}
}


static void
setdefaults(void)			/* set default values */
{
	extern char	*atos();
	int	decades;
	char	buf[256];

	if (vdef(ANIMATE)) {
		vval(OCTREE) = NULL;
		vdef(OCTREE) = 0;
	} else if (!vdef(OCTREE)) {
		fprintf(stderr, "%s: either %s or %s must be defined\n",
				progname, vnam(OCTREE), vnam(ANIMATE));
		quit(1);
	}
	if (!vdef(VIEWFILE)) {
		fprintf(stderr, "%s: %s undefined\n", progname, vnam(VIEWFILE));
		quit(1);
	}
	if (!vdef(HOST)) {
		vval(HOST) = LHOSTNAME;
		vdef(HOST)++;
	}
	if (!vdef(START)) {
		vval(START) = "1";
		vdef(START)++;
	}
	if (!vdef(END)) {
		sprintf(buf, "%d", countviews()+vint(START)-1);
		vval(END) = savqstr(buf);
		vdef(END)++;
	}
	if (vint(END) < vint(START)) {
		fprintf(stderr, "%s: ending frame less than starting frame\n",
				progname);
		quit(1);
	}
	if (!vdef(BASENAME)) {
		decades = (int)log10((double)vint(END)) + 1;
		if (decades < 3) decades = 3;
		sprintf(buf, "%s/frame%%0%dd", vval(DIRECTORY), decades);
		vval(BASENAME) = savqstr(buf);
		vdef(BASENAME)++;
	}
	if (!vdef(RESOLUTION)) {
		vval(RESOLUTION) = "640";
		vdef(RESOLUTION)++;
	}
	if (!vdef(OVERSAMP)) {
		vval(OVERSAMP) = "2";
		vdef(OVERSAMP)++;
	}
	if (!vdef(INTERP)) {
		vval(INTERP) = "0";
		vdef(INTERP)++;
	}
	if (!vdef(MBLUR)) {
		vval(MBLUR) = "0";
		vdef(MBLUR)++;
	}
	if (!vdef(RTRACE)) {
		vval(RTRACE) = "F";
		vdef(RTRACE)++;
	}
	if (!vdef(DISKSPACE)) {
		if (!nowarn)
			fprintf(stderr,
		"%s: warning - no %s setting, assuming 100 Mbytes available\n",
					progname, vnam(DISKSPACE));
		vval(DISKSPACE) = "100";
		vdef(DISKSPACE)++;
	}
	if (!vdef(RSH)) {
		vval(RSH) = REMSH;
		vdef(RSH)++;
	}
				/* locate remote shell program */
	atos(buf, sizeof(buf), vval(RSH));
	if ((remsh = getpath(buf, getenv("PATH"), X_OK)) != NULL)
		remsh = savqstr(remsh);
	else
		remsh = vval(RSH);	/* will generate error if used */

				/* append rendering options */
	if (vdef(RENDER))
		sprintf(rendopt+strlen(rendopt), " %s", vval(RENDER));
}


static void
sethosts(void)			/* set up process servers */
{
	extern char	*iskip();
	char	buf[256], *dir, *uname;
	int	np;
	register char	*cp;
	int	i;

	npslots = 0;
	if (noaction)
		return;
	for (i = 0; i < vdef(HOST); i++) {	/* add each host */
		dir = uname = NULL;
		np = 1;
		strcpy(cp=buf, nvalue(HOST, i));	/* copy to buffer */
		cp = sskip(cp);				/* skip host name */
		while (isspace(*cp))
			*cp++ = '\0';
		if (*cp) {				/* has # processes? */
			np = atoi(cp);
			if ((cp = iskip(cp)) == NULL || (*cp && !isspace(*cp)))
				badvalue(HOST);
			while (isspace(*cp))
				cp++;
			if (*cp) {			/* has directory? */
				dir = cp;
				cp = sskip(cp);			/* skip dir. */
				while (isspace(*cp))
					*cp++ = '\0';
				if (*cp) {			/* has user? */
					uname = cp;
					if (*sskip(cp))
						badvalue(HOST);
				}
			}
		}
		if (addpserver(buf, dir, uname, np) == NULL) {
			if (!nowarn)
				fprintf(stderr,
					"%s: cannot execute on host \"%s\"\n",
						progname, buf);
		} else
			npslots += np;
	}
	if (npslots == 0) {
		fprintf(stderr, "%s: no working process servers\n", progname);
		quit(1);
	}
	pslot = (struct pslot *)calloc(npslots, sizeof(struct pslot));
	if (pslot == NULL) {
		perror("malloc");
		quit(1);
	}
}

static void
getradfile(char *rfargs)		/* run rad and get needed variables */
{
	static short	mvar[] = {OCTREE,PFILT,RESOLUTION,EXPOSURE,-1};
	char	combuf[256];
	register int	i;
	register char	*cp;
	char	*pippt = NULL;
					/* create rad command */
	sprintf(rendopt, " @%s/render.opt", vval(DIRECTORY));
	sprintf(combuf,
	"rad -v 0 -s -e -w %s OPTFILE=%s | egrep '^[ \t]*(NOMATCH",
			rfargs, rendopt+2);
	cp = combuf;
	while (*cp) {
		if (*cp == '|') pippt = cp;
		cp++;
	}				/* match unset variables */
	for (i = 0; mvar[i] >= 0; i++)
		if (!vdef(mvar[i])) {
			*cp++ = '|';
			strcpy(cp, vnam(mvar[i]));
			while (*cp) cp++;
			pippt = NULL;
		}
	if (pippt != NULL)
		strcpy(pippt, "> " NULL_DEVICE);	/* nothing to match */
	else {
		sprintf(cp, ")[ \t]*=' > %s/radset.var", vval(DIRECTORY));
		cp += 11;		/* point to file name */
	}
	system(combuf);			/* ignore exit code */
	if (pippt == NULL) {		/* load variables and remove file */
		loadvars(cp);
		unlink(cp);
	}
}


static void
animate(void)			/* run animation */
{
	int	xres, yres;
	float	pa, mult;
	int	frames_batch;
	register int	i;
	double	d1, d2;
					/* compute rpict resolution */
	i = sscanf(vval(RESOLUTION), "%d %d %f", &xres, &yres, &pa);
	mult = vflt(OVERSAMP);
	if (i == 3) {
		sprintf(rresopt, "-x %d -y %d -pa %.3f", (int)(mult*xres),
				(int)(mult*yres), pa);
		sprintf(fresopt, "-x %d -y %d -pa %.3f", xres, yres, pa);
	} else if (i) {
		if (i == 1) yres = xres;
		sprintf(rresopt, "-x %d -y %d", (int)(mult*xres),
				(int)(mult*yres));
		sprintf(fresopt, "-x %d -y %d -pa 1", xres, yres);
	} else
		badvalue(RESOLUTION);
					/* consistency checks */
	if (vdef(ANIMATE)) {
		if (vint(INTERP)) {
			if (!nowarn)
				fprintf(stderr,
					"%s: resetting %s=0 for animation\n",
						progname, vnam(INTERP));
			vval(INTERP) = "0";
		}
		if (strcmp(vval(MBLUR),"0")) {	/* can't handle this */
			if (!nowarn)
				fprintf(stderr,
					"%s: resetting %s=0 for animation\n",
						progname, vnam(MBLUR));
			vval(MBLUR) = "0";
		}
	}
					/* figure # frames per batch */
	d1 = mult*xres*mult*yres*4;		/* space for orig. picture */
	if ((i=vint(INTERP)) || getblur(NULL, NULL) > 1)
		d1 += mult*xres*mult*yres*sizeof(float);	/* Z-buffer */
	d2 = xres*yres*4;			/* space for final picture */
	frames_batch = (i+1)*(vflt(DISKSPACE)*1048576.-d1)/(d1+i*d2);
	if (frames_batch < i+2) {
		fprintf(stderr, "%s: insufficient disk space allocated\n",
				progname);
		quit(1);
	}
					/* initialize archive argument list */
	i = vdef(ARCHIVE) ? strlen(vval(ARCHIVE))+132 : 132;
	arcnext = arcfirst = arcargs + i;
					/* initialize status file */
	if (astat.rnext == 0)
		astat.rnext = astat.fnext = astat.tnext = vint(START);
	putastat();
					/* render in batches */
	while (astat.tnext <= vint(END)) {
		renderframes(frames_batch);
		filterframes();
		transferframes();
	}
					/* mark status as finished */
	astat.pid = 0;
	putastat();
					/* close open files */
	getview(0);
	getexp(0);
}


static void
renderframes(int nframes)		/* render next nframes frames */
{
	static char	vendbuf[16];
	VIEW	*vp;
	FILE	*fp = NULL;
	char	vfname[128];
	int	lastframe;
	register int	i;

	if (astat.tnext < astat.rnext)	/* other work to do first */
		return;
					/* create batch view file */
	if (!vdef(ANIMATE)) {
		sprintf(vfname, "%s/anim.vf", vval(DIRECTORY));
		if ((fp = fopen(vfname, "w")) == NULL) {
			perror(vfname);
			quit(1);
		}
	}
					/* bound batch properly */
	lastframe = astat.rnext + nframes - 1;
	if ((lastframe-1) % (vint(INTERP)+1))	/* need even interval */
		lastframe += vint(INTERP)+1 - ((lastframe-1)%(vint(INTERP)+1));
	if (lastframe > vint(END))		/* check for end */
		lastframe = vint(END);
					/* render each view */
	for (i = astat.rnext; i <= lastframe; i++) {
		if ((vp = getview(i)) == NULL) {
			if (!nowarn)
				fprintf(stderr,
				"%s: ran out of views before last frame\n",
					progname);
			sprintf(vval(END)=vendbuf, "%d", i-1);
			lastframe = i - 1;
			break;
		}
		if (vdef(ANIMATE))		/* animate frame */
			animrend(i, vp);
		else {				/* else record it */
			fputs(VIEWSTR, fp);
			fprintview(vp, fp);
			putc('\n', fp);
		}
	}
	if (vdef(ANIMATE))		/* wait for renderings to finish */
		bwait(0);
	else {				/* else if walk-through */
		fclose(fp);		/* close view file */
		walkwait(astat.rnext, lastframe, vfname);	/* walk it */
		unlink(vfname);		/* remove view file */
	}
	astat.rnext = i;		/* update status */
	putastat();
}


static void
filterframes(void)				/* catch up with filtering */
{
	register int	i;

	if (astat.tnext < astat.fnext)	/* other work to do first */
		return;
					/* filter each view */
	for (i = astat.fnext; i < astat.rnext; i++)
		dofilt(i, 0);

	bwait(0);			/* wait for filter processes */
	archive();			/* archive originals */
	astat.fnext = i;		/* update status */
	putastat();
}


static void
transferframes(void)			/* catch up with picture transfers */
{
	char	combuf[10240], *fbase;
	register char	*cp;
	register int	i;

	if (astat.tnext >= astat.fnext)	/* nothing to do, yet */
		return;
	if (!vdef(TRANSFER)) {		/* no transfer function -- leave 'em */
		astat.tnext = astat.fnext;
		putastat();		/* update status */
		return;
	}
	strcpy(combuf, "cd ");		/* start transfer command */
	fbase = dirfile(cp = combuf+3, vval(BASENAME));
	if (*cp) {
		while (*++cp) ;
		*cp++ = ';'; *cp++ = ' ';
	} else
		cp = combuf;
	strcpy(cp, vval(TRANSFER));
	while (*cp) cp++;
					/* make argument list */
	for (i = astat.tnext; i < astat.fnext; i++) {
		*cp++ = ' ';
		sprintf(cp, fbase, i);
		while (*cp) cp++;
		strcpy(cp, ".pic");
		cp += 4;
	}
	if (runcom(combuf)) {		/* transfer frames */
		fprintf(stderr, "%s: error running transfer command\n",
				progname);
		quit(1);
	}
	astat.tnext = i;		/* update status */
	putastat();
}


static void
animrend(			/* start animation frame */
int	frame,
VIEW	*vp
)
{
	char	combuf[2048];
	char	fname[128];

	sprintf(fname, vval(BASENAME), frame);
	strcat(fname, ".unf");
	if (access(fname, F_OK) == 0)
		return;
	sprintf(combuf, "%s %d | rpict%s%s -w0 %s > %s", vval(ANIMATE), frame,
			rendopt, viewopt(vp), rresopt, fname);
	bruncom(combuf, frame, recover);	/* run in background */
}


static void
walkwait(		/* walk-through frames */
int	first,
int last,
char	*vfn
)
{
	double	mblurf, dblurf;
	int	nblur = getblur(&mblurf, &dblurf);
	char	combuf[2048];
	register char	*inspoint;
	register int	i;

	if (!noaction && vint(INTERP))		/* create dummy frames */
		for (i = first; i <= last; i++)
			if (i < vint(END) && (i-1) % (vint(INTERP)+1)) {
				sprintf(combuf, vval(BASENAME), i);
				strcat(combuf, ".unf");
				close(open(combuf, O_RDONLY|O_CREAT, 0666));
			}
					/* create command */
	sprintf(combuf, "rpict%s%s -w0", rendopt,
			viewopt(getview(first>1 ? first-1 : 1)));
	inspoint = combuf;
	while (*inspoint) inspoint++;
	if (nblur) {
		sprintf(inspoint, " -pm %.3f", mblurf/nblur);
		while (*inspoint) inspoint++;
		sprintf(inspoint, " -pd %.3f", dblurf/nblur);
		while (*inspoint) inspoint++;
	}
	if (nblur > 1 || vint(INTERP)) {
		sprintf(inspoint, " -z %s.zbf", vval(BASENAME));
		while (*inspoint) inspoint++;
	}
	sprintf(inspoint, " -o %s.unf %s -S %d",
			vval(BASENAME), rresopt, first);
	while (*inspoint) inspoint++;
	sprintf(inspoint, " %s < %s", vval(OCTREE), vfn);
					/* run in parallel */
	i = (last-first+1)/(vint(INTERP)+1);
	if (i < 1) i = 1;
	if (pruncom(combuf, inspoint, i)) {
		fprintf(stderr, "%s: error rendering frames %d through %d\n",
				progname, first, last);
		quit(1);
	}
	if (!noaction && vint(INTERP))		/* remove dummy frames */
		for (i = first; i <= last; i++)
			if (i < vint(END) && (i-1) % (vint(INTERP)+1)) {
				sprintf(combuf, vval(BASENAME), i);
				strcat(combuf, ".unf");
				unlink(combuf);
			}
}


static int
recover(int frame)			/* recover the specified frame */
{
	static int	*rfrm;		/* list of recovered frames */
	static int	nrfrms = 0;
	double	mblurf, dblurf;
	int	nblur = getblur(&mblurf, &dblurf);
	char	combuf[2048];
	char	fname[128];
	register char	*cp;
	register int	i;
					/* check to see if recovered already */
	for (i = nrfrms; i--; )
		if (rfrm[i] == frame)
			return(0);
					/* build command */
	sprintf(fname, vval(BASENAME), frame);
	if (vdef(ANIMATE))
		sprintf(combuf, "%s %d | rpict%s -w0",
				vval(ANIMATE), frame, rendopt);
	else
		sprintf(combuf, "rpict%s -w0", rendopt);
	cp = combuf;
	while (*cp) cp++;
	if (nblur) {
		sprintf(cp, " -pm %.3f", mblurf/nblur);
		while (*cp) cp++;
		sprintf(cp, " -pd %.3f", dblurf/nblur);
		while (*cp) cp++;
	}
	if (nblur > 1 || vint(INTERP)) {
		sprintf(cp, " -z %s.zbf", fname);
		while (*cp) cp++;
	}
	sprintf(cp, " -ro %s.unf", fname);
	while (*cp) cp++;
	if (!vdef(ANIMATE)) {
		*cp++ = ' ';
		strcpy(cp, vval(OCTREE));
	}
	if (runcom(combuf))		/* run command */
		return(1);
					/* add frame to recovered list */
	if (nrfrms)
		rfrm = (int *)realloc((void *)rfrm, (nrfrms+1)*sizeof(int));
	else
		rfrm = (int *)malloc(sizeof(int));
	if (rfrm == NULL) {
		perror("malloc");
		quit(1);
	}
	rfrm[nrfrms++] = frame;
	return(0);
}


static int
frecover(int frame)				/* recover filtered frame */
{
	if (dofilt(frame, 2) && dofilt(frame, 1))
		return(1);
	return(0);
}


static void
archive(void)			/* archive and remove renderings */
{
#define RMCOML	(sizeof(rmcom)-1)
	static char	rmcom[] = "rm -f";
	char	basedir[128];
	int	dlen, alen;
	register int	j;

	if (arcnext == arcfirst)
		return;				/* nothing to do */
	dirfile(basedir, vval(BASENAME));
	dlen = strlen(basedir);
	if (vdef(ARCHIVE)) {			/* run archive command */
		alen = strlen(vval(ARCHIVE));
		if (dlen) {
			j = alen + dlen + 5;
			strncpy(arcfirst-j, "cd ", 3);
			strncpy(arcfirst-j+3, basedir, dlen);
			(arcfirst-j)[dlen+3] = ';'; (arcfirst-j)[dlen+4] = ' ';
		} else
			j = alen;
		strncpy(arcfirst-alen, vval(ARCHIVE), alen);
		if (runcom(arcfirst-j)) {
			fprintf(stderr, "%s: error running archive command\n",
					progname);
			quit(1);
		}
	}
	if (dlen) {
		j = RMCOML + dlen + 5;
		strncpy(arcfirst-j, "cd ", 3);
		strncpy(arcfirst-j+3, basedir, dlen);
		(arcfirst-j)[dlen+3] = ';'; (arcfirst-j)[dlen+4] = ' ';
	} else
		j = RMCOML;
						/* run remove command */
	strncpy(arcfirst-RMCOML, rmcom, RMCOML);
	runcom(arcfirst-j);
	arcnext = arcfirst;			/* reset argument list */
#undef RMCOML
}


static int
dofilt(				/* filter frame */
int	frame,
int	rvr
)
{
	static int	iter = 0;
	double	mblurf, dblurf;
	int	nblur = getblur(&mblurf, &dblurf);
	VIEW	*vp = getview(frame);
	char	*ep = getexp(frame);
	char	fnbefore[128], fnafter[128], *fbase;
	char	combuf[1024], fname0[128], fname1[128];
	int	usepinterp, usepfilt, nora_rgbe;
	int	frseq[2];
						/* check what is needed */
	if (vp == NULL) {
		 fprintf(stderr,
			"%s: unexpected error reading view for frame %d\n",
		                          progname, frame);
		quit(1);
	}
	usepinterp = (nblur > 1);
	usepfilt = pfiltalways | (ep==NULL);
	if (ep != NULL && !strcmp(ep, "1"))
		ep = "+0";
	nora_rgbe = strcmp(vval(OVERSAMP),"1") || ep==NULL ||
			*ep != '+' || *ep != '-' || !isint(ep);
						/* compute rendered views */
	frseq[0] = frame - ((frame-1) % (vint(INTERP)+1));
	frseq[1] = frseq[0] + vint(INTERP) + 1;
	fbase = dirfile(NULL, vval(BASENAME));
	if (frseq[1] > vint(END))
		frseq[1] = vint(END);
	if (frseq[1] == frame) {			/* pfilt only */
		frseq[0] = frseq[1];
		usepinterp = 0;			/* update what's needed */
		usepfilt |= nora_rgbe;
	} else if (frseq[0] == frame) {		/* no interpolation needed */
		if (!rvr && frame > 1+vint(INTERP)) {	/* archive previous */
			*arcnext++ = ' ';
			sprintf(arcnext, fbase, frame-vint(INTERP)-1);
			while (*arcnext) arcnext++;
			strcpy(arcnext, ".unf");
			arcnext += 4;
			if (usepinterp || vint(INTERP)) {	/* and Z-buf */
				*arcnext++ = ' ';
				sprintf(arcnext, fbase, frame-vint(INTERP)-1);
				while (*arcnext) arcnext++;
				strcpy(arcnext, ".zbf");
				arcnext += 4;
			}
		}
		if (!usepinterp)		/* update what's needed */
			usepfilt |= nora_rgbe;
	} else					/* interpolation needed */
		usepinterp++;
	if (frseq[1] >= astat.rnext)		/* next batch unavailable */
		frseq[1] = frseq[0];
	sprintf(fnbefore, vval(BASENAME), frseq[0]);
	sprintf(fnafter, vval(BASENAME), frseq[1]);
	if (rvr == 1 && recover(frseq[0]))	/* recover before frame? */
		return(1);
						/* generate command */
	if (usepinterp) {			/* using pinterp */
		if (rvr == 2 && recover(frseq[1]))	/* recover after? */
			return(1);
		if (nblur > 1) {		/* with pdmblur */
			sprintf(fname0, "%s/vw0%c", vval(DIRECTORY),
					'a'+(iter%26));
			sprintf(fname1, "%s/vw1%c", vval(DIRECTORY),
					'a'+(iter%26));
			if (!noaction) {
				FILE	*fp;		/* motion blurring */
				if ((fp = fopen(fname0, "w")) == NULL) {
					perror(fname0); quit(1);
				}
				fputs(VIEWSTR, fp);
				fprintview(vp, fp);
				putc('\n', fp); fclose(fp);
				if ((vp = getview(frame+1)) == NULL) {
					fprintf(stderr,
			"%s: unexpected error reading view for frame %d\n",
							progname, frame+1);
					quit(1);
				}
				if ((fp = fopen(fname1, "w")) == NULL) {
					perror(fname1); quit(1);
				}
				fputs(VIEWSTR, fp);
				fprintview(vp, fp);
				putc('\n', fp); fclose(fp);
			}
			sprintf(combuf,
			"(pmdblur %.3f %.3f %d %s %s; rm -f %s %s) | pinterp -B -a",
					mblurf, dblurf, nblur,
					fname0, fname1, fname0, fname1);
			iter++;
		} else				/* no blurring */
			strcpy(combuf, "pinterp");
		strcat(combuf, viewopt(vp));
		if (vbool(RTRACE))
			sprintf(combuf+strlen(combuf), " -ff -fr '%s -w0 %s'",
					rendopt+1, vval(OCTREE));
		if (vdef(PINTERP))
			sprintf(combuf+strlen(combuf), " %s", vval(PINTERP));
		if (usepfilt)
			sprintf(combuf+strlen(combuf), " %s", rresopt);
		else
			sprintf(combuf+strlen(combuf), " -a %s -e %s",
					fresopt, ep);
		sprintf(combuf+strlen(combuf), " %s.unf %s.zbf",
				fnbefore, fnbefore);
		if (frseq[1] != frseq[0])
			 sprintf(combuf+strlen(combuf), " %s.unf %s.zbf",
					fnafter, fnafter);
		if (usepfilt) {			/* also pfilt */
			if (vdef(PFILT))
				sprintf(combuf+strlen(combuf), " | pfilt %s",
						vval(PFILT));
			else
				strcat(combuf, " | pfilt");
			if (ep != NULL)
				sprintf(combuf+strlen(combuf), " -1 -e %s %s",
						ep, fresopt);
			else
				sprintf(combuf+strlen(combuf), " %s", fresopt);
		}
	} else if (usepfilt) {			/* pfilt only */
		if (rvr == 2)
			return(1);
		if (vdef(PFILT))
			sprintf(combuf, "pfilt %s", vval(PFILT));
		else
			strcpy(combuf, "pfilt");
		if (ep != NULL)
			sprintf(combuf+strlen(combuf), " -1 -e %s %s %s.unf",
				ep, fresopt, fnbefore);
		else
			sprintf(combuf+strlen(combuf), " %s %s.unf",
					fresopt, fnbefore);
	} else {				/* else just check it */
		if (rvr == 2)
			return(1);
		sprintf(combuf, "ra_rgbe -e %s -r %s.unf", ep, fnbefore);
	}
						/* output file name */
	sprintf(fname0, vval(BASENAME), frame);
	sprintf(combuf+strlen(combuf), " > %s.pic", fname0);
	if (rvr)				/* in recovery */
		return(runcom(combuf));
	bruncom(combuf, frame, frecover);	/* else run in background */
	return(0);
}


static VIEW *
getview(int n)			/* get view number n */
{
	static FILE	*viewfp = NULL;		/* view file pointer */
	static int	viewnum = 0;		/* current view number */
	static VIEW	curview = STDVIEW;	/* current view */
	char	linebuf[256];

	if (n == 0) {			/* signal to close file and clean up */
		if (viewfp != NULL) {
			fclose(viewfp);
			viewfp = NULL;
			viewnum = 0;
			curview = stdview;
		}
		return(NULL);
	}
	if (viewfp == NULL) {			/* open file */
		if ((viewfp = fopen(vval(VIEWFILE), "r")) == NULL) {
			perror(vval(VIEWFILE));
			quit(1);
		}
	} else if (n > 0 && n < viewnum) {	/* rewind file */
		if (viewnum == 1 && feof(viewfp))
			return(&curview);		/* just one view */
		if (fseek(viewfp, 0L, 0) == EOF) {
			perror(vval(VIEWFILE));
			quit(1);
		}
		curview = stdview;
		viewnum = 0;
	}
	if (n < 0) {				/* get next view */
		register int	c = getc(viewfp);
		if (c == EOF)
			return((VIEW *)NULL);		/* that's it */
		ungetc(c, viewfp);
		n = viewnum + 1;
	}
	while (n > viewnum) {		/* scan to desired view */
		if (fgets(linebuf, sizeof(linebuf), viewfp) == NULL)
			return(viewnum==1 ? &curview : (VIEW *)NULL);
		if (isview(linebuf) && sscanview(&curview, linebuf) > 0)
			viewnum++;
	}
	return(&curview);		/* return it */
}


static int
countviews(void)			/* count views in view file */
{
	int	n;

	if (getview(n=1) == NULL)
		return(0);
	while (getview(-1) != NULL)
		n++;
	return(n);
}


static char *
getexp(int n)			/* get exposure for nth frame */
{
	extern char	*fskip();
	static char	expval[32];
	static FILE	*expfp = NULL;
	static long	*exppos;
	static int	curfrm;
	register char	*cp;

	if (n == 0) {				/* signal to close file */
		if (expfp != NULL) {
			fclose(expfp);
			free((void *)exppos);
			expfp = NULL;
		}
		return(NULL);
	} else if (n > vint(END))		/* request past end (error?) */
		return(NULL);
	if (!vdef(EXPOSURE))			/* no setting (auto) */
		return(NULL);
	if (isflt(vval(EXPOSURE)))		/* always the same */
		return(vval(EXPOSURE));
	if (expfp == NULL) {			/* open exposure file */
		if ((expfp = fopen(vval(EXPOSURE), "r")) == NULL) {
			fprintf(stderr,
			"%s: cannot open exposure file \"%s\"\n",
					progname, vval(EXPOSURE));
			quit(1);
		}
		curfrm = vint(END) + 1;		/* init lookup tab. */
		exppos = (long *)malloc(curfrm*sizeof(long *));
		if (exppos == NULL) {
			perror(progname);
			quit(1);
		}
		while (curfrm--)
			exppos[curfrm] = -1L;
		curfrm = 0;
	}
						/* find position in file */
	if (n-1 != curfrm && n != curfrm && exppos[n-1] >= 0 &&
				fseek(expfp, exppos[curfrm=n-1], 0) == EOF) {
		fprintf(stderr, "%s: seek error on exposure file\n", progname);
		quit(1);
	}
	while (n > curfrm) {			/* read exposure */
		if (exppos[curfrm] < 0)
			exppos[curfrm] = ftell(expfp);
		if (fgets(expval, sizeof(expval), expfp) == NULL) {
			fprintf(stderr, "%s: too few exposures\n",
					vval(EXPOSURE));
			quit(1);
		}
		curfrm++;
		cp = fskip(expval);			/* check format */
		if (cp != NULL)
			while (isspace(*cp))
				*cp++ = '\0';
		if (cp == NULL || *cp) {
			fprintf(stderr,
				"%s: exposure format error on line %d\n",
					vval(EXPOSURE), curfrm);
			quit(1);
		}
	}
	return(expval);				/* return value */
}


static struct pslot *
findpslot(RT_PID pid)			/* find or allocate a process slot */
{
	register struct pslot	*psempty = NULL;
	register int	i;

	for (i = 0; i < npslots; i++) {		/* look for match */
		if (pslot[i].pid == pid)
			return(pslot+i);
		if (psempty == NULL && pslot[i].pid == 0)
			psempty = pslot+i;
	}
	return(psempty);		/* return emtpy slot (error if NULL) */
}


static int
donecom(		/* clean up after finished process */
	PSERVER	*ps,
	int	pn,
	int	status
)
{
	register NETPROC	*pp;
	register struct pslot	*psl;

	pp = ps->proc + pn;
	if (pp->elen) {			/* pass errors */
		if (ps->hostname[0])
			fprintf(stderr, "%s: ", ps->hostname);
		fprintf(stderr, "Error output from: %s\n", pp->com);
		fputs(pp->errs, stderr);
		fflush(stderr);
		if (ps->hostname[0])
			status = 1;	/* because rsh doesn't return status */
	}
	lastpserver = NULL;
	psl = findpslot(pp->pid);	/* check for bruncom() slot */
	if (psl->pid) {
		if (status) {
			if (psl->rcvf != NULL)	/* attempt recovery */
				status = (*psl->rcvf)(psl->fout);
			if (status) {
				fprintf(stderr,
					"%s: error rendering frame %d\n",
						progname, psl->fout);
				quit(1);
			}
			lastpserver = ps;
		}
		psl->pid = 0;			/* free process slot */
	} else if (status)
		lastpserver = ps;
	freestr(pp->com);		/* free command string */
	return(status);
}


static int
serverdown(void)			/* check status of last process server */
{
	if (lastpserver == NULL || !lastpserver->hostname[0])
		return(0);
	if (pserverOK(lastpserver))	/* server still up? */
		return(0);
	delpserver(lastpserver);	/* else delete it */
	if (pslist == NULL) {
		fprintf(stderr, "%s: all process servers are down\n",
				progname);
		quit(1);
	}
	return(1);
}


static RT_PID
bruncom(		/* run a command in the background */
char	*com,
int	fout,
int	(*rf)()
)
{
	RT_PID	pid;
	register struct pslot	*psl;

	if (noaction) {
		if (!silent)
			printf("\t%s\n", com);	/* echo command */
		return(0);
	}
	com = savestr(com);		/* else start it when we can */
	while ((pid = startjob(NULL, com, donecom)) == -1)
		bwait(1);
	if (!silent) {				/* echo command */
		PSERVER	*ps;
		RT_PID	psn = pid;
		ps = findjob(&psn);
		printf("\t%s\n", com);
		printf("\tProcess started on %s\n", phostname(ps));
		fflush(stdout);
	}
	psl = findpslot(pid);		/* record info. in appropriate slot */
	psl->pid = pid;
	psl->fout = fout;
	psl->rcvf = rf;
	return(pid);
}


static void
bwait(int ncoms)				/* wait for batch job(s) to finish */
{
	int	status;

	if (noaction)
		return;
	while ((status = wait4job(NULL, -1)) != -1) {
		serverdown();		/* update server status */
		if (--ncoms == 0)
			break;		/* done enough */
	}
}


static int
pruncom(	/* run a command in parallel over network */
char	*com,
char	*ppins,
int	maxcopies
)
{
	int	retstatus = 0;
	int	hostcopies;
	char	buf[10240], *com1, *s;
	int	status;
	int	pfd;
	register int	n;
	register PSERVER	*ps;

	if (!silent)
		printf("\t%s\n", com);	/* echo command */
	if (noaction)
		return(0);
	fflush(stdout);
					/* start jobs on each server */
	for (ps = pslist; ps != NULL; ps = ps->next) {
		hostcopies = 0;
		if (maxcopies > 1 && ps->nprocs > 1 && ppins != NULL) {
			strcpy(com1=buf, com);	/* build -PP command */
			sprintf(com1+(ppins-com), " -PP %s/%s.persist",
					vval(DIRECTORY), phostname(ps));
			unlink(com1+(ppins-com)+5);
			strcat(com1, ppins);
		} else
			com1 = com;
		while (maxcopies > 0) {
			s = savestr(com1);
			if (startjob(ps, s, donecom) != -1) {
				sleep(20);
				hostcopies++;
				maxcopies--;
			} else {
				freestr(s);
				break;
			}
		}
		if (!silent && hostcopies) {
			if (hostcopies > 1)
				printf("\t%d duplicate processes", hostcopies);
			else
				printf("\tProcess");
			printf(" started on %s\n", phostname(ps));
			fflush(stdout);
		}
	}
					/* wait for jobs to finish */
	while ((status = wait4job(NULL, -1)) != -1)
		retstatus += status && !serverdown();
					/* terminate parallel rpict's */
	for (ps = pslist; ps != NULL; ps = ps->next) {
		sprintf(buf, "%s/%s.persist", vval(DIRECTORY), phostname(ps));
		if ((pfd = open(buf, O_RDONLY)) >= 0) {
			n = read(pfd, buf, sizeof(buf)-1);	/* get PID */
			buf[n] = '\0';
			close(pfd);
			for (n = 0; buf[n] && !isspace(buf[n]); n++)
				;
								/* terminate */
			sprintf(buf, "kill -ALRM %d", atoi(buf+n));
			wait4job(ps, startjob(ps, buf, NULL));
		}
	}
	return(retstatus);
}


static int
runcom(char *cs)			/* run a command locally and wait for it */
{
	if (!silent)		/* echo it */
		printf("\t%s\n", cs);
	if (noaction)
		return(0);
	fflush(stdout);		/* flush output and pass to shell */
	return(system(cs));
}


static int
rmfile(char *fn)			/* remove a file */
{
	if (!silent)
#ifdef _WIN32
		printf("\tdel %s\n", fn);
#else
		printf("\trm -f %s\n", fn);
#endif
	if (noaction)
		return(0);
	return(unlink(fn));
}


static void
badvalue(int vc)			/* report bad variable value and exit */
{
	fprintf(stderr, "%s: bad value for variable '%s'\n",
			progname, vnam(vc));
	quit(1);
}


static char *
dirfile(		/* separate path into directory and file */
char	*df,
register char	*path
)
{
	register int	i;
	int	psep;

	for (i = 0, psep = -1; path[i]; i++)
		if (path[i] == '/')
			psep = i;
	if (df != NULL) {
		if (psep == 0) {
			df[0] = '/';
			df[1] = '\0';
		} else if (psep > 0) {
			strncpy(df, path, psep);
			df[psep] = '\0';
		} else
			df[0] = '\0';
	}
	return(path+psep+1);
}


static int
getblur(double *mbf, double *dbf)	/* get # blur samples (and fraction) */
{
	double	mblurf, dblurf;
	int	nmblur, ndblur;
	char	*s;
					/* get motion blur */
	if (!vdef(MBLUR) || (mblurf = atof(vval(MBLUR))) < 0.0)
		mblurf = 0.0;
	if (mbf != NULL)
		*mbf = mblurf;
	if (mblurf <= FTINY)
		nmblur = 0;
	else if (!*(s = sskip(vval(MBLUR))))
		nmblur = DEF_NBLUR;
	else if ((nmblur = atoi(s)) <= 0)
		nmblur = 1;
					/* get depth-of-field blur */
	if (!vdef(DBLUR) || (dblurf = atof(vval(DBLUR))) < 0.0)
		dblurf = 0.0;
	if (dbf != NULL)
		*dbf = dblurf;
	if (dblurf <= FTINY)
		ndblur = 0;
	else if (!*(s = sskip(vval(DBLUR))))
		ndblur = DEF_NBLUR;
	else if ((ndblur = atoi(s)) <= 0)
		ndblur = 1;
	if ((nmblur == 1) & (ndblur == 1))
		return(1);
					/* return combined samples */
	return(nmblur + ndblur);
}
