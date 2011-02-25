#ifndef lint
static const char	RCSid[] = "$Id: rad.c,v 2.105 2011/02/25 06:51:51 greg Exp $";
#endif
/*
 * Executive program for oconv, rpict and pfilt
 */

#include "standard.h"

#include <ctype.h>
#include <time.h>

#include "platform.h"
#include "rtprocess.h"
#include "view.h"
#include "paths.h"
#include "vars.h"

#ifdef _WIN32
  #define DELCMD "del"
  #define RENAMECMD "rename"
#else
  #define DELCMD "rm -f"
  #define RENAMECMD "mv"
  #include <sys/types.h>
  #include <sys/wait.h>
#endif

				/* variables (alphabetical by name) */
#define AMBFILE		0		/* ambient file name */
#define DETAIL		1		/* level of scene detail */
#define EXPOSURE	2		/* picture exposure setting */
#define EYESEP		3		/* interocular distance */
#define ILLUM		4		/* mkillum input files */
#define INDIRECT	5		/* indirection in lighting */
#define MATERIAL	6		/* material files */
#define MKILLUM		7		/* mkillum options */
#define OBJECT		8		/* object files */
#define OCONV		9		/* oconv options */
#define OCTREE		10		/* octree file name */
#define OPTFILE		11		/* rendering options file */
#define PENUMBRAS	12		/* shadow penumbras are desired */
#define PFILT		13		/* pfilt options */
#define PICTURE		14		/* picture file root name */
#define QUALITY		15		/* desired rendering quality */
#define RAWFILE		16		/* raw picture file root name */
#define RENDER		17		/* rendering options */
#define REPORT		18		/* report frequency and errfile */
#define RESOLUTION	19		/* maximum picture resolution */
#define RPICT		20		/* rpict parameters */
#define RVU		21		/* rvu parameters */
#define SCENE		22		/* scene files */
#define UP		23		/* view up (X, Y or Z) */
#define VARIABILITY	24		/* level of light variability */
#define VIEWS		25		/* view(s) for picture(s) */
#define ZFILE		26		/* distance file root name */
#define ZONE		27		/* simulation zone */
				/* total number of variables */
int NVARS = 28;

VARIABLE	vv[] = {		/* variable-value pairs */
	{"AMBFILE",	3,	0,	NULL,	onevalue},
	{"DETAIL",	3,	0,	NULL,	qualvalue},
	{"EXPOSURE",	3,	0,	NULL,	fltvalue},
	{"EYESEP",	3,	0,	NULL,	fltvalue},
	{"illum",	3,	0,	NULL,	catvalues},
	{"INDIRECT",	3,	0,	NULL,	intvalue},
	{"materials",	3,	0,	NULL,	catvalues},
	{"mkillum",	3,	0,	NULL,	catvalues},
	{"objects",	3,	0,	NULL,	catvalues},
	{"oconv",	3,	0,	NULL,	catvalues},
	{"OCTREE",	3,	0,	NULL,	onevalue},
	{"OPTFILE",	3,	0,	NULL,	onevalue},
	{"PENUMBRAS",	3,	0,	NULL,	boolvalue},
	{"pfilt",	2,	0,	NULL,	catvalues},
	{"PICTURE",	3,	0,	NULL,	onevalue},
	{"QUALITY",	3,	0,	NULL,	qualvalue},
	{"RAWFILE",	3,	0,	NULL,	onevalue},
	{"render",	3,	0,	NULL,	catvalues},
	{"REPORT",	3,	0,	NULL,	onevalue},
	{"RESOLUTION",	3,	0,	NULL,	onevalue},
	{"rpict",	3,	0,	NULL,	catvalues},
	{"rvu",		3,	0,	NULL,	catvalues},
	{"scene",	3,	0,	NULL,	catvalues},
	{"UP",		2,	0,	NULL,	onevalue},
	{"VARIABILITY",	3,	0,	NULL,	qualvalue},
	{"view",	2,	0,	NULL,	NULL},
	{"ZFILE",	2,	0,	NULL,	onevalue},
	{"ZONE",	2,	0,	NULL,	onevalue},
};

				/* overture calculation file */
#ifdef NULL_DEVICE
char	overfile[] = NULL_DEVICE;
#else
char	overfile[] = "overture.unf";
#endif


time_t	scenedate;		/* date of latest scene or object file */
time_t	octreedate;		/* date of octree */
time_t	matdate;		/* date of latest material file */
time_t	illumdate;		/* date of last illum file */

char	*oct0name;		/* name of pre-mkillum octree */
time_t	oct0date;		/* date of pre-mkillum octree */
char	*oct1name;		/* name of post-mkillum octree */
time_t	oct1date;		/* date of post-mkillum octree (>= matdate) */

int	nowarn = 0;		/* no warnings */
int	explicate = 0;		/* explicate variables */
int	silent = 0;		/* do work silently */
int	touchonly = 0;		/* touch files only */
int	nprocs = 1;		/* maximum executing processes */
int	sayview = 0;		/* print view out */
char	*rvdevice = NULL;	/* rvu output device */
char	*viewselect = NULL;	/* specific view only */

#define DEF_RPICT_PATH	"rpict"		/* default rpict path */

				/* command paths */
char	c_oconv[256] = "oconv";
char	c_mkillum[256] = "mkillum";
char	c_rvu[256] = "rvu";
char	c_rpict[256] = DEF_RPICT_PATH;
char	c_rpiece[] = "rpiece";
char	c_pfilt[256] = "pfilt";

int	overture = 0;		/* overture calculation needed */

int	children_running = 0;	/* set negative in children */

char	*progname;		/* global argv[0] */
char	*rifname;		/* global rad input file name */

char	radname[PATH_MAX];	/* root Radiance file name */

#define	inchild()	(children_running < 0)

static void rootname(char	*rn, char	*fn);
static time_t checklast(char	*fnames);
static char * newfname(char	*orig, int	pred);
static void checkfiles(void);
static void getoctcube(double	org[3], double	*sizp);
static void setdefaults(void);
static void oconv(void);
static char * addarg(char	*op, char	*arg);
static void oconvopts(char	*oo);
static void mkillumopts(char	*mo);
static void checkambfile(void);
static double ambval(void);
static void renderopts(char	*op, char	*po);
static void lowqopts(char	*op, char	*po);
static void medqopts(char	*op, char	*po);
static void hiqopts(char	*op, char	*po);
static void xferopts(char	*ro);
static void pfiltopts(char	*po);
static int matchword(char	*s1, char	*s2);
static char * specview(char	*vs);
static char * getview(int	n, char	*vn);
static int myprintview(char	*vopts, FILE	*fp);
static void rvu(char	*opts, char	*po);
static void rpict(char	*opts, char	*po);
static int touch(char	*fn);
static int runcom(char	*cs);
static int rmfile(char	*fn);
static int mvfile(char	*fold, char	*fnew);
static int next_process(int	reserve);
static void wait_process(int	all);
static void finish_process(void);
static void badvalue(int	vc);
static void syserr(char	*s);


int
main(
	int	argc,
	char	*argv[]
)
{
	char	ropts[512];
	char	popts[64];
	int	i;

	progname = argv[0];
				/* get options */
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 's':
			silent++;
			break;
		case 'n':
			nprocs = 0;
			break;
		case 'N':
			nprocs = atoi(argv[++i]);
			if (nprocs < 0)
				nprocs = 0;
			break;
		case 't':
			touchonly++;
			break;
		case 'e':
			explicate++;
			break;
		case 'o':
			rvdevice = argv[++i];
			break;
		case 'V':
			sayview++;
			break;
		case 'v':
			viewselect = argv[++i];
			break;
		case 'w':
			nowarn++;
			break;
		default:
			goto userr;
		}
	if (i >= argc)
		goto userr;
	rifname = argv[i];
				/* assign Radiance root file name */
	rootname(radname, rifname);
				/* load variable values */
	loadvars(rifname);
				/* get any additional assignments */
	for (i++; i < argc; i++)
		if (setvariable(argv[i], matchvar) < 0) {
			fprintf(stderr, "%s: unknown variable: %s\n",
					progname, argv[i]);
			quit(1);
		}
				/* check assignments */
	checkvalues();
				/* check files and dates */
	checkfiles();
				/* set default values as necessary */
	setdefaults();
				/* print all values if requested */
	if (explicate)
		printvars(stdout);
				/* build octree (and run mkillum) */
	oconv();
				/* check date on ambient file */
	checkambfile();
				/* run simulation */
	renderopts(ropts, popts);
	xferopts(ropts);
	if (rvdevice != NULL)
		rvu(ropts, popts);
	else
		rpict(ropts, popts);
	quit(0);
userr:
	fprintf(stderr,
"Usage: %s [-w][-s][-n|-N npr][-t][-e][-V][-v view][-o dev] rfile [VAR=value ..]\n",
			progname);
	quit(1);
	return 1; /* pro forma return */
}


static void
rootname(		/* remove tail from end of fn */
	register char	*rn,
	register char	*fn
)
{
	char	*tp, *dp;

	for (tp = NULL, dp = rn; (*rn = *fn++); rn++)
		if (ISDIRSEP(*rn))
			dp = rn;
		else if (*rn == '.')
			tp = rn;
	if (tp != NULL && tp > dp)
		*tp = '\0';
}


static time_t
checklast(			/* check files and find most recent */
	register char	*fnames
)
{
	char	thisfile[PATH_MAX];
	time_t	thisdate, lastdate = 0;

	if (fnames == NULL)
		return(0);
	while ((fnames = nextword(thisfile, PATH_MAX, fnames)) != NULL) {
		if (thisfile[0] == '!' ||
				(thisfile[0] == '\\' && thisfile[1] == '!')) {
			if (!lastdate)
				lastdate = 1;
			continue;
		}
		if (!(thisdate = fdate(thisfile)))
			syserr(thisfile);
		if (thisdate > lastdate)
			lastdate = thisdate;
	}
	return(lastdate);
}


static char *
newfname(		/* create modified file name */
	char	*orig,
	int	pred
)
{
	register char	*cp;
	register int	n;
	int	suffix;

	n = 0; cp = orig; suffix = -1;		/* suffix position, length */
	while (*cp) {
		if (*cp == '.') suffix = n;
		else if (ISDIRSEP(*cp)) suffix = -1;
		cp++; n++;
	}
	if (suffix == -1) suffix = n;
	if ((cp = bmalloc(n+2)) == NULL)
		syserr(progname);
	strncpy(cp, orig, suffix);
	cp[suffix] = pred;			/* root name + pred + suffix */
	strcpy(cp+suffix+1, orig+suffix);
	return(cp);
}


static void
checkfiles(void)			/* check for existence and modified times */
{
	time_t	objdate;

	if (!vdef(OCTREE)) {
		if ((vval(OCTREE) = bmalloc(strlen(radname)+5)) == NULL)
			syserr(progname);
		sprintf(vval(OCTREE), "%s.oct", radname);
		vdef(OCTREE)++;
	} else if (vval(OCTREE)[0] == '!') {
		fprintf(stderr, "%s: illegal '%s' specification\n",
				progname, vnam(OCTREE));
		quit(1);
	}
	octreedate = fdate(vval(OCTREE));
	if (vdef(ILLUM)) {		/* illum requires secondary octrees */
		oct0name = newfname(vval(OCTREE), '0');
		oct1name = newfname(vval(OCTREE), '1');
		oct0date = fdate(oct0name);
		oct1date = fdate(oct1name);
	} else
		oct0name = oct1name = vval(OCTREE);
	if ((scenedate = checklast(vval(SCENE))) &&
			(objdate = checklast(vval(OBJECT))) > scenedate)
		scenedate = objdate;
	illumdate = checklast(vval(ILLUM));
	if (!octreedate & !scenedate & !illumdate) {
		fprintf(stderr, "%s: need '%s' or '%s' or '%s'\n", progname,
				vnam(OCTREE), vnam(SCENE), vnam(ILLUM));
		quit(1);
	}
	matdate = checklast(vval(MATERIAL));
}	


static void
getoctcube(		/* get octree bounding cube */
	double	org[3],
	double	*sizp
)
{
	static double	oorg[3], osiz = 0.;
	double	min[3], max[3];
	char	buf[1024];
	FILE	*fp;
	register int	i;

	if (osiz <= FTINY) {
		if (!nprocs && fdate(oct1name) <
				(scenedate>illumdate?scenedate:illumdate)) {
							/* run getbbox */
			sprintf(buf, "getbbox -w -h %s",
				vdef(SCENE) ? vval(SCENE) : vval(ILLUM));
			if ((fp = popen(buf, "r")) == NULL)
				syserr("getbbox");
			if (fscanf(fp, "%lf %lf %lf %lf %lf %lf",
					&min[0], &max[0], &min[1], &max[1],
					&min[2], &max[2]) != 6) {
				fprintf(stderr,
			"%s: error reading bounding box from getbbox\n",
						progname);
				quit(1);
			}
			for (i = 0; i < 3; i++)
				if (max[i] - min[i] > osiz)
					osiz = max[i] - min[i];
			for (i = 0; i < 3; i++)
				oorg[i] = (max[i]+min[i]-osiz)*.5;
			pclose(fp);
		} else {				/* from octree */
			oconv();	/* does nothing if done already */
			sprintf(buf, "getinfo -d < %s", oct1name);
			if ((fp = popen(buf, "r")) == NULL)
				syserr("getinfo");
			if (fscanf(fp, "%lf %lf %lf %lf", &oorg[0], &oorg[1],
					&oorg[2], &osiz) != 4) {
				fprintf(stderr,
			"%s: error reading bounding cube from getinfo\n",
						progname);
				quit(1);
			}
			pclose(fp);
		}
	}
	org[0] = oorg[0]; org[1] = oorg[1]; org[2] = oorg[2]; *sizp = osiz;
}


static void
setdefaults(void)			/* set default values for unassigned var's */
{
	double	org[3], lim[3], size;
	char	buf[128];

	if (!vdef(ZONE)) {
		getoctcube(org, &size);
		sprintf(buf, "E %g %g %g %g %g %g", org[0], org[0]+size,
				org[1], org[1]+size, org[2], org[2]+size);
		vval(ZONE) = savqstr(buf);
		vdef(ZONE)++;
	}
	if (!vdef(EYESEP)) {
		if (sscanf(vval(ZONE), "%*s %lf %lf %lf %lf %lf %lf",
				&org[0], &lim[0], &org[1], &lim[1],
				&org[2], &lim[2]) != 6)
			badvalue(ZONE);
		sprintf(buf, "%f",
			0.01*(lim[0]-org[0]+lim[1]-org[1]+lim[2]-org[2]));
		vval(EYESEP) = savqstr(buf);
		vdef(EYESEP)++;
	}
	if (!vdef(INDIRECT)) {
		vval(INDIRECT) = "0";
		vdef(INDIRECT)++;
	}
	if (!vdef(QUALITY)) {
		vval(QUALITY) = "L";
		vdef(QUALITY)++;
	}
	if (!vdef(RESOLUTION)) {
		vval(RESOLUTION) = "512";
		vdef(RESOLUTION)++;
	}
	if (!vdef(PICTURE)) {
		vval(PICTURE) = radname;
		vdef(PICTURE)++;
	}
	if (!vdef(VIEWS)) {
		vval(VIEWS) = "X";
		vdef(VIEWS)++;
	}
	if (!vdef(DETAIL)) {
		vval(DETAIL) = "M";
		vdef(DETAIL)++;
	}
	if (!vdef(PENUMBRAS)) {
		vval(PENUMBRAS) = "F";
		vdef(PENUMBRAS)++;
	}
	if (!vdef(VARIABILITY)) {
		vval(VARIABILITY) = "L";
		vdef(VARIABILITY)++;
	}
}


static void
oconv(void)				/* run oconv and mkillum if necessary */
{
	static char	illumtmp[] = "ilXXXXXX";
	char	combuf[PATH_MAX], ocopts[64], mkopts[1024];

	oconvopts(ocopts);		/* get options */
	if (octreedate < scenedate) {	/* check date on original octree */
		if (touchonly && octreedate)
			touch(vval(OCTREE));
		else {				/* build command */
			if (vdef(MATERIAL))
				sprintf(combuf, "%s%s %s %s > %s", c_oconv,
						ocopts, vval(MATERIAL),
						vval(SCENE), vval(OCTREE));
			else
				sprintf(combuf, "%s%s %s > %s", c_oconv, ocopts,
						vval(SCENE), vval(OCTREE));
			
			if (runcom(combuf)) {		/* run it */
				fprintf(stderr,
				"%s: error generating octree\n\t%s removed\n",
						progname, vval(OCTREE));
				unlink(vval(OCTREE));
				quit(1);
			}
		}
		octreedate = time((time_t *)NULL);
		if (octreedate < scenedate)	/* in case clock is off */
			octreedate = scenedate;
	}
	if (oct1name == vval(OCTREE))		/* no mkillum? */
		oct1date = octreedate > matdate ? octreedate : matdate;
	if ((oct1date >= octreedate) & (oct1date >= matdate)
			& (oct1date >= illumdate))	/* all done */
		return;
						/* make octree0 */
	if ((oct0date < scenedate) | (oct0date < illumdate)) {
		if (touchonly && oct0date)
			touch(oct0name);
		else {				/* build command */
			if (octreedate)
				sprintf(combuf, "%s%s -i %s %s > %s", c_oconv,
					ocopts, vval(OCTREE),
					vval(ILLUM), oct0name);
			else if (vdef(MATERIAL))
				sprintf(combuf, "%s%s %s %s > %s", c_oconv,
					ocopts, vval(MATERIAL),
					vval(ILLUM), oct0name);
			else
				sprintf(combuf, "%s%s %s > %s", c_oconv,
					ocopts, vval(ILLUM), oct0name);
			if (runcom(combuf)) {		/* run it */
				fprintf(stderr,
				"%s: error generating octree\n\t%s removed\n",
						progname, oct0name);
				unlink(oct0name);
				quit(1);
			}
		}
		oct0date = time((time_t *)NULL);
		if (oct0date < octreedate)	/* in case clock is off */
			oct0date = octreedate;
		if (oct0date < illumdate)	/* ditto */
			oct0date = illumdate;
		}
	if (touchonly && oct1date)
		touch(oct1name);
	else {
		mkillumopts(mkopts);		/* build mkillum command */
		mktemp(illumtmp);
		sprintf(combuf, "%s%s %s \"<\" %s > %s", c_mkillum, mkopts,
				oct0name, vval(ILLUM), illumtmp);
		if (runcom(combuf)) {			/* run it */
			fprintf(stderr, "%s: error running %s\n",
					progname, c_mkillum);
			unlink(illumtmp);
			quit(1);
		}
						/* make octree1 (frozen) */
		if (octreedate)
			sprintf(combuf, "%s%s -f -i %s %s > %s", c_oconv,
				ocopts, vval(OCTREE), illumtmp, oct1name);
		else if (vdef(MATERIAL))
			sprintf(combuf, "%s%s -f %s %s > %s", c_oconv,
				ocopts, vval(MATERIAL), illumtmp, oct1name);
		else
			sprintf(combuf, "%s%s -f %s > %s", c_oconv, ocopts,
				illumtmp, oct1name);
		if (runcom(combuf)) {		/* run it */
			fprintf(stderr,
				"%s: error generating octree\n\t%s removed\n",
					progname, oct1name);
			unlink(oct1name);
			unlink(illumtmp);
			quit(1);
		}
		rmfile(illumtmp);
	}
	oct1date = time((time_t *)NULL);
	if (oct1date < oct0date)	/* in case clock is off */
		oct1date = oct0date;
}


static char *
addarg(				/* append argument and advance pointer */
register char	*op,
register char	*arg
)
{
	while (*op)
		op++;
	*op = ' ';
	while ( (*++op = *arg++) )
		;
	return(op);
}


static void
oconvopts(				/* get oconv options */
	register char	*oo
)
{
	/* BEWARE:  This may be called via setdefaults(), so no assumptions */

	*oo = '\0';
	if (vdef(OCONV))
		if (vval(OCONV)[0] != '-') {
			atos(c_oconv, sizeof(c_oconv), vval(OCONV));
			oo = addarg(oo, sskip2(vval(OCONV), 1));
		} else
			oo = addarg(oo, vval(OCONV));
}


static void
mkillumopts(				/* get mkillum options */
	char	*mo
)
{
	/* BEWARE:  This may be called via setdefaults(), so no assumptions */

	if (nprocs > 1)
		sprintf(mo, " -n %d", nprocs);
	else
		*mo = '\0';
	if (vdef(MKILLUM))
		if (vval(MKILLUM)[0] != '-') {
			atos(c_mkillum, sizeof(c_mkillum), vval(MKILLUM));
			mo = addarg(mo, sskip2(vval(MKILLUM), 1));
		} else
			mo = addarg(mo, vval(MKILLUM));
}


static void
checkambfile(void)			/* check date on ambient file */
{
	time_t	afdate;

	if (!vdef(AMBFILE))
		return;
	if (!(afdate = fdate(vval(AMBFILE))))
		return;
	if (oct1date > afdate) {
		if (touchonly)
			touch(vval(AMBFILE));
		else
			rmfile(vval(AMBFILE));
	}
}


static double
ambval(void)				/* compute ambient value */
{
	if (vdef(EXPOSURE)) {
		if (vval(EXPOSURE)[0] == '+' || vval(EXPOSURE)[0] == '-')
			return(.5/pow(2.,vflt(EXPOSURE)));
		return(.5/vflt(EXPOSURE));
	}
	if (vlet(ZONE) == 'E')
		return(10.);
	if (vlet(ZONE) == 'I')
		return(.01);
	badvalue(ZONE);
	return 0; /* pro forma return */
}


static void
renderopts(			/* set rendering options */
	char	*op,
	char	*po
)
{
	switch(vscale(QUALITY)) {
	case LOW:
		lowqopts(op, po);
		break;
	case MEDIUM:
		medqopts(op, po);
		break;
	case HIGH:
		hiqopts(op, po);
		break;
	}
	if (vdef(RENDER))
		op = addarg(op, vval(RENDER));
	if (rvdevice != NULL) {
		if (vdef(RVU))
			if (vval(RVU)[0] != '-') {
				atos(c_rvu, sizeof(c_rvu), vval(RVU));
				po = addarg(po, sskip2(vval(RVU), 1));
			} else
				po = addarg(po, vval(RVU));
	} else {
		if (vdef(RPICT))
			if (vval(RPICT)[0] != '-') {
				atos(c_rpict, sizeof(c_rpict), vval(RPICT));
				po = addarg(po, sskip2(vval(RPICT), 1));
			} else
				po = addarg(po, vval(RPICT));
	}
}


static void
lowqopts(			/* low quality rendering options */
	register char	*op,
	char	*po
)
{
	double	d, org[3], siz[3];

	*op = '\0';
	*po = '\0';
	if (sscanf(vval(ZONE), "%*s %lf %lf %lf %lf %lf %lf", &org[0],
			&siz[0], &org[1], &siz[1], &org[2], &siz[2]) != 6)
		badvalue(ZONE);
	siz[0] -= org[0]; siz[1] -= org[1]; siz[2] -= org[2];
	if ((siz[0] <= FTINY) | (siz[1] <= FTINY) | (siz[2] <= FTINY))
		badvalue(ZONE);
	getoctcube(org, &d);
	d *= 3./(siz[0]+siz[1]+siz[2]);
	switch (vscale(DETAIL)) {
	case LOW:
		po = addarg(po, "-ps 16");
		op = addarg(op, "-dp 64");
		sprintf(op, " -ar %d", (int)(8*d));
		op += strlen(op);
		break;
	case MEDIUM:
		po = addarg(po, "-ps 8");
		op = addarg(op, "-dp 128");
		sprintf(op, " -ar %d", (int)(16*d));
		op += strlen(op);
		break;
	case HIGH:
		po = addarg(po, "-ps 4");
		op = addarg(op, "-dp 256");
		sprintf(op, " -ar %d", (int)(32*d));
		op += strlen(op);
		break;
	}
	po = addarg(po, "-pt .16");
	if (vbool(PENUMBRAS))
		op = addarg(op, "-ds .4");
	else
		op = addarg(op, "-ds 0");
	op = addarg(op, "-dt .2 -dc .25 -dr 0 -ss 0 -st .5");
	if (vdef(AMBFILE)) {
		sprintf(op, " -af %s", vval(AMBFILE));
		op += strlen(op);
	} else
		overture = 0;
	switch (vscale(VARIABILITY)) {
	case LOW:
		op = addarg(op, "-aa .3 -ad 256");
		break;
	case MEDIUM:
		op = addarg(op, "-aa .25 -ad 512");
		break;
	case HIGH:
		op = addarg(op, "-aa .2 -ad 1024");
		break;
	}
	op = addarg(op, "-as 0");
	d = ambval();
	sprintf(op, " -av %.2g %.2g %.2g", d, d, d);
	op += strlen(op);
	op = addarg(op, "-lr 6 -lw .003");
}


static void
medqopts(			/* medium quality rendering options */
	register char	*op,
	char	*po
)
{
	double	d, org[3], siz[3], asz;

	*op = '\0';
	*po = '\0';
	if (sscanf(vval(ZONE), "%*s %lf %lf %lf %lf %lf %lf", &org[0],
			&siz[0], &org[1], &siz[1], &org[2], &siz[2]) != 6)
		badvalue(ZONE);
	siz[0] -= org[0]; siz[1] -= org[1]; siz[2] -= org[2];
	if ((siz[0] <= FTINY) | (siz[1] <= FTINY) | (siz[2] <= FTINY))
		badvalue(ZONE);
	getoctcube(org, &d);
	asz = (siz[0]+siz[1]+siz[2])/3.;
	d /= asz;
	switch (vscale(DETAIL)) {
	case LOW:
		po = addarg(po, vbool(PENUMBRAS) ? "-ps 4" : "-ps 8");
		op = addarg(op, "-dp 256");
		sprintf(op, " -ar %d", (int)(16*d));
		op += strlen(op);
		sprintf(op, " -ms %.2g", asz/20.);
		op += strlen(op);
		break;
	case MEDIUM:
		po = addarg(po, vbool(PENUMBRAS) ? "-ps 3" : "-ps 6");
		op = addarg(op, "-dp 512");
		sprintf(op, " -ar %d", (int)(32*d));
		op += strlen(op);
		sprintf(op, " -ms %.2g", asz/40.);
		op += strlen(op);
		break;
	case HIGH:
		po = addarg(po, vbool(PENUMBRAS) ? "-ps 2" : "-ps 4");
		op = addarg(op, "-dp 1024");
		sprintf(op, " -ar %d", (int)(64*d));
		op += strlen(op);
		sprintf(op, " -ms %.2g", asz/80.);
		op += strlen(op);
		break;
	}
	po = addarg(po, "-pt .08");
	if (vbool(PENUMBRAS))
		op = addarg(op, "-ds .2 -dj .9");
	else
		op = addarg(op, "-ds .3");
	op = addarg(op, "-dt .1 -dc .5 -dr 1 -ss 1 -st .1");
	if ( (overture = vint(INDIRECT)) ) {
		sprintf(op, " -ab %d", overture);
		op += strlen(op);
	}
	if (vdef(AMBFILE)) {
		sprintf(op, " -af %s", vval(AMBFILE));
		op += strlen(op);
	} else
		overture = 0;
	switch (vscale(VARIABILITY)) {
	case LOW:
		op = addarg(op, "-aa .2 -ad 329 -as 42");
		break;
	case MEDIUM:
		op = addarg(op, "-aa .15 -ad 800 -as 128");
		break;
	case HIGH:
		op = addarg(op, "-aa .1 -ad 1536 -as 392");
		break;
	}
	d = ambval();
	sprintf(op, " -av %.2g %.2g %.2g", d, d, d);
	op += strlen(op);
	op = addarg(op, "-lr 8 -lw 1e-4");
}


static void
hiqopts(				/* high quality rendering options */
	register char	*op,
	char	*po
)
{
	double	d, org[3], siz[3], asz;

	*op = '\0';
	*po = '\0';
	if (sscanf(vval(ZONE), "%*s %lf %lf %lf %lf %lf %lf", &org[0],
			&siz[0], &org[1], &siz[1], &org[2], &siz[2]) != 6)
		badvalue(ZONE);
	siz[0] -= org[0]; siz[1] -= org[1]; siz[2] -= org[2];
	if ((siz[0] <= FTINY) | (siz[1] <= FTINY) | (siz[2] <= FTINY))
		badvalue(ZONE);
	getoctcube(org, &d);
	asz = (siz[0]+siz[1]+siz[2])/3.;
	d /= asz;
	switch (vscale(DETAIL)) {
	case LOW:
		po = addarg(po, vbool(PENUMBRAS) ? "-ps 1" : "-ps 8");
		op = addarg(op, "-dp 1024");
		sprintf(op, " -ar %d", (int)(32*d));
		op += strlen(op);
		sprintf(op, " -ms %.2g", asz/40.);
		op += strlen(op);
		break;
	case MEDIUM:
		po = addarg(po, vbool(PENUMBRAS) ? "-ps 1" : "-ps 5");
		op = addarg(op, "-dp 2048");
		sprintf(op, " -ar %d", (int)(64*d));
		op += strlen(op);
		sprintf(op, " -ms %.2g", asz/80.);
		op += strlen(op);
		break;
	case HIGH:
		po = addarg(po, vbool(PENUMBRAS) ? "-ps 1" : "-ps 3");
		op = addarg(op, "-dp 4096");
		sprintf(op, " -ar %d", (int)(128*d));
		op += strlen(op);
		sprintf(op, " -ms %.2g", asz/160.);
		op += strlen(op);
		break;
	}
	po = addarg(po, "-pt .04");
	if (vbool(PENUMBRAS))
		op = addarg(op, "-ds .1 -dj .9");
	else
		op = addarg(op, "-ds .2");
	op = addarg(op, "-dt .05 -dc .75 -dr 3 -ss 16 -st .01");
	sprintf(op, " -ab %d", overture=vint(INDIRECT)+1);
	op += strlen(op);
	if (vdef(AMBFILE)) {
		sprintf(op, " -af %s", vval(AMBFILE));
		op += strlen(op);
	} else
		overture = 0;
	switch (vscale(VARIABILITY)) {
	case LOW:
		op = addarg(op, "-aa .125 -ad 512 -as 64");
		break;
	case MEDIUM:
		op = addarg(op, "-aa .1 -ad 1536 -as 768");
		break;
	case HIGH:
		op = addarg(op, "-aa .075 -ad 4096 -as 2048");
		break;
	}
	d = ambval();
	sprintf(op, " -av %.2g %.2g %.2g", d, d, d);
	op += strlen(op);
	op = addarg(op, "-lr 12 -lw 1e-5");
}


static void
xferopts(				/* transfer options if indicated */
	char	*ro
)
{
	int	fd, n;
	register char	*cp;
	
	n = strlen(ro);
	if (n < 2)
		return;
	if (vdef(OPTFILE)) {
		for (cp = ro; cp[1]; cp++)
			if (isspace(cp[1]) && (cp[2] == '@' ||
					(cp[2] == '-' && isalpha(cp[3]))))
				*cp = '\n';
			else
				*cp = cp[1];
		*cp = '\n';
		fd = open(vval(OPTFILE), O_WRONLY|O_CREAT|O_TRUNC, 0666);
		if (fd < 0 || write(fd, ro, n) != n || close(fd) < 0)
			syserr(vval(OPTFILE));
		sprintf(ro, " @%s", vval(OPTFILE));
	}
#ifdef _WIN32
	else if (n > 50) {
		setenv("ROPT", ro+1);
		strcpy(ro, " $ROPT");
	}
#endif
}


static void
pfiltopts(				/* get pfilt options */
	register char	*po
)
{
	*po = '\0';
	if (vdef(EXPOSURE)) {
		po = addarg(po, "-1 -e");
		po = addarg(po, vval(EXPOSURE));
	}
	switch (vscale(QUALITY)) {
	case MEDIUM:
		po = addarg(po, "-r .6");
		break;
	case HIGH:
		po = addarg(po, "-m .25");
		break;
	}
	if (vdef(PFILT))
		if (vval(PFILT)[0] != '-') {
			atos(c_pfilt, sizeof(c_pfilt), vval(PFILT));
			po = addarg(po, sskip2(vval(PFILT), 1));
		} else
			po = addarg(po, vval(PFILT));
}


static int
matchword(			/* match white-delimited words */
	register char	*s1,
	register char	*s2
)
{
	while (isspace(*s1)) s1++;
	while (isspace(*s2)) s2++;
	while (*s1 && !isspace(*s1))
		if (*s1++ != *s2++)
			return(0);
	return(!*s2 || isspace(*s2));
}


static char *
specview(				/* get proper view spec from vs */
	register char	*vs
)
{
	static char	vup[7][12] = {"-vu 0 0 -1","-vu 0 -1 0","-vu -1 0 0",
			"-vu 0 0 1", "-vu 1 0 0","-vu 0 1 0","-vu 0 0 1"};
	static char	viewopts[128];
	register char	*cp;
	int	xpos, ypos, zpos, viewtype, upax;
	register int	i;
	double	cent[3], dim[3], mult, d;

	if (vs == NULL || *vs == '-')
		return(vs);
	upax = 0;			/* get the up vector */
	if (vdef(UP)) {
		if (vval(UP)[0] == '-' || vval(UP)[0] == '+')
			upax = 1-'X'+UPPER(vval(UP)[1]);
		else
			upax = 1-'X'+vlet(UP);
		if ((upax < 1) | (upax > 3))
			badvalue(UP);
		if (vval(UP)[0] == '-')
			upax = -upax;
	}
					/* check standard view names */
	xpos = ypos = zpos = 0;
	if (*vs == 'X') {
		xpos = 1; vs++;
	} else if (*vs == 'x') {
		xpos = -1; vs++;
	}
	if (*vs == 'Y') {
		ypos = 1; vs++;
	} else if (*vs == 'y') {
		ypos = -1; vs++;
	}
	if (*vs == 'Z') {
		zpos = 1; vs++;
	} else if (*vs == 'z') {
		zpos = -1; vs++;
	}
	switch (*vs) {
	case VT_PER:
	case VT_PAR:
	case VT_ANG:
	case VT_HEM:
	case VT_PLS:
	case VT_CYL:
		viewtype = *vs++;
		break;
	default:
		viewtype = VT_PER;
		break;
	}
	cp = viewopts;
	if ((!*vs || isspace(*vs)) && (xpos|ypos|zpos)) {	/* got one! */
		*cp++ = '-'; *cp++ = 'v'; *cp++ = 't'; *cp++ = viewtype;
		if (sscanf(vval(ZONE), "%*s %lf %lf %lf %lf %lf %lf",
				&cent[0], &dim[0], &cent[1], &dim[1],
				&cent[2], &dim[2]) != 6)
			badvalue(ZONE);
		for (i = 0; i < 3; i++) {
			dim[i] -= cent[i];
			if (dim[i] <= FTINY)
				badvalue(ZONE);
			cent[i] += .5*dim[i];
		}
		mult = vlet(ZONE)=='E' ? 2. : .45 ;
		sprintf(cp, " -vp %.2g %.2g %.2g -vd %.2g %.2g %.2g",
				cent[0]+xpos*mult*dim[0],
				cent[1]+ypos*mult*dim[1],
				cent[2]+zpos*mult*dim[2],
				-xpos*dim[0], -ypos*dim[1], -zpos*dim[2]);
		cp += strlen(cp);
					/* redirect up axis if necessary */
		switch (upax) {
		case 3:			/* plus or minus Z axis */
		case -3:
		case 0:
			if (!(xpos|ypos))
				upax = 2;
			break;
		case 2:			/* plus or minus Y axis */
		case -2:
			if (!(xpos|zpos))
				upax = 1;
			break;
		case 1:			/* plus or minus X axis */
		case -1:
			if (!(ypos|zpos))
				upax = 3;
			break;
		}
		cp = addarg(cp, vup[upax+3]);
		switch (viewtype) {
		case VT_PER:
			cp = addarg(cp, "-vh 45 -vv 45");
			break;
		case VT_PAR:
			d = sqrt(dim[0]*dim[0]+dim[1]*dim[1]+dim[2]*dim[2]);
			sprintf(cp, " -vh %.2g -vv %.2g", d, d);
			cp += strlen(cp);
			break;
		case VT_ANG:
		case VT_HEM:
		case VT_PLS:
			cp = addarg(cp, "-vh 180 -vv 180");
			break;
		case VT_CYL:
			cp = addarg(cp, "-vh 180 -vv 90");
			break;
		}
	} else {
		while (!isspace(*vs))		/* else skip id */
			if (!*vs++)
				return(NULL);
		if (upax) {			/* prepend up vector */
			strcpy(cp, vup[upax+3]);
			cp += strlen(cp);
		}
	}
	if (cp == viewopts)		/* append any additional options */
		vs++;		/* skip prefixed space if unneeded */
	strcpy(cp, vs);
#ifdef _WIN32
	if (strlen(viewopts) > 40) {
		setenv("VIEW", viewopts);
		return("$VIEW");
	}
#endif
	return(viewopts);
}


static char *
getview(				/* get view n, or NULL if none */
	int	n,
	char	*vn		/* returned view name */
)
{
	register char	*mv;

	if (viewselect != NULL) {		/* command-line selected */
		if (n)				/* only do one */
			return(NULL);
					
		if (isint(viewselect)) {	/* view number? */
			n = atoi(viewselect)-1;
			goto numview;
		}
		if (viewselect[0] == '-') {	/* already specified */
			if (vn != NULL)
				strcpy(vn, "0");
			return(viewselect);
		}
		if (vn != NULL) {
			for (mv = viewselect; *mv && !isspace(*mv);
					*vn++ = *mv++)
				;
			*vn = '\0';
		}
						/* check list */
		while ((mv = nvalue(VIEWS, n++)) != NULL)
			if (matchword(viewselect, mv))
				return(specview(mv));

		return(specview(viewselect));	/* standard view? */
	}
numview:
	mv = nvalue(VIEWS, n);		/* use view n */
	if ((vn != NULL) & (mv != NULL))
		if (*mv != '-') {
			register char	*mv2 = mv;
			while (*mv2 && !isspace(*mv2))
				*vn++ = *mv2++;
			*vn = '\0';
		} else
			sprintf(vn, "%d", n+1);

	return(specview(mv));
}


static int
myprintview(			/* print out selected view */
	register char	*vopts,
	FILE	*fp
)
{
	VIEW	vwr;
	char	buf[128];
	register char	*cp;
#ifdef _WIN32
/* XXX Should we allow something like this for all platforms? */
/* XXX Or is it still required at all? */
again:
#endif
	if (vopts == NULL)
		return(-1);
#ifdef _WIN32
	if (vopts[0] == '$') {
		vopts = getenv(vopts+1);
		goto again;
	}
#endif
	vwr = stdview;
	sscanview(&vwr, cp=vopts);		/* set initial options */
	while ((cp = strstr(cp, "-vf ")) != NULL &&
			*atos(buf, sizeof(buf), cp += 4)) {
		viewfile(buf, &vwr, NULL);	/* load -vf file */
		sscanview(&vwr, cp);		/* reset tail */
	}
	fputs(VIEWSTR, fp);
	fprintview(&vwr, fp);			/* print full spec. */
	fputc('\n', fp);
	return(0);
}


static void
rvu(				/* run rvu with first view */
	char	*opts,
	char	*po
)
{
	char	*vw;
	char	combuf[PATH_MAX];
					/* build command */
	if (touchonly || (vw = getview(0, NULL)) == NULL)
		return;
	if (sayview)
		myprintview(vw, stdout);
	sprintf(combuf, "%s %s%s%s -R %s ", c_rvu, vw, opts, po, rifname);
	if (nprocs > 1)
		sprintf(combuf+strlen(combuf), "-n %d ", nprocs);
	if (rvdevice != NULL)
		sprintf(combuf+strlen(combuf), "-o %s ", rvdevice);
	if (vdef(EXPOSURE))
		sprintf(combuf+strlen(combuf), "-pe %s ", vval(EXPOSURE));
	strcat(combuf, oct1name);
	if (runcom(combuf)) {		/* run it */
		fprintf(stderr, "%s: error running %s\n", progname, c_rvu);
		quit(1);
	}
}


static void
rpict(				/* run rpict and pfilt for each view */
	char	*opts,
	char	*po
)
{
#define do_rpiece	(sfile[0]!='\0')
	char	combuf[5*PATH_MAX+512];
	char	rawfile[PATH_MAX], picfile[PATH_MAX];
	char	zopt[PATH_MAX+4], rep[PATH_MAX+16], res[32];
	char	rppopt[32], sfile[PATH_MAX], *pfile = NULL;
	char	pfopts[128];
	char	vs[32], *vw;
	int	vn, mult;
	FILE	*fp;
	time_t	rfdt, pfdt;
	int	xres, yres;
	double	aspect;
	int	n;
					/* get pfilt options */
	pfiltopts(pfopts);
					/* get resolution, reporting */
	switch (vscale(QUALITY)) {
	case LOW:
		mult = 1;
		break;
	case MEDIUM:
		mult = 2;
		break;
	case HIGH:
		mult = 3;
		break;
	}
	n = sscanf(vval(RESOLUTION), "%d %d %lf", &xres, &yres, &aspect);
	if (n == 3)
		sprintf(res, "-x %d -y %d -pa %.3f",
				mult*xres, mult*yres, aspect);
	else if (n) {
		aspect = 1.;
		if (n == 1) yres = xres;
		sprintf(res, "-x %d -y %d", mult*xres, mult*yres);
	} else
		badvalue(RESOLUTION);
	rep[0] = '\0';
	if (vdef(REPORT)) {
		double	minutes;
		n = sscanf(vval(REPORT), "%lf %s", &minutes, rawfile);
		if (n == 2)
			sprintf(rep, " -t %d -e %s", (int)(minutes*60), rawfile);
		else if (n == 1)
			sprintf(rep, " -t %d", (int)(minutes*60));
		else
			badvalue(REPORT);
	}
					/* set up parallel rendering */
	sfile[0] = '\0';
	if ((nprocs > 1) & !touchonly & !vdef(ZFILE) &&
					getview(0, vs) != NULL) {
		if (!strcmp(c_rpict, DEF_RPICT_PATH) &&
				getview(1, NULL) == NULL) {
			sprintf(sfile, "%s_%s_rpsync.txt", 
				vdef(RAWFILE) ? vval(RAWFILE) : vval(PICTURE),
					vs);
			strcpy(rppopt, "-PP pfXXXXXX");
		} else {
			strcpy(rppopt, "-S 1 -PP pfXXXXXX");
		}
		pfile = rppopt + strlen(rppopt) - 8;
		if (mktemp(pfile) == NULL) {
			if (do_rpiece) {
				fprintf(stderr, "%s: cannot create\n", pfile);
				quit(1);
			}
			pfile = NULL;
		}
	}
	vn = 0;					/* do each view */
	while ((vw = getview(vn++, vs)) != NULL) {
		if (sayview)
			myprintview(vw, stdout);
		sprintf(picfile, "%s_%s.hdr", vval(PICTURE), vs);
		if (vdef(ZFILE))
			sprintf(zopt, " -z %s_%s.zbf", vval(ZFILE), vs);
		else
			zopt[0] = '\0';
						/* check date on picture */
		pfdt = fdate(picfile);
		if (pfdt >= oct1date)
			continue;
						/* get raw file name */
		sprintf(rawfile, "%s_%s.unf",
			vdef(RAWFILE) ? vval(RAWFILE) : vval(PICTURE), vs);
		rfdt = fdate(rawfile);
		if (touchonly) {		/* update times only */
			if (rfdt) {
				if (rfdt < oct1date)
					touch(rawfile);
			} else if (pfdt && pfdt < oct1date)
				touch(picfile);
			continue;
		}
						/* parallel running? */
		if (do_rpiece) {
			if (rfdt < oct1date || !fdate(sfile)) {
				int	xdiv = 8+nprocs/3, ydiv = 8+nprocs/3;
				if (rfdt >= oct1date) {
					fprintf(stderr,
		"%s: partial output not created with rpiece\n", rawfile);
					quit(1);
				}
				if (rfdt) {	/* start fresh */
					rmfile(rawfile);
					rfdt = 0;
				}
				if (!silent)
					printf("\techo %d %d > %s\n",
							xdiv, ydiv, sfile);
				if ((fp = fopen(sfile, "w")) == NULL) {
					fprintf(stderr, "%s: cannot create\n",
							sfile);
					quit(1);
				}
				fprintf(fp, "%d %d\n", xdiv, ydiv);
				fclose(fp);
			}
		} else if (next_process(0)) {
			if (pfile != NULL)
				sleep(10);
			continue;
		} else if (!inchild())
			pfile = NULL;
		/* XXX Remember to call finish_process() */
						/* build rpict command */
		if (rfdt >= oct1date) {		/* already in progress */
			if (do_rpiece) {
				sprintf(combuf, "%s -R %s %s%s %s %s%s%s -o %s %s",
						c_rpiece, sfile, rppopt, rep, vw,
						res, opts, po, rawfile, oct1name);
				while (next_process(1)) {
					sleep(10);
					combuf[strlen(c_rpiece)+2] = 'F';
				}
			} else
				sprintf(combuf, "%s%s%s%s%s -ro %s %s", c_rpict,
					rep, opts, po, zopt, rawfile, oct1name);
			if (runcom(combuf))	/* run rpict/rpiece */
				goto rperror;
		} else {
			if (overture) {		/* run overture calculation */
				sprintf(combuf,
					"%s%s %s%s -x 64 -y 64 -ps 1 %s > %s",
						c_rpict, rep, vw, opts,
						oct1name, overfile);
				if (!do_rpiece || !next_process(0)) {
					if (runcom(combuf)) {
						fprintf(stderr,
					"%s: error in overture for view %s\n",
							progname, vs);
						quit(1);
					}
#ifndef NULL_DEVICE
					rmfile(overfile);
#endif
				} else if (do_rpiece)
					sleep(20);
			}
			if (do_rpiece) {
				sprintf(combuf, "%s -F %s %s%s %s %s%s%s -o %s %s",
						c_rpiece, sfile, rppopt, rep, vw,
						res, opts, po, rawfile, oct1name);
				while (next_process(1))
					sleep(10);
			} else {
				sprintf(combuf, "%s%s %s %s%s%s%s %s > %s",
						c_rpict, rep, vw, res, opts, po,
						zopt, oct1name, rawfile);
			}
			if ((pfile != NULL) & !do_rpiece) {
				if (!silent)
					printf("\t%s\n", combuf);
				fflush(stdout);
				sprintf(combuf, "%s%s %s %s%s%s %s > %s",
						c_rpict, rep, rppopt, res,
						opts, po, oct1name, rawfile);
				fp = popen(combuf, "w");
				if (fp == NULL)
					goto rperror;
				myprintview(vw, fp);
				if (pclose(fp))
					goto rperror;
			} else if (runcom(combuf))
				goto rperror;
		}
		if (do_rpiece) {		/* need to finish raw, first */
			finish_process();
			wait_process(1);
			/* XXX should check sync file to see if really done? */
		}
		if (!vdef(RAWFILE) || strcmp(vval(RAWFILE),vval(PICTURE))) {
						/* build pfilt command */
			if (do_rpiece)
				sprintf(combuf,
					"%s%s -x %d -y %d -p %.3f %s > %s",
					c_pfilt, pfopts, xres, yres, aspect,
					rawfile, picfile);
			else if (mult > 1)
				sprintf(combuf, "%s%s -x /%d -y /%d %s > %s",
					c_pfilt, pfopts, mult, mult,
					rawfile, picfile);
			else
				sprintf(combuf, "%s%s %s > %s", c_pfilt,
					pfopts, rawfile, picfile);
			if (runcom(combuf)) {	/* run pfilt */
				fprintf(stderr,
				"%s: error filtering view %s\n\t%s removed\n",
						progname, vs, picfile);
				unlink(picfile);
				quit(1);
			}
		}
						/* remove/rename raw file */
		if (vdef(RAWFILE)) {
			sprintf(combuf, "%s_%s.hdr", vval(RAWFILE), vs);
			mvfile(rawfile, combuf);
		} else
			rmfile(rawfile);
		if (do_rpiece)			/* done with sync file */
			rmfile(sfile);
		else
			finish_process();	/* exit if child */
	}
	wait_process(1);		/* wait for children to finish */
	if (pfile != NULL) {		/* clean up persistent rpict */
		RT_PID	pid;
		fp = fopen(pfile, "r");
		if (fp != NULL) {
			if (fscanf(fp, "%*s %d", &pid) != 1 ||
					kill(pid, 1) < 0)
				unlink(pfile);
			fclose(fp);
		}
	}
	return;
rperror:
	fprintf(stderr, "%s: error rendering view %s\n", progname, vs);
	quit(1);
#undef do_rpiece
}


static int
touch(			/* update a file */
	char	*fn
)
{
	if (!silent)
		printf("\ttouch %s\n", fn);
	if (!nprocs)
		return(0);
	return(setfdate(fn, time((time_t *)NULL)));
}


static int
runcom(			/* run command */
	char	*cs
)
{
	if (!silent)		/* echo it */
		printf("\t%s\n", cs);
	if (!nprocs)
		return(0);
	fflush(NULL);		/* flush output and pass to shell */
	return(system(cs));
}


static int
rmfile(			/* remove a file */
	char	*fn
)
{
	if (!silent)
		printf("\t%s %s\n", DELCMD, fn);
	if (!nprocs)
		return(0);
	return(unlink(fn));
}


static int
mvfile(		/* move a file */
	char	*fold,
	char	*fnew
)
{
	if (!silent)
		printf("\t%s %s %s\n", RENAMECMD, fold, fnew);
	if (!nprocs)
		return(0);
	return(rename(fold, fnew));
}


#ifdef RHAS_FORK_EXEC
static int
next_process(int reserve)		/* fork the next process */
{
	RT_PID	child_pid;

	if (nprocs <= 1)
		return(0);		/* it's us or no one */
	if (inchild()) {
		fprintf(stderr, "%s: internal error 1 in next_process()\n",
				progname);
		quit(1);
	}
	if (reserve > 0 && children_running >= nprocs-reserve)
		return(0);		/* caller holding back process(es) */
	if (children_running >= nprocs)
		wait_process(0);	/* wait for someone to finish */
	fflush(NULL);			/* flush output */
	child_pid = fork();		/* split process */
	if (child_pid == 0) {		/* we're the child */
		children_running = -1;
		nprocs = 1;
		return(0);
	}
	if (child_pid > 0) {		/* we're the parent */
		++children_running;
		return(1);
	}
	fprintf(stderr, "%s: warning -- fork() failed\n", progname);
	return(0);
}

static void
wait_process(			/* wait for process(es) to finish */
	int	all
)
{
	int	ourstatus = 0, status;
	RT_PID	pid;

	if (all)
		all = children_running;
	else if (children_running > 0)
		all = 1;
	while (all-- > 0) {
		pid = wait(&status);
		if (pid < 0)
			syserr(progname);
		status = status>>8 & 0xff;
		--children_running;
		if (status != 0) {	/* child's problem is our problem */
			if ((ourstatus == 0) & (children_running > 0))
				fprintf(stderr, "%s: waiting for remaining processes\n",
						progname);
			ourstatus = status;
			all = children_running;
		}
	}
	if (ourstatus != 0)
		quit(ourstatus);	/* bad status from child */
}
#else	/* ! RHAS_FORK_EXEC */
static int
next_process(int reserve)
{
	return(0);			/* cannot start new process */
}
static void
wait_process(all)
int	all;
{
	(void)all;			/* no one to wait for */
}
int
kill(pid, sig) /* win|unix_process.c should also wait and kill */
RT_PID pid;
int sig;
{
	return 0;
}
#endif	/* ! RHAS_FORK_EXEC */

static void
finish_process(void)			/* exit a child process */
{
	if (!inchild())
		return;			/* in parent -- noop */
	exit(0);
}

#ifdef _WIN32
setenv(vname, value)		/* set an environment variable */
char	*vname, *value;
{
	register char	*evp;

	evp = bmalloc(strlen(vname)+strlen(value)+2);
	if (evp == NULL)
		syserr(progname);
	sprintf(evp, "%s=%s", vname, value);
	if (putenv(evp) != 0) {
		fprintf(stderr, "%s: out of environment space\n", progname);
		quit(1);
	}
	if (!silent)
		printf("set %s\n", evp);
}
#endif


static void
badvalue(			/* report bad variable value and exit */
	int	vc
)
{
	fprintf(stderr, "%s: bad value for variable '%s'\n",
			progname, vnam(vc));
	quit(1);
}


static void
syserr(			/* report a system error and exit */
	char	*s
)
{
	perror(s);
	quit(1);
}


void
quit(ec)			/* exit program */
int	ec;
{
	exit(ec);
}
