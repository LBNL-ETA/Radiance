#ifndef lint
static const char RCSid[] = "$Id";
#endif
/*
 *  Radiance object animation program
 *
 *  Main program and control file handling.
 *
 *  See ranimove.h and the ranimove(1) man page for details.
 */

#include "copyright.h"

#include <time.h>
#ifdef _WIN32
  #include <winsock.h> /* struct timeval. XXX find a replacement? */
#else
  #include <sys/time.h>
#endif
#include <ctype.h>
#include <string.h>

#include "paths.h"
#include "ranimove.h"

int		NVARS = NV_INIT; /* total number of variables */

VARIABLE	vv[] = VV_INIT;	/* variable-value pairs */

extern int	nowarn;		/* don't report warnings? */
int		silent = 0;	/* run silently? */

int		quickstart = 0;	/* time initial frame as well? */

int		nprocs = 1;	/* number of rendering processes */

int		rtperfrm = 60;	/* seconds to spend per frame */

double		ndthresh = 2.;	/* noticeable difference threshold */
int		ndtset = 0;	/* user threshold -- stop when reached? */

int		fbeg = 1;	/* starting frame */
int		fend = 0;	/* ending frame */
int		fcur;		/* current frame being rendered */

char		lorendoptf[32];	/* low quality options file */
RAYPARAMS	lorendparams;	/* low quality rendering parameters */
char		hirendoptf[32];	/* high quality options file */
RAYPARAMS	hirendparams;	/* high quality rendering parameters */
RAYPARAMS	*curparams;	/* current parameter settings */
int		twolevels;	/* low and high quality differ */

double		mblur;		/* vflt(MBLUR) */
double		rate;		/* vflt(RATE) */

char		objtmpf[32];	/* object temporary file */

struct ObjMove	*obj_move;	/* object movements */

int		haveprio = 0;	/* high-level saliency specified */

int		gargc;		/* global argc for printargs */
char		**gargv;	/* global argv for printargs */

static void setdefaults(void);
static void setmove(struct ObjMove	*om, char	*ms);
static void setrendparams(char		*optf, char		*qval);
static void getradfile(char	*rfargs);
static void animate(void);
static int countviews(void); /* XXX duplicated function */
static char * getobjname(struct ObjMove	*om);
static char * getxf(struct ObjMove	*om, int	n);


int
main(
	int	argc,
	char	*argv[]
)
{
	int	explicate = 0;
	char	*cfname;
	int	i;

	progname = argv[0];			/* get arguments */
	gargc = argc;
	gargv = argv;
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 't':			/* seconds per frame */
			rtperfrm = atoi(argv[++i]);
			break;
		case 'd':			/* noticeable difference */
			ndthresh = atof(argv[++i]);
			ndtset = 1;
			break;
		case 'e':			/* print variables */
			explicate++;
			break;
		case 's':			/* silent running */
			silent++;
			break;
		case 'q':			/* start quickly */
			quickstart++;
			break;
		case 'w':			/* turn off warnings */
			nowarn++;
			break;
		case 'f':			/* frame range */
			switch (sscanf(argv[++i], "%d,%d", &fbeg, &fend)) {
			case 2:
				if ((fbeg <= 0) | (fend < fbeg))
					goto userr;
				break;
			case 1:
				if (fbeg <= 0)
					goto userr;
				fend = 0;
				break;
			default:
				goto userr;
			}
			break;
		case 'n':			/* number of processes */
			nprocs = atoi(argv[++i]);
			break;
		default:
			goto userr;
		}
	if (rtperfrm <= 0) {
		if (!ndtset)
			error(USER, "specify -d jnd with -t 0");
		rtperfrm = 7*24*3600;
	}
	if (i != argc-1)
		goto userr;
	cfname = argv[i];
						/* load variables */
	loadvars(cfname);
						/* check variables */
	checkvalues();
						/* load RIF if any */
	if (vdef(RIF))
		getradfile(vval(RIF));
						/* set defaults */
	setdefaults();
						/* print variables */
	if (explicate)
		printvars(stdout);
						/* run animation */
	if (nprocs > 0)
		animate();
						/* all done */
	if (lorendoptf[0])
		unlink(lorendoptf);
	if (hirendoptf[0] && strcmp(hirendoptf, lorendoptf))
		unlink(hirendoptf);
	if (objtmpf[0])
		unlink(objtmpf);
	return(0);
userr:
	fprintf(stderr,
"Usage: %s [-n nprocs][-f beg,end][-t sec][-d jnd][-s][-w][-e] anim_file\n",
			progname);
	quit(1);
	return 1; /* pro forma return */
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
setdefaults(void)			/* set default values */
{
	int	nviews;
	int	decades;
	char	buf[256];
	int	i;

	if (!vdef(OCTREEF)) {
		sprintf(errmsg, "%s or %s must be defined",
				vnam(OCTREEF), vnam(RIF));
		error(USER, errmsg);
	}
	if (!vdef(VIEWFILE)) {
		sprintf(errmsg, "%s must be defined", vnam(VIEWFILE));
		error(USER, errmsg);
	}
	nviews = countviews();
	if (!nviews)
		error(USER, "no views in view file");
	if (!vdef(END)) {
		if (nviews == 1) {
			sprintf(errmsg, "%s must be defined for single view",
					vnam(END));
			error(USER, errmsg);
		}
		sprintf(buf, "%d", nviews);
		vval(END) = savqstr(buf);
		vdef(END)++;
	}
	if (!fend)
		fend = vint(END);
	if (fbeg > fend)
		error(USER, "begin after end");
	if (!vdef(BASENAME)) {
		decades = (int)log10((double)vint(END)) + 1;
		if (decades < 3) decades = 3;
		sprintf(buf, "frame%%0%dd", decades);
		vval(BASENAME) = savqstr(buf);
		vdef(BASENAME)++;
	}
	if (!vdef(RATE)) {
		vval(RATE) = "8";
		vdef(RATE)++;
	}
	rate = vflt(RATE);
	if (!vdef(RESOLUTION)) {
		vval(RESOLUTION) = "640";
		vdef(RESOLUTION)++;
	}
	if (!vdef(MBLUR)) {
		vval(MBLUR) = "0";
		vdef(MBLUR)++;
	}
	mblur = vflt(MBLUR);
	if (mblur > 1.)
		mblur = 1.;
				/* set up objects */
	if (vdef(MOVE)) {
		obj_move = (struct ObjMove *)malloc(
				sizeof(struct ObjMove)*vdef(MOVE) );
		if (obj_move == NULL)
			error(SYSTEM, "out of memory in setdefaults");
		for (i = 0; i < vdef(MOVE); i++)
			setmove(&obj_move[i], nvalue(MOVE, i));
	}
				/* set up high quality options */
	setrendparams(hirendoptf, vval(HIGHQ));
	ray_save(&hirendparams);
				/* set up low quality options */
	setrendparams(lorendoptf, vval(LOWQ));
	ray_save(&lorendparams);
	curparams = &lorendparams;
	twolevels = memcmp(&lorendparams, &hirendparams, sizeof(RAYPARAMS));
}


static void
setmove(			/* assign a move object from spec. */
	struct ObjMove	*om,
	char	*ms
)
{
	char	parname[128];
	char	*cp;
	
	ms = nextword(parname, sizeof(parname), ms);
	if (!parname[0])
		goto badspec;
	for (cp = parname; *cp; cp++)
		if (isspace(*cp))
			*cp = '_';
	for (om->parent = (om - obj_move); om->parent--; )
		if (!strcmp(parname, obj_move[om->parent].name))
			break;
	if (om->parent < 0 &&
			strcmp(parname, ".") && strcmp(parname, VOIDID)) {
		sprintf(errmsg, "undefined parent object '%s'", parname);
		error(USER, errmsg);
	}
	ms = nextword(om->name, sizeof(om->name), ms);
	if (!om->name[0])
		goto badspec;
	for (cp = om->name; *cp; cp++)
		if (isspace(*cp))
			*cp = '_';
	ms = nextword(om->xf_file, sizeof(om->xf_file), ms);
	if (!strcmp(om->xf_file, "."))
		om->xf_file[0] = '\0';
	if (!om->xf_file[0]) {
		om->xfs[0] = '\0';
	} else if (om->xf_file[0] == '-') {
		strcpy(om->xfs, om->xf_file);
		om->xf_file[0] = '\0';
	}
	ms = nextword(om->spec, sizeof(om->spec), ms);
	if (om->spec[0]) {
		if (!strcmp(om->spec, ".") || !strcmp(om->spec, VOIDID))
			om->spec[0] = '\0';
		ms = nextword(om->prio_file, sizeof(om->prio_file), ms);
	} else
		om->prio_file[0] = '\0';
	if (om->prio_file[0]) {
		if (isflt(om->prio_file)) {
			om->prio = atof(om->prio_file);
			om->prio_file[0] = '\0';
			haveprio |= ((om->prio < 0.95) | (om->prio > 1.05));
		} else
			haveprio = 1;
	} else
		om->prio = 1.;
	om->cfm = -1;
	return;
badspec:
	error(USER, "bad object specification");
}


static void
setrendparams(	/* set global rendering parameters */
	char		*optf,
	char		*qval
)
{
	char	*argv[1024];
	char	**av = argv;
	int	ac, i, rval;

	av[ac=0] = NULL;
				/* load options from file, first */
	if (optf != NULL && *optf) {
		ac = wordfile(av, optf);
		if (ac < 0) {
			sprintf(errmsg, "cannot load options file \"%s\"",
					optf);
			error(SYSTEM, errmsg);
		}
	}
				/* then from options string */
	if (qval != NULL && qval[0] == '-')
		ac += wordstring(av+ac, qval);

				/* restore default parameters */
	ray_restore(NULL);
				/* set what we have */
	for (i = 0; i < ac; i++) {
		while ((rval = expandarg(&ac, &av, i)) > 0)
			;
		if (rval < 0) {
			sprintf(errmsg, "cannot expand '%s'", av[i]);
			error(SYSTEM, errmsg);
		}
		if (!strcmp(av[i], "-w")) {
			nowarn++;
			continue;
		}
		rval = getrenderopt(ac-i, av+i);
		if (rval < 0) {
			sprintf(errmsg, "bad render option at '%s'", av[i]);
			error(USER, errmsg);
		}
		i += rval;
	}
}


static void
getradfile(		/* run rad and get needed variables */
	char	*rfargs
)
{
	static short	mvar[] = {OCONV,OCTREEF,RESOLUTION,EXPOSURE,-1};
	char	combuf[256];
	register int	i;
	register char	*cp;
	char	*pippt = NULL;
					/* create rad command */
	strcpy(lorendoptf, "ranim0.opt");
	sprintf(combuf,
	"rad -v 0 -s -e -w %s %s oconv=-f OPTFILE=%s | egrep '^[ \t]*(NOMATCH",
			rfargs,
		(vdef(LOWQ) && vval(LOWQ)[0]!='-') ? vval(LOWQ) : "",
			lorendoptf);
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
		strcpy(cp, ")[ \t]*=' > ranimove.var");
		cp += 11;		/* point to file name */
	}
	system(combuf);			/* ignore exit code */
	if (pippt == NULL) {		/* load variables and remove file */
		loadvars(cp);
		unlink(cp);
	}
	if (!vdef(HIGHQ) || vval(HIGHQ)[0]=='-') {
		strcpy(hirendoptf, lorendoptf);
		return;
	}
					/* get high quality options */
	strcpy(hirendoptf, "ranim1.opt");
	sprintf(combuf, "rad -v 0 -s -w %s %s OPTFILE=%s",
			rfargs, vval(HIGHQ), hirendoptf);
	system(combuf);
}


static void
animate(void)			/* run through animation */
{
	int	rpass;

	if (fbeg > 1)				/* synchronize transforms */
		getoctspec(fbeg-1);

	for (fcur = fbeg; fcur <= fend; fcur++) {
		if (!silent)
			printf("Frame %d:\n", fcur);
						/* base rendering */
		init_frame();
						/* refinement */
		for (rpass = 0; refine_frame(rpass); rpass++)
			;
						/* final filter pass */
		filter_frame();
						/* output frame */
		send_frame();
	}
					/* free resources */
	free_frame();
	if (nprocs > 1)
		ray_pdone(1);
	else
		ray_done(1);
	getview(0);
	getexp(0);
}


extern VIEW *
getview(			/* get view number n */
	int	n
)
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
			return(NULL);			/* that's it */
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


extern char *
getexp(			/* get exposure for nth frame */
	int	n
)
{
	extern char	*fskip();
	static char	expval[32];
	static FILE	*expfp = NULL;
	static int	curfrm = 0;
	register char	*cp;

	if (n == 0) {				/* signal to close file */
		if (expfp != NULL) {
			fclose(expfp);
			expfp = NULL;
		}
		return(NULL);
	}
	if (!vdef(EXPOSURE))			/* no setting (auto) */
		return(NULL);
	if (isflt(vval(EXPOSURE)))		/* always the same */
		return(vval(EXPOSURE));
	if (expfp == NULL) {			/* open exposure file */
		if ((expfp = fopen(vval(EXPOSURE), "r")) == NULL) {
			sprintf(errmsg, "cannot open exposure file \"%s\"",
					vval(EXPOSURE));
			error(SYSTEM, errmsg);
		}
		curfrm = 0;
	}
	if (curfrm > n) {			/* rewind if necessary */
		rewind(expfp);
		curfrm = 0;
	}
	while (n > curfrm) {			/* read to line */
		if (fgets(expval, sizeof(expval), expfp) == NULL) {
			sprintf(errmsg, "%s: too few exposures",
					vval(EXPOSURE));
			error(USER, errmsg);
		}
		if (strlen(expval) == sizeof(expval)-1)
			goto formerr;
		curfrm++;
	}
	cp = fskip(expval);			/* check format */
	if (cp != NULL)
		while (isspace(*cp))
			*cp++ = '\0';
	if (cp == NULL || *cp)
		goto formerr;
	return(expval);				/* return value */
formerr:
	sprintf(errmsg, "%s: exposure format error on line %d",
			vval(EXPOSURE), curfrm);
	error(USER, errmsg);
	return NULL; /* pro forma return */
}


extern double
expspec_val(			/* get exposure value from spec. */
	char	*s
)
{
	double	expval;

	if (s == NULL || !*s)
		return(1.0);

	expval = atof(s);
	if ((s[0] == '+') | (s[0] == '-'))
		return(pow(2.0, expval));
	return(expval);
}


extern char *
getoctspec(			/* get octree for the given frame */
	int	n
)
{
	static char	combuf[1024];
	int		cfm = 0;
	int	uses_inline;
	FILE	*fp;
	int	i;
					/* is octree static? */
	if (!vdef(MOVE))
		return(vval(OCTREEF));
					/* done already */
	if (n == cfm)
		return(combuf);
					/* else create object file */
	strcpy(objtmpf, "movinobj.rad");
	fp = fopen(objtmpf, "w");
	if (fp == NULL) {
		sprintf(errmsg, "cannot write to moving objects file '%s'",
				objtmpf);
		error(SYSTEM, errmsg);
	}
	uses_inline = 0;
	for (i = 0; i < vdef(MOVE); i++) {
		int	inlc = (obj_move[i].spec[0] == '!');
		if (!obj_move[i].spec[0])
			continue;
		if (inlc)
			fprintf(fp, "%s %d \\\n\t| xform",
					obj_move[i].spec, n);
		else
			fputs("!xform", fp);
		fprintf(fp, " -n %s", getobjname(&obj_move[i]));
		fputs(getxf(&obj_move[i], n), fp);
		if (!inlc)
			fprintf(fp, " %s\n", obj_move[i].spec);
		else
			fputc('\n', fp);
		uses_inline |= inlc;
	}
	if (fclose(fp) == EOF)
		error(SYSTEM, "error writing moving objects file");
	if (uses_inline)
		sprintf(combuf, "!oconv %s -f -i '%s' %s",
				vdef(OCONV) ? vval(OCONV) : "",
				vval(OCTREEF), objtmpf);
	else
		sprintf(combuf, "!xform -f %s | oconv -f -i '%s' -",
				objtmpf, vval(OCTREEF));
	return(combuf);
}


static char *
getobjname(			/* get fully qualified object name */
	register struct ObjMove	*om
)
{
	static char	objName[512];
	register char	*cp = objName;
	
	strcpy(cp, om->name);
	while (om->parent >= 0) {
		while (*cp) cp++;
		*cp++ = '@';
		om = &obj_move[om->parent];
		strcpy(cp, om->name);
	}
	return(objName);
}


static char *
getxf(			/* get total transform for object */
	register struct ObjMove	*om,
	int	n
)
{
	static char	xfsbuf[4096];
	char		*xfp;
	int		framestep = 0;
	XF		oxf;
	FILE		*fp;
	char		abuf[512];
	char		*av[64];
	int		ac;
	int		i;
	register char	*cp;
					/* get parent transform, first */
	if (om->parent >= 0)
		xfp = getxf(&obj_move[om->parent], n);
	else
		*(xfp = xfsbuf + sizeof(xfsbuf)-1) = '\0';
					/* get transform spec. & priority */
	if (om->cfm != n) {
		if (om->xf_file[0]) {
			fp = fopen(om->xf_file, "r");
			if (fp == NULL) {
				sprintf(errmsg,
					"cannot open transform file '%s'",
						om->xf_file);
				error(SYSTEM, errmsg);
			}
			for (i = 0; i < n; i++)
				if (fgetline(om->xfs, sizeof(om->xfs), fp)
						== NULL) {
					sprintf(errmsg,
					"too few transforms in file '%s'",
							om->xf_file);
					error(USER, errmsg);
				}
			fclose(fp);
		}
		strcpy(cp=abuf, om->xfs);
		ac = 0;
		for ( ; ; ) {
			while (isspace(*cp))
				*cp++ = '\0';
			if (!*cp)
				break;
			av[ac++] = cp;
			while (*++cp && !isspace(*cp))
				;
		}
		av[ac] = NULL;
		if (xf(&oxf, ac, av) != ac ||
				fabs(oxf.sca) <= FTINY) {
			sprintf(errmsg, "bad transform args: %s",
					om->xfs);
			error(USER, errmsg);
		}
		copymat4(om->xfm, oxf.xfm);
		if (om->prio_file[0]) {
			fp = fopen(om->prio_file, "r");
			if (fp == NULL) {
				sprintf(errmsg,
					"cannot open priority file '%s'",
						om->prio_file);
				error(SYSTEM, errmsg);
			}
			for (i = 0; i < n; i++)
				if (fgets(abuf, sizeof(abuf), fp)
						== NULL) {
					sprintf(errmsg,
					"too few priorities in file '%s'",
							om->prio_file);
					error(USER, errmsg);
				}
			fclose(fp);
			cp = fskip(abuf);
			if (cp != NULL)
				while (isspace(*cp))
					*cp++ = '\0';
			if (cp == NULL || *cp) {
				sprintf(errmsg,
				"priority '%s' in file '%s' not a number",
						abuf, om->prio_file);
				error(USER, errmsg);
			}
			om->prio = atof(abuf);
		}
		framestep = (n == om->cfm + 1);
		om->cfm = n;
	}
					/* prepend to parent transform */
	if (om->xfs[0]) {
		i = strlen(om->xfs);
		if (xfp - i <= xfsbuf)
			error(INTERNAL, "transform too long in getxf");
		cp = om->xfs + i;
		while (i--)
			*--xfp = *--cp;
		*--xfp = ' ';
	}
	if (framestep)
		copymat4(oxf.xfm, om->cxfm);
	if (om->parent >= 0) {
		multmat4(om->cxfm, om->xfm, obj_move[om->parent].cxfm);
		om->cprio = obj_move[om->parent].cprio * om->prio;
	} else {
		copymat4(om->cxfm, om->xfm);
		om->cprio = om->prio;
	}
					/* XXX bxfm relies on call order */
	if (framestep) {
		if (invmat4(om->bxfm, om->cxfm))
			multmat4(om->bxfm, om->bxfm, oxf.xfm);
		else
			setident4(om->bxfm);
	}
					/* all done */
	return(xfp);
}


extern int
getmove(				/* find matching move object */
	OBJECT	obj
)
{
	static int	lasti;
	static OBJECT	lasto = OVOID;
	char	*onm, *objnm;
	int	len, len2;
	register int	i;

	if (obj == OVOID)
		return(-1);
	if (obj == lasto)
		return(lasti);
					/* look for matching object */
	onm = objptr(obj)->oname;
	for (i = vdef(MOVE); i--; ) {
		objnm = obj_move[i].name;
		len = strlen(objnm);
		if (!strncmp(onm, objnm, len)) {
			if ((obj_move[i].parent < 0) & (onm[len] == '.'))
				break;
			objnm = getobjname(&obj_move[i]) + len;
			len2 = strlen(objnm);
			if (!strncmp(onm+len, objnm, len2) && onm[len+len2] == '.')
				break;
		}
	}
	lasto = obj;			/* cache what we found */
	return(lasti = i);
}


extern double
obj_prio(			/* return priority for object */
	OBJECT	obj
)
{
	int	moi;
	
	if (obj == OVOID || (moi = getmove(obj)) < 0)
		return(1.0);
	return(obj_move[moi].cprio);
}


extern double
getTime(void)			/* get current time (CPU or real) */
{
	struct timeval	time_now;
					/* return CPU time if one process */
	if (nprocs == 1)
		return((double)clock()*(1.0/(double)CLOCKS_PER_SEC));
					/* otherwise, return wall time */
	gettimeofday(&time_now, NULL);
	return((double)time_now.tv_sec + 1e-6*(double)time_now.tv_usec);
}
