#ifndef lint
static const char RCSid[] = "$Id: rpict.c,v 2.94 2020/07/20 15:52:30 greg Exp $";
#endif
/*
 *  rpict.c - routines and variables for picture generation.
 */

#include "copyright.h"

#include  <sys/types.h>

#include  "platform.h"
#ifdef NON_POSIX
 #ifdef MINGW
  #include  <sys/time.h>
 #endif
#else
 #ifdef BSD
  #include  <sys/time.h>
  #include  <sys/resource.h>
 #else
  #include  <sys/times.h>
  #include  <unistd.h>
 #endif
#endif

#include  <time.h>
#include  <signal.h>

#include  "ray.h"
#include  "paths.h"
#include  "ambient.h"
#include  "view.h"
#include  "random.h"
#include  "paths.h"
#include  "hilbert.h"
#include  "pmapbias.h"
#include  "pmapdiag.h"

#define	 RFTEMPLATE	"rfXXXXXX"

#ifndef SIGCONT
#ifdef SIGIO     /* XXX can we live without this? */
#define SIGCONT		SIGIO
#endif
#endif

CUBE  thescene;				/* our scene */
OBJECT	nsceneobjs;			/* number of objects in our scene */

int  dimlist[MAXDIM];			/* sampling dimensions */
int  ndims = 0;				/* number of sampling dimensions */
int  samplendx;				/* sample index number */

void  (*addobjnotify[])() = {ambnotify, NULL};

VIEW  ourview = STDVIEW;		/* view parameters */
int  hresolu = 512;			/* horizontal resolution */
int  vresolu = 512;			/* vertical resolution */
double	pixaspect = 1.0;		/* pixel aspect ratio */

int  psample = 4;			/* pixel sample size */
double	maxdiff = .05;			/* max. difference for interpolation */
double	dstrpix = 0.67;			/* square pixel distribution */

double  mblur = 0.;			/* motion blur parameter */

double  dblur = 0.;			/* depth-of-field blur parameter */

void  (*trace)() = NULL;		/* trace call */

int  do_irrad = 0;			/* compute irradiance? */

int  rand_samp = 0;			/* pure Monte Carlo sampling? */

double	dstrsrc = 0.0;			/* square source distribution */
double	shadthresh = .05;		/* shadow threshold */
double	shadcert = .5;			/* shadow certainty */
int  directrelay = 1;			/* number of source relays */
int  vspretest = 512;			/* virtual source pretest density */
int  directvis = 1;			/* sources visible? */
double	srcsizerat = .25;		/* maximum ratio source size/dist. */

COLOR  cextinction = BLKCOLOR;		/* global extinction coefficient */
COLOR  salbedo = BLKCOLOR;		/* global scattering albedo */
double  seccg = 0.;			/* global scattering eccentricity */
double  ssampdist = 0.;			/* scatter sampling distance */

double	specthresh = .15;		/* specular sampling threshold */
double	specjitter = 1.;		/* specular sampling jitter */

int  backvis = 1;			/* back face visibility */

int  maxdepth = 7;			/* maximum recursion depth */
double	minweight = 1e-3;		/* minimum ray weight */

char  *ambfile = NULL;			/* ambient file name */
COLOR  ambval = BLKCOLOR;		/* ambient value */
int  ambvwt = 0;			/* initial weight for ambient value */
double	ambacc = 0.2;			/* ambient accuracy */
int  ambres = 64;			/* ambient resolution */
int  ambdiv = 512;			/* ambient divisions */
int  ambssamp = 128;			/* ambient super-samples */
int  ambounce = 0;			/* ambient bounces */
char  *amblist[AMBLLEN];		/* ambient include/exclude list */
int  ambincl = -1;			/* include == 1, exclude == 0 */

int  ralrm = 0;				/* seconds between reports */

double	pctdone = 0.0;			/* percentage done */
time_t  tlastrept = 0L;			/* time at last report */
time_t  tstart;				/* starting time */

#define	 MAXDIV		16		/* maximum sample size */

#define	 pixjitter()	(.5+dstrpix*(.5-frandom()))

int  hres, vres;			/* resolution for this frame */

static VIEW	lastview;		/* the previous view input */

static void report(int);
static int nextview(FILE *fp);
static void render(char *zfile, char *oldfile);
static void fillscanline(COLOR *scanline, float *zline, char *sd, int xres,
		int y, int xstep);
static void fillscanbar(COLOR *scanbar[], float *zbar[], int xres,
		int y, int ysize);
static int fillsample(COLOR *colline, float *zline, int x, int y,
		int xlen, int ylen, int b);
static double pixvalue(COLOR  col, int  x, int  y);
static int salvage(char  *oldfile);
static int pixnumber(int  x, int  y, int  xres, int  yres);



void
quit(int code)			/* quit program */
{
	if (code)			/* report status */
		report(0);
#ifndef NON_POSIX
	headclean();			/* delete header file */
	pfclean();			/* clean up persist files */
#endif
	exit(code);
}


#ifndef NON_POSIX
static void
report(int dummy)		/* report progress */
{
	char			bcStat [128];
	double		u, s;
#ifdef BSD
	struct rusage	rubuf;
#else
	double		period = 1.0 / 60.0;
	struct tms	tbuf;
#endif

	tlastrept = time((time_t *)NULL);
#ifdef BSD
	getrusage(RUSAGE_SELF, &rubuf);
	u = rubuf.ru_utime.tv_sec + rubuf.ru_utime.tv_usec*1e-6;
	s = rubuf.ru_stime.tv_sec + rubuf.ru_stime.tv_usec*1e-6;
	getrusage(RUSAGE_CHILDREN, &rubuf);
	u += rubuf.ru_utime.tv_sec + rubuf.ru_utime.tv_usec*1e-6;
	s += rubuf.ru_stime.tv_sec + rubuf.ru_stime.tv_usec*1e-6;
#else
	times(&tbuf);
#ifdef _SC_CLK_TCK
	period = 1.0 / sysconf(_SC_CLK_TCK);
#endif
	u = ( tbuf.tms_utime + tbuf.tms_cutime ) * period;
	s = ( tbuf.tms_stime + tbuf.tms_cstime ) * period;
#endif

	/* PMAP: Get photon map bias compensation statistics */
	pmapBiasCompReport(bcStat);
	
	sprintf(errmsg,
			"%lu rays, %s%4.2f%% after %.3fu %.3fs %.3fr hours on %s (PID %d)\n",
			nrays, bcStat, pctdone, u*(1./3600.), s*(1./3600.),
			(tlastrept-tstart)*(1./3600.), myhostname(), getpid());
	eputs(errmsg);
#ifdef SIGCONT
	signal(SIGCONT, report);
#endif
}
#else
static void
report(int dummy)		/* report progress */
{
	char	bcStat [128];
	
	tlastrept = time((time_t *)NULL);

	/* PMAP: Get photon map bias compensation statistics */
	pmapBiasCompReport(bcStat);
	
	sprintf(errmsg, "%lu rays, %s%4.2f%% after %5.4f hours\n",
			nrays, bcStat, pctdone, (tlastrept-tstart)/3600.0);
	eputs(errmsg);
}
#endif


void
rpict(			/* generate image(s) */
	int  seq,
	char  *pout,
	char  *zout,
	char  *prvr
)
/*
 * If seq is greater than zero, then we will render a sequence of
 * images based on view parameter strings read from the standard input.
 * If pout is NULL, then all images will be sent to the standard ouput.
 * If seq is greater than zero and prvr is an integer, then it is the
 * frame number at which rendering should begin.  Preceeding view parameter
 * strings will be skipped in the input.
 * If pout and prvr are the same, prvr is renamed to avoid overwriting.
 * Note that pout and zout should contain %d format specifications for
 * sequenced file naming.
 */
{
	char  fbuf[128], fbuf2[128];
	int  npicts;
	char  *cp;
	RESOLU	rs;
	double	pa;
					/* check sampling */
	if (psample < 1)
		psample = 1;
	else if (psample > MAXDIV) {
		sprintf(errmsg, "pixel sampling reduced from %d to %d",
				psample, MAXDIV);
		error(WARNING, errmsg);
		psample = MAXDIV;
	}
					/* get starting frame */
	if (seq <= 0)
		seq = 0;
	else if (prvr != NULL && isint(prvr)) {
		int  rn;		/* skip to specified view */
		if ((rn = atoi(prvr)) < seq)
			error(USER, "recover frame less than start frame");
		if (pout == NULL)
			error(USER, "missing output file specification");
		for ( ; seq < rn; seq++)
			if (nextview(stdin) == EOF)
				error(USER, "unexpected EOF on view input");
		setview(&ourview);
		prvr = fbuf;			/* mark for renaming */
	}
	if ((pout != NULL) & (prvr != NULL)) {
		sprintf(fbuf, pout, seq);
		if (!strcmp(prvr, fbuf)) {	/* rename */
			strcpy(fbuf2, fbuf);
			for (cp = fbuf2; *cp; cp++)
				;
			while (cp > fbuf2 && !ISDIRSEP(cp[-1]))
				cp--;
			strcpy(cp, RFTEMPLATE);
			prvr = mktemp(fbuf2);
			if (rename(fbuf, prvr) < 0) {
				if (errno == ENOENT) {	/* ghost file */
					sprintf(errmsg,
						"new output file \"%s\"",
						fbuf);
					error(WARNING, errmsg);
					prvr = NULL;
				} else {		/* serious error */
					sprintf(errmsg,
					"cannot rename \"%s\" to \"%s\"",
						fbuf, prvr);
					error(SYSTEM, errmsg);
				}
			}
		}
	}
	npicts = 0;			/* render sequence */
	do {
		if (seq && nextview(stdin) == EOF)
			break;
		pctdone = 0.0;
		if (pout != NULL) {
			int	myfd;
			close(1);			/* reassign stdout */
			sprintf(fbuf, pout, seq);
		tryagain:
			errno = 0;			/* exclusive open */
			if ((myfd = open(fbuf, O_WRONLY|O_CREAT|O_EXCL, 0666)) < 0) {
				if ((errno != EEXIST) | (prvr != NULL) ||
						!strcmp(fbuf, pout)) {
					sprintf(errmsg,
						"cannot open output file \"%s\"",
						fbuf);
					error(SYSTEM, errmsg);
				}
				setview(&ourview);
				continue;		/* don't clobber */
			}
			if (myfd != 1) {
				unlink(fbuf);
				goto tryagain;		/* leave it open */
			}
			SET_FD_BINARY(1);
			dupheader();
		}
		hres = hresolu; vres = vresolu; pa = pixaspect;
		if (prvr != NULL) {
			if (viewfile(prvr, &ourview, &rs) <= 0) {
				sprintf(errmsg,
			"cannot recover view parameters from \"%s\"", prvr);
				error(WARNING, errmsg);
			} else {
				pa = 0.0;
				hres = scanlen(&rs);
				vres = numscans(&rs);
			}
		}
		if ((cp = setview(&ourview)) != NULL)
			error(USER, cp);
		normaspect(viewaspect(&ourview), &pa, &hres, &vres);
		if (seq) {
			if (ralrm > 0) {
				fprintf(stderr, "FRAME %d:", seq);
				fprintview(&ourview, stderr);
				putc('\n', stderr);
				fflush(stderr);
			}
			printf("FRAME=%d\n", seq);
		}
		fputs(VIEWSTR, stdout);
		fprintview(&ourview, stdout);
		putchar('\n');
		if ((pa < .99) | (pa > 1.01))
			fputaspect(pa, stdout);
		fputnow(stdout);
		fputformat(COLRFMT, stdout);
		putchar('\n');
		if (zout != NULL)
			sprintf(cp=fbuf, zout, seq);
		else
			cp = NULL;
		render(cp, prvr);
		prvr = NULL;
		npicts++;
	} while (seq++);
					/* check that we did something */
	if (npicts == 0)
		error(WARNING, "no output produced");
}


static int
nextview(				/* get next view from fp */
	FILE  *fp
)
{
	char  linebuf[256];

	lastview = ourview;
	while (fgets(linebuf, sizeof(linebuf), fp) != NULL)
		if (isview(linebuf) && sscanview(&ourview, linebuf) > 0)
			return(0);
	return(EOF);
}	


static void
render(				/* render the scene */
	char  *zfile,
	char  *oldfile
)
{
	COLOR  *scanbar[MAXDIV+1];	/* scanline arrays of pixel values */
	float  *zbar[MAXDIV+1];		/* z values */
	char  *sampdens;		/* previous sample density */
	int  ypos;			/* current scanline */
	int  ystep;			/* current y step size */
	int  hstep;			/* h step size */
	int  zfd;
	COLOR  *colptr;
	float  *zptr;
	int  i;
					/* check for empty image */
	if ((hres <= 0) | (vres <= 0)) {
		error(WARNING, "empty output picture");
		fprtresolu(0, 0, stdout);
		return;
	}
					/* allocate scanlines */
	for (i = 0; i <= psample; i++) {
		scanbar[i] = (COLOR *)malloc(hres*sizeof(COLOR));
		if (scanbar[i] == NULL)
			goto memerr;
	}
	hstep = (psample*140+49)/99;		/* quincunx sampling */
	ystep = (psample*99+70)/140;
	if (hstep > 2) {
		i = hres/hstep + 2;
		if ((sampdens = (char *)malloc(i)) == NULL)
			goto memerr;
		while (i--)
			sampdens[i] = hstep;
	} else
		sampdens = NULL;
					/* open z-file */
	if (zfile != NULL) {
		if ((zfd = open(zfile, O_WRONLY|O_CREAT, 0666)) == -1) {
			sprintf(errmsg, "cannot open z-file \"%s\"", zfile);
			error(SYSTEM, errmsg);
		}
		SET_FD_BINARY(zfd);
		for (i = 0; i <= psample; i++) {
			zbar[i] = (float *)malloc(hres*sizeof(float));
			if (zbar[i] == NULL)
				goto memerr;
		}
	} else {
		zfd = -1;
		for (i = 0; i <= psample; i++)
			zbar[i] = NULL;
	}
					/* write out boundaries */
	fprtresolu(hres, vres, stdout);
					/* recover file and compute first */
	i = salvage(oldfile);
	if (i >= vres)
		goto alldone;
	if ((zfd != -1) & (i > 0) &&
			lseek(zfd, (off_t)i*hres*sizeof(float), SEEK_SET) < 0)
		error(SYSTEM, "z-file seek error in render");
	pctdone = 100.0*i/vres;
	if (ralrm > 0)			/* report init stats */
		report(0);
#ifdef SIGCONT
	else
	signal(SIGCONT, report);
#endif
	ypos = vres-1 - i;			/* initialize sampling */
	if (directvis)
		init_drawsources(psample);
	fillscanline(scanbar[0], zbar[0], sampdens, hres, ypos, hstep);
						/* compute scanlines */
	for (ypos -= ystep; ypos > -ystep; ypos -= ystep) {
							/* bottom adjust? */
		if (ypos < 0) {
			ystep += ypos;
			ypos = 0;
		}
		colptr = scanbar[ystep];		/* move base to top */
		scanbar[ystep] = scanbar[0];
		scanbar[0] = colptr;
		zptr = zbar[ystep];
		zbar[ystep] = zbar[0];
		zbar[0] = zptr;
							/* fill base line */
		fillscanline(scanbar[0], zbar[0], sampdens,
				hres, ypos, hstep);
							/* fill bar */
		fillscanbar(scanbar, zbar, hres, ypos, ystep);
		if (directvis)				/* add bitty sources */
			drawsources(scanbar, zbar, 0, hres, ypos, ystep);
							/* write it out */
#ifdef SIGCONT
		signal(SIGCONT, SIG_IGN);	/* don't interrupt writes */
#endif
		for (i = ystep; i > 0; i--) {
			if (zfd != -1 && write(zfd, (char *)zbar[i],
					hres*sizeof(float))
					< hres*sizeof(float))
				goto writerr;
			if (fwritescan(scanbar[i], hres, stdout) < 0)
				goto writerr;
		}
		if (fflush(stdout) == EOF)
			goto writerr;
							/* record progress */
		pctdone = 100.0*(vres-1-ypos)/vres;
		if (ralrm > 0 && time((time_t *)NULL) >= tlastrept+ralrm)
			report(0);
#ifdef SIGCONT
		else
			signal(SIGCONT, report);
#endif
	}
						/* clean up */
#ifdef SIGCONT
	signal(SIGCONT, SIG_IGN);
#endif
	if (zfd != -1 && write(zfd, (char *)zbar[0], hres*sizeof(float))
				< hres*sizeof(float))
		goto writerr;
	fwritescan(scanbar[0], hres, stdout);
	if (fflush(stdout) == EOF)
		goto writerr;
alldone:
	if (zfd != -1) {
		if (close(zfd) == -1)
			goto writerr;
		for (i = 0; i <= psample; i++)
			free((void *)zbar[i]);
	}
	for (i = 0; i <= psample; i++)
		free((void *)scanbar[i]);
	if (sampdens != NULL)
		free(sampdens);
	pctdone = 100.0;
	if (ralrm > 0)
		report(0);
#ifdef SIGCONT
	signal(SIGCONT, SIG_DFL);
#endif
	return;
writerr:
	error(SYSTEM, "write error in render");
memerr:
	error(SYSTEM, "out of memory in render");
}


static void
fillscanline(	/* fill scan at y */
	COLOR	*scanline,
	float	*zline,
	char  *sd,
	int  xres,
	int  y,
	int  xstep
)
{
	static int  nc = 0;		/* number of calls */
	int  bl = xstep, b = xstep;
	double	z;
	int  i;

	z = pixvalue(scanline[0], 0, y);
	if (zline) zline[0] = z;
				/* zig-zag start for quincunx pattern */
	for (i = ++nc & 1 ? xstep : xstep/2; i < xres-1+xstep; i += xstep) {
		if (i >= xres) {
			xstep += xres-1-i;
			i = xres-1;
		}
		z = pixvalue(scanline[i], i, y);
		if (zline) zline[i] = z;
		if (sd) b = sd[0] > sd[1] ? sd[0] : sd[1];
		if (i <= xstep)
			b = fillsample(scanline, zline, 0, y, i, 0, b/2);
		else
			b = fillsample(scanline+i-xstep,
					zline ? zline+i-xstep : (float *)NULL,
					i-xstep, y, xstep, 0, b/2);
		if (sd) *sd++ = nc & 1 ? bl : b;
		bl = b;
	}
	if (sd && nc & 1) *sd = bl;
}


static void
fillscanbar(	/* fill interior */
	COLOR	*scanbar[],
	float	*zbar[],
	int  xres,
	int  y,
	int  ysize
)
{
	COLOR  vline[MAXDIV+1];
	float  zline[MAXDIV+1];
	int  b = ysize;
	int  i, j;

	for (i = 0; i < xres; i++) {
		copycolor(vline[0], scanbar[0][i]);
		copycolor(vline[ysize], scanbar[ysize][i]);
		if (zbar[0]) {
			zline[0] = zbar[0][i];
			zline[ysize] = zbar[ysize][i];
		}
		b = fillsample(vline, zbar[0] ? zline : (float *)NULL,
				i, y, 0, ysize, b/2);

		for (j = 1; j < ysize; j++)
			copycolor(scanbar[j][i], vline[j]);
		if (zbar[0])
			for (j = 1; j < ysize; j++)
				zbar[j][i] = zline[j];
	}
}


static int
fillsample( /* fill interior points */
	COLOR	*colline,
	float	*zline,
	int  x,
	int  y,
	int  xlen,
	int  ylen,
	int  b
)
{
	double	ratio;
	double	z;
	COLOR  ctmp;
	int  ncut;
	int  len;

	if (xlen > 0)			/* x or y length is zero */
		len = xlen;
	else
		len = ylen;

	if (len <= 1)			/* limit recursion */
		return(0);

	if (b > 0 ||
	(zline && 2.*fabs(zline[0]-zline[len]) > maxdiff*(zline[0]+zline[len]))
			|| bigdiff(colline[0], colline[len], maxdiff)) {

		z = pixvalue(colline[len>>1], x + (xlen>>1), y + (ylen>>1));
		if (zline) zline[len>>1] = z;
		ncut = 1;
	} else {					/* interpolate */
		copycolor(colline[len>>1], colline[len]);
		ratio = (double)(len>>1) / len;
		scalecolor(colline[len>>1], ratio);
		if (zline) zline[len>>1] = zline[len] * ratio;
		ratio = 1.0 - ratio;
		copycolor(ctmp, colline[0]);
		scalecolor(ctmp, ratio);
		addcolor(colline[len>>1], ctmp);
		if (zline) zline[len>>1] += zline[0] * ratio;
		ncut = 0;
	}
							/* recurse */
	ncut += fillsample(colline, zline, x, y, xlen>>1, ylen>>1, (b-1)/2);

	ncut += fillsample(colline+(len>>1),
			zline ? zline+(len>>1) : (float *)NULL,
			x+(xlen>>1), y+(ylen>>1),
			xlen-(xlen>>1), ylen-(ylen>>1), b/2);

	return(ncut);
}


static double
pixvalue(		/* compute pixel value */
	COLOR  col,			/* returned color */
	int  x,			/* pixel position */
	int  y
)
{
	extern void  SDsquare2disk(double ds[2], double seedx, double seedy);
	RAY  thisray;
	FVECT	lorg, ldir;
	double	hpos, vpos, vdist, lmax;
	int	i;
						/* compute view ray */
	setcolor(col, 0.0, 0.0, 0.0);
	hpos = (x+pixjitter())/hres;
	vpos = (y+pixjitter())/vres;
	if ((thisray.rmax = viewray(thisray.rorg, thisray.rdir,
					&ourview, hpos, vpos)) < -FTINY)
		return(0.0);

	vdist = ourview.vdist;
						/* set pixel index */
	samplendx = pixnumber(x,y,hres,vres);
						/* optional motion blur */
	if (lastview.type && mblur > FTINY && (lmax = viewray(lorg, ldir,
					&lastview, hpos, vpos)) >= -FTINY) {
		double  d = mblur*(.5-urand(281+samplendx));

		thisray.rmax = (1.-d)*thisray.rmax + d*lmax;
		for (i = 3; i--; ) {
			thisray.rorg[i] = (1.-d)*thisray.rorg[i] + d*lorg[i];
			thisray.rdir[i] = (1.-d)*thisray.rdir[i] + d*ldir[i];
		}
		if (normalize(thisray.rdir) == 0.0)
			return(0.0);
		vdist = (1.-d)*vdist + d*lastview.vdist;
	}
						/* optional depth-of-field */
	if (dblur > FTINY) {
		double  vc, df[2];
						/* random point on disk */
		SDsquare2disk(df, frandom(), frandom());
		df[0] *= .5*dblur;
		df[1] *= .5*dblur;
		if ((ourview.type == VT_PER) | (ourview.type == VT_PAR)) {
			double	adj = 1.0;
			if (ourview.type == VT_PER)
				adj /= DOT(thisray.rdir, ourview.vdir);
			df[0] /= sqrt(ourview.hn2);
			df[1] /= sqrt(ourview.vn2);
			for (i = 3; i--; ) {
				vc = ourview.vp[i] + adj*vdist*thisray.rdir[i];
				thisray.rorg[i] += df[0]*ourview.hvec[i] +
							df[1]*ourview.vvec[i] ;
				thisray.rdir[i] = vc - thisray.rorg[i];
			}
		} else {			/* non-standard view case */
			double	dfd = PI/4.*dblur*(.5 - frandom());
			if ((ourview.type != VT_ANG) & (ourview.type != VT_PLS)) {
				if (ourview.type != VT_CYL)
					df[0] /= sqrt(ourview.hn2);
				df[1] /= sqrt(ourview.vn2);
			}
			for (i = 3; i--; ) {
				vc = ourview.vp[i] + vdist*thisray.rdir[i];
				thisray.rorg[i] += df[0]*ourview.hvec[i] +
							df[1]*ourview.vvec[i] +
							dfd*ourview.vdir[i] ;
				thisray.rdir[i] = vc - thisray.rorg[i];
			}
		}
		if (normalize(thisray.rdir) == 0.0)
			return(0.0);
	}

	rayorigin(&thisray, PRIMARY, NULL, NULL);

	rayvalue(&thisray);			/* trace ray */

	copycolor(col, thisray.rcol);		/* return color */

	return(raydistance(&thisray));		/* return distance */
}


static int
salvage(		/* salvage scanlines from killed program */
	char  *oldfile
)
{
	COLR  *scanline;
	FILE  *fp;
	int  x, y;

	if (oldfile == NULL)
		goto gotzip;

	if ((fp = fopen(oldfile, "r")) == NULL) {
		sprintf(errmsg, "cannot open recover file \"%s\"", oldfile);
		error(WARNING, errmsg);
		goto gotzip;
	}
	SET_FILE_BINARY(fp);
				/* discard header */
	getheader(fp, NULL, NULL);
				/* get picture size */
	if (!fscnresolu(&x, &y, fp)) {
		sprintf(errmsg, "bad recover file \"%s\" - not removed",
				oldfile);
		error(WARNING, errmsg);
		fclose(fp);
		goto gotzip;
	}

	if ((x != hres) | (y != vres)) {
		sprintf(errmsg, "resolution mismatch in recover file \"%s\"",
				oldfile);
		error(USER, errmsg);
	}

	scanline = (COLR *)malloc(hres*sizeof(COLR));
	if (scanline == NULL)
		error(SYSTEM, "out of memory in salvage");
	for (y = 0; y < vres; y++) {
		if (freadcolrs(scanline, hres, fp) < 0)
			break;
		if (fwritecolrs(scanline, hres, stdout) < 0)
			goto writerr;
	}
	if (fflush(stdout) == EOF)
		goto writerr;
	free((void *)scanline);
	fclose(fp);
	unlink(oldfile);
	return(y);
gotzip:
	if (fflush(stdout) == EOF)
		error(SYSTEM, "error writing picture header");
	return(0);
writerr:
	sprintf(errmsg, "write error during recovery of \"%s\"", oldfile);
	error(SYSTEM, errmsg);
	return -1; /* pro forma return */
}

static int
pixnumber(		/* compute pixel index (screen door) */
	int  x,
	int  y,
	int  xres,
	int  yres
)
{
	unsigned	nbits = 0;
	bitmask_t	coord[2];

	if (xres < yres) xres = yres;
	while (xres > 0) {
		xres >>= 1;
		++nbits;
	}
	coord[0] = x; coord[1] = y;
	return ((int)hilbert_c2i(2, nbits, coord));
}
