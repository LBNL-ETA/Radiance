#ifndef lint
static const char	RCSid[] = "$Id: pinterp.c,v 2.43 2012/02/08 13:43:21 greg Exp $";
#endif
/*
 * Interpolate and extrapolate pictures with different view parameters.
 *
 *	Greg Ward	09Dec89
 */

#include "copyright.h"

#include <ctype.h>
#include <string.h>

#include "platform.h"
#include "standard.h"
#include "rtprocess.h" /* Windows: must come before color.h */
#include "view.h"
#include "color.h"

#define LOG2		0.69314718055994530942

#define pscan(y)	(ourpict+(y)*hresolu)
#define sscan(y)	(ourspict+(y)*hresolu)
#define wscan(y)	(ourweigh+(y)*hresolu)
#define zscan(y)	(ourzbuf+(y)*hresolu)
#define bscan(y)	(ourbpict+(y)*hresolu)
#define averaging	(ourweigh != NULL)
#define blurring	(ourbpict != NULL)
#define usematrix	(hasmatrix & !averaging)
#define zisnorm		((!usematrix) | (ourview.type != VT_PER))

#define MAXWT		1000.		/* maximum pixel weight (averaging) */

#define F_FORE		1		/* fill foreground */
#define F_BACK		2		/* fill background */

#define PACKSIZ		256		/* max. calculation packet size */

#define RTCOM		"rtrace -h- -ovl -fff -ld- -i- -I- "

#define ABS(x)		((x)>0?(x):-(x))

struct position {int x,y; float z;};

#define NSTEPS		64		/* number steps in overlap prescan */
#define MINSTEP		4		/* minimum worthwhile preview step */

struct bound {int min,max;};

VIEW	ourview = STDVIEW;		/* desired view */
int	hresolu = 512;			/* horizontal resolution */
int	vresolu = 512;			/* vertical resolution */
double	pixaspect = 1.0;		/* pixel aspect ratio */

double	zeps = .02;			/* allowed z epsilon */

COLR	*ourpict;			/* output picture (COLR's) */
COLOR	*ourspict;			/* output pixel sums (averaging) */
float	*ourweigh = NULL;		/* output pixel weights (averaging) */
float	*ourzbuf;			/* corresponding z-buffer */
COLOR	*ourbpict = NULL;		/* blurred picture (view averaging) */

VIEW	avgview;			/* average view for -B option */
int	nvavg;				/* number of views averaged */

char	*progname;

int	fillo = F_FORE|F_BACK;		/* selected fill options */
int	fillsamp = 0;			/* sample separation (0 == inf) */
COLR	backcolr = BLKCOLR;		/* background color */
COLOR	backcolor = BLKCOLOR;		/* background color (float) */
double	backz = 0.0;			/* background z value */
int	normdist = 1;			/* i/o normalized distance? */
char	ourfmt[LPICFMT+1] = PICFMT;	/* original picture format */
double	ourexp = -1;			/* original picture exposure */
int	expadj = 0;			/* exposure adjustment (f-stops) */
double	rexpadj = 1;			/* real exposure adjustment */

VIEW	theirview;			/* input view */
int	gotview;			/* got input view? */
int	wrongformat = 0;		/* input in another format? */
RESOLU	tresolu;			/* input resolution */
double	theirexp;			/* input picture exposure */
MAT4	theirs2ours;			/* transformation matrix */
int	hasmatrix = 0;			/* has transformation matrix */

static SUBPROC PDesc = SP_INACTIVE; /* rtrace process descriptor */
unsigned short	queue[PACKSIZ][2];	/* pending pixels */
int	packsiz;			/* actual packet size */
int	queuesiz = 0;			/* number of pixels pending */

typedef void fillfunc_t(int x, int y);

static gethfunc headline;
static int nextview(FILE *fp);
static void compavgview(void);
static void addpicture(char *pfile, char *zspec);
static int pixform(MAT4 xfmat, VIEW *vw1, VIEW *vw2);
static void addscanline(struct bound *xl, int y,
	COLR *pline, float *zline, struct position *lasty);
static void addpixel(struct position *p0, struct position *p1,
	struct position *p2, COLR pix, double w, double z);
static double movepixel(FVECT pos);
static int getperim(struct bound *xl, struct bound *yl, float *zline, int zfd);
static void backpicture(fillfunc_t *fill, int samp);
static void fillpicture(fillfunc_t *fill);
static int clipaft(void);
static int addblur(void);
static void writepicture(void);
static void writedistance(char *fname);
static fillfunc_t backfill;
static fillfunc_t rcalfill;
static void calstart(char *prog, char *args);
static void caldone(void);
static void clearqueue(void);
static void syserror(char *s);

fillfunc_t *fillfunc = backfill;	/* selected fill function */

int
main(			/* interpolate pictures */
	int	argc,
	char	*argv[]
)
{
#define  check(ol,al)		if (argv[an][ol] || \
				badarg(argc-an-1,argv+an+1,al)) \
				goto badopt
	int	gotvfile = 0;
	int	doavg = -1;
	int	doblur = 0;
	char	*zfile = NULL;
	char	*expcomp = NULL;
	int	i, an, rval;

        SET_DEFAULT_BINARY();
        SET_FILE_BINARY(stdout);

	progname = argv[0];

	for (an = 1; an < argc && argv[an][0] == '-'; an++) {
		rval = getviewopt(&ourview, argc-an, argv+an);
		if (rval >= 0) {
			an += rval;
			continue;
		}
		switch (argv[an][1]) {
		case 'e':				/* exposure */
			check(2,"f");
			expcomp = argv[++an];
			break;
		case 't':				/* threshold */
			check(2,"f");
			zeps = atof(argv[++an]);
			break;
		case 'a':				/* average */
			check(2,NULL);
			doavg = 1;
			break;
		case 'B':				/* blur views */
			check(2,NULL);
			doblur = 1;
			break;
		case 'q':				/* quick (no avg.) */
			check(2,NULL);
			doavg = 0;
			break;
		case 'n':				/* dist. normalized? */
			check(2,NULL);
			normdist = !normdist;
			break;
		case 'f':				/* fill type */
			switch (argv[an][2]) {
			case '0':				/* none */
				check(3,NULL);
				fillo = 0;
				break;
			case 'f':				/* foreground */
				check(3,NULL);
				fillo = F_FORE;
				break;
			case 'b':				/* background */
				check(3,NULL);
				fillo = F_BACK;
				break;
			case 'a':				/* all */
				check(3,NULL);
				fillo = F_FORE|F_BACK;
				break;
			case 's':				/* sample */
				check(3,"i");
				fillsamp = atoi(argv[++an]);
				break;
			case 'c':				/* color */
				check(3,"fff");
				fillfunc = backfill;
				setcolor(backcolor, atof(argv[an+1]),
					atof(argv[an+2]), atof(argv[an+3]));
				setcolr(backcolr, colval(backcolor,RED),
						colval(backcolor,GRN),
						colval(backcolor,BLU));
				an += 3;
				break;
			case 'z':				/* z value */
				check(3,"f");
				fillfunc = backfill;
				backz = atof(argv[++an]);
				break;
			case 'r':				/* rtrace */
				check(3,"s");
				fillfunc = rcalfill;
				calstart(RTCOM, argv[++an]);
				break;
			default:
				goto badopt;
			}
			break;
		case 'z':				/* z file */
			check(2,"s");
			zfile = argv[++an];
			break;
		case 'x':				/* x resolution */
			check(2,"i");
			hresolu = atoi(argv[++an]);
			break;
		case 'y':				/* y resolution */
			check(2,"i");
			vresolu = atoi(argv[++an]);
			break;
		case 'p':				/* pixel aspect */
			if (argv[an][2] != 'a')
				goto badopt;
			check(3,"f");
			pixaspect = atof(argv[++an]);
			break;
		case 'v':				/* view file */
			if (argv[an][2] != 'f')
				goto badopt;
			check(3,"s");
			gotvfile = viewfile(argv[++an], &ourview, NULL);
			if (gotvfile < 0)
				syserror(argv[an]);
			else if (gotvfile == 0) {
				fprintf(stderr, "%s: bad view file\n",
						argv[an]);
				exit(1);
			}
			break;
		default:
		badopt:
			fprintf(stderr, "%s: command line error at '%s'\n",
					progname, argv[an]);
			goto userr;
		}
	}
						/* check arguments */
	if ((argc-an)%2)
		goto userr;
	if (fillsamp == 1)
		fillo &= ~F_BACK;
	if (doavg < 0)
		doavg = (argc-an) > 2;
	if (expcomp != NULL) {
		if ((expcomp[0] == '+') | (expcomp[0] == '-')) {
			expadj = atof(expcomp) + (expcomp[0]=='+' ? .5 : -.5);
			if (doavg | doblur)
				rexpadj = pow(2.0, atof(expcomp));
			else
				rexpadj = pow(2.0, (double)expadj);
		} else {
			if (!isflt(expcomp))
				goto userr;
			rexpadj = atof(expcomp);
			expadj = log(rexpadj)/LOG2 + (rexpadj>1 ? .5 : -.5);
			if (!(doavg | doblur))
				rexpadj = pow(2.0, (double)expadj);
		}
	}
						/* set view */
	if (nextview(doblur ? stdin : (FILE *)NULL) == EOF) {
		fprintf(stderr, "%s: no view on standard input!\n",
				progname);
		exit(1);
	}
	normaspect(viewaspect(&ourview), &pixaspect, &hresolu, &vresolu);
						/* allocate frame */
	if (doavg) {
		ourspict = (COLOR *)bmalloc(hresolu*vresolu*sizeof(COLOR));
		ourweigh = (float *)bmalloc(hresolu*vresolu*sizeof(float));
		if ((ourspict == NULL) | (ourweigh == NULL))
			syserror(progname);
	} else {
		ourpict = (COLR *)bmalloc(hresolu*vresolu*sizeof(COLR));
		if (ourpict == NULL)
			syserror(progname);
	}
	if (doblur) {
		ourbpict = (COLOR *)bmalloc(hresolu*vresolu*sizeof(COLOR));
		if (ourbpict == NULL)
			syserror(progname);
	}
	ourzbuf = (float *)bmalloc(hresolu*vresolu*sizeof(float));
	if (ourzbuf == NULL)
		syserror(progname);
							/* new header */
	newheader("RADIANCE", stdout);
	fputnow(stdout);
							/* run pictures */
	do {
		memset((char *)ourzbuf, '\0', hresolu*vresolu*sizeof(float));
		for (i = an; i < argc; i += 2)
			addpicture(argv[i], argv[i+1]);
		if (fillo&F_BACK)			/* fill in spaces */
			backpicture(fillfunc, fillsamp);
		else
			fillpicture(fillfunc);
							/* aft clipping */
		clipaft();
	} while (addblur() && nextview(stdin) != EOF);
							/* close calculation */
	caldone();
							/* add to header */
	printargs(argc, argv, stdout);
	compavgview();
	if (doblur | gotvfile) {
		fputs(VIEWSTR, stdout);
		fprintview(&avgview, stdout);
		putc('\n', stdout);
	}
	if ((pixaspect < .99) | (pixaspect > 1.01))
		fputaspect(pixaspect, stdout);
	if (ourexp > 0)
		ourexp *= rexpadj;
	else
		ourexp = rexpadj;
	if ((ourexp < .995) | (ourexp > 1.005))
		fputexpos(ourexp, stdout);
	if (strcmp(ourfmt, PICFMT))		/* print format if known */
		fputformat(ourfmt, stdout);
	putc('\n', stdout);
							/* write picture */
	writepicture();
							/* write z file */
	if (zfile != NULL)
		writedistance(zfile);

	exit(0);
userr:
	fprintf(stderr,
	"Usage: %s [view opts][-t eps][-z zout][-e spec][-B][-a|-q][-fT][-n] pfile zspec ..\n",
			progname);
	exit(1);
#undef check
}


static int
headline(				/* process header string */
	char	*s,
	void	*p
)
{
	char	fmt[32];

	if (isheadid(s))
		return(0);
	if (formatval(fmt, s)) {
		if (globmatch(ourfmt, fmt)) {
			wrongformat = 0;
			strcpy(ourfmt, fmt);
		} else
			wrongformat = 1;
		return(0);
	}
	if (nvavg < 2) {
		putc('\t', stdout);
		fputs(s, stdout);
	}
	if (isexpos(s)) {
		theirexp *= exposval(s);
		return(0);
	}
	if (isview(s) && sscanview(&theirview, s) > 0)
		gotview++;
	return(0);
}


static int
nextview(				/* get and set next view */
	FILE	*fp
)
{
	char	linebuf[256];
	char	*err;
	register int	i;

	if (fp != NULL) {
		do			/* get new view */
			if (fgets(linebuf, sizeof(linebuf), fp) == NULL)
				return(EOF);
		while (!isview(linebuf) || !sscanview(&ourview, linebuf));
	}
					/* set new view */
	if ((err = setview(&ourview)) != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	if (!nvavg) {			/* first view */
		avgview = ourview;
		return(nvavg++);
	}
					/* add to average view */
	for (i = 0; i < 3; i++) {
		avgview.vp[i] += ourview.vp[i];
		avgview.vdir[i] += ourview.vdir[i];
		avgview.vup[i] += ourview.vup[i];
	}
	avgview.vdist += ourview.vdist;
	avgview.horiz += ourview.horiz;
	avgview.vert += ourview.vert;
	avgview.hoff += ourview.hoff;
	avgview.voff += ourview.voff;
	avgview.vfore += ourview.vfore;
	avgview.vaft += ourview.vaft;
	return(nvavg++);
}


static void
compavgview(void)				/* compute average view */
{
	register int	i;
	double	f;

	if (nvavg < 2)
		return;
	f = 1.0/nvavg;
	for (i = 0; i < 3; i++) {
		avgview.vp[i] *= f;
		avgview.vdir[i] *= f;
		avgview.vup[i] *= f;
	}
	avgview.vdist *= f;
	avgview.horiz *= f;
	avgview.vert *= f;
	avgview.hoff *= f;
	avgview.voff *= f;
	avgview.vfore *= f;
	avgview.vaft *= f;
	if (setview(&avgview) != NULL)		/* in case of emergency... */
		avgview = ourview;
	pixaspect = viewaspect(&avgview) * hresolu / vresolu;
}


static void
addpicture(		/* add picture to output */
	char	*pfile,
	char	*zspec
)
{
	FILE	*pfp;
	int	zfd;
	char	*err;
	COLR	*scanin;
	float	*zin;
	struct position	*plast;
	struct bound	*xlim, ylim;
	int	y;
					/* open picture file */
	if ((pfp = fopen(pfile, "r")) == NULL)
		syserror(pfile);
					/* get header with exposure and view */
	theirexp = 1.0;
	theirview = stdview;
	gotview = 0;
	if (nvavg < 2)
		printf("%s:\n", pfile);
	getheader(pfp, headline, NULL);
	if (wrongformat || !gotview || !fgetsresolu(&tresolu, pfp)) {
		fprintf(stderr, "%s: picture format error\n", pfile);
		exit(1);
	}
	if (ourexp <= 0)
		ourexp = theirexp;
	else if (ABS(theirexp-ourexp) > .01*ourexp)
		fprintf(stderr, "%s: different exposure (warning)\n", pfile);
	if ( (err = setview(&theirview)) ) {
		fprintf(stderr, "%s: %s\n", pfile, err);
		exit(1);
	}
					/* compute transformation */
	hasmatrix = pixform(theirs2ours, &theirview, &ourview);
					/* get z specification or file */
	zin = (float *)malloc(scanlen(&tresolu)*sizeof(float));
	if (zin == NULL)
		syserror(progname);
	if ((zfd = open(zspec, O_RDONLY)) == -1) {
		double	zvalue;
		register int	x;
		if (!isflt(zspec) || (zvalue = atof(zspec)) <= 0.0)
			syserror(zspec);
		for (x = scanlen(&tresolu); x-- > 0; )
			zin[x] = zvalue;
	}
					/* compute transferrable perimeter */
	xlim = (struct bound *)malloc(numscans(&tresolu)*sizeof(struct bound));
	if (xlim == NULL)
		syserror(progname);
	if (!getperim(xlim, &ylim, zin, zfd)) {	/* overlapping area? */
		free((void *)zin);
		free((void *)xlim);
		if (zfd != -1)
			close(zfd);
		fclose(pfp);
		return;
	}
					/* allocate scanlines */
	scanin = (COLR *)malloc(scanlen(&tresolu)*sizeof(COLR));
	plast = (struct position *)calloc(scanlen(&tresolu),
			sizeof(struct position));
	if ((scanin == NULL) | (plast == NULL))
		syserror(progname);
					/* skip to starting point */
	for (y = 0; y < ylim.min; y++)
		if (freadcolrs(scanin, scanlen(&tresolu), pfp) < 0) {
			fprintf(stderr, "%s: read error\n", pfile);
			exit(1);
		}
	if (zfd != -1 && lseek(zfd,
			(off_t)ylim.min*scanlen(&tresolu)*sizeof(float),
			SEEK_SET) < 0)
		syserror(zspec);
					/* load image */
	for (y = ylim.min; y <= ylim.max; y++) {
		if (freadcolrs(scanin, scanlen(&tresolu), pfp) < 0) {
			fprintf(stderr, "%s: read error\n", pfile);
			exit(1);
		}
		if (zfd != -1 && read(zfd, (char *)zin,
				scanlen(&tresolu)*sizeof(float))
				< scanlen(&tresolu)*sizeof(float))
			syserror(zspec);
		addscanline(xlim+y, y, scanin, zin, plast);
	}
					/* clean up */
	free((void *)xlim);
	free((void *)scanin);
	free((void *)zin);
	free((void *)plast);
	fclose(pfp);
	if (zfd != -1)
		close(zfd);
}


static int
pixform(		/* compute view1 to view2 matrix */
	register MAT4	xfmat,
	register VIEW	*vw1,
	register VIEW	*vw2
)
{
	double	m4t[4][4];

	if ((vw1->type != VT_PER) & (vw1->type != VT_PAR))
		return(0);
	if ((vw2->type != VT_PER) & (vw2->type != VT_PAR))
		return(0);
	setident4(xfmat);
	xfmat[0][0] = vw1->hvec[0];
	xfmat[0][1] = vw1->hvec[1];
	xfmat[0][2] = vw1->hvec[2];
	xfmat[1][0] = vw1->vvec[0];
	xfmat[1][1] = vw1->vvec[1];
	xfmat[1][2] = vw1->vvec[2];
	xfmat[2][0] = vw1->vdir[0];
	xfmat[2][1] = vw1->vdir[1];
	xfmat[2][2] = vw1->vdir[2];
	xfmat[3][0] = vw1->vp[0];
	xfmat[3][1] = vw1->vp[1];
	xfmat[3][2] = vw1->vp[2];
	setident4(m4t);
	m4t[0][0] = vw2->hvec[0]/vw2->hn2;
	m4t[1][0] = vw2->hvec[1]/vw2->hn2;
	m4t[2][0] = vw2->hvec[2]/vw2->hn2;
	m4t[3][0] = -DOT(vw2->vp,vw2->hvec)/vw2->hn2;
	m4t[0][1] = vw2->vvec[0]/vw2->vn2;
	m4t[1][1] = vw2->vvec[1]/vw2->vn2;
	m4t[2][1] = vw2->vvec[2]/vw2->vn2;
	m4t[3][1] = -DOT(vw2->vp,vw2->vvec)/vw2->vn2;
	m4t[0][2] = vw2->vdir[0];
	m4t[1][2] = vw2->vdir[1];
	m4t[2][2] = vw2->vdir[2];
	m4t[3][2] = -DOT(vw2->vp,vw2->vdir);
	multmat4(xfmat, xfmat, m4t);
	return(1);
}


static void
addscanline(	/* add scanline to output */
	struct bound	*xl,
	int	y,
	COLR	*pline,
	float	*zline,
	struct position	*lasty		/* input/output */
)
{
	FVECT	pos;
	struct position	lastx, newpos;
	double	wt;
	register int	x;

	lastx.z = 0;
	for (x = xl->max; x >= xl->min; x--) {
		pix2loc(pos, &tresolu, x, y);
		pos[2] = zline[x];
		if ((wt = movepixel(pos)) <= FTINY) {
			lasty[x].z = lastx.z = 0;	/* mark invalid */
			continue;
		}
					/* add pixel to our image */
		newpos.x = pos[0] * hresolu;
		newpos.y = pos[1] * vresolu;
		newpos.z = zline[x];
		addpixel(&newpos, &lastx, &lasty[x], pline[x], wt, pos[2]);
		lasty[x].x = lastx.x = newpos.x;
		lasty[x].y = lastx.y = newpos.y;
		lasty[x].z = lastx.z = newpos.z;
	}
}


static void
addpixel(		/* fill in pixel parallelogram */
	struct position	*p0,
	struct position	*p1,
	struct position	*p2,
	COLR	pix,
	double	w,
	double	z
)
{
	double	zt = 2.*zeps*p0->z;		/* threshold */
	COLOR	pval;				/* converted+weighted pixel */
	int	s1x, s1y, s2x, s2y;		/* step sizes */
	int	l1, l2, c1, c2;			/* side lengths and counters */
	int	p1isy;				/* p0p1 along y? */
	int	x1, y1;				/* p1 position */
	register int	x, y;			/* final position */

					/* compute vector p0p1 */
	if (fillo&F_FORE && ABS(p1->z-p0->z) <= zt) {
		s1x = p1->x - p0->x;
		s1y = p1->y - p0->y;
		l1 = ABS(s1x);
		if ( (p1isy = (ABS(s1y) > l1)) )
			l1 = ABS(s1y);
		else if (l1 < 1)
			l1 = 1;
	} else {
		l1 = s1x = s1y = 1;
		p1isy = -1;
	}
					/* compute vector p0p2 */
	if (fillo&F_FORE && ABS(p2->z-p0->z) <= zt) {
		s2x = p2->x - p0->x;
		s2y = p2->y - p0->y;
		if (p1isy == 1)
			l2 = ABS(s2x);
		else {
			l2 = ABS(s2y);
			if (p1isy != 0 && ABS(s2x) > l2)
				l2 = ABS(s2x);
		}
		if (l2 < 1)
			l2 = 1;
	} else
		l2 = s2x = s2y = 1;
					/* fill the parallelogram */
	if (averaging) {
		colr_color(pval, pix);
		scalecolor(pval, w);
	}
	for (c1 = l1; c1-- > 0; ) {
		x1 = p0->x + c1*s1x/l1;
		y1 = p0->y + c1*s1y/l1;
		for (c2 = l2; c2-- > 0; ) {
			x = x1 + c2*s2x/l2;
			if ((x < 0) | (x >= hresolu))
				continue;
			y = y1 + c2*s2y/l2;
			if ((y < 0) | (y >= vresolu))
				continue;
			if (averaging) {
				if (zscan(y)[x] <= 0 || zscan(y)[x]-z
						> zeps*zscan(y)[x]) {
					copycolor(sscan(y)[x], pval);
					wscan(y)[x] = w;
					zscan(y)[x] = z;
				} else if (z-zscan(y)[x] <= zeps*zscan(y)[x]) {
					addcolor(sscan(y)[x], pval);
					wscan(y)[x] += w;
				}
			} else if (zscan(y)[x] <= 0 || zscan(y)[x]-z
						> zeps*zscan(y)[x]) {
				copycolr(pscan(y)[x], pix);
				zscan(y)[x] = z;
			}
		}
	}
}


static double
movepixel(				/* reposition image point */
	register FVECT	pos
)
{
	FVECT	pt, tdir, odir;
	double	d;
	
	if (pos[2] <= 0)		/* empty pixel */
		return(0);
	if (usematrix) {
		pos[0] += theirview.hoff - .5;
		pos[1] += theirview.voff - .5;
		if (normdist & (theirview.type == VT_PER))
			d = sqrt(1. + pos[0]*pos[0]*theirview.hn2
					+ pos[1]*pos[1]*theirview.vn2);
		else
			d = 1.;
		pos[2] += d*theirview.vfore;
		if (theirview.type == VT_PER) {
			pos[2] /= d;
			pos[0] *= pos[2];
			pos[1] *= pos[2];
		}
		multp3(pos, pos, theirs2ours);
		if (pos[2] <= ourview.vfore)
			return(0);
		if (ourview.type == VT_PER) {
			pos[0] /= pos[2];
			pos[1] /= pos[2];
		}
		pos[0] += .5 - ourview.hoff;
		pos[1] += .5 - ourview.voff;
		pos[2] -= ourview.vfore;
	} else {
		if (viewray(pt, tdir, &theirview, pos[0], pos[1]) < -FTINY)
			return(0);
		if ((!normdist) & (theirview.type == VT_PER))	/* adjust */
			pos[2] *= sqrt(1. + pos[0]*pos[0]*theirview.hn2
					+ pos[1]*pos[1]*theirview.vn2);
		pt[0] += tdir[0]*pos[2];
		pt[1] += tdir[1]*pos[2];
		pt[2] += tdir[2]*pos[2];
		viewloc(pos, &ourview, pt);
		if (pos[2] <= 0)
			return(0);
	}
	if ((pos[0] < 0) | (pos[0] >= 1-FTINY) | (pos[1] < 0) | (pos[1] >= 1-FTINY))
		return(0);
	if (!averaging)
		return(1);
						/* compute pixel weight */
	if (ourview.type == VT_PAR) {
		d = DOT(ourview.vdir,tdir);
		d = 1. - d*d;
	} else {
		VSUB(odir, pt, ourview.vp);
		d = DOT(odir,tdir);
		d = 1. - d*d/DOT(odir,odir);
	}
	if (d <= 1./MAXWT/MAXWT)
		return(MAXWT);		/* clip to maximum weight */
	return(1./sqrt(d));
}


static int
getperim(		/* compute overlapping image area */
	register struct bound	*xl,
	struct bound	*yl,
	float	*zline,
	int	zfd
)
{
	int	step;
	FVECT	pos;
	register int	x, y;
						/* set up step size */
	if (scanlen(&tresolu) < numscans(&tresolu))
		step = scanlen(&tresolu)/NSTEPS;
	else
		step = numscans(&tresolu)/NSTEPS;
	if (step < MINSTEP) {			/* not worth cropping? */
		yl->min = 0;
		yl->max = numscans(&tresolu) - 1;
		x = scanlen(&tresolu) - 1;
		for (y = numscans(&tresolu); y--; ) {
			xl[y].min = 0;
			xl[y].max = x;
		}
		return(1);
	}
	yl->min = 32000; yl->max = 0;		/* search for points on image */
	for (y = step - 1; y < numscans(&tresolu); y += step) {
		if (zfd != -1) {
			if (lseek(zfd, (off_t)y*scanlen(&tresolu)*sizeof(float),
					SEEK_SET) < 0)
				syserror("lseek");
			if (read(zfd, (char *)zline,
					scanlen(&tresolu)*sizeof(float))
					< scanlen(&tresolu)*sizeof(float))
				syserror("read");
		}
		xl[y].min = 32000; xl[y].max = 0;		/* x max */
		for (x = scanlen(&tresolu); (x -= step) > 0; ) {
			pix2loc(pos, &tresolu, x, y);
			pos[2] = zline[x];
			if (movepixel(pos) > FTINY) {
				xl[y].max = x + step - 1;
				xl[y].min = x - step + 1;	/* x min */
				if (xl[y].min < 0)
					xl[y].min = 0;
				for (x = step - 1; x < xl[y].max; x += step) {
					pix2loc(pos, &tresolu, x, y);
					pos[2] = zline[x];
					if (movepixel(pos) > FTINY) {
						xl[y].min = x - step + 1;
						break;
					}
				}
				if (y < yl->min)		/* y limits */
					yl->min = y - step + 1;
				yl->max = y + step - 1;
				break;
			}
		}
							/* fill in between */
		if (y < step) {
			xl[y-1].min = xl[y].min;
			xl[y-1].max = xl[y].max;
		} else {
			if (xl[y].min < xl[y-step].min)
				xl[y-1].min = xl[y].min;
			else
				xl[y-1].min = xl[y-step].min;
			if (xl[y].max > xl[y-step].max)
				xl[y-1].max = xl[y].max;
			else
				xl[y-1].max = xl[y-step].max;
		}
		for (x = 2; x < step; x++)
			*(xl+y-x) = *(xl+y-1);
	}
	if (yl->max >= numscans(&tresolu))
		yl->max = numscans(&tresolu) - 1;
	y -= step;
	for (x = numscans(&tresolu) - 1; x > y; x--)	/* fill bottom rows */
		*(xl+x) = *(xl+y);
	return(yl->max >= yl->min);
}


static void
backpicture(			/* background fill algorithm */
	fillfunc_t *fill,
	int	samp
)
{
	int	*yback, xback;
	int	y;
	register int	x, i;
							/* get back buffer */
	yback = (int *)malloc(hresolu*sizeof(int));
	if (yback == NULL)
		syserror(progname);
	for (x = 0; x < hresolu; x++)
		yback[x] = -2;
	/*
	 * Xback and yback are the pixel locations of suitable
	 * background values in each direction.
	 * A value of -2 means unassigned, and -1 means
	 * that there is no suitable background in this direction.
	 */
							/* fill image */
	for (y = 0; y < vresolu; y++) {
		xback = -2;
		for (x = 0; x < hresolu; x++)
			if (zscan(y)[x] <= 0) {		/* empty pixel */
				/*
				 * First, find background from above or below.
				 * (farthest assigned pixel)
				 */
				if (yback[x] == -2) {
					for (i = y+1; i < vresolu; i++)
						if (zscan(i)[x] > 0)
							break;
					if (i < vresolu
				&& (y <= 0 || zscan(y-1)[x] < zscan(i)[x]))
						yback[x] = i;
					else
						yback[x] = y-1;
				}
				/*
				 * Next, find background from left or right.
				 */
				if (xback == -2) {
					for (i = x+1; i < hresolu; i++)
						if (zscan(y)[i] > 0)
							break;
					if (i < hresolu
				&& (x <= 0 || zscan(y)[x-1] < zscan(y)[i]))
						xback = i;
					else
						xback = x-1;
				}
				/*
				 * If we have no background for this pixel,
				 * use the given fill function.
				 */
				if (xback < 0 && yback[x] < 0)
					goto fillit;
				/*
				 * Compare, and use the background that is
				 * farther, unless one of them is next to us.
				 * If the background is too distant, call
				 * the fill function.
				 */
				if ( yback[x] < 0
					|| (xback >= 0 && ABS(x-xback) <= 1)
					|| ( ABS(y-yback[x]) > 1
						&& zscan(yback[x])[x]
						< zscan(y)[xback] ) ) {
					if (samp > 0 && ABS(x-xback) >= samp)
						goto fillit;
					if (averaging) {
						copycolor(sscan(y)[x],
							sscan(y)[xback]);
						wscan(y)[x] = wscan(y)[xback];
					} else
						copycolr(pscan(y)[x],
							pscan(y)[xback]);
					zscan(y)[x] = zscan(y)[xback];
				} else {
					if (samp > 0 && ABS(y-yback[x]) > samp)
						goto fillit;
					if (averaging) {
						copycolor(sscan(y)[x],
							sscan(yback[x])[x]);
						wscan(y)[x] =
							wscan(yback[x])[x];
					} else
						copycolr(pscan(y)[x],
							pscan(yback[x])[x]);
					zscan(y)[x] = zscan(yback[x])[x];
				}
				continue;
			fillit:
				(*fill)(x,y);
				if (fill == rcalfill) {		/* use it */
					clearqueue();
					xback = x;
					yback[x] = y;
				}
			} else {				/* full pixel */
				yback[x] = -2;
				xback = -2;
			}
	}
	free((void *)yback);
}


static void
fillpicture(		/* paint in empty pixels using fill */
	fillfunc_t *fill
)
{
	register int	x, y;

	for (y = 0; y < vresolu; y++)
		for (x = 0; x < hresolu; x++)
			if (zscan(y)[x] <= 0)
				(*fill)(x,y);
	if (fill == rcalfill)
		clearqueue();
}


static int
clipaft(void)			/* perform aft clipping as indicated */
{
	register int	x, y;
	int	adjtest = (ourview.type == VT_PER) & zisnorm;
	double	tstdist;
	double	yzn2, vx;

	if (ourview.vaft <= FTINY)
		return(0);
	tstdist = ourview.vaft - ourview.vfore;
	for (y = 0; y < vresolu; y++) {
		if (adjtest) {				/* adjust test */
			yzn2 = (y+.5)/vresolu + ourview.voff - .5;
			yzn2 = 1. + yzn2*yzn2*ourview.vn2;
			tstdist = (ourview.vaft - ourview.vfore)*sqrt(yzn2);
		}
		for (x = 0; x < hresolu; x++)
			if (zscan(y)[x] > tstdist) {
				if (adjtest) {
					vx = (x+.5)/hresolu + ourview.hoff - .5;
					if (zscan(y)[x] <= (ourview.vaft -
							ourview.vfore) *
						sqrt(vx*vx*ourview.hn2 + yzn2))
						continue;
				}
				if (averaging)
					memset(sscan(y)[x], '\0', sizeof(COLOR));
				else
					memset(pscan(y)[x], '\0', sizeof(COLR));
				zscan(y)[x] = 0.0;
			}
	}
	return(1);
}


static int
addblur(void)				/* add to blurred picture */
{
	COLOR	cval;
	double	d;
	register int	i;

	if (!blurring)
		return(0);
	i = hresolu*vresolu;
	if (nvavg < 2)
		if (averaging)
			while (i--) {
				copycolor(ourbpict[i], ourspict[i]);
				d = 1.0/ourweigh[i];
				scalecolor(ourbpict[i], d);
			}
		else
			while (i--)
				colr_color(ourbpict[i], ourpict[i]);
	else
		if (averaging)
			while (i--) {
				copycolor(cval, ourspict[i]);
				d = 1.0/ourweigh[i];
				scalecolor(cval, d);
				addcolor(ourbpict[i], cval);
			}
		else
			while (i--) {
				colr_color(cval, ourpict[i]);
				addcolor(ourbpict[i], cval);
			}
				/* print view */
	printf("VIEW%d:", nvavg);
	fprintview(&ourview, stdout);
	putchar('\n');
	return(1);
}


static void
writepicture(void)				/* write out picture (alters buffer) */
{
	int	y;
	register int	x;
	double	d;

	fprtresolu(hresolu, vresolu, stdout);
	for (y = vresolu-1; y >= 0; y--)
		if (blurring) {
			for (x = 0; x < hresolu; x++) {	/* compute avg. */
				d = rexpadj/nvavg;
				scalecolor(bscan(y)[x], d);
			}
			if (fwritescan(bscan(y), hresolu, stdout) < 0)
				syserror(progname);
		} else if (averaging) {
			for (x = 0; x < hresolu; x++) {	/* average pixels */
				d = rexpadj/wscan(y)[x];
				scalecolor(sscan(y)[x], d);
			}
			if (fwritescan(sscan(y), hresolu, stdout) < 0)
				syserror(progname);
		} else {
			if (expadj)
				shiftcolrs(pscan(y), hresolu, expadj);
			if (fwritecolrs(pscan(y), hresolu, stdout) < 0)
				syserror(progname);
		}
}


static void
writedistance(			/* write out z file (alters buffer) */
	char	*fname
)
{
	int	donorm = normdist & !zisnorm ? 1 :
			(ourview.type == VT_PER) & !normdist & zisnorm ? -1 : 0;
	int	fd;
	int	y;

	if ((fd = open(fname, O_WRONLY|O_CREAT|O_TRUNC, 0666)) == -1)
		syserror(fname);
	for (y = vresolu-1; y >= 0; y--) {
		if (donorm) {
			double	vx, yzn2, d;
			register int	x;
			yzn2 = (y+.5)/vresolu + ourview.voff - .5;
			yzn2 = 1. + yzn2*yzn2*ourview.vn2;
			for (x = 0; x < hresolu; x++) {
				vx = (x+.5)/hresolu + ourview.hoff - .5;
				d = sqrt(vx*vx*ourview.hn2 + yzn2);
				if (donorm > 0)
					zscan(y)[x] *= d;
				else
					zscan(y)[x] /= d;
			}
		}
		if (write(fd, (char *)zscan(y), hresolu*sizeof(float))
				< hresolu*sizeof(float))
			syserror(fname);
	}
	close(fd);
}


static void
backfill(				/* fill pixel with background */
	int	x,
	int	y
)
{
	if (averaging) {
		copycolor(sscan(y)[x], backcolor);
		wscan(y)[x] = 1;
	} else
		copycolr(pscan(y)[x], backcolr);
	zscan(y)[x] = backz;
}


static void
calstart(                    /* start fill calculation */
	char	*prog,
	char	*args
)
{
	char	combuf[512];
	char	*argv[64];
	int	rval;
	register char	**wp, *cp;

	if (PDesc.running) {
		fprintf(stderr, "%s: too many calculations\n", progname);
		exit(1);
	}
	strcpy(combuf, prog);
	strcat(combuf, args);
	cp = combuf;
	wp = argv;
	for ( ; ; ) {
		while (isspace(*cp))	/* nullify spaces */
			*cp++ = '\0';
		if (!*cp)		/* all done? */
			break;
		*wp++ = cp;		/* add argument to list */
		while (*++cp && !isspace(*cp))
			;
	}
	*wp = NULL;
						/* start process */
	if ((rval = open_process(&PDesc, argv)) < 0)
		syserror(progname);
	if (rval == 0) {
		fprintf(stderr, "%s: command not found\n", argv[0]);
		exit(1);
	}
	packsiz = rval/(6*sizeof(float)) - 1;
	if (packsiz > PACKSIZ)
		packsiz = PACKSIZ;
	queuesiz = 0;
}


static void
caldone(void)                               /* done with calculation */
{
	if (!PDesc.running)
		return;
	clearqueue();
	close_process(&PDesc);
}


static void
rcalfill(				/* fill with ray-calculated pixel */
	int	x,
	int	y
)
{
	if (queuesiz >= packsiz)	/* flush queue if needed */
		clearqueue();
					/* add position to queue */
	queue[queuesiz][0] = x;
	queue[queuesiz][1] = y;
	queuesiz++;
}


static void
clearqueue(void)				/* process queue */
{
	FVECT	orig, dir;
	float	fbuf[6*(PACKSIZ+1)];
	register float	*fbp;
	register int	i;
	double	vx, vy;

	if (queuesiz == 0)
		return;
	fbp = fbuf;
	for (i = 0; i < queuesiz; i++) {
		viewray(orig, dir, &ourview,
				(queue[i][0]+.5)/hresolu,
				(queue[i][1]+.5)/vresolu);
		*fbp++ = orig[0]; *fbp++ = orig[1]; *fbp++ = orig[2];
		*fbp++ = dir[0]; *fbp++ = dir[1]; *fbp++ = dir[2];
	}
					/* mark end and get results */
	memset((char *)fbp, '\0', 6*sizeof(float));
	if (process(&PDesc, (char *)fbuf, (char *)fbuf,
			4*sizeof(float)*(queuesiz+1),
			6*sizeof(float)*(queuesiz+1)) !=
			4*sizeof(float)*(queuesiz+1)) {
		fprintf(stderr, "%s: error reading from rtrace process\n",
				progname);
		exit(1);
	}
	fbp = fbuf;
	for (i = 0; i < queuesiz; i++) {
		if (ourexp > 0 && ourexp != 1.0) {
			fbp[0] *= ourexp;
			fbp[1] *= ourexp;
			fbp[2] *= ourexp;
		}
		if (averaging) {
			setcolor(sscan(queue[i][1])[queue[i][0]],
					fbp[0], fbp[1], fbp[2]);
			wscan(queue[i][1])[queue[i][0]] = 1;
		} else
			setcolr(pscan(queue[i][1])[queue[i][0]],
					fbp[0], fbp[1], fbp[2]);
		if (zisnorm)
			zscan(queue[i][1])[queue[i][0]] = fbp[3];
		else {
			vx = (queue[i][0]+.5)/hresolu + ourview.hoff - .5;
			vy = (queue[i][1]+.5)/vresolu + ourview.voff - .5;
			zscan(queue[i][1])[queue[i][0]] = fbp[3] / sqrt(1. +
					vx*vx*ourview.hn2 + vy*vy*ourview.vn2);
		}
		fbp += 4;
	}
	queuesiz = 0;
}


static void
syserror(			/* report error and exit */
	char	*s
)
{
	perror(s);
	exit(1);
}
