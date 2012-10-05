#ifndef lint
static const char RCSid[] = "$Id: pmblur2.c,v 2.3 2012/10/05 15:25:12 greg Exp $";
#endif
/*
 *  pmblur2.c - program to computer better motion blur from ranimove frames.
 *  
 *  Created by Greg Ward on 10/3/12.
 */

#include "copyright.h"

#include "standard.h"
#include "platform.h"
#include "color.h"
#include "view.h"

#define ZGONE	(2.*FHUGE)	/* non-existent depth value */

const char	*hdrspec;	/* HDR image specification string */
const char	*zbfspec;	/* depth buffer specification string */
const char	*mvospec;	/* motion vector specification string */

RESOLU		imres;		/* image size and orientation */
double		pixaspect = 1.;	/* pixel aspect ratio */

COLOR		*imsum;		/* blurred image in progress */
VIEW		vwsum;		/* averaged view for blurred image */
char		imfmt[32] = PICFMT;	/* image format */

int		nsamps;		/* number of time samples */

COLOR		*imbuf;		/* accumulation buffer for this time sample */
float		*zbuf;		/* associated depth buffer */

COLR		*imprev;	/* previous loaded image frame */
float		*zprev;		/* loaded depth buffer */
double		exprev;		/* previous image exposure */
VIEW		vwprev = STDVIEW;	/* previous view */
int		fnprev;		/* frame number for previous view */

char		*progname;	/* global argv[0] */

typedef struct {
	VIEW	*vp;		/* view pointer */
	int	gotview;	/* got view parameters? */
	double	ev;		/* exposure value */
} IMGHEAD;	/* holder for individual header information */


/* Process line from picture header */
static int
headline(char *s, void *p)
{
	IMGHEAD	*ip = (IMGHEAD *)p;
	char	fmt[32];

	if (isview(s)) {
		ip->gotview += (sscanview(ip->vp, s) > 0);
		return(1);
	}
	if (isexpos(s)) {
		ip->ev *= exposval(s);
		return(1);
	}
	if (formatval(fmt, s)) {
		if (globmatch(imfmt, fmt)) {
			strcpy(imfmt, fmt);
			return(1);
		}
		error(WARNING, "wrong picture format");
		return(-1);
	}
	if (isaspect(s)) {
		pixaspect = aspectval(s);
		return(1);
	}
	return(0);
}


/* Load previous frame (allocates arrays on first call) */
static void
loadprev(int fno)
{
	char	*err;
	IMGHEAD	ih;
	RESOLU	rs;
	char	fname[256];
	FILE	*fp;
	int	y, hres;
	int	fd;

	if (fno == fnprev)		/* already loaded? */
		return;
	sprintf(fname, hdrspec, fno);	/* open picture file */
	if ((fp = fopen(fname, "r")) == NULL) {
		sprintf(errmsg, "cannot open picture file \"%s\"", fname);
		error(SYSTEM, errmsg);
	}
	ih.vp = &vwprev;
	ih.gotview = 0;
	ih.ev = 1.;
	if (getheader(fp, headline, &ih) < 0)
		goto readerr;
	if (!ih.gotview) {
		sprintf(errmsg, "missing view in picture \"%s\"", fname);
		error(USER, errmsg);
	}
	if ((err = setview(&vwprev)) != NULL)
		error(USER, err);
	exprev = ih.ev;
	if (!fgetsresolu(&rs, fp))
		goto readerr;
	if (rs.rt != PIXSTANDARD) {
		sprintf(errmsg, "unsupported orientation in picture \"%s\"", fname);
		error(USER, errmsg);
	}
	if (!imres.xr) {		/* allocate buffers */
		imres = rs;
		imsum = (COLOR *)ecalloc(imres.xr*imres.yr, sizeof(COLOR));
		imbuf = (COLOR *)emalloc(sizeof(COLOR)*imres.xr*imres.yr);
		zbuf = (float *)emalloc(sizeof(float)*imres.xr*imres.yr);
		imprev = (COLR *)emalloc(sizeof(COLR)*rs.xr*rs.yr);
		zprev = (float *)emalloc(sizeof(float)*rs.xr*rs.yr);
	} else if ((imres.xr != rs.xr) | (imres.yr != rs.yr)) {
		sprintf(errmsg, "resolution mismatch in picture \"%s\"", fname);
		error(USER, errmsg);
	}
	hres = scanlen(&rs);		/* load picture */
	for (y = 0; y < numscans(&rs); y++)
		if (freadcolrs(imprev + y*hres, hres, fp) < 0)
			goto readerr;
	fclose(fp);
	sprintf(fname, zbfspec, fno);	/* load depth buffer */
	if ((fd = open(fname, O_RDONLY)) < 0) {
		sprintf(errmsg, "cannot open depth buffer \"%s\"", fname);
		error(SYSTEM, errmsg);
	}
	if (read(fd, zprev, sizeof(float)*rs.xr*rs.yr) != sizeof(float)*rs.xr*rs.yr)
		goto readerr;
	close(fd);
	fnprev = fno;
	return;
readerr:
	sprintf(errmsg, "error loading file \"%s\"", fname);
	error(USER, errmsg);
}


/* Interpolate between two image pixels */
/* XXX skipping expensive ray vector calculations for now */
static void
interp_pixel(COLOR con, const VIEW *vwn, int xn, int yn, double zval,
			double pos, int xp, int yp, int interpOK)
{
	const int	hres = scanlen(&imres);
	RREAL		ipos;
	FVECT		wprev, wcoor, rdir;
	int		np, xd, yd, nd;
	COLOR		cpr;
	double		sf;
					/* check if off image */
	if ((xp < 0) | (xp >= hres))
		return;
	if ((yp < 0) | (yp >= numscans(&imres)))
		return;
	np = yp*hres + xp;		/* get indexed destination */
	xd = (1.-pos)*xp + pos*xn + .5;
	yd = (1.-pos)*yp + pos*yn + .5;
	nd = yd*hres + xd;
					/* check depth */
	if (interpOK)
		zval = (1.-pos)*zprev[np] + pos*zval;
	if (zval >= zbuf[nd])
		return;
	zbuf[nd] = zval;		/* assign new depth */
	copycolor(imbuf[nd], con);	/* assign new color */
	if (!interpOK)
		return;
	scalecolor(imbuf[nd], pos);
	sf = (1.-pos)/exprev;
	colr_color(cpr, imprev[np]);
	scalecolor(cpr, sf);
	addcolor(imbuf[nd], cpr);
}


/* Find neighbor with minimum depth from buffer */
static int
neigh_zmin(const float *zb, int n)
{
	const int	hres = scanlen(&imres);
	const int	vc = n / hres;
	const int	hc = n - vc*hres;
	float		zbest = ZGONE;
	int		nbest = 0;

	if (hc > 0 && zb[n-1] < zbest)
		zbest = zb[nbest = n-1];
	if (hc < hres-1 && zb[n+1] < zbest)
		zbest = zb[nbest = n+1];
	if (vc > 0 && zb[n-hres] < zbest)
		zbest = zb[nbest = n-hres];
	if (vc < numscans(&imres)-1 && zb[n+hres] < zbest)
		return(n+hres);
	return(nbest);
}


/* Fill in missing pixels from immediate neighbors */
static void
fill_missing(void)
{
	int	n, m;

	for (n = imres.xr*imres.yr; n--; )
		if (zbuf[n] >= .9*FHUGE &&
				zbuf[m = neigh_zmin(zbuf,n)] < .9*FHUGE)
			copycolor(imbuf[n], imbuf[m]);
}


/* Load a subframe */
static void
addframe(double fpos)
{
	COLOR		backg;
	char		*err;
	VIEW		fvw;
	RESOLU		rs;
	IMGHEAD		ihead;
	char		fname[256];
	FILE		*fpimg;
	int		fdzbf, fdmvo;
	COLR		*iscan;
	float		*zscan;
	uint16		*mscan;
	int		hres, n, h;
	
	loadprev((int)fpos);		/* make sure we have previous */
	fpos -= floor(fpos);
	if (fpos <= .02) {		/* special case */
		COLOR	col;
		for (n = imres.xr*imres.yr; n--; ) {
			colr_color(col, imprev[n]);
			addcolor(imsum[n], col);
		}
		return;
	}
					/* open input files */
	sprintf(fname, hdrspec, fnprev+1);
	if ((fpimg = fopen(fname, "r")) == NULL)
		goto openerr;
	fvw = vwprev;
	ihead.vp = &fvw;
	ihead.gotview = 0;
	ihead.ev = 1.;
	if (getheader(fpimg, headline, &ihead) < 0)
		goto readerr;
	if (!ihead.gotview) {
		sprintf(errmsg, "missing view in picture \"%s\"", fname);
		error(USER, errmsg);
	}
	if ((err = setview(&fvw)) != NULL)
		error(USER, err);
	if (!fgetsresolu(&rs, fpimg))
		goto readerr;
	if (rs.rt != PIXSTANDARD) {
		sprintf(errmsg, "unsupported orientation in picture \"%s\"", fname);
		error(USER, errmsg);
	}
	if ((rs.xr != imres.xr) | (rs.yr != imres.yr)) {
		sprintf(errmsg, "resolution mismatch for picture \"%s\"", fname);
		error(USER, errmsg);
	}
	vwsum.type = fvw.type;		/* average in view */
	for (n = 3; n--; ) {
		vwsum.vp[n] += (1.-fpos)*vwprev.vp[n] + fpos*fvw.vp[n];
		vwsum.vdir[n] += (1.-fpos)*vwprev.vdir[n] + fpos*fvw.vdir[n];
		vwsum.vup[n] += (1.-fpos)*vwprev.vdir[n] + fpos*fvw.vdir[n];
	}
	vwsum.vdist += (1.-fpos)*vwprev.vdist + fpos*fvw.vdist;
	vwsum.horiz += (1.-fpos)*vwprev.horiz + fpos*fvw.horiz;
	vwsum.vert += (1.-fpos)*vwprev.vert + fpos*fvw.vert;
	vwsum.voff += (1.-fpos)*vwprev.voff + fpos*fvw.voff;
	vwsum.hoff += (1.-fpos)*vwprev.hoff + fpos*fvw.hoff;
	vwsum.vfore += (1.-fpos)*vwprev.vfore + fpos*fvw.vfore;
	vwsum.vaft += (1.-fpos)*vwprev.vaft + fpos*fvw.vaft;
	sprintf(fname, zbfspec, fnprev+1);
	if ((fdzbf = open(fname, O_RDONLY)) < 0)
		goto openerr;
	sprintf(fname, mvospec, fnprev+1);
	if ((fdmvo = open(fname, O_RDONLY)) < 0)
		goto openerr;
	setcolor(backg, .0, .0, .0);	/* clear image & depth buffers */
	for (n = rs.xr*rs.yr; n--; ) {
		if (zprev[n] >= .9*FHUGE &&
				zprev[neigh_zmin(zprev,n)] >= .9*FHUGE)
			colr_color(backg, imprev[n]);
		copycolor(imbuf[n], backg);
		zbuf[n] = ZGONE;
	}
	hres = scanlen(&rs);		/* process a scanline at a time */
	iscan = (COLR *)emalloc(sizeof(COLR)*hres);
	zscan = (float *)emalloc(sizeof(float)*hres);
	mscan = (uint16 *)emalloc(sizeof(uint16)*3*hres);
	for (n = 0; n < numscans(&rs); n++) {
		const double	sf = 1./ihead.ev;
		COLOR		cval;
		if (freadcolrs(iscan, hres, fpimg) < 0)
			goto readerr;
		if (read(fdzbf, zscan, sizeof(float)*hres) !=
					sizeof(float)*hres)
			goto readerr;
		if (read(fdmvo, mscan, sizeof(uint16)*3*hres) !=
					sizeof(uint16)*3*hres)
			goto readerr;
		for (h = hres; h--; ) {
			if (!mscan[3*h])
				continue;	/* unknown motion */
			colr_color(cval, iscan[h]);
			scalecolor(cval, sf);
			interp_pixel(cval, &fvw, h, n, zscan[h], fpos,
					h + (int)mscan[3*h] - 32768,
					n - (int)mscan[3*h+1] + 32768,
					mscan[3*h+2]);
		}
	}
					/* fill in missing pixels */
	fill_missing();
					/* add in results */
	for (n = imres.xr*imres.yr; n--; )
		addcolor(imsum[n], imbuf[n]);
					/* clean up */
	free(mscan);
	free(zscan);
	free(iscan);
	close(fdmvo);
	close(fdzbf);
	fclose(fpimg);
	return;
openerr:
	sprintf(errmsg, "cannot open input \"%s\"", fname);
	error(SYSTEM, errmsg);
readerr:
	error(SYSTEM, "error reading input file");
}


/* Average and write our image out to the given file */
static void
write_average(FILE *fp)
{
	const int	hres = scanlen(&imres);
	const double	sf = 1./(double)nsamps;
	const float	csf = exprev*sf;
	int		n, h;
					/* normalize view & image */
	for (n = 3; n--; ) {
		vwsum.vp[n] *= sf;
		vwsum.vdir[n] *= sf;
		vwsum.vup[n] *= sf;
	}
	vwsum.vdist *= sf;
	vwsum.horiz *= sf;
	vwsum.vert *= sf;
	vwsum.voff *= sf;
	vwsum.hoff *= sf;
	vwsum.vfore *= sf;
	vwsum.vaft *= sf;
	for (n = imres.xr*imres.yr; n--; )
		scalecolor(imsum[n], csf);
					/* write it out */
	fputs(VIEWSTR, fp);
	fprintview(&vwsum, fp);
	fputc('\n', fp);
	if ((pixaspect < .98) | (pixaspect > 1.02))
		fputaspect(pixaspect, fp);
	fputexpos(exprev, fp);
	if (strcmp(imfmt, PICFMT))
		fputformat(imfmt, fp);
	fputc('\n', fp);		/* end of info. header */
	fputsresolu(&imres, fp);
	for (n = 0; n < numscans(&imres); n++)
		if (fwritescan(imsum + n*hres, hres, fp) < 0)
			error(SYSTEM, "write error on picture output");
}

/* Load frame sequence and compute motion blur output */
int
main(int argc, char *argv[])
{
	double	fstart, fend, fstep, fcur;

	progname = argv[0];
	SET_DEFAULT_BINARY();
	if (argc != 5)
		goto userr;
				/* get frame range & sampling */
	switch (sscanf(argv[1], "%lf,%lf/%d", &fstart, &fend, &nsamps)) {
	case 1:
		nsamps = 0;
		fend = fstart;
		break;
	case 2:
		nsamps = 0;
		/* fall through */
	case 3:
		if (fend < fstart)
			goto userr;
		break;
	default:
		goto userr;
	}
	if (fstart < 1)
		goto userr;
	if (nsamps <= 0)
		nsamps = (fend - fstart)*17.;
	if (nsamps) {
		if (nsamps > 100) nsamps = 100;
		fstep = (fend - fstart)/nsamps;
	} else {
		fstep = 1.;
		fstart -= .5;
		nsamps = 1;
	}
	hdrspec = argv[2];
	zbfspec = argv[3];
	mvospec = argv[4];
				/* load and filter each subframe */
	for (fcur = fstart + .5*fstep; fcur <= fend; fcur += fstep)
		addframe(fcur);
				/* write results to stdout */
	SET_FILE_BINARY(stdout);
	newheader("RADIANCE", stdout);
	printargs(argc, argv, stdout);
	fputnow(stdout);
	write_average(stdout);
	return(fflush(stdout) == EOF);
userr:
	fprintf(stderr, "Usage: %s f0,f1 HDRspec ZBUFspec MVOspec\n", progname);
	return(1);
}
