#ifndef lint
static const char	RCSid[] = "$Id: pcond4.c,v 3.16 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 * Routines for veiling glare and loss of acuity.
 */

#include "pcond.h"

/************** VEILING STUFF *****************/

#define	VADAPT		0.08		/* fraction of adaptation from veil */

static COLOR	*veilimg = NULL;	/* veiling image */

#define veilscan(y)	(veilimg+(y)*fvxr)

static float	(*raydir)[3] = NULL;	/* ray direction for each pixel */

#define rdirscan(y)	(raydir+(y)*fvxr)


compraydir()				/* compute ray directions */
{
	FVECT	rorg, rdir;
	double	h, v;
	register int	x, y;

	if (raydir != NULL)		/* already done? */
		return;
	raydir = (float (*)[3])malloc(fvxr*fvyr*3*sizeof(float));
	if (raydir == NULL)
		syserror("malloc");

	for (y = 0; y < fvyr; y++) {
		switch (inpres.rt) {
		case YMAJOR: case YMAJOR|XDECR:
			v = (y+.5)/fvyr; break;
		case YMAJOR|YDECR: case YMAJOR|YDECR|XDECR:
			v = 1. - (y+.5)/fvyr; break;
		case 0: case YDECR:
			h = (y+.5)/fvyr; break;
		case XDECR: case XDECR|YDECR:
			h = 1. - (y+.5)/fvyr; break;
		}
		for (x = 0; x < fvxr; x++) {
			switch (inpres.rt) {
			case YMAJOR: case YMAJOR|YDECR:
				h = (x+.5)/fvxr; break;
			case YMAJOR|XDECR: case YMAJOR|XDECR|YDECR:
				h = 1. - (x+.5)/fvxr; break;
			case 0: case XDECR:
				v = (x+.5)/fvxr; break;
			case YDECR: case YDECR|XDECR:
				v = 1. - (x+.5)/fvxr; break;
			}
			if (viewray(rorg, rdir, &ourview, h, v)
					>= -FTINY) {
				rdirscan(y)[x][0] = rdir[0];
				rdirscan(y)[x][1] = rdir[1];
				rdirscan(y)[x][2] = rdir[2];
			} else {
				rdirscan(y)[x][0] =
				rdirscan(y)[x][1] =
				rdirscan(y)[x][2] = 0.0;
			}
		}
	}
}


compveil()				/* compute veiling image */
{
	double	t2, t2sum;
	COLOR	ctmp, vsum;
	int	px, py;
	register int	x, y;

	if (veilimg != NULL)		/* already done? */
		return;
					/* compute ray directions */
	compraydir();
					/* compute veil image */
	veilimg = (COLOR *)malloc(fvxr*fvyr*sizeof(COLOR));
	if (veilimg == NULL)
		syserror("malloc");
	for (py = 0; py < fvyr; py++)
		for (px = 0; px < fvxr; px++) {
			t2sum = 0.;
			setcolor(vsum, 0., 0., 0.);
			for (y = 0; y < fvyr; y++)
				for (x = 0; x < fvxr; x++) {
					if (x == px && y == py) continue;
					t2 = DOT(rdirscan(py)[px],
							rdirscan(y)[x]);
					if (t2 <= FTINY) continue;
					/*	use approximation instead
					t3 = acos(t2);
					t2 = t2/(t3*t3);
					*/
					t2 *= .5 / (1. - t2);
					copycolor(ctmp, fovscan(y)[x]);
					scalecolor(ctmp, t2);
					addcolor(vsum, ctmp);
					t2sum += t2;
				}
			/* VADAPT of original is subtracted in addveil() */
			if (t2sum > FTINY)
				scalecolor(vsum, VADAPT/t2sum);
			copycolor(veilscan(py)[px], vsum);
		}
					/* modify FOV sample image */
	for (y = 0; y < fvyr; y++)
		for (x = 0; x < fvxr; x++) {
			scalecolor(fovscan(y)[x], 1.-VADAPT);
			addcolor(fovscan(y)[x], veilscan(y)[x]);
		}
	comphist();			/* recompute histogram */
}


#if ADJ_VEIL
/*
 * The following veil adjustment was added to compensate for
 * the fact that contrast reduction gets confused with veil
 * in the human visual system.  Therefore, we reduce the
 * veil in portions of the image where our mapping has
 * already reduced contrast below the target value.
 * This gets called after the intial veil has been computed
 * and added to the foveal image, and the mapping has been
 * determined.
 */
adjveil()				/* adjust veil image */
{
	float	*crfptr = crfimg;
	COLOR	*fovptr = fovimg;
	COLOR	*veilptr = veilimg;
	double	s2nits = 1./inpexp;
	double	vl, vl2, fovl, vlsum;
	double	deltavc[3];
	int	i, j;

	if (lumf == rgblum)
		s2nits *= WHTEFFICACY;

	for (i = fvxr*fvyr; i--; crfptr++, fovptr++, veilptr++) {
		if (crfptr[0] >= 0.95)
			continue;
		vl = plum(veilptr[0]);
		fovl = (plum(fovptr[0]) - vl) * (1./(1.-VADAPT));
		if (vl <= 0.05*fovl)
			continue;
		vlsum = vl;
		for (j = 2; j < 11; j++) {
			vlsum += crfptr[0]*vl - (1.0 - crfptr[0])*fovl;
			vl2 = vlsum / (double)j;
			if (vl2 < 0.0)
				vl2 = 0.0;
			crfptr[0] = crfactor(fovl + vl2);
		}
		/* desaturation code causes color fringes at this level */
		for (j = 3; j--; ) {
			double	vc = colval(veilptr[0],j);
			double	fovc = (colval(fovptr[0],j) - vc) *
						(1./(1.-VADAPT));
			deltavc[j] = (1.-crfptr[0])*(fovl/s2nits - fovc);
			if (vc + deltavc[j] < 0.0)
				break;
		}
		if (j < 0)
			addcolor(veilptr[0], deltavc);
		else
			scalecolor(veilptr[0], vl2/vl);
	}
	smoothveil();			/* smooth our result */
}


smoothveil()				/* smooth veil image */
{
	COLOR 	*nveilimg;
	COLOR	*ovptr, *nvptr;
	int	x, y, i;

	nveilimg = (COLOR *)malloc(fvxr*fvyr*sizeof(COLOR));
	if (nveilimg == NULL)
		return;
	for (y = 1; y < fvyr-1; y++) {
		ovptr = veilimg + y*fvxr + 1;
		nvptr = nveilimg + y*fvxr + 1;
		for (x = 1; x < fvxr-1; x++, ovptr++, nvptr++)
			for (i = 3; i--; )
				nvptr[0][i] = 0.5 * ovptr[0][i]
					+ (1./12.) *
				(ovptr[-1][i] + ovptr[-fvxr][i] +
					ovptr[1][i] + ovptr[fvxr][i])
					+ (1./24.) *
				(ovptr[-fvxr-1][i] + ovptr[-fvxr+1][i] +
					ovptr[fvxr-1][i] + ovptr[fvxr+1][i]);
	}
	ovptr = veilimg + 1;
	nvptr = nveilimg + 1;
	for (x = 1; x < fvxr-1; x++, ovptr++, nvptr++)
		for (i = 3; i--; )
			nvptr[0][i] = 0.5 * ovptr[0][i]
				+ (1./9.) *
			(ovptr[-1][i] + ovptr[1][i] + ovptr[fvxr][i])
				+ (1./12.) *
			(ovptr[fvxr-1][i] + ovptr[fvxr+1][i]);
	ovptr = veilimg + (fvyr-1)*fvxr + 1;
	nvptr = nveilimg + (fvyr-1)*fvxr + 1;
	for (x = 1; x < fvxr-1; x++, ovptr++, nvptr++)
		for (i = 3; i--; )
			nvptr[0][i] = 0.5 * ovptr[0][i]
				+ (1./9.) *
			(ovptr[-1][i] + ovptr[1][i] + ovptr[-fvxr][i])
				+ (1./12.) *
			(ovptr[-fvxr-1][i] + ovptr[-fvxr+1][i]);
	ovptr = veilimg + fvxr;
	nvptr = nveilimg + fvxr;
	for (y = 1; y < fvyr-1; y++, ovptr += fvxr, nvptr += fvxr)
		for (i = 3; i--; )
			nvptr[0][i] = 0.5 * ovptr[0][i]
				+ (1./9.) *
			(ovptr[-fvxr][i] + ovptr[1][i] + ovptr[fvxr][i])
				+ (1./12.) *
			(ovptr[-fvxr+1][i] + ovptr[fvxr+1][i]);
	ovptr = veilimg + fvxr - 1;
	nvptr = nveilimg + fvxr - 1;
	for (y = 1; y < fvyr-1; y++, ovptr += fvxr, nvptr += fvxr)
		for (i = 3; i--; )
			nvptr[0][i] = 0.5 * ovptr[0][i]
				+ (1./9.) *
			(ovptr[-fvxr][i] + ovptr[-1][i] + ovptr[fvxr][i])
				+ (1./12.) *
			(ovptr[-fvxr-1][i] + ovptr[fvxr-1][i]);
	for (i = 3; i--; ) {
		nveilimg[0][i] = veilimg[0][i];
		nveilimg[fvxr-1][i] = veilimg[fvxr-1][i];
		nveilimg[(fvyr-1)*fvxr][i] = veilimg[(fvyr-1)*fvxr][i];
		nveilimg[fvyr*fvxr-1][i] = veilimg[fvyr*fvxr-1][i];
	}
	free((void *)veilimg);
	veilimg = nveilimg;
}
#endif

addveil(sl, y)				/* add veil to scanline */
COLOR	*sl;
int	y;
{
	int	vx, vy;
	double	dx, dy;
	double	lv, uv;
	register int	x, i;

	vy = dy = (y+.5)/numscans(&inpres)*fvyr - .5;
	while (vy >= fvyr-1) vy--;
	dy -= (double)vy;
	for (x = 0; x < scanlen(&inpres); x++) {
		vx = dx = (x+.5)/scanlen(&inpres)*fvxr - .5;
		while (vx >= fvxr-1) vx--;
		dx -= (double)vx;
		for (i = 0; i < 3; i++) {
			lv = (1.-dy)*colval(veilscan(vy)[vx],i) +
					dy*colval(veilscan(vy+1)[vx],i);
			uv = (1.-dy)*colval(veilscan(vy)[vx+1],i) +
					dy*colval(veilscan(vy+1)[vx+1],i);
			colval(sl[x],i) = (1.-VADAPT)*colval(sl[x],i) +
					(1.-dx)*lv + dx*uv;
		}
	}
}


/****************** ACUITY STUFF *******************/

typedef struct {
	short	sampe;		/* sample area size (exponent of 2) */
	short	nscans;		/* number of scanlines in this bar */
	int	len;		/* individual scanline length */
	int	nread;		/* number of scanlines loaded */
	COLOR	*sdata;		/* scanbar data */
} SCANBAR;

#define bscan(sb,y)	((COLOR *)(sb)->sdata+((y)%(sb)->nscans)*(sb)->len)

SCANBAR	*rootbar;		/* root scan bar (lowest resolution) */

float	*inpacuD;		/* input acuity data (cycles/degree) */

#define tsampr(x,y)	inpacuD[(y)*fvxr+(x)]


double
hacuity(La)			/* return visual acuity in cycles/degree */
double	La;
{
					/* functional fit */
	return(17.25*atan(1.4*log10(La) + 0.35) + 25.72);
}


COLOR *
getascan(sb, y)			/* find/read scanline y for scanbar sb */
register SCANBAR	*sb;
int	y;
{
	register COLOR	*sl0, *sl1, *mysl;
	register int	i;

	if (y < sb->nread - sb->nscans)			/* too far back? */
		return(NULL);
	for ( ; y >= sb->nread; sb->nread++) {		/* read as necessary */
		mysl = bscan(sb, sb->nread);
		if (sb->sampe == 0) {
			if (freadscan(mysl, sb->len, infp) < 0) {
				fprintf(stderr, "%s: %s: scanline read error\n",
						progname, infn);
				exit(1);
			}
		} else {
			sl0 = getascan(sb+1, 2*y);
			if (sl0 == NULL)
				return(NULL);
			sl1 = getascan(sb+1, 2*y+1);
			for (i = 0; i < sb->len; i++) {
				copycolor(mysl[i], sl0[2*i]);
				addcolor(mysl[i], sl0[2*i+1]);
				addcolor(mysl[i], sl1[2*i]);
				addcolor(mysl[i], sl1[2*i+1]);
				scalecolor(mysl[i], 0.25);
			}
		}
	}
	return(bscan(sb, y));
}


acuscan(scln, y)		/* get acuity-sampled scanline */
COLOR	*scln;
int	y;
{
	double	sr;
	double	dx, dy;
	int	ix, iy;
	register int	x;
					/* compute foveal y position */
	iy = dy = (y+.5)/numscans(&inpres)*fvyr - .5;
	while (iy >= fvyr-1) iy--;
	dy -= (double)iy;
	for (x = 0; x < scanlen(&inpres); x++) {
					/* compute foveal x position */
		ix = dx = (x+.5)/scanlen(&inpres)*fvxr - .5;
		while (ix >= fvxr-1) ix--;
		dx -= (double)ix;
					/* interpolate sample rate */
		sr = (1.-dy)*((1.-dx)*tsampr(ix,iy) + dx*tsampr(ix+1,iy)) +
			dy*((1.-dx)*tsampr(ix,iy+1) + dx*tsampr(ix+1,iy+1));

		acusample(scln[x], x, y, sr);	/* compute sample */
	}
}


acusample(col, x, y, sr)	/* interpolate sample at (x,y) using rate sr */
COLOR	col;
int	x, y;
double	sr;
{
	COLOR	c1;
	double	d;
	register SCANBAR	*sb0;

	for (sb0 = rootbar; sb0->sampe != 0 && 1<<sb0[1].sampe > sr; sb0++)
		;
	ascanval(col, x, y, sb0);
	if (sb0->sampe == 0)		/* don't extrapolate highest */
		return;
	ascanval(c1, x, y, sb0+1);
	d = ((1<<sb0->sampe) - sr)/(1<<sb0[1].sampe);
	scalecolor(col, 1.-d);
	scalecolor(c1, d);
	addcolor(col, c1);
}


ascanval(col, x, y, sb)		/* interpolate scanbar at orig. coords (x,y) */
COLOR	col;
int	x, y;
SCANBAR	*sb;
{
	COLOR	*sl0, *sl1, c1, c1y;
	double	dx, dy;
	int	ix, iy;

	if (sb->sampe == 0) {		/* no need to interpolate */
		sl0 = getascan(sb, y);
		copycolor(col, sl0[x]);
		return;
	}
					/* compute coordinates for sb */
	ix = dx = (x+.5)/(1<<sb->sampe) - .5;
	while (ix >= sb->len-1) ix--;
	dx -= (double)ix;
	iy = dy = (y+.5)/(1<<sb->sampe) - .5;
	while (iy >= (numscans(&inpres)>>sb->sampe)-1) iy--;
	dy -= (double)iy;
					/* get scanlines */
	sl0 = getascan(sb, iy);
#ifdef DEBUG
	if (sl0 == NULL)
		error(INTERNAL, "cannot backspace in ascanval");
#endif
	sl1 = getascan(sb, iy+1);
					/* 2D linear interpolation */
	copycolor(col, sl0[ix]);
	scalecolor(col, 1.-dx);
	copycolor(c1, sl0[ix+1]);
	scalecolor(c1, dx);
	addcolor(col, c1);
	copycolor(c1y, sl1[ix]);
	scalecolor(c1y, 1.-dx);
	copycolor(c1, sl1[ix+1]);
	scalecolor(c1, dx);
	addcolor(c1y, c1);
	scalecolor(col, 1.-dy);
	scalecolor(c1y, dy);
	addcolor(col, c1y);
	for (ix = 0; ix < 3; ix++)	/* make sure no negative */
		if (colval(col,ix) < 0.)
			colval(col,ix) = 0.;
}


SCANBAR	*
sballoc(se, ns, sl)		/* allocate scanbar */
int	se;		/* sampling rate exponent */
int	ns;		/* number of scanlines */
int	sl;		/* original scanline length */
{
	SCANBAR	*sbarr;
	register SCANBAR	*sb;

	sbarr = sb = (SCANBAR *)malloc((se+1)*sizeof(SCANBAR));
	if (sb == NULL)
		syserror("malloc");
	do {
		sb->len = sl>>se;
		if (sb->len <= 0)
			continue;
		sb->sampe = se;
		sb->nscans = ns;
		sb->sdata = (COLOR *)malloc(sb->len*ns*sizeof(COLOR));
		if (sb->sdata == NULL)
			syserror("malloc");
		sb->nread = 0;
		ns <<= 1;
		sb++;
	} while (--se >= 0);
	return(sbarr);
}


initacuity()			/* initialize variable acuity sampling */
{
	FVECT	diffx, diffy, cp;
	double	omega, maxsr;
	register int	x, y, i;

	compraydir();			/* compute ray directions */

	inpacuD = (float *)malloc(fvxr*fvyr*sizeof(float));
	if (inpacuD == NULL)
		syserror("malloc");
	maxsr = 1.;			/* compute internal sample rates */
	for (y = 1; y < fvyr-1; y++)
		for (x = 1; x < fvxr-1; x++) {
			for (i = 0; i < 3; i++) {
				diffx[i] = 0.5*fvxr/scanlen(&inpres) *
						(rdirscan(y)[x+1][i] -
						rdirscan(y)[x-1][i]);
				diffy[i] = 0.5*fvyr/numscans(&inpres) *
						(rdirscan(y+1)[x][i] -
						rdirscan(y-1)[x][i]);
			}
			fcross(cp, diffx, diffy);
			omega = 0.5 * sqrt(DOT(cp,cp));
			if (omega <= FTINY*FTINY)
				tsampr(x,y) = 1.;
			else if ((tsampr(x,y) = PI/180. / sqrt(omega) /
					hacuity(plum(fovscan(y)[x]))) > maxsr)
				maxsr = tsampr(x,y);
		}
					/* copy perimeter (easier) */
	for (x = 1; x < fvxr-1; x++) {
		tsampr(x,0) = tsampr(x,1);
		tsampr(x,fvyr-1) = tsampr(x,fvyr-2);
	}
	for (y = 0; y < fvyr; y++) {
		tsampr(0,y) = tsampr(1,y);
		tsampr(fvxr-1,y) = tsampr(fvxr-2,y);
	}
					/* initialize with next power of two */
	rootbar = sballoc((int)(log(maxsr)/log(2.))+1, 2, scanlen(&inpres));
}
