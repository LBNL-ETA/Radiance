#ifndef lint
static const char	RCSid[] = "$Id: pcond2.c,v 3.10 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 * Input and output conditioning routines for pcond.
 *  Added white-balance adjustment 10/01 (GW).
 */

#include "pcond.h"
#include "warp3d.h"


RGBPRIMP	inprims = stdprims;	/* input primaries */
COLORMAT	inrgb2xyz;		/* convert input RGB to XYZ */
RGBPRIMP	outprims = stdprims;	/* output primaries */

double	(*lumf)() = rgblum;		/* input luminance function */
double	inpexp = 1.0;			/* input exposure value */

char	*mbcalfile = NULL;		/* macbethcal mapping file */
char	*cwarpfile = NULL;		/* color space warping file */

static struct mbc {
	COLORMAT	cmat;
	float	xa[3][6], ya[3][6];
	COLOR	cmin, cmax;
}	mbcond;				/* macbethcal conditioning struct */

static WARP3D	*cwarp;			/* color warping structure */

static COLOR	*scanbuf;		/* scanline processing buffer */
static int	nread;			/* number of scanlines processed */


double
rgblum(clr, scotopic)		/* compute (scotopic) luminance of RGB color */
COLOR	clr;
int	scotopic;
{
	if (scotopic)		/* approximate */
		return( WHTSEFFICACY * (colval(clr,RED)*.062 +
			colval(clr,GRN)*.608 + colval(clr,BLU)*.330) );
	return( WHTEFFICACY * (colval(clr,RED)*inrgb2xyz[1][0] +
			colval(clr,GRN)*inrgb2xyz[1][1] +
			colval(clr,BLU)*inrgb2xyz[1][2]) );
}


double
cielum(xyz, scotopic)		/* compute (scotopic) luminance of CIE color */
COLOR	xyz;
int	scotopic;
{
	if (scotopic)		/* approximate */
		return(colval(xyz,CIEY) *
			(1.33*(1. + (colval(xyz,CIEY)+colval(xyz,CIEZ))/
			colval(xyz,CIEX)) - 1.68));
	return(colval(xyz,CIEY));
}


COLOR *
nextscan()				/* read and condition next scanline */
{
	if (nread >= numscans(&inpres)) {
		if (cwarpfile != NULL)
			free3dw(cwarp);
		free((void *)scanbuf);
		return(scanbuf = NULL);
	}
	if (what2do&DO_ACUITY)
		acuscan(scanbuf, nread);
	else if (freadscan(scanbuf, scanlen(&inpres), infp) < 0) {
		fprintf(stderr, "%s: %s: scanline read error\n",
				progname, infn);
		exit(1);
	}
	if (what2do&DO_VEIL)			/* add veiling */
		addveil(scanbuf, nread);
	if (what2do&DO_COLOR)			/* scotopic color loss */
		scotscan(scanbuf, scanlen(&inpres));
	if (what2do&DO_LINEAR)			/* map luminances */
		sfscan(scanbuf, scanlen(&inpres), scalef);
	else
		mapscan(scanbuf, scanlen(&inpres));
	if (mbcalfile != NULL)			/* device color correction */
		mbscan(scanbuf, scanlen(&inpres), &mbcond);
	else if (cwarpfile != NULL)		/* device color space warp */
		cwscan(scanbuf, scanlen(&inpres), cwarp);
	else if (lumf == cielum | inprims != outprims)
		matscan(scanbuf, scanlen(&inpres), mbcond.cmat);
	nread++;
	return(scanbuf);
}


COLOR *
firstscan()				/* return first processed scanline */
{
	if (mbcalfile != NULL)		/* load macbethcal file */
		getmbcalfile(mbcalfile, &mbcond);
	else if (cwarpfile != NULL) {
		if ((cwarp = load3dw(cwarpfile, NULL)) == NULL)
			syserror(cwarpfile);
	} else
		if (lumf == rgblum)
			comprgb2rgbWBmat(mbcond.cmat, inprims, outprims);
		else
			compxyz2rgbWBmat(mbcond.cmat, outprims);
	if (what2do&DO_ACUITY)
		initacuity();
	scanbuf = (COLOR *)malloc(scanlen(&inpres)*sizeof(COLOR));
	if (scanbuf == NULL)
		syserror("malloc");
	nread = 0;
	return(nextscan());
}


sfscan(sl, len, sf)			/* apply scalefactor to scanline */
register COLOR	*sl;
int	len;
double	sf;
{
	while (len--) {
		scalecolor(sl[0], sf);
		sl++;
	}
}


matscan(sl, len, mat)			/* apply color matrix to scaline */
register COLOR	*sl;
int	len;
COLORMAT	mat;
{
	while (len--) {
		colortrans(sl[0], mat, sl[0]);
		clipgamut(sl[0], bright(sl[0]), CGAMUT, cblack, cwhite);
		sl++;
	}
}


mbscan(sl, len, mb)			/* apply macbethcal adj. to scaline */
COLOR	*sl;
int	len;
register struct mbc	*mb;
{
	double	d;
	register int	i, j;

	while (len--) {
		colortrans(sl[0], mb->cmat, sl[0]);
		clipgamut(sl[0], bright(sl[0]), CGAMUT, mb->cmin, mb->cmax);
		for (i = 0; i < 3; i++) {
			d = colval(sl[0],i);
			for (j = 0; j < 4 && mb->xa[i][j+1] <= d; j++)
				;
			colval(sl[0],i) = ( (mb->xa[i][j+1] - d)*mb->ya[i][j] +
					(d - mb->xa[i][j])*mb->ya[i][j+1] ) /
					(mb->xa[i][j+1] - mb->xa[i][j]);
		}
		sl++;
	}
}


cwscan(sl, len, wp)			/* apply color space warp to scaline */
COLOR	*sl;
int	len;
WARP3D	*wp;
{
	int	rval;

	while (len--) {
		rval = warp3d(sl[0], sl[0], wp);
		if (rval & W3ERROR)
			syserror("warp3d");
		if (rval & W3BADMAP) {
			fprintf(stderr, "%s: %s: bad color space map\n",
					progname, cwarpfile);
			exit(1);
		}
		clipgamut(sl[0], bright(sl[0]), CGAMUT, cblack, cwhite);
		sl++;
	}
}


getmbcalfile(fn, mb)			/* load macbethcal file */
char	*fn;
register struct mbc	*mb;
{
	extern char	*fgets();
	char	buf[128];
	FILE	*fp;
	int	inpflags = 0;
	register int	i;

	if ((fp = fopen(fn, "r")) == NULL)
		syserror(fn);
	while (fgets(buf, sizeof(buf), fp) != NULL) {
		if (!(inpflags & 01) &&
				sscanf(buf,
				"rxa(i) : select(i,%f,%f,%f,%f,%f,%f)",
				&mb->xa[0][0], &mb->xa[0][1],
				&mb->xa[0][2], &mb->xa[0][3],
				&mb->xa[0][4], &mb->xa[0][5]
				) == 6)
			inpflags |= 01;
		else if (!(inpflags & 02) &&
				sscanf(buf,
				"rya(i) : select(i,%f,%f,%f,%f,%f,%f)",
				&mb->ya[0][0], &mb->ya[0][1],
				&mb->ya[0][2], &mb->ya[0][3],
				&mb->ya[0][4], &mb->ya[0][5]
				) == 6)
			inpflags |= 02;
		else if (!(inpflags & 04) &&
				sscanf(buf,
				"gxa(i) : select(i,%f,%f,%f,%f,%f,%f)",
				&mb->xa[1][0], &mb->xa[1][1],
				&mb->xa[1][2], &mb->xa[1][3],
				&mb->xa[1][4], &mb->xa[1][5]
				) == 6)
			inpflags |= 04;
		else if (!(inpflags & 010) &&
				sscanf(buf,
				"gya(i) : select(i,%f,%f,%f,%f,%f,%f)",
				&mb->ya[1][0], &mb->ya[1][1],
				&mb->ya[1][2], &mb->ya[1][3],
				&mb->ya[1][4], &mb->ya[1][5]
				) == 6)
			inpflags |= 010;
		else if (!(inpflags & 020) &&
				sscanf(buf,
				"bxa(i) : select(i,%f,%f,%f,%f,%f,%f)",
				&mb->xa[2][0], &mb->xa[2][1],
				&mb->xa[2][2], &mb->xa[2][3],
				&mb->xa[2][4], &mb->xa[2][5]
				) == 6)
			inpflags |= 020;
		else if (!(inpflags & 040) &&
				sscanf(buf,
				"bya(i) : select(i,%f,%f,%f,%f,%f,%f)",
				&mb->ya[2][0], &mb->ya[2][1],
				&mb->ya[2][2], &mb->ya[2][3],
				&mb->ya[2][4], &mb->ya[2][5]
				) == 6)
			inpflags |= 040;
		else if (!(inpflags & 0100) &&
				sscanf(buf,
				"ro = %f*rn + %f*gn + %f*bn",
				&mb->cmat[0][0], &mb->cmat[0][1],
				&mb->cmat[0][2]) == 3)
			inpflags |= 0100;
		else if (!(inpflags & 0200) &&
				sscanf(buf,
				"go = %f*rn + %f*gn + %f*bn",
				&mb->cmat[1][0], &mb->cmat[1][1],
				&mb->cmat[1][2]) == 3)
			inpflags |= 0200;
		else if (!(inpflags & 0400) &&
				sscanf(buf,
				"bo = %f*rn + %f*gn + %f*bn",
				&mb->cmat[2][0], &mb->cmat[2][1],
				&mb->cmat[2][2]) == 3)
			inpflags |= 0400;
	}
	if (inpflags != 0777) {
		fprintf(stderr,
		"%s: cannot grok macbethcal file \"%s\" (inpflags==0%o)\n",
			progname, fn, inpflags);
		exit(1);
	}
	fclose(fp);
					/* compute gamut */
	for (i = 0; i < 3; i++) {
		colval(mb->cmin,i) = mb->xa[i][0] -
				mb->ya[i][0] *
				(mb->xa[i][1]-mb->xa[i][0]) /
				(mb->ya[i][1]-mb->ya[i][0]);
		colval(mb->cmax,i) = mb->xa[i][4] +
				(1.-mb->ya[i][4]) *
				(mb->xa[i][5] - mb->xa[i][4]) /
				(mb->ya[i][5] - mb->ya[i][4]);
	}
}
