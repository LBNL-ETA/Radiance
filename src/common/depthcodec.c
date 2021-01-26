#ifndef lint
static const char RCSid[] = "$Id: depthcodec.c,v 2.11 2021/01/26 18:47:25 greg Exp $";
#endif
/*
 * Routines to encode/decoded 16-bit depths
 */

#include "copyright.h"

#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "rtio.h"
#include "fvect.h"
#include "depthcodec.h"


#ifndef depth2code
/* Encode depth as 16-bit signed integer */
int
depth2code(double d, double dref)
{
	if (d <= .0)
		return -32768;

	if (d > dref)
		return (int)(32768.001 - 32768.*dref/d) - 1;

	return (int)(32767.*d/dref - 32768.999);
}
#endif


#ifndef code2depth
/* Decode depth from 16-bit signed integer */
double
code2depth(int c, double dref)
{
	if (c <= -32768)
		return .0;

	if (c >= 32767)
		return FHUGE;

	if (c >= -1)
		return dref*32768./(32766.5 - c);

	return dref*(32768.5 + c)*(1./32767.);
}
#endif


/* Set codec defaults */
void
set_dc_defaults(DEPTHCODEC *dcp)
{
	memset(dcp, 0, sizeof(DEPTHCODEC));
	dcp->finp = stdin;
	dcp->inpname = "<stdin>";
	dcp->format = 'a';
	dcp->refdepth = 1.;
	dcp->depth_unit[0] = '1';
	dcp->vw = stdview;
	dcp->res.rt = PIXSTANDARD;
	if (!progname) progname = "depth_codec";
}


/* process header line */
static int
headline(char *s, void *p)
{
	DEPTHCODEC	*dcp = (DEPTHCODEC *)p;
	int		rv;

	if (formatval(dcp->inpfmt, s))	/* don't pass format */
		return 0;

	if ((rv = isbigendian(s)) >= 0) {
		dcp->swapped = (nativebigendian() != rv);
		return 0;
	}
					/* check for reference depth */
	if (!strncmp(s, DEPTHSTR, LDEPTHSTR)) {
		char	*cp;
		strlcpy(dcp->depth_unit, s+LDEPTHSTR, sizeof(dcp->depth_unit));
		cp = dcp->depth_unit;
		while (*cp) cp++;
		while (cp > dcp->depth_unit && isspace(cp[-1])) cp--;
		*cp = '\0';
		dcp->refdepth = atof(dcp->depth_unit);
		if (dcp->refdepth <= .0) {
			if (dcp->hdrflags & HF_STDERR) {
				fputs(dcp->inpname, stderr);
				fputs(": bad reference depth in input header\n",
						stderr);
			}
			return -1;
		}
		if (dcp->hdrflags & HF_ENCODE)
			return 0;	/* will add this later */
	} else if (!strncmp(s, "SAMP360=", 8))
		dcp->gotview--;
	else if (isview(s))		/* get view params */
		dcp->gotview += (sscanview(&dcp->vw, s) > 0);
	if (dcp->hdrflags & HF_HEADOUT)
		fputs(s, stdout);	/* copy to standard output */
	return 1;
}


/* Load/copy header */
int
process_dc_header(DEPTHCODEC *dcp, int ac, char *av[])
{
	if (dcp->hdrflags & HF_HEADIN &&
			getheader(dcp->finp, headline, dcp) < 0) {
		if (dcp->hdrflags & HF_STDERR) {
			fputs(dcp->inpname, stderr);
			fputs(": bad header\n", stderr);
		}
		return 0;
	}
	dcp->gotview *= (dcp->gotview > 0);
	if (dcp->hdrflags & HF_HEADOUT) {	/* finish header */
		if (!(dcp->hdrflags & HF_HEADIN))
			newheader("RADIANCE", stdout);
		if (ac > 0)
			printargs(ac, av, stdout);
		if (dcp->hdrflags & HF_ENCODE) {
			fputs(DEPTHSTR, stdout);
			fputs(dcp->depth_unit, stdout);
			fputc('\n', stdout);
			fputformat(DEPTH16FMT, stdout);
		} else
			switch (dcp->format) {
			case 'a':
				fputformat("ascii", stdout);
				break;
			case 'f':
				fputendian(stdout);
				fputformat("float", stdout);
				break;
			case 'd':
				fputendian(stdout);
				fputformat("double", stdout);
				break;
			}
		fputc('\n', stdout);
	}
					/* get/put resolution string */
	if (dcp->hdrflags & HF_RESIN && !fgetsresolu(&dcp->res, dcp->finp)) {
		if (dcp->hdrflags & HF_STDERR) {
			fputs(dcp->inpname, stderr);
			fputs(": bad resolution string\n", stderr);
		}
		return 0;
	}
	if (dcp->hdrflags & HF_RESOUT)
		fputsresolu(&dcp->res, stdout);

	dcp->dstart = dcp->curpos = ftell(dcp->finp);
	return 1;
}


/* Check that we have what we need to decode depths */
int
check_decode_depths(DEPTHCODEC *dcp)
{
	if (dcp->hdrflags & HF_ENCODE) {
		if (dcp->hdrflags & HF_STDERR) {
			fputs(progname, stderr);
			fputs(": wrong header mode for decode\n", stderr);
		}
		return 0;
	}
	if (dcp->inpfmt[0] && strcmp(dcp->inpfmt, DEPTH16FMT)) {
		if (dcp->hdrflags & HF_STDERR) {
			fputs(dcp->inpname, stderr);
			fputs(": unexpected input format: ", stderr);
			fputs(dcp->inpfmt, stderr);
			fputc('\n', stderr);
		}
		return 0;
	}
	return 1;
}


/* Check that we have what we need to decode world positions */
int
check_decode_worldpos(DEPTHCODEC *dcp)
{
	char	*err;

	if (!check_decode_depths(dcp))
		return 0;
	if ((dcp->res.xr <= 0) | (dcp->res.yr <= 0)) {
		if (dcp->hdrflags & HF_STDERR) {
			fputs(progname, stderr);
			fputs(": missing map resolution\n", stderr);
		}
		return 0;
	}
	if (!dcp->gotview) {
		if (dcp->hdrflags & HF_STDERR) {
			fputs(dcp->inpname, stderr);
			fputs(": missing view\n", stderr);
		}
		return 0;
	}
	if ((err = setview(&dcp->vw)) != NULL) {
		if (dcp->hdrflags & HF_STDERR) {
			fputs(dcp->inpname, stderr);
			fputs(": input view error: ", stderr);
			fputs(err, stderr);
			fputc('\n', stderr);
		}
		return 0;
	}
	return 1;	
}


/* Decode next depth pixel */
double
decode_depth_next(DEPTHCODEC *dcp)
{
	int	c;

	if (dcp->use_last) {
		dcp->use_last = 0;
		return code2depth(dcp->last_dc, dcp->refdepth);
	}
	c = getint(2, dcp->finp);

	if (c == EOF && feof(dcp->finp))
		return -1.;

	dcp->last_dc = c;
	dcp->curpos += 2;

	return code2depth(c, dcp->refdepth);
}


/* Compute world position from depth */
int
compute_worldpos(FVECT wpos, DEPTHCODEC *dcp, int x, int y, double d)
{
	RREAL	loc[2];
	FVECT	rdir;

	if (d >= FHUGE*.99)
		goto badval;

	pix2loc(loc, &dcp->res, x, y);

	if (viewray(wpos, rdir, &dcp->vw, loc[0], loc[1]) >= -FTINY) {
		VSUM(wpos, wpos, rdir, d);
		return 1;
	}
badval:
	VCOPY(wpos, dcp->vw.vp);
	return 0;
}


/* Decode the next world position */
int
decode_worldpos_next(FVECT wpos, DEPTHCODEC *dcp)
{
	const long	n = (dcp->curpos - dcp->dstart)>>1;
	int		x, y;
	double		d;

	d = decode_depth_next(dcp);
	if (d < -FTINY)
		return -1;

	x = scanlen(&dcp->res);
	y = n / x;
	x = n - (long)y*x;

	return compute_worldpos(wpos, dcp, x, y, d);
}


/* Seek to the indicated pixel position */
int
seek_dc_pix(DEPTHCODEC *dcp, int x, int y)
{
	long	seekpos;

	if ((dcp->res.xr <= 0) | (dcp->res.yr <= 0)) {
		if (dcp->hdrflags & HF_STDERR) {
			fputs(progname, stderr);
			fputs(": need map resolution to seek\n", stderr);
		}
		return -1;
	}
	if ((x < 0) | (y < 0) ||
			(x >= scanlen(&dcp->res)) | (y >= numscans(&dcp->res))) {
		if (dcp->hdrflags & HF_STDERR) {
			fputs(dcp->inpname, stderr);
			fputs(": warning - pixel index off map\n", stderr);
		}
		return 0;
	}
	seekpos = dcp->dstart + 2*((long)y*scanlen(&dcp->res) + x);

	if (seekpos == dcp->curpos-2) {
		dcp->use_last++;	/* avoids seek/read */
		return 1;
	}
	if (seekpos != dcp->curpos &&
			fseek(dcp->finp, seekpos, SEEK_SET) == EOF) {
		if (dcp->hdrflags & HF_STDERR) {
			fputs(dcp->inpname, stderr);
			fputs(": seek error\n", stderr);
		}
		return -1;
	}
	dcp->curpos = seekpos;
	dcp->use_last = 0;
	return 1;
}


/* Read and decode depth for the given pixel */
double
decode_depth_pix(DEPTHCODEC *dcp, int x, int y)
{
	int	rval = seek_dc_pix(dcp, x, y);

	if (rval < 0)
		return -1.;
	if (!rval)
		return .0;

	return decode_depth_next(dcp);
}


/* Read and decode the world position at the given pixel */
int
get_worldpos_pix(FVECT wpos, DEPTHCODEC *dcp, int x, int y)
{
	double	d = decode_depth_pix(dcp, x, y);

	if (d < -FTINY)
		return -1;
	
	return compute_worldpos(wpos, dcp, x, y, d);
}
