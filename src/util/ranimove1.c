#ifndef lint
static const char	RCSid[] = "$Id: ranimove1.c,v 3.21 2012/10/05 00:59:38 greg Exp $";
#endif
/*
 *  ranimove1.c
 *
 *  Basic frame rendering routines for ranimove(1).
 *
 *  Created by Gregory Ward on Wed Jan 08 2003.
 */

#include "copyright.h"

#include <string.h>

#include "platform.h"
#include "ranimove.h"
#include "otypes.h"
#include "source.h"
#include "random.h"

double		acctab[256];	/* accuracy value table */

int		hres, vres;	/* frame resolution (fcur) */
double		pixaspect;	/* pixel aspect ratio */

VIEW		vw;		/* view for this frame */
COLOR		*cbuffer;	/* color at each pixel */
float		*zbuffer;	/* depth at each pixel */
OBJECT		*obuffer;	/* object id at each pixel */
short		*xmbuffer;	/* x motion at each pixel */
short		*ymbuffer;	/* y motion at each pixel */
uby8		*abuffer;	/* accuracy at each pixel */
uby8		*sbuffer;	/* sample count per pixel */
COLOR		*outbuffer;	/* output buffer (may equal cbuffer) */

VIEW		vwprev;		/* last frame's view */
COLOR		*cprev;		/* last frame colors */
float		*zprev;		/* last frame depth */
OBJECT		*oprev;		/* last frame objects */
uby8		*aprev;		/* last frame accuracy */

float		*cerrmap;	/* conspicuous error map */
COLOR		*val2map;	/* value-squared map for variance */

double		frm_stop;	/* when to stop rendering this frame */

double		hlsmax;		/* maximum high-level saliency this frame */


static void next_frame(void);
static int sample_here(int x, int y);
static int offset_cmp(const void *p1, const void *p2);
static void setmotion(int n, FVECT wpos);
static void init_frame_sample(void);


#if 0
void
write_map(		/* write out float map (debugging) */
	float	*mp,
	char	*fn
)
{
	FILE	*fp = fopen(fn, "w");
	COLOR	scanbuf[2048];
	int	x, y;

	if (fp == NULL)
		return;
	newheader("RADIANCE", fp);
	fputformat(COLRFMT, fp);
	fputc('\n', fp);		/* end header */
	fprtresolu(hres, vres, fp);
	for (y = vres; y--; ) {		/* write scanlines */
		float	*bp = mp + (y+1)*hres - 1;
		for (x = hres; x--; bp--)
			setcolor(scanbuf[x], *bp, *bp, *bp);
		if (fwritescan(scanbuf, hres, fp) < 0)
			break;
	}
	fclose(fp);
}
#endif


static void
next_frame(void)			/* prepare next frame buffer */
{
	VIEW	*fv;
	char	*err;
					/* get previous view */
	if (vw.type != 0)
		vwprev = vw;
	else if (fcur > 1 && (fv = getview(fcur-1)) != NULL) {
		vwprev = *fv;
		if (setview(&vwprev) != NULL)
			vwprev.type = 0;
	}
					/* get current view */
	if ((fv = getview(fcur)) == NULL) {
		sprintf(errmsg, "cannot get view for frame %d", fcur);
		error(USER, errmsg);
	}
	vw = *fv;
	if ((err = setview(&vw)) != NULL) {
		sprintf(errmsg, "view error at frame %d: %s", fcur, err);
		error(USER, errmsg);
	}
	if (cbuffer == NULL) {
		int	n;		/* compute resolution and allocate */
		switch (sscanf(vval(RESOLUTION), "%d %d %lf",
				&hres, &vres, &pixaspect)) {
		case 1:
			vres = hres;
			/* fall through */
		case 2:
			pixaspect = 1.;
			/* fall through */
		case 3:
			if ((hres > 0) & (vres > 0))
				break;
			/* fall through */
		default:
			sprintf(errmsg, "bad %s value", vnam(RESOLUTION));
			error(USER, errmsg);
		}
		normaspect(viewaspect(&vw), &pixaspect, &hres, &vres);
		cbuffer = (COLOR *)malloc(sizeof(COLOR)*hres*vres);
		zbuffer = (float *)malloc(sizeof(float)*hres*vres);
		obuffer = (OBJECT *)malloc(sizeof(OBJECT)*hres*vres);
		xmbuffer = (short *)malloc(sizeof(short)*hres*vres);
		ymbuffer = (short *)malloc(sizeof(short)*hres*vres);
		abuffer = (uby8 *)calloc(hres*vres, sizeof(uby8));
		sbuffer = (uby8 *)calloc(hres*vres, sizeof(uby8));
		cprev = (COLOR *)malloc(sizeof(COLOR)*hres*vres);
		zprev = (float *)malloc(sizeof(float)*hres*vres);
		oprev = (OBJECT *)malloc(sizeof(OBJECT)*hres*vres);
		aprev = (uby8 *)malloc(sizeof(uby8)*hres*vres);
		if (mblur > .02)
			outbuffer = (COLOR *)malloc(sizeof(COLOR)*hres*vres);
		else
			outbuffer = cbuffer;
		if ((cbuffer==NULL) | (zbuffer==NULL) | (obuffer==NULL) |
				(xmbuffer==NULL) | (ymbuffer==NULL) |
				(abuffer==NULL) | (sbuffer==NULL) |
				(cprev==NULL) | (zprev == NULL) |
				(oprev==NULL) | (aprev==NULL) |
				(outbuffer==NULL))
			error(SYSTEM, "out of memory in init_frame");
		for (n = hres*vres; n--; ) {
			zprev[n] = -1.f;
 			oprev[n] = OVOID;
		}
		frm_stop = getTime() + rtperfrm;
	} else {			/* else just swap buffers */
		float	*fp;
		OBJECT	*op;
		uby8	*bp;
		if (outbuffer != cbuffer) {
			COLOR	*cp = cprev;
			cprev = cbuffer; cbuffer = cp;
		} else {
			outbuffer = cprev; cprev = cbuffer;
			cbuffer = outbuffer;
		}
		fp = zprev; zprev = zbuffer; zbuffer = fp;
		op = oprev; oprev = obuffer; obuffer = op;
		bp = aprev; aprev = abuffer; abuffer = bp;
		memset(abuffer, 0, sizeof(uby8)*hres*vres);
		memset(sbuffer, 0, sizeof(uby8)*hres*vres);
		frm_stop += rtperfrm;
	}
	cerrmap = NULL;
	val2map = NULL;
}


#define SAMPDIST	3	/* Maximum distance to neighbor sample */
#define	SAMPDIST2	(SAMPDIST*SAMPDIST)


static int
sample_here(		/* 4x4 quincunx sample at this pixel? */
	int	x,
	int	y
)
{
	if (y & 0x1)		/* every other row has samples */
		return(0);
	if (y & 0x3)		/* every fourth row is offset */
		x += 2;
	return((x & 0x3) == 0);	/* every fourth column is sampled */
}


void
sample_pos(	/* compute jittered sample position */
	double	hv[2],
	int	x,
	int	y,
	int	sn
)
{
	int	hl[2];

	hl[0] = x; hl[1] = y;
	multisamp(hv, 2, urand(ilhash(hl,2) + sn));
	hv[0] = ((double)x + hv[0]) / (double)hres;
	hv[1] = ((double)y + hv[1]) / (double)vres;
}


double
sample_wt(		/* compute interpolant sample weight */
	int	xo,
	int yo
)
{
	static double	etab[400];
	/* we can't use the name rad2 here, for some reason Visual C
	   thinks that is a constant (compiler bug?) */
	int	rad_2 = xo*xo + yo*yo;
	int	i;

	if (etab[0] <= FTINY)		/* initialize exponent table */
		for (i = 400; i--; )
			etab[i] = exp(-0.1*i);

					/* look up Gaussian */
	i = (int)((10.*3./(double)SAMPDIST2)*rad_2 + .5);
	if (i >= 400)
		return(0.0);
	return(etab[i]);
}


static int
offset_cmp(		/* compare offset distances */
	const void	*p1,
	const void	*p2
)
{
	return(*(const int *)p1 - *(const int *)p2);
}


int
getclosest(	/* get nc closest neighbors on same object */
	int	*iarr,
	int	nc,
	int	x,
	int	y
)
{
#define	NSCHECK		((2*SAMPDIST+1)*(2*SAMPDIST+1))
	static int	hro, vro;
	static int	ioffs[NSCHECK];
	OBJECT	myobj;
	int	i0, nf;
	int	i, j;
					/* get our object number */
	myobj = obuffer[fndx(x, y)];
					/* special case for borders */
	if ((x < SAMPDIST) | (x >= hres-SAMPDIST) |
			(y < SAMPDIST) | (y >= vres-SAMPDIST)) {
		int	tndx[NSCHECK][2];
		nf = 0;
		for (j = y - SAMPDIST; j <= y + SAMPDIST; j++) {
		    if (j >= vres) break;
		    if (j < 0) j = 0;
		    for (i = x - SAMPDIST; i <= x + SAMPDIST; i++) {
			if (i >= hres) break;
			if (i < 0) i = 0;
			i0 = fndx(i, j);
			if (!sbuffer[i0])
				continue;
			if ((myobj != OVOID) & (obuffer[i0] != myobj))
				continue;
			tndx[nf][0] = (i-x)*(i-x) + (j-y)*(j-y);
			tndx[nf][1] = i0;
			nf++;
		    }
		}
		qsort((void *)tndx, nf, 2*sizeof(int), offset_cmp);
		if (nf > nc)
			nf = nc;
		for (i = nf; i--; )
			iarr[i] = tndx[i][1];
		return(nf);
	}
					/* initialize offset array */
	if ((hres != hro) | (vres != vro)) {
		int	toffs[NSCHECK][2];
		i0 = fndx(SAMPDIST, SAMPDIST);
		nf = 0;
		for (i = 0; i <= 2*SAMPDIST; i++)
		    for (j = 0; j <= 2*SAMPDIST; j++) {
			toffs[nf][0] = (i-SAMPDIST)*(i-SAMPDIST) +
					(j-SAMPDIST)*(j-SAMPDIST);
			toffs[nf][1] = fndx(i, j) - i0;
			nf++;
		    }
		qsort((void *)toffs, nf, 2*sizeof(int), offset_cmp);
		for (i = NSCHECK; i--; )
			ioffs[i] = toffs[i][1];
		hro = hres;
		vro = vres;
	}
					/* find up to nc neighbors */
	i0 = fndx(x, y);
	for (j = 0, nf = 0; (j < NSCHECK) & (nf < nc); j++) {
		i = i0 + ioffs[j];
		if (sbuffer[i] && (myobj == OVOID) | (obuffer[i] == myobj))
			iarr[nf++] = i;
	}
					/* return number found */
	return(nf);
#undef NSCHECK
}


static void
setmotion(		/* compute motion vector for this pixel */
		int	n,
		FVECT	wpos
)
{
	FVECT	ovp;
	int	moi;
	int	xp, yp;
					/* ID object and update maximum HLS */
	moi = getmove(obuffer[n]);
	if (moi >= 0 && obj_move[moi].cprio > hlsmax)
		hlsmax = obj_move[moi].cprio;
	if (vwprev.type == 0)		/* return leaves MO_UNK */
		return;
	if (moi >= 0) {			/* move object point back */
		multp3(ovp, wpos, obj_move[moi].bxfm);
		wpos = ovp;
	}
	viewloc(ovp, &vwprev, wpos);
	if (ovp[2] <= FTINY)
		return;
	xp = (int)(ovp[0]*hres);
	yp = (int)(ovp[1]*vres);
	xmbuffer[n] = xp - (n % hres);
	ymbuffer[n] = yp - (n / hres);
	if ((xp < 0) | (xp >= hres))
		return;
	if ((yp < 0) | (yp >= vres))
		return;
	n = fndx(xp, yp);
	if ((zprev[n] < 0.97*ovp[2]) | (zprev[n] > 1.03*ovp[2]))
		oprev[n] = OVOID;	/* assume it's a bad match */
}


static void
init_frame_sample(void)		/* sample our initial frame */
{
	RAY	ir;
	int	x, y;
	int	n;

	if (!silent) {
		printf("\tComputing initial samples...");
		fflush(stdout);
	}
	hlsmax = CSF_SMN;
	for (y = vres; y--; )
	    for (x = hres; x--; ) {
		double	hv[2];
		n = fndx(x, y);
		xmbuffer[n] = ymbuffer[n] = MO_UNK;
		sample_pos(hv, x, y, 0);
		ir.rmax = viewray(ir.rorg, ir.rdir, &vw, hv[0], hv[1]);
		if (ir.rmax < -FTINY) {
			setcolor(cbuffer[n], 0., 0., 0.);
			zbuffer[n] = FHUGE;
			obuffer[n] = OVOID;
			abuffer[n] = ADISTANT;
			continue;
		}
		if (!sample_here(x, y)) {	/* just cast */
			rayorigin(&ir, PRIMARY, NULL, NULL);
			if (!localhit(&ir, &thescene)) {
				if (ir.ro != &Aftplane && sourcehit(&ir)) {
					rayshade(&ir, ir.ro->omod);
					rayparticipate(&ir);
				}
				copycolor(cbuffer[n], ir.rcol);
				zbuffer[n] = ir.rot;
				obuffer[n] = ir.robj;
				abuffer[n] = ADISTANT;
				sbuffer[n] = 1;
			} else {
				zbuffer[n] = ir.rot;
				obuffer[n] = ir.robj;
				setmotion(n, ir.rop);
			}	
			continue;
		}
		if (nprocs > 1) {		/* get sample */
			int	rval;
			rayorigin(&ir, PRIMARY, NULL, NULL);
			ir.rno = n;
			rval = ray_pqueue(&ir);
			if (!rval)
				continue;
			if (rval < 0)
				quit(1);
			n = ir.rno;
		} else
			ray_trace(&ir);
		copycolor(cbuffer[n], ir.rcol);
		zbuffer[n] = ir.rot;
		obuffer[n] = ir.robj;
		sbuffer[n] = 1;
		if (ir.rot >= 0.99*FHUGE)
			abuffer[n] = ADISTANT;
		else {
			abuffer[n] = ALOWQ;
			setmotion(n, ir.rop);
		}
	    }
	if (nprocs > 1)			/* get stragglers */
		while (ray_presult(&ir, 0)) {
			n = ir.rno;
			copycolor(cbuffer[n], ir.rcol);
			zbuffer[n] = ir.rot;
			obuffer[n] = ir.robj;
			sbuffer[n] = 1;
			if (ir.rot >= FHUGE)
				abuffer[n] = ADISTANT;
			else {
				abuffer[n] = ALOWQ;
				setmotion(n, ir.rop);
			}
		}
					/* ambiguate object boundaries */
	for (y = vres-1; y--; )
	    for (x = hres-1; x--; ) {
		OBJECT	obj;
		n = fndx(x, y);
		if ((obj = obuffer[n]) == OVOID)
			continue;
		if ((obuffer[n+1] != OVOID) & (obuffer[n+1] != obj)) {
			obuffer[n] = OVOID;
			obuffer[n+1] = OVOID;
		}
		if ((obuffer[n+hres] != OVOID) & (obuffer[n+hres] != obj)) {
			obuffer[n] = OVOID;
			obuffer[n+hres] = OVOID;
		}
	    }

	if (!silent)
		printf("done\n");
}


int
getambcolor(		/* get ambient color for object if we can */
		COLOR	clr,
		int	obj
)
{
	OBJREC	*op;

	if (obj == OVOID)
		return(0);
	op = objptr(obj);		/* search for material */
	if (op->omod == OVOID)
		return(0);
	op = findmaterial(objptr(op->omod));
	if (op == NULL)
		return(0);
	/*
	 * Since this routine is called to compute the difference
	 * from rendering with and without interreflections,
	 * we don't want to return colors for materials that are
	 * explicitly excluded from the HQ ambient calculation.
	 */
	if (hirendparams.ambincl >= 0) {
		int	i;
		char	*lv;
		for (i = 0; (lv = rpambmod(&hirendparams,i)) != NULL; i++)
			if (lv[0] == op->oname[0] &&
					!strcmp(lv+1, op->oname+1))
				break;
		if ((lv != NULL) ^ hirendparams.ambincl)
			return(0);
	}
	switch (op->otype) {
	case MAT_PLASTIC:
	case MAT_METAL:
	case MAT_PLASTIC2:
	case MAT_METAL2:
	case MAT_PFUNC:
	case MAT_MFUNC:
	case MAT_PDATA:
	case MAT_MDATA:
	case MAT_TRANS:
	case MAT_TRANS2:
	case MAT_TFUNC:
	case MAT_TDATA:
		if (op->oargs.nfargs < 3)
			return(0);
		setcolor(clr, op->oargs.farg[0], op->oargs.farg[1],
				op->oargs.farg[2]);
		return(1);
	case MAT_BRTDF:
		if (op->oargs.nfargs < 6)
			return(0);
		setcolor(clr, op->oargs.farg[0]+op->oargs.farg[3],
				op->oargs.farg[1]+op->oargs.farg[4],
				op->oargs.farg[2]+op->oargs.farg[5]);
		scalecolor(clr, 0.5);
		return(1);
	case MAT_LIGHT:
	case MAT_GLOW:
	case MAT_ILLUM:
		setcolor(clr, 0., 0., 0.);
		return(1);
	}
	return(0);
}


double
estimaterr(		/* estimate relative error from samples */
		COLOR	cs,
		COLOR	cs2,
		int	ns,
		int ns0
)
{
	double	d, d2, brt;

	if (ns <= 1 || (brt = bright(cs)/ns) < 1e-14)
		return(1.0);
					/* use largest of RGB std. dev. */
	d2 = colval(cs2,RED) - colval(cs,RED)*colval(cs,RED)/ns;
	d = colval(cs2,GRN) - colval(cs,GRN)*colval(cs,GRN)/ns;
	if (d > d2) d2 = d;
	d = colval(cs2,BLU) - colval(cs,BLU)*colval(cs,BLU)/ns;
	if (d > d2) d2 = d;
					/* use s.d. if <= 1 central sample */
	if (ns0 <= 1)
		return(sqrt(d2/(ns-1))/brt);
					/* use s.d./sqrt(ns0) otherwise */
	return(sqrt(d2/((ns-1)*ns0))/brt);
}


double
comperr(		/* estimate relative error in neighborhood */
		int	*neigh,
		int	nc,
		int	ns0
)
{
	COLOR	csum, csum2;
	COLOR	ctmp;
	int	i;
	int	ns;
	int	n;
					/* add together samples */
	setcolor(csum, 0., 0., 0.);
	setcolor(csum2, 0., 0., 0.);
	for (i = 0, ns = 0; (i < nc) & (ns < NSAMPOK); i++) {
		n = neigh[i];
		addcolor(csum, cbuffer[n]);
		if (val2map != NULL) {
			addcolor(csum2, val2map[n]);
			ns += sbuffer[n];
			continue;
		}
		if (sbuffer[n] != 1)
			error(CONSISTENCY, "bad count in comperr");
		copycolor(ctmp, cbuffer[n]);
		multcolor(ctmp, ctmp);
		addcolor(csum2, ctmp);
		ns++;
	}
	return(estimaterr(csum, csum2, ns, ns0));
}


void
comp_frame_error(void)		/* initialize frame error values */
{
	uby8	*edone = NULL;
	COLOR	objamb;
	double	eest;
	int	neigh[NSAMPOK];
	int	nc;
	int	x, y, i;
	int	n;

	if (!silent) {
		printf("\tComputing error map\n");
		fflush(stdout);
	}
	if (acctab[0] <= FTINY)		/* initialize accuracy table */
		for (i = 256; i--; )
			acctab[i] = errorf(i);
					/* estimate sample error */
	if (!curparams->ambounce && hirendparams.ambounce) {
		/*
		 * Our error estimate for the initial value is based
		 * on the assumption that most of it comes from the
		 * lack of an interreflection calculation.  The relative
		 * error should be less than the ambient value divided
		 * by the returned ray value -- we take half of this.
		 */
		edone = (uby8 *)calloc(hres*vres, sizeof(uby8));
		for (n = hres*vres; n--; ) {
			if ((abuffer[n] != ALOWQ) | (obuffer[n] == OVOID))
				continue;
			if (!getambcolor(objamb, obuffer[n]))
				continue;
			multcolor(objamb, ambval);
			if ((eest = bright(cbuffer[n])) <= FTINY)
				continue;
			eest = bright(objamb) / eest;
			if (eest > 1.)	/* should we report this? */
				continue;
			eest *= 0.50;	/* use 50% ambient error */
			i = errori(eest);
			if (i < AMIN) i = AMIN;
			else if (i >= ADISTANT/2) i = ADISTANT/2-1;
			abuffer[n] = i;
			edone[n] = 1;
		}
	}
					/* final statistical estimate */
	for (y = vres; y--; )
	    for (x = hres; x--; ) {
		n = fndx(x, y);
		if (abuffer[n] == ADISTANT)
			continue;	/* don't update these */
		if (edone != NULL && edone[n])
			continue;	/* already done this */
		if (sbuffer[n] >= 255) {
			abuffer[n] = ADISTANT;
			continue;	/* can't take any more */
		}
		nc = getclosest(neigh, NSAMPOK, x, y);
		if (nc <= 0) {
			abuffer[n] = ANOVAL;
			continue;	/* no clue what to do for him */
		}
		i = errori(comperr(neigh, nc, sbuffer[n]));
		if (i < AMIN) i = AMIN;
		else if (i >= ADISTANT) i = ADISTANT-1;
		abuffer[n] = i;
					/* can't be better than closest */
		if (i < abuffer[neigh[0]] && abuffer[neigh[0]] >= AMIN)
			abuffer[n] = abuffer[neigh[0]];
	    }
	if (edone != NULL)
		free((void *)edone);
}


void
init_frame(void)			/* render base (low quality) frame */
{
	int	restart;
					/* allocate/swap buffers */
	next_frame();
					/* check rendering status */
	restart = (!nobjects || vdef(MOVE));
	if (!restart && curparams != &lorendparams && nprocs > 1)
		restart = -1;
					/* post low quality parameters */
	if (curparams != &lorendparams)
		ray_restore(curparams = &lorendparams);
	if (restart > 0) {		/* load new octree */
		char	*oct = getoctspec(fcur);
		if (oct == NULL) {
			sprintf(errmsg, "cannot get scene for frame %d", fcur);
			error(USER, errmsg);
		}
		if (!silent) {
			printf("\tLoading octree...");
			fflush(stdout);
		}
		if (nprocs > 1)
			ray_pinit(oct, nprocs);
		else
			ray_init(oct);
	} else if (restart < 0) {	/* update children */
		if (!silent) {
			printf("\tRestarting %d processes...", nprocs);
			fflush(stdout);
		}
		ray_pclose(0);
		ray_popen(nprocs);
	}
	if (restart && !silent)
		printf("done\n");
					/* sample frame buffer */
	init_frame_sample();
					/* initialize frame error */
	comp_frame_error();
#if 0
{
	float	*ebuf = (float *)malloc(sizeof(float)*hres*vres);
	char	fnm[256];
	int	n;
	for (n = hres*vres; n--; )
		ebuf[n] = acctab[abuffer[n]];
	sprintf(fnm, vval(BASENAME), fcur);
	strcat(fnm, "_inerr.pic");
	write_map(ebuf, fnm);
	free((void *)ebuf);
}
#endif
}


void
filter_frame(void)			/* interpolation, motion-blur, and exposure */
{
	const double	expval = expspec_val(getexp(fcur));
	int		x, y;
	int		neigh[NPINTERP];
	int		nc;
	COLOR		cval;
	double		w, wsum;
	int		n;

#if 0
	/* XXX TEMPORARY!! */
	conspicuity();
	write_map(cerrmap, "outcmap.pic");
	{
		float	*ebuf = (float *)malloc(sizeof(float)*hres*vres);
		for (n = hres*vres; n--; )
			ebuf[n] = acctab[abuffer[n]];
		write_map(ebuf, "outerr.pic");
		free((void *)ebuf);
	}
#endif
	if (!silent) {
		printf("\tFiltering frame\n");
		fflush(stdout);
	}
					/* normalize samples */
	for (n = hres*vres; n--; ) {
		if (sbuffer[n] <= 1)
			continue;
		w = 1.0/(double)sbuffer[n];
		scalecolor(cbuffer[n], w);
	}
					/* interpolate samples */
	for (y = vres; y--; )
	    for (x = hres; x--; ) {
		n = fndx(x, y);
		if (sbuffer[n])
			continue;
		nc = getclosest(neigh, NPINTERP, x, y);
		setcolor(cbuffer[n], 0., 0., 0.);
		if (nc <= 0) {		/* no acceptable neighbors */
			if (y < vres-1)
				nc = fndx(x, y+1);
			else if (x < hres-1)
				nc = fndx(x+1, y);
			else
				continue;
			copycolor(cbuffer[n], cbuffer[nc]);
			continue;
		}
		wsum = 0.;
		while (nc-- > 0) {
			copycolor(cval, cbuffer[neigh[nc]]);
			w = sample_wt((neigh[nc]%hres) - x,
					(neigh[nc]/hres) - y);
			scalecolor(cval, w);
			addcolor(cbuffer[n], cval);
			wsum += w;
		}
		w = 1.0/wsum;
		scalecolor(cbuffer[n], w);
	    }
					/* motion blur if requested */
	if (outbuffer != cbuffer) {
		int	xs, ys, xl, yl;
		int	rise, run;
		long	rise2, run2;
		int	n2;
		int	cnt;
fprintf(stderr, "outbuffer=0x%lx, cbuffer=0x%lx, cprev=0x%lx\n",
(size_t)outbuffer, (size_t)cbuffer, (size_t)cprev);
					/* sum in motion streaks */
		memset(outbuffer, 0, sizeof(COLOR)*hres*vres);
		memset(wbuffer, 0, sizeof(float)*hres*vres);
		for (y = vres; y--; )
		    for (x = hres; x--; ) {
			n = fndx(x, y);
			if (xmbuffer[n] == MO_UNK) {
				run = rise = 0;
			} else {
				run = (int)(mblur*xmbuffer[n]);
				rise = (int)(mblur*ymbuffer[n]);
			}
			if (!(run | rise)) {
				addcolor(outbuffer[n], cbuffer[n]);
				wbuffer[n] += 1.;
				continue;
			}
			xl = x; yl = y;
			if (run < 0) { xs = -1; run = -run; }
			else xs = 1;
			if (rise < 0) { ys = -1; rise = -rise; }
			else ys = 1;
			rise2 = run2 = 0L;
			if (rise > run) {
				cnt = rise + 1;
				w = 1./cnt;
				copycolor(cval, cbuffer[n]);
				scalecolor(cval, w);
				while (cnt)
					if (rise2 >= run2) {
						if ((xl >= 0) & (xl < hres) &
						    (yl >= 0) & (yl < vres)) {
							n2 = fndx(xl, yl);
							addcolor(outbuffer[n2],
									cval);
							wbuffer[n2] += w;
						}
						yl += ys;
						run2 += run;
						cnt--;
					} else {
						xl += xs;
						rise2 += rise;
					}
			} else {
				cnt = run + 1;
				w = 1./cnt;
				copycolor(cval, cbuffer[n]);
				scalecolor(cval, w);
				while (cnt)
					if (run2 >= rise2) {
						if ((xl >= 0) & (xl < hres) &
						    (yl >= 0) & (yl < vres)) {
							n2 = fndx(xl, yl);
							addcolor(outbuffer[n2],
									cval);
							wbuffer[n2] += w;
						}
						xl += xs;
						rise2 += rise;
						cnt--;
					} else {
						yl += ys;
						run2 += run;
					}
			}
		    }
					/* compute final results */
		for (n = hres*vres; n--; )
			if (wbuffer[n] > 1.02) {
				w = 1./wbuffer[n];
				scalecolor(outbuffer[n], w);
			} else if (wbuffer[n] < 0.98) {
				w = 1.-wbuffer[n];
				copycolor(cval, cprev[n]);
				scalecolor(cval, w);
				addcolor(outbuffer[n], cval);
			}
	}
					/* adjust exposure */
	for (n = ((expval < 0.98) | (expval > 1.02))*hres*vres; n--; )
		scalecolor(outbuffer[n], expval);
	/*
	   for (n = hres*vres; n--; )
		   if (!sbuffer[n])
			   setcolor(outbuffer[n], 0., 0., 0.);
	 */
#if 0
	{
		float	*sbuf = (float *)malloc(sizeof(float)*hres*vres);
		char	fnm[256];
		sprintf(fnm, vval(BASENAME), fcur);
		strcat(fnm, "_outsamp.pic");
		for (n = hres*vres; n--; )
			sbuf[n] = (float)sbuffer[n];
		write_map(sbuf, fnm);
		free((void *)sbuf);
	}
#endif
}


void
send_frame(void)			/* send frame to destination */
{
	char	fname[1024];
	double	d;
	FILE	*fp;
	int	y;
					/* open output picture */
	sprintf(fname, vval(BASENAME), fcur);
	strcat(fname, ".hdr");
	fp = fopen(fname, "w");
	if (fp == NULL) {
		sprintf(errmsg, "cannot open output frame \"%s\"", fname);
		error(SYSTEM, errmsg);
	}
	SET_FILE_BINARY(fp);
	if (!silent) {
		printf("\tWriting to \"%s\"\n", fname);
		fflush(stdout);
	}
					/* write header */
	newheader("RADIANCE", fp);
	printargs(gargc, gargv, fp);
	fprintf(fp, "SOFTWARE= %s\n", VersionID);
	fprintf(fp, "FRAME=%d\n", fcur);
	fputnow(fp);
	fputs(VIEWSTR, fp); fprintview(&vw, fp); fputc('\n', fp);
	d = expspec_val(getexp(fcur));
	if ((d < 0.99) | (d > 1.01))
		fputexpos(d, fp);
	d = viewaspect(&vw) * hres / vres;
	if ((d < 0.99) | (d > 1.01))
		fputaspect(d, fp);
	fputformat(COLRFMT, fp);
	fputc('\n', fp);		/* end header */
	fprtresolu(hres, vres, fp);
	if (fflush(fp) == EOF)
		goto writerr;
#if (PIXSTANDARD != (YMAJOR|YDECR))
	error(CONSISTENCY, "bad code in send_frame");
#endif
	for (y = vres; y--; )		/* write scanlines */
		if (fwritescan(outbuffer+y*hres, hres, fp) < 0)
			goto writerr;
	if (fclose(fp) == EOF)
		goto writerr;
	if (vdef(ZNAME)) {		/* output z-buffer */
		sprintf(fname, vval(ZNAME), fcur);
		strcat(fname, ".zbf");
		fp = fopen(fname, "w");
		if (fp == NULL) {
			sprintf(errmsg, "cannot open z-file \"%s\"\n", fname);
			error(SYSTEM, errmsg);
		}
		SET_FILE_BINARY(fp);
		if (!silent) {
			printf("\tWriting depths to \"%s\"\n", fname);
			fflush(stdout);
		}
		for (y = vres; y--; )
			if (fwrite(zbuffer+y*hres, sizeof(float),
						hres, fp) != hres)
				goto writerr;
		if (fclose(fp) == EOF)
			goto writerr;
	}
	if (vdef(MNAME)) {		/* output motion buffer */
		unsigned short	*mbuffer = (unsigned short *)malloc(
						sizeof(unsigned short)*3*hres);
		int		x, n;
		if (mbuffer == NULL)
			error(SYSTEM, "out of memory in send_frame");
		sprintf(fname, vval(MNAME), fcur);
		strcat(fname, ".mvo");
		fp = fopen(fname, "w");
		if (fp == NULL) {
			sprintf(errmsg, "cannot open motion file \"%s\"\n",
							fname);
			error(SYSTEM, errmsg);
		}
		SET_FILE_BINARY(fp);
		if (!silent) {
			printf("\tWriting motion vectors to \"%s\"\n", fname);
			fflush(stdout);
		}
		for (y = vres; y--; ) {
			for (x = hres; x--; ) {
				n = fndx(x,y);
				mbuffer[3*x] = xmbuffer[n] + 0x8000;
				mbuffer[3*x+1] = ymbuffer[n] + 0x8000;
				mbuffer[3*x+2] = (oprev[n]!=OVOID)*0x8000;
			}
			if (fwrite(mbuffer, sizeof(*mbuffer),
						3*hres, fp) != 3*hres)
				goto writerr;
		}
		free((void *)mbuffer);
		if (fclose(fp) == EOF)
			goto writerr;
	}
	return;				/* all is well */
writerr:
	sprintf(errmsg, "error writing file \"%s\"", fname);
	error(SYSTEM, errmsg);
}


void
free_frame(void)			/* free frame allocation */
{
	if (cbuffer == NULL)
		return;
	if (outbuffer != cbuffer)
		free((void *)outbuffer);
	outbuffer = NULL;
	free((void *)cbuffer); cbuffer = NULL;
	free((void *)zbuffer); zbuffer = NULL;
	free((void *)obuffer); obuffer = NULL;
	free((void *)xmbuffer); xmbuffer = NULL;
	free((void *)ymbuffer); ymbuffer = NULL;
	free((void *)cprev); cprev = NULL;
	free((void *)zprev); zprev = NULL;
	free((void *)oprev); oprev = NULL;
	cerrmap = NULL;
	val2map = NULL;
	hres = vres = 0;
	vw.type = vwprev.type = 0;
	frm_stop = 0;
}
