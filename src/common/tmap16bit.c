#ifndef lint
static const char RCSid[] = "$Id: tmap16bit.c,v 1.2 2003/07/15 00:23:36 greg Exp $";
#endif
/*
 * Routines for tone-mapping 16-bit/primary pixels
 *
 * Externals declared in tonemap.h
 */

#include "copyright.h"

#include <stdio.h>
#include <math.h>
#include "tmprivat.h"

#define LOGTABBITS	11	/* log table is 1<<LOGTABBITS long */
#define GAMTABBITS	9	/* gamma table is 1<<GAMTABBITS long */

static float	logtab[1<<LOGTABBITS];
static float	gamtab[1<<GAMTABBITS];
static float	gammul[16];
static double	cur_gam = 0.;

#define imultpow2(i,s)	((s)>=0 ? (i)<<(s) : (i)>>-(s))


/* Fill our log table */
static void
mkLogTable()
{	
	int	i;

	if (logtab[0] != 0.f)
		return;
	for (i = 1<<LOGTABBITS; i--; )
		logtab[i] = log(imultpow2(i,15-LOGTABBITS)*(1./(1L<<16))
					+ .5);
}


/* Fill our gamma table */
static void
mkGamTable(double gv)
{
	int	i;

	if (gv == cur_gam)
		return;
	for (i = 1<<GAMTABBITS; i--; )
		gamtab[i] = pow((i+.5)*(1./(1<<GAMTABBITS)), gv);
	for (i = 16; i--; )
		gammul[i] = pow((double)(1L<<i), -gv);
	cur_gam = gv;
}


/* Find normalizing shift value for a 2-byte unsigned integer */
static int
normShift16(int i)
{
	int	s = 0;
	
	if (!i)
		return(-1);
	while (!(i & 0x8000)) {
		i <<= 1;
		++s;
	}
	return(s);
}


/* Find common normalizing shift for 3 2-byte unsigned integers */
static int
normShift48(uint16 clr48[3])
{
	int	imax = (clr48[1] > clr48[0] ? clr48[1] : clr48[0]);
	if (clr48[2] > imax)
		imax = clr48[2];
	return(normShift16(imax));
}


/* convert at 48-bit tristimulus value to a COLOR */
static void
rgb48_color(COLOR col, uint16 clr48[3], double gv)
{
	int	nshft;

	if (gv == 1.) {				/* linear case */
		col[0] = clr48[0]*(1./(1L<<16));
		col[1] = clr48[1]*(1./(1L<<16));
		col[2] = clr48[2]*(1./(1L<<16));
		return;
	}
						/* non-linear case */
	/* XXX Uncomment if routine is made public
	if (gv != cur_gam)
		mkGamTable(gv);
	*/
	nshft = normShift48(clr48);
	if (nshft < 0) {
		col[0] = col[1] = col[2] = .0f;
		return;
	}
	col[0] = gamtab[imultpow2(clr48[0],GAMTABBITS-16+nshft)] * 
			gammul[nshft];
	col[1] = gamtab[imultpow2(clr48[1],GAMTABBITS-16+nshft)] *
			gammul[nshft];
	col[2] = gamtab[imultpow2(clr48[2],GAMTABBITS-16+nshft)] *
			gammul[nshft];
}


/* Convert 16-bit gray scanline to encoded luminance */
int
tmCvGray16(TMbright *ls, uint16 *scan, int len, double gv)
{
	static char	funcName[] = "tmCvGray16";
	static double	cur_inpsf = 1.;
	static double	log_inpsf = 0.;
	int		nshft;
	double		d;

	if (tmTop == NULL)
		returnErr(TM_E_TMINVAL);
	if (ls == NULL | scan == NULL | len < 0)
		returnErr(TM_E_ILLEGAL);
	if (gv <= 0.)
		gv = DEFGAM;
						/* initialize log table */
	if (logtab[0] == 0.f)
		mkLogTable();
	if (cur_inpsf != tmTop->inpsf)
		log_inpsf = log(cur_inpsf = tmTop->inpsf);
						/* convert 16-bit grays */
	while (len--) {
		nshft = normShift16(*scan);
		if (nshft < 0) {		/* bogus value */
			*ls++ = TM_NOBRT;
			scan++;
			continue;
		}
		d = logtab[ imultpow2(*scan,LOGTABBITS-15+nshft) &
					((1L<<LOGTABBITS)-1) ]
				- M_LN2*nshft;
		d = (double)TM_BRTSCALE * (gv*d + log_inpsf);
		*ls++ = (d>0. ? d+.5 : d-.5);
		scan++;
	}
	returnOK;
}

/* Convert a 48-bit RGB scanline to encoded luminance/chrominance */
int
tmCvRGB48(TMbright *ls, BYTE *cs, uint16 (*scan)[3], int len, double gv)
{
	static char	funcName[] = "tmCvRGB48";
	static double	cur_inpsf = 1.;
	static double	log_inpsf = 0.;
	int		i;

	if (tmTop == NULL)
		returnErr(TM_E_TMINVAL);
	if (ls == NULL | scan == NULL | len < 0)
		returnErr(TM_E_ILLEGAL);
	if (gv <= 0.)
		gv = DEFGAM;
						/* update gamma table */
	if (gv != 1. & gv != cur_gam)
		mkGamTable(gv);
#if 0
	if (tmNeedMatrix(tmTop)) {		/* need floating point */
#else
	{
#endif
		COLOR	*newscan;
		newscan = (COLOR *)tempbuffer(len*sizeof(COLOR));
		if (newscan == NULL)
			returnErr(TM_E_NOMEM);
		for (i = len; i--; )
			rgb48_color(newscan[i], scan[i], gv);
		return(tmCvColors(ls, cs, newscan, len));
	}
#if 0
						/* initialize log table */
	if (logtab[0] == 0.f)
		mkLogTable();
	if (cur_inpsf != tmTop->inpsf)
		log_inpsf = log(cur_inpsf = tmTop->inpsf);
						/* convert scanline */
	for (i = len; i--; ) {
		copycolr(cmon, scan[i]);
							/* world luminance */
		li =  ( cd->clfb[RED]*cmon[RED] +
			cd->clfb[GRN]*cmon[GRN] +
			cd->clfb[BLU]*cmon[BLU] ) >> 8;
		bi = BRT2SCALE(cmon[EXP]-COLXS) +
				logi[li] + cd->inpsfb;
		if (li <= 0) {
			bi = TM_NOBRT;			/* bogus value */
			li = 1;				/* avoid li==0 */
		}
		ls[i] = bi;
		if (cs == TM_NOCHROM)			/* no color? */
			continue;
							/* mesopic adj. */
		if (tmTop->flags & TM_F_MESOPIC && bi < BMESUPPER) {
			int	pf, sli = normscot(cmon);
			if (bi < BMESLOWER)
				cmon[RED] = cmon[GRN] = cmon[BLU] = sli;
			else {
				if (tmTop->flags & TM_F_BW)
					cmon[RED] = cmon[GRN] = cmon[BLU] = li;
				pf = tmMesofact[bi-BMESLOWER];
				sli *= 256 - pf;
				cmon[RED] = ( sli + pf*cmon[RED] ) >> 8;
				cmon[GRN] = ( sli + pf*cmon[GRN] ) >> 8;
				cmon[BLU] = ( sli + pf*cmon[BLU] ) >> 8;
			}
		} else if (tmTop->flags & TM_F_BW) {
			cmon[RED] = cmon[GRN] = cmon[BLU] = li;
		}
		bi = ( (int32)GAMTSZ*cd->clfb[RED]*cmon[RED]/li ) >> 8;
		cs[3*i  ] = bi>=GAMTSZ ? 255 : cd->gamb[bi];
		bi = ( (int32)GAMTSZ*cd->clfb[GRN]*cmon[GRN]/li ) >> 8;
		cs[3*i+1] = bi>=GAMTSZ ? 255 : cd->gamb[bi];
		bi = ( (int32)GAMTSZ*cd->clfb[BLU]*cmon[BLU]/li ) >> 8;
		cs[3*i+2] = bi>=GAMTSZ ? 255 : cd->gamb[bi];
	}
	returnOK;
#endif
}
