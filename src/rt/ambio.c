#ifndef lint
static const char	RCSid[] = "$Id: ambio.c,v 2.15 2023/11/15 18:02:52 greg Exp $";
#endif
/*
 * Read and write portable ambient values
 *
 *  Declarations of external symbols in ambient.h
 */

#include "copyright.h"

#include "ray.h"
#include "ambient.h"


int	*AMB_CNDX = CNDX;	/* open ambient file RGBE indices */
float	*AMB_WLPART = WLPART;	/* open ambient file limits+partitions (nm) */


#define  badflt(x)	(((x) < -FHUGE) | ((x) > FHUGE))

#define  badvec(v)	(badflt((v)[0]) | badflt((v)[1]) | badflt((v)[2]))


int
amb_headline(		/* check ambient headline (copy to file if *p) */
	char *hl,
	void *p
)
{
	static int	ambcndx[4];
	static float	ambwlpart[4];
	FILE		*fout = (FILE *)p;
	char		fmt[MAXFMTLEN];

	if (formatval(fmt, hl)) {
		if (strcmp(fmt, AMBFMT))
			return(-1);
		return(0);
	}
	if (isncomp(hl)) {
		ambcndx[3] = ncompval(hl);
		if ((ambcndx[3] < 3) | (ambcndx[3] > MAXCSAMP))
			return(-1);
		AMB_CNDX = ambcndx;
		return(1);
	}
	if (iswlsplit(hl)) {
		if (!wlsplitval(ambwlpart, hl))
			return(-1);
		AMB_WLPART = ambwlpart;
		return(1);
	}
	if (fout)
		fputs(hl, fout);
	return(0);
}


void
putambmagic(fp)			/* write out ambient value magic number */
FILE  *fp;
{
	putint(AMBMAGIC, 2, fp);
}


int
hasambmagic(fp)			/* read in and check validity of magic # */
FILE  *fp;
{
	int  magic;

	magic = getint(2, fp);
	if (feof(fp))
		return(0);
	return(magic == AMBMAGIC);
}


#define  putpos(v,fp)	putflt((v)[0],fp);putflt((v)[1],fp);putflt((v)[2],fp)

#define  getpos(v,fp)	(v)[0]=getflt(fp);(v)[1]=getflt(fp);(v)[2]=getflt(fp)

#define  putv2(v2,fp)	putflt((v2)[0],fp);putflt((v2)[1],fp)

#define  getv2(v2,fp)	(v2)[0]=getflt(fp);(v2)[1]=getflt(fp)

int
writambval(			/* write ambient value to stream */
	AMBVAL  *av,
	FILE  *fp
)
{
	SCOLR  sclr;

	putint(av->lvl, 1, fp);
	putflt(av->weight, fp);
	putpos(av->pos, fp);
	putint(av->ndir, sizeof(av->ndir), fp);
	putint(av->udir, sizeof(av->udir), fp);
	putv2(av->rad, fp);
	putv2(av->gpos, fp);
	putv2(av->gdir, fp);
	putint(av->corral, sizeof(av->corral), fp);
	if ((AMB_CNDX == CNDX) & (AMB_WLPART == WLPART)) {
		scolor_scolr(sclr, av->val);
	} else {
		SCOLOR	scol;
		convertscolor(scol, AMB_CNDX[3], AMB_WLPART[0], AMB_WLPART[3],
				av->val, NCSAMP, WLPART[0], WLPART[3]);
		scolor2scolr(sclr, scol, AMB_CNDX[3]);
	}
	putbinary(sclr, AMB_CNDX[3]+1, 1, fp);

	return(ferror(fp) ? -1 : 0);
}


int
ambvalOK(			/* check consistency of ambient value */
	AMBVAL  *av
)
{
	int	i;

	if (badvec(av->pos)) return(0);
	if (!av->ndir | !av->udir) return(0);
	if ((av->weight <= 0.) | (av->weight > 1.)) return(0);
	if ((av->rad[0] <= 0.) | (av->rad[0] >= FHUGE)) return(0);
	if ((av->rad[1] <= 0.) | (av->rad[1] >= FHUGE)) return(0);
	if (av->rad[0] > av->rad[1]+FTINY) return(0);
	if (badflt(av->gpos[0]) | badflt(av->gpos[1])) return(0);
	if (badflt(av->gdir[0]) | badflt(av->gdir[1])) return(0);
	for (i = AMB_CNDX[3]; i-- > 0; )
		if ((av->val[i] < 0.) | (av->val[i] >= FHUGE))
			return(0);
	return(1);
}


int
readambval(			/* read ambient value from stream */
	AMBVAL  *av,
	FILE  *fp
)
{
	SCOLR  sclr;

	av->lvl = getint(1, fp) & 0xff;
	if (feof(fp))
		return(0);
	av->weight = getflt(fp);
	getpos(av->pos, fp);
	av->ndir = getint(sizeof(av->ndir), fp);
	av->udir = getint(sizeof(av->udir), fp);
	getv2(av->rad, fp);
	getv2(av->gpos, fp);
	getv2(av->gdir, fp);
	av->corral = (uint32)getint(sizeof(av->corral), fp);
	if (getbinary(sclr, AMB_CNDX[3]+1, 1, fp) != 1)
		return(0);
	if ((AMB_CNDX == CNDX) & (AMB_WLPART == WLPART)) {
		scolr_scolor(av->val, sclr);
	} else {
		SCOLOR	scol;
		scolr2scolor(scol, sclr, AMB_CNDX[3]);
		convertscolor(av->val, NCSAMP, WLPART[0], WLPART[3],
				scol, AMB_CNDX[3], AMB_WLPART[0], AMB_WLPART[3]);
	}
	return(feof(fp) ? 0 : ambvalOK(av));
}
