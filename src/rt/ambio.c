#ifndef lint
static const char	RCSid[] = "$Id: ambio.c,v 2.14 2020/11/12 21:25:20 greg Exp $";
#endif
/*
 * Read and write portable ambient values
 *
 *  Declarations of external symbols in ambient.h
 */

#include "copyright.h"

#include "ray.h"
#include "ambient.h"


#define  badflt(x)	(((x) < -FHUGE) | ((x) > FHUGE))

#define  badvec(v)	(badflt((v)[0]) | badflt((v)[1]) | badflt((v)[2]))


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
	COLR  clr;

	putint(av->lvl, 1, fp);
	putflt(av->weight, fp);
	putpos(av->pos, fp);
	putint(av->ndir, sizeof(av->ndir), fp);
	putint(av->udir, sizeof(av->udir), fp);
	setcolr(clr, colval(av->val,RED),
			colval(av->val,GRN), colval(av->val,BLU));
	putbinary((char *)clr, sizeof(clr), 1, fp);
	putv2(av->rad, fp);
	putv2(av->gpos, fp);
	putv2(av->gdir, fp);
	putint(av->corral, sizeof(av->corral), fp);
	return(ferror(fp) ? -1 : 0);
}


int
ambvalOK(			/* check consistency of ambient value */
	AMBVAL  *av
)
{
	double	d;

	if (badvec(av->pos)) return(0);
	if (!av->ndir | !av->udir) return(0);
	if ((av->weight <= 0.) | (av->weight > 1.)) return(0);
	if ((av->rad[0] <= 0.) | (av->rad[0] >= FHUGE)) return(0);
	if ((av->rad[1] <= 0.) | (av->rad[1] >= FHUGE)) return(0);
	if (av->rad[0] > av->rad[1]+FTINY) return(0);
	if ((colval(av->val,RED) < 0.) |
			(colval(av->val,RED) >= FHUGE) |
			(colval(av->val,GRN) < 0.) |
			(colval(av->val,GRN) >= FHUGE) |
			(colval(av->val,BLU) < 0.) |
			(colval(av->val,BLU) >= FHUGE)) return(0);
	if (badflt(av->gpos[0]) | badflt(av->gpos[1])) return(0);
	if (badflt(av->gdir[0]) | badflt(av->gdir[1])) return(0);
	return(1);
}


int
readambval(			/* read ambient value from stream */
	AMBVAL  *av,
	FILE  *fp
)
{
	COLR  clr;

	av->lvl = getint(1, fp) & 0xff;
	if (feof(fp))
		return(0);
	av->weight = getflt(fp);
	getpos(av->pos, fp);
	av->ndir = getint(sizeof(av->ndir), fp);
	av->udir = getint(sizeof(av->udir), fp);
	if (getbinary((char *)clr, sizeof(clr), 1, fp) != 1)
		return(0);
	colr_color(av->val, clr);
	getv2(av->rad, fp);
	getv2(av->gpos, fp);
	getv2(av->gdir, fp);
	av->corral = (uint32)getint(sizeof(av->corral), fp);
	return(feof(fp) ? 0 : ambvalOK(av));
}
