#ifndef lint
static const char	RCSid[] = "$Id: ambio.c,v 2.5 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 * Read and write portable ambient values
 *
 *  Declarations of external symbols in ambient.h
 */

#include "copyright.h"

#include "ray.h"

#include "ambient.h"


#define  putvec(v,fp)	putflt((v)[0],fp);putflt((v)[1],fp);putflt((v)[2],fp)

#define  getvec(v,fp)	(v)[0]=getflt(fp);(v)[1]=getflt(fp);(v)[2]=getflt(fp)

#define  badflt(x)	((x) < -FHUGE || (x) > FHUGE)

#define  badvec(v)	(badflt((v)[0]) || badflt((v)[1]) || badflt((v)[2]))



void
putambmagic(fp)			/* write out ambient value magic number */
FILE  *fp;
{
	putint((long)AMBMAGIC, 2, fp);
}


int
hasambmagic(fp)			/* read in and check validity of magic # */
FILE  *fp;
{
	register int  magic;

	magic = getint(2, fp);
	if (feof(fp))
		return(0);
	return(magic == AMBMAGIC);
}


int
writambval(av, fp)		/* write ambient value to stream */
register AMBVAL  *av;
FILE  *fp;
{
	COLR  col;

	putint((long)av->lvl, 1, fp);
	putflt(av->weight, fp);
	putvec(av->pos, fp);
	putvec(av->dir, fp);
	setcolr(col, colval(av->val,RED),
			colval(av->val,GRN), colval(av->val,BLU));
	fwrite((char *)col, sizeof(col), 1, fp);
	putflt(av->rad, fp);
	putvec(av->gpos, fp);
	putvec(av->gdir, fp);
	return(ferror(fp) ? -1 : 0);
}


int
ambvalOK(av)			/* check consistency of ambient value */
register AMBVAL  *av;
{
	double	d;

	if (badvec(av->pos)) return(0);
	if (badvec(av->dir)) return(0);
	d = DOT(av->dir,av->dir);
	if (d < 0.9999 || d > 1.0001) return(0);
	if (av->lvl < 0 || av->lvl > 100) return(0);
	if (av->weight <= 0. || av->weight > 1.) return(0);
	if (av->rad <= 0. || av->rad >= FHUGE) return(0);
	if (colval(av->val,RED) < 0. ||
			colval(av->val,RED) > FHUGE ||
			colval(av->val,GRN) < 0. ||
			colval(av->val,GRN) > FHUGE ||
			colval(av->val,BLU) < 0. ||
			colval(av->val,BLU) > FHUGE) return(0);
	if (badvec(av->gpos)) return(0);
	if (badvec(av->gdir)) return(0);
	return(1);
}


int
readambval(av, fp)		/* read ambient value from stream */
register AMBVAL  *av;
FILE  *fp;
{
	COLR  col;

	av->lvl = getint(1, fp);
	if (feof(fp))
		return(0);
	av->weight = getflt(fp);
	getvec(av->pos, fp);
	getvec(av->dir, fp);
	if (fread((char *)col, sizeof(col), 1, fp) != 1)
		return(0);
	colr_color(av->val, col);
	av->rad = getflt(fp);
	getvec(av->gpos, fp);
	getvec(av->gdir, fp);
	return(feof(fp) ? 0 : ambvalOK(av));
}
