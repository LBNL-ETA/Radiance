#ifndef lint
static const char	RCSid[] = "$Id: ambio.c,v 2.4 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 * Read and write portable ambient values
 *
 *  Declarations of external symbols in ambient.h
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

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
