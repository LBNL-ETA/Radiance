/* RCSid: $Id: tmprivat.h,v 3.6 2003/02/22 02:07:22 greg Exp $ */
/*
 * Private header file for tone mapping routines.
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

#ifndef	MEM_PTR
#define	MEM_PTR		void *
#endif
#include	"color.h"
#include	"tonemap.h"

#ifdef __cplusplus
extern "C" {
#endif

				/* required constants */
#ifndef M_LN2
#define M_LN2		0.69314718055994530942
#endif
#ifndef M_LN10
#define M_LN10		2.30258509299404568402
#endif
				/* minimum values and defaults */
#define	MINGAM		0.75
#define DEFGAM		2.2
#define	MINLDDYN	2.
#define DEFLDDYN	32.
#define	MINLDMAX	1.
#define	DEFLDMAX	100.

#define BRT2SCALE(l2)	(int)(M_LN2*TM_BRTSCALE*(l2) + ((l2)>0 ? .5 : -.5))

#define HISTEP		8		/* steps in BRTSCALE for each bin */

#define MINBRT		(-16*TM_BRTSCALE)	/* minimum usable brightness */
#define MINLUM		(1.125352e-7)		/* tmLuminance(MINBRT) */

#define LMESLOWER	(5.62e-3)		/* lower mesopic limit */
#define	LMESUPPER	(5.62)			/* upper mesopic limit */
#if	(TM_BRTSCALE==128)
#define BMESLOWER	(-663)			/* encoded LMESLOWER */
#define BMESUPPER	(221)			/* encoded LMESUPPER */
#else
#define BMESLOWER	((int)(-5.18*TM_BRTSCALE-.5))
#define BMESUPPER	((int)(1.73*TM_BRTSCALE+.5))
#endif

#ifndef	int4
#define	int4		int			/* assume 32-bit integers */
#endif
						/* approximate scotopic lum. */
#define	SCO_rf		0.062
#define SCO_gf		0.608
#define SCO_bf		0.330
#define scotlum(c)	(SCO_rf*(c)[RED] + SCO_gf*(c)[GRN] + SCO_bf*(c)[BLU])
#define normscot(c)	( (	(int4)(SCO_rf*256.+.5)*(c)[RED] + \
				(int4)(SCO_gf*256.+.5)*(c)[GRN] + \
				(int4)(SCO_bf*256.+.5)*(c)[BLU]	) >> 8 )

#ifndef	malloc
MEM_PTR		malloc();
MEM_PTR		calloc();
#endif
extern int	tmErrorReturn();

						/* lookup for mesopic scaling */
extern BYTE	tmMesofact[BMESUPPER-BMESLOWER];

extern void	tmMkMesofact();			/* build tmMesofact */

#define	returnErr(code)	return(tmErrorReturn(funcName,code))
#define returnOK	return(TM_E_OK)

#define	FEQ(a,b)	((a) < (b)+1e-5 && (b) < (a)+1e-5)

#define	PRIMEQ(p1,p2)	(FEQ((p1)[0][0],(p2)[0][0])&&FEQ((p1)[0][1],(p2)[0][1])\
			&&FEQ((p1)[1][0],(p2)[1][0])&&FEQ((p1)[1][1],(p2)[1][1])\
			&&FEQ((p1)[2][0],(p2)[2][0])&&FEQ((p1)[2][1],(p2)[2][1])\
			&&FEQ((p1)[3][0],(p2)[3][0])&&FEQ((p1)[3][1],(p2)[3][1]))

#ifdef __cplusplus
}
#endif
