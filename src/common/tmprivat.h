/* RCSid $Id: tmprivat.h,v 3.23 2011/03/07 20:49:19 greg Exp $ */
/*
 * Private header file for tone mapping routines.
 */
#ifndef _RAD_TMPRIVAT_H_
#define _RAD_TMPRIVAT_H_

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

#define BRT2SCALE(l2)	(int)(M_LN2*TM_BRTSCALE*(l2) + .5 - ((l2) < 0))

#define HISTEP		16		/* steps in BRTSCALE for each bin */

#define MINBRT		(-16*TM_BRTSCALE)	/* minimum usable brightness */
#define MINLUM		(1.125352e-7)		/* tmLuminance(MINBRT) */

#define HISTI(li)	(((li)-MINBRT)/HISTEP)
#define HISTV(i)	(MINBRT + HISTEP/2 + (i)*HISTEP)

#define LMESLOWER	(5.62e-3)		/* lower mesopic limit */
#define	LMESUPPER	(5.62)			/* upper mesopic limit */
#define BMESLOWER	((int)(-5.18*TM_BRTSCALE-.5))
#define BMESUPPER	((int)(1.73*TM_BRTSCALE+.5))

						/* approximate scotopic lum. */
#define	SCO_rf		0.062
#define SCO_gf		0.608
#define SCO_bf		0.330
#define scotlum(c)	(SCO_rf*(c)[RED] + SCO_gf*(c)[GRN] + SCO_bf*(c)[BLU])
#define normscot(c)	( (	(int32)(SCO_rf*256.+.5)*(c)[RED] + \
				(int32)(SCO_gf*256.+.5)*(c)[GRN] + \
				(int32)(SCO_bf*256.+.5)*(c)[BLU]	) >> 8 )

extern int	tmNewMap(TMstruct *tms);	/* allocate new tone-mapping */

extern int	tmErrorReturn(const char *, TMstruct *, int);

						/* lookup for mesopic scaling */
extern BYTE	tmMesofact[BMESUPPER-BMESLOWER];

extern void	tmMkMesofact(void);			/* build tmMesofact */

#define	returnErr(code)	return(tmErrorReturn(funcName,tms,code))
#define returnOK	return(TM_E_OK)

#define	FEQ(a,b)	((a) < (b)+1e-5 && (b) < (a)+1e-5)

#define	PRIMEQ(p1,p2)	(FEQ((p1)[0][0],(p2)[0][0])&&FEQ((p1)[0][1],(p2)[0][1])\
			&&FEQ((p1)[1][0],(p2)[1][0])&&FEQ((p1)[1][1],(p2)[1][1])\
			&&FEQ((p1)[2][0],(p2)[2][0])&&FEQ((p1)[2][1],(p2)[2][1])\
			&&FEQ((p1)[3][0],(p2)[3][0])&&FEQ((p1)[3][1],(p2)[3][1]))

#ifdef __cplusplus
}
#endif
#endif /* _RAD_TMPRIVAT_H_ */

