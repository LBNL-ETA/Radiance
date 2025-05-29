/* RCSid $Id: otspecial.h,v 2.12 2025/05/29 16:42:28 greg Exp $ */
/*
 * Special type flags for objects used in rendering.
 * Depends on definitions in otypes.h
 */
#ifndef _RAD_OTSPECIAL_H_
#define _RAD_OTSPECIAL_H_

#ifdef __cplusplus
extern "C" {
#endif
		/* test for nominally transparent materials */
#define  T_TRANSP	T_SP1
#define  istransp(m)	(ofun[(m)->otype].flags & T_TRANSP || \
			  (((m)->otype==MAT_WGMDF) & ((m)->oargs.nsargs > 5) \
				&& strcmp((m)->oargs.sarg[5], "0")) || \
			  (((m)->otype==MAT_BRTDF) & ((m)->oargs.nsargs > 5) \
				&& strcmp((m)->oargs.sarg[3], "0") | \
				   strcmp((m)->oargs.sarg[4], "0") | \
				   strcmp((m)->oargs.sarg[5], "0")))

		/* test for opaque (SHADOW) materials */
#define  T_OPAQUE       T_SP2
#define  isopaque(m)    (ofun[(m)->otype].flags & T_OPAQUE || \
			  (((m)->otype==MAT_WGMDF) & ((m)->oargs.nsargs > 5) \
				&& !strcmp((m)->oargs.sarg[5], "0")) || \
			  (((m)->otype==MAT_BRTDF) & ((m)->oargs.nsargs > 5) \
				&& !strcmp((m)->oargs.sarg[3], "0") & \
				   !strcmp((m)->oargs.sarg[4], "0") & \
				   !strcmp((m)->oargs.sarg[5], "0")))

		/* test if we have a BSDF proxy surface */
#define isBSDFproxy(m)	(((m)->otype==MAT_BSDF) & ((m)->oargs.nsargs > 0) \
				&& strcmp((m)->oargs.sarg[0], "0"))

		/* defined in initotypes.c */
extern OBJREC   *findmaterial(OBJREC *o);

#ifdef __cplusplus
}
#endif

#endif /* _RAD_OTSPECIAL_H_ */
