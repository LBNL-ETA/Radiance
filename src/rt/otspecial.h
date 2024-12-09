/* RCSid $Id: otspecial.h,v 2.11 2024/12/09 00:44:29 greg Exp $ */
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
				&& strcmp((m)->oargs.sarg[5], "0")))

		/* test for completely opaque materials */
#define  T_OPAQUE       T_SP2
#define  isopaque(m)    (ofun[(m)->otype].flags & T_OPAQUE || \
			  (((m)->otype==MAT_WGMDF) & ((m)->oargs.nsargs > 5) \
				&& !strcmp((m)->oargs.sarg[5], "0")))

		/* test if we have a BSDF proxy surface */
#define isBSDFproxy(m)	(((m)->otype==MAT_BSDF) & ((m)->oargs.nsargs > 0) \
				&& strcmp((m)->oargs.sarg[0], "0"))

		/* defined in initotypes.c */
extern OBJREC   *findmaterial(OBJREC *o);

#ifdef __cplusplus
}
#endif

#endif /* _RAD_OTSPECIAL_H_ */
