/* RCSid $Id: otspecial.h,v 2.5 2003/06/27 06:53:22 greg Exp $ */
/*
 * Special type flags for objects used in rendering.
 * Depends on definitions in otypes.h
 */
#ifndef _RAD_OTSPECIAL_H_
#define _RAD_OTSPECIAL_H_
#ifdef __cplusplus
extern "C" {
#endif

		/* flag for materials to ignore during irradiance comp. */
#define  T_IRR_IGN	T_SP1

#define  irr_ignore(t)	(ofun[t].flags & T_IRR_IGN)


#ifdef __cplusplus
}
#endif
#endif /* _RAD_OTSPECIAL_H_ */

