/* RCSid $Id: otspecial.h,v 2.6 2003/12/31 01:50:02 greg Exp $ */
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

		/* flag for completely opaque materials */
#define  T_OPAQUE       T_SP2

#define  irr_ignore(t)	(ofun[t].flags & T_IRR_IGN)


#ifdef __cplusplus
}
#endif
#endif /* _RAD_OTSPECIAL_H_ */

