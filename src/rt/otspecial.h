/* RCSid $Id: otspecial.h,v 2.4 2003/06/07 00:54:58 schorsch Exp $ */
/*
 * Special type flags for objects used in rendering.
 * Depends on definitions in otypes.h
 */
#ifndef _RAD_OTSPECIAL_H_
#define _RAD_OTSPECIAL_H_
#ifdef __cplusplus
extern "C" {
#endif


#include "copyright.h"

		/* flag for materials to ignore during irradiance comp. */
#define  T_IRR_IGN	T_SP1

#define  irr_ignore(t)	(ofun[t].flags & T_IRR_IGN)


#ifdef __cplusplus
}
#endif
#endif /* _RAD_OTSPECIAL_H_ */

