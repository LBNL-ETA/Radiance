/* RCSid $Id: otspecial.h,v 2.3 2003/02/25 02:47:23 greg Exp $ */
/*
 * Special type flags for objects used in rendering.
 * Depends on definitions in otypes.h
 */

#include "copyright.h"

		/* flag for materials to ignore during irradiance comp. */
#define  T_IRR_IGN	T_SP1

#define  irr_ignore(t)	(ofun[t].flags & T_IRR_IGN)
