#ifndef lint
static const char	RCSid[] = "$Id: devtable.c,v 2.8 2011/10/05 17:20:55 greg Exp $";
#endif
/*
 *  devtable.c - device table for rview.
 */

#include "copyright.h"

#include  <stdio.h>

#include  "color.h"
#include  "driver.h"

#if !defined(HAS_X11) && !defined(HAS_QT)
#define HAS_X11
#endif

#ifdef HAS_X11
extern dr_initf_t x11_init;
char  dev_default[] = "x11";
#else
char  dev_default[] = "qt";
#endif

#ifdef HAS_QT
extern dr_initf_t qt_init;
#endif

struct device  devtable[] = {			/* supported devices */
	{"slave", "Slave driver", slave_init},
#ifdef HAS_X11
	{"x11", "X11 color or greyscale display", x11_init},
	{"x11d", "X11 display using stdin/stdout", x11_init},
#endif
#ifdef HAS_QT
	{"qt", "QT display", qt_init},
#endif
	{0}					/* terminator */
};
