#ifndef lint
static const char	RCSid[] = "$Id: devtable.c,v 2.10 2018/04/06 17:52:33 greg Exp $";
#endif
/*
 *  devtable.c - device table for rview.
 */

#include "copyright.h"

#include  "driver.h"

#if !defined(HAS_X11) && !defined(HAS_QT) && !defined(WIN_RVIEW)
/* weird logic ... */
#define HAS_X11 1
#endif

#if HAS_X11
extern dr_initf_t x11_init;
char  dev_default[] = "x11";
#elif defined(HAS_QT)
char  dev_default[] = "qt";
#elif defined(WIN_RVIEW)
char  dev_default[] = "win";
#endif

#ifdef HAS_QT
extern dr_initf_t qt_init;
#endif

#ifdef WIN_RVIEW
extern dr_initf_t win_rvudev_init;
#endif

struct device  devtable[] = {			/* supported devices */
	{"slave", "Slave driver", slave_init},
#if HAS_X11
	{"x11", "X11 color or greyscale display", x11_init},
	{"x11d", "X11 display using stdin/stdout", x11_init},
#endif
#ifdef HAS_QT
	{"qt", "QT display", qt_init},
#endif
#ifdef WIN_RVIEW
	{"win", "Windows display", win_rvudev_init},
#endif
	{0}					/* terminator */
};
