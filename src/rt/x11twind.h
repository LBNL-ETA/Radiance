/* RCSid $Id: x11twind.h,v 2.5 2003/06/07 00:54:58 schorsch Exp $ */
/*
 *  xtwind.h - header for X text window routines.
 *
 *  Written by G. Ward
 *	10/30/87
 *
 *  Modified for X11 B. V. Smith
 *	9/26/88
 */
#ifndef _RAD_X11TWIND_H_
#define _RAD_X11TWIND_H_
#ifdef __cplusplus
extern "C" {
#endif


#include "copyright.h"

#define LEFTMAR	2			/* left margin width */

typedef struct {
	Display  *dpy;			/* the display */
	Window  w;			/* window */
	XFontStruct  *f;		/* font information */
	GC	gc;			/* graphics context */
	short  nc, nr;			/* text number of columns and rows */
	char  **lp;			/* null-terminated lines */
	short  c, r;			/* current position */
	short	cursor;			/* cursor type */
}  TEXTWIND;			/* a text window */

#define TNOCURS		0
#define TBLKCURS	1

extern TEXTWIND	*xt_open(Display *dpy, Window parent,
			int x, int y, int width, int height,
			int bw, unsigned long fore, unsigned long back,
			char *fontname);
extern void	xt_puts(char *s, TEXTWIND *t);
extern void	xt_putc(int c, TEXTWIND *t);
extern void	xt_delete(TEXTWIND *t, int r);
extern void	xt_insert(TEXTWIND *t, int r);
extern void	xt_redraw(TEXTWIND *t);
extern void	xt_clear(TEXTWIND *t);
extern void	xt_move(TEXTWIND *t, int r, int c);
extern int	xt_cursor(TEXTWIND *t, int curs);
extern void	xt_close(TEXTWIND *t);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_X11TWIND_H_ */

