/* RCSid: $Id: x11twind.h,v 2.2 2003/02/22 02:07:29 greg Exp $ */
/*
 *  xtwind.h - header for X text window routines.
 *
 *  Written by G. Ward
 *	10/30/87
 *
 *  Modified for X11 B. V. Smith
 *	9/26/88
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

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

#ifdef NOPROTO

extern TEXTWIND  *xt_open();
extern void	xt_puts();
extern void	xt_putc();
extern void	xt_delete();
extern void	xt_insert();
extern void	xt_redraw();
extern void	xt_clear();
extern void	xt_move();
extern int	xt_cursor();
extern void	xt_close();

#else

extern TEXTWIND	*xt_open(Display *dpy, Window parent,
			int x, int y, int width, int height,
			int bw, unsigned long fore, unsigned long back,
			char *fontname);
extern void	xt_puts(char *s, TEXTWIND *t);
extern void	xt_putc(char c, TEXTWIND *t);
extern void	xt_delete(TEXTWIND *t, int r);
extern void	xt_insert(TEXTWIND *t, int r);
extern void	xt_redraw(TEXTWIND *t);
extern void	xt_clear(TEXTWIND *t);
extern void	xt_move(TEXTWIND *t, int r, int c);
extern int	xt_cursor(TEXTWIND *t, int curs);
extern void	xt_close(TEXTWIND *t);

#endif
