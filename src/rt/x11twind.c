#ifndef lint
static const char	RCSid[] = "$Id: x11twind.c,v 2.8 2003/04/11 19:32:00 greg Exp $";
#endif
/*
 *  x11twind.c - routines for X11 text windows.
 *
 *  Written by G. Ward
 *	10/30/87
 *
 *  Modified for X11 by B. V. Smith
 *	9/26/88
 */

#include "copyright.h"

#include  <stdio.h>

#include  <stdlib.h>

#include  <X11/Xlib.h>

#include  "x11twind.h"

#ifndef  BSD
#define  bzero(d,n)		(void)memset(d,0,n)
#endif

#define checkcurs(t)		if ((t)->cursor) togglecurs(t)

#define restorecurs		checkcurs

/* width/height of a character in fontstruct f */
#define	Width(f)		((f)->max_bounds.rbearing - (f)->min_bounds.lbearing)
#define	Height(f)		((f)->ascent + (f)->descent)
#define YStart(f)		((f)->ascent)

static void  togglecurs();


TEXTWIND *
xt_open(dpy, parent, x, y, width, height, bw, fore, back, fontname)
Display  *dpy;
Window  parent;
int  x, y;
int  width, height;
int  bw;
unsigned long  fore, back;
char  *fontname;
{
	register int  i;
	register TEXTWIND  *t;

	if ((t = (TEXTWIND *)malloc(sizeof(TEXTWIND))) == NULL)
		return(NULL);

	t->dpy = dpy;
	t->w = XCreateSimpleWindow(dpy, parent, x, y, width, height,
					bw, fore, back);
	if (t->w == 0)
		return(NULL);
	XMapWindow(dpy, t->w);

	if ((t->f = XLoadQueryFont(dpy, fontname)) == 0)
		return(NULL);

	/* if (!t->f.fixedwidth)   check for fixedwidth later 
		return(NULL); */

	t->gc = XCreateGC(dpy,t->w,0,NULL);
	XSetState(dpy, t->gc, fore, back, GXcopy, AllPlanes);
	XSetFont(dpy, t->gc, t->f->fid);

	t->nc = (width - LEFTMAR) / 
		Width(t->f);		/* number of columns */
	t->nr = height /
		Height(t->f);		/* number of rows */
	if (t->nc < 1 || t->nr < 1)
		return(NULL);
	if ((t->lp = (char **)calloc(t->nr, sizeof(char *))) == NULL)
		return(NULL);
	for (i = 0; i < t->nr; i++)
		if ((t->lp[i] = calloc(t->nc+1, 1)) == NULL)
			return(NULL);
	t->r = t->c = 0;
	t->cursor = TNOCURS;
	return(t);

}


void
xt_putc(c, t)				/* output a character */
int  c;
register TEXTWIND  *t;
{
	char	ch[2];

	checkcurs(t);
	switch (c) {
	case '\n':
		if (t->r >= t->nr - 1)
			xt_delete(t, 0);	/* scroll up 1 line */
		else if (t->r < t->nr - 1)
			t->r++;
	/* fall through */
	case '\r':
		t->c = 0;
		break;
	case '\b':
		while (t->c < 1 && t->r > 0)
			t->c = strlen(t->lp[--t->r]);
		if (t->c > 0)
			t->c--;
		break;
	default:
		if (t->c >= t->nc)
			xt_putc('\n', t);
		ch[0] = c; ch[1] = '\0';
		XDrawImageString(t->dpy, t->w, t->gc, LEFTMAR+t->c*Width(t->f), 
			YStart(t->f)+t->r*Height(t->f), ch, 1);
		t->lp[t->r][t->c++] = c;
		break;
	}
	restorecurs(t);
}


void
xt_puts(s, t)				/* output a string */
register char  *s;
TEXTWIND  *t;
{
	int	oldcurs;

	oldcurs = xt_cursor(t, TNOCURS);	/* for efficiency */
	while (*s)
		xt_putc(*s++, t);
	xt_cursor(t, oldcurs);
}


void
xt_delete(t, r)				/* delete a line */
register TEXTWIND  *t;
int  r;
{
	char  *cp;
	register int  i;

	if (r < 0 || r >= t->nr)
		return;
	checkcurs(t);
					/* move lines */
	XCopyArea(t->dpy, t->w, t->w, t->gc, LEFTMAR, (r+1)*Height(t->f),
			t->nc*Width(t->f), (t->nr-1-r)*Height(t->f),
			LEFTMAR, r*Height(t->f));
	cp = t->lp[r];
	for (i = r; i < t->nr-1; i++)
		t->lp[i] = t->lp[i+1];
	t->lp[t->nr-1] = cp;
					/* clear bottom */
	XClearArea(t->dpy, t->w, LEFTMAR, (t->nr-1)*Height(t->f),
			t->nc*Width(t->f), Height(t->f),(Bool) 0);

	bzero(cp, t->nc);
	restorecurs(t);			/* should we reposition cursor? */
}


void
xt_insert(t, r)				/* insert a line */
register TEXTWIND  *t;
int  r;
{
	char  *cp;
	register int  i;

	if (r < 0 || r >= t->nr)
		return;
	checkcurs(t);
					/* move lines */
	XCopyArea(t->dpy, t->w, t->w, t->gc, LEFTMAR, r*Height(t->f), LEFTMAR,
			(r+1)*Height(t->f), t->nc*Width(t->f),
			(t->nr-1-r)*Height(t->f));
	cp = t->lp[t->nr-1];
	for (i = t->nr-1; i > r; i--)
		t->lp[i] = t->lp[i-1];
	t->lp[r] = cp;
					/* clear new line */
	XClearArea(t->dpy, t->w, LEFTMAR, r*Height(t->f),
			t->nc*Width(t->f), Height(t->f), (Bool) 0);
	bzero(cp, t->nc);
	restorecurs(t);			/* should we reposition cursor? */
}


void
xt_redraw(t)				/* redraw text window */
register TEXTWIND  *t;
{
	register int  i;

	XClearWindow(t->dpy, t->w);
	for (i = 0; i < t->nr; i++)
	    if (strlen(t->lp[i]) > 0)
		XDrawImageString(t->dpy, t->w, t->gc, LEFTMAR, 
				YStart(t->f)+i*Height(t->f),
				t->lp[i], strlen(t->lp[i]));
	restorecurs(t);
}


void
xt_clear(t)				/* clear text window */
register TEXTWIND  *t;
{
	register int  i;

	XClearWindow(t->dpy, t->w);
	for (i = 0; i < t->nr; i++)
		bzero(t->lp[i], t->nc);
	t->r = t->c = 0;
	restorecurs(t);
}


void
xt_move(t, r, c)			/* move to new position */
register TEXTWIND  *t;
int  r, c;
{
	if (r < 0 || c < 0 || r >= t->nr || c >= t->nc)
		return;
	checkcurs(t);
	t->r = r;
	t->c = c;
	restorecurs(t);
}


int
xt_cursor(t, curs)			/* change cursor */
register TEXTWIND  *t;
register int  curs;
{
	register int	oldcurs;

	if (curs != TNOCURS && curs != TBLKCURS)
		return(-1);
	oldcurs = t->cursor;
	if (curs != oldcurs)
		togglecurs(t);
	t->cursor = curs;
	return(oldcurs);
}


void
xt_close(t)				/* close text window */
register TEXTWIND  *t;
{
	register int  i;

	XFreeFont(t->dpy, t->f);
	XFreeGC(t->dpy,t->gc);
	XDestroyWindow(t->dpy, t->w);
	for (i = 0; i < t->nr; i++)
		free(t->lp[i]);
	free((void *)t->lp);
	free((void *)t);
}


static void
togglecurs(t)
register TEXTWIND  *t;
{
	XSetFunction(t->dpy, t->gc, GXinvert);
	XSetPlaneMask(t->dpy, t->gc, 1L);
	XFillRectangle(t->dpy, t->w, t->gc, 
			t->c*Width(t->f)+LEFTMAR, t->r*Height(t->f),
			Width(t->f), Height(t->f));
	XSetFunction(t->dpy, t->gc, GXcopy);
	XSetPlaneMask(t->dpy, t->gc, ~0L);
}
