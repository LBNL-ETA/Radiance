#ifndef lint
static const char	RCSid[] = "$Id: x11findwind.c,v 3.2 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 * find a window by its name under X
 */

#include "copyright.h"

#include <stdio.h>
#include <X11/Xlib.h>


Window
xfindwind(dpy, win, name, depth)
Display	*dpy;
Window	win;
char	*name;
int	depth;
{
	char	*nr;
	Window	rr, pr, *cl;
	Window	wr;
	unsigned int	nc;
	register int	i;

	if (depth == 0)		/* negative depths search all */
		return(None);
	if (!XQueryTree(dpy, win, &rr, &pr, &cl, &nc) || nc == 0)
		return(None);
	wr = None;		/* breadth first search */
	for (i = 0; wr == None && i < nc; i++)
		if (XFetchName(dpy, cl[i], &nr)) {
			if (!strcmp(nr, name))
				wr = cl[i];
			free(nr);
		}
	for (i = 0; wr == None && i < nc; i++)
		wr = xfindwind(dpy, cl[i], name, depth-1);
	XFree((char *)cl);
	return(wr);
}
