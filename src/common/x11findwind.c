#ifndef lint
static const char	RCSid[] = "$Id: x11findwind.c,v 2.6 2004/07/04 12:08:47 schorsch Exp $";
#endif
/*
 * find a window by its name under X
 */

#include "copyright.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <X11/Xlib.h>


Window
xfindwind(
	Display	*dpy,
	Window	win,
	char	*name,
	int	depth
)
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
