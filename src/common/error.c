#ifndef lint
static const char	RCSid[] = "$Id: error.c,v 2.12 2023/02/09 18:34:31 greg Exp $";
#endif
/*
 *  error.c - standard error reporting function
 *
 *  External symbols declared in standard.h
 */

#include "copyright.h"

#include  <stdio.h>
#include  <stdlib.h>

#include  "rterror.h"

extern char	*strerror(int);
				/* global list of error actions */
struct erract	erract[NERRS] = ERRACT_INIT;

char  errmsg[2048];		/* global error message buffer */


void
error(int etype, const char *emsg)	/* report error, quit if necessary */
{
	struct erract	*ep;

	if ((etype < 0) | (etype >= NERRS))
		return;
	ep = erract + etype;
	if (ep->pf != NULL) {
		if (ep->pre[0]) (*ep->pf)(ep->pre);
		if (emsg != NULL && emsg[0]) (*ep->pf)(emsg);
		if (etype == SYSTEM && errno > 0) {
			(*ep->pf)(": ");
			(*ep->pf)(strerror(errno));
		}
		(*ep->pf)("\n");
	}
	if (!ep->ec)		/* non-fatal */
		return;
	if (ep->ec < 0)		/* dump core */
		abort();
	quit(ep->ec);		/* quit calls exit after cleanup */
}
