#ifndef lint
static const char	RCSid[] = "$Id: duphead.c,v 2.5 2003/06/05 19:29:34 schorsch Exp $";
#endif
/*
 * Duplicate header on stdout.
 *
 * Externals declared in ray.h
 */

#include "copyright.h"

#include  "standard.h"
#include  "paths.h"


int  headismine = 1;		/* true if header file belongs to me */

static char  *headfname = NULL;	/* temp file name */
static FILE  *headfp = NULL;	/* temp file pointer */


void
headclean()			/* remove header temp file (if one) */
{
	if (headfname == NULL)
		return;
	if (headfp != NULL)
		fclose(headfp);
	if (headismine)
		unlink(headfname);
}


void
openheader()			/* save standard output to header file */
{
	static char  template[] = TEMPLATE;

	headfname = mktemp(template);
	if (freopen(headfname, "w", stdout) == NULL) {
		sprintf(errmsg, "cannot open header file \"%s\"", headfname);
		error(SYSTEM, errmsg);
	}
}


void
dupheader()			/* repeat header on standard output */
{
	register int  c;

	if (headfp == NULL) {
		if ((headfp = fopen(headfname, "r")) == NULL)
			error(SYSTEM, "error reopening header file");
		SET_FILE_BINARY(headfp);
	} else if (fseek(headfp, 0L, 0) < 0)
		error(SYSTEM, "seek error on header file");
	while ((c = getc(headfp)) != EOF)
		putchar(c);
}
