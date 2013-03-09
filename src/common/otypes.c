#ifndef lint
static const char RCSid[] = "$Id: otypes.c,v 2.6 2013/03/09 19:20:31 greg Exp $";
#endif
/*
 * Object type lookup and error reporting
 *
 *  External symbols declared in object.h
 */

#include "copyright.h"

#include  "standard.h"

#include  "object.h"

#include  "otypes.h"


int
otype(				/* get object function number from its name */
	char  *ofname
)
{
	int  i;

	for (i = 0; i < NUMOTYPE; i++)
		if (ofun[i].funame[0] == ofname[0] &&
				!strcmp(ofun[i].funame, ofname))
			return(i);

	return(-1);		/* not found */
}


void
objerror(			/* report error related to object */
	OBJREC  *o,
	int  etyp,
	char  *msg
)
{
	char  msgbuf[512];

	sprintf(msgbuf, "%s for %s \"%s\"",
			msg, ofun[o->otype].funame,
			o->oname!=NULL ? o->oname : "(NULL)");
	error(etyp, msgbuf);
}
