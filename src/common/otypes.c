#ifndef lint
static const char RCSid[] = "$Id: otypes.c,v 2.4 2003/03/10 17:13:29 greg Exp $";
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
otype(ofname)			/* get object function number from its name */
register char  *ofname;
{
	register int  i;

	for (i = 0; i < NUMOTYPE; i++)
		if (!strcmp(ofun[i].funame, ofname))
			return(i);

	return(-1);		/* not found */
}


void
objerror(o, etyp, msg)		/* report error related to object */
OBJREC  *o;
int  etyp;
char  *msg;
{
	char  msgbuf[128];

	sprintf(msgbuf, "%s for %s \"%s\"",
			msg, ofun[o->otype].funame,
			o->oname!=NULL ? o->oname : "(NULL)");
	error(etyp, msgbuf);
}
