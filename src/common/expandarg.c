#ifndef lint
static const char	RCSid[] = "$Id: expandarg.c,v 2.8 2003/06/30 14:59:11 schorsch Exp $";
#endif
/*
 * Get additional command arguments from file or environment.
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include <errno.h>
#include <string.h>

#include "rtio.h"
#include "rtmisc.h"


#define MAXARGEXP	512		/* maximum argument expansion */

			/* set the following to suit, -1 to disable */
int	envexpchr = '$';		/* environment expansion character */
int	filexpchr = '@';		/* file expansion character */


int
expandarg(acp, avp, n)		/* expand list at argument n */
int	*acp;
register char	***avp;
int	n;
{
	int	ace;
	char	*ave[MAXARGEXP];
	char	**newav;
					/* check argument */
	if (n >= *acp)
		return(0);
	errno = 0;	
	if ((*avp)[n][0] == filexpchr) {		/* file name */
		ace = wordfile(ave, (*avp)[n]+1);
		if (ace < 0)
			return(-1);	/* no such file */
	} else if ((*avp)[n][0] == envexpchr) {		/* env. variable */
		ace = wordstring(ave, getenv((*avp)[n]+1));
		if (ace < 0)
			return(-1);	/* no such variable */
	} else						/* regular argument */
		return(0);
					/* allocate new pointer list */
	newav = (char **)bmalloc((*acp+ace)*sizeof(char *));
	if (newav == NULL)
		return(-1);
					/* copy preceeding arguments */
	memcpy((void *)newav, (void *)*avp, n*sizeof(char *));
					/* copy expanded argument */
	memcpy((void *)(newav+n), (void *)ave, ace*sizeof(char *));
					/* copy trailing arguments + NULL */
	memcpy((void *)(newav+n+ace), (void *)(*avp+n+1), (*acp-n)*sizeof(char *));
					/* free old list */
	bfree((char *)*avp, (*acp+1)*sizeof(char *));
					/* assign new list */
	*acp += ace-1;
	*avp = newav;
	return(1);			/* return success */
}
