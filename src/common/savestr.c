#ifndef lint
static const char	RCSid[] = "$Id: savestr.c,v 2.7 2003/06/07 12:50:20 schorsch Exp $";
#endif
/*
 *  savestr.c - routines for efficient string storage.
 *
 *	Savestr(s) stores a shared read-only string.
 *  Freestr(s) indicates a client is finished with a string.
 *	All strings must be null-terminated.  There is
 *  no imposed length limit.
 *	Strings stored with savestr(s) can be equated
 *  reliably using their pointer values.
 *	Calls to savestr(s) and freestr(s) should be
 *  balanced (obviously).  The last call to freestr(s)
 *  frees memory associated with the string; it should
 *  never be referenced again.
 *
 *  External symbols declared in standard.h
 */

#include "copyright.h"

#include <string.h>
#include <stdlib.h>

#include "standard.h"

#ifndef  NHASH
#define  NHASH		509		/* hash table size (prime!) */
#endif

typedef struct s_head {
	struct s_head  *next;		/* next in hash list */
	int  nl;			/* links count */
}  S_HEAD;				/* followed by the string itself */

static S_HEAD  *stab[NHASH];

#define  hash(s)	(shash(s)%NHASH)

#define  string(sp)	((char *)((sp)+1))

#define  salloc(str)	(S_HEAD *)malloc(sizeof(S_HEAD)+1+strlen(str))

#define  sfree(sp)	free((void *)(sp))


char *
savestr(str)				/* save a string */
char  *str;
{
	register int  hval;
	register S_HEAD  *sp;

	if (str == NULL)
		return(NULL);
	hval = hash(str);
	for (sp = stab[hval]; sp != NULL; sp = sp->next)
		if (!strcmp(str, string(sp))) {
			sp->nl++;
			return(string(sp));
		}
	if ((sp = salloc(str)) == NULL) {
		eputs("Out of memory in savestr\n");
		quit(1);
	}
	strcpy(string(sp), str);
	sp->nl = 1;
	sp->next = stab[hval];
	stab[hval] = sp;
	return(string(sp));
}


void
freestr(s)				/* free a string */
char  *s;
{
	int  hval;
	register S_HEAD  *spl, *sp;

	if (s == NULL)
		return;
	hval = hash(s);
	for (spl = NULL, sp = stab[hval]; sp != NULL; spl = sp, sp = sp->next)
		if (s == string(sp)) {
			if (--sp->nl > 0)
				return;
			if (spl != NULL)
				spl->next = sp->next;
			else
				stab[hval] = sp->next;
			sfree(sp);
			return;
		}
}


int
shash(s)
register char  *s;
{
	register int  h = 0;

	while (*s)
		h = (h<<1 & 0x7fff) ^ (*s++ & 0xff);
	return(h);
}
