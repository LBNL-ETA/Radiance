#ifndef lint
static const char	RCSid[] = "$Id: misc.c,v 1.3 2003/07/14 16:05:45 greg Exp $";
#endif
/*
 *   Miscellaneous functions for meta-files
 */


#include  "meta.h"


char  coms[] = COML;

char  errmsg[128];



int
comndx(c)		/* return index for command */

register int  c;

{
 register char  *cp;

 if (!isalpha(c))
    return(-1);

 for (cp = coms; *cp; cp++)
    if(*cp == c)
       return(cp-coms);

 return(-1);
 }






PRIMITIVE  *
pop(pl)			/* pop top off plist */

register PLIST  *pl;

{
 register PRIMITIVE  *p;

 if ((p = pl->ptop) != NULL)  {
    if ((pl->ptop = p->pnext) == NULL)
       pl->pbot = NULL;
    p->pnext = NULL;
    }

 return(p);
 }




push(p, pl)		/* push primitive onto plist */

register PRIMITIVE  *p;
register PLIST  *pl;

{

 if ((p->pnext = pl->ptop) == NULL)
    pl->pbot = p;
 pl->ptop = p;

 }




add(p, pl)		/* add primitive to plist */

register PRIMITIVE  *p;
register PLIST  *pl;

{

 if (pl->ptop == NULL)
    pl->ptop = p;
 else
    pl->pbot->pnext = p;
 p->pnext = NULL;
 pl->pbot = p;

 }




append(pl1, pl2)		/* append pl1 to the end of pl2 */

register PLIST  *pl1, *pl2;

{

    if (pl1->ptop != NULL) {
        if (pl2->ptop != NULL)
            pl2->pbot->pnext = pl1->ptop;
        else
            pl2->ptop = pl1->ptop;
        pl2->pbot = pl1->pbot;
    }

}




fargs(p)		/* free any arguments p has */

register PRIMITIVE  *p;

{

 if (p->args != NULL) {
    freestr(p->args);
    p->args = NULL;
    }
    
 }



char *
nextscan(start, format, result)		/* scan and advance through string */

register char  *start;
char  *format;
char  *result;

{

    if (start == NULL) return(NULL);
    
    while (isspace(*start)) start++;
    
    if (sscanf(start, format, result) != 1) return(NULL);
    
    while (*start && !isspace(*start)) start++;
    
    return(start);
}



mcopy(p1, p2, n)	/* copy p2 into p1 size n */

register char  *p1, *p2;
register int  n;

{

    while (n--)
	*p1++ = *p2++;

}

