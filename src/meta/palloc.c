#ifndef lint
static const char	RCSid[] = "$Id: palloc.c,v 1.2 2003/06/08 12:03:10 schorsch Exp $";
#endif
/*
 *   Limited dynamic storage allocator for primitives
 */


#define  FBSIZE  72		/* size of malloc call */


#include  "meta.h"


extern int  maxalloc;		/* number of prims to allocate */
int  nalloc = 0;		/* number allocated so far */

static PLIST  freelist = {NULL, NULL};


static int morefree(void);


PRIMITIVE *
palloc(void)		/* allocate a primitive */

{
 register PRIMITIVE  *p;

 if (maxalloc > 0 && nalloc >= maxalloc)
    return(NULL);

 if ((p = pop(&freelist)) == NULL)
    if (morefree())
       p = pop(&freelist);
    else {
       sprintf(errmsg, "out of memory in palloc (nalloc = %d)", nalloc);
       error(SYSTEM, errmsg);
       }

 nalloc++;
 return(p);
}



void
pfree(		/* free a primitive */
register PRIMITIVE  *p
)
{

 if (p->args != NULL) {
    freestr(p->args);
    p->args = NULL;
    } 
 push(p, &freelist);
 nalloc--;

 }



void
plfree(		/* free a primitive list */
register PLIST  *pl
)
{
    register PRIMITIVE  *p;
    
    for (p = pl->ptop; p != NULL; p = p->pnext) {
        if (p->args != NULL) {
           freestr(p->args);
           p->args = NULL;
           }
        nalloc--;
        }

    append(pl, &freelist);
    pl->ptop = pl->pbot = NULL;

}
    


static int
morefree(void)		/* get more free space */

{
 register PRIMITIVE  *p;
 register int  i;
 int rnu;

 if (maxalloc > 0 && (i = maxalloc-nalloc) < FBSIZE)
    rnu = i;
 else
    rnu = i = FBSIZE;

 p = (PRIMITIVE *)malloc((unsigned)i * sizeof(PRIMITIVE));

 if (p == NULL)
    return(0);

 while (i--) {
    p->args = NULL;
    push(p, &freelist);
    p++;
    }

 return(rnu);
 }
