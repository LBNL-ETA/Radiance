#ifndef lint
static const char	RCSid[] = "$Id: segment.c,v 1.2 2003/11/15 02:13:37 schorsch Exp $";
#endif
/*
 *   Routines for segment storage and inclusion for meta-files
 *
 *   12/18/84
 */

#include <string.h>

#include  "rtio.h"
#include  "meta.h"


#define  NODEC  0		/* null declaration (must be 0) */

#define  MAXDEC  2048		/* maximum number of declarations */

#define  NHASH  101		/* hash size (prime) */


static struct declaration {
		    char	*sname;		/* segment name */
		    PLIST	sprims;		/* segment primitives */
		    int		context;	/* declaration containing us */
		    int		nexthash;	/* next in hash list */
		} dectabl[MAXDEC];

static int	hashtabl[NHASH];

static int	dectop = 0;		/* top of declaration table */

static int	curdec = NODEC;		/* current declaration */


static int hash(register char	*s);
static int lookup(register char	*name);



static int
hash(			/* hash a string */
	register char	*s
)

{
    register int	hval = 0;

    while (*s)
	hval += *s++;

    return(hval % NHASH);
}





static int
lookup(		/* find name in declaration table */
	register char	*name
)

{
    register int	curd;

    if (name == NULL)
	error(USER, "illegal null segment name in lookup");

    for (curd = hashtabl[hash(name)]; curd != NODEC; curd = dectabl[curd].nexthash)
	if (strcmp(name, dectabl[curd].sname) == 0)
	    break;

    return(curd);
}



int
inseg(void)			/* return TRUE if currently in a segment */

{

    return(curdec != NODEC);
}



void
openseg(			/* open a new segment */
	char	*name
)

{
    register int  olddec;

    if ((olddec = lookup(name)) != NODEC &&
    		dectabl[olddec].context == curdec)	/* redefined segment */
        plfree(&dectabl[olddec].sprims);

    if (++dectop >= MAXDEC)
	error(SYSTEM, "too many segments in openseg");

    dectabl[dectop].sname = savestr(name);
    				/* save previous context */
    dectabl[dectop].context = curdec;
    curdec = dectop;
    
}




void
segprim(		/* store primitive in current segment */
	register PRIMITIVE  *p
)

{
    register PRIMITIVE  *newp;

    if (!inseg())
        error(SYSTEM, "illegal call to segprim");

    switch (p->com) {

	case POPEN:
	    openseg(p->args);
	    break;
	    
	case PSEG:
	    segment(p, segprim);
	    break;

	case PCLOSE:
	    closeseg();
	    break;

	default:
	    if ((newp = palloc()) == NULL)
		error(SYSTEM, "memory limit exceeded in segprim");
	    mcopy((char *)newp, (char *)p, sizeof(PRIMITIVE));
	    newp->args = savestr(p->args);
	    add(newp, &dectabl[curdec].sprims);
	    break;
    }

}


void
closeseg(void)		/* close the current segment */

{
    register int  i;
    
    if (!inseg())
        error(SYSTEM, "illegal call to closeseg");
        
				/* undefine internal declarations */
    for (i = dectop; i > curdec; i--) {
	hashtabl[hash(dectabl[i].sname)] = dectabl[i].nexthash;
	freestr(dectabl[i].sname);
	plfree(&dectabl[i].sprims);
    }
    dectop = curdec;
				/* define this declaration */
    i = hash(dectabl[curdec].sname);
    dectabl[curdec].nexthash = hashtabl[i];
    hashtabl[i] = curdec;
    				/* return context */
    curdec = dectabl[curdec].context;
}



void
segment(			/* expand segment p */
	PRIMITIVE	*p,
	void	(*funcp)(PRIMITIVE *p)
)

{
    int		decln;
    PRIMITIVE	curp;
    register PRIMITIVE	*sp;

    if ((decln = lookup(p->args)) == NODEC)  {
	sprintf(errmsg, "reference to undefined segment \"%s\" in segment",
						p->args);
	error(USER, errmsg);
    }

    for (sp = dectabl[decln].sprims.ptop; sp != NULL; sp = sp->pnext)

	if (isglob(sp->com))

	    (*funcp)(sp);

	else {

	    switch (curp.com = sp->com) {

		case PSEG:
		case PVSTR:
		case PTFILL:
		case PPFILL:
		    curp.arg0 = (sp->arg0 & 0100) |
				(((sp->arg0 & 060) + p->arg0) & 060);
		    break;

		case PLSEG:
		    if (p->arg0 & 020)
			curp.arg0 = (~sp->arg0 & 0100) | (sp->arg0 & 060);
		    else
			curp.arg0 = sp->arg0 & 0160;
		    break;

		default:
		    curp.arg0 = sp->arg0 & 0160;
		    break;
	    }
	    if (p->arg0 & 014)
		curp.arg0 |= p->arg0 & 014;
	    else
		curp.arg0 |= sp->arg0 & 014;
	    if (p->arg0 & 03)
		curp.arg0 |= p->arg0 & 03;
	    else
		curp.arg0 |= sp->arg0 & 03;

	    curp.xy[XMN] = xlate(XMN, sp, p);
	    curp.xy[YMN] = xlate(YMN, sp, p);
	    curp.xy[XMX] = xlate(XMX, sp, p);
	    curp.xy[YMX] = xlate(YMX, sp, p);
	    curp.args = sp->args;

            (*funcp)(&curp);
	}

}





int
xlate(		/* return extrema from p through px */
	short			extrema,
	PRIMITIVE		*p,
	register PRIMITIVE	*px
)

{
    short	oldex;
    int		val;

    if (((oldex = (extrema + 4 - ((px->arg0 >> 4) & 03))%4) ^ extrema) & 02)
	val = (XYSIZE-1) - p->xy[oldex];
    else
	val = p->xy[oldex];
    
    if (extrema & 01)
	return(CONV(val, px->xy[YMX] - px->xy[YMN]) + px->xy[YMN]);
    else
	return(CONV(val, px->xy[XMX] - px->xy[XMN]) + px->xy[XMN]);
}
