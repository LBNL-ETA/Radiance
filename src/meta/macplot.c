#ifndef lint
static const char	RCSid[] = "$Id: macplot.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  Plotting routines for MacIntosh
 */


#include  "meta.h"

#include  "macplot.h"


#define  NPATS  9		/* number of fill patterns */


int  pati[4] = {0, 1, 2, 3};	/* pattern indices */

/*
 *  Patterns are 8 X 8, ordered left column to right, one byte per column.
 */

Pattern  macpat[NPATS+4] = {
	{0377,0377,0377,0377,0377,0377,0377,0377},	/* solid */
	{0217,0307,0343,0361,0370,0174,076,037},	/* thick \\\ (dashed) */
	{0231,0314,0146,063,0231,0314,0146,063},	/* thin \\\ (dotted) */
	{0201,0300,0140,060,030,014,06,03},		/* sparse \\\ */
	{037,076,0174,0370,0361,0343,0307,0217},	/* thick /// (dashed) */
	{063,0146,0314,0231,063,0146,0314,0231},	/* thin /// (dotted) */
	{03,06,014,030,060,0140,0300,0201},		/* sparse /// */
	{0201,0102,044,030,030,044,0102,0201},		/* crisscross */
	{0201,0101,041,021,011,05,03,0377}		/* web */
    };

static int  context = 0;

struct setting {
	int  cnx;		/* setting context */
	char  *val;		/* attribute value */
	struct setting  *snext;	/* next setting in list */
};

static struct setting  *sets[0200];




set(attrib, value)		/* set attribute to value */

int  attrib;
char  *value;

{
    struct setting  *newset;

    switch (attrib) {
	case SALL:
	    context++;
	    return;
	case SPAT0:
	case SPAT1:
	case SPAT2:
	case SPAT3:
	    if (!spat(attrib, value)) {
		sprintf(errmsg, "Bad pattern set value: %s", value);
		error(WARNING, errmsg);
		return;
	    }
	    break;
	default:
	    sprintf(errmsg, "Can't set attribute: %o", attrib);
	    error(WARNING, errmsg);
	    return;
    }
    newset = (struct setting *)malloc(sizeof(struct setting));
    newset->cnx = context;
    newset->val = savestr(value);
    newset->snext = sets[attrib];
    sets[attrib] = newset;
}




unset(attrib)			/* return attribute to previous setting */

int  attrib;

{
    register int  i;

    if (attrib == SALL) {
	if (context == 0)
	    reset(SALL);
	else {
	    context--;
	    for (i = 0; i < 0200; i++)
		while (sets[i] != NULL && sets[i]->cnx > context)
		    spop(i);
	}
	return;
    }
    spop(attrib);
    if (sets[attrib] == NULL)
	reset(attrib);

    switch (attrib) {
	case SPAT0:
	case SPAT1:
	case SPAT2:
	case SPAT3:
	    spat(attrib, sets[attrib]->val);
	    break;
	default:
	    sprintf(errmsg, "Can't unset attribute: %o", attrib);
	    error(WARNING, errmsg);
	    return;
    }
}





reset(attrib)			/* return attribute to default setting */

int  attrib;

{
    switch (attrib) {
	case SALL:
	    reset(SPAT0);
	    reset(SPAT1);
	    reset(SPAT2);
	    reset(SPAT3);
	    context = 0;
	    return;
	case SPAT0:
	    spat(SPAT0, "P0");
	    break;
	case SPAT1:
	    spat(SPAT1, "P1");
	    break;
	case SPAT2:
	    spat(SPAT2, "P2");
	    break;
	case SPAT3:
	    spat(SPAT3, "P3");
	    break;
	default:
	    sprintf(errmsg, "Can't reset attribute: %o", attrib);
	    error(WARNING, errmsg);
	    return;
    }
    while (sets[attrib] != NULL)
	spop(attrib);
}




static
spop(attrib)			/* pop top off attrib settings list */

int  attrib;

{

    if (sets[attrib] != NULL) {
	if (sets[attrib]->val != NULL)
	    free(sets[attrib]->val);
	free((char *)sets[attrib]);
	sets[attrib] = sets[attrib]->snext;
    }

}




static int
spat(pat, patval)			/* set a pattern */

int  pat;
char  *patval;

{
    int  n, v;
    char  *nextscan();
    register char  *cp;
    register int  j;

    if (patval == NULL) return(FALSE);

    if (patval[0] == 'P' || patval[0] == 'p') {
	if (nextscan(patval+1, "%d", &n) == NULL || n < 0 || n >= NPATS)
	    return(FALSE);
    } else {
	n = NPATS + pat - SPAT0;
	cp = patval;
	for (j = 0; j < 8; j++) {
	    if ((cp = nextscan(cp, "%o", &v)) == NULL || v < 0 || v > 0377)
		return(FALSE);
	    macpat[n][j] = v;
	}
    }
    pati[pat-SPAT0] = n;
    return(TRUE);
}
