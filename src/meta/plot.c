#ifndef lint
static const char	RCSid[] = "$Id: plot.c,v 1.2 2003/06/08 12:03:10 schorsch Exp $";
#endif
/*
 *  Plotting routines for meta-files
 */


#include  "meta.h"

#include  "plot.h"


int  pati[4] = {0, 1, 2, 3};	/* pattern indices */

/*
 *  Patterns are 16 X 16, ordered left to right, bottom to top.
 *  Each byte represents a column in an 8-bit high row.
 */

unsigned char	pattern[NPATS+4][PATSIZE/8][PATSIZE] = {
    {						/* solid */
	{
	    0377,0377,0377,0377,0377,0377,0377,0377,
	    0377,0377,0377,0377,0377,0377,0377,0377
	}, {
	    0377,0377,0377,0377,0377,0377,0377,0377,
	    0377,0377,0377,0377,0377,0377,0377,0377
	}
    }, {					/* thick \\\ (dashed) */
	{
	    0377,0177,077,037,017,0207,0303,0341,
	    0360,0370,0374,0376,0377,0377,0377,0377
	}, {
	    0360,0370,0374,0376,0377,0377,0377,0377,
	    0377,0177,077,037,017,0207,0303,0341
	}
    }, {					/* thin \\\ (dotted) */
	{
	    0314,0146,063,0231,0314,0146,063,0231,
	    0314,0146,063,0231,0314,0146,063,0231
	}, {
	    0314,0146,063,0231,0314,0146,063,0231,
	    0314,0146,063,0231,0314,0146,063,0231
	}
    }, {					/* mix \\\ (dotted-dashed) */
	{
	    0376,0177,077,037,0217,0307,0343,0161,
	    070,034,0216,0307,0343,0361,0370,0374
	}, {
	    070,034,0216,0307,0343,0361,0370,0374,
	    0376,0177,077,037,0217,0307,0343,0161
	}
    }, {					/* thick /// (dashed) */
	{
	    0377,0377,0377,0377,0376,0374,0370,0360,
	    0341,0303,0207,017,037,077,0177,0377
	}, {
	    0341,0303,0207,017,037,077,0177,0377,
	    0377,0377,0377,0377,0376,0374,0370,0360
	}
    }, {					/* thin /// (dotted) */
	{
	    0231,063,0146,0314,0231,063,0146,0314,
	    0231,063,0146,0314,0231,063,0146,0314
	}, {
	    0231,063,0146,0314,0231,063,0146,0314,
	    0231,063,0146,0314,0231,063,0146,0314
	}
    }, {					/* mix /// (dotted-dashed) */
	{
	    0374,0370,0361,0343,0307,0216,034,070,
	    0161,0343,0307,0217,037,077,0177,0376
	}, {
	    0161,0343,0307,0217,037,077,0177,0376,
	    0374,0370,0361,0343,0307,0216,034,070
	}
    }, {					/* crisscross */
	{
	    0201,0102,044,030,030,044,0102,0201,
	    0201,0102,044,030,030,044,0102,0201
	}, {
	    0201,0102,044,030,030,044,0102,0201,
	    0201,0102,044,030,030,044,0102,0201
	}
    }, {					/* web */
	{
	    0377,0300,0240,0220,0210,0204,0202,0201,
	    0201,0202,0204,0210,0220,0240,0300,0200
	}, {
	    0377,02,04,010,020,040,0100,0200,
	    0200,0100,040,020,010,04,02,01
	}
    }
};

static int  context = 0;

struct setting {
	int  cnx;		/* setting context */
	char  *val;		/* attribute value */
	struct setting  *snext;	/* next setting in list */
};

static struct setting  *sets[0200];

static void spop(int attrib);
static int spat(int  pat, char  *patval);

void
set(		/* set attribute to value */
int  attrib,
char  *value
)
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



void
unset(			/* return attribute to previous setting */
int  attrib
)
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



void
reset(			/* return attribute to default setting */
int  attrib
)
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



static void
spop(			/* pop top off attrib settings list */
int  attrib
)
{

    if (sets[attrib] != NULL) {
	if (sets[attrib]->val != NULL)
	    freestr(sets[attrib]->val);
	free((char *)sets[attrib]);
	sets[attrib] = sets[attrib]->snext;
    }

}



static int
spat(			/* set a pattern */
int  pat,
char  *patval
)
{
    int  n, i, j, v;
    char  *nextscan();
    register char  *cp;

    if (patval == NULL) return(FALSE);

    if (patval[0] == 'P' || patval[0] == 'p') {
	if (nextscan(patval+1, "%d", &n) == NULL || n < 0 || n >= NPATS)
	    return(FALSE);
    } else {
	n = NPATS + pat - SPAT0;
	cp = patval;
	for (i = 0; i < PATSIZE>>3; i++)
	    for (j = 0; j < PATSIZE; j++) {
		if ((cp = nextscan(cp, "%o", &v)) == NULL || v < 0 || v > 0377)
		    return(FALSE);
		pattern[n][i][j] = v;
	    }
    }
    pati[pat-SPAT0] = n;
    return(TRUE);
}
