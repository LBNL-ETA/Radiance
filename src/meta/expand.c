#ifndef lint
static const char	RCSid[] = "$Id: expand.c,v 1.4 2003/11/15 02:13:36 schorsch Exp $";
#endif
/*
 *   Expansion routine for command implementation
 *   on simple drivers
 */

#include  <string.h>

#include  "meta.h"

				/* vector characters file */
#define  VINPUT  "vchars.mta"

#define  FIRSTCHAR  ' '			/* first vector character */

#define  LASTCHAR  '~'			/* last vector character */

#define  PEPS  5			/* minimum size for fill area */

#define  fillok(p)  ((p)->xy[XMX] - (p)->xy[XMN] > PEPS && \
			(p)->xy[YMX] - (p)->xy[YMN] > PEPS)



static short  *xlist;			/* pointer to expansion list */

static void exfile(FILE  *fp);
static void exprim(PRIMITIVE  *p);
static void sendmstr(register PRIMITIVE  *p);
static void sendvstr(register PRIMITIVE	*p);
static void include(int  code, char  *fname);
static void polyfill(PRIMITIVE  *p);
static void polyedge(PRIMITIVE  *p, int  x1, int y1, int x2, int y2);


void
expand(		/* expand requested commands */
FILE  *infp,
short  *exlist
)
{
 static PRIMITIVE  pincl = {PINCL, 02, {-1, -1, -1, -1}, VINPUT, NULL};
 static short  vcloaded = FALSE;

 xlist = exlist;

 if (exlist[comndx(PVSTR)] == 1 && !vcloaded) {
    vcloaded = TRUE;
    exprim(&pincl);
    }

 exfile(infp);

 if (inseg())
    error(USER, "unclosed segment in expand");
    
 }




static void
exfile(			/* expand the given file */
register FILE  *fp
)
{
    PRIMITIVE  curp;

    while (readp(&curp, fp)) {
	exprim(&curp);
	fargs(&curp);
    }

}




static void
exprim(			/* expand primitive */
register PRIMITIVE  *p
)
{
    int  xflag = xlist[comndx(p->com)];
				/* 1==expand, 0==pass, -1==discard */
    if (xflag != -1) {
       if (xflag == 1)
          switch (p->com)  {

	     case POPEN:
	        openseg(p->args);
		break;

	     case PSEG:
		segment(p, exprim);
		break;

	     case PCLOSE:
		closeseg();
		break;
	        
	     case PINCL:
		include(p->arg0, p->args);
	        break;

	     case PMSTR:
		sendmstr(p);
		break;
	        
	     case PVSTR:
		sendvstr(p);
		break;

	     case PPFILL:
	        polyfill(p);
	        break;
	        
	     default:
		sprintf(errmsg, "unsupported command '%c' in exprim", p->com);
		error(WARNING, errmsg);
	        break;
	     }
       else if (inseg())
          segprim(p);
       else
          writep(p, stdout);
    }
}




static void
sendmstr(			/* expand a matrix string */
register PRIMITIVE  *p
)
{
    PRIMITIVE  pout;
    int  cheight, cwidth, cthick, ccol;

    cheight = 350;
    if (p->arg0 & 010)
	cheight *= 2;
    cwidth = (6 - ((p->arg0 >> 4) & 03)) * 35;
    if (p->arg0 & 04)
	cwidth *= 2;
    cthick = (p->arg0 & 0100) ? 1 : 0;
    ccol = p->arg0 & 03;

    pout.com = PVSTR;
    pout.arg0 = (cthick << 2) | ccol;
    pout.xy[XMN] = p->xy[XMN];
    pout.xy[YMN] = p->xy[YMX] - cheight/2;
    pout.xy[XMX] = p->xy[XMN] + strlen(p->args)*cwidth;
    if (pout.xy[XMX] >= XYSIZE)
	pout.xy[XMX] = XYSIZE-1;
    pout.xy[YMX] = p->xy[YMX] + cheight/2;
    pout.args = p->args;

    exprim(&pout);

}



static void
sendvstr(			/* expand a vector string */
register PRIMITIVE	*p
)
{
    PRIMITIVE  pout;
    int  xadv = 0, yadv = 0;
    char  s[3];
    register char  *cp;

    if (p->args == NULL)
	error(USER, "illegal empty string in sendvstr");

    pout.com = PSEG;
    pout.arg0 = p->arg0;
    switch (p->arg0 & 060) {
	case 0:			/* right */
	    xadv = (p->xy[XMX] - p->xy[XMN])/strlen(p->args);
	    pout.xy[XMN] = p->xy[XMN];
	    pout.xy[XMX] = p->xy[XMN] + xadv;
	    yadv = 0;
	    pout.xy[YMN] = p->xy[YMN];
	    pout.xy[YMX] = p->xy[YMX];
	    break;
	case 020:		/* up */
	    xadv = 0;
	    pout.xy[XMN] = p->xy[XMN];
	    pout.xy[XMX] = p->xy[XMX];
	    yadv = (p->xy[YMX] - p->xy[YMN])/strlen(p->args);
	    pout.xy[YMN] = p->xy[YMN];
	    pout.xy[YMX] = p->xy[YMN] + yadv;
	    break;
	case 040:		/* left */
	    xadv = -(p->xy[XMX] - p->xy[XMN])/strlen(p->args);
	    pout.xy[XMN] = p->xy[XMX] + xadv;
	    pout.xy[XMX] = p->xy[XMX];
	    yadv = 0;
	    pout.xy[YMN] = p->xy[YMN];
	    pout.xy[YMX] = p->xy[YMX];
	    break;
	case 060:		/* down */
	    xadv = 0;
	    pout.xy[XMN] = p->xy[XMN];
	    pout.xy[XMX] = p->xy[XMX];
	    yadv = -(p->xy[YMX] - p->xy[YMN])/strlen(p->args);
	    pout.xy[YMN] = p->xy[YMX] + yadv;
	    pout.xy[YMX] = p->xy[YMX];
	    break;
    }

    pout.args = s;
    s[1] = '\'';
    s[2] = '\0';
    for (cp = p->args; *cp; cp++)
	if (*cp < FIRSTCHAR || *cp > LASTCHAR) {
	    sprintf(errmsg, "unknown character (%d) in sendvstr", *cp);
	    error(WARNING, errmsg);
	}
	else {
	    s[0] = *cp;
	    exprim(&pout);
	    pout.xy[XMN] += xadv;
	    pout.xy[XMX] += xadv;
	    pout.xy[YMN] += yadv;
	    pout.xy[YMX] += yadv;
	}

}



static void
include(			/* load an include file */
int  code,
char  *fname
)
{
    register FILE  *fp = NULL;

    if (fname == NULL)
	error(USER, "missing include file name in include");
    
    if (code == 2 || (fp = fopen(fname, "r")) == NULL) {
	if (code != 0)
	    fp = mfopen(fname, "r");
	else {
	    sprintf(errmsg, "cannot open user include file \"%s\"", fname);
	    error(USER, errmsg);
	}
    }
    
    exfile(fp);
    fclose(fp);

}



static void
polyfill(			/* expand polygon fill command */
PRIMITIVE  *p
)
{
    char  firstx, firsty, x, y;
	int lastx, lasty;
    register char  *cp;
    
    if ((cp=nextscan(nextscan(p->args,"%d",&firstx),"%d",&firsty)) == NULL) {
        sprintf(errmsg, "illegal polygon spec \"%s\" in polyfill", p->args);
        error(WARNING, errmsg);
        return;
    }
    
    lastx = firstx;
    lasty = firsty;
    
    while ((cp=nextscan(nextscan(cp,"%d",&x),"%d",&y)) != NULL) {
    
        polyedge(p, lastx, lasty, x, y);
        lastx = x;
        lasty = y;
    }
    
    polyedge(p, lastx, lasty, firstx, firsty);
}



static void
polyedge(		/* expand edge of polygon */
PRIMITIVE  *p,
int  x1, int y1, int x2, int y2
)
{
    int  reverse;
    PRIMITIVE  pin, pout;
    
    if (x1 < x2) {
        pin.xy[XMN] = x1;
        pin.xy[XMX] = x2;
        reverse = FALSE;
    } else {
        pin.xy[XMN] = x2;
        pin.xy[XMX] = x1;
        reverse = TRUE;
    }
    
    if (y1 < y2) {
        pin.xy[YMN] = y1;
        pin.xy[YMX] = y2;
    } else {
        pin.xy[YMN] = y2;
        pin.xy[YMX] = y1;
        reverse = y1 > y2 && !reverse;
    }
    
    pout.xy[XMN] = xlate(XMN, &pin, p);
    pout.xy[XMX] = xlate(XMX, &pin, p);
    pout.xy[YMN] = xlate(YMN, &pin, p);
    pout.xy[YMX] = xlate(YMX, &pin, p);
    pout.com = PTFILL;
    pout.arg0 = 0100 | (reverse << 4) | (p->arg0 & 017);
    pout.args = NULL;
    if (fillok(&pout))
	exprim(&pout);
    
    if (p->arg0 & 0100) {		/* draw border */
        pout.com = PLSEG;
        pout.arg0 = reverse << 6;
        exprim(&pout);
    }
    
    pout.com = PRFILL;
    pout.arg0 = 0100 | (p->arg0 & 017);
    pout.xy[XMN] = pout.xy[XMX];
    pout.xy[XMX] = p->xy[XMX];
    if (fillok(&pout))
	exprim(&pout);
    
}
