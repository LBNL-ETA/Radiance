#ifndef lint
static const char	RCSid[] = "$Id: mplot.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *   Plotting routines for meta-files to line-at-a-time printers
 */


#include  <fcntl.h>

#include  "meta.h"

#include  "plot.h"

#include  "span.h"



int  minwidth = 0;

static PLIST  inqueue = {NULL, NULL};

static PRIMITIVE  nextp;



plot(infp)		/* plot meta-file */

FILE  *infp;

{

    do {
	readp(&nextp, infp);
	initplot();
	while (inqueue.ptop != NULL || isprim(nextp.com))
	    plotspan(infp);
	doglobal(&nextp);
	fargs(&nextp);
    } while (nextp.com != PEOF);

}





initplot()			/* initialize this plot */

{

    thispage();
    outspan.xleft = 0;
    outspan.xright = dxsize - 1;
    outspan.ytop = dysize + linhite - 1;
    outspan.ybot = dysize;

}





doglobal(g)			/* execute a global command */

PRIMITIVE  *g;

{
    char  c;
    int  tty;

    switch (g->com) {

	case PEOF:
	    break;

	case PDRAW:
	    fflush(stdout);
	    break;

	case PEOP:
	    if (g->arg0 & 0200)		/* advance page */
	    	nextpage();
	    else if (g->arg0 == 3)	/* continue down */
	        contpage();
	    else
	        error(USER, "illegal continue direction in doglobal");
	    break;

	case PPAUSE:
	    fflush(stdout);
	    tty = open(TTY, O_RDWR);
	    if (g->args != NULL) {
		write(tty, g->args, strlen(g->args));
		write(tty, " - (hit return to continue)", 27);
	    } else
		write(tty, "\007", 1);
	    do {
	        c = '\n';
	        read(tty, &c, 1);
	    } while (c != '\n');
	    close(tty);
	    break;

	case PSET:
	    set(g->arg0, g->args);
	    break;

	case PUNSET:
	    unset(g->arg0);
	    break;

	case PRESET:
	    reset(g->arg0);
	    break;

	default:
	    sprintf(errmsg, "unknown command '%c' in doglobal", g->com);
	    error(WARNING, errmsg);
	    break;
	}

}





plotspan(infp)			/* plot next span */

FILE  *infp;

{
    PLIST  lastinq;
    register PRIMITIVE  *p;

					/* clear span */
    nextspan();
    					/* plot from queue */
    lastinq.ptop = inqueue.ptop;
    lastinq.pbot = inqueue.pbot;
    inqueue.ptop = inqueue.pbot = NULL;
    while ((p = pop(&lastinq)) != NULL) {
	doprim(p);
	pfree(p);
    }
					/* plot from file */
    while (isprim(nextp.com) && CONV(nextp.xy[YMX],dysize) >= outspan.ybot) {
        doprim(&nextp);
	fargs(&nextp);
	readp(&nextp, infp);
    }
    					/* print out span */
    outputspan();

}





nextspan()		/* prepare next span */

{
    register int  i;
    register char  *colp, *tcolp;
	
    if (spanmin <= spanmax) {			/* clear span */

	i = nrows*dxsize;
	colp = outspan.cols;
	tcolp = outspan.tcols;
	while (i--)
	    *colp++ = *tcolp++ = '\0';
    }
            
    outspan.ytop -= linhite;			/* advance to next */
    outspan.ybot -= linhite;
    spanmin = dxsize;
    spanmax = 0;

}



outputspan()		/* output span to printer */
{
    register int  i;
    register char  *colp, *tcolp;

    if (spanmin <= spanmax) {			/* overlay spans */

	i = nrows*dxsize;
	colp = outspan.cols;
	tcolp = outspan.tcols;
	while (i--)
	    *colp++ |= *tcolp++;
    }
    printspan();			/* print span */
}



doprim(p)		/* plot primitive */

register PRIMITIVE  *p;

{
    register PRIMITIVE  *newp;
    
    switch (p->com) {

	case PLSEG:
	    plotlseg(p);
	    break;

	case PRFILL:
	    fill((p->arg0&0103) | (pati[(p->arg0>>2)&03]<<2),
			CONV(p->xy[XMN],dxsize),CONV(p->xy[YMN],dysize),
			CONV(p->xy[XMX],dxsize)+(p->arg0&0100?-1:0),
			CONV(p->xy[YMX],dysize)+(p->arg0&0100?-1:0));
	    break;

	case PTFILL:
	    tfill(p);
	    break;
		
	case PMSTR:
	    printstr(p);
	    break;

	default:
	    sprintf(errmsg, "unknown command '%c' in doprim", p->com);
	    error(WARNING, errmsg);
	    return;
    }

    if (CONV(p->xy[YMN],dysize) < outspan.ybot) {	/* save for next time */
        if ((newp = palloc()) == NULL)
            error(SYSTEM, "memory limit exceeded in doprim");
    	mcopy((char *)newp, (char *)p, sizeof(PRIMITIVE));
	newp->args = savestr(p->args);
        add(newp, &inqueue);
    }
        
}





plotlseg(p)		/* plot a line segment */

register PRIMITIVE  *p;

{
    register int  ti;
    long  run2 = 0L, rise2 = 0L;
    int  x, y, run, rise, xstop, ystop, hrad, vrad, lpat, n;

						/* compute line pattern */
    lpat = (p->arg0 >> 2) & 014;
    if (p->arg0 & 0100 && lpat != 0)
	lpat += 014;
    lpat |= p->arg0 & 03;

    ti = (p->arg0 >> 2) & 03;			/* compute line radius */
    ti = WIDTH(ti) / 2;
    hrad = CONV(ti, dxsize);
    vrad = CONV(ti, dysize);
    if (hrad < minwidth)
	hrad = minwidth;
    if (vrad < minwidth)
	vrad = minwidth;

    x = CONV(p->xy[XMX], dxsize);		/* start at top */
    y = CONV(p->xy[YMX], dysize);
    run = CONV(p->xy[XMN], dxsize) - x;
    rise = CONV(p->xy[YMN], dysize) - y;

    if (p->arg0 & 0100)				/* slope < 0; reverse x */
	x -= (run = -run);

    xstop = x + run;				/* compute end point */
    ystop = y + rise;

    if ((ti = outspan.ytop+vrad+1-y) < 0) {	/* adjust to top of span */
	run2 = rise2 = (long)ti*run;
	x += rise2/rise;
	y += ti;
    }

    if ((ti = outspan.ybot-vrad-1-ystop) > 0) {	/* adjust to bottom of span */
        xstop += (long)ti*run/rise;
	ystop += ti;
    }

    if (abs(run) > -rise)
	n = abs(xstop - x);
    else
	n = y - ystop;
    
    paintline(x, y, run, rise, hrad, vrad, lpat, run2, rise2, n);

}



/*
 *  This routine paints a line with calls to fill().  The line can
 *    start and end at arbitrary points on a longer line segment.
 */

paintline(x, y, run, rise, hrad, vrad, lpat, run2, rise2, n)

register int  x, y;
int  run, rise;
int  hrad, vrad;
int  lpat;
long  run2, rise2;
int  n;

{
    int  xstep, ystep;

    if (run >= 0)
	xstep = 1;
    else {
	xstep = -1;
	run = -run;
    }
    if (rise >= 0)
	ystep = 1;
    else {
	ystep = -1;
	rise = -rise;
    }
    if (n < 0)
	n = max(run, rise);

    if (run > rise)
	while (n >= 0)
	    if (run2 >= rise2) {
		fill(lpat, x-hrad, y-vrad, x+hrad, y+vrad);
		n--;
		x += xstep;
		rise2 += rise;
	    } else {
		y += ystep;
		run2 += run;
	    }
    else
	while (n >= 0)
	    if (rise2 >= run2) {
		fill(lpat, x-hrad, y-vrad, x+hrad, y+vrad);
		n--;
		y += ystep;
		run2 += run;
	    } else {
		x += xstep;
		rise2 += rise;
	    }

}



tfill(p)			/* fill a triangle */
register PRIMITIVE  *p;
{
    register int  x, txmin, txmax;	
    int  xmn, ymn, tpat;
    long  xsz, ysz;

    xmn = CONV(p->xy[XMN], dxsize);
    xsz = CONV(p->xy[XMX], dxsize) - xmn;
    ymn = CONV(p->xy[YMN], dysize);
    ysz = CONV(p->xy[YMX], dysize) - ymn;
    if (xsz <= 0 || ysz <= 0)
	return;	
    txmin = (outspan.ybot - ymn)*xsz/ysz;
    txmax = (outspan.ytop - ymn)*xsz/ysz;
    if (p->arg0 & 020) {				/* up or down */
	x = txmin;
	txmin = xsz - txmax;
	txmax = xsz - x;
    }
    txmin += xmn;
    txmax += xmn;
    txmin = max(txmin, xmn);
    txmax = min(txmax, xmn + (int)xsz - 1);
    tpat = (p->arg0&0103) | (pati[(p->arg0>>2)&03]<<2);

    if (p->arg0 & 040) {				/* left or down */
	fill(tpat, xmn, ymn, txmin - 1, ymn + (int)ysz - 1);
	for (x = txmin; x <= txmax; x++)
	    if (p->arg0 & 020)			/* down */
		fill(tpat, x, ymn, x, (int)(ysz-(x-xmn)*ysz/xsz) + ymn - 1);
	    else				/* left */
		fill(tpat, x, (int)((x-xmn)*ysz/xsz) + ymn, x, ymn + (int)ysz - 1);
    } else {						/* right or up */
	for (x = txmin; x <= txmax; x++)
	    if (p->arg0 & 020)			/* up */
		fill(tpat, x, (int)(ysz-(x-xmn)*ysz/xsz) + ymn, x, ymn + (int)ysz - 1);
	    else				/* right */
		fill(tpat, x, ymn, x, (int)((x-xmn)*ysz/xsz) + ymn - 1);
	fill(tpat, txmax + 1, ymn, xmn + (int)xsz - 1, ymn + (int)ysz - 1);
    }
}




fill(attrib, xmin, ymin, xmax, ymax)	/* fill rectangle with attribute */

int  attrib;
int  xmin, ymin, xmax, ymax;

{
    int  filpat;
    int  rpos;
    unsigned char  *pattr;
    register char  *colp;
    register int  i;

    xmin = max(xmin, outspan.xleft) - outspan.xleft;
    ymin = max(ymin, outspan.ybot) - outspan.ybot;
    xmax = min(xmax, outspan.xright) - outspan.xleft;
    ymax = min(ymax, outspan.ytop) - outspan.ybot;

    for (rpos = 0; rpos < nrows; rpos++)

	if (rpos >= ymin >> 3 && rpos <= ymax >> 3) {

	    filpat = 0377;
	    if (rpos == ymin >> 3) {
		i = ymin & 07;
		filpat = (filpat >> i) << i;
	    }
	    if (rpos == ymax >> 3) {
		i = ~ymax & 07;
		filpat = ((filpat << i) & 0377) >> i;
	    }

	    pattr = pattern[(attrib&074)>>2]
				[((outspan.ybot>>3)+rpos)%(PATSIZE>>3)];

	    if (attrib & 0100) {
		colp = &outspan.tcols[rpos*dxsize + xmin];
		for (i = xmin; i <= xmax; i++)
		    *colp++ ^= filpat & pattr[i%PATSIZE];
	    } else {
		colp = &outspan.cols[rpos*dxsize + xmin];
		for (i = xmin; i <= xmax; i++)
		    *colp++ |= filpat & pattr[i%PATSIZE];
	    }

	    spanmin = min(xmin, spanmin);
	    spanmax = max(xmax, spanmax);
	}

}
