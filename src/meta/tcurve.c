#ifndef lint
static const char	RCSid[] = "$Id: tcurve.c,v 1.5 2003/11/15 02:13:37 schorsch Exp $";
#endif
/*
 *       PROGRAM TO PLOT TEL-A-GRAF CURVES TO METAFILE
 *
 *              Greg Ward
 *              12/12/84
 *
 *      cc -o ../tcurve tcurve.c tgraph.o primout.o mfio.o syscalls.o misc.o -lm
 */

#include "rtprocess.h" /* getpid() */
#include  "tgraph.h"
#include  "plot.h"

#define  XLEGEND  (XBEG+XSIZ+4*TSIZ)	/* x start of legend */

#define  YLEGEND  (YBEG+2*YSIZ/3)	/* y start of legend */

short  usecurve[NCUR];		/* booleans for curve usage */

double  xmin, ymin, xmax, ymax;		/* domain */

double  xsize, ysize;			/* axis dimensions */

double  xmnset = -FHUGE, xmxset = FHUGE,	/* domain settings */
        ymnset = -FHUGE, ymxset = FHUGE;

short  logx = FALSE, logy = FALSE;	/* flags for log plots */

short  polar = FALSE;			/* flag for polar plot */
       
short  grid = FALSE;			/* flag for grid */

int    curtype[NCUR] = {0, 020, 040, 060, 01, 021, 041, 061,
			 02, 022, 042, 062, 03, 023, 043, 063};

char   *sym[NCUR] = {"ex", "triangle", "square", "triangle2", "diamond",
			"cross", "octagon", "crosssquare", "exsquare",
			"trianglesquare", "triangle2square", "crossdiamond",
			"crossoctagon", "exoctagon", "block", "bullet"};
		       
int    ncurves;

int    xlegend,
       ylegend;				/* current legend position */

int    symrad = SYMRAD;			/* symbol radius */ 

char  *progname;


static void limline(int a0, double x,double y,double xout,double yout, char *s);


int
main(
	int  argc,
	char  **argv
)
/*
 *     Take Tel-A-Graf runnable files and convert them to
 *  metafile primitives to send to standard output
 */

{
 char  tfname[MAXFNAME];
 FILE  *fp;
 int  axflag;

 progname = *argv++;
 argc--;

 initialize();

 for (; argc && (**argv == '-' || **argv == '+'); argc--, argv++)
    option(*argv);

 if (polar)			/* avoid dumb choices */
    logx = FALSE;

 axflag = XTICS|XNUMS|YTICS|YNUMS;
 if (grid)
    axflag |= XGRID|YGRID;
 if (polar)
    axflag |= ORIGIN;
 else
    axflag |= BOX;

 pglob(PINCL, 2, "symbols.mta");

 if (argc)

    for ( ; argc--; argv++)  {

       fp = efopen(*argv, "r");
       normalize(fp, NULL);
       makeaxis(axflag);
       fseek(fp, 0L, 0);
       plot(fp);
       fclose(fp);
       }
 else  {

    sprintf(tfname, "%stc%d", TDIR, getpid());
    fp = efopen(tfname, "w+");
    normalize(stdin, fp);
    makeaxis(axflag);
    fseek(fp, 0L, 0);
    plot(fp);
    fclose(fp);
    unlink(tfname);
    }

 pglob(PEOF, 0200, NULL);

 return(0);
 }




void
plot(			/* read file and generate plot */
	FILE  *fp
)
{
 int  ncur = 0; 		/* curves seen so far */
 int  cur = 0;			/* current curve pattern */
 char  line[255], *s;
 double  x, y;
 double lastx = 0, lasty = 0;
 short  oobounds = FALSE, firstpoint = TRUE;

 xlegend = XLEGEND;
 ylegend = YLEGEND;

 if (ncurves > 0) {
    pprim(PMSTR, 0100, xlegend, ylegend+800, xlegend, ylegend+800, "Legend:");
    pprim(PMSTR, 0100, xlegend, ylegend+800, xlegend, ylegend+800, "______");
    }

 while (fgets(line, sizeof line, fp) != NULL)

    if (istitle(line)) {
       s = snagquo(line);
       boxstring(0, 0, YBEG+YSIZ+1000, XYSIZE-1, YBEG+YSIZ+1500, s);
       }

    else if (isxlabel(line)) {
       s = snagquo(line);
       boxstring(0, XBEG, YBEG-1250, XBEG+XSIZ, YBEG-900, s);
       }

    else if (isylabel(line)) {
       s = snagquo(line);
       boxstring(020, XBEG-1900, YBEG, XBEG-1550, YBEG+YSIZ, s);
       }

    else if (islabel(line)) {
       if (++ncur < NCUR && usecurve[ncur]) {
	  oobounds = FALSE;
	  firstpoint = TRUE;
	  cur++;
	  plseg(curtype[cur], xlegend, ylegend, xlegend+2000, ylegend);
	  symout(curtype[cur] & 03, xlegend, ylegend, sym[ncur]);
	  pprim(PMSTR, 020, xlegend+400, ylegend+200,
	 	      xlegend+400, ylegend+200, snagquo(line));
	  ylegend -= 500;
	  }
       }
       
    else if (isdata(line) && usecurve[ncur])  {

       if (getdata(line, &x, &y) >= 0) {

	  if (oobounds) {
	     oobounds = FALSE;
	     limline(curtype[cur], x, y, lastx, lasty, sym[ncur]);
	     }
	  else if (firstpoint) {
	     firstpoint = FALSE;
	     limline(curtype[cur], x, y, x, y, sym[ncur]);
	     }
	  else
	     plseg(curtype[cur], XCONV(lastx), YCONV(lasty),
				 XCONV(x), YCONV(y));

	  lastx = x;
	  lasty = y;

          }
       else {
	  if (!oobounds) {
	     oobounds = TRUE;
	     if (firstpoint)
		firstpoint = FALSE;
	     else
		limline(curtype[cur], lastx, lasty, x, y, NULL);
	     }
	  lastx = x;
	  lasty = y;
	  }

       }

 pglob(PEOP, 0200, NULL);

 }



void
limline(	/* print line from/to out of bounds */
	int  a0,
	double  x,
	double  y,
	double  xout,
	double  yout,
	char  *s
)
{

    for ( ; ; )
	if (xout < xmin) {
	    yout = (yout - y)*(xmin - x)/(xout - x) + y;
	    xout = xmin;
	} else if (yout < ymin) {
	    xout = (xout - x)*(ymin - y)/(yout - y) + x;
	    yout = ymin;
	} else if (xout > xmax) {
	    yout = (yout - y)*(xmax - x)/(xout - x) + y;
	    xout = xmax;
	} else if (yout > ymax) {
	    xout = (xout - x)*(ymax - y)/(yout - y) + x;
	    yout = ymax;
	} else {
	    if (s != NULL)
		symout(a0 & 03, XCONV(xout), YCONV(yout), s);
	    plseg(a0, XCONV(x), YCONV(y), XCONV(xout), YCONV(yout));
	    return;
	}        

}
