#ifndef lint
static const char	RCSid[] = "$Id: tbar.c,v 1.4 2003/08/01 14:14:24 schorsch Exp $";
#endif
/*
 *       PROGRAM TO PLOT TEL-A-GRAF POINTS TO METAFILE
 *
 *              Greg Ward
 *              12/12/84
 *
 *      cc -o ../tbar tbar.c tgraph.o primout.o mfio.o syscalls.o misc.o -lm
 */

#include "rtprocess.h" /* getpid() */
#include "tgraph.h"


#define  XLEGEND  (XBEG+XSIZ+4*TSIZ)	/* x start of legend */

#define  YLEGEND  (YBEG+2*YSIZ/3)	/* y start of legend */

#define  MAXBAR  1000			/* maximum bar width */

short  usecurve[NCUR];			/* booleans for curve usage */

double  xmin, ymin, xmax, ymax;		/* domain */

double  xsize, ysize;			/* axis dimensions */

double  xmnset = -FHUGE, xmxset = FHUGE,	/* domain settings */
        ymnset = -FHUGE, ymxset = FHUGE;

short  logx = FALSE, logy = FALSE;	/* flags for log plots */

short  polar = FALSE;			/* flag for polar plots */

short  grid = FALSE;			/* flag for grid */

int    curtype[NCUR] = {0, 04, 010, 014, 01, 05, 011, 015,
			 02, 06, 012, 016, 03, 07, 013, 017};

int    ncurves;

int    xlegend,
       ylegend;				/* current legend position */

int    symrad = SYMRAD;			/* symbol radius */ 

char  *progname;





main(argc, argv)

int  argc;
char  **argv;

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

 polar = FALSE;			/* override stupid choices */
 logx = FALSE;
 axflag = BOX|ORIGIN|YTICS|YNUMS;
 if (grid)
    axflag |= YGRID;

 if (argc)

    for ( ; argc--; argv++)  {

       fp = efopen(*argv, "r");
       normalize(fp, NULL);
       if (ymin > 0)
          ymin = 0;
       else if (ymax < 0)
	  ymax = 0;
       xmax++;
       makeaxis(axflag);
       fseek(fp, 0L, 0);
       plot(fp);
       fclose(fp);
       }
 else  {

    sprintf(tfname, "%stb%d", TDIR, getpid());
    fp = efopen(tfname, "w+");
    normalize(stdin, fp);
    if (ymin > 0)
       ymin = 0;
    else if (ymax < 0)
       ymax = 0;
    xmax++;
    makeaxis(axflag);
    fseek(fp, 0L, 0);
    plot(fp);
    fclose(fp);
    unlink(tfname);
    }

 pglob(PEOF, 0200, NULL);

 return(0);
 }





plot(fp)			/* read file and generate plot */

FILE  *fp;

{
 int  ncur = 0; 		/* curves seen so far */
 char  line[255], *s;
 double  x, y;

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

    else if (isdivlab(line)) {
       s = line;
       x = 1;
       while ((s = snagquo(s)) != NULL) {
          if (x >= xmin && x <= xmax)
	     boxstring(0, XCONV(x), YBEG-750, XCONV(x+0.66), YBEG-400, s);
	  s += strlen(s)+1;
	  x++;
	  }
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
          boxout(curtype[ncur], xlegend-200, ylegend, xlegend+200, ylegend+250);
	  pprim(PMSTR, 020, xlegend+400, ylegend+200,
		      xlegend+400, ylegend+200, snagquo(line));
	  ylegend -= 500;
	  }
       }
       
    else if (usecurve[ncur] && isdata(line)) {

       if (getdata(line, &x, &y) >= 0)
          barout(ncur, x, y);

       }

 pglob(PEOP, 0200, NULL);

 }



barout(cn, x, y)		/* output bar for curve cn, value (x,y) */

int  cn;
double  x, y;

{
    int  barleft, barwidth;
    int  barlower, barheight;
    
    barwidth = XSIZ/xsize/(ncurves+1)*0.66;
    if (barwidth > MAXBAR)
	barwidth = MAXBAR;

    barleft = XCONV(x) + cn*barwidth;
    
    if (y < 0.0) {
	barlower = YCONV(y);
	barheight = YCONV(0.0) - barlower;
    } else {
	barlower = YCONV(0.0);
	barheight = YCONV(y) - barlower;
    }
    boxout(curtype[cn], barleft, barlower,
		barleft+barwidth, barlower+barheight);
    
}



boxout(a0, xmn, ymn, xmx, ymx)			/* output a box */

int  a0;
int  xmn, ymn, xmx, ymx;

{
    
    pprim(PRFILL, a0, xmn, ymn, xmx, ymx, NULL);
    plseg(0, xmn, ymn, xmx, ymn);
    plseg(0, xmn, ymn, xmn, ymx);
    plseg(0, xmn, ymx, xmx, ymx);
    plseg(0, xmx, ymn, xmx, ymx);
    
}
