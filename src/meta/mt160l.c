#ifndef lint
static const char	RCSid[] = "$Id: mt160l.c,v 1.2 2003/08/01 14:14:24 schorsch Exp $";
#endif
/*
 *  Program to print meta-files on a dot-matrix printer
 *
 *  cc -o mt160l mt160l.c mplot.o plot.o mfio.o syscalls.o palloc.o misc.o
 *
 *  Mannesman Tally MT160L high-density
 */


#define  MAXALLOC  5000

#define  DXSIZE  1064		/* x resolution */

#define  DYSIZE  1024		/* y resolution */

#define  LINWIDT  DXSIZE	/* line width */

#define  LINHITE  16		/* line height */

#define  NLINES  (DYSIZE/LINHITE)	/* number of lines */

#define  CHARWIDTH  8

#define  PNORM  "\033[0y\033[6w"

#define  PINIT  "\033[6~\033[7z\033[0y\033[6w"

#define  PUNINIT  "\033[6~"

#define  DBLON  "\033[=z"

#define  DBLOFF  "\033[>z"

#define  OUTHIGH  "\033[9~\033%6"

#define  OUTLOW  "\033[8~\033%6"

#define  XCOM  "pexpand +vOCIsp %s | psort -Y +x"




#include  "meta.h"

#include  "plot.h"

#include  "span.h"




char  *progname;

struct span  outspan;

int  dxsize = DXSIZE, dysize = DYSIZE,
     linwidt = LINWIDT, linhite = LINHITE,
     nrows = (LINHITE-1)/8+1;

int  maxalloc = MAXALLOC;

int  spanmin = 0, spanmax = LINWIDT-1;

int  charwidth = CHARWIDTH;

static char  chrtype[16][9] = {
				"\033[1y\033[4y",
				"\033[0y\033[0w",
				"\033[1y\033[4y",
				"\033[0y\033[0w",
				"\033[1y\033[5y",
				"\033[0y\033[1w",
				"\033[1y\033[5y",
				"\033[0y\033[1w",
				"\033[0y\033[6w",
				"\033[0y\033[2w",
				"\033[0y\033[6w",
				"\033[0y\033[2w",
				"\033[0y\033[7w",
				"\033[0y\033[3w",
				"\033[0y\033[7w",
				"\033[0y\033[3w"
				};

static int  lineno = 0;

static short  condonly = FALSE,
	      conditioned = FALSE;


main(argc, argv)

int  argc;
char  **argv;

{
 FILE  *fp;
 FILE  *popen();
 char  comargs[200], command[300];

 progname = *argv++;
 argc--;
 
 condonly = FALSE;
 conditioned = FALSE;

 while (argc && **argv == '-')  {
    switch (*(*argv+1))  {
       case 'c':
	  condonly = TRUE;
	  break;
       case 'r':
	  conditioned = TRUE;
	  break;
       default:
	  error(WARNING, "unknown option");
	  break;
       }
    argv++;
    argc--;
    }

 if (conditioned) {
    fputs(PINIT, stdout);
    if (argc)
       while (argc)  {
	  fp = efopen(*argv, "r");
	  plot(fp);
	  fclose(fp);
	  argv++;
	  argc--;
	  }
    else
       plot(stdin);
    if (lineno)
       nextpage();
    fputs(PUNINIT, stdout);
 } else  {
    comargs[0] = '\0';
    while (argc)  {
       strcat(comargs, " ");
       strcat(comargs, *argv);
       argv++;
       argc--;
       }
    sprintf(command, XCOM, comargs);
    if (condonly)
       return(system(command));
    else  {
       fputs(PINIT, stdout);
       if ((fp = popen(command, "r")) == NULL)
          error(SYSTEM, "cannot execute input filter");
       plot(fp);
       pclose(fp);
       if (lineno)
	  nextpage();
       fputs(PUNINIT, stdout);
       }
    }

 return(0);
 }







thispage()		/* rewind and initialize current page */

{

    if (lineno)
	error(USER, "cannot restart page in thispage");

}




nextpage()		/* advance to next page */

{

    fputs("\f\r", stdout);

    lineno = 0;

}




contpage()		/* continue new plot on current page */

{

    while (lineno++ < NLINES)
        putc('\n', stdout);
    
    lineno = 0;
    
}



printspan()		/* output span to printer */

{
    register unsigned shiftreg;
    unsigned short  outc;
    register  k;
    unsigned short offset, mask;
    int  i,j;

    if (spanmin <= spanmax)
	for (offset = 0; offset < 2; offset++) {

	    k = spanmin/charwidth;
	    while (k--)
		putc(' ', stdout);

	    k = spanmin%charwidth;
	    if (offset)
		fputs(OUTHIGH, stdout);
	    else
		fputs(OUTLOW, stdout);
	    putc((spanmax-spanmin+k+1)%256, stdout);
	    putc((spanmax-spanmin+k+1)/256, stdout);

	    while (k--)
		putc('\0', stdout);

	    for (i = spanmin; i <= spanmax; i++) {
		outc = '\0';
		for (j = 0; j < 2; j++) {
		    shiftreg = outspan.cols[j*dxsize + i] >> offset;
		    if (j) {
			shiftreg <<= 4;
			mask = 1 << 4;
		    } else
			mask = 1;
		    for (k = 0; k < 4; k++) {
			outc |= shiftreg & (mask << k);
			shiftreg >>= 1;
		    }
		}
		putc(outc, stdout);
	    }
	    putc('\r', stdout);
	}
    putc('\n', stdout);
    lineno++;

}





printstr(p)		/* output a string to the printer */

PRIMITIVE  *p;

{
    int  i;

    i = CONV(p->xy[XMN], dxsize)/charwidth;
    while (i--)
	putc(' ', stdout);

    if (p->arg0 & 0100)		/* double strike */
    	fputs(DBLON, stdout);
    else
    	fputs(DBLOFF, stdout);
    	
    fputs(chrtype[(p->arg0 >> 2) & 017], stdout);
    fputs(p->args, stdout);
    fputs(PNORM, stdout);
    putc('\r', stdout);

}
