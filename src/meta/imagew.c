#ifndef lint
static const char	RCSid[] = "$Id: imagew.c,v 1.2 2003/08/01 14:14:24 schorsch Exp $";
#endif
/*
 *  Program to print meta-files on a dot-matrix printer
 *
 *  cc -o image image.c mplot.o plot.o mfio.o syscalls.o palloc.o misc.o
 *
 *  Apple Imagewriter
 */


#define  MAXALLOC  5000

#define  DXSIZE  576		/* x resolution */

#define  DYSIZE  576		/* y resolution */

#define  LINWIDT  DXSIZE	/* line width */

#define  LINHITE  8		/* line height */

#define  NLINES  (DYSIZE/LINHITE)

#define  CHARWIDTH  8

#define  PNORM  "\033n\017\033\""

#define  PCLEAR  "\033c"

#define  PINIT  "\033>\033n\033T16"

#define  LFNORM  "\033f"

#define  LFREV  "\033r"

#define  DBLON  "\033!"

#define  DBLOFF  ""

#define  PTAB  "\033F"

#define  OUTPUT  "\033g"

#define  XCOM  "pexpand +vtOCIps %s | psort -Y +x"




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

static char  chrtype[16][5] = {
				"\033N",
				"\033N\016",
				"\033N",
				"\033N\016",
				"\033E",
				"\033E\016",
				"\033E",
				"\033E\016",
				"\033q",
				"\033q\016",
				"\033q",
				"\033q\016",
				"\033Q",
				"\033Q\016",
				"\033Q",
				"\033Q\016"
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
    pinit();
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
    if (lineno > 0)
	nextpage();
    puninit();
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
       pinit();
       if ((fp = popen(command, "r")) == NULL)
          error(SYSTEM, "cannot execute input filter");
       plot(fp);
       pclose(fp);
       if (lineno > 0)
	   nextpage();
       puninit();
       }
    }

 return(0);
 }






thispage()		/* rewind and initialize current page */

{

    if (lineno != 0) {
        fputs(LFREV, stdout);
	while (lineno) {
	    putc('\n', stdout);
	    lineno--;
	}
	fputs(LFNORM, stdout);
    }

}




nextpage()		/* advance to next page */

{

    fputs("\f\r", stdout);

    lineno = 0;

}



contpage()		/* continue on this page */

{

    while (lineno++ < NLINES)
    	putc('\n', stdout);

    lineno = 0;

}



printspan()		/* output span to printer */

{
    register int  i;

    if (spanmin <= spanmax) {

	fprintf(stdout, "%s%4d", PTAB, spanmin);
	fprintf(stdout, "%s%3d", OUTPUT, (spanmax-spanmin)/8 + 1);

	for (i = spanmin; i <= spanmax; i++)
	    putc(flipbyte(outspan.cols[i]), stdout);

	i = 7 - (spanmax-spanmin)%8;
	while (i--)
	    putc('\0', stdout);

	putc('\r', stdout);
    }

    putc('\n', stdout);
    lineno++;

}




int
flipbyte(b)		/* flip an 8-bit byte end-to-end */

register int  b;

{
    register int  i, a = 0;

    if (b)
	for (i = 0; i < 8; i++) {
	    a <<= 1;
	    a |= b & 01;
	    b >>= 1;
	}

    return(a);
}




printstr(p)		/* output a string to the printer */

PRIMITIVE  *p;

{

    fprintf(stdout, "%s%4d", PTAB, CONV(p->xy[XMN], dxsize));

    if (p->arg0 & 0100)		/* double strike */
    	fputs(DBLON, stdout);
    else
    	fputs(DBLOFF, stdout);
    	
    fputs(chrtype[(p->arg0 >> 2) & 017], stdout);
    fputs(p->args, stdout);
    fputs(PNORM, stdout);
    putc('\r', stdout);

}



pinit()				/* initialize printer for output */
{
    pclear();
				/* get eighth bit */
    fputs("\033Z", stdout);
    putc('\0', stdout);
    putc(' ', stdout);
    fputs(PINIT, stdout);
}


puninit()			/* uninitialize printer */
{
    pclear();
}


pclear()			/* clear printer and sleep */

{
    register int  i, j;

    fputs(PCLEAR, stdout);
    fflush(stdout);
    for (i = 0; i < 1000; i++)
	for (j = 0; j < 500; j++)
	    ;

}
