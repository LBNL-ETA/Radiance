#ifndef lint
static const char	RCSid[] = "$Id: mx80.c,v 1.2 2003/08/01 14:14:24 schorsch Exp $";
#endif
/*
 *  Program to print meta-files on a dot-matrix printer
 *
 *  cc -o mx80 mx80.c mplot.o plot.o mfio.o syscalls.o palloc.o misc.o
 *
 *  Epson MX-80 printer
 */


#define  MAXALLOC  5000

#define  DXSIZE  480		/* x resolution */

#define  DYSIZE  616		/* y resolution */

#define  LINWIDT  DXSIZE	/* line width */

#define  LINHITE  8		/* line height */

#define  NLINES  (DYSIZE/LINHITE)	/* number of lines */

#define  CHARWIDTH  6

#define  PNORM  ""

#define  PINIT  "\033@\033#\0331\033U1"

#define  PUNINIT  "\033@\033="

#define  OUTPUT  "\033K"

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
    register int  k;

    if (spanmin <= spanmax) {

	k = spanmin/charwidth;
	while (k--)
	    putc(' ', stdout);

	k = spanmin%charwidth;
	fputs(OUTPUT, stdout);
	putc((spanmax-spanmin+k+1)%256, stdout);
	putc((spanmax-spanmin+k+1)/256, stdout);
	while (k--)
	    putc('\0', stdout);

	for (k = spanmin; k <= spanmax; k++)
	    putc(outspan.cols[k], stdout);

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

    fputs(p->args, stdout);
    fputs(PNORM, stdout);
    putc('\r', stdout);

}
