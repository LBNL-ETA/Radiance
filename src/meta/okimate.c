#ifndef lint
static const char	RCSid[] = "$Id: okimate.c,v 1.3 2003/10/27 10:28:59 schorsch Exp $";
#endif
/*
 *  Program to print meta-files on a dot-matrix printer
 *
 *  cc -o okimate okimate.c mplot.o plot.o mfio.o syscalls.o palloc.o misc.o
 *
 *  Okimate 20 printer (version for black ribbon)
 */


#define  MAXALLOC  5000

#define  DXSIZE  960		/* x resolution */

#define  DYSIZE  1152		/* y resolution */

#define  LINWIDT  DXSIZE	/* line width */

#define  LINHITE  24		/* line height */

#define  NLINES  (DYSIZE/LINHITE)	/* number of lines */

#define  CHARWIDTH  12

#define  PNORM  "\024\022\033%H"

#define  PINIT  "\033O\033T\033I\002\033%H\033T\022\033A\014\0332"

#define  PUNINIT  ""

#define  DBLON  "\033%G"

#define  GRPH	"\033%O"

#define  XCOM  "pexpand +vOCIsp %s | psort -Y +x"


#include  "rtprocess.h"
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

static char  chrtype[16][4] = {
		"",
		"\016",
		"",
		"\016",
		"\033:",
		"\033:\016",
		"\033:",
		"\033:\016",
		"\017",
		"\017\016",
		"\017",
		"\017\016",
		"\017",
		"",
		"\017",
		""
				};

static int  lineno = 0;

static short  condonly = FALSE,
	      conditioned = FALSE;


main(argc, argv)

int  argc;
char  **argv;

{
 FILE  *fp;
 char  comargs[200], command[300];

 progname = *argv++;
 argc--;
 
 condonly = FALSE;
 conditioned = FALSE;

 minwidth = 1;			/* so lines aren't invisible */

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

    fputs("\r\f", stdout);

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
    register int  i,j;

    if (spanmin <= spanmax) {
	    j = spanmin/charwidth;
	    while (j--)
		putc(' ', stdout);

	    fputs(GRPH, stdout);
	    j = spanmin%charwidth;
	    putc((spanmax-spanmin+j+1)%256, stdout);
	    putc((spanmax-spanmin+j+1)/256, stdout);

	    j *= 3;
	    while (j--)
		putc('\0', stdout);

	    for (i = spanmin; i <= spanmax; i++)
		for (j = 2; j >= 0; j--)
		    putc(outspan.cols[j*dxsize + i], stdout);
    }
    fputs("\r\n", stdout);
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
    	
    fputs(chrtype[(p->arg0 >> 2) & 017], stdout);
    fputs(p->args, stdout);
    fputs(PNORM, stdout);
    putc('\r', stdout);

}
