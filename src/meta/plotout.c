#ifndef lint
static const char	RCSid[] = "$Id: plotout.c,v 1.2 2003/07/01 21:21:40 greg Exp $";
#endif
/*
 *  Program to send meta-files to plot(3X) drivers
 *
 *  cc -o plotout plotout.c mfio.o syscalls.o misc.o -lplot
 *
 *  Plot drivers
 */


#ifdef  FORTEK
#define  XCOM  "pexpand +vOCIsm -rtpSUR %s"
#else
#define  XCOM  "pexpand +vOCIs -rtpSUR %s"
#endif




#include  <fcntl.h>

#include  "meta.h"




char  *progname;

static short  newpage = TRUE; 

static int  curx = -1,		/* current position */
	    cury = -1;

static short  curmod = -1;	/* current line drawing mode */

static char  lmode[4][16] = {"solid", "shortdashed",
			     "dotted", "dotdashed"};

static PRIMITIVE  nextp;



main(argc, argv)

int  argc;
char  **argv;

{
 FILE  *fp, *popen();
 char  comargs[200], command[300];
 short  condonly = FALSE, conditioned = FALSE;

 progname = *argv++;
 argc--;

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

 if (conditioned)
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
 else  {
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
       if ((fp = popen(command, "r")) == NULL)
          error(SYSTEM, "cannot execute input filter");
       plot(fp);
       pclose(fp);
       }
    }

 return(0);
 }





plot(infp)		/* plot meta-file */

FILE  *infp;

{

    openpl();
    space(0, 0, XYSIZE, XYSIZE);

    do {
	readp(&nextp, infp);
	while (isprim(nextp.com)) {
	    doprim(&nextp);
	    fargs(&nextp);
	    readp(&nextp, infp);
	}
	doglobal(&nextp);
	fargs(&nextp);
    } while (nextp.com != PEOF);

    closepl();

}





doglobal(g)			/* execute a global command */

PRIMITIVE  *g;

{
    int  tty;
    char  c;

    switch (g->com) {

	case PEOF:
	    break;

	case PDRAW:
	    fflush(stdout);
	    break;

	case PEOP:
	    newpage = TRUE;
	    if (!isatty(fileno(stdout)))
		break;
	/* fall through */

	case PPAUS:
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
	    } while (c != '\n' && c != '\r');
	    close(tty);
	    break;

	default:
	    sprintf(errmsg, "unknown command '%c' in doglobal", g->com);
	    error(WARNING, errmsg);
	    break;
	}

}




doprim(p)		/* plot primitive */

register PRIMITIVE  *p;

{

    if (newpage) {
	erase();
	newpage = FALSE;
    }
    
    switch (p->com) {

	case PLSEG:
	    plotlseg(p);
	    break;

	case PMSTR:
	    printstr(p);
	    break;

	default:
	    sprintf(errmsg, "unknown command '%c' in doprim", p->com);
	    error(WARNING, errmsg);
	    return;
    }

}






printstr(p)		/* output a string */

register PRIMITIVE  *p;

{

    move(p->xy[XMN], p->xy[YMX]);
    label(p->args);
    curx = -1;
    cury = -1;

}





plotlseg(p)		/* plot a line segment */

register PRIMITIVE  *p;

{
    static short  right = FALSE;
    int  y1, y2;
    short  lm = (p->arg0 >> 4) & 03;

    if (p->arg0 & 0100) {
	y1 = p->xy[YMX];
	y2 = p->xy[YMN];
    } else {
	y1 = p->xy[YMN];
	y2 = p->xy[YMX];
    }

    if (lm != curmod) {
	linemod(lmode[lm]);
	curmod = lm;
    }

    if (p->xy[XMN] == curx && y1 == cury) {
	cont(p->xy[XMX], y2);
	curx = p->xy[XMX];
	cury = y2;
    } else if (p->xy[XMX] == curx && y2 == cury) {
	cont(p->xy[XMN], y1);
	curx = p->xy[XMN];
	cury = y1;
    } else if (right = !right) {
	line(p->xy[XMN], y1, p->xy[XMX], y2);
	curx = p->xy[XMX];
	cury = y2;
    } else {
	line(p->xy[XMX], y2, p->xy[XMN], y1);
	curx = p->xy[XMN];
	cury = y1;
    }
}
