#ifndef lint
static const char	RCSid[] = "$Id: macprt.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  Program to send meta-files to MacIntosh printer
 *
 *  cc -o macprt macprt.c macplot.o mfio.o syscalls.o misc.o -lplot
 */


#include  <fcntl.h>

#include  "meta.h"

#include  "macplot.h"

#include  <print.h>



char  *progname;

int  dxsize, dysize;		/* plot dimensions */

static short  newpage = TRUE; 

static PRIMITIVE  nextp;

static GrafPort  ourPort;



main(argc, argv)

int  argc;
char  **argv;

{
 FILE  *fp;

 progname = *argv++;
 argc--;

 while (argc && **argv == '-')  {
    switch (*(*argv+1))  {
       default:
	  error(WARNING, "unknown option");
	  break;
       }
    argv++;
    argc--;
    }

 InitGraf(&thePort);

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

 return(0);
}




plot(infp)		/* plot meta-file */

register FILE  *infp;

{

    OpenPort(&ourPort);
    TextFont(4);		/* 9-point Monaco */
    TextSize(9);

    dxsize = dysize = min(ourPort.portRect.right - ourPort.portRect.left,
			ourPort.portRect.bottom - ourPort.portRect.top);

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

    ClosePort(&ourPort);

}





doglobal(g)			/* execute a global command */

register PRIMITIVE  *g;

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
		/* fall through */

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
	    } while (c != '\n' && c != '\r');
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




doprim(p)		/* plot primitive */

register PRIMITIVE  *p;

{

    if (newpage) {
	EraseRgn(ourPort.visRgn);
	newpage = FALSE;
    }
    
    switch (p->com) {

	case PLSEG:
	    plotlseg(p);
	    break;

	case PMSTR:
	    printstr(p);
	    break;

	case PRFILL:
	    fillrect(p);
	    break;

	case PTFILL:
	    filltri(p);
	    break;

	case PPFILL:
	    fillpoly(p);
	    break;

	default:
	    sprintf(errmsg, "unknown command '%c' in doprim", p->com);
	    error(WARNING, errmsg);
	    return;
    }

}
