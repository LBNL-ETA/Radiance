#ifndef lint
static const char	RCSid[] = "$Id: impress.c,v 1.3 2003/10/27 10:28:59 schorsch Exp $";
#endif
/*
 *  Program to convert meta-files to imPress format
 *
 *
 *     1/2/86
 */


#include  "rtprocess.h"
#include  "meta.h"
#include  "imPcodes.h"
#include  "imPfuncs.h"


#define  XCOM  "pexpand +OCIsv -P %s"


char  *progname;

FILE  *imout;

static short  newpage = TRUE; 



main(argc, argv)

int  argc;
char  **argv;

{
 FILE  *fp;
 short  condonly, conditioned;
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

 imInit();

 if (conditioned) {
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
       if ((fp = popen(command, "r")) == NULL)
          error(SYSTEM, "cannot execute input filter");
       plot(fp);
       pclose(fp);
       }
    }

 if (!newpage)
     imEndPage();
 imEof();

 return(0);
 }



plot(infp)		/* plot meta-file */

register FILE  *infp;

{
    PRIMITIVE  nextp;

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

}






doglobal(g)			/* execute a global command */

register PRIMITIVE  *g;

{
    char  c;

    switch (g->com) {

	case PEOF:
	    break;

	case PDRAW:
	    fflush(stdout);
	    break;

	case PEOP:
	    if (!newpage)
		imEndPage();		/* don't waste paper */
	    newpage = TRUE;
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

    switch (p->com) {

	case PMSTR:
	    printstr(p);
	    break;

	case PLSEG:
	    plotlseg(p);
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

    newpage = FALSE;

}
