#ifndef lint
static const char	RCSid[] = "$Id: impress.c,v 1.4 2003/11/15 02:13:37 schorsch Exp $";
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
#include  "plot.h"


#define  XCOM  "pexpand +OCIsv -P %s"


char  *progname;

FILE  *imout;

static short  newpage = TRUE; 

static void doprim(register PRIMITIVE  *p);
static void doglobal(register PRIMITIVE  *g);



int
main(
	int  argc,
	char  **argv
)

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


void
plot(		/* plot meta-file */
	register FILE  *infp
)

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



void
doglobal(			/* execute a global command */
	register PRIMITIVE  *g
)

{
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



void
doprim(		/* plot primitive */
	register PRIMITIVE  *p
)

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
