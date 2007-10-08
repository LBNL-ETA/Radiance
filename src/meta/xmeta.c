#ifndef lint
static const char	RCSid[] = "$Id: xmeta.c,v 1.4 2007/10/08 18:07:57 greg Exp $";
#endif
/*
 *  Program to output meta-files to X window system.
 *
 *  cc -o Xmeta Xmeta.c Xplot.o plot.o mfio.o syscalls.o misc.o
 *
 *     2/26/86
 */

#include  "rtprocess.h"
#include  "string.h"
#include  "meta.h"
#include  "plot.h"


#define  overlap(p,xmn,ymn,xmx,ymx)  (ov((p)->xy[XMN],(p)->xy[XMX],xmn,xmx) \
				     &&ov((p)->xy[YMN],(p)->xy[YMX],ymn,ymx))

#define  ov(mn1,mx1,mn2,mx2)  ((mn1)<(mx2)&&(mn2)<(mx1))

#define  XCOM  "pexpand +OCIsv -P %s"


char  *progname;

static short  newpage = TRUE; 

static PLIST  recording;

int maxalloc = 0;		/* no limit */

extern void	pXFlush(void);


void
save(p)				/* record primitive */
PRIMITIVE  *p;
{
    register PRIMITIVE  *pnew;

    if ((pnew = palloc()) == NULL)
	    error(SYSTEM, "out of memory in save");
    mcopy((char *)pnew, (char *)p, sizeof(PRIMITIVE));
    add(pnew, &recording);
}


void
doglobal(g, sf)			/* execute a global command */
register PRIMITIVE  *g;
int  sf;
{

    switch (g->com) {

	case PEOF:
	    return;

	case PDRAW:
	    pXFlush();
	    break;

	case PEOP:
	    endpage();
	    plfree(&recording);
	    set(SALL, NULL);
	    newpage = TRUE;
	    return;

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
	    return;
	}
	if (sf)
	    save(g);

}


void
doprim(p, sf)		/* plot primitive */
register PRIMITIVE  *p;
int  sf;
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
    if (sf) {
	save(p);
	newpage = FALSE;
    }

}


void
plot(infp)		/* plot meta-file */
FILE  *infp;
{
    PRIMITIVE  nextp;

    set(SALL, NULL);
    do {
	readp(&nextp, infp);
	while (isprim(nextp.com)) {
	    doprim(&nextp, 1);
	    readp(&nextp, infp);
	}
	doglobal(&nextp, 1);
    } while (nextp.com != PEOF);

}


void
replay(xmin, ymin, xmax, ymax)		/* play back region */
int  xmin, ymin, xmax, ymax;
{
    register PRIMITIVE  *p;

    unset(SALL);
    set(SALL, NULL);
    for (p = recording.ptop; p != NULL; p = p->pnext)
	if (isprim(p->com)) {
	    if (overlap(p, xmin, ymin, xmax, ymax))
		doprim(p, 0);
	} else
	    doglobal(p, 0);

}


int
main(argc, argv)
int  argc;
char  **argv;
{
 FILE  *fp;
 char  *geometry = NULL;
 short  condonly, conditioned;
 char  comargs[500], command[600];

 progname = *argv++;
 argc--;

 condonly = FALSE;
 conditioned = FALSE;

 for ( ; argc; argv++, argc--)
    if (!strcmp(*argv, "-c"))
	condonly = TRUE;
    else if (!strcmp(*argv, "-r"))
	conditioned = TRUE;
    else if (**argv == '=')
	geometry = *argv;
    else
	break;

 if (conditioned) {
    init(progname, geometry);
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
       init(progname, geometry);
       if ((fp = popen(command, "r")) == NULL)
          error(SYSTEM, "cannot execute input filter");
       plot(fp);
       pclose(fp);
       }
    }

 if (!newpage)
     endpage();

 return(0);
 }
