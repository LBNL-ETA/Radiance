#ifndef lint
static const char	RCSid[] = "$Id: psmeta.c,v 1.4 2003/11/15 02:13:37 schorsch Exp $";
#endif
/*
*  Program to convert meta-files to PostScript.
*
*     9/23/88
*/


#include  "meta.h"
#include  "plot.h"
#include  "psplot.h"


char  *progname;

static short  newpage = TRUE;

static void doprim(PRIMITIVE  *p);
static int doglobal(PRIMITIVE  *g);

int
main(
	int  argc,
	char  **argv
)
{
	FILE  *fp;

	progname = *argv++;
	argc--;

	init(progname);
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

	if (!newpage)
		endpage();

	done();
	return(0);
}



void
plot(		/* plot meta-file */
	FILE  *infp
)
{
	PRIMITIVE  nextp;

	do {
		readp(&nextp, infp);
		while (isprim(nextp.com)) {
			doprim(&nextp);
			readp(&nextp, infp);
		}
	} while (doglobal(&nextp));

}


int
doglobal(			/* execute a global command */
	PRIMITIVE  *g
)
{
	FILE  *fp = NULL;

	switch (g->com) {

	case PEOF:
		return(0);

	case PPAUS:
		break;

	case PINCL:
		if (g->args == NULL)
		    error(USER, "missing include file name in include");
		if (g->arg0 == 2 || (fp = fopen(g->args, "r")) == NULL) {
		    if (g->arg0 != 0)
			fp = mfopen(g->args, "r");
		    else {
			sprintf(errmsg, "cannot open user include file \"%s\"",
					g->args);
			error(USER, errmsg);
		    }
		}
		plot(fp);
		fclose(fp);
		break;

	case PDRAW:
		fflush(stdout);
		break;

	case PEOP:
		endpage();
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

	case POPEN:
		segopen(g->args);
		break;

	case PCLOSE:
		segclose();
		break;

	default:
		sprintf(errmsg, "unknown command '%c' in doglobal", g->com);
		error(WARNING, errmsg);
		break;
	}

	return(1);
}



void
doprim(		/* plot primitive */
	PRIMITIVE  *p
)
{

	switch (p->com) {

	case PMSTR:
		printstr(p);
		break;

	case PVSTR:
		plotvstr(p);
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

	case PSEG:
		doseg(p);
		break;

	default:
		sprintf(errmsg, "unknown command '%c' in doprim", p->com);
		error(WARNING, errmsg);
		return;
	}
	newpage = FALSE;
}

