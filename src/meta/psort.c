#ifndef lint
static const char	RCSid[] = "$Id: psort.c,v 1.3 2003/11/15 02:13:37 schorsch Exp $";
#endif
/*
 *   Program to sort a meta-file
 *
 *   cc psort.c sort.o palloc.o mfio.o syscalls.o misc.o
 */


#define  MAXALLOC  1100		/* must be >= PBSIZE in sort.c */


#include  "meta.h"


char  *progname;

int  maxalloc = MAXALLOC;

static int  val[5], extrema[5];


static int pcompare(PRIMITIVE  **pp1, PRIMITIVE  **pp2);


int
main(
	int  argc,
	char  **argv
)
{
 FILE  *fp;
 int  i;

 progname = *argv++;
 argc--;

 for (i = 0; i < 4 && argc && (**argv == '+' || **argv == '-'); i++)  {
    val[i] = (**argv == '+') ? 1 : -1;
    switch (*(*argv+1))  {
       case 'x':
	  extrema[i] = XMN;
	  break;
       case 'y':
	  extrema[i] = YMN;
	  break;
       case 'X':
	  extrema[i] = XMX;
	  break;
       case 'Y':
	  extrema[i] = YMX;
	  break;
       default:
          sprintf(errmsg, "unknown option \"%s\"", *argv);
	  error(USER, errmsg);
	  break;
       }
    argv++;
    argc--;
    }

 val[i] = 0;

 if (argc)
    while (argc)  {
       fp = efopen(*argv, "r");
       sort(fp, pcompare);
       fclose(fp);
       argv++;
       argc--;
       }
 else
    sort(stdin, pcompare);

 writeof(stdout);

 return(0);
 }




int
pcompare(
	PRIMITIVE  **pp1,
	PRIMITIVE  **pp2
)
{
 register PRIMITIVE  *p1 = *pp1, *p2 = *pp2;
 register int  i;

 for (i = 0; val[i]; i++)
    if (p1->xy[extrema[i]] > p2->xy[extrema[i]])
       return(val[i]);
    else if (p1->xy[extrema[i]] < p2->xy[extrema[i]])
       return(-val[i]);

 return(0);
 }
