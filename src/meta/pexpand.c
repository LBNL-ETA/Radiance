#ifndef lint
static const char	RCSid[] = "$Id: pexpand.c,v 1.3 2003/11/15 02:13:37 schorsch Exp $";
#endif
/*
 *   Program to expand meta-file commands
 *
 *   cc pexpand.c expand.o mfio.o segment.o palloc.o syscalls.o misc.o
 */


#include  "meta.h"


char  *progname;

int  maxalloc = 0;		/* no limit */


int
main(
	int  argc,
	char  **argv
)
{
 FILE  *fp;
 int  i;
 char  *cp;
 int  com;
 short  exlist[NCOMMANDS];	/* 1==expand, 0==pass, -1==discard */

 progname = *argv++;
 argc--;

 for (i = 0; i < NCOMMANDS; i++)
    exlist[i] = 0;

 while (argc && (**argv == '+' || **argv == '-'))  {
    i = (**argv == '+') ? 1 : -1;
    for (cp = *argv+1; *cp ; cp++)  {
       if ((com = comndx(*cp)) == -1 || *cp == PEOF) {
          sprintf(errmsg, "unknown option '%c'", *cp);
	  error(WARNING, errmsg);
	  }
       else
	  exlist[com] = i;
       }
    argv++;
    argc--;
    }

 if (argc)
    while (argc)  {
       fp = efopen(*argv, "r");
       expand(fp, exlist);
       fclose(fp);
       argv++;
       argc--;
       }
 else
    expand(stdin, exlist);

 writeof(stdout);

 return(0);
 }
