#ifndef lint
static const char	RCSid[] = "$Id: cv.c,v 1.2 2003/08/01 14:14:24 schorsch Exp $";
#endif
/*
 *    Human-readable file I/O conversion program
 *
 *    cc -o ../cv cv.c mfio.o cvhfio.o syscalls.o misc.o
 */


#include  "meta.h"



char  *progname;


main(argc, argv)

int  argc;
char  **argv;

{
 FILE  *fp;
 PRIMITIVE  curp;
 short  htom = TRUE;

 progname = *argv++;
 argc--;

 if (argc && **argv == '-') {
    htom = FALSE;
    argv++;
    argc--;
    }

 if (argc)
    for (; argc; argc--, argv++) {
       fp = efopen(*argv, "r");
       if (htom)
          while (scanp(&curp, fp)) {
             writep(&curp, stdout);
	     fargs(&curp);
             }
       else
          while (readp(&curp, fp)) {
             printp(&curp, stdout);
	     fargs(&curp);
             }
       fclose(fp);
       }
 else
    if (htom)
       while (scanp(&curp, stdin)) {
          writep(&curp, stdout);
	  fargs(&curp);
          }
    else
       while (readp(&curp, stdin)) {
          printp(&curp, stdout);
	  fargs(&curp);
          }

 if (htom)
    writeof(stdout);
    
 return(0);
 }
