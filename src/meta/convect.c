#ifndef lint
static const char	RCSid[] = "$Id: convect.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  Convert vector file of type A to segments
 *
 *   1/16/85
 *
 *   cc convect.c primout.o hfio.o syscalls.o misc.o
 */


#include  "meta.h"

#define  MAXX  32

#define  MAXY  32

#define  MINX  0

#define  MINY  0

#define  CVX(x)  ((int)((long)(x-MINX)*(XYSIZE-1)/(MAXX-MINX)))

#define  CVY(y)  ((int)((long)(y-MINY)*(XYSIZE-1)/(MAXY-MINY)))

char  *progname;


main(argc, argv)	/* single argument is input file name */

int  argc;
char  **argv;

{
    char  stemp[16];
    int  i, npairs;
    FILE  *fp;
    PRIMITIVE  curp;
    int  x, y, lastx, lasty;

    progname = *argv++;
    argc--;

    if (argc != 1)
	error(USER, "arg count");
    
    fp = efopen(*argv, "r");

    while (fscanf(fp, "%s", stemp) == 1) {

	fscanf(fp, "%*d %*d %*d %*d");
	curp.com = POPEN;
	curp.arg0 = 0;
	curp.xy[XMN] = curp.xy[YMN] = curp.xy[XMX] = curp.xy[YMX] = -1;
	curp.args = stemp;
	writep(&curp, stdout);

	for ( ; ; ) {
	    fscanf(fp, "%d", &i);
	    if (i == 0)
		break;
	    fscanf(fp, "%d", &npairs);
	    for (i = 0; i < npairs; i++) {
	        fscanf(fp, "%d %d", &x, &y);
		x = CVX(x);
		y = CVY(y);
		if (i)
		    plseg(0, lastx, lasty, x, y, stdout);
		lastx = x;
		lasty = y;
	    }
	}
	curp.com = PCLOSE;
	curp.arg0 = 0200;
	curp.xy[XMN] = curp.xy[YMN] = curp.xy[XMX] = curp.xy[YMX] = -1;
	curp.args = NULL;
	writep(&curp, stdout);
    }

    writeof(stdout);
    fclose(fp);

    return(0);
}
