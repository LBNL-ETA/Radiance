#ifndef lint
static const char	RCSid[] = "$Id: plot4.c,v 1.4 2003/11/15 02:13:37 schorsch Exp $";
#endif
/*
 *  plot4.c - program to put four metafile pages onto one.
 *
 *	Greg Ward
 *	7/10/86
 */

#include  "rtprocess.h"
#include  "meta.h"
#include  "plot.h"


#define  OUTFILT  "pexpand +OCIms"	/* output filter */

#define  SEGNAME  "plot4seg"		/* segment name */


extern FILE  *pout;			/* the output stream */

char  *progname;

static void plot4(FILE  *fp);
static void doseg(int  n);


int
main(
	int  argc,
	char  *argv[]
)
{
	FILE  *fp;
	int  i;

	progname = argv[0];

	pout = popen(OUTFILT, "w");

	if (argc > 1)
		for (i = 1; i < argc; i++) {
			fp = efopen(argv[i], "r");
			plot4(fp);
			fclose(fp);
		}
	else
		plot4(stdin);
	
	pglob(PEOF, 0200, NULL);

	return(pclose(pout));
}


void
plot4(			/* put a file into its place on page */
	FILE  *fp
)
{
	static int  nplts = 0;
	PRIMITIVE  curp;

	pglob(POPEN, 0, SEGNAME);

	while (readp(&curp, fp))
		if (curp.com == PEOP) {
			pglob(PCLOSE, 0200, NULL);
			doseg(nplts++ % 4);
			pglob(POPEN, 0, SEGNAME);
		} else
			writep(&curp, pout);
	
	pglob(PCLOSE, 0200, NULL);
}


void
doseg(		/* do segment number n */
	int  n
)
{
	switch (n) {
	case 0:				/* upper left */
		pprim(PSEG, 0, 0, XYSIZE/2, XYSIZE/2-1, XYSIZE-1, SEGNAME);
		break;
	case 1:				/* upper right */
		pprim(PSEG, 0, XYSIZE/2, XYSIZE/2, XYSIZE-1, XYSIZE-1, SEGNAME);
		break;
	case 2:				/* lower left */
		pprim(PSEG, 0, 0, 0, XYSIZE/2-1, XYSIZE/2-1, SEGNAME);
		break;
	case 3:				/* lower right, end page */
		pprim(PSEG, 0, XYSIZE/2, 0, XYSIZE-1, XYSIZE/2-1, SEGNAME);
		pglob(PEOP, 0200, NULL);
		break;
	}
}
