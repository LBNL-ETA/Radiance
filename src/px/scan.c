#ifndef lint
static const char	RCSid[] = "$Id: scan.c,v 1.4 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 *  writescan.c - fortran interface to picture output routines.
 *
 *	4/26/88
 */

#include <stdio.h>


static FILE	*outfp;			/* output file pointer */
static char	outfile[128];		/* output file name */


initscan_(fname, width, height)		/* initialize output file */
char	*fname;
int	*width, *height;
{
	extern char	*strcpy();

	if (fname == NULL || fname[0] == '\0') {
		outfp = stdout;
		strcpy(outfile, "<stdout>");
	} else {
		if ((outfp = fopen(fname, "w")) == NULL) {
			perror(fname);
			exit(1);
		}
		strcpy(outfile, fname);
	}
	fprintf(outfp, "%dx%d picture\n\n-Y %d +X %d\n",
			*width, *height, *height, *width);
}


writescan_(scan, width)			/* output scanline */
float	*scan;
int	*width;
{
	if (fwritescan(scan, *width, outfp) < 0) {
		perror(outfile);
		exit(1);
	}
}


donescan_()				/* clean up */
{
	if (fclose(outfp) < 0) {
		perror(outfile);
		exit(1);
	}
}
