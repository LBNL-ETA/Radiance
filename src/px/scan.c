#ifndef lint
static const char	RCSid[] = "$Id: scan.c,v 1.5 2004/01/02 12:45:36 schorsch Exp $";
#endif
/*
 *  writescan.c - fortran interface to picture output routines.
 *
 *	4/26/88
 */

#include <stdio.h>
#include <string.h>


static FILE	*outfp;			/* output file pointer */
static char	outfile[128];		/* output file name */


void
initscan_(		/* initialize output file */
	char	*fname,
	int	*width,
	int	*height
)
{
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


void
writescan_(			/* output scanline */
	float	*scan,
	int	*width
)
{
	if (fwritescan(scan, *width, outfp) < 0) {
		perror(outfile);
		exit(1);
	}
}


void
donescan_(void)				/* clean up */
{
	if (fclose(outfp) < 0) {
		perror(outfile);
		exit(1);
	}
}
