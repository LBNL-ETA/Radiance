#ifndef lint
static const char	RCSid[] = "$Id: panim.c,v 2.4 2003/10/27 10:24:51 schorsch Exp $";
#endif
/*
 *  Send pictures to PC animation system.
 *
 *	6/20/88
 */

#include <stdio.h>
#include <string.h>

#include "rtprocess.h"
#include "random.h"
#include "color.h"
#include "clntrpc.h"
#include "scan.h"

#define GAMMA		2.0		/* gamma correction factor */


char	*pcom = NULL;			/* pipe command */

BYTE	gammamap[256];			/* gamma correction table */


main(argc, argv)
int	argc;
char	*argv[];
{
	char	*progname;
	int	nframes;
	int	port = PCPROGRAM;
					/* initialize */
	compgamma();
	for (progname = *argv++; --argc; argv++)
		if (!strcmp(*argv, "-p") && argv[1]) {
			port = atoi(*++argv); argc--;
		} else if (!strcmp(*argv, "-u") && argv[1]) {
			pcom = *++argv; argc--;
		} else
			break;
	if (!argc) {
		fputs("Usage: ", stderr);
		fputs(progname, stderr);
		fputs(" [-p port] [-u uncompress] hostname [-c copies][-r record] [frame] ..\n",
				stderr);
		exit(1);
	}
	scry_open(*argv++, TCP, port);
	scry_set_compress(NONE);
	scry_set_image(TARGA_IMAGE);
	scry_set_record(PREVIEW);
	scry_set_copy_num(3);
					/* send frames */
	nframes = 0;
	for ( ; --argc; argv++)
		if (!strcmp(*argv, "-r") && argv[1]) {
			scry_set_record(atoi(*++argv)); argc--;
		} else if (!strcmp(*argv, "-c") && argv[1]) {
			scry_set_copy_num(atoi(*++argv)); argc--;
		} else {
			sendframe(*argv);
			nframes++;
		}
	if (nframes == 0)		/* get stdin */
		sendframe(NULL);
					/* clean up */
	scry_close();
	exit(0);
}


sendframe(file)			/* convert and send a frame */
char	*file;
{
	char	command[128];
	COLR	scanin[SCANLINE];
	int	xres, yres;
	int	xbeg, ybeg;
	FILE	*fp;
	int	y;
	register int	x;
						/* open file */
	if (file == NULL) {
		if (pcom != NULL)
			fp = popen(pcom, "r");
		else
			fp = stdin;
		file = "<stdin>";
	} else {
		if (pcom != NULL) {
			sprintf(command, "( %s ) < %s", pcom, file);
			fp = popen(command, "r");
		} else
			fp = fopen(file, "r");
	}
	if (fp == NULL) {
		perror(file);
		exit(1);
	}
						/* get dimensions */
	getheader(fp, NULL);
	if (checkheader(fp, COLRFMT, NULL) < 0) {
		fputs(file, stderr);
		fputs(": not a Radiance picture\n", stderr);
		exit(1);
	}
	if (fgetresolu(&xres, &yres, fp) != (YMAJOR|YDECR) ||
			xres > SCANLINE || yres > NUMSCANS) {
		fputs(file, stderr);
		fputs(": bad picture size\n", stderr);
		exit(1);
	}
						/* compute center */
	xbeg = (SCANLINE-xres)/2;
	ybeg = (NUMSCANS-yres)/2;
						/* clear output */
	memset(sc_frame_arr, '\0', sizeof(sc_frame_arr));
						/* get frame */
	for (y = yres-1; y >= 0; y--) {
		if (freadcolrs(scanin, xres, fp) < 0) {
			fputs(file, stderr);
			fputs(": read error\n", stderr);
			exit(1);
		}
		normcolrs(scanin, xres, 0);	/* normalize */
		for (x = 0; x < xres; x++)	/* convert */
			sc_frame_arr[y+ybeg][x+xbeg] =
			((gammamap[scanin[x][RED]]+(random()&7))&0xf8)<<7
			| ((gammamap[scanin[x][GRN]]+(random()&7))&0xf8)<<2
			| (gammamap[scanin[x][BLU]]+(random()&7))>>3;
	}
						/* send frame */
	scry_send_frame();
						/* close file */
	if (pcom != NULL)
		pclose(fp);
	else
		fclose(fp);
}


compgamma()				/* compute gamma correction map */
{
	extern double  pow();
	register int  i, val;

	for (i = 0; i < 256; i++) {
		val = pow((i+0.5)/256.0, 1.0/GAMMA) * 256.0;
		if (val > 248) val = 248;
		gammamap[i] = val;
	}
}
