#ifndef lint
static const char	RCSid[] = "$Id: pranim.c,v 1.3 2003/10/27 10:32:06 schorsch Exp $";
#endif
/*
 *  Send Sun rasterfiles to PC animation system.
 *
 *	6/20/88
 */

#include <stdio.h>
#include <rasterfile.h>

#include "rtprocess.h"
#include "client/clnt.h"


char	*pcom = NULL;			/* uncompress command */


main(argc, argv)
int	argc;
char	**argv;
{
	int	startrec = 0;
	char	*progname;
	int	prognum = PCPROGRAM;
					/* initialize */
	for (progname = *argv++; --argc; argv++)
		if (!strcmp(*argv, "-p") && argv[1]) {
			prognum = atoi(*++argv); argc--;
		} else if (!strcmp(*argv, "-u") && argv[1]) {
			pcom = *++argv; argc--;
		} else if (!strcmp(*argv, "-r") && argv[1]) {
			startrec = atoi(*++argv); argc--;
		} else
			break;
	if (!argc) {
		fputs("Usage: ", stderr);
		fputs(progname, stderr);
		fputs(" [-p prognum] [-u uncompress] [-r record] hostname [frame..]\n",
				stderr);
		exit(1);
	}
        scry_open(*argv++, prognum); argc--;
	scry_get_info(1) ;
	scry_set_compress(NONE);
	if (startrec > 0)
		scry_init_record(startrec, argc>0?argc:1);
	else
		scry_init_record(PREVIEW,0);
					/* send frames */
	if (argc <= 0)
		sendframe(NULL);
	else
		for ( ; argc > 0; argc--, argv++)
			sendframe(*argv);
					/* clean up */
	scry_close();
	exit(0);
}


sendframe(file)			/* convert and send a frame */
char	*file;
{
	static struct rasterfile	head;
	char	command[128];
	unsigned char	cmap[3][256];
	unsigned char   dmap[256][3] ;
	static unsigned short	targa[256];
	static int	xbeg, ybeg;
	register FILE	*fp;
	register int	i, j;
	register int	c;
	/* DAVIDR */
	unsigned char *image_read, *newimage, *compdata ;
	int new_ht, new_wd, new_depth ;
	/* end DAVIDR */
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
						/* check format */
	if (fread(&head, sizeof(head), 1, fp) != 1)
		goto rasfmt;
	if (head.ras_magic != RAS_MAGIC || head.ras_type != RT_STANDARD)
		goto rasfmt;
	if (head.ras_maptype != RMT_EQUAL_RGB || head.ras_maplength != 3*256)
		goto rasfmt;
						/* get color table */
	if (fread(cmap, 256, 3, fp) != 3)
		goto readerr;
						/* convert table */
#ifdef LATER
	for (j = 0; j < 253; j++)
		targa[j+1] = (cmap[0][j] & 0xf8)<<7 |
				(cmap[1][j] & 0xf8)<<2 | cmap[2][j]>>3;
#endif
        for (j = 0 ; j < 256 ; j++)
	{
	    dmap[j][0] = cmap[0][j] ;
	    dmap[j][1] = cmap[1][j] ;
	    dmap[j][2] = cmap[2][j] ;
	}
						/* send table */
	scry_set_map(254, dmap);
	image_read = (unsigned char *) malloc(head.ras_width*head.ras_height) ;
	for (j = 0 ; j < head.ras_height ; j++)
	    fread (&(image_read[(head.ras_height-j-1)*head.ras_width]),head.ras_width,1,fp) ;
	fread(image_read,head.ras_width*head.ras_height,1,fp) ;
	newimage = NULL ;
						/* send frame */
	scry_conv_vis (image_read,&newimage,head.ras_height,head.ras_width,&new_ht,&new_wd,1,&new_depth) ;
	free (image_read) ;
	compdata = NULL ;
	scry_compress(newimage,&compdata,new_ht,new_wd,new_depth) ;
	scry_send_frame(compdata);
	free(compdata) ;
	free (newimage) ;
						/* close file */
	if (pcom != NULL)
		pclose(fp);
	else
		fclose(fp);
	return;
rasfmt:
	fputs(file, stderr);
	fputs(": wrong rasterfile format\n", stderr);
	exit(1);
readerr:
	fputs(file, stderr);
	fputs(": read error\n", stderr);
	exit(1);
}
