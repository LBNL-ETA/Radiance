#ifndef lint
static const char	RCSid[] = "$Id: t16anim.c,v 1.3 2003/10/27 10:32:06 schorsch Exp $";
#endif
/*
 *  Send Targa 16-bit files to PC animation system.
 *
 *	6/20/88
 */

#include <stdio.h>

#include "rtprocess.h"
#include "client/clnt.h"
#include  "targa.h"

#define  goodpic(h)	(((h)->dataType==IM_RGB || (h)->dataType==IM_CRGB) \
				&& (h)->dataBits==16)


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
	scry_get_info(2) ;
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
	static struct hdStruct	head;
	char	command[128];
	register FILE	*fp;
	register int	i, j;
	register int	c;
	/* DAVIDR */
	unsigned char *image_read,*newimage,*compdata ;
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
	if (getthead(&head, NULL, fp) == -1 || !goodpic(&head))
		goto badfmt;
	/* DAVIDR */
	image_read = (unsigned char *) malloc((head.dataBits)/8 * head.x * head.y) ;
						/* get frame */
	for (i = 0 ; i < head.y ; i++)
	    fread(&(image_read[(head.y-i-1)*head.x*(head.dataBits/8)]),(head.dataBits)/8,head.x,fp) ;
	newimage = NULL ;
	/* convert to server size and visual */
	scry_conv_vis (image_read,&newimage,head.y,head.x,&new_ht,&new_wd,(head.dataBits)/8,&new_depth) ;
	free (image_read) ;
	compdata = NULL ;
	/* compress if necessary */
	scry_compress(newimage,&compdata,new_ht,new_wd,new_depth) ;
						/* send frame */
	scry_send_frame(compdata);
	free(newimage) ;
	free(compdata) ;
	/* end DAVIDR */
						/* close file */
	if (pcom != NULL)
		pclose(fp);
	else
		fclose(fp);
	return;
badfmt:
	fputs(file, stderr);
	fputs(": wrong targa format\n", stderr);
	exit(1);
readerr:
	fputs(file, stderr);
	fputs(": read error\n", stderr);
	exit(1);
}


getthead(hp, ip, fp)		/* read header from input */
struct hdStruct  *hp;
char  *ip;
register FILE  *fp;
{
	int	nidbytes;

	if ((nidbytes = getc(fp)) == EOF)
		return(-1);
	hp->mapType = getc(fp);
	hp->dataType = getc(fp);
	hp->mapOrig = getint2(fp);
	hp->mapLength = getint2(fp);
	hp->CMapBits = getc(fp);
	hp->XOffset = getint2(fp);
	hp->YOffset = getint2(fp);
	hp->x = getint2(fp);
	hp->y = getint2(fp);
	hp->dataBits = getc(fp);
	hp->imType = getc(fp);

	if (ip != NULL)
		if (nidbytes)
			fread(ip, nidbytes, 1, fp);
		else
			*ip = '\0';
	else if (nidbytes)
		fseek(fp, (long)nidbytes, 1);

	return(feof(fp) || ferror(fp) ? -1 : 0);
}




int
getint2(fp)			/* get a 2-byte positive integer */
register FILE	*fp;
{
	register int b1, b2;

	if ((b1 = getc(fp)) == EOF || (b2 = getc(fp)) == EOF) {
		fputs("Unexpected EOF\n", stderr);
		exit(1);
	}
	return(b1 | b2<<8);
}


