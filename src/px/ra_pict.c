#ifndef lint
static const char	RCSid[] = "$Id: ra_pict.c,v 2.8 2003/02/22 02:07:27 greg Exp $";
#endif
/*		Convert an Radiance image to APPLE pict format.
 *
 *			Orginally Iris to PICT by	Paul Haeberli - 1990
 *			Hacked into Rad to PICT by Russell Street 1990
 *
 *	History:
 *	    V 1			-- Does basic conversion
 *	    V 1.1 (2/11/91)	-- Added options for Gamma
 *				-- verbose option
 *				-- man page
 *				-- better about allocating buffers
 *	    V 1.2 (19/11/91)	-- Revised to handle opening "The Radiance Way"
 *				-- Added exposure level adjustment
 */

#include <stdio.h>
#include <math.h>
#ifdef MSDOS
#include <fcntl.h>
#endif
#include  <time.h>

#include "pict.h"
#include "color.h"
#include "resolu.h"

int	outbytes;		    /* This had better be 32 bits! */
char	*progname;
int	verbose = 0;
float	gamcor = 2.0;
int	bradj = 0;

	/* First some utility routines */

putrect(xorg,yorg,xsize,ysize)
int xorg, yorg, xsize, ysize;
{
    putashort(yorg);
    putashort(xorg);
    putashort(ysize);
    putashort(xsize);
}

putfprect(xorg,yorg,xsize,ysize)
int xorg, yorg, xsize, ysize;
{
    putalong(yorg<<16);
    putalong(xorg<<16);
    putalong(ysize<<16);
    putalong(xsize<<16);
}

putalong(l)
long l;
{
    putbyte((l>>24)&0xff);
    putbyte((l>>16)&0xff);
    putbyte((l>>8)&0xff);
    putbyte((l>>0)&0xff);
}

putashort(s)
short s;
{
    putbyte((s>>8)&0xff);
    putbyte((s>>0)&0xff);
}

putbyte(b)
int b;
{
    if (putc(b,stdout) == EOF && ferror(stdout)) {
	fprintf(stderr,"%s: error on write\n", progname);
	exit(1);
    }
    outbytes++;
}

putbytes(buf,n)
unsigned char *buf;
int n;
{
    if(!fwrite(buf,n,1,stdout)) {
	fprintf(stderr,"%s: error on write\n", progname);
	exit(1);
    }
    outbytes+=n;
}

main(argc,argv)
int argc;
char **argv;
{
    int xsize, ysize;
    int i, picsize;
    int ssizepos, lsizepos;
#ifdef MSDOS
    extern int	_fmode;
    _fmode = O_BINARY;
    setmode(fileno(stdin), O_BINARY);
    setmode(fileno(stdout), O_BINARY);
#endif
    progname = argv[0];

    for (i = 1; i < argc ; i++)
	if (argv[i][0] ==  '-')
	    switch (argv[i][1]) {
		case 'g':	gamcor = atof(argv[++i]);
				break;

		case 'e':	if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
				   usage();
				else
				    bradj = atoi(argv[++i]);
				break;

		case 'v':	;
				verbose = 1;
				break;

		case 'r':	fprintf(stderr, "Sorry. Get a Macintosh :-)\n");
				exit(1);

		case '-':	i++;
				goto outofparse;
				break;		    /* NOTREACHED */

		otherwise:	usage();
				break;
		 }
	else
	    break;

outofparse:

    if (i < argc - 2)
	usage();

    if (i <= argc - 1 && freopen(argv[i], "r", stdin) == NULL) {
	fprintf(stderr, "%s: can not open input \"%s\"\n",
	    progname, argv[i]);
	exit(1);
    }

    if (i <= argc - 2 && freopen(argv[i+1], "w", stdout) == NULL) {
	fprintf(stderr, "%s: can not open input \"%s\"\n",
	    progname, argv[i+1]);
	exit(1);
    }
	
#ifdef DEBUG
	fprintf(stderr, "Input file: %s\n", i <= argc - 1 ? argv[i] : "stdin");
	fprintf(stderr, "Outut file: %s\n", i <= argc - 2 ? argv[i+1] : "stdout" );
	fprintf(stderr, "Gamma: %f\n", gamcor);
	fprintf(stderr, "Brightness adjust: %d\n", bradj);
	fprintf(stderr, "Verbose: %s\n", verbose ? "on" : "off");
#endif


	     /* OK. Now we read the size of the Radiance picture */
    if (checkheader(stdin, COLRFMT, NULL) < 0 ||
	    fgetresolu(&xsize, &ysize, stdin) < 0 /* != (YMAJOR|YDECR) */ ) {
	fprintf(stderr, "%s: not a radiance picture\n", progname);
	exit(1);
	}

	    /* Set the gamma correction */

    setcolrgam(gamcor);

    for(i=0; i<HEADER_SIZE; i++) 
	putbyte(0);

    ssizepos = outbytes;
    putashort(0);		/* low 16 bits of file size less HEADER_SIZE */
    putrect(0,0,xsize,ysize);	/* bounding box of picture */
    putashort(PICT_picVersion);
    putashort(0x02ff);		/* version 2 pict */
    putashort(PICT_reservedHeader);	/* reserved header opcode */

    lsizepos = outbytes;
    putalong(0);		/* full size of the file */
    putfprect(0,0,xsize,ysize); /* fixed point bounding box of picture */
    putalong(0);		/* reserved */

    putashort(PICT_clipRgn);	/* the clip region */
    putashort(10);
    putrect(0,0,xsize,ysize);

    if (verbose)
	fprintf(stderr, "%s: The picture is %d by %d, with a gamma of %f\n",
	    progname, xsize, ysize, gamcor);


    putpict(xsize, ysize);	/* Here is where all the work is done */

    putashort(PICT_EndOfPicture); /* end of pict */

    picsize = outbytes-HEADER_SIZE;
    fseek(stdout,ssizepos,0);
    putashort(picsize&0xffff);
    fseek(stdout,lsizepos,0);
    putalong(picsize);

    fclose(stdout);
    fclose(stdin);
    
    exit(0);
    return 0;	    /* lint fodder */
}

putpict(xsize, ysize)
int xsize;
int ysize;
{
    int	    y;
    int	    nbytes, rowbytes;
    char    *cbuf, *pbuf;

    cbuf = (char *)malloc(4 * xsize);

    if (cbuf == NULL) {
	fprintf(stderr, "%s: not enough memory\n", progname);
	exit(1);
	}

    pbuf = (char *)malloc(4 * xsize);

    if (pbuf == NULL) {
	fprintf(stderr, "%s: not enough memory\n", progname);
	exit(1);
	}

    putashort(PICT_Pack32BitsRect); /* 32 bit rgb */
    rowbytes = 4*xsize;
    putalong(0x000000ff);		/* base address */


    if(rowbytes&1)
	rowbytes++;
    putashort(rowbytes|0x8000); /* rowbytes */
    putrect(0,0,xsize,ysize);	/* bounds */
    putashort(0);		/* version */

    putashort(4);	/* packtype */
    putalong(0);	/* packsize */
    putalong(72<<16);	/* hres */
    putalong(72<<16);	/* vres */

    putashort(16);	/* pixeltype */
    putashort(32);	/* pixelsize */
    putashort(3);	/* cmpcount */


    putashort(8);	/* cmpsize */
    putalong(0);	/* planebytes */
    putalong(0);	/* pmtable */
    putalong(0);	/* pmreserved */


    putrect(0,0,xsize,ysize);	/* scr rect */
    putrect(0,0,xsize,ysize);	/* dest rect */


    putashort(0x40);	/* transfer mode */

    for(y=0; y<ysize; y++) {
	getrow(stdin, cbuf, xsize);

	nbytes = packbits(cbuf,pbuf,24*xsize);
	if(rowbytes>250) 
	    putashort(nbytes);
	else
	    putbyte(nbytes);
	putbytes(pbuf,nbytes);
	}

    if(outbytes&1) 
	putbyte(0);

    free(cbuf);
    free(pbuf);
}

int getrow(in, cbuf, xsize)
FILE *in;
char *cbuf;
int xsize;
{
    extern char *tempbuffer();		/* defined in color.c */
    COLR    *scanin = NULL;
    int	    x;

    if ((scanin = (COLR *)tempbuffer(xsize*sizeof(COLR))) == NULL) {
	fprintf(stderr, "%s: not enough memory\n", progname);
	exit(1);
    }

    if (freadcolrs(scanin, xsize, in) < 0) {
	fprintf(stderr, "%s: read error\n", progname);
	exit(1);
	}

    if (bradj)	    /* Adjust exposure level */
	shiftcolrs(scanin, xsize, bradj);


    colrs_gambs(scanin, xsize);	    /* Gamma correct it */
    
    for (x = 0; x < xsize; x++) {
	cbuf[x] = scanin[x][RED];
	cbuf[xsize + x] = scanin[x][GRN];
	cbuf[2*xsize + x] = scanin[x][BLU];
	}

}


packbits(ibits,pbits,nbits)
unsigned char *ibits, *pbits;
int nbits;
{
    int bytes;			    /* UNUSED */
    unsigned char *sptr;
    unsigned char *ibitsend;
    unsigned char *optr = pbits;
    int nbytes, todo, cc, count;

    nbytes = ((nbits-1)/8)+1;
    ibitsend = ibits+nbytes;
    while(ibits<ibitsend) {
	sptr = ibits;
	ibits += 2;
	while((ibits<ibitsend)&&((ibits[-2]!=ibits[-1])||(ibits[-1]!=ibits[0])))
	    ibits++;
	if(ibits != ibitsend) {
	    ibits -= 2;
	}
	count = ibits-sptr;
	while(count) {
	    todo = count>127 ? 127:count;
	    count -= todo;
	    *optr++ = todo-1;
	    while(todo--)
		*optr++ = *sptr++;
	}
	if(ibits == ibitsend)
	    break;
	sptr = ibits;
	cc = *ibits++;
	while( (ibits<ibitsend) && (*ibits == cc) )
	    ibits++;
	count = ibits-sptr;
	while(count) {
	    todo = count>128 ? 128:count;
	    count -= todo;
	    *optr++ = 257-todo;
	    *optr++ = cc;
	}
    }
    return optr-pbits;
}

usage()
{
    fprintf(stderr, "Usage: %s [-v] [-g gamma] [infile [outfile]]\n",
	progname);
    exit(2);
}
