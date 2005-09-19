#ifndef lint
static const char RCSid[] = "$Id: meta2bmp.c,v 1.2 2005/09/19 12:48:09 schorsch Exp $";
#endif
/*
 *  Program to convert meta-files to BMP 8-bit color-mapped format
 */

#include  "copyright.h"

#include  "rtprocess.h"
#include  "meta.h"
#include  "plot.h"
#include  "rast.h"
#include  "bmpfile.h"
#include  "targa.h"

#define  MAXALLOC  30000
#define  DXSIZE  400		/* default x resolution */
#define  DYSIZE  400		/* default y resolution */
#define  XCOM  "pexpand +vOCImsp -DP %s | psort +y"


char  *progname;

SCANBLOCK	outblock;

int  dxsize = DXSIZE, dysize = DYSIZE;

int  maxalloc = MAXALLOC;

int  ydown = 0;

static char  outname[64];
static char  *outtack = NULL;

static BMPWriter  *bmpw = NULL;

static int  lineno = 0;

static short  condonly = FALSE,
	      conditioned = FALSE;




char *
findtack(s)			/* find place to tack on suffix */
register char *s;
{
	while (*s && *s != '.')
		s++;
	return(s);
}


int
main(
	int  argc,
	char  **argv
)

{
 FILE  *fp;
 char  comargs[200], command[300];

 progname = *argv++;
 argc--;

 condonly = FALSE;
 conditioned = FALSE;
 
 while (argc && **argv == '-')  {
    switch (*(*argv+1))  {
       case 'c':
	  condonly = TRUE;
	  break;
       case 'r':
	  conditioned = TRUE;
	  break;
       case 'm':
	  minwidth = atoi(*++argv);
	  argc--;
	  break;
       case 'x':
	  dxsize = atoi(*++argv);
	  argc--;
	  break;
       case 'y':
	  dysize = atoi(*++argv);
	  argc--;
	  break;
       case 'o':
	  strcpy(outname, *++argv);
	  outtack = findtack(outname);
	  argc--;
	  break;
       default:
	  error(WARNING, "unknown option");
	  break;
       }
    argv++;
    argc--;
    }

 if (conditioned) {
    if (argc)
       while (argc)  {
	  fp = efopen(*argv, "r");
	  plot(fp);
	  fclose(fp);
	  argv++;
	  argc--;
	  }
    else
       plot(stdin);
    if (lineno)
	nextpage();
 } else  {
    comargs[0] = '\0';
    while (argc)  {
       strcat(comargs, " ");
       strcat(comargs, *argv);
       argv++;
       argc--;
       }
    sprintf(command, XCOM, comargs);
    if (condonly)
       return(system(command));
    else  {
       if ((fp = popen(command, "r")) == NULL)
          error(SYSTEM, "cannot execute input filter");
       plot(fp);
       pclose(fp);
       if (lineno)
	  nextpage();
       }
    }

 return(0);
 }



void
thispage(void)		/* rewind current file */
{
    if (lineno)
	error(USER, "cannot restart page in thispage");
}


void
initfile(void)		/* initialize this file */
{
    /*
    static unsigned char  cmap[24] = {255,255,255, 255,152,0, 0,188,0, 0,0,255,
			179,179,0, 255,0,255, 0,200,200, 0,0,0};
     */
    static const unsigned char  cmap[8][3] = {{0,0,0}, {0,0,255}, {0,188,0},
		{255,152,0}, {0,200,200}, {255,0,255}, {179,179,0}, {255,255,255}};
    static int  filenum = 0;
    BMPHeader  *hp;
    register int  i;

    hp = BMPmappedHeader(dxsize, dysize, 0, 256);
    hp->compr = BI_RLE8;
    if (hp == NULL)
	error(USER, "Illegal image parameter(s)");
    for (i = 0; i < 8; i++) {
	hp->palette[i].r = cmap[i][2];
	hp->palette[i].g = cmap[i][1];
	hp->palette[i].b = cmap[i][0];
    }
    hp->impColors = 8;
    if (outtack != NULL) {
	sprintf(outtack, "%d.bmp", ++filenum);
	bmpw = BMPopenOutputFile(outname, hp);
    } else {
	bmpw = BMPopenOutputStream(stdout, hp);
    }
    if (bmpw == NULL)
	error(USER, "Cannot create output");
}



void
nextpage(void)		/* advance to next page */

{

    if (lineno == 0)
	return;
    if (bmpw != NULL) {
	while (lineno < dysize) {
	    nextblock();
	    outputblock();
	}
	BMPcloseOutput(bmpw);
	bmpw = NULL;
    }
    lineno = 0;

}



#define MINRUN	4

extern void
printblock(void)		/* output scanline block to file */

{
    int  i;
    register unsigned char  *scanline;

    if (lineno == 0)
	initfile();
    for (i = outblock.ybot; i <= outblock.ytop && i < dysize; i++) {
	scanline = outblock.cols[i-outblock.ybot] + outblock.xleft;
	memcpy((void *)bmpw->scanline, (void *)scanline,
			sizeof(uint8)*(outblock.xright-outblock.xleft+1));
	if (BMPwriteScanline(bmpw) != BIR_OK)
		error(SYSTEM, "Error writing BMP file");
	lineno++;
    }
    
}
