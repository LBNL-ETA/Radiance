#ifndef lint
static const char	RCSid[] = "$Id: meta2tga.c,v 1.4 2003/10/27 10:28:59 schorsch Exp $";
#endif
/*
 *  Program to convert meta-files to Targa 8-bit color-mapped format
 */

#include  "copyright.h"

#include  "rtprocess.h"
#include  "meta.h"
#include  "plot.h"
#include  "rast.h"
#include  "targa.h"

#define  MAXALLOC  5000
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

static FILE  *fout;

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


main(argc, argv)

int  argc;
char  **argv;

{
 FILE  *fp;
 char  comargs[200], command[300];

  fout = stdout;
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






thispage()		/* rewind current file */
{
    if (lineno)
	error(USER, "cannot restart page in thispage");
}



initfile()		/* initialize this file */
{
    static int  filenum = 0;
    /*
    static unsigned char  cmap[24] = {255,255,255, 255,152,0, 0,188,0, 0,0,255,
			179,179,0, 255,0,255, 0,200,200, 0,0,0};
     */
    static unsigned char  cmap[24] = {0,0,0, 0,0,255, 0,188,0, 255,152,0,
			0,200,200, 255,0,255, 179,179,0, 255,255,255};
    struct hdStruct  thead;
    register int  i;

    if (outtack != NULL) {
	sprintf(outtack, "%d.tga", ++filenum);
	fout = efopen(outname, "w");
    }
    if (fout == NULL)
	error(USER, "no output file");
    thead.mapType = CM_HASMAP;
    thead.dataType = IM_CCMAP;
    thead.mapOrig = 0;
    thead.mapLength = 256;
    thead.CMapBits = 24;
    thead.XOffset = 0;
    thead.YOffset = 0;
    thead.x = dxsize;
    thead.y = dysize;
    thead.dataBits = 8;
    thead.imType = 0;
    putthead(&thead, NULL, fout);
    for (i = 0; i < 8*3; i++)
	putc(cmap[i], fout);
    while (i++ < 256*3)
	putc(0, fout);
}




nextpage()		/* advance to next page */

{

    if (lineno == 0)
	return;
    if (fout != NULL) {
	while (lineno < dysize) {
	    nextblock();
	    outputblock();
	}
	fclose(fout);
	fout = NULL;
    }
    lineno = 0;

}



#define MINRUN	4


printblock()		/* output scanline block to file */

{
    int  i, c2;
    register unsigned char  *scanline;
    register int  j, beg, cnt;

    if (lineno == 0)
	initfile();
    for (i = outblock.ybot; i <= outblock.ytop && i < dysize; i++) {
	scanline = outblock.cols[i-outblock.ybot];
	for (j = outblock.xleft; j <= outblock.xright; j += cnt) {
	    for (beg = j; beg <= outblock.xright; beg += cnt) {
		for (cnt = 1; cnt < 128 && beg+cnt <= outblock.xright &&
			scanline[beg+cnt] == scanline[beg]; cnt++)
		    ;
		if (cnt >= MINRUN)
		    break;			/* long enough */
	    }
	    while (j < beg) {		/* write out non-run */
		if ((c2 = beg-j) > 128) c2 = 128;
		putc(c2-1, fout);
		while (c2--)
		    putc(scanline[j++], fout);
	    }
	    if (cnt >= MINRUN) {		/* write out run */
		putc(127+cnt, fout);
		putc(scanline[beg], fout);
	    } else
		cnt = 0;
	}
	lineno++;
    }
    
}


putint2(i, fp)			/* put a 2-byte positive integer */
register int  i;
register FILE	*fp;
{
	putc(i&0xff, fp);
	putc(i>>8&0xff, fp);
}


putthead(hp, ip, fp)		/* write header to output */
struct hdStruct  *hp;
char  *ip;
register FILE  *fp;
{
	if (ip != NULL)
		putc(strlen(ip), fp);
	else
		putc(0, fp);
	putc(hp->mapType, fp);
	putc(hp->dataType, fp);
	putint2(hp->mapOrig, fp);
	putint2(hp->mapLength, fp);
	putc(hp->CMapBits, fp);
	putint2(hp->XOffset, fp);
	putint2(hp->YOffset, fp);
	putint2(hp->x, fp);
	putint2(hp->y, fp);
	putc(hp->dataBits, fp);
	putc(hp->imType, fp);

	if (ip != NULL)
		fputs(ip, fp);

	return(ferror(fp) ? -1 : 0);
}
