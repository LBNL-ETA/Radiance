#ifndef lint
static const char	RCSid[] = "$Id: ra_pixar.c,v 2.4 2004/01/02 10:25:13 schorsch Exp $";
#endif
/* ra_pixar.c */
/*
 * convert from RADIANCE image file to PIXAR image file.  (or vice versa)
 *
 *	David G. Jones	July 1989
 */
/* External functions called:
 *
 * header.c:
 *	getheader(fp,func);
 *	printargs(argc,argv,fp);
 * color.c:
 *	freadscan(scanline,len,fp);
 *	fwritescan(scanline,len,fp);
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include  <time.h>

/* PIXAR */
#include <picio.h>
#include <pixeldef.h>
#include <chad.h>

/* RADIANCE (color.h conflicts with PIXAR includes) ugh!! */
#ifdef undef
#include "color.h"
#else
#define			RED			0
#define			GRN			1
#define			BLU			2
#define			EXP			3
#define			colval(col,pri)		((col)[pri])
#define			setcolor(col,r,g,b)	\
	((col)[RED]=(r),(col)[GRN]=(g),(col)[BLU]=(b))
typedef float  COLOR[3];	/* red, green, blue */
#endif
#include "resolu.h"


char			*ProgramName;
int			global_argc;
char			**global_argv;

/* RADIANCE global variables */
FILE			*radiance_fp		= NULL;
double			radiance_pixel();

usage()
{
   fprintf(stderr,"usage:  %s radiance.pic pixar.pix\n",ProgramName);
   fprintf(stderr,"  or    %s -r pixar.pix radiance.pic\n",ProgramName);
   exit(1);
}


main(argc,argv)
int			argc;
char			*argv[];
{
   char			*infile;
   char			*outfile;
   int			i;
   int			reverse			= 0;

   ProgramName=argv[0];
   if (argc < 3)
      usage();
   infile=argv[argc-2];
   outfile=argv[argc-1];
   for (i=1; i < argc-2 ; ++i)
      if (!strcmp(argv[i],"-r"))
	 reverse=1;
      else
	 usage();

   if (reverse)
   {
      global_argc=argc;
      global_argv=argv;
      pix_ra(infile,outfile);
   }
   else
      ra_pix(infile,outfile);
}

ra_pix(infile,outfile)
char			*infile;
char			*outfile;
{
/* PIXAR variables */
   RGBAPixelType	*pixar_scanline;
   PFILE		*out;
/* RADIANCE variables */
   COLOR		*radiance_scanline;
/* regular variables */
   int			width;
   int			height;
   register int		x;
   register int		y;

/* get width, height from radiance file */
   radiance_getsize(infile,&width,&height);
   radiance_scanline=(COLOR *)malloc(width*sizeof(COLOR));
   pixar_scanline=(RGBAPixelType *)malloc(width*sizeof(RGBAPixelType));
   if (!radiance_scanline || !pixar_scanline)
   {
      fprintf(stderr,"not enough memory?\n");
      perror("malloc");
      exit(1);
   }
   memset(pixar_scanline, '\0', width*sizeof(RGBAPixelType));

   PicSetForce(1);
   PicSetPsize(width,height);
   PicSetTsize(width,height);
   PicSetPformat(PF_RGB);
   PicSetPstorage(PS_12DUMP);
   PicSetPmatting(PM_NONE);
   PicSetOffset(0,0);
   PicSetLabel("converted from RADIANCE format by ra_pixar");
   if (!(out=PicCreat(outfile,0666)))
   {
      fprintf(stderr,"can't open PIXAR image file `%s'\n",outfile);
      perror("open");
      exit(1);
   }

   picPreEncodeScanline(out,0);
   for (y=0; y < height ; ++y)
   {
      radiance_readscanline(radiance_scanline,width);
      for (x=0; x < width ; ++x)
      {
	 SetRGBColor(&pixar_scanline[x],
	    DBL2PXL(radiance_pixel(radiance_scanline[x],RED)),
	    DBL2PXL(radiance_pixel(radiance_scanline[x],GRN)),
	    DBL2PXL(radiance_pixel(radiance_scanline[x],BLU)));
      }
      picEncodeScanline(out,pixar_scanline);
   }
   picPostEncodeScanline(out);
   fclose(radiance_fp);
   PicClose(out);
}


radiance_getsize(filename,w,h)
char		*filename;
int		*w;
int		*h;
{
   if (!(radiance_fp=fopen(filename,"r")))
   {
      fprintf(stderr,"can't open RADIANCE image file `%s'\n",filename);
      perror("open");
      exit(1);
   }
   getheader(radiance_fp,NULL,NULL);
   if (fgetresolu(w, h, radiance_fp) < 0)
   {
      fprintf(stderr,"bad RADIANCE format\n");
      exit(1);
   }
}


radiance_readscanline(buf,x)
COLOR			*buf;
int			x;
{
   if (freadscan(buf,x,radiance_fp) < 0)
   {
      fprintf(stderr,"read error?\n");
      perror("fread");
      exit(1);
   }
}


double radiance_pixel(pixel,i)
COLOR			pixel;
int			i;
{
   double		value;

   value=colval(pixel,i);
   if (value < 0.0)
      return 0.0;
   else if (value> 1.0)
      return 1.0;
   return value;
}

pix_ra(infile,outfile)
char			*infile;
char			*outfile;
{
/* PIXAR variables */
   RGBAPixelType	*pixar_scanline;
   PFILE		*in;
/* RADIANCE variables */
   COLOR		*radiance_scanline;
   FILE			*out;
/* regular variables */
   int			width;
   int			height;
   register int		x;
   register int		y;

   if (!(in=PicOpen(infile,"r")))
   {
      fprintf(stderr,"can't open PIXAR image file `%s'\n",infile);
      perror("open");
      exit(1);
   }
   if (!(out=fopen(outfile,"w")))
   {
      fprintf(stderr,"can't open RADIANCE image file `%s'\n",outfile);
      perror("open");
      exit(1);
   }

/* get width, height from PIXAR file */
   width=in->Pwidth;
   height=in->Pheight;

/* allocate scan line space */
   radiance_scanline=(COLOR *)malloc(width*sizeof(COLOR));
   pixar_scanline=(RGBAPixelType *)malloc(width*sizeof(RGBAPixelType));
   if (!radiance_scanline || !pixar_scanline)
   {
      fprintf(stderr,"not enough memory?\n");
      perror("malloc");
      exit(1);
   }

/* write out the RADIANCE header */
   radiance_writeheader(out,width,height);

   picPreDecodeScanline(in,0);
   for (y=0; y < height ; ++y)
   {
      picDecodeScanline(in,pixar_scanline);
      for (x=0; x < width ; ++x)
      {
	 setcolor(radiance_scanline[x],
	    PXL2DBL(pixar_scanline[x].Red),
	    PXL2DBL(pixar_scanline[x].Green),
	    PXL2DBL(pixar_scanline[x].Blue));
      }
      radiance_writescanline(out,radiance_scanline,width);
   }
   picPostDecodeScanline(in);
   PicClose(in);
   fclose(out);
}


radiance_writeheader(fp,x,y)
FILE		*fp;
int		x;
int		y;
{
   printargs(global_argc,global_argv,fp);
   fputc('\n',fp);
   fprtresolu(x, y, fp);
   fflush(fp);
}


radiance_writescanline(fp,buf,x)
FILE			*fp;
COLOR			*buf;
int			x;
{
   if (fwritescan(buf,x,fp) < 0)
   {
      fprintf(stderr,"write error?\n");
      perror("fwrite");
      exit(1);
   }
}
