#ifndef lint
static const char	RCSid[] = "$Id: ra_gif.c,v 2.9 2003/07/21 22:30:19 schorsch Exp $";
#endif
/*
 * Convert from Radiance picture file to Compuserve GIF.
 * Currently, we don't know how to get back.
 */

#include  <stdio.h>
#include  <time.h>
#include  <math.h>

#include  "platform.h"
#include  "color.h"
#include  "resolu.h"

#define MAXCOLORS		256

int rmap[MAXCOLORS];
int gmap[MAXCOLORS];
int bmap[MAXCOLORS];

int currow;

extern long  ftell();

long  picstart;

BYTE  clrtab[256][3];

extern int  samplefac;

extern int  getgifpix();

COLR	*scanln;
BYTE	*pixscan;

int  xmax, ymax;			/* picture size */

double	gamv = 2.2;			/* gamma correction */

int  greyscale = 0;			/* convert to B&W? */

int  dither = 1;			/* dither colors? */

int  bradj = 0;				/* brightness adjustment */

int  ncolors = 0;			/* number of colors requested */

char  *progname;


main(argc, argv)
int  argc;
char  *argv[];
{
	int  bitsperpix;
	int  i;
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	progname = argv[0];
	samplefac = 0;

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'g':
				gamv = atof(argv[++i]);
				break;
			case 'b':
				greyscale = 1;
				break;
			case 'd':
				dither = !dither;
				break;
			case 'c':
				ncolors = atoi(argv[++i]);
				break;
			case 'e':
				if (argv[i+1][0] != '+' && argv[i+1][0] != '-')
					goto userr;
				bradj = atoi(argv[++i]);
				break;
			case 'n':
				samplefac = atoi(argv[++i]);
				break;
			default:
				goto userr;
			}
		else
			break;

	if (i < argc-2)
		goto userr;
	if (i <= argc-1 && freopen(argv[i], "r", stdin) == NULL) {
		fprintf(stderr, "%s: cannot open input \"%s\"\n",
				progname, argv[i]);
		exit(1);
	}
	if (i == argc-2 && freopen(argv[i+1], "w", stdout) == NULL) {
		fprintf(stderr, "%s: cannot open output \"%s\"\n",
				progname, argv[i+1]);
		exit(1);
	}
	if (checkheader(stdin, COLRFMT, NULL) < 0 ||
			fgetresolu(&xmax, &ymax, stdin) < 0) {
		fprintf(stderr, "%s: bad picture format\n", progname);
		exit(1);
	}
	picstart = ftell(stdin);
	currow = -1;
	scanln = (COLR *)malloc(xmax*sizeof(COLR));
	if (scanln == NULL) {
		fprintf(stderr, "%s: out of memory\n", progname);
		exit(1);
	}
				/* set up gamma correction */
	setcolrgam(gamv);
				/* figure out the bits per pixel */
	if (ncolors < 2 | ncolors > MAXCOLORS)
		ncolors = MAXCOLORS;
	for (bitsperpix = 1; ncolors > 1<<bitsperpix; bitsperpix++)
		;
				/* compute color map */
	if (greyscale)
		mkgrymap(ncolors);
	else
		mkclrmap(ncolors);
				/* convert image */
	GIFEncode(stdout, xmax, ymax, 0, 0, bitsperpix,
			rmap, gmap, bmap, getgifpix);
	exit(0);
userr:
	fprintf(stderr,
	"Usage: %s [-b][-d][-n samp][-c ncolors][-g gamv][-e +/-stops] input [output]\n",
			progname);
	exit(1);
}


getrow(y)			/* get the specified row from our image */
int  y;
{
	if (y == currow)
		return;
	if (y < currow) {
		fseek(stdin, picstart, 0);
		currow = -1;
	}
	do
		if (freadcolrs(scanln, xmax, stdin) < 0) {
			fprintf(stderr, "%s: error reading picture (y==%d)\n",
					progname, ymax-1-y);
			exit(1);
		}
	while (++currow < y);
	if (bradj)
		shiftcolrs(scanln, xmax, bradj);
	colrs_gambs(scanln, xmax);
	if (pixscan != NULL) {
		if (samplefac)
			neu_dith_colrs(pixscan, scanln, xmax);
		else
			dith_colrs(pixscan, scanln, xmax);
	}
}


mkclrmap(nc)			/* make our color map */
int	nc;
{
	register int	i;

	if ((samplefac ? neu_init(xmax*ymax) : new_histo(xmax*ymax)) == -1)
		goto memerr;
	for (i = 0; i < ymax; i++) {
		getrow(i);
		if (samplefac)
			neu_colrs(scanln, xmax);
		else
			cnt_colrs(scanln, xmax);
	}
	if (samplefac)
		neu_clrtab(nc);
	else
		new_clrtab(nc);
	for (i = 0; i < nc; i++) {
		rmap[i] = clrtab[i][RED];
		gmap[i] = clrtab[i][GRN];
		bmap[i] = clrtab[i][BLU];
	}
	if (dither && (pixscan = (BYTE *)malloc(xmax)) == NULL)
		goto memerr;
	return;
memerr:
	fprintf(stderr, "%s: out of memory\n", progname);
	exit(1);
}


mkgrymap(nc)
int	nc;
{
	register int	i;

	for (i = 0; i < nc; i++) {
		rmap[i] =
		gmap[i] =
		bmap[i] = ((unsigned)i<<8)/nc;
	}
}


int
getgifpix(x, y)			/* get a single pixel from our picture */
int  x, y;
{
	int pix;

	getrow(y);
	if (greyscale)
		return((normbright(scanln[x])*ncolors)>>8);
	if (pixscan != NULL)
		return(pixscan[x]);
	return(samplefac ? neu_map_pixel(scanln[x]) : map_pixel(scanln[x]));
}


/*
 * SCARY GIF code follows . . . . sorry.
 *
 * Based on GIFENCOD by David Rowley <mgardi@watdscu.waterloo.edu>.A
 * Lempel-Zim compression based on "compress".
 *
 */

/*****************************************************************************
 *
 * GIFENCODE.C    - GIF Image compression interface
 *
 * GIFEncode( FName, GHeight, GWidth, GInterlace, Background,
 *            BitsPerPixel, Red, Green, Blue, GetPixel )
 *
 *****************************************************************************/
typedef int (* ifunptr)();

#define TRUE 1
#define FALSE 0

int Width, Height;
int curx, cury;
long CountDown;
int Pass;
int Interlace;
unsigned long cur_accum = 0;
int cur_bits = 0;

/*
 * Bump the 'curx' and 'cury' to point to the next pixel
 */
BumpPixel()
{
    curx++;
    if( curx == Width ) {
	curx = 0;
	if( !Interlace ) {
	    cury++;
	} else {
	    switch( Pass ) {
	     	case 0:
		    cury += 8;
		    if( cury >= Height ) {
			Pass++;
			cury = 4;
		    } 
		    break;
		case 1:
		    cury += 8;
		        if( cury >= Height ) {
			    Pass++;
			    cury = 2;
		      	}
		      	break;
	       	case 2:
		    cury += 4;
		    if( cury >= Height ) {
		     	Pass++;
			cury = 1;
		    }
		    break;
	       	case 3:
		    cury += 2;
		    break;
	    }
	}
    }
}

/*
 * Return the next pixel from the image
 */
GIFNextPixel( getpixel )
ifunptr getpixel;
{
    int r;

    if( CountDown == 0 )
	return EOF;
    CountDown--;
    r = (*getpixel)( curx, cury );
    BumpPixel();
    return r;
}

/*
 * public GIFEncode
 */
GIFEncode( fp, GWidth, GHeight, GInterlace, Background,
           BitsPerPixel, Red, Green, Blue, GetPixel )
FILE *fp;
int GWidth, GHeight;
int GInterlace;
int Background;
int BitsPerPixel;
int Red[], Green[], Blue[];
ifunptr GetPixel;
{
    int B;
    int RWidth, RHeight;
    int LeftOfs, TopOfs;
    int Resolution;
    int ColorMapSize;
    int InitCodeSize;
    int i;

    long cur_accum = 0;
    cur_bits = 0;

    Interlace = GInterlace;
    ColorMapSize = 1 << BitsPerPixel;
    RWidth = Width = GWidth;
    RHeight = Height = GHeight;
    LeftOfs = TopOfs = 0;
    Resolution = BitsPerPixel;

    CountDown = (long)Width * (long)Height;
    Pass = 0;
    if( BitsPerPixel <= 1 )
	InitCodeSize = 2;
    else
	InitCodeSize = BitsPerPixel;
    curx = cury = 0;
    fwrite( "GIF87a", 1, 6, fp );
    Putword( RWidth, fp );
    Putword( RHeight, fp );
    B = 0x80;       /* Yes, there is a color map */
    B |= (Resolution - 1) << 5;
    B |= (BitsPerPixel - 1);
    fputc( B, fp );
    fputc( Background, fp );
    fputc( 0, fp );
    for( i=0; i<ColorMapSize; i++ ) {
	fputc( Red[i], fp );
	fputc( Green[i], fp );
	fputc( Blue[i], fp );
    }
    fputc( ',', fp );
    Putword( LeftOfs, fp );
    Putword( TopOfs, fp );
    Putword( Width, fp );
    Putword( Height, fp );
    if( Interlace )
	fputc( 0x40, fp );
    else
	fputc( 0x00, fp );
    fputc( InitCodeSize, fp );
    compress( InitCodeSize+1, fp, GetPixel );
    fputc( 0, fp );
    fputc( ';', fp );
    fclose( fp );
}

/*
 * Write out a word to the GIF file
 */
Putword( w, fp )
int w;
FILE *fp;
{
    fputc( w & 0xff, fp );
    fputc( (w/256) & 0xff, fp );
}


/***************************************************************************
 *
 *  GIFCOMPR.C       - GIF Image compression routines
 *
 *  Lempel-Ziv compression based on 'compress'.  GIF modifications by
 *  David Rowley (mgardi@watdcsu.waterloo.edu)
 *
 ***************************************************************************/

#define CBITS    12
#define HSIZE  5003            /* 80% occupancy */

/*
 * a code_int must be able to hold 2**CBITS values of type int, and also -1
 */
typedef int             code_int;
typedef long int        count_int;
typedef unsigned char  	char_type;

/*
 *
 * GIF Image compression - modified 'compress'
 *
 * Based on: compress.c - File compression ala IEEE Computer, June 1984.
 *
 * By Authors:  Spencer W. Thomas       (decvax!harpo!utah-cs!utah-gr!thomas)
 *              Jim McKie               (decvax!mcvax!jim)
 *              Steve Davies            (decvax!vax135!petsd!peora!srd)
 *              Ken Turkowski           (decvax!decwrl!turtlevax!ken)
 *              James A. Woods          (decvax!ihnp4!ames!jaw)
 *              Joe Orost               (decvax!vax135!petsd!joe)
 *
 */
#include <ctype.h>

#define ARGVAL() (*++(*argv) || (--argc && *++argv))

int n_bits;                        /* number of bits/code */
int maxbits = CBITS;                /* user settable max # bits/code */
code_int maxcode;                  /* maximum code, given n_bits */
code_int maxmaxcode = (code_int)1 << CBITS; /* should NEVER generate this code */
# define MAXCODE(n_bits)        (((code_int) 1 << (n_bits)) - 1)

count_int htab [HSIZE];
unsigned short codetab [HSIZE];
#define HashTabOf(i)       htab[i]
#define CodeTabOf(i)    codetab[i]

code_int hsize = HSIZE;                 /* for dynamic table sizing */

/*
 * To save much memory, we overlay the table used by compress() with those
 * used by decompress().  The tab_prefix table is the same size and type
 * as the codetab.  The tab_suffix table needs 2**CBITS characters.  We
 * get this from the beginning of htab.  The output stack uses the rest
 * of htab, and contains characters.  There is plenty of room for any
 * possible stack (stack used to be 8000 characters).
 */
#define tab_prefixof(i) CodeTabOf(i)
#define tab_suffixof(i)        ((char_type *)(htab))[i]
#define de_stack               ((char_type *)&tab_suffixof((code_int)1<<CBITS))

code_int free_ent = 0;                  /* first unused entry */

/*
 * block compression parameters -- after all codes are used up,
 * and compression rate changes, start over.
 */
int clear_flg = 0;
int offset;
long int in_count = 1;            /* length of input */
long int out_count = 0;           /* # of codes output (for debugging) */

/*
 * compress stdin to stdout
 *
 * Algorithm:  use open addressing double hashing (no chaining) on the
 * prefix code / next character combination.  We do a variant of Knuth's
 * algorithm D (vol. 3, sec. 6.4) along with G. Knott's relatively-prime
 * secondary probe.  Here, the modular division first probe is gives way
 * to a faster exclusive-or manipulation.  Also do block compression with
 * an adaptive reset, whereby the code table is cleared when the compression
 * ratio decreases, but after the table fills.  The variable-length output
 * codes are re-sized at this point, and a special CLEAR code is generated
 * for the decompressor.  Late addition:  construct the table according to
 * file size for noticeable speed improvement on small files.  Please direct
 * questions about this implementation to ames!jaw.
 */

int g_init_bits;
FILE *g_outfile;
int ClearCode;
int EOFCode;

compress( init_bits, outfile, ReadValue )
int init_bits;
FILE *outfile;
ifunptr ReadValue;
{
    register long fcode;
    register code_int i = 0;
    register int c;
    register code_int ent;
    register code_int disp;
    register code_int hsize_reg;
    register int hshift;

    /*
     * Set up the globals:  g_init_bits - initial number of bits
     *                      g_outfile   - pointer to output file
     */
    g_init_bits = init_bits;
    g_outfile = outfile;
    /*
     * Set up the necessary values
     */
    offset = 0;
    out_count = 0;
    clear_flg = 0;
    in_count = 1;
    maxcode = MAXCODE(n_bits = g_init_bits);
    ClearCode = (1 << (init_bits - 1));
    EOFCode = ClearCode + 1;
    free_ent = ClearCode + 2;
    char_init();
    ent = GIFNextPixel( ReadValue );
    hshift = 0;
    for ( fcode = (long) hsize;  fcode < 65536L; fcode *= 2L )
        hshift++;
    hshift = 8 - hshift;                /* set hash code range bound */
    hsize_reg = hsize;
    cl_hash( (count_int) hsize_reg);            /* clear hash table */
    output( (code_int)ClearCode );
    while ( (c = GIFNextPixel( ReadValue )) != EOF ) {
        in_count++;
        fcode = (long) (((long) c << maxbits) + ent);
        /* i = (((code_int)c << hshift) ~ ent);    /* xor hashing */
        i = (((code_int)c << hshift) ^ ent);    /* xor hashing */
        if ( HashTabOf (i) == fcode ) {
            ent = CodeTabOf (i);
            continue;
        } else if ( (long)HashTabOf (i) < 0 )      /* empty slot */
            goto nomatch;
        disp = hsize_reg - i;           /* secondary hash (after G. Knott) */
        if ( i == 0 )
            disp = 1;
probe:
        if ( (i -= disp) < 0 )
            i += hsize_reg;
        if ( HashTabOf (i) == fcode ) {
            ent = CodeTabOf (i);
            continue;
        }
        if ( (long)HashTabOf (i) > 0 )
            goto probe;
nomatch:
        output ( (code_int) ent );
        out_count++;
        ent = c;
        if ( free_ent < maxmaxcode ) {
            CodeTabOf (i) = free_ent++; /* code -> hashtable */
            HashTabOf (i) = fcode;
        } else
                cl_block();
    }
    /*
     * Put out the final code.
     */
    output( (code_int)ent );
    out_count++;
    output( (code_int) EOFCode );
    return;
}

/*****************************************************************
 * TAG( output )
 *
 * Output the given code.
 * Inputs:
 *      code:   A n_bits-bit integer.  If == -1, then EOF.  This assumes
 *              that n_bits =< (long)wordsize - 1.
 * Outputs:
 *      Outputs code to the file.
 * Assumptions:
 *      Chars are 8 bits long.
 * Algorithm:
 *      Maintain a CBITS character long buffer (so that 8 codes will
 * fit in it exactly).  Use the VAX insv instruction to insert each
 * code in turn.  When the buffer fills up empty it and start over.
 */

unsigned long masks[] = { 0x0000, 0x0001, 0x0003, 0x0007, 0x000F,
                                  0x001F, 0x003F, 0x007F, 0x00FF,
                                  0x01FF, 0x03FF, 0x07FF, 0x0FFF,
                                  0x1FFF, 0x3FFF, 0x7FFF, 0xFFFF };

output( code )
code_int  code;
{
    cur_accum &= masks[ cur_bits ];
    if( cur_bits > 0 )
        cur_accum |= ((long)code << cur_bits);
    else
        cur_accum = code;
    cur_bits += n_bits;
    while( cur_bits >= 8 ) {
        char_out( (unsigned int)(cur_accum & 0xff) );
        cur_accum >>= 8;
        cur_bits -= 8;
    }

    /*
     * If the next entry is going to be too big for the code size,
     * then increase it, if possible.
     */
   if ( free_ent > maxcode || clear_flg ) {
            if( clear_flg ) {
                maxcode = MAXCODE (n_bits = g_init_bits);
                clear_flg = 0;
            } else {
                n_bits++;
                if ( n_bits == maxbits )
                    maxcode = maxmaxcode;
                else
                    maxcode = MAXCODE(n_bits);
            }
    }
    if( code == EOFCode ) {
        /*
         * At EOF, write the rest of the buffer.
         */
        while( cur_bits > 0 ) {
                char_out( (unsigned int)(cur_accum & 0xff) );
                cur_accum >>= 8;
                cur_bits -= 8;
        }
        flush_char();
        fflush( g_outfile );
        if( ferror( g_outfile ) )
                writeerr();
    }
}

/*
 * Clear out the hash table
 */
cl_block ()             /* table clear for block compress */
{
        cl_hash ( (count_int) hsize );
        free_ent = ClearCode + 2;
        clear_flg = 1;
        output( (code_int)ClearCode );
}

cl_hash(hsize)          /* reset code table */
register count_int hsize;
{
        register count_int *htab_p = htab+hsize;
        register long i;
        register long m1 = -1;

        i = hsize - 16;
        do {                            /* might use Sys V memset(3) here */
                *(htab_p-16) = m1;
                *(htab_p-15) = m1;
                *(htab_p-14) = m1;
                *(htab_p-13) = m1;
                *(htab_p-12) = m1;
                *(htab_p-11) = m1;
                *(htab_p-10) = m1;
                *(htab_p-9) = m1;
                *(htab_p-8) = m1;
                *(htab_p-7) = m1;
                *(htab_p-6) = m1;
                *(htab_p-5) = m1;
                *(htab_p-4) = m1;
                *(htab_p-3) = m1;
                *(htab_p-2) = m1;
                *(htab_p-1) = m1;
                htab_p -= 16;
        } while ((i -= 16) >= 0);
        for ( i += 16; i > 0; i-- )
                *--htab_p = m1;
}

writeerr()
{
        printf( "error writing output file\n" );
        exit(1);
}

/******************************************************************************
 *
 * GIF Specific routines
 *
 ******************************************************************************/

/*
 * Number of characters so far in this 'packet'
 */
int a_count;

/*
 * Set up the 'byte output' routine
 */
char_init()
{
        a_count = 0;
}

/*
 * Define the storage for the packet accumulator
 */
char accum[256];

/*
 * Add a character to the end of the current packet, and if it is 254
 * characters, flush the packet to disk.
 */
char_out( c )
int c;
{
        accum[ a_count++ ] = c;
        if( a_count >= 254 )
                flush_char();
}

/*
 * Flush the packet to disk, and reset the accumulator
 */
flush_char()
{
        if( a_count > 0 ) {
                fputc( a_count, g_outfile );
                fwrite( accum, 1, a_count, g_outfile );
                a_count = 0;
        }
}

static float curgamma;
static short gamtab[256];

gammawarp(sbuf,gam,n)
short *sbuf;
float gam;
int n;
{
    int i;
    float f;

    if(gam!=curgamma) {
	for(i=0; i<256; i++) 
	    gamtab[i] = 255*pow(i/255.0,gam)+0.5;
	curgamma = gam;
    }
    while(n--) {
	*sbuf = gamtab[*sbuf];
	sbuf++;
    }
}
