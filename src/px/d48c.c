#ifndef lint
static const char	RCSid[] = "$Id: d48c.c,v 1.3 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 *  d48c.c - driver for dicomed D48 film recorder w/ color optics.
 *
 *     3/26/87
 */

#include  <stdio.h>

#include  <signal.h>

#include  "color.h"


#define  HSIZE		4096		/* device width */
#define  VSIZE		4096		/* device height */

#define  MAXMULT	15		/* maximum element size */

				/* dicomed commands */
#define  O_ECS		0		/* encoded command select */
#define  O_CONF		03400		/* configuration code */
#define  O_FES		020000		/* function element select */
#define  O_VPA		040000		/* vector or position absolute */
#define  O_ICS		0100000		/* initial condition select */
#define  O_PES		0120000		/* point element select */

				/* ECS commands */
#define  ECS_FA		07		/* film advance */
#define  ECS_RED	011		/* select red filter */
#define  ECS_GRN	012		/* select green filter */
#define  ECS_BLU	014		/* select blue filter */
#define  ECS_NEU	017		/* select neutral filter */
#define  ECS_RAS	020		/* select raster mode */
#define  ECS_COL	022		/* select color mode */
#define  ECS_BW		023		/* select black and white mode */
#define  ECS_BT		024		/* bypass translation tables */

				/* FES functions */
#define  FES_RE		030000		/* raster element */
#define  FES_RL		0170000		/* run length */

				/* ICS fields */
#define  ICS_HVINT	040		/* interchange h and v axes */
#define  ICS_VREV	0100		/* reverse v axis direction */
#define  ICS_HREV	0200		/* reverse h axis direction */

int  nrows;				/* number of rows for output */
int  ncols;				/* number of columns for output */

int  mult;				/* number of h/v points per elem. */

int  hhome;				/* offset for h */
int  vhome;				/* offset for v */

				/* dicomed configurations */
#define  C_POLAROID	0		/* polaroid print */
#define  C_35FRES	1		/* 35mm full resolution */
#define  C_35FFRA	2		/* 35mm full frame */
#define  C_SHEET	3		/* 4X5 sheet film */

#define  NCONF		4		/* number of configurations */

struct  config {
	char  *name;			/* configuration name */
	int  ccode;			/* configuration code */
	int  icsfield;			/* initial condition select */
	int  hbeg, vbeg;		/* coord. origin */
	int  hsiz, vsiz;		/* dimensions */
}  ctab[NCONF] = {
	{ "polaroid", 11, ICS_HVINT|ICS_VREV, 0, 0, 4096, 4096 },
	{ "35r", 5, 0, 0, 0, 4096, 4096 },
	{ "35f", 7, 0, 39, 692, 4018, 2712 },
	{ "sheet", 10, ICS_HVINT|ICS_VREV, 0, 0, 4096, 4096 },
};

#define  NEU		3		/* neutral color */

struct  filter {
	float  cv[3];			/* filter color vector */
	char  *fn;			/* file name for this color */
	FILE  *fp;			/* stream for this color */
}  ftab[4] = {
	{ {1.0, 0.0, 0.0}, NULL, NULL },	/* red */
	{ {0.0, 1.0, 0.0}, NULL, NULL },	/* green */
	{ {0.0, 0.0, 1.0}, NULL, NULL },	/* blue */
	{ {0.3, 0.59, 0.11}, NULL, NULL },	/* neutral */
};

#define  f_val(col,f)	( colval(col,0)*ftab[f].cv[0] + \
			  colval(col,1)*ftab[f].cv[1] + \
			  colval(col,2)*ftab[f].cv[2] )

float  neumap[] = {			/* neutral map, log w/ E6 */
	.0,		.0,
	.0004,		.059,
	.0008,		.122,
	.0018,		.184,
	.0033,		.247,
	.0061,		.310,
	.0119,		.373,
	.0202,		.435,
	.0373,		.498,
	.065,		.561,
	.113,		.624,
	.199,		.686,
	.300,		.749,
	.432,		.812,
	.596,		.875,
	.767,		.937,
	1.000, 		1.000,
};


#define  TEMPLATE	"/tmp/dXXXXXX"

int  docolor = 1;		/* color? */

int  cftype = C_35FFRA;		/* configuration type */

char  *progname;		/* program name */


main(argc, argv)
int  argc;
char  *argv[];
{
	int  quit();
	int  i;

	if (signal(SIGINT, quit) == SIG_IGN)
		signal(SIGINT, SIG_IGN);
	if (signal(SIGHUP, quit) == SIG_IGN)
		signal(SIGINT, SIG_IGN);
	signal(SIGTERM, quit);
	signal(SIGPIPE, quit);
	signal(SIGXCPU, quit);
	signal(SIGXFSZ, quit);

	progname = argv[0];

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'c':			/* color */
			docolor = 1;
			break;
		case 'b':			/* black and white */
			docolor = 0;
			break;
		case 't':			/* type */
			settype(argv[++i]);
			break;
		default:
			fprintf(stderr, "%s: Unknown option: %s\n",
					progname, argv[i]);
			quit(1);
			break;
		}
			
	if (i < argc)
		for ( ; i < argc; i++)
			dofile(argv[i]);
	else
		dofile(NULL);

	quit(0);
}


settype(cname)			/* set configuration type */
char  *cname;
{
	for (cftype = 0; cftype < NCONF; cftype++)
		if (!strcmp(cname, ctab[cftype].name))
			return;
	
	fprintf(stderr, "%s: Unknown type: %s\n", progname, cname);
	quit(1);
}


void
quit(code)			/* quit program */
int  code;
{
	int  i;

	for (i = 0; i < 4; i++)
		if (ftab[i].fn != NULL)
			unlink(ftab[i].fn);
	exit(code);
}


dofile(fname)			/* convert file to dicomed format */
char  *fname;
{
	FILE  *fin;
	char  sbuf[128];
	COLOR  scanline[HSIZE];
	int  i;

	if (fname == NULL) {
		fin = stdin;
		fname = "<stdin>";
	} else if ((fin = fopen(fname, "r")) == NULL) {
		fprintf(stderr, "%s: Cannot open: %s\n", progname, fname);
		quit(1);
	}
				/* discard header */
	while (fgets(sbuf, sizeof(sbuf), fin) != NULL && sbuf[0] != '\n')
		;
				/* get picture size */
	if (fgets(sbuf, sizeof(sbuf), fin) == NULL ||
			sscanf(sbuf, "-Y %d +X %d\n", &nrows, &ncols) != 2) {
		fprintf(stderr, "%s: Bad picture size: %s\n", progname, fname);
		quit(1);
	}

	mult = ctab[cftype].hsiz / ncols;
	if (ctab[cftype].vsiz / nrows < mult)
		mult = ctab[cftype].vsiz / nrows;

	if (mult < 1) {
		fprintf(stderr, "%s: Resolution mismatch: %s\n",
				progname, fname);
		quit(1);
	}
	if (mult > MAXMULT)
		mult = MAXMULT;		/* maximum element size */

	hhome = ctab[cftype].hbeg + (ctab[cftype].hsiz - ncols*mult)/2;
	vhome = ctab[cftype].vbeg + (ctab[cftype].vsiz - nrows*mult)/2;

	d_init();

	for (i = nrows-1; i >= 0; i--) {
		if (freadscan(scanline, ncols, fin) < 0) {
			fprintf(stderr, "%s: Read error in row %d: %s\n",
					progname, i, fname);
			quit(1);
		}
		scanout(scanline, i);
	}

	d_done();

	fclose(fin);
}


d_init()			/* set up dicomed, files */
{
	char  *mktemp();
					/* initial condition select */
	putw(O_ICS|ctab[cftype].icsfield, stdout);
					/* configuration code */
	putw(O_CONF|ctab[cftype].ccode, stdout);
					/* select raster mode */
	putw(O_ECS|ECS_RAS, stdout);
					/* color or black and white */
	putw(O_ECS|(docolor?ECS_COL:ECS_BW), stdout);
					/* bypass translation tables */
	putw(O_ECS|ECS_BT, stdout);
					/* point element select */
	putw(O_PES|mult<<9|010, stdout);
	putw(mult<<9|010, stdout);
	putw(mult<<4|mult, stdout);

	if (docolor) {			/* color output */
						/* set up red file */
		if (ftab[RED].fn == NULL)
			ftab[RED].fn = mktemp(TEMPLATE);
		if ((ftab[RED].fp = fopen(ftab[RED].fn, "w+")) == NULL) {
			fprintf(stderr, "%s: Cannot open: %s\n",
					progname, ftab[RED].fn);
			quit(1);
		}
		putw(O_ECS|ECS_RED, ftab[RED].fp);	/* red filter */
		putw(O_VPA, ftab[RED].fp);		/* home */
		putw(hhome<<3, ftab[RED].fp);
		putw(vhome<<3, ftab[RED].fp);
						/* set up green file */
		ftab[GRN].fp = stdout;
		putw(O_ECS|ECS_GRN, ftab[GRN].fp);	/* green filter */
		putw(O_VPA, ftab[GRN].fp);		/* home */
		putw(hhome<<3, ftab[GRN].fp);
		putw(vhome<<3, ftab[GRN].fp);
						/* set up blue file */
		if (ftab[BLU].fn == NULL)
			ftab[BLU].fn = mktemp(TEMPLATE);
		if ((ftab[BLU].fp = fopen(ftab[BLU].fn, "w+")) == NULL) {
			fprintf(stderr, "%s: Cannot open: %s\n",
					progname, ftab[BLU].fn);
			quit(1);
		}
		putw(O_ECS|ECS_BLU, ftab[BLU].fp);	/* blue filter */
		putw(O_VPA, ftab[BLU].fp);		/* home */
		putw(hhome<<3, ftab[BLU].fp);
		putw(vhome<<3, ftab[BLU].fp);

	} else {			/* black and white */
		ftab[NEU].fp = stdout;
		putw(O_ECS|ECS_NEU, ftab[NEU].fp);	/* neutral filter */
		putw(O_VPA, ftab[NEU].fp);		/* home */
		putw(hhome<<3, ftab[NEU].fp);
		putw(vhome<<3, ftab[NEU].fp);
	}
}


scanout(scan, y)			/* output scan line */
COLOR  *scan;
int  y;
{
	int  i;
						/* uses horiz. flyback */
	if (docolor)
		for (i = 0; i < 3; i++)
			runlength(scan, i);
	else
		runlength(scan, NEU);
}


d_done()				/* flush files, finish frame */
{
	if (docolor) {
		transfer(ftab[RED].fp, stdout);
		transfer(ftab[BLU].fp, stdout);
	}
	putw(O_ECS|ECS_FA, stdout);		/* film advance */
	fflush(stdout);
	if (ferror(stdout)) {
		fprintf(stderr, "%s: Write error: <stdout>\n", progname);
		quit(1);
	}
}


runlength(scan, f)			/* do scanline for filter f */
COLOR  *scan;
int  f;
{
	double  mapfilter();
	BYTE  ebuf[2*HSIZE];
	register int  j;
	register BYTE  *ep;

	ep = ebuf;
	ep[0] = mapfilter(scan[0], f) * 255.9999;
	ep[1] = 1;
	ep += 2;
	for (j = 1; j < ncols; j++) {
		ep[0] = mapfilter(scan[j], f) * 255.9999;
		if (ep[0] == ep[-2] && ep[-1] < 255)
			ep[-1]++;
		else {
			ep[1] = 1;
			ep += 2;
		}
	}
	j = ep - ebuf;
	putw(O_FES, ftab[f].fp);
	putw(FES_RL|j>>1, ftab[f].fp);
	for (ep = ebuf; j-- > 0; ep++)
		putc(*ep, ftab[f].fp);

	if (ferror(ftab[f].fp)) {
		fprintf(stderr, "%s: Write error: %s\n", progname,
				ftab[f].fn==NULL ? "<stdout>" : ftab[f].fn);
		quit(1);
	}
}


double
mapfilter(col, f)			/* map filter value */
COLOR  col;
register int  f;
{
	static float  *mp[4] = {neumap, neumap, neumap, neumap};
	double  x, y;

	y = f_val(col,f);

	if (y >= 1.0)
		return(1.0);

	while (y >= mp[f][2])
		mp[f] += 2;
	while (y < mp[f][0])
		mp[f] -= 2;

	x = (y - mp[f][0]) / (mp[f][2] - mp[f][0]);

	return((1.0-x)*mp[f][1] + x*mp[f][3]);
}


putw(w, fp)			/* output a (short) word */
int  w;
register FILE  *fp;
{
	putc(w&0377, fp);		/* in proper order! */
	putc(w>>8, fp);
}


transfer(fin, fout)		/* transfer fin contents to fout, close fin */
register FILE  *fin, *fout;
{
	register int  c;

	fseek(fin, 0L, 0);			/* rewind input */
	
	while ((c = getc(fin)) != EOF)		/* transfer file */
		putc(c, fout);

	fclose(fin);				/* close input */
}
