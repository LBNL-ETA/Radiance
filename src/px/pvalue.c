#ifndef lint
static const char RCSid[] = "$Id: pvalue.c,v 2.19 2003/06/08 12:03:10 schorsch Exp $";
#endif
/*
 *  pvalue.c - program to print pixel values.
 *
 *     4/23/86
 */

#include  <time.h>

#include  "standard.h"
#include  "platform.h"
#include  "color.h"
#include  "resolu.h"

typedef	unsigned short uint16;	/* sizeof (uint16) must == 2 */

#define	 min(a,b)		((a)<(b)?(a):(b))

				/* what to put out (also RED, GRN, BLU) */
#define  ALL		3
#define  BRIGHT		4

RESOLU	picres;			/* resolution of picture */

int  uniq = 0;			/* print only unique values? */

int  doexposure = 0;		/* exposure change? (>100 to print) */

int  dataonly = 0;		/* data only format? */

int  putprim = ALL;		/* what to put out */

int  reverse = 0;		/* reverse conversion? */

int  format = 'a';		/* input/output format */
char  *fmtid = "ascii";		/* format identifier for header */

int  header = 1;		/* do header? */

long  skipbytes = 0;		/* skip bytes in input? */

int  swapbytes = 0;		/* swap bytes in 16-bit words? */

int  interleave = 1;		/* file is interleaved? */

int  resolution = 1;		/* put/get resolution string? */

int  original = 0;		/* convert to original values? */

int  wrongformat = 0;		/* wrong input format? */

double	gamcor = 1.0;		/* gamma correction */

int  ord[3] = {RED, GRN, BLU};	/* RGB ordering */
int  rord[4];			/* reverse ordering */

COLOR  exposure = WHTCOLOR;

char  *progname;

FILE  *fin;
FILE  *fin2 = NULL, *fin3 = NULL;	/* for other color channels */

int  (*getval)(), (*putval)();

double
rgb_bright(clr)
COLOR  clr;
{
	return(bright(clr));
}

double
xyz_bright(clr)
COLOR  clr;
{
	return(clr[CIEY]);
}

double  (*mybright)() = &rgb_bright;


main(argc, argv)
int  argc;
char  **argv;
{
	extern int  checkhead();
	extern long  atol();
	double  d, expval = 1.0;
	int  i;

	progname = argv[0];

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-' || argv[i][0] == '+')
			switch (argv[i][1]) {
			case 'h':		/* header */
				header = argv[i][0] == '+';
				break;
			case 'H':		/* resolution string */
				resolution = argv[i][0] == '+';
				break;
			case 's':		/* skip bytes in header */
				skipbytes = atol(argv[++i]);
				break;
			case 'u':		/* unique values */
				uniq = argv[i][0] == '-';
				break;
			case 'o':		/* original values */
				original = argv[i][0] == '-';
				break;
			case 'g':		/* gamma correction */
				gamcor = atof(argv[i+1]);
				if (argv[i][0] == '+')
					gamcor = 1.0/gamcor;
				i++;
				break;
			case 'e':		/* exposure correction */
				d = atof(argv[i+1]);
				if (argv[i+1][0] == '-' || argv[i+1][0] == '+')
					d = pow(2.0, d);
				if (argv[i][0] == '-')
					expval *= d;
				scalecolor(exposure, d);
				doexposure++;
				i++;
				break;
			case 'R':		/* reverse byte sequence */
				if (argv[i][0] == '-') {
					ord[0]=BLU; ord[1]=GRN; ord[2]=RED;
				} else {
					ord[0]=RED; ord[1]=GRN; ord[2]=BLU;
				}
				break;
			case 'r':		/* reverse conversion */
				reverse = argv[i][0] == '-';
				break;
			case 'n':		/* non-interleaved RGB */
				interleave = argv[i][0] == '+';
				break;
			case 'b':		/* brightness values */
				putprim = argv[i][0] == '-' ? BRIGHT : ALL;
				break;
			case 'p':		/* put primary */
				switch (argv[i][2]) {
				case 'r': case 'R': putprim = RED; break;
				case 'g': case 'G': putprim = GRN; break;
				case 'b': case 'B': putprim = BLU; break;
				default: goto unkopt;
				}
				break;
			case 'd':		/* data only (no indices) */
				dataonly = argv[i][0] == '-';
				switch (argv[i][2]) {
				case '\0':
				case 'a':		/* ascii */
					format = 'a';
					fmtid = "ascii";
					break;
				case 'i':		/* integer */
					format = 'i';
					fmtid = "ascii";
					break;
				case 'b':		/* byte */
					dataonly = 1;
					format = 'b';
					fmtid = "byte";
					break;
				case 'W':		/* 16-bit swapped */
					swapbytes = 1;
				case 'w':		/* 16-bit */
					dataonly = 1;
					format = 'w';
					fmtid = "16-bit";
					break;
				case 'f':		/* float */
					dataonly = 1;
					format = 'f';
					fmtid = "float";
					break;
				case 'd':		/* double */
					dataonly = 1;
					format = 'd';
					fmtid = "double";
					break;
				default:
					goto unkopt;
				}
				break;
			case 'x':		/* x resolution */
			case 'X':		/* x resolution */
				resolution = 0;
				if (argv[i][0] == '-')
					picres.rt |= XDECR;
				picres.xr = atoi(argv[++i]);
				break;
			case 'y':		/* y resolution */
			case 'Y':		/* y resolution */
				resolution = 0;
				if (argv[i][0] == '-')
					picres.rt |= YDECR;
				if (picres.xr == 0)
					picres.rt |= YMAJOR;
				picres.yr = atoi(argv[++i]);
				break;
			default:
unkopt:
				fprintf(stderr, "%s: unknown option: %s\n",
						progname, argv[i]);
				quit(1);
				break;
			}
		else
			break;
					/* recognize special formats */
	if (dataonly && format == 'b')
		if (putprim == ALL)
			fmtid = "24-bit_rgb";
		else
			fmtid = "8-bit_grey";
	if (dataonly && format == 'w')
		if (putprim == ALL)
			fmtid = "48-bit_rgb";
		else
			fmtid = "16-bit_grey";
					/* assign reverse ordering */
	rord[ord[0]] = 0;
	rord[ord[1]] = 1;
	rord[ord[2]] = 2;
					/* get input */
	if (i == argc) {
		fin = stdin;
	} else if (i < argc) {
		if ((fin = fopen(argv[i], "r")) == NULL) {
			fprintf(stderr, "%s: can't open file \"%s\"\n",
						progname, argv[i]);
			quit(1);
		}
		if (reverse && putprim != BRIGHT && i == argc-3) {
			if ((fin2 = fopen(argv[i+1], "r")) == NULL) {
				fprintf(stderr, "%s: can't open file \"%s\"\n",
						progname, argv[i+1]);
				quit(1);
			}
			if ((fin3 = fopen(argv[i+2], "r")) == NULL) {
				fprintf(stderr, "%s: can't open file \"%s\"\n",
						progname, argv[i+2]);
				quit(1);
			}
			interleave = -1;
		} else if (i != argc-1)
			fin = NULL;
		if (reverse && putprim != BRIGHT && !interleave) {
			fin2 = fopen(argv[i], "r");
			fin3 = fopen(argv[i], "r");
		}
		if (skipbytes && (fseek(fin, skipbytes, 0) || (fin2 != NULL &&
				(fseek(fin2, skipbytes, 0) ||
				fseek(fin3, skipbytes, 0))))) {
			fprintf(stderr, "%s: cannot skip %ld bytes on input\n",
					progname, skipbytes);
			quit(1);
		}
	}
	if (fin == NULL) {
		fprintf(stderr, "%s: bad # file arguments\n", progname);
		quit(1);
	}

	if (reverse) {
#ifdef _WIN32
		SET_FILE_BINARY(stdout);
		if (format != 'a' && format != 'i')
			SET_FILE_BINARY(fin);
#endif
					/* get header */
		if (header) {
			if (checkheader(fin, fmtid, stdout) < 0) {
				fprintf(stderr, "%s: wrong input format\n",
						progname);
				quit(1);
			}
			if (fin2 != NULL) {
				getheader(fin2, NULL, NULL);
				getheader(fin3, NULL, NULL);
			}
		} else
			newheader("RADIANCE", stdout);
					/* get resolution */
		if ((resolution && !fgetsresolu(&picres, fin)) ||
				picres.xr <= 0 || picres.yr <= 0) {
			fprintf(stderr, "%s: missing resolution\n", progname);
			quit(1);
		}
		if (resolution && fin2 != NULL) {
			RESOLU  pres2;
			if (!fgetsresolu(&pres2, fin2) ||
					pres2.rt != picres.rt ||
					pres2.xr != picres.xr ||
					pres2.yr != picres.yr ||
					!fgetsresolu(&pres2, fin3) ||
					pres2.rt != picres.rt ||
					pres2.xr != picres.xr ||
					pres2.yr != picres.yr) {
				fprintf(stderr, "%s: resolution mismatch\n",
						progname);
				quit(1);
			}
		}
						/* add to header */
		printargs(i, argv, stdout);
		if (expval < .99 || expval > 1.01)
			fputexpos(expval, stdout);
		fputformat(COLRFMT, stdout);
		putchar('\n');
		fputsresolu(&picres, stdout);	/* always put resolution */
		valtopix();
	} else {
#ifdef _WIN32
		SET_FILE_BINARY(fin);
		if (format != 'a' && format != 'i')
			SET_FILE_BINARY(stdout);
#endif
						/* get header */
		getheader(fin, checkhead, NULL);
		if (wrongformat) {
			fprintf(stderr,
				"%s: input not a Radiance RGBE picture\n",
					progname);
			quit(1);
		}
		if (!fgetsresolu(&picres, fin)) {
			fprintf(stderr, "%s: missing resolution\n", progname);
			quit(1);
		}
		if (header) {
			printargs(i, argv, stdout);
			if (expval < .99 || expval > 1.01)
				fputexpos(expval, stdout);
			fputformat(fmtid, stdout);
			putchar('\n');
		}
		if (resolution)			/* put resolution */
			fputsresolu(&picres, stdout);
		pixtoval();
	}

	quit(0);
}


int
checkhead(line)				/* deal with line from header */
char  *line;
{
	char	fmt[32];
	double	d;
	COLOR	ctmp;

	if (formatval(fmt, line)) {
		if (!strcmp(fmt, CIEFMT)) {
			mybright = &xyz_bright;
			if (original) {
				scalecolor(exposure, 1./WHTEFFICACY);
				doexposure++;
			}
		} else if (!strcmp(fmt, COLRFMT))
			mybright = &rgb_bright;
		else
			wrongformat++;
	} else if (original && isexpos(line)) {
		d = 1.0/exposval(line);
		scalecolor(exposure, d);
		doexposure++;
	} else if (original && iscolcor(line)) {
		colcorval(ctmp, line);
		setcolor(exposure, colval(exposure,RED)/colval(ctmp,RED),
				colval(exposure,GRN)/colval(ctmp,GRN),
				colval(exposure,BLU)/colval(ctmp,BLU));
		doexposure++;
	} else if (header)
		fputs(line, stdout);
	return(0);
}


pixtoval()				/* convert picture to values */
{
	register COLOR	*scanln;
	int  dogamma;
	COLOR  lastc;
	FLOAT  hv[2];
	int  startprim, endprim;
	long  startpos;
	int  y;
	register int  x;

	scanln = (COLOR *)malloc(scanlen(&picres)*sizeof(COLOR));
	if (scanln == NULL) {
		fprintf(stderr, "%s: out of memory\n", progname);
		quit(1);
	}
	dogamma = gamcor < .95 || gamcor > 1.05;
	if (putprim == ALL && !interleave) {
		startprim = RED; endprim = BLU;
		startpos = ftell(fin);
	} else {
		startprim = putprim; endprim = putprim;
	}
	for (putprim = startprim; putprim <= endprim; putprim++) {
		if (putprim != startprim && fseek(fin, startpos, 0)) {
			fprintf(stderr, "%s: seek error on input file\n",
					progname);
			quit(1);
		}
		set_io();
		setcolor(lastc, 0.0, 0.0, 0.0);
		for (y = 0; y < numscans(&picres); y++) {
			if (freadscan(scanln, scanlen(&picres), fin) < 0) {
				fprintf(stderr, "%s: read error\n", progname);
				quit(1);
			}
			for (x = 0; x < scanlen(&picres); x++) {
				if (uniq)
					if (	colval(scanln[x],RED) ==
							colval(lastc,RED) &&
						colval(scanln[x],GRN) ==
							colval(lastc,GRN) &&
						colval(scanln[x],BLU) ==
							colval(lastc,BLU)	)
						continue;
					else
						copycolor(lastc, scanln[x]);
				if (doexposure)
					multcolor(scanln[x], exposure);
				if (dogamma)
					setcolor(scanln[x],
					pow(colval(scanln[x],RED), 1.0/gamcor),
					pow(colval(scanln[x],GRN), 1.0/gamcor),
					pow(colval(scanln[x],BLU), 1.0/gamcor));
				if (!dataonly) {
					pix2loc(hv, &picres, x, y);
					printf("%7d %7d ",
							(int)(hv[0]*picres.xr),
							(int)(hv[1]*picres.yr));
				}
				if ((*putval)(scanln[x]) < 0) {
					fprintf(stderr, "%s: write error\n",
							progname);
					quit(1);
				}
			}
		}
	}
	free((void *)scanln);
}


valtopix()			/* convert values to a pixel file */
{
	int  dogamma;
	register COLOR	*scanln;
	int  y;
	register int  x;

	scanln = (COLOR *)malloc(scanlen(&picres)*sizeof(COLOR));
	if (scanln == NULL) {
		fprintf(stderr, "%s: out of memory\n", progname);
		quit(1);
	}
	dogamma = gamcor < .95 || gamcor > 1.05;
	set_io();
	for (y = 0; y < numscans(&picres); y++) {
		for (x = 0; x < scanlen(&picres); x++) {
			if (!dataonly) {
				fscanf(fin, "%*d %*d");
				if (fin2 != NULL) {
					fscanf(fin2, "%*d %*d");
					fscanf(fin3, "%*d %*d");
				}
			}
			if ((*getval)(scanln[x]) < 0) {
				fprintf(stderr, "%s: read error\n", progname);
				quit(1);
			}
			if (dogamma)
				setcolor(scanln[x],
					pow(colval(scanln[x],RED), gamcor),
					pow(colval(scanln[x],GRN), gamcor),
					pow(colval(scanln[x],BLU), gamcor));
			if (doexposure)
				multcolor(scanln[x], exposure);
		}
		if (fwritescan(scanln, scanlen(&picres), stdout) < 0) {
			fprintf(stderr, "%s: write error\n", progname);
			quit(1);
		}
	}
	free((void *)scanln);
}


void
quit(code)
int  code;
{
	exit(code);
}


swap16(wp, n)		/* swap n 16-bit words */
register uint16  *wp;
int  n;
{
	while (n-- > 0) {
		*wp = *wp << 8 | ((*wp >> 8) & 0xff);
		wp++;
	}
}

getcascii(col)		/* get an ascii color value from stream(s) */
COLOR  col;
{
	double	vd[3];

	if (fin2 == NULL) {
		if (fscanf(fin, "%lf %lf %lf", &vd[0], &vd[1], &vd[2]) != 3)
			return(-1);
	} else {
		if (fscanf(fin, "%lf", &vd[0]) != 1 ||
				fscanf(fin2, "%lf", &vd[1]) != 1 ||
				fscanf(fin3, "%lf", &vd[2]) != 1)
			return(-1);
	}
	setcolor(col, vd[rord[RED]], vd[rord[GRN]], vd[rord[BLU]]);
	return(0);
}


getcdouble(col)		/* get a double color value from stream(s) */
COLOR  col;
{
	double	vd[3];

	if (fin2 == NULL) {
		if (fread((char *)vd, sizeof(double), 3, fin) != 3)
			return(-1);
	} else {
		if (fread((char *)vd, sizeof(double), 1, fin) != 1 ||
			fread((char *)(vd+1), sizeof(double), 1, fin2) != 1 ||
			fread((char *)(vd+2), sizeof(double), 1, fin3) != 1)
			return(-1);
	}
	setcolor(col, vd[rord[RED]], vd[rord[GRN]], vd[rord[BLU]]);
	return(0);
}


getcfloat(col)		/* get a float color value from stream(s) */
COLOR  col;
{
	float  vf[3];

	if (fin2 == NULL) {
		if (fread((char *)vf, sizeof(float), 3, fin) != 3)
			return(-1);
	} else {
		if (fread((char *)vf, sizeof(float), 1, fin) != 1 ||
			fread((char *)(vf+1), sizeof(float), 1, fin2) != 1 ||
			fread((char *)(vf+2), sizeof(float), 1, fin3) != 1)
			return(-1);
	}
	setcolor(col, vf[rord[RED]], vf[rord[GRN]], vf[rord[BLU]]);
	return(0);
}


getcint(col)		/* get an int color value from stream(s) */
COLOR  col;
{
	int  vi[3];

	if (fin2 == NULL) {
		if (fscanf(fin, "%d %d %d", &vi[0], &vi[1], &vi[2]) != 3)
			return(-1);
	} else {
		if (fscanf(fin, "%d", &vi[0]) != 1 ||
				fscanf(fin2, "%d", &vi[1]) != 1 ||
				fscanf(fin3, "%d", &vi[2]) != 1)
			return(-1);
	}
	setcolor(col, (vi[rord[RED]]+.5)/256.,
			(vi[rord[GRN]]+.5)/256., (vi[rord[BLU]]+.5)/256.);
	return(0);
}


getcbyte(col)		/* get a byte color value from stream(s) */
COLOR  col;
{
	BYTE  vb[3];

	if (fin2 == NULL) {
		if (fread((char *)vb, sizeof(BYTE), 3, fin) != 3)
			return(-1);
	} else {
		if (fread((char *)vb, sizeof(BYTE), 1, fin) != 1 ||
			fread((char *)(vb+1), sizeof(BYTE), 1, fin2) != 1 ||
			fread((char *)(vb+2), sizeof(BYTE), 1, fin3) != 1)
			return(-1);
	}
	setcolor(col, (vb[rord[RED]]+.5)/256.,
			(vb[rord[GRN]]+.5)/256., (vb[rord[BLU]]+.5)/256.);
	return(0);
}


getcword(col)		/* get a 16-bit color value from stream(s) */
COLOR  col;
{
	uint16  vw[3];

	if (fin2 == NULL) {
		if (fread((char *)vw, sizeof(uint16), 3, fin) != 3)
			return(-1);
	} else {
		if (fread((char *)vw, sizeof(uint16), 1, fin) != 1 ||
			fread((char *)(vw+1), sizeof(uint16), 1, fin2) != 1 ||
			fread((char *)(vw+2), sizeof(uint16), 1, fin3) != 1)
			return(-1);
	}
	if (swapbytes)
		swap16(vw, 3);
	setcolor(col, (vw[rord[RED]]+.5)/65536.,
			(vw[rord[GRN]]+.5)/65536., (vw[rord[BLU]]+.5)/65536.);
	return(0);
}


getbascii(col)		/* get an ascii brightness value from fin */
COLOR  col;
{
	double	vd;

	if (fscanf(fin, "%lf", &vd) != 1)
		return(-1);
	setcolor(col, vd, vd, vd);
	return(0);
}


getbdouble(col)		/* get a double brightness value from fin */
COLOR  col;
{
	double	vd;

	if (fread((char *)&vd, sizeof(double), 1, fin) != 1)
		return(-1);
	setcolor(col, vd, vd, vd);
	return(0);
}


getbfloat(col)		/* get a float brightness value from fin */
COLOR  col;
{
	float  vf;

	if (fread((char *)&vf, sizeof(float), 1, fin) != 1)
		return(-1);
	setcolor(col, vf, vf, vf);
	return(0);
}


getbint(col)		/* get an int brightness value from fin */
COLOR  col;
{
	int  vi;
	double	d;

	if (fscanf(fin, "%d", &vi) != 1)
		return(-1);
	d = (vi+.5)/256.;
	setcolor(col, d, d, d);
	return(0);
}


getbbyte(col)		/* get a byte brightness value from fin */
COLOR  col;
{
	BYTE  vb;
	double	d;

	if (fread((char *)&vb, sizeof(BYTE), 1, fin) != 1)
		return(-1);
	d = (vb+.5)/256.;
	setcolor(col, d, d, d);
	return(0);
}


getbword(col)		/* get a 16-bit brightness value from fin */
COLOR  col;
{
	uint16  vw;
	double	d;

	if (fread((char *)&vw, sizeof(uint16), 1, fin) != 1)
		return(-1);
	if (swapbytes)
		swap16(&vw, 1);
	d = (vw+.5)/65536.;
	setcolor(col, d, d, d);
	return(0);
}


putcascii(col)			/* put an ascii color to stdout */
COLOR  col;
{
	fprintf(stdout, "%15.3e %15.3e %15.3e\n",
			colval(col,ord[0]),
			colval(col,ord[1]),
			colval(col,ord[2]));

	return(ferror(stdout) ? -1 : 0);
}


putcfloat(col)			/* put a float color to stdout */
COLOR  col;
{
	float  vf[3];

	vf[0] = colval(col,ord[0]);
	vf[1] = colval(col,ord[1]);
	vf[2] = colval(col,ord[2]);
	fwrite((char *)vf, sizeof(float), 3, stdout);

	return(ferror(stdout) ? -1 : 0);
}


putcdouble(col)			/* put a double color to stdout */
COLOR  col;
{
	double	vd[3];

	vd[0] = colval(col,ord[0]);
	vd[1] = colval(col,ord[1]);
	vd[2] = colval(col,ord[2]);
	fwrite((char *)vd, sizeof(double), 3, stdout);

	return(ferror(stdout) ? -1 : 0);
}


putcint(col)			/* put an int color to stdout */
COLOR  col;
{
	fprintf(stdout, "%d %d %d\n",
			(int)(colval(col,ord[0])*256.),
			(int)(colval(col,ord[1])*256.),
			(int)(colval(col,ord[2])*256.));

	return(ferror(stdout) ? -1 : 0);
}


putcbyte(col)			/* put a byte color to stdout */
COLOR  col;
{
	long  i;
	BYTE  vb[3];

	i = colval(col,ord[0])*256.;
	vb[0] = min(i,255);
	i = colval(col,ord[1])*256.;
	vb[1] = min(i,255);
	i = colval(col,ord[2])*256.;
	vb[2] = min(i,255);
	fwrite((char *)vb, sizeof(BYTE), 3, stdout);

	return(ferror(stdout) ? -1 : 0);
}


putcword(col)			/* put a 16-bit color to stdout */
COLOR  col;
{
	long  i;
	uint16  vw[3];

	i = colval(col,ord[0])*65536.;
	vw[0] = min(i,65535);
	i = colval(col,ord[1])*65536.;
	vw[1] = min(i,65535);
	i = colval(col,ord[2])*65536.;
	vw[2] = min(i,65535);
	if (swapbytes)
		swap16(vw, 3);
	fwrite((char *)vw, sizeof(uint16), 3, stdout);

	return(ferror(stdout) ? -1 : 0);
}


putbascii(col)			/* put an ascii brightness to stdout */
COLOR  col;
{
	fprintf(stdout, "%15.3e\n", (*mybright)(col));

	return(ferror(stdout) ? -1 : 0);
}


putbfloat(col)			/* put a float brightness to stdout */
COLOR  col;
{
	float  vf;

	vf = (*mybright)(col);
	fwrite((char *)&vf, sizeof(float), 1, stdout);

	return(ferror(stdout) ? -1 : 0);
}


putbdouble(col)			/* put a double brightness to stdout */
COLOR  col;
{
	double	vd;

	vd = (*mybright)(col);
	fwrite((char *)&vd, sizeof(double), 1, stdout);

	return(ferror(stdout) ? -1 : 0);
}


putbint(col)			/* put an int brightness to stdout */
COLOR  col;
{
	fprintf(stdout, "%d\n", (int)((*mybright)(col)*256.));

	return(ferror(stdout) ? -1 : 0);
}


putbbyte(col)			/* put a byte brightness to stdout */
COLOR  col;
{
	register int  i;
	BYTE  vb;

	i = (*mybright)(col)*256.;
	vb = min(i,255);
	fwrite((char *)&vb, sizeof(BYTE), 1, stdout);

	return(ferror(stdout) ? -1 : 0);
}


putbword(col)			/* put a 16-bit brightness to stdout */
COLOR  col;
{
	long  i;
	uint16  vw;

	i = (*mybright)(col)*65536.;
	vw = min(i,65535);
	if (swapbytes)
		swap16(&vw, 1);
	fwrite((char *)&vw, sizeof(uint16), 1, stdout);

	return(ferror(stdout) ? -1 : 0);
}


putpascii(col)			/* put an ascii primary to stdout */
COLOR  col;
{
	fprintf(stdout, "%15.3e\n", colval(col,putprim));

	return(ferror(stdout) ? -1 : 0);
}


putpfloat(col)			/* put a float primary to stdout */
COLOR  col;
{
	float  vf;

	vf = colval(col,putprim);
	fwrite((char *)&vf, sizeof(float), 1, stdout);

	return(ferror(stdout) ? -1 : 0);
}


putpdouble(col)			/* put a double primary to stdout */
COLOR  col;
{
	double	vd;

	vd = colval(col,putprim);
	fwrite((char *)&vd, sizeof(double), 1, stdout);

	return(ferror(stdout) ? -1 : 0);
}


putpint(col)			/* put an int primary to stdout */
COLOR  col;
{
	fprintf(stdout, "%d\n", (int)(colval(col,putprim)*256.));

	return(ferror(stdout) ? -1 : 0);
}


putpbyte(col)			/* put a byte primary to stdout */
COLOR  col;
{
	long  i;
	BYTE  vb;

	i = colval(col,putprim)*256.;
	vb = min(i,255);
	fwrite((char *)&vb, sizeof(BYTE), 1, stdout);

	return(ferror(stdout) ? -1 : 0);
}


putpword(col)			/* put a 16-bit primary to stdout */
COLOR  col;
{
	long  i;
	uint16  vw;

	i = colval(col,putprim)*65536.;
	vw = min(i,65535);
	if (swapbytes)
		swap16(&vw, 1);
	fwrite((char *)&vw, sizeof(uint16), 1, stdout);

	return(ferror(stdout) ? -1 : 0);
}


set_io()			/* set put and get functions */
{
	switch (format) {
	case 'a':					/* ascii */
		if (putprim == BRIGHT) {
			getval = getbascii;
			putval = putbascii;
		} else if (putprim != ALL) {
			getval = getbascii;
			putval = putpascii;
		} else {
			getval = getcascii;
			putval = putcascii;
			if (reverse && !interleave) {
				fprintf(stderr,
				"%s: ASCII input files must be interleaved\n",
						progname);
				quit(1);
			}
		}
		return;
	case 'f':					/* binary float */
		if (putprim == BRIGHT) {
			getval = getbfloat;
			putval = putbfloat;
		} else if (putprim != ALL) {
			getval = getbfloat;
			putval = putpfloat;
		} else {
			getval = getcfloat;
			putval = putcfloat;
			if (reverse && !interleave) {
				if (fin2 == NULL)
					goto namerr;
				if (fseek(fin2,
				(long)sizeof(float)*picres.xr*picres.yr, 1))
					goto seekerr;
				if (fseek(fin3,
				(long)sizeof(float)*2*picres.xr*picres.yr, 1))
					goto seekerr;
			}
		}
		return;
	case 'd':					/* binary double */
		if (putprim == BRIGHT) {
			getval = getbdouble;
			putval = putbdouble;
		} else if (putprim != ALL) {
			getval = getbdouble;
			putval = putpdouble;
		} else {
			getval = getcdouble;
			putval = putcdouble;
			if (reverse && !interleave) {
				if (fin2 == NULL)
					goto namerr;
				if (fseek(fin2,
				(long)sizeof(double)*picres.xr*picres.yr, 1))
					goto seekerr;
				if (fseek(fin3,
				(long)sizeof(double)*2*picres.xr*picres.yr, 1))
					goto seekerr;
			}
		}
		return;
	case 'i':					/* integer */
		if (putprim == BRIGHT) {
			getval = getbint;
			putval = putbint;
		} else if (putprim != ALL) {
			getval = getbint;
			putval = putpint;
		} else {
			getval = getcint;
			putval = putcint;
			if (reverse && !interleave) {
				fprintf(stderr,
				"%s: integer input files must be interleaved\n",
						progname);
				quit(1);
			}
		}
		return;
	case 'b':					/* byte */
		if (putprim == BRIGHT) {
			getval = getbbyte;
			putval = putbbyte;
		} else if (putprim != ALL) {
			getval = getbbyte;
			putval = putpbyte;
		} else {
			getval = getcbyte;
			putval = putcbyte;
			if (reverse && !interleave) {
				if (fin2 == NULL)
					goto namerr;
				if (fseek(fin2,
				(long)sizeof(BYTE)*picres.xr*picres.yr, 1))
					goto seekerr;
				if (fseek(fin3,
				(long)sizeof(BYTE)*2*picres.xr*picres.yr, 1))
					goto seekerr;
			}
		}
		return;
	case 'w':					/* 16-bit */
		if (putprim == BRIGHT) {
			getval = getbword;
			putval = putbword;
		} else if (putprim != ALL) {
			getval = getbword;
			putval = putpword;
		} else {
			getval = getcword;
			putval = putcword;
			if (reverse && !interleave) {
				if (fin2 == NULL)
					goto namerr;
				if (fseek(fin2,
				(long)sizeof(uint16)*picres.xr*picres.yr, 1))
					goto seekerr;
				if (fseek(fin3,
				(long)sizeof(uint16)*2*picres.xr*picres.yr, 1))
					goto seekerr;
			}
		}
		return;
	}
badopt:
	fprintf(stderr, "%s: botched file type\n", progname);
	quit(1);
namerr:
	fprintf(stderr, "%s: non-interleaved file(s) must be named\n",
			progname);
	quit(1);
seekerr:
	fprintf(stderr, "%s: cannot seek on interleaved input file\n",
			progname);
	quit(1);
}
