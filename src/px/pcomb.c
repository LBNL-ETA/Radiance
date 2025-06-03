#ifndef lint
static const char	RCSid[] = "$Id: pcomb.c,v 2.61 2025/06/03 21:31:51 greg Exp $";
#endif
/*
 *  Combine picture files according to calcomp functions.
 *
 *	1/4/89
 */

#include "platform.h"
#include "standard.h"
#include "paths.h"
#include "color.h"
#include "calcomp.h"
#include "view.h"

#define MAXINP		1024		/* maximum number of input files */
#define WINSIZ		127		/* scanline window size */
#define MIDSCN		((WINSIZ-1)/2+1)

struct {
	char	*name;		/* file or command name */
	FILE	*fp;		/* stream pointer */
	VIEW	vw;		/* view for picture */
	RESOLU	rs;		/* image resolution and orientation */
	int	infloat;	/* input is floating point (#comp)? */
	float	pa;		/* pixel aspect ratio */
	COLOR	*scan[WINSIZ];	/* input scanline window */
	COLOR	coef;		/* coefficient */
	COLOR	expos;		/* recorded exposure */
}	input[MAXINP];			/* input pictures */

int	nfiles;				/* number of input files */

VIEW	*commvp = NULL;			/* common view parameters */

char	ourfmt[MAXFMTLEN] = PICFMT;	/* input picture format */
int	outfloat = 0;			/* #component float output? */

char	StandardInput[] = "<stdin>";
char	Command[] = "<Command>";
char	vcolin[3][4] = {"ri", "gi", "bi"};
char	vcolout[3][4] = {"ro", "go", "bo"};
char	vbrtin[] = "li";
char	vbrtout[] = "lo";
char	vcolexp[3][4] = {"re", "ge", "be"};
char	vbrtexp[] = "le";
char	vpixaspect[] = "pa";

char	vray[7][4] = {"Ox", "Oy", "Oz", "Dx", "Dy", "Dz", "T"};

char	vpsize[] = "S";

char	vnfiles[] = "nfiles";
char	vwhteff[] = "WE";
char	vxmax[] = "xmax";
char	vymax[] = "ymax";
char	vxres[] = "xres";
char	vyres[] = "yres";
char	vxpos[] = "x";
char	vypos[] = "y";

int	nowarn = 0;			/* no warning messages? */

int	xmax = 0, ymax = 0;		/* input resolution */

int	xscan, yscan;			/* input position */

int	xres, yres;			/* output resolution */

int	xpos, ypos;			/* output position */

int	echoheader = 1;
int	wrongformat = 0;
int	gotview;


static gethfunc headline;
static void checkfile(int orig);
static double rgb_bright(COLOR  clr);
static double xyz_bright(COLOR  clr);
static void init(void);
static void combine(void);
static void advance(void);
static double l_expos(char	*nam);
static double l_pixaspect(char *nm);
static double l_colin(char	*nam);
static double l_ray(char	*nam);
static double l_psize(char *nm);


int
main(
	int	argc,
	char	*argv[]
)
{
	int	original;
	double	f;
	int	a;

	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	fixargv0(argv[0]);		/* sets global progname */
	esupport |= E_VARIABLE|E_FUNCTION|E_RCONST;
	esupport &= ~(E_OUTCHAN|E_INCHAN);
						/* scan options */
	for (a = 1; a < argc; a++) {
		if (argv[a][0] == '-')
			switch (argv[a][1]) {
			case 'x':
			case 'y':
				a++;
				continue;
			case 'w':
				nowarn = !nowarn;
				continue;
			case 'h':
				echoheader = !echoheader;
				continue;
			case 'f':
				if (argv[a][2]) {
					outfloat = (argv[a][2] == 'f');
					continue;
				}
				/* fall through */
			case 'e':
				a++;
				continue;
			}
		break;
	}
	newheader("RADIANCE", stdout);	/* start header */
	fputnow(stdout);
					/* process files */
	for (nfiles = 0; nfiles < MAXINP; nfiles++) {
		setcolor(input[nfiles].coef, 1.0, 1.0, 1.0);
		setcolor(input[nfiles].expos, 1.0, 1.0, 1.0);
		input[nfiles].vw = stdview;
		input[nfiles].pa = 1.0;
	}
	nfiles = 0;
	original = 0;
	for ( ; a < argc; a++) {
		if (nfiles >= MAXINP) {
			eputs(argv[0]);
			eputs(": too many picture files\n");
			quit(1);
		}
		if (argv[a][0] == '-')
			switch (argv[a][1]) {
			case '\0':
				input[nfiles].name = StandardInput;
				input[nfiles].fp = stdin;
				break;
			case 'o':
				original++;
				continue;
			case 's':
				f = atof(argv[++a]);
				scalecolor(input[nfiles].coef, f);
				continue;
			case 'c':
				colval(input[nfiles].coef,RED)*=atof(argv[++a]);
				colval(input[nfiles].coef,GRN)*=atof(argv[++a]);
				colval(input[nfiles].coef,BLU)*=atof(argv[++a]);
				continue;
			default:
				goto usage;
			}
		else {
			if (argv[a][0] == '!') {
				input[nfiles].name = Command;
				input[nfiles].fp = popen(argv[a]+1, "r");
			} else {
				input[nfiles].name = argv[a];
				input[nfiles].fp = fopen(argv[a], "r");
			}
			if (input[nfiles].fp == NULL) {
				perror(argv[a]);
				quit(1);
			}
		}
		checkfile(original);
		nfiles++;
		original = 0;
	}
	init();				/* set constants */
					/* go back and get expressions */
	for (a = 1; a < argc; a++) {
		char	*fpath;
		if (argv[a][0] == '-')
			switch (argv[a][1]) {
			case 'x':
				varset(vxres, ':', eval(argv[++a]));
				continue;
			case 'y':
				varset(vyres, ':', eval(argv[++a]));
				continue;
			case 'w':
				continue;
			case 'h':
				continue;
			case 'f':
				if (argv[a][2])
					continue;
				fpath = getpath(argv[++a], getrlibpath(), 0);
				if (fpath == NULL) {
					eputs(argv[0]);
					eputs(": cannot find file '");
					eputs(argv[a]);
					eputs("'\n");
					quit(1);
				}
				fcompile(fpath);
				continue;
			case 'e':
				scompile(argv[++a], NULL, 0);
				continue;
			}
		break;
	}
						/* get output resolution */
	xres = varvalue(vxres) + .5;
	yres = varvalue(vyres) + .5;
	if (xres <= 0 || yres <= 0) {
		eputs(argv[0]);
		eputs(": illegal output resolution\n");
		quit(1);
	}
	if (!vardefined(vbrtout))		/* single or 3-channel? */
		outfloat *= 3;

	printargs(argc, argv, stdout);		/* complete header */
	if (commvp != NULL) {
		fputs(VIEWSTR, stdout);
		fprintview(commvp, stdout);
		fputc('\n', stdout);
	}
	if (outfloat) {				/* print output format */
		printf("NROWS=%d\nNCOLS=%d\n", yres, xres);
		fputncomp(outfloat, stdout);
		fputendian(stdout);
		fputformat("float", stdout);
	} else if (strcmp(ourfmt, PICFMT))
		fputformat(ourfmt, stdout);
	else
		fputformat(COLRFMT, stdout);
	fputc('\n', stdout);			/* end header */
	if (!outfloat)
		fprtresolu(xres, yres, stdout);

	doptimize(1);				/* optimize definitions */
	combine();				/* combine pictures */
	quit(0);
usage:
	eputs("Usage: ");
	eputs(argv[0]);
	eputs(
" [-w][-h][-ff][-x xr][-y yr][-e expr][-f file] [ [-o][-s f][-c r g b] hdr ..]\n");
	quit(1);
	return 1; /* pro forma return */
}


static int
headline(			/* check header line & echo if requested */
	char	*s,
	void	*p
)
{
	int	orig = *(int *)p;
	char	fmt[MAXFMTLEN];
	double	d;
	int	bigend;
	COLOR	ctmp;

	if (isheadid(s))			/* header id */
		return(0);	/* don't echo */
	if (formatval(fmt, s)) {		/* check format */
		if (globmatch(ourfmt, fmt)) {
			wrongformat = 0;
			strcpy(ourfmt, fmt);
		} else if (!strcmp("float", fmt))
			input[nfiles].infloat *= -1;
		else
			wrongformat = globmatch(PICFMT, fmt) ? 1 : -1;
		return(0);	/* don't echo */
	}
	if ((bigend = isbigendian(s)) >= 0) {
		if (bigend != nativebigendian()) {
			eputs(input[nfiles].name);
			eputs(": unsupported input byte ordering\n");
			quit(1);
		}
		return(0);	/* don't echo */
	}
	if (!strncmp("NROWS=", s, 6)) {		/* X-resolution */
		input[nfiles].rs.yr = atoi(s+6);
		return(0);	/* don't echo */
	}
	if (!strncmp("NCOLS=", s, 6)) {		/* Y-resolution */
		input[nfiles].rs.xr = atoi(s+6);
		return(0);	/* don't echo */
	}
	if (isncomp(s))				/* # components */
		input[nfiles].infloat *= ncompval(s);

	if (orig) {				/* undo exposure? */
		if (isexpos(s)) {
			d = 1./exposval(s);
			scalecolor(input[nfiles].coef, d);
			return(0);	/* don't echo */
		}
		if (iscolcor(s)) {
			int	i;
			colcorval(ctmp, s);
			for (i = 3; i--; )
				colval(input[nfiles].coef,i) /= colval(ctmp,i);
			return(0);	/* don't echo */
		}
	} else if (isexpos(s)) {		/* record exposure? */
		d = exposval(s);
		scalecolor(input[nfiles].expos, d);
	} else if (iscolcor(s)) {
		colcorval(ctmp, s);
		multcolor(input[nfiles].expos, ctmp);
	}
	if (isaspect(s))
		input[nfiles].pa *= aspectval(s);
	else if (isview(s) && sscanview(&input[nfiles].vw, s) > 0)
		gotview++;

	if (echoheader) {			/* echo line */
		putchar('\t');
		return(fputs(s, stdout));
	}
	return(0);
}


static void
checkfile(int orig)			/* ready a file */
{
	int	i;
					/* process header */
	gotview = 0;
	input[nfiles].infloat = -1;
	input[nfiles].rs.rt = PIXSTANDARD;
	input[nfiles].rs.xr = input[nfiles].rs.yr = 0;
	if (echoheader) {
		fputs(input[nfiles].name, stdout);
		fputs(":\n", stdout);
	}
	getheader(input[nfiles].fp, headline, &orig);

	if (input[nfiles].infloat <= 0)
		input[nfiles].infloat = 0;
	else if ((input[nfiles].infloat != 3) &
				(input[nfiles].infloat != 1)) {
		eputs(input[nfiles].name);
		eputs(": unsupported number of components\n");
		quit(1);
	}
	if (wrongformat < 0) {
		eputs(input[nfiles].name);
		eputs(": not a Radiance picture or float matrix\n");
		quit(1);
	}
	if (wrongformat > 0) {
		wputs(input[nfiles].name);
		wputs(": warning -- incompatible picture format\n");
	}
	if (!gotview || setview(&input[nfiles].vw) != NULL)
		input[nfiles].vw.type = 0;
	else if (commvp == NULL)
		commvp = &input[nfiles].vw;
	if ((input[nfiles].rs.xr <= 0) | (input[nfiles].rs.yr <= 0) &&
			!fgetsresolu(&input[nfiles].rs, input[nfiles].fp)) {
		eputs(input[nfiles].name);
		eputs(": bad picture size\n");
		quit(1);
	}
	if (xmax == 0 && ymax == 0) {
		xmax = scanlen(&input[nfiles].rs);
		ymax = numscans(&input[nfiles].rs);
	} else if (scanlen(&input[nfiles].rs) != xmax ||
			numscans(&input[nfiles].rs) != ymax) {
		eputs(input[nfiles].name);
		eputs(": resolution mismatch\n");
		quit(1);
	}
					/* allocate scanlines */
	for (i = 0; i < WINSIZ; i++)
		input[nfiles].scan[i] = (COLOR *)emalloc(xmax*sizeof(COLOR));
}


static double
rgb_bright(
	COLOR  clr
)
{
	return(bright(clr));
}


static double
xyz_bright(
	COLOR  clr
)
{
	return(clr[CIEY]);
}


double	(*ourbright)() = rgb_bright;


static void
init(void)					/* perform final setup */
{
	int	i;
						/* define constants */
	varset("PI", ':', PI);
	varset(vnfiles, ':', (double)nfiles);
	varset(vxmax, ':', (double)xmax);
	varset(vymax, ':', (double)ymax);
						/* set functions */
	for (i = 0; i < 3; i++) {
		funset(vcolexp[i], 1, ':', l_expos);
		funset(vcolin[i], 1, '=', l_colin);
	}
	funset(vbrtexp, 1, ':', l_expos);
	funset(vbrtin, 1, '=', l_colin);
	funset(vpixaspect, 1, ':', l_pixaspect);
	for (i = 0; i < 7; i++)
		funset(vray[i], 1, '=', l_ray);
	funset(vpsize, 1, '=', l_psize);
						/* set brightness function */
	if (!strcmp(ourfmt, CIEFMT)) {
		varset(vwhteff, ':', 1.0);
		ourbright = xyz_bright;
	} else
		varset(vwhteff, ':', WHTEFFICACY);
						/* these may be overridden */
	varset(vxres, ':', (double)xmax);
	varset(vyres, ':', (double)ymax);
}


static void
combine(void)			/* combine pictures */
{
	EPNODE	*coldef[3], *brtdef;
	int	set_x, set_y;
	COLOR	*scanout;
	double	d;
	int	i, j;
						/* check defined variables */
	for (j = 0; j < 3; j++) {
		if (vardefined(vcolout[j]))
			coldef[j] = eparse(vcolout[j]);
		else
			coldef[j] = NULL;
	}
	if (vardefined(vbrtout))
		brtdef = eparse(vbrtout);
	else
		brtdef = NULL;
						/* what to set */
	set_x = varlookup(vxpos) != NULL && !vardefined(vxpos);
	set_y = varlookup(vypos) != NULL && !vardefined(vypos);
						/* allocate scanline */
	scanout = (COLOR *)emalloc(xres*sizeof(COLOR));
						/* set input position */
	yscan = ymax+MIDSCN;
						/* combine files */
	for (ypos = yres-1; ypos >= 0; ypos--) {
	    advance();
	    if (set_y) varset(vypos, '=', (double)ypos);
	    for (xpos = 0; xpos < xres; xpos++) {
		xscan = (xpos+.5)*xmax/xres;
		if (set_x) varset(vxpos, '=', (double)xpos);
		eclock++;
		if (brtdef != NULL) {
		    d = evalue(brtdef);
		    d *= (outfloat > 0) | (d >= 0);
		    setcolor(scanout[xpos], d, d, d);
		} else {
		    for (j = 0; j < 3; j++) {
			if (coldef[j] != NULL) {
				d = evalue(coldef[j]);
			} else {
			    d = 0.0;
			    for (i = 0; i < nfiles; i++)
				d += colval(input[i].scan[MIDSCN][xscan],j);
			}
			d *= (outfloat > 0) | (d >= 0);
			colval(scanout[xpos],j) = d;
		    }
		}
	    }
	    switch (outfloat) {
	    case 3:			/* writing out float triplets */
	    	if (fwrite(scanout, sizeof(float)*3, xres, stdout) != xres)
		    goto writerr;
	    	break;
	    case 1:			/* writing out gray float values */
	    	for (xpos = 0; xpos < xres; xpos++)
		    if (putbinary(&scanout[xpos][CIEY],
		    		sizeof(float), 1, stdout) != 1)
		    	goto writerr;
	    	break;
	    default:			/* writing out Radiance picture */
	   	if (fwritescan(scanout, xres, stdout) < 0)
	    	    goto writerr;
		break;
	    }
	}
	efree((char *)scanout);
	return;
writerr:
    perror("write error");
    quit(1);
}


static void
advance(void)			/* read in data for next scanline */
{
	int	ytarget;
	COLOR	*st;
	int	i, j;

	for (ytarget = (ypos+.5)*ymax/yres; yscan > ytarget; yscan--)
		for (i = 0; i < nfiles; i++) {
			st = input[i].scan[WINSIZ-1];
			for (j = WINSIZ-1; j > 0; j--)	/* rotate window */
				input[i].scan[j] = input[i].scan[j-1];
			input[i].scan[0] = st;
			if (yscan <= MIDSCN)		/* hit bottom? */
				continue;
				 			 /* read scanline */
			if (input[i].infloat) {
				if (fread(st, sizeof(float)*input[i].infloat,
						xmax, input[i].fp) != xmax)
					goto readerr;
				if (input[i].infloat == 1) {
					const COLORV	*f = st[0] + xmax;
					int		x = xmax;
					while (x-- > 0)
						st[x][BLU] = st[x][GRN] =
							st[x][RED] = *--f;
				}
			} else if (freadscan(st, xmax, input[i].fp) < 0)
				goto readerr;
			for (j = 0; j < xmax; j++)	/* adjust color */
				multcolor(st[j], input[i].coef);
		}
	return;
readerr:
	eputs(input[i].name);
	eputs(": read error\n");
	quit(1);
}


static double
l_expos(			/* return picture exposure */
	char	*nam
)
{
	int	fn, n;
	double	d;

	d = argument(1);
	if (d <= -0.5 || d >= nfiles+0.5) {
		errno = EDOM;
		return(0.0);
	}
	if (d < 0.5)
		return((double)nfiles);
	fn = d - 0.5;
	if (nam == vbrtexp)
		return((*ourbright)(input[fn].expos));
	n = 3;
	while (n--)
		if (nam == vcolexp[n])
			return(colval(input[fn].expos,n));
	eputs("Bad call to l_expos()!\n");
	quit(1);
	return 1; /* pro forma return */
}


static double
l_pixaspect(char *nm)		/* return pixel aspect ratio */
{
	int	fn;
	double	d;

	d = argument(1);
	if (d <= -0.5 || d >= nfiles+0.5) {
		errno = EDOM;
		return(0.0);
	}
	if (d < 0.5)
		return((double)nfiles);
	fn = d - 0.5;
	return(input[fn].pa);
}


static double
l_colin(			/* return color value for picture */
	char	*nam
)
{
	int	fn;
	int	n, xoff, yoff;
	double	d;

	d = argument(1);
	if (d <= -0.5 || d >= nfiles+0.5) {
		errno = EDOM;
		return(0.0);
	}
	if (d < 0.5)
		return((double)nfiles);
	fn = d - 0.5;
	xoff = yoff = 0;
	n = nargum();
	if (n >= 2) {
		d = argument(2);
		if (d < 0.0) {
			xoff = d-.5;
			if (xscan+xoff < 0)
				xoff = -xscan;
		} else {
			xoff = d+.5;
			if (xscan+xoff >= xmax)
				xoff = xmax-1-xscan;
		}
	}
	if (n >= 3) {
		d = argument(3);
		if (d < 0.0) {
			yoff = d-.5;
			if (yoff+MIDSCN < 0)
				yoff = -MIDSCN;
			if (yscan+yoff < 0)
				yoff = -yscan;
		} else {
			yoff = d+.5;
			if (yoff+MIDSCN >= WINSIZ)
				yoff = WINSIZ-1-MIDSCN;
			if (yscan+yoff >= ymax)
				yoff = ymax-1-yscan;
		}
	}
	if (nam == vbrtin)
		return((*ourbright)(input[fn].scan[MIDSCN+yoff][xscan+xoff]));
	n = 3;
	while (n--)
	    if (nam == vcolin[n])
		return(colval(input[fn].scan[MIDSCN+yoff][xscan+xoff],n));
	eputs("Bad call to l_colin()!\n");
	quit(1);
	return 1; /* pro forma return */
}


static double
l_ray(		/* return ray origin or direction */
	char	*nam
)
{
	static unsigned long	ltick[MAXINP];
	static FVECT	lorg[MAXINP], ldir[MAXINP];
	static double	ldist[MAXINP];
	RREAL	loc[2];
	double	d;
	int	fn;
	int	i;

	d = argument(1);
	if (d <= -0.5 || d >= nfiles+0.5) {
		errno = EDOM;
		return(0.0);
	}
	if (d < 0.5)
		return((double)nfiles);
	fn = d - 0.5;
	if (ltick[fn] != eclock) {		/* need to compute? */
		if (input[fn].vw.type == 0) {
			lorg[fn][0] = lorg[fn][1] = lorg[fn][2] = 0.0;
			ldir[fn][0] = ldir[fn][1] = ldir[fn][2] = 0.0;
			ldist[fn] = -1.0;
			errno = EDOM;
		} else {
			pix2loc(loc, &input[fn].rs, xscan, ymax-1-yscan);
			ldist[fn] = viewray(lorg[fn], ldir[fn],
					&input[fn].vw, loc[0], loc[1]);
			if (ldist[fn] < -FTINY)
				errno = EDOM;
		}
		ltick[fn] = eclock;
	}
	if (nam == vray[i=6])
		return(ldist[fn]);
	while (i--)
		if (nam == vray[i])
			return(i < 3 ? lorg[fn][i] : ldir[fn][i-3]);
	eputs("Bad call to l_ray()!\n");
	quit(1);
	return 1; /* pro forma return */
}


static double
l_psize(char *nm)		/* compute pixel size in steradians */
{
	static unsigned long	ltick[MAXINP];
	static double	psize[MAXINP];
	FVECT	dir0, org, dirx, diry;
	RREAL	locx[2], locy[2];
	double	d;
	int	fn;
	int	i;

	d = argument(1);
	if (d <= -0.5 || d >= nfiles+0.5) {
		errno = EDOM;
		return(0.0);
	}
	if (d < 0.5)
		return((double)nfiles);
	fn = d - 0.5;
	if (ltick[fn] != eclock) {		/* need to compute? */
		psize[fn] = 0.0;
		if (input[fn].vw.type == 0)
			errno = EDOM;
		else if (input[fn].vw.type != VT_PAR &&
				funvalue(vray[6], 1, &d) >= -FTINY) {
			for (i = 0; i < 3; i++)
				dir0[i] = funvalue(vray[3+i], 1, &d);
			pix2loc(locx, &input[fn].rs, xscan+1, ymax-1-yscan);
			pix2loc(locy, &input[fn].rs, xscan, ymax-yscan);
			if (viewray(org, dirx, &input[fn].vw,
					locx[0], locx[1]) >= -FTINY &&
					viewray(org, diry, &input[fn].vw,
					locy[0], locy[1]) >= -FTINY) {
						/* approximate solid angle */
				for (i = 0; i < 3; i++) {
					dirx[i] -= dir0[i];
					diry[i] -= dir0[i];
				}
				fcross(dir0, dirx, diry);
				psize[fn] = VLEN(dir0);
			}
		}
		ltick[fn] = eclock;
	}
	return(psize[fn]);
}


extern void
wputs(const char *msg)
{
	if (!nowarn)
		eputs(msg);
}


extern void
eputs(const char *msg)
{
	fputs(msg, stderr);
}


extern void
quit(int code)		/* exit gracefully */
{
	int  i;
				/* close input files */
	for (i = 0; i < nfiles; i++)
		if (input[i].name == Command)
			pclose(input[i].fp);
		else
			fclose(input[i].fp);
	exit(code);
}
