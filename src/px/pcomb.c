#ifndef lint
static const char	RCSid[] = "$Id: pcomb.c,v 2.28 2003/10/22 02:06:35 greg Exp $";
#endif
/*
 *  Combine picture files according to calcomp functions.
 *
 *	1/4/89
 */

#include "platform.h"
#include "standard.h"
#include "color.h"
#include "calcomp.h"
#include "view.h"

#define MAXINP		32		/* maximum number of input files */
#define WINSIZ		64		/* scanline window size */
#define MIDSCN		((WINSIZ-1)/2+1)

struct {
	char	*name;		/* file or command name */
	FILE	*fp;		/* stream pointer */
	VIEW	vw;		/* view for picture */
	RESOLU	rs;		/* image resolution and orientation */
	float	pa;		/* pixel aspect ratio */
	COLOR	*scan[WINSIZ];	/* input scanline window */
	COLOR	coef;		/* coefficient */
	COLOR	expos;		/* recorded exposure */
}	input[MAXINP];			/* input pictures */

int	nfiles;				/* number of input files */

char	ourfmt[LPICFMT+1] = PICFMT;	/* input picture format */

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

char	*progname;			/* global argv[0] */

int	wrongformat = 0;
int	gotview;

FILE	*popen();

extern char	*emalloc();


main(argc, argv)
int	argc;
char	*argv[];
{
	int	original;
	double	f;
	int	a, i;
	SET_DEFAULT_BINARY();
	SET_FILE_BINARY(stdin);
	SET_FILE_BINARY(stdout);
	progname = argv[0];
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
			case 'f':
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
				input[nfiles].name = "<stdin>";
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
		checkfile();
		if (original) {
			colval(input[nfiles].coef,RED) /=
					colval(input[nfiles].expos,RED);
			colval(input[nfiles].coef,GRN) /=
					colval(input[nfiles].expos,GRN);
			colval(input[nfiles].coef,BLU) /=
					colval(input[nfiles].expos,BLU);
			setcolor(input[nfiles].expos, 1.0, 1.0, 1.0);
		}
		nfiles++;
		original = 0;
	}
	init();				/* set constants */
					/* go back and get expressions */
	for (a = 1; a < argc; a++) {
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
			case 'f':
				fcompile(argv[++a]);
				continue;
			case 'e':
				scompile(argv[++a], NULL, 0);
				continue;
			}
		break;
	}
						/* set/get output resolution */
	if (!vardefined(vxres))
		varset(vxres, ':', (double)xmax);
	if (!vardefined(vyres))
		varset(vyres, ':', (double)ymax);
	xres = varvalue(vxres) + .5;
	yres = varvalue(vyres) + .5;
	if (xres <= 0 || yres <= 0) {
		eputs(argv[0]);
		eputs(": illegal output resolution\n");
		quit(1);
	}
						/* complete header */
	printargs(argc, argv, stdout);
	if (strcmp(ourfmt, PICFMT))
		fputformat(ourfmt, stdout);	/* print format if known */
	putchar('\n');
	fprtresolu(xres, yres, stdout);
						/* combine pictures */
	combine();
	quit(0);
usage:
	eputs("Usage: ");
	eputs(argv[0]);
	eputs(
" [-w][-x xr][-y yr][-e expr][-f file] [ [-o][-s f][-c r g b] pic ..]\n");
	quit(1);
}


tabputs(s)			/* put out string preceded by a tab */
char	*s;
{
	char	fmt[32];
	double	d;
	COLOR	ctmp;

	if (isheadid(s))			/* header id */
		return(0);	/* don't echo */
	if (formatval(fmt, s)) {		/* check format */
		if (globmatch(ourfmt, fmt)) {
			wrongformat = 0;
			strcpy(ourfmt, fmt);
		} else
			wrongformat = globmatch(PICFMT, fmt) ? 1 : -1;
		return(0);	/* don't echo */
	}
	if (isexpos(s)) {			/* exposure */
		d = exposval(s);
		scalecolor(input[nfiles].expos, d);
	} else if (iscolcor(s)) {		/* color correction */
		colcorval(ctmp, s);
		multcolor(input[nfiles].expos, ctmp);
	} else if (isaspect(s))
		input[nfiles].pa *= aspectval(s);
	else if (isview(s) && sscanview(&input[nfiles].vw, s) > 0)
		gotview++;
						/* echo line */
	putchar('\t');
	return(fputs(s, stdout));
}


checkfile()			/* ready a file */
{
	register int	i;
					/* process header */
	gotview = 0;
	fputs(input[nfiles].name, stdout);
	fputs(":\n", stdout);
	getheader(input[nfiles].fp, tabputs, NULL);
	if (wrongformat < 0) {
		eputs(input[nfiles].name);
		eputs(": not a Radiance picture\n");
		quit(1);
	}
	if (wrongformat > 0) {
		wputs(input[nfiles].name);
		wputs(": warning -- incompatible picture format\n");
	}
	if (!gotview || setview(&input[nfiles].vw) != NULL)
		input[nfiles].vw.type = 0;
	if (!fgetsresolu(&input[nfiles].rs, input[nfiles].fp)) {
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


double	(*ourbright)() = rgb_bright;


init()					/* perform final setup */
{
	double	l_colin(char *), l_expos(char *), l_pixaspect(char *),
			l_ray(char *), l_psize(char *);
	register int	i;
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
}


combine()			/* combine pictures */
{
	EPNODE	*coldef[3], *brtdef;
	COLOR	*scanout;
	double	d;
	register int	i, j;
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
						/* allocate scanline */
	scanout = (COLOR *)emalloc(xres*sizeof(COLOR));
						/* set input position */
	yscan = ymax+MIDSCN;
						/* combine files */
	for (ypos = yres-1; ypos >= 0; ypos--) {
	    advance();
	    varset(vypos, '=', (double)ypos);
	    for (xpos = 0; xpos < xres; xpos++) {
		xscan = (long)xpos*xmax/xres;
		varset(vxpos, '=', (double)xpos);
		eclock++;
		if (brtdef != NULL) {
		    d = evalue(brtdef);
		    if (d < 0.0)
			d = 0.0;
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
			if (d < 0.0)
			    d = 0.0;
			colval(scanout[xpos],j) = d;
		    }
		}
	    }
	    if (fwritescan(scanout, xres, stdout) < 0) {
		    perror("write error");
		    quit(1);
	    }
	}
	efree((char *)scanout);
}


advance()			/* read in data for next scanline */
{
	int	ytarget;
	register COLOR	*st;
	register int	i, j;

	for (ytarget = (long)ypos*ymax/yres; yscan > ytarget; yscan--)
		for (i = 0; i < nfiles; i++) {
			st = input[i].scan[WINSIZ-1];
			for (j = WINSIZ-1; j > 0; j--)	/* rotate window */
				input[i].scan[j] = input[i].scan[j-1];
			input[i].scan[0] = st;
			if (yscan <= MIDSCN)		/* hit bottom? */
				continue;
			if (freadscan(st, xmax, input[i].fp) < 0) {  /* read */
				eputs(input[i].name);
				eputs(": read error\n");
				quit(1);
			}
			if (fabs(colval(input[i].coef,RED)-1.0) > 1e-3 ||
				fabs(colval(input[i].coef,GRN)-1.0) > 1e-3 ||
				fabs(colval(input[i].coef,BLU)-1.0) > 1e-3)
				for (j = 0; j < xmax; j++)  /* adjust color */
					multcolor(st[j], input[i].coef);
		}
}


double
l_expos(nam)			/* return picture exposure */
register char	*nam;
{
	register int	fn, n;

	fn = argument(1) - .5;
	if (fn < 0 || fn >= nfiles)
		return(1.0);
	if (nam == vbrtexp)
		return((*ourbright)(input[fn].expos));
	n = 3;
	while (n--)
		if (nam == vcolexp[n])
			return(colval(input[fn].expos,n));
	eputs("Bad call to l_expos()!\n");
	quit(1);
}


double
l_pixaspect(char *nm)		/* return pixel aspect ratio */
{
	register int	fn;

	fn = argument(1) - .5;
	if (fn < 0 || fn >= nfiles)
		return(1.0);
	return(input[fn].pa);
}


double
l_colin(nam)			/* return color value for picture */
register char	*nam;
{
	int	fn;
	register int	n, xoff, yoff;
	double	d;

	fn = argument(1) - .5;
	if (fn < 0 || fn >= nfiles) {
		errno = EDOM;
		return(0.0);
	}
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
}


double
l_ray(nam)		/* return ray origin or direction */
register char	*nam;
{
	static unsigned long	ltick[MAXINP];
	static FVECT	lorg[MAXINP], ldir[MAXINP];
	static double	ldist[MAXINP];
	RREAL	loc[2];
	int	fn;
	register int	i;

	fn = argument(1) - .5;
	if (fn < 0 || fn >= nfiles) {
		errno = EDOM;
		return(0.0);
	}
	if (ltick[fn] != eclock) {		/* need to compute? */
		lorg[fn][0] = lorg[fn][1] = lorg[fn][2] = 0.0;
		ldir[fn][0] = ldir[fn][1] = ldir[fn][2] = 0.0;
		ldist[fn] = -1.0;
		if (input[fn].vw.type == 0)
			errno = EDOM;
		else {
			pix2loc(loc, &input[fn].rs, xscan, ymax-1-yscan);
			ldist[fn] = viewray(lorg[fn], ldir[fn],
					&input[fn].vw, loc[0], loc[1]);
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
}


double
l_psize(char *nm)		/* compute pixel size in steradians */
{
	static unsigned long	ltick[MAXINP];
	static double	psize[MAXINP];
	FVECT	dir0, org, dirx, diry;
	RREAL	locx[2], locy[2];
	double	d;
	int	fn;
	register int	i;

	d = argument(1);
	if (d < .5 || d >= nfiles+.5) {
		errno = EDOM;
		return(0.0);
	}
	fn = d - .5;
	if (ltick[fn] != eclock) {		/* need to compute? */
		psize[fn] = 0.0;
		if (input[fn].vw.type == 0)
			errno = EDOM;
		else if (input[fn].vw.type != VT_PAR &&
				funvalue(vray[6], 1, &d) >= 0) {
			for (i = 0; i < 3; i++)
				dir0[i] = funvalue(vray[3+i], 1, &d);
			pix2loc(locx, &input[fn].rs, xscan+1, ymax-1-yscan);
			pix2loc(locy, &input[fn].rs, xscan, ymax-yscan);
			if (viewray(org, dirx, &input[fn].vw,
					locx[0], locx[1]) >= 0 &&
					viewray(org, diry, &input[fn].vw,
					locy[0], locy[1]) >= 0) {
						/* approximate solid angle */
				for (i = 0; i < 3; i++) {
					dirx[i] -= dir0[i];
					diry[i] -= dir0[i];
				}
				fcross(dir0, dirx, diry);
				psize[fn] = sqrt(DOT(dir0,dir0));
			}
		}
		ltick[fn] = eclock;
	}
	return(psize[fn]);
}


void
wputs(msg)
char	*msg;
{
	if (!nowarn)
		eputs(msg);
}


void
eputs(msg)
char	*msg;
{
	fputs(msg, stderr);
}


void
quit(code)		/* exit gracefully */
int  code;
{
	register int  i;
				/* close input files */
	for (i = 0; i < nfiles; i++)
		if (input[i].name == Command)
			pclose(input[i].fp);
		else
			fclose(input[i].fp);
	exit(code);
}
