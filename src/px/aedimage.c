#ifndef lint
static const char	RCSid[] = "$Id: aedimage.c,v 2.5 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 *  aedimage.c - RADIANCE driver for AED 512 terminal.
 *
 *     3/20/86
 *     3/3/88	Added calls to Paul Heckbert's ciq routines
 */

#include  <stdio.h>

#include  <math.h>

#include  <signal.h>

#include  <sys/ioctl.h>

#include  "pic.h"

#include  "color.h"


	/* AED command characters */

#define AEDFMT	"1888N"	/* Format string to send to AED */
#define CSTAT	0100	/* Console status: lower case */
#define SCS	96	/* Set Console Status */
#define SCT	75	/* Set color lookup table */
#define SEC	67	/* Set color for vector drawing */
#define MOV	81	/* Set CAP (current access pointer) to x, y */
#define MVR	105	/* Set CAP relative */
#define DVA	65	/* Draw vector to new CAP x, y */
#define DVR	108	/* Draw vector relative */
#define WHS	88	/* Write horizontal scan */
#define WHR	92	/* Write horizontal runs */
#define OPT	40	/* Miscellaneous terminal control */
#define SEN	71	/* Set encoding types */
#define RST	48	/* Reset terminal */
#define FFD	12	/* Clear screen */
#define SCS	96	/* Set status */
#define DAI	114	/* Define area of interest */
#define EJC	85	/* Enable Joystick cursor positioning */
#define DJC	100	/* Disable Joystick cursor positioning */
#define SCC	99	/* Set Cursor Colors */
#define SCP	93	/* Set Cursor Parameters */
#define RCP	106	/* Read Cursor Position */
#define MAXRLEN	255	/* Maximum runlength for a span */
#define NUL	0	/* Null terminates run sequences */
#define ESC	27	/* Escape starts a command sequence */
#define END	1	/* Ctrl-A used to terminate command mode */

#define BLK	0	/* color table entry for black */
#define WHT	7	/* color table entry for white */

#define command(c)	(putc(ESC, stdout), putc(c, stdout))
#define byte(b)		putc(b, stdout)
#define flush()		fflush(stdout)

#define  GAMMA		2.5		/* gamma value used in correction */

#define  MINCOLOR	8		/* start of device color table */

#define  NCOLORS	248		/* available color table size */

#define  NROWS		512		/* maximum # rows for output */
#define  NCOLS		512		/* maximum # columns for output */

int  dither = 0;			/* dither colors? */

int  greyscale = 0;			/* grey only? */

int  initialized = 0;
struct	sgttyb	oldflags;
int  oldlmode;

int  intrcnt = 0;			/* interrupt well */

#define  disintr()	(intrcnt--)	/* disable interrupts */
#define  enaintr()	onintr()	/* enable interrupts */

char  *progname;

char  errmsg[128];

COLR  scanline[NCOLS];

colormap  colrtbl;

int  cury;

int  xmax, ymax;

FILE  *fin;

extern long  ftell();
long  scanpos[NROWS];

double  exposure = 1.0;
int  wrong_fmt = 0;


main(argc, argv)
int  argc;
char  *argv[];
{
	int  onintr(), checkhead();
	char  sbuf[256];
	register int  i;
	
	progname = argv[0];

#ifdef  SIGTSTP
	signal(SIGTSTP, SIG_IGN);
#endif
	if (signal(SIGINT, onintr) == SIG_IGN)
		signal(SIGINT, SIG_IGN);

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'd':
				dither = !dither;
				break;
			case 'b':
				greyscale = !greyscale;
				break;
			default:
				goto userr;
			}
		else
			break;

	if (i+1 != argc)
		goto userr;

	if ((fin = fopen(argv[i], "r")) == NULL) {
		sprintf(errmsg, "can't open file \"%s\"", argv[i]);
		quitmsg(errmsg);
	}
				/* get header */
	getheader(fin, checkhead, NULL);
	if (wrong_fmt)
		quitmsg("input must be a Radiance picture");
				/* get picture dimensions */
	if (fgetresolu(&xmax, &ymax, fin) < 0)
		quitmsg("bad picture size");
	if (xmax > NCOLS || ymax > NROWS)
		quitmsg("resolution mismatch");

	init();			/* initialize terminal for output */

	for (i = 0; i < NROWS; i++)
		scanpos[i] = -1;

	if (greyscale)		/* map colors and display */
		biq(dither,NCOLORS,1,colrtbl);
	else
		ciq(dither,NCOLORS,1,colrtbl);

	loopcom();

	quitmsg(NULL);
userr:
	fprintf(stderr, "Usage: %s [-d][-b] input\n", progname);
	quit(1);
}


int
checkhead(line)				/* deal with line from header */
char  *line;
{
	char	fmt[32];

	if (isexpos(line))
		exposure *= exposval(line);
	else if (isformat(line)) {
		formatval(fmt, line);
		wrong_fmt = strcmp(fmt, COLRFMT);
	}
	return(0);
}


init()			/* initialize terminal */
{
	struct sgttyb  flags;
	int  lmode;

	if (isatty(1)) {
				/* Set literal mode so AED gets bit 8 */
		ioctl(1, TIOCLGET, &oldlmode);
		lmode = oldlmode | LLITOUT | LNOFLSH;
		ioctl(1, TIOCLSET, &lmode);
		ioctl(1, TIOCLGET, &lmode);	/* 4.2 bug */
				/* Turn terminal echoing off */
		ioctl(1, TIOCGETP, &oldflags);
		ioctl(1, TIOCGETP, &flags);
		flags.sg_flags &= ~ECHO;
		ioctl(1, TIOCSETP, &flags);
		initialized = 1;
	} else
		quitmsg("output must be to a terminal");

		/* Reset AED and tell it the data format */
	command(RST);
	longwait(2);
	command(OPT);
	byte(6); byte(1);			/* AED command set */
	command(SEN); fputs(AEDFMT, stdout);
	command(SCS); byte(CSTAT);		/* console status */
	command(FFD);				/* clear screen */
	flush();
}


onintr()		/* die gracefully */
{
	if (++intrcnt > 0)
		quitmsg("interrupt");
}


quitmsg(err)		/* uninitialize terminal and exit */
char  *err;
{
	if (initialized) {		/* close AED */
		command(DJC);		/* disable cursor */
		command(SEC);
		byte(WHT);		/* Set color to white */
		aedsetcap(0, 0);	/* Put cursor in a good place */
		byte(END);
		flush();
		ioctl(1, TIOCSETP, &oldflags);
		ioctl(1, TIOCLSET, &oldlmode);
		putchar('\n');
	}
	if (err != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	exit(0);
}


void
eputs(s)
char  *s;
{
	fputs(s, stderr);
}


void
quit(status)
int  status;
{
	quitmsg(status ? "aborted" : NULL);
}


loopcom()				/* print pixel values interactively */
{
	int  coffset[3];		/* color table offset */
	struct sgttyb  flags;
	COLOR  cval, ctmp;
	int  rept, x, y, j;
	register int  com, i;
						/* Set raw mode on input */
	ioctl(0, TIOCGETP, &flags);
	flags.sg_flags |= RAW;
	flags.sg_flags &= ~ECHO;
	ioctl(0, TIOCSETP, &flags);

	coffset[0] = coffset[1] = coffset[2] = 0;
						/* set cursor type */
	command(SCC);
	byte(BLK); byte(WHT); byte(15);
	command(SCP);
	byte('+'); byte(0); byte(1);
						/* loop on command */
	for ( ; ; ) {
		command(EJC);			/* enable cursor */
		flush();
		rept = 0;			/* get command */
		while ((com = getc(stdin)) >= '0' && com <= '9')
			rept = rept*10 + com - '0';
		if (rept == 0) rept = 1;
		command(DJC);			/* disable cursor */
		switch (com) {
		case '\r':
		case '\n':
		case 'l':
		case 'L':
		case 'c':
		case 'C':
		case 'p':
		case 'P':
		case '=':
			aedgetcap(&x, &y);	/* get cursor position */
			y = ymax-1-y;
			if (y < 0 || scanpos[y] == -1)	/* off picture? */
				break;
						/* get value */
			setcolor(cval, 0.0, 0.0, 0.0);
			for (j = rept-1; j >= 0; j--) {
				if (y+j >= NROWS || scanpos[y+j] == -1)
					continue;
				getscan(y+j);
				for (i = 0; i < rept; i++)
					if (x+i < xmax) {
						colr_color(ctmp, scanline[x+i]);
						addcolor(cval, ctmp);
					}
			}
			scalecolor(cval, 1.0/(rept*rept));
						/* print value */
			command(SEC);
			byte(scanline[(x+20)%xmax][EXP] >= COLXS ? BLK : WHT);
			byte(END);
			switch (com) {
			case '\r':
			case '\n':
				printf("%-3g", intens(cval)/exposure);
				break;
			case 'l':
			case 'L':
				printf("%-3gL", luminance(cval)/exposure);
				break;
			case 'c':
			case 'C':
				printf("(%.2f,%.2f,%.2f)",
						colval(cval,RED),
						colval(cval,GRN),
						colval(cval,BLU));
				break;
			case 'p':
			case 'P':
				printf("(%d,%d)", x, y);
				break;
			case '=':
				printf("(%d,%d,%d)", coffset[0],
						coffset[1], coffset[2]);
				break;
			}
			break;
		case 'w':
			rept = -rept;
		case 'W':
			coffset[0] += rept;
			coffset[1] += rept;
			coffset[2] += rept;
			movecolor(coffset);
			break;
		case 'r':
			rept = -rept;
		case 'R':
			coffset[0] += rept;
			movecolor(coffset);
			break;
		case 'g':
			rept = -rept;
		case 'G':
			coffset[1] += rept;
			movecolor(coffset);
			break;
		case 'b':
			rept = -rept;
		case 'B':
			coffset[2] += rept;
			movecolor(coffset);
			break;
		case '!':
			coffset[0] = coffset[1] = coffset[2] = 0;
			movecolor(coffset);
			break;
		case 'C'-'@':
		case 'q':
		case 'Q':
			return;
		}
	}
}


getscan(y)
int  y;
{
	if (y != cury) {
		if (scanpos[y] == -1)
			quitmsg("cannot seek in getscan");
		if (fseek(fin, scanpos[y], 0) == -1)
			quitmsg("fseek error");
		cury = y;
	} else
		scanpos[y] = ftell(fin);

	if (freadcolrs(scanline, xmax, fin) < 0)
		quitmsg("read error");

	cury++;
}


picreadline3(y, l3)			/* read in 3-byte scanline */
int  y;
register rgbpixel  *l3;
{
	register int	i;
							/* read scanline */
	getscan(y);
							/* convert scanline */
	normcolrs(scanline, xmax, 0);
	for (i = 0; i < xmax; i++) {
		l3[i].r = scanline[i][RED];
		l3[i].g = scanline[i][GRN];
		l3[i].b = scanline[i][BLU];
	}
}


picwriteline(y, l)		/* output scan line of pixel values */
int  y;
register pixel  *l;
{
	int  curpix;
	register int  n, ncols;
	
	disintr();			/* disable interrupts */
	aedsetcap(0, ymax-1-y);
	command(DAI);
	aedcoord(xmax-1, ymax-1-y);
	ncols = xmax;
	if (dither) {		/* pixel run lengths short */
		command(WHS);
		while (ncols-- > 0)
			byte(MINCOLOR + *l++);
	} else {			/* pixel run lengths long */
		command(WHR);
		while (ncols-- > 0) {
			curpix = *l++;
			for (n = 1; n < MAXRLEN; n++) {
				if (ncols <= 0 || *l != curpix)
					break;
				l++; ncols--;
			}
			byte(n);
			byte(MINCOLOR + curpix);
		}
		byte(NUL);
	}
	flush();
	enaintr();			/* enable interrupts */
}


/*
 * aedsetcap - sets AED's current access pointer to (x, y). Used
 * before writing a span to put it in the right place.
 */

aedsetcap(x, y)
register int  x, y;
{
	command(MOV);
	aedcoord(x, y);
}

/*
 * aedcoord - puts out an (x, y) coordinate in AED 8 bit format.
 */

aedcoord(x, y)
register int  x, y;
{
	putc(((x >> 4) & 0x30) | ((y >> 8) & 0x3), stdout);
	putc(x & 0xff, stdout);
	putc(y & 0xff, stdout);
}


aedgetcap(xp, yp)		/* get cursor postion */
int  *xp, *yp;
{
	register int  c;
				/* get cursor postition */
	command(RCP);
	flush();
	c = getc(stdin);
	*xp = (c & 0x30) << 4;
	*yp = (c & 0x3) << 8;
	*xp |= getc(stdin);
	*yp |= getc(stdin);
}


movecolor(offset)		/* move our color table up or down */
int  offset[3];
{
	register int  i, j, c;

	disintr();
	command(SCT);
	byte(MINCOLOR);
	byte(NCOLORS);
	for (i = 0; i < NCOLORS; i++)
		for (j = 0; j < 3; j++) {
			c = colrtbl[j][i] + offset[j];
			if (c < 0)
				byte(0);
			else if (c > 255)
				byte(255);
			else
				byte(c);
		}
	shortwait(20);
	enaintr();
}


picreadcm(map)			/* do gamma correction */
colormap  map;
{
	register int  i, val;

	for (i = 0; i < 256; i++) {
		val = pow((i+0.5)/256.0, 1.0/GAMMA) * 256.0;
		map[0][i] = map[1][i] = map[2][i] = val;
	}
}


picwritecm(map)			/* write out color map */
colormap  map;
{
	register int  i, j;

	disintr();
	command(SCT);
	byte(MINCOLOR);
	byte(NCOLORS);
	for (i = 0; i < NCOLORS; i++)
		for (j = 0; j < 3; j++)
			byte(map[j][i]);
	shortwait(20);
	enaintr();
}


longwait(t)		/* longer wait */
int  t;
{
	flush();
	sleep(t);
}


shortwait(t)		/* shorter wait */
int  t;
{
	register long  l = t*1000;

	flush();
	while (l--)
		;
}
