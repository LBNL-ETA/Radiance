#ifndef lint
static const char	RCSid[] = "$Id: aed5.c,v 1.2 2003/07/01 21:21:40 greg Exp $";
#endif
/*
 *  aed.c - driver for AED 512 terminal.
 *
 *     7/8/86
 */

#include  <fcntl.h>

#include  "meta.h"


	/* AED command characters */

#define AEDFMT	"1888N"	/* Format string to send to AED */
#define SCT	75	/* Set color lookup table */
#define SEC	67	/* Set color for vector drawing */
#define SBC	91	/* Set background color */
#define MOV	81	/* Set CAP (current access pointer) to x, y */
#define MVR	105	/* Set CAP relative */
#define DVA	65	/* Draw vector to new CAP x, y */
#define DVR	108	/* Draw vector relative */
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
#define SAP	94	/* Set alphanumeric parameters */
#define SLS	49	/* Set Line Style */
#define DFR	111	/* Draw Filled Rectangle */
#define MAXRLEN	255	/* Maximum runlength for a span */
#define NUL	0	/* Null terminates run sequences */
#define ESC	27	/* Escape starts a command sequence */
#define END	1	/* Ctrl-A used to terminate command mode */

#define  command(c)	(putc(ESC, stdout), putc(c, stdout))
#define  byte(b)	putc(b, stdout)
#define  ascii(c)	putc(c, stdout)

#define BKGCOL	7	/* Background color (white) */

#define  NROWS		512		/* # rows for output */
#define  NCOLS		512		/* # columns for output */

#define  xconv(x)	CONV(x, NCOLS)
#define  yconv(y)	CONV(y, NROWS)

#define  XCOM  "pexpand +vOCIs -tpSUR %s"

char  *progname;

static short  newpage = TRUE; 

static int  curx = -1,		/* current position */
	    cury = -1;

static int  curcol = -1;	/* current color */

static short  curmod = -1;	/* current line drawing mode */

static short  cmode[4] = {0, 1, 2, 4};		/* color map */

static short  lmode[4] = {255, 15, 85, 39};	/* line map */

static PRIMITIVE  nextp;



main(argc, argv)
int  argc;
char  **argv;
{
	FILE  *fp, *popen();
	char  comargs[200], shcom[300];
	short  condonly = FALSE, conditioned = FALSE;

	progname = *argv++;
	argc--;

	while (argc && **argv == '-')  {
		switch (*(*argv+1))  {
		case 'c':
			condonly = TRUE;
			break;
		case 'r':
			conditioned = TRUE;
			break;
		default:
			error(WARNING, "unknown option");
			break;
		}
		argv++;
		argc--;
	}

	if (conditioned) {
		init();
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
	} else  {
		comargs[0] = '\0';
		while (argc)  {
			strcat(comargs, " ");
			strcat(comargs, *argv);
			argv++;
			argc--;
		}
		sprintf(shcom, XCOM, comargs);
		if (condonly)
			return(system(shcom));
		else  {
			init();
			if ((fp = popen(shcom, "r")) == NULL)
				error(SYSTEM, "cannot execute input filter");
			plot(fp);
			return(pclose(fp));
		}
	}
	return(0);
}





plot(infp)		/* plot meta-file */

FILE  *infp;

{

	do {
		readp(&nextp, infp);
		while (isprim(nextp.com)) {
			doprim(&nextp);
			fargs(&nextp);
			readp(&nextp, infp);
		}
		doglobal(&nextp);
		fargs(&nextp);
	} 
	while (nextp.com != PEOF);

}





doglobal(g)			/* execute a global command */

PRIMITIVE  *g;

{
	int  tty;
	char  c;

	switch (g->com) {

	case PEOF:
		break;

	case PDRAW:
		fflush(stdout);
		break;

	case PEOP:
		newpage = TRUE;
		if (!isatty(fileno(stdout)))
			break;
		/* fall through */

	case PPAUS:
		fflush(stdout);
		tty = open(TTY, O_RDWR);
		if (g->args != NULL) {
			write(tty, g->args, strlen(g->args));
			write(tty, " - (hit return to continue)", 27);
		} 
		else
			write(tty, "\007", 1);
		do {
			c = '\n';
			read(tty, &c, 1);
		} 
		while (c != '\n' && c != '\r');
		close(tty);
		break;

	default:
		sprintf(errmsg, "unknown command '%c' in doglobal", g->com);
		error(WARNING, errmsg);
		break;
	}

}




doprim(p)		/* plot primitive */

register PRIMITIVE  *p;

{

	if (newpage) {				/* clear screen */
		command(FFD);
		shortwait(50);
		newpage = FALSE;
	}

	switch (p->com) {

	case PLSEG:
		plotlseg(p);
		break;

	case PMSTR:
		printstr(p);
		break;

	case PRFILL:
		fillrect(p);
		break;

	default:
		sprintf(errmsg, "unknown command '%c' in doprim", p->com);
		error(WARNING, errmsg);
		return;
	}

}






printstr(p)		/* output a string */

register PRIMITIVE  *p;

{
	static int  hsp[4] = {6, 5, 4, 3};
	char  font, size;
	int  hspace, vspace = 18;

	hspace = hsp[(p->arg0 >> 4) & 03];

	if (p->arg0 & 04)
		hspace *= 2;

	if (hspace > 16) {
		font = '7'; size = '2';
	} else if (hspace > 12) {
		font = '5'; size = '2';
	} else if (hspace > 8) {
		font = '7'; size = '1';
	} else {
		font = '5'; size = '1';
	}

	command(SAP);
	ascii(size); ascii(font);
	byte(hspace); byte(vspace);
	ascii('L');

	setcolor(p->arg0 & 03);

	move(p->xy[XMN], p->xy[YMX]);

	putc(END, stdout);
	fputs(p->args, stdout);

	curx = -1;
	cury = -1;
}





plotlseg(p)		/* plot a line segment */

register PRIMITIVE  *p;

{
	static short  right = FALSE;
	int  y1, y2;
	short  lm = (p->arg0 >> 4) & 03;

	if (p->arg0 & 0100) {
		y1 = p->xy[YMX];
		y2 = p->xy[YMN];
	} 
	else {
		y1 = p->xy[YMN];
		y2 = p->xy[YMX];
	}

	setcolor(p->arg0 & 03);

	if (lm != curmod) {
		command(SLS);
		byte(lmode[lm]); byte(85);
		curmod = lm;
	}

	if (p->xy[XMN] == curx && y1 == cury)
		draw(p->xy[XMX], y2);
	else if (p->xy[XMX] == curx && y2 == cury)
		draw(p->xy[XMN], y1);
	else if (right = !right) {
		move(p->xy[XMN], y1);
		draw(p->xy[XMX], y2);
	} else {
		move(p->xy[XMX], y2);
		draw(p->xy[XMN], y1);
	}
}


fillrect(p)

register PRIMITIVE  *p;

{

	setcolor(p->arg0 & 03);
	move(p->xy[XMN], p->xy[YMN]);
	curx = xconv(p->xy[XMX]);
	cury = yconv(p->xy[YMX]);
	command(DFR);
	aedcoord(curx, cury);

}

init()			/* initialize terminal */
{
		/* Reset AED and tell it the data format */
	command(RST);
	longwait(2);
	command(OPT);
	byte(6); byte(1);			/* AED command set */
	command(SEN); fputs(AEDFMT, stdout);
	command(SBC); byte(BKGCOL);
}


setcolor(cn)		/* set color */
int  cn;
{
	if (cn != curcol) {
		command(SEC);
		byte(cmode[cn]);
		curcol = cn;
	}
}


move(x, y)		/* move to coordinate (x,y) */
int  x, y;
{
	command(MOV);
	aedcoord(xconv(x), yconv(y));
	curx = x;
	cury = y;
}


draw(x, y)		/* draw vector from CAP to (x,y) */
int  x, y;
{
	command(DVA);
	aedcoord(xconv(x), yconv(y));
	curx = x;
	cury = y;
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


longwait(t)		/* longer wait */
int  t;
{
	fflush(stdout);
	sleep(t);
}


shortwait(t)		/* shorter wait */
int  t;
{
	register long  l = t*1000;

	fflush(stdout);
	while (l--)
		;
}
