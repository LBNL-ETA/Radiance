#ifndef lint
static const char	RCSid[] = "$Id: aed.c,v 2.6 2004/03/30 16:13:00 schorsch Exp $";
#endif
/*
 *  aed.c - driver for AED 512 terminal.
 */

#include "copyright.h"

#include  <stdio.h>

#include  "rterror.h"
#include  "color.h"
#include  "driver.h"


	/* AED command characters */

#define AEDFMT	"1888N"	/* Format string to send to AED */
#define CSTAT	0100	/* Console status: lower case */
#define SCT	75	/* Set color lookup table */
#define SEC	67	/* Set color for vector drawing */
#define MOV	81	/* Set CAP (current access pointer) to x, y */
#define OPT	40	/* Miscellaneous terminal control */
#define SEN	71	/* Set encoding types */
#define RST	48	/* Reset terminal */
#define FFD	12	/* Clear screen */
#define EJC	85	/* Enable Joystick cursor positioning */
#define DJC	100	/* Disable Joystick cursor positioning */
#define SCC	99	/* Set Cursor Colors */
#define SCP	93	/* Set Cursor Parameters */
#define DFR	111	/* Draw Filled Rectangle */
#define RCP	106	/* Read Cursor Position */
#define ESC	27	/* Escape starts a command sequence */
#define SCS	96	/* Set Console Status */
#define END	1	/* Ctrl-A used to terminate command mode */

#define BLK	0	/* color table entry for black */
#define WHT	7	/* color table entry for white */

#define command(c)	(putc(ESC, stdout), putc(c, stdout))
#define byte(b)		putc(b, stdout)
#define string(s)	fputs(s, stdout)
#define flush()		fflush(stdout)

#define  GAMMA		2.5		/* exponent for color correction */

#define  NCOLORS	248		/* our color table size */
#define  MINPIX		8		/* minimum hardware color */

#define  NCOLS		512		/* maximum # columns for output */
#define  NROWS		483-COMHT	/* maximum # rows for output */
#define  COMHT		16		/* height of command line */
#define  COMCW		63		/* maximum chars on command line */

static void aed_errout(char*);
static void longwait(int t);
static void aedgetcap(int  *xp, int *yp);
static void aedsetcap(int  x, int y);
static dr_newcolrf_t anewcolr;
static void aedcoord(int  x, int y);

static dr_closef_t aed_close;
static dr_clearf_t aed_clear;
static dr_paintrf_t aed_paintr;
static dr_getcurf_t aed_getcur;
static dr_comoutf_t comout;

static struct driver  aed_driver = {
	aed_close, aed_clear, aed_paintr, aed_getcur,
	aed_comout, NULL, NULL,
	1.0, NCOLS, NROWS
};


struct driver *
aed_init(name, id)			/* open AED */
char  *name, *id;
{
	if (ttyset(&aed_driver, fileno(stdin)) < 0) {	/* set tty driver */
		eputs("cannot access terminal\n");
		return(NULL);
	}
	command(RST);					/* reset AED */
	longwait(2);
	command(OPT);
	byte(6); byte(1);				/* AED command set */
	command(SEN);
	string(AEDFMT);
	command(SCS);					/* set console status */
	byte(CSTAT);
	command(SCC);					/* set cursor type */
	byte(BLK); byte(WHT); byte(15);
	command(SCP);
	byte('+'); byte(0); byte(1);
	make_gmap(GAMMA);				/* make color map */
	erract[USER].pf =				/* set error vector */
	erract[SYSTEM].pf =
	erract[INTERNAL].pf =
	erract[CONSISTENCY].pf = aed_errout;
	erract[COMMAND].pf = aed_errout;
	if (erract[WARNING].pf != NULL)
		erract[WARNING].pf = aed_errout;
	return(&aed_driver);
}


static void
aed_close(void)					/* close AED */
{
	erract[USER].pf = 			/* reset error vector */
	erract[SYSTEM].pf =
	erract[INTERNAL].pf =
	erract[CONSISTENCY].pf = eputs;
	erract[COMMAND].pf = NULL;
	if (erract[WARNING].pf != NULL)
		erract[WARNING].pf = wputs;
	aedsetcap(0, 0);			/* go to bottom */
	command(SEC);
	byte(WHT);				/* white text */
	byte(END);
	byte('\n');				/* scroll */
	flush();
	ttydone();
}


static
aed_clear(x, y)					/* clear AED */
int  x, y;
{
	command(FFD);
	new_ctab(NCOLORS);			/* init color table */
	flush();
}


static void
aed_paintr(col, xmin, ymin, xmax, ymax)		/* paint a rectangle */
COLOR  col;
int  xmin, ymin, xmax, ymax;
{
	int  ndx;

	ndx = get_pixel(col, anewcolr);		/* calls anewcolr() */
	command(SEC);				/* draw rectangle */
	byte(ndx+MINPIX);
	aedsetcap(xmin, ymin+COMHT);
	command(DFR);
	aedcoord(xmax-1, ymax+(-1+COMHT));
	flush();
}


static int
aed_getcur(xp, yp)			/* get cursor position */
int  *xp, *yp;
{
	int  com;

	command(EJC);			/* enable cursor */
	com = getch();			/* get key */
	command(DJC);			/* disable cursor */
	aedgetcap(xp, yp);		/* get cursor position */
	*yp -= COMHT;			/* convert coordinates */
	return(com);			/* return key hit */
}


static void
aed_comout(out)				/* output to command line */
register char  *out;
{
	static int  newl = 1;
	static int  inmsg = 0;

	while (*out) {
		if (newl) {			/* new line */
			command(SEC);
			byte(BLK);
			aedsetcap(0, 0);
			command(DFR);
			aedcoord(NCOLS-1, COMHT-1);
			aedsetcap(1, 4);
			command(SEC);		/* white text */
			byte(WHT);
			byte(END);
			newl = 0;
		}
		switch (*out) {
		case '\n':			/* end of line */
			newl = 1;
		/* fall through */
		case '\r':
			inmsg = 0;
			byte('\r');
			break;
		case '\b':			/* back space */
			if (inmsg > 0 && --inmsg < COMCW)
				byte('\b');
			break;
		default:			/* character */
			if (inmsg++ < COMCW)
				byte(*out);
			break;
		}
		out++;
	}
	flush();
}


static void
aed_errout(msg)				/* print an error message */
char  *msg;
{
	aed_comout(msg);
	if (*msg && msg[strlen(msg)-1] == '\n')
		longwait(2);
}


/*
 * aedsetcap - sets AED's current access pointer to (x, y).
 */

static void
aedsetcap(x, y)
register int  x, y;
{
	command(MOV);
	aedcoord(x, y);
}

/*
 * aedcoord - puts out an (x, y) coordinate in AED 8 bit format.
 */

static void
aedcoord(x, y)
register int  x, y;
{
	byte(((x >> 4) & 0x30) | ((y >> 8) & 0x3));
	byte(x & 0xff);
	byte(y & 0xff);
}


static void
aedgetcap(xp, yp)		/* get cursor postion */
int  *xp, *yp;
{
	register int  c;
				/* get cursor postition */
	command(RCP);
	flush();
	c = getch();
	*xp = (c & 0x30) << 4;
	*yp = (c & 0x3) << 8;
	*xp |= getch();
	*yp |= getch();
}


static void
anewcolr(index, r, g, b)		/* enter a color into our table */
int  index;
int  r, g, b;
{
	command(SCT);
	byte((index+MINPIX)&255);
	byte(1);
	byte(r);
	byte(g);
	byte(b);
	flush();
}


static void
longwait(t)		/* longer wait */
int  t;
{
	flush();
	sleep(t);
}
