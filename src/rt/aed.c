#ifndef lint
static const char	RCSid[] = "$Id: aed.c,v 2.3 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 *  aed.c - driver for AED 512 terminal.
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

#include  <stdio.h>

#include  "driver.h"

#include  "color.h"


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

int  aed_close(), aed_clear(), aed_paintr(),
		aed_getcur(), aed_comout(), aed_errout();

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


static
aed_close()					/* close AED */
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


static
aed_paintr(col, xmin, ymin, xmax, ymax)		/* paint a rectangle */
COLOR  col;
int  xmin, ymin, xmax, ymax;
{
	extern int  anewcolr();
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


static
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


static
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

static
aedsetcap(x, y)
register int  x, y;
{
	command(MOV);
	aedcoord(x, y);
}

/*
 * aedcoord - puts out an (x, y) coordinate in AED 8 bit format.
 */

static
aedcoord(x, y)
register int  x, y;
{
	byte(((x >> 4) & 0x30) | ((y >> 8) & 0x3));
	byte(x & 0xff);
	byte(y & 0xff);
}


static
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


static
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


static
longwait(t)		/* longer wait */
int  t;
{
	flush();
	sleep(t);
}
