#ifndef lint
static const char	RCSid[] = "$Id: devmain.c,v 2.3 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 *  devmain.c - main for independent drivers.
 *
 *	Redefine your initialization routine to dinit.
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

#include "standard.h"

#include "color.h"

#include "driver.h"

struct driver	*dev = NULL;			/* output device */

FILE	*devin, *devout;			/* communications channels */

char	*progname;				/* driver name */

int	r_clear(), r_paintr(), r_getcur(), r_comout(), r_comin(), r_flush();

int	(*dev_func[NREQUESTS])() = {		/* request handlers */
		r_clear, r_paintr,
		r_getcur, r_comout,
		r_comin, r_flush
	};


main(argc, argv)		/* set up communications and main loop */
int	argc;
char	*argv[];
{
	extern struct driver	*dinit();
	int  com;
						/* set up I/O */
	progname = argv[0];
	if (argc < 3) {
		eputs("arg count\n");
		quit(1);
	}
	devin = fdopen(atoi(argv[1]), "r");
	devout = fdopen(atoi(argv[2]), "w");
	if (devin == NULL || devout == NULL || getw(devin) != COM_SENDM) {
		eputs("connection failure\n");
		quit(1);
	}
						/* open device */
	if ((dev = dinit(argv[0], argv[3])) == NULL)
		quit(1);
	putw(COM_RECVM, devout);		/* verify initialization */
	sendstate();
	fflush(devout);
						/* loop on requests */
	while ((com = getc(devin)) != EOF) {
		if (com >= NREQUESTS || dev_func[com] == NULL) {
			eputs("invalid request\n");
			quit(1);
		}
		(*dev_func[com])();		/* process request */
	}
	quit(0);			/* all done, clean up and exit */
}


void
quit(code)				/* close device and exit */
int	code;
{
	if (dev != NULL)
		(*dev->close)();
	exit(code);
}


r_clear()				/* clear screen */
{
	int	xres, yres;

	xres = getw(devin);
	yres = getw(devin);
	(*dev->clear)(xres, yres);
}


r_paintr()				/* paint a rectangle */
{
	COLOR	col;
	int	xmin, ymin, xmax, ymax;

	fread((char *)col, sizeof(COLOR), 1, devin);
	xmin = getw(devin); ymin = getw(devin);
	xmax = getw(devin); ymax = getw(devin);
	(*dev->paintr)(col, xmin, ymin, xmax, ymax);
}


r_flush()				/* flush output */
{
	if (dev->flush != NULL)
		(*dev->flush)();
	putc(COM_FLUSH, devout);
	sendstate();
	fflush(devout);
}


r_getcur()				/* get and return cursor position */
{
	int	c;
	int	x, y;
					/* get it if we can */
	if (dev->getcur == NULL) {
		c = ABORT;
		x = y = 0;
	} else
		c = (*dev->getcur)(&x, &y);
					/* reply */
	putc(COM_GETCUR, devout);
	putc(c, devout);
	putw(x, devout);
	putw(y, devout);
	fflush(devout);
}


r_comout()				/* print string to command line */
{
	char	str[256];

	mygets(str, devin);
	(*dev->comout)(str);
}


r_comin()				/* read string from command line */
{
	char	buf[256], *prompt;
					/* get prompt */
	if (getc(devin)) {
		mygets(buf, devin);
		prompt = buf;
	} else
		prompt = NULL;
					/* get string */
	(*dev->comin)(buf, prompt);
					/* reply */
	putc(COM_COMIN, devout);
	myputs(buf, devout);
	sendstate();
	fflush(devout);
}


mygets(s, fp)				/* get string from file (with nul) */
register char	*s;
register FILE	*fp;
{
	register int	c;

	while ((c = getc(fp)) != EOF)
		if ((*s++ = c) == '\0')
			return;
	*s = '\0';
}


myputs(s, fp)				/* put string to file (with nul) */
register char	*s;
register FILE	*fp;
{
	do
		putc(*s, fp);
	while (*s++);
}


void
eputs(s)				/* put string to stderr */
register char  *s;
{
	static int  midline = 0;

	if (!*s)
		return;
	if (!midline++) {
		fputs(progname, stderr);
		fputs(": ", stderr);
	}
	fputs(s, stderr);
	if (s[strlen(s)-1] == '\n') {
		fflush(stderr);
		midline = 0;
	}
}


sendstate()				/* send driver state variables */
{
	fwrite((char *)&dev->pixaspect, sizeof(dev->pixaspect), 1, devout);
	putw(dev->xsiz, devout);
	putw(dev->ysiz, devout);
	putw(dev->inpready, devout);
}
