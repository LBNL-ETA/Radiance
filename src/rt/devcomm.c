#ifndef lint
static const char	RCSid[] = "$Id: devcomm.c,v 2.6 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 *  devcomm.c - communication routines for separate drivers.
 *
 *  External symbols declared in driver.h
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

#include "vfork.h"

#ifndef DEVPATH
#define DEVPATH		getenv("PATH")	/* device search path */
#endif

static int	comm_getcur();
static void	comm_close(), comm_clear(), comm_paintr(),
		comm_comin(), comm_comout(), comm_flush();

struct driver	comm_driver = {
	comm_close, comm_clear, comm_paintr, comm_getcur,
	comm_comout, comm_comin, comm_flush
};

static void	mygets(), myputs(), reply_error(), getstate();

FILE	*devin, *devout;

int	devchild;


static struct driver *
final_connect()				/* verify and initialize connection */
{
	putw(COM_SENDM, devout);
	fflush(devout);
	if (getw(devin) != COM_RECVM)
		return(NULL);
						/* get driver parameters */
	getstate();
						/* set error vectors */
	erract[COMMAND].pf = comm_comout;
	if (erract[WARNING].pf != NULL)
		erract[WARNING].pf = comm_comout;
	return(&comm_driver);
}


struct driver *
slave_init(dname, id)			/* run rview in slave mode */
char	*dname, *id;
{
	devchild = -1;				/* we're the slave here */
	devout = stdout;			/* use standard input */
	devin = stdin;				/* and standard output */
	return(final_connect());		/* verify initialization */
}


struct driver *
comm_init(dname, id)			/* set up and execute driver */
char	*dname, *id;
{
	char	*dvcname;
	int	p1[2], p2[2];
	char	pin[16], pout[16];
						/* find driver program */
	if ((dvcname = getpath(dname, DEVPATH, X_OK)) == NULL) {
		eputs(dname);
		eputs(": not found\n");
		return(NULL);
	}
						/* open communication pipes */
	if (pipe(p1) == -1 || pipe(p2) == -1)
		goto syserr;
	if ((devchild = vfork()) == 0) {	/* fork driver process */
		close(p1[1]);
		close(p2[0]);
		sprintf(pin, "%d", p1[0]);
		sprintf(pout, "%d", p2[1]);
		execl(dvcname, dname, pin, pout, id, 0);
		perror(dvcname);
		_exit(127);
	}
	if (devchild == -1)
		goto syserr;
	close(p1[0]);
	close(p2[1]);
	if ((devout = fdopen(p1[1], "w")) == NULL)
		goto syserr;
	if ((devin = fdopen(p2[0], "r")) == NULL)
		goto syserr;
	return(final_connect());		/* verify initialization */
syserr:
	perror(dname);
	return(NULL);
}


static void
comm_close()			/* done with driver */
{
	int	pid;

	erract[COMMAND].pf = NULL;		/* reset error vectors */
	if (erract[WARNING].pf != NULL)
		erract[WARNING].pf = wputs;
	fclose(devout);
	fclose(devin);
	if (devchild < 0)
		return;
	while ((pid = wait(0)) != -1 && pid != devchild)
		;
}


static void
comm_clear(xres, yres)				/* clear screen */
int	xres, yres;
{
	putc(COM_CLEAR, devout);
	putw(xres, devout);
	putw(yres, devout);
	fflush(devout);
}


static void
comm_paintr(col, xmin, ymin, xmax, ymax)	/* paint a rectangle */
COLOR	col;
int	xmin, ymin, xmax, ymax;
{
	putc(COM_PAINTR, devout);
	fwrite((char *)col, sizeof(COLOR), 1, devout);
	putw(xmin, devout);
	putw(ymin, devout);
	putw(xmax, devout);
	putw(ymax, devout);
}


static void
comm_flush()				/* flush output to driver */
{
	putc(COM_FLUSH, devout);
	fflush(devout);
	if (getc(devin) != COM_FLUSH)
		reply_error("flush");
	getstate();
}


static int
comm_getcur(xp, yp)			/* get and return cursor position */
int	*xp, *yp;
{
	int	c;

	putc(COM_GETCUR, devout);
	fflush(devout);
	if (getc(devin) != COM_GETCUR)
		reply_error("getcur");
	c = getc(devin);
	*xp = getw(devin);
	*yp = getw(devin);
	return(c);
}


static void
comm_comout(str)			/* print string to command line */
char	*str;
{
	putc(COM_COMOUT, devout);
	myputs(str, devout);
	if (str[strlen(str)-1] == '\n')
		fflush(devout);
}


static void
comm_comin(buf, prompt)			/* read string from command line */
char	*buf;
char	*prompt;
{
	putc(COM_COMIN, devout);
	if (prompt == NULL)
		putc(0, devout);
	else {
		putc(1, devout);
		myputs(prompt, devout);
	}
	fflush(devout);
	if (getc(devin) != COM_COMIN)
		reply_error("comin");
	mygets(buf, devin);
	getstate();
}


static void
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


static void
myputs(s, fp)				/* put string to file (with nul) */
register char	*s;
register FILE	*fp;
{
	do
		putc(*s, fp);
	while (*s++);
}


static void
reply_error(routine)			/* what should we do here? */
char	*routine;
{
	eputs(routine);
	eputs(": driver reply error\n");
	quit(1);
}


static void
getstate()				/* get driver state variables */
{
	fread((char *)&comm_driver.pixaspect,
			sizeof(comm_driver.pixaspect), 1, devin);
	comm_driver.xsiz = getw(devin);
	comm_driver.ysiz = getw(devin);
	comm_driver.inpready = getw(devin);
}
