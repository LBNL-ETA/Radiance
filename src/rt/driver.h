/* RCSid: $Id: driver.h,v 2.4 2003/02/22 02:07:28 greg Exp $ */
/*
 *  driver.h - header file for interactive device drivers.
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

struct driver {				/* driver functions */
	void  (*close)();			/* close device */
	void  (*clear)();			/* clear device */
	void  (*paintr)();			/* paint rectangle */
	int  (*getcur)();			/* get cursor position */
	void  (*comout)();			/* command line output */
	void  (*comin)();			/* command line input */
	void  (*flush)();			/* flush output */
	double  pixaspect;			/* pixel aspect ratio */
	int  xsiz, ysiz;			/* device size */
	int  inpready;				/* input ready on device */
};
					/* magic numbers for verification */
#define COM_SENDM		0x6f37
#define COM_RECVM		0x51da
					/* stream commands */
#define COM_CLEAR		0
#define COM_PAINTR		1
#define COM_GETCUR		2
#define COM_COMOUT		3
#define COM_COMIN		4
#define COM_FLUSH		5
#define NREQUESTS		6	/* number of valid requests */

extern struct device {			/* interactive device */
	char  *name;				/* device name */
	char  *descrip;				/* description */
	struct driver  *(*init)();		/* initialize device */
}  devtable[];				/* supported devices */

extern char  dev_default[];		/* default device name */

#define  MB1		('\n')		/* mouse button 1 */
#define  MB2		('\r')		/* mouse button 2 */
#define  MB3		(' ')		/* mouse button 3 */
#define  ABORT		('C'-'@')	/* abort key */

/*
 *  struct driver *
 *  dname_init(name, id)
 *  char  *name, *id;
 *  {
 *	Initialize device and return pointer to driver
 *	functions.  Returns NULL if an error occurred.
 *	The name string identifies the driver,
 *	and the id string identifies the client.
 *	A device can be open by at most one client.
 *	Be verbose in error reports; call eputs().
 *	If device has its own error output, set erract.
 *  }
 *  (*dev->close)()
 *  {
 *	Close the device.  Reset error vectors.
 *  }
 *  (*dev->clear)(xres, yres)
 *  int  xres, yres;
 *  {
 *	Clear the device for xres by yres output.  This call will
 *	be made prior to any output.
 *  }
 *  (*dev->paintr)(col, xmin, ymin, xmax, ymax)
 *  COLOR  col;
 *  int  xmin, ymin, xmax, ymax;
 *  {
 *	Paint a half-open rectangle from (xmin,ymin) to (xmax,ymax)
 *	with the color col.
 *  }
 *  (*dev->getcur)(xp, yp)
 *  int  *xp, *yp;
 *  {
 *	Get the cursor position entered by the user via mouse,
 *	joystick, etc.  Return the key hit by the user (usually
 *	MB1 or MB2).  Return ABORT to cancel.
 *	Can be NULL for devices without this capability.
 *  }
 *  (*dev->comout)(out)
 *  char  *out;
 *  {
 *	Print the string out on the device command line.  If the
 *	string ends with '\n', the message is considered complete,
 *	and the next call can erase it.
 *  }
 *  (*dev->comin)(in, prompt)
 *  char  *in, *prompt;
 *  {
 *	Print a prompt then read an edited input command line
 *	assuming the in buffer is big enough.  Unless prompt is NULL,
 *	the driver may substitute its own rview command.  This is
 *	the most reliable way to repaint areas of the screen.
 *	If the user enters an unrecognized control character,
 *	terminate input and return the string with only that character.
 *	The input string should not contain a newline.  The routines in
 *	editline.c may be useful.  Comin must work in consort with comout.
 *  }
 *  (*dev->flush)()
 *  {
 *	Flush output to the display.  This is guaranteed to be called
 *	frequently enough to keep the display up to date.
 *	This is an ideal time to check for device input.
 *	This function can be NULL for devices that don't need it.
 *  }
 *  xsiz, ysiz
 *	The maximum allowable x and y dimensions.  If any
 *	size is allowable, these should be set to MAXRES.
 *  inpready
 *	This variable should be made positive asynchronously
 *	when characters are ready on the input.  (Often easiest
 *	to check for input during calls to paintr.)
 */

#ifdef NOPROTO

extern void	editline();
extern void	tocombuf();
extern int	fromcombuf();
extern struct driver	*slave_init();
extern struct driver	*comm_init();
extern int	new_ctab();
extern int	get_pixel();
extern void	make_gmap();
extern void	set_cmap();
extern void	map_color();

#else
					/* defined in editline.c */
extern void	editline(char *buf, int (*c_get)(), void (*s_put)());
extern void	tocombuf(char *b, struct driver *d);
extern int	fromcombuf(char *b, struct driver *d);
					/* defined in devcomm.c */
extern struct driver	*slave_init(char *dname, char *id);
extern struct driver	*comm_init(char *dname, char *id);
					/* defined in colortab.c */
extern int	new_ctab(int ncolors);
extern int	get_pixel(COLOR col, void (*set_pixel)());
extern void	make_gmap(double gam);
extern void	set_cmap(BYTE *rmap, BYTE *gmap, BYTE *bmap);
extern void	map_color(BYTE rgb[3], COLOR col);

#endif
