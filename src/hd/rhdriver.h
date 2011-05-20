/* RCSid: $Id: rhdriver.h,v 3.26 2011/05/20 02:06:39 greg Exp $ */
/*
 * Header file for holodeck device driver routines.
 */
#ifndef _RAD_RHDRIVER_H_
#define _RAD_RHDRIVER_H_

#include "view.h"

#ifdef __cplusplus
extern "C" {
#endif

extern struct driver {
	char	*name;		/* holodeck name or title */
	VIEW	v;		/* base view parameters */
	int	hres, vres;	/* base view resolution */
	int	ifd;		/* input file descriptor (for select) */
	int	firstuse;	/* non-zero if driver can't recycle samples */
	int	inpready;	/* number of unprocessed input events */
} odev;			/* our open device */

extern char	odev_args[];	/* command arguments, if any */

extern int	imm_mode;	/* bundles are being delivered immediately */

extern double	eyesepdist;	/* world eye separation distance */

				/* user commands */
#define	DC_SETVIEW	0		/* set the base view */
#define	DC_GETVIEW	1		/* print the current base view */
#define	DC_LASTVIEW	2		/* restore previous view */
#define DC_FOCUS	3		/* view focus */
#define	DC_PAUSE	4		/* pause the current calculation */
#define	DC_RESUME	5		/* resume the calculation */
#define	DC_REDRAW	6		/* redraw from server */
#define	DC_KILL		7		/* kill rtrace process(es) */
#define	DC_RESTART	8		/* restart rtrace process(es) */
#define DC_CLOBBER	9		/* clobber holodeck file */
#define	DC_QUIT		10		/* quit the program */

#define	DC_NCMDS	11		/* number of commands */

				/* dev_input() returns flags from above */
#define DFL(dc)		(1<<(dc))

#define	CTRL(c)		((c)-'@')
				/* commands entered in display window */
#define DV_INIT		{'\0','v','l','f','p','\r',CTRL('L'),'K','R','C','q'}
				/* commands entered on stdin */
#define	DC_INIT		{"VIEW=","where","last","frame","pause","resume",\
				"redraw","kill","restart","clobber","quit"}


/************************************************************************
 * Driver routines (all are required):


void
dev_open(nm)		: prepare the device
char	*nm;		: appropriate title bar annotation

Sets fields of odev structure and prepares the display for i/o.
The base view type, horizontal and vertical view angles and other
default parameters in odev.v should also be assigned.


int
dev_view(nv)		: set base view parameters
VIEW	*nv;		: the new view

Updates the display for the given base view change.
Look for nv==&odev.v when making view current after dev_input()
returns DEV_NEWVIEW flag.  Return 1 on success, or 0 if the
new view would conflict with device requirements.  In the latter
case, reset parameters in nv to make it more agreeable, calling
error(COMMAND, "appropriate user warning").


void
dev_clear()		: clear device memory

Clear the device memory in preparation for fresh data.  Clearing
the screen is optional.


void
dev_value(c, d, p)	: register new point of light
COLR	c;		: pixel color (RGBE)
FVECT	d;		: ray direction vector
FVECT	p;		: world intersection point

Add the given color point to the display output queue.  If imm_mode is
non-zero, then values are being sent in rapid succession.  If p is NULL,
then the point is at infinity.


int
dev_flush()		: flush the output and prepare for select call

Updates display, taking any pending action required before select(2) call.
Returns non-zero if there is device input available, setting odev.inpready.


int
dev_input()		: process pending display input

Called when odev struct file descriptor shows input is ready.
Returns flags indicating actions to take in the control process.
If the DC_VIEW or DC_RESIZE flag is returned, the odev
structure must be updated beforehand.  If the DC_FOCUS
flag is set, then odev_args contains the frame dimensions.
No events shall remain when this function returns,
and odev.inpready will therefore be 0.

void
dev_auxcom(cmd, args)	: process auxiliary command
char	*cmd, *args;	: command name and argument string

Execute an auxiliary command (not one of those listed at the head of
this file).  The cmd argument points to the command name itself, and
the args argument points to a string with the rest of the input line.
If the command isn't known or there ARE no auxiliary commands, print
an appropriate COMMAND error message and return.


VIEW *
dev_auxview(n, hv)	: return nth auxiliary view
int	n;		: auxiliary view number
int	hv[2];		: returned horiz. and vert. image resolution

Return the nth auxiliary view associated with the current base view.
The hv entries are assigned the horizontal and vertical view resolution,
respectively.  Function returns NULL if there are no more auxiliary
views.  The zeroeth auxiliary view is the base view itself.


void
dev_section(gf, pf)	: add geometry and ports for rendering
char	*gf;		: geometry file name
char	*pf;		: portal file name(s)

Add the given geometry file to the list of geometry to render for
intermediate views if direct geometry rendering is available.  The
second argument gives the name(s) of any portal geometry files
associated with this section.  Multiple portal file names are separated
by spaces.  A single octree file will be given for the geometry, ending
in the ".oct" suffix.  Portal files will be given as zero or more
Radiance scene description file names.  If no portals are given for
this section, the string may be NULL.  The character strings are
guaranteed to be static (or permanently allocated) such that they may
be safely stored as a pointer.  The same pointers or file lists may be
(and often are) given repeatedly.  If a given geometry file does not
exist, the call should be silently ignored.  If gf is NULL, then the
last section has been given, and the display can be updated with the
new information.


void
dev_close()		: close the display

Close the display device and free up resources in preparation for exit.
Set odev.v.type=0 and odev.hres=odev.vres=0 when done.


 ************************************************************************/


/*
extern VIEW	*dev_auxview();
*/
extern int16	*beam_view(VIEW *vn, int hr, int vr);


extern void dev_open(char  *id);
extern void dev_close(void);
extern void dev_clear(void);
extern int dev_view(register VIEW *nv);
extern void dev_section(char *gfn, char *pfn);
extern void dev_auxcom(char *cmd, char *args);
extern VIEW *dev_auxview(int n, int hvres[2]);
extern int dev_input(void);
extern void dev_value(COLR c, FVECT d, FVECT p);
extern int dev_flush(void);
extern void dev_paintr(uby8 rgb[3], int xmin, int ymin, int xmax, int ymax);
extern void dev_cone(uby8 rgb[3], FVECT ip, double rad);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RHDRIVER_H_ */

