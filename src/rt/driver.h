/* RCSid $Id: driver.h,v 2.10 2011/05/20 02:06:39 greg Exp $ */
/*
 *  driver.h - header file for interactive device drivers.
 */
#ifndef _RAD_DRIVER_H_
#define _RAD_DRIVER_H_

#include "color.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct driver *dr_initf_t(char *dname, char *id);
typedef int dr_getchf_t(void);
typedef void dr_newcolrf_t(int ndx, int r, int g, int b);

typedef void dr_closef_t(void);
typedef void dr_clearf_t(int, int);
typedef void dr_paintrf_t(COLOR col, int xmin, int ymin, int xmax, int ymax);
typedef int  dr_getcurf_t(int*,int*);
typedef void dr_comoutf_t(char*);
typedef void dr_cominf_t(char*,char*);
typedef void dr_flushf_t(void);

struct driver {				/* driver functions */
	dr_closef_t *close;			/* close device */
	dr_clearf_t *clear;			/* clear device */
	dr_paintrf_t *paintr;			/* paint rectangle */
	dr_getcurf_t *getcur;			/* get cursor position */
	dr_comoutf_t *comout;			/* command line output */
	dr_cominf_t *comin;			/* command line input */
	dr_flushf_t *flush;			/* flush output */
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
	dr_initf_t *init;		/* initialize device */
}  devtable[];			/* supported devices */

extern char  dev_default[];		/* default device name */

#define  MB1		('\n')		/* mouse button 1 */
#define  MB2		('\r')		/* mouse button 2 */
#define  MB3		(' ')		/* mouse button 3 */
#define  ABORT		('C'-'@')	/* abort key */

/*
 *  How to write an interactive display driver for rview.
 *  ----------------------------------------------------
 *
 *  static struct driver dname_driver;
 *
 *  extern dr_initf_t dname_init;
 *
 *  extern struct driver *
 *  dname_init(
 *    char  *name,
 *    char  *id
 *  )
 *  {
 *	Initialize device and return pointer to driver
 *	dname_driver. Return NULL if an error occurred.
 *	The name string identifies the driver,
 *	and the id string identifies the client.
 *	A device can be open by at most one client.
 *	Be verbose in error reports; call eputs().
 *	If device has its own error output, set erract.
 *  This function then needs to be inserted into the
 *  device table in devtable.c.
 *  }
 *
 *
 *  static dr_closef_t dname_close;
 *  
 *  dname_driver.close = dname_close;
 *
 *  static void
 *  dname_close(void)
 *  {
 *	Close the device.  Reset error vectors.
 *  }
 *
 *
 *  static dr_clearf_t dname_clear;
 *  
 *  dname_driver.clear = dname_clear;
 *
 *  static void
 *  dname_clear(
 *    int  xres,
 *    int  yres;
 *  )
 *  {
 *	Clear the device for xres by yres output.  This call will
 *	be made prior to any output.
 *  }
 *
 *
 *  static dr_paintrf_t dname_paintr;
 *  
 *  dname_driver.paintr = dname_paintr;
 *
 *  static void
 *  dname_paintr(
 *    COLOR  col,
 *    int  xmin,
 *    int  ymin,
 *    int  xmax,
 *    int  ymax
 *  )
 *  {
 *	Paint a half-open rectangle from (xmin,ymin) to (xmax,ymax)
 *	with the color col.
 *  }
 *
 *
 *  static dr_getcurf_t dname_getcur;
 *  
 *  dname_driver.getcur = dname_getcur;
 *
 *  static int
 *  dname_getcur(
 *    int  *xp,
 *    int  *yp
 *  )
 *  {
 *	Get the cursor position entered by the user via mouse,
 *	joystick, etc.  Return the key hit by the user (usually
 *	MB1 or MB2).  Return ABORT to cancel.
 *	Can be NULL for devices without this capability.
 *  }
 *
 *
 *  static dr_comoutf_t dname_comout;
 *  
 *  dname_driver.comout = dname_comout;
 *
 *  static void
 *  dname_comout(
 *    char  *out
 *  )
 *  {
 *	Print the string out on the device command line.  If the
 *	string ends with '\n', the message is considered complete,
 *	and the next call can erase it.
 *  }
 *
 *
 *  static dr_cominf_t dname_comin;
 *  
 *  dname_driver.comin = dname_comin;
 *
 *  static void
 *  dname_comin(
 *    char  *in,
 *    char  *prompt
 *  )
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
 *
 *
 *  static dr_flushf_t dname_flush;
 *  
 *  dname_driver.flush = dname_flush;
 *
 *  static void
 *  dname_flush(void)
 *  {
 *	Flush output to the display.  This is guaranteed to be called
 *	frequently enough to keep the display up to date.
 *	This is an ideal time to check for device input.
 *	This function can be NULL for devices that don't need it.
 *  }
 *
 *
 *  dname_driver.xsiz
 *  dname_driver.ysiz
 *
 *	The maximum allowable x and y dimensions.  If any
 *	size is allowable, these should be set to MAXRES.
 *
 *
 *  dname_driver.inpready
 *
 *	This variable should be made positive asynchronously
 *	when characters are ready on the input.  (Often easiest
 *	to check for input during calls to paintr.)
 */

					/* defined in editline.c */
extern void	editline(char *buf, dr_getchf_t *c_get, dr_comoutf_t *s_put);
extern void	tocombuf(char *b, struct driver *d);
extern int	fromcombuf(char *b, struct driver *d);

					/* defined in devcomm.c */
extern dr_initf_t slave_init; /* XXX should probably be in a seperate file */
extern struct driver	*comm_init(char *dname, char *id);

					/* defined in colortab.c */
extern int	new_ctab(int ncolors);
extern int	get_pixel(COLOR col, dr_newcolrf_t *newcolr);
extern void	make_gmap(double gam);
extern void	set_cmap(uby8 *rmap, uby8 *gmap, uby8 *bmap);
extern void	map_color(uby8 rgb[3], COLOR col);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_DRIVER_H_ */

