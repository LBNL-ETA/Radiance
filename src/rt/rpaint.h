/* RCSid $Id: rpaint.h,v 2.3 2003/02/25 02:47:23 greg Exp $ */
/*
 *  rpaint.h - header file for image painting.
 */

#include "copyright.h"

#include  "driver.h"

#include  "view.h"

typedef short  COORD;		/* an image coordinate */

typedef struct pnode {
	struct pnode  *kid;		/* children */
	COORD  x, y;			/* position */
	COLOR  v;			/* value */
}  PNODE;			/* a paint node */

				/* child ordering */
#define  DL		0		/* down left */
#define  DR		1		/* down right */
#define  UL		2		/* up left */
#define  UR		3		/* up right */

#define  newptree()	(PNODE *)calloc(4, sizeof(PNODE))

typedef struct {
	COORD  l, d, r, u;		/* left, down, right, up */
}  RECT;			/* a rectangle */

extern PNODE  ptrunk;		/* the base of the image tree */

extern VIEW  ourview;		/* current view parameters */
extern VIEW  oldview;		/* previous view parameters */
extern int  hresolu, vresolu;	/* image resolution */

extern int  greyscale;		/* map colors to brightness? */

extern int  pdepth;		/* image depth in current frame */
extern RECT  pframe;		/* current frame rectangle */

extern double  exposure;	/* exposure for scene */

extern struct driver  *dev;	/* driver functions */

#ifdef NOPROTO

extern void	devopen();
extern void	devclose();
extern void	printdevices();
extern void	fillreserves();
extern void	freereserves();
extern void	command();
extern void	rsample();
extern int	refine();
extern void	getframe();
extern void	getrepaint();
extern void	getview();
extern void	lastview();
extern void	saveview();
extern void	loadview();
extern void	getaim();
extern void	getmove();
extern void	getrotate();
extern void	getpivot();
extern void	getexposure();
extern int	getparam();
extern void	setparam();
extern void	traceray();
extern void	writepict();
extern int	getrect();
extern int	getinterest();
extern float	*greyof();
extern void	paint();
extern void	newimage();
extern void	redraw();
extern void	repaint();
extern void	paintrec();
extern PNODE	*findrect();
extern void	scalepict();
extern void	getpictcolrs();
extern void	freepkids();
extern void	newview();
extern void	moveview();
extern void	pcopy();
extern void	zoomview();

#else
				/* defined in rview.c */
extern void	devopen(char *dname);
extern void	devclose(void);
extern void	printdevices(void);
extern void	fillreserves(void);
extern void	freereserves(void);
extern void	command(char *prompt);
extern void	rsample(void);
extern int	refine(PNODE *p, int xmin, int ymin,
				int xmax, int ymax, int pd);
				/* defined in rv2.c */
extern void	getframe(char *s);
extern void	getrepaint(char *s);
extern void	getview(char *s);
extern void	lastview(char *s);
extern void	saveview(char *s);
extern void	loadview(char *s);
extern void	getaim(char *s);
extern void	getmove(char *s);
extern void	getrotate(char *s);
extern void	getpivot(char *s);
extern void	getexposure(char *s);
extern int	getparam(char *str, char *dsc, int typ, void *p);
extern void	setparam(char *s);
extern void	traceray(char *s);
extern void	writepict(char *s);
				/* defined in rv3.c */
extern int	getrect(char *s, RECT *r);
extern int	getinterest(char *s, int direc, FVECT vec, double *mp);
extern float	*greyof(COLOR col);
extern void	paint(PNODE *p, int xmin, int ymin, int xmax, int ymax);
extern void	newimage(void);
extern void	redraw(void);
extern void	repaint(int xmin, int ymin, int xmax, int ymax);
extern void	paintrect(PNODE *p, int xmin, int ymin,
				int xmax, int ymax, RECT *r);
extern PNODE	*findrect(int x, int y, PNODE *p, RECT *r, int pd);
extern void	scalepict(PNODE *p, double sf);
extern void	getpictcolrs(int yoff, COLR *scan, PNODE *p,
			int xsiz, int ysiz);
extern void	freepkids(PNODE *p);
extern void	newview(VIEW *vp);
extern void	moveview(double angle, double elev, double mag, FVECT vc);
extern void	pcopy(PNODE *p1, PNODE *p2);
extern void	zoomview(VIEW *vp, double zf);

#endif
