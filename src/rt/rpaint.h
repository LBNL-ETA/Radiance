/* RCSid: $Id: rpaint.h,v 2.2 2003/02/22 02:07:29 greg Exp $ */
/*
 *  rpaint.h - header file for image painting.
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
