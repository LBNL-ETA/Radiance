/* RCSid: $Id: view.h,v 2.7 2003/02/22 02:07:23 greg Exp $ */
/*
 *  view.h - header file for image generation.
 *
 *  Include after fvect.h
 *  Includes resolu.h
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

#include  "resolu.h"

				/* view types */
#define  VT_PER		'v'		/* perspective */
#define  VT_PAR		'l'		/* parallel */
#define  VT_ANG		'a'		/* angular fisheye */
#define  VT_HEM		'h'		/* hemispherical fisheye */
#define  VT_CYL		'c'		/* cylindrical panorama */

typedef struct {
	int  type;		/* view type */
	FVECT  vp;		/* view origin */
	FVECT  vdir;		/* view direction */
	FVECT  vup;		/* view up */
	double  horiz;		/* horizontal view size */
	double  vert;		/* vertical view size */
	double  hoff;		/* horizontal image offset */
	double  voff;		/* vertical image offset */
	double  vfore;		/* fore clipping plane */
	double  vaft;		/* aft clipping plane (<=0 for inf) */
	FVECT  hvec;		/* computed horizontal image vector */
	FVECT  vvec;		/* computed vertical image vector */
	double  hn2;		/* DOT(hvec,hvec) */
	double  vn2;		/* DOT(vvec,vvec) */
} VIEW;			/* view parameters */

extern VIEW  stdview;

#define  viewaspect(v)	sqrt((v)->vn2/(v)->hn2)

#define  STDVIEW	{VT_PER,{0.,0.,0.},{0.,1.,0.},{0.,0.,1.}, \
				45.,45.,0.,0.,0.,0., \
				{0.,0.,0.},{0.,0.,0.},0.,0.}

#define  VIEWSTR	"VIEW="
#define  VIEWSTRL	5

#ifdef NOPROTO

extern char	*setview();
extern void	normaspect();
extern double	viewray();
extern void	viewloc();
extern void	pix2loc();
extern void	loc2pix();
extern int	getviewopt();
extern int	sscanview();
extern void	fprintview();
extern char	*viewopt();
extern int	isview();
extern int	viewfile();

#else

extern char	*setview(VIEW *v);
extern void	normaspect(double va, double *ap, int *xp, int *yp);
extern double	viewray(FVECT orig, FVECT direc, VIEW *v, double x, double y);
extern void	viewloc(FVECT ip, VIEW *v, FVECT p);
extern void	pix2loc(FLOAT loc[2], RESOLU *rp, int px, int py);
extern void	loc2pix(int pp[2], RESOLU *rp, double lx, double ly);
extern int	getviewopt(VIEW *v, int ac, char *av[]);
extern int	sscanview(VIEW *vp, char *s);
extern void	fprintview(VIEW *vp, FILE *fp);
extern char	*viewopt(VIEW *vp);
extern int	isview(char *s);
extern int	viewfile(char *fname, VIEW *vp, RESOLU *rp);

#endif
