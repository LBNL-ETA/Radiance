/* RCSid: $Id: ray.h,v 2.12 2003/02/22 02:07:29 greg Exp $ */
/*
 *  ray.h - header file for routines using rays.
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

#include  "standard.h"

#include  "octree.h"

#include  "object.h"

#include  "color.h"

#define  MAXDIM		32	/* maximum number of dimensions */

				/* ray type flags */
#define  PRIMARY	01		/* original ray */
#define  SHADOW		02		/* ray to light source */
#define  REFLECTED	04		/* reflected ray */
#define  REFRACTED	010		/* refracted (bent) ray */
#define  TRANS		020		/* transmitted/transferred ray */
#define  AMBIENT	040		/* ray scattered for interreflection */
#define  SPECULAR	0100		/* ray scattered for specular */

				/* reflected ray types */
#define  RAYREFL	(SHADOW|REFLECTED|AMBIENT|SPECULAR)

typedef struct ray {
	unsigned long  rno;	/* unique ray number */
	int	rlvl;		/* number of reflections for this ray */
	float	rweight;	/* cumulative weight of this ray */
	short	rtype;		/* ray type */
	short	crtype;		/* cumulative ray type */
	struct ray  *parent;	/* ray this originated from */
	FVECT	rorg;		/* origin of ray */
	FVECT	rdir;		/* normalized direction of ray */
	double	rmax;		/* maximum distance (aft clipping plane) */
	int	rsrc;		/* source we're aiming for */
	OBJECT	*clipset;	/* set of objects currently clipped */
	OBJECT	*newcset;	/* next clipset, used for transmission */
	void	(*revf)();	/* evaluation function for this ray */
	OBJECT	robj;		/* intersected object number */
 	OBJREC	*ro;		/* intersected object (one with material) */
	double	rot;		/* distance to object */
	FVECT	rop;		/* intersection point */
	FVECT	ron;		/* intersection surface normal */
	double	rod;		/* -DOT(rdir, ron) */
	FULLXF	*rox;		/* object transformation */
	FVECT	pert;		/* surface normal perturbation */
	COLOR	pcol;		/* pattern color */
	COLOR	rcol;		/* returned ray value */
	double	rt;		/* returned effective ray length */
	COLOR	cext;		/* medium extinction coefficient */
	COLOR	albedo;		/* medium scattering albedo */
	float	gecc;		/* scattering eccentricity coefficient */
	int	*slights;	/* list of lights to test for scattering */
}  RAY;

#define  rayvalue(r)	(*(r)->revf)(r)

extern char  VersionID[];	/* Radiance version ID string */

extern CUBE	thescene;	/* our scene */
extern OBJECT	nsceneobjs;	/* number of objects in our scene */

extern unsigned long	raynum;	/* next ray ID */
extern unsigned long	nrays;	/* total rays traced so far */

extern OBJREC  Lamb;		/* a Lambertian surface */
extern OBJREC  Aftplane;	/* aft clipping object */

extern void	(*trace)();	/* global trace reporting callback */

extern int	dimlist[];	/* dimension list for distribution */
extern int	ndims;		/* number of dimensions so far */
extern int	samplendx;	/* index for this sample */

extern int	ray_savesiz;	/* size of parameter save buffer */

extern int	do_irrad;	/* compute irradiance? */

extern double	dstrsrc;	/* square source distribution */
extern double	shadthresh;	/* shadow threshold */
extern double	shadcert;	/* shadow testing certainty */
extern int	directrelay;	/* number of source relays */
extern int	vspretest;	/* virtual source pretest density */
extern int	directvis;	/* light sources visible to eye? */
extern double	srcsizerat;	/* maximum source size/dist. ratio */

extern double	specthresh;	/* specular sampling threshold */
extern double	specjitter;	/* specular sampling jitter */

extern COLOR	cextinction;	/* global extinction coefficient */
extern COLOR	salbedo;	/* global scattering albedo */
extern double	seccg;		/* global scattering eccentricity */
extern double	ssampdist;	/* scatter sampling distance */

extern int	backvis;	/* back face visibility */

extern int	maxdepth;	/* maximum recursion depth */
extern double	minweight;	/* minimum ray weight */

extern char	*ambfile;	/* ambient file name */
extern COLOR	ambval;		/* ambient value */
extern int	ambvwt;		/* initial weight for ambient value */
extern double	ambacc;		/* ambient accuracy */
extern int	ambres;		/* ambient resolution */
extern int	ambdiv;		/* ambient divisions */
extern int	ambssamp;	/* ambient super-samples */
extern int	ambounce;	/* ambient bounces */
extern char	*amblist[];	/* ambient include/exclude list */
extern int	ambincl;	/* include == 1, exclude == 0 */

extern int	ray_pnprocs;	/* number of child processes */
extern int	ray_pnidle;	/* number of idle processes */

#define AMBLLEN		128	/* max. ambient list length */
#define AMBWORD		8	/* average word length */

typedef struct {		/* rendering parameter holder */
	int	do_irrad;
	double	dstrsrc;
	double	shadthresh;
	double	shadcert;
	int	directrelay;
	int	vspretest;
	int	directvis;
	double	srcsizerat;
	COLOR	cextinction;
	COLOR	salbedo;
	double	seccg;
	double	ssampdist;
	double	specthresh;
	double	specjitter;
	int	backvis;
	int	maxdepth;
	double	minweight;
	char	ambfile[512];
	COLOR	ambval;
	int	ambvwt;
	double	ambacc;
	int	ambres;
	int	ambdiv;
	int	ambssamp;
	int	ambounce;
	int	ambincl;
	short	amblndx[AMBLLEN+1];
	char	amblval[AMBLLEN*AMBWORD];
} RAYPARAMS;

#define rpambmod(p,i)	( (i)>=AMBLLEN||(p)->amblndx[i]<0 ? \
			  (char *)NULL : (p)->amblval+(p)->amblndx[i] )

#ifdef NOPROTO

extern void	headclean();
extern void	openheader();
extern void	dupheader();
extern void	pfdetach();
extern void	pfclean();
extern void	pflock();
extern void	pfhold();
extern void	io_process();
extern int	free_objs();
extern void	free_objmem();
extern int	load_os();
extern void	preload_objs();
extern void	ray_init();
extern void	ray_trace();
extern void	ray_done();
extern void	ray_save();
extern void	ray_restore();
extern void	ray_defaults();
extern void	ray_pinit();
extern void	ray_psend();
extern int	ray_pqueue();
extern int	ray_presult();
extern void	ray_pdone();
extern void	ray_popen();
extern void	ray_pclose();
extern int	rayorigin();
extern void	rayclear();
extern void	raytrace();
extern void	raycont();
extern void	raytrans();
extern int	rayshade();
extern void	rayparticipate();
extern int	raymixture();
extern double	raydist();
extern double	raynormal();
extern void	newrayxf();
extern void	flipsurface();
extern int	localhit();
extern int	getrenderopt();
extern void	print_rdefaults();
extern void	drawsources();
extern void	rtrace();
extern void	rview();
extern void	rpict();

#else
					/* defined in duphead.c */
extern void	headclean(void);
extern void	openheader(void);
extern void	dupheader(void);
extern void	pfdetach(void);
extern void	pfclean(void);
extern void	pflock(int lf);
extern void	pfhold(void);
extern void	io_process(void);
					/* defined in freeobjmem.c */
extern int	free_objs(OBJECT on, OBJECT no);
extern void	free_objmem(void);
					/* defined in preload.c */
extern int	load_os(OBJREC *op);
extern void	preload_objs(void);
					/* defined in raycalls.c */
extern void	ray_init(char *otnm);
extern void	ray_trace(RAY *r);
extern void	ray_done(int freall);
extern void	ray_save(RAYPARAMS *rp);
extern void	ray_restore(RAYPARAMS *rp);
extern void	ray_defaults(RAYPARAMS *rp);
					/* defined in raypcalls.c */
extern void	ray_pinit(char *otnm, int nproc);
extern void	ray_psend(RAY *r);
extern int	ray_pqueue(RAY *r);
extern int	ray_presult(RAY *r, int poll);
extern void	ray_pdone(int freall);
extern void	ray_popen(int nadd);
extern void	ray_pclose(int nsub);
					/* defined in raytrace.c */
extern int	rayorigin(RAY *r, RAY *ro, int rt, double rw);
extern void	rayclear(RAY *r);
extern void	raytrace(RAY *r);
extern void	raycont(RAY *r);
extern void	raytrans(RAY *r);
extern int	rayshade(RAY *r, int mod);
extern void	rayparticipate(RAY *r);
extern int	raymixture(RAY *r, OBJECT fore, OBJECT back, double coef);
extern double	raydist(RAY *r, int flags);
extern double	raynormal(FVECT norm, RAY *r);
extern void	newrayxf(RAY *r);
extern void	flipsurface(RAY *r);
extern int	localhit(RAY *r, CUBE *scene);
					/* defined in renderopts.c */
extern int	getrenderopt(int ac, char *av[]);
extern void	print_rdefaults(void);
					/* defined in srcdraw.c */
extern void	drawsources(COLOR *pic[], float *zbf[],
			int x0, int xsiz, int y0, int ysiz);
					/* module main procedures */
extern void	rtrace(char *fname);
extern void	rview(void);
extern void	rpict(int seq, char *pout, char *zout, char *prvr);

#endif
