/* RCSid $Id: ray.h,v 2.23 2004/11/05 17:36:55 greg Exp $ */
/*
 *  ray.h - header file for routines using rays.
 */
#ifndef _RAD_RAY_H_
#define _RAD_RAY_H_

#include  "standard.h"
#include  "octree.h"
#include  "object.h"
#include  "color.h"

#ifdef __cplusplus
extern "C" {
#endif

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
	void	(*revf)(struct ray *);	/* evaluation function for this ray */
	void	(*hitf)(OBJECT *, struct ray *);	/* custom hit test */
	OBJECT	robj;		/* intersected object number */
 	OBJREC	*ro;		/* intersected object (one with material) */
	double	rot;		/* distance to object */
	FVECT	rop;		/* intersection point */
	FVECT	ron;		/* intersection surface normal */
	double	rod;		/* -DOT(rdir, ron) */
	FULLXF	*rox;		/* object transformation */
	RREAL	uv[2];		/* local coordinates */
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

#ifndef AMBLLEN
#define AMBLLEN		512	/* max. ambient list length */
#endif
#define AMBWORD		12	/* average word length */

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

					/* defined in duphead.c */
extern void	headclean(void);
extern void	openheader(void);
extern void	dupheader(void);
					/* defined in persist.c */
extern void persistfile(char *pfn);
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
extern void	rayhit(OBJECT *oset, RAY *r);
extern void	raycont(RAY *r);
extern void	raytrans(RAY *r);
extern int	rayshade(RAY *r, int mod);
extern void	rayparticipate(RAY *r);
extern void	raytexture(RAY *r, OBJECT mod);
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
extern void init_drawsources(int rad);
					/* defined in rt/initotypes.c */
extern void initotypes(void);
					/* module main procedures */
extern void	rtrace(char *fname);
extern char * formstr(int  f);
extern void	rview(void);
extern void	rpict(int seq, char *pout, char *zout, char *prvr);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_RAY_H_ */

