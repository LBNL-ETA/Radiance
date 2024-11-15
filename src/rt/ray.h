/* RCSid $Id: ray.h,v 2.57 2024/11/15 20:47:42 greg Exp $ */
/*
 *  ray.h - header file for routines using rays.
 */
#ifndef _RAD_RAY_H_
#define _RAD_RAY_H_

#include  "standard.h"
#include  "octree.h"
#include  "object.h"
#include  "color.h"
#include  "pmapparm.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef RNUMBER
#define RNUMBER		size_t		/* ray counter (>= sizeof pointer) */
#endif

#define  MAXDIM		32	/* maximum number of sampling dimensions */

				/* ray type flags */
#define  PRIMARY	01		/* original ray */
#define  RSHADOW	02		/* reflected ray to light source */
#define  REFLECTED	04		/* reflected ray */
#define  REFRACTED	010		/* refracted (bent) ray */
#define  TRANS		020		/* transmitted/transferred ray */
#define  RAMBIENT	040		/* reflected diffuse interreflection */
#define  RSPECULAR	0100		/* reflected specular */
#define  TSHADOW	0200		/* transmitted shadow */
#define  TAMBIENT	0400		/* transmitted ambient */
#define  TSPECULAR	01000		/* transmitted specular */
#define  SHADOW		(RSHADOW|TSHADOW)
#define  AMBIENT	(RAMBIENT|TAMBIENT)
#define  SPECULAR	(RSPECULAR|TSPECULAR)

				/* reflected ray types */
#define  RAYREFL	(RSHADOW|REFLECTED|RAMBIENT|RSPECULAR)

/* Arrange so double's come first for optimal alignment */
/* Pointers and long's come second for 64-bit mode */
/* Int's next (unknown length), then floats, followed by short's & char's */
typedef struct ray {
	FVECT	rorg;		/* origin of ray */
	FVECT	rdir;		/* normalized direction of ray */
	RREAL	rmax;		/* maximum distance (aft clipping plane) */
	RREAL	rot;		/* distance to object */
	FVECT	rop;		/* intersection point */
	FVECT	ron;		/* intersection surface normal */
	RREAL	rod;		/* -DOT(rdir, ron) */
	RREAL	uv[2];		/* local coordinates */
	FVECT	pert;		/* surface normal perturbation */
	RREAL	rmt;		/* returned mirrored ray length */
	RREAL	rxt;		/* returned unmirrored ray length */
	const struct ray  *parent;	/* ray this originated from */
	OBJECT	*clipset;	/* set of objects currently clipped */
	OBJECT	*newcset;	/* next clipset, used for transmission */
	void	(*revf)(struct ray *);	/* ray evaluation function */
	void	(*hitf)(OBJECT *, struct ray *);	/* custom hit test */
 	OBJREC	*ro;		/* intersected object (one with material) */
	FULLXF	*rox;		/* object transformation */
	int	*slights;	/* list of lights to test for scattering */
	RNUMBER	rno;		/* unique ray number */
	OBJECT	robj;		/* intersected object number */
	int	rsrc;		/* source we're aiming for (or ones to skip) */
#ifdef SSKIPOPT
	float	scorr;		/* correction factor for included sources */
#endif
	float	rweight;	/* cumulative weight (for termination) */
	float	gecc;		/* scattering eccentricity coefficient */
	SCOLOR	rcoef;		/* contribution coefficient w.r.t. parent */
	SCOLOR	pcol;		/* pattern color */
	SCOLOR	mcol;		/* mirrored color contribution */
	SCOLOR	rcol;		/* returned radiance value */
	COLOR	cext;		/* medium extinction coefficient */
	COLOR	albedo;		/* medium scattering albedo */
	short	rflips;		/* surface orientation has been reversed */
	short	rlvl;		/* number of reflections for this ray */
	short	rtype;		/* ray type */
	short	crtype;		/* cumulative ray type */
}  RAY;

#define  rayvalue(r)	(*(r)->revf)(r)

#define  thrudir(r,v)	((r)->rod > 0 ^ DOT((r)->ron,v) > 0)

#define  raydistance(r)	(pbright((r)->mcol) > 0.5*pbright((r)->rcol) ? \
				(r)->rmt : (r)->rxt)

#define  rayreorient(r)	if ((r)->rflips & 1) flipsurface(r); else

extern char	VersionID[];	/* Radiance version ID string */
extern char	RFeatureList[];	/* newline-separated feature list */

extern CUBE	thescene;	/* our scene */
extern OBJECT	nsceneobjs;	/* number of objects in our scene */

extern RNUMBER	raynum;		/* next ray ID */
extern RNUMBER	nrays;		/* total rays traced so far */

extern OBJREC  Lamb;		/* a Lambertian surface */
extern OBJREC  Aftplane;	/* aft clipping object */

extern void	(*trace)(RAY*);	/* global trace reporting callback */

extern int	dimlist[];	/* dimension list for distribution */
extern int	ndims;		/* number of dimensions so far */
extern int	samplendx;	/* index for this sample */

extern int	do_irrad;	/* compute irradiance? */

extern int	rand_samp;	/* pure Monte Carlo sampling? */

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
	int	rand_samp;
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
	
	/* PMAP: photon mapping parameters */
	PhotonMapParams pmapParams [NUM_PMAP_TYPES];
} RAYPARAMS;

#define rpambmod(p,i)	( (i)>=AMBLLEN||(p)->amblndx[i]<0 ? \
			  (char *)NULL : (p)->amblval+(p)->amblndx[i] )

					/* defined in duphead.c */
extern void	headclean(void);
extern void	openheader(void);
extern void	dupheader(void);
					/* defined in persist.c */
extern void	persistfile(char *pfn);
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
extern char	*shm_boundary;
extern void	cow_memshare(void);
extern void	cow_doneshare(void);
					/* defined in raycalls.c */
extern void	ray_init(char *otnm);
extern void	ray_trace(RAY *r);
extern void	ray_done(int freall);
extern void	ray_save(RAYPARAMS *rp);
extern void	ray_restore(RAYPARAMS *rp);
extern void	ray_defaults(RAYPARAMS *rp);
					/* defined in raypcalls.c */
extern void	ray_pinit(char *otnm, int nproc);
extern int	ray_psend(RAY *r);
extern int	ray_pqueue(RAY *r);
extern int	ray_presult(RAY *r, int poll);
extern void	ray_pdone(int freall);
extern void	ray_popen(int nadd);
extern void	ray_pclose(int nsub);
					/* defined in ray_fifo.c */
extern int	(*ray_fifo_out)(RAY *r);
extern int	ray_fifo_in(RAY *r);
extern int	ray_fifo_flush(void);
					/* defined in raytrace.c */
extern int	rayorigin(RAY *r, int rt, const RAY *ro, const SCOLOR rc);
extern void	rayclear(RAY *r);
extern void	raytrace(RAY *r);
extern int	rayreject(OBJREC *o, RAY *r, double t, double rod);
extern void	rayhit(OBJECT *oset, RAY *r);
extern void	raycont(RAY *r);
extern void	raytrans(RAY *r);
extern int	raytirrad(OBJREC *m, RAY *r);
extern int	rayshade(RAY *r, int mod);
extern void	rayparticipate(RAY *r);
extern void	raytexture(RAY *r, OBJECT mod);
extern int	raymixture(RAY *r, OBJECT fore, OBJECT back, double coef);
extern void	raycontrib(SCOLOR rc, const RAY *r, int flags);
extern double	raydist(const RAY *r, int flags);
extern double	raynormal(FVECT norm, RAY *r);
extern void	newrayxf(RAY *r);
extern void	flipsurface(RAY *r);
extern int	localhit(RAY *r, CUBE *scene);
					/* defined in renderopts.c */
extern int	feature_status(int ac, char *av[]);
extern int	getrenderopt(int ac, char *av[]);
extern void	print_rdefaults(void);
					/* defined in srcdraw.c */
extern void	init_drawsources(int rad);
extern void	drawsources(COLORV *pic[], RGBPRIMP primp, float *zbf[],
			int x0, int xsiz, int y0, int ysiz);
					/* defined in rt/initotypes.c */
extern void	initotypes(void);
					/* module main procedures */
extern void	rtrace(char *fname, int nproc);
extern const char	*formstr(int  f);
extern void	rview(void);
extern void	rpict(int seq, char *pout, char *zout, char *prvr);

#ifdef __FAST_MATH__
#define	checknorm(vn)	(void)normalize(vn)
#else
#define checknorm(vn)
#endif

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RAY_H_ */

