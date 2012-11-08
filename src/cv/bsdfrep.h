/* RCSid $Id: bsdfrep.h,v 2.6 2012/11/08 00:31:17 greg Exp $ */
/*
 * Definitions for BSDF representation used to interpolate measured data.
 *
 *	G. Ward
 */

#include "bsdf.h"

#define DEBUG		1

#ifndef GRIDRES
#define GRIDRES		200		/* grid resolution per side */
#endif
					/* convert to/from coded radians */
#define ANG2R(r)	(int)((r)*((1<<16)/M_PI))
#define R2ANG(c)	(((c)+.5)*(M_PI/(1<<16)))

typedef struct {
	float		vsum;		/* DSF sum */
	unsigned short	nval;		/* number of values in sum */
	unsigned short	crad;		/* radius (coded angle) */
} GRIDVAL;			/* grid value */

typedef struct {
	float		peak;		/* lobe value at peak */
	unsigned short	crad;		/* radius (coded angle) */
	unsigned char	gx, gy;		/* grid position */
} RBFVAL;			/* radial basis function value */

struct s_rbfnode;		/* forward declaration of RBF struct */

typedef struct s_migration {
	struct s_migration	*next;		/* next in global edge list */
	struct s_rbfnode	*rbfv[2];	/* from,to vertex */
	struct s_migration	*enxt[2];	/* next from,to sibling */
	float			mtx[1];		/* matrix (extends struct) */
} MIGRATION;			/* migration link (winged edge structure) */

typedef struct s_rbfnode {
	int			ord;		/* ordinal position in list */
	struct s_rbfnode	*next;		/* next in global RBF list */
	MIGRATION		*ejl;		/* edge list for this vertex */
	FVECT			invec;		/* incident vector direction */
	double			vtotal;		/* volume for normalization */
	int			nrbf;		/* number of RBFs */
	RBFVAL			rbfa[1];	/* RBF array (extends struct) */
} RBFNODE;			/* RBF representation of DSF @ 1 incidence */

				/* symmetry operations */
#define MIRROR_X	1		/* mirror(ed) x-coordinate */
#define MIRROR_Y	2		/* mirror(ed) y-coordinate */

				/* represented incident quadrants */
#define INP_QUAD1	1		/* 0-90 degree quadrant */
#define	INP_QUAD2	2		/* 90-180 degree quadrant */
#define INP_QUAD3	4		/* 180-270 degree quadrant */
#define INP_QUAD4	8		/* 270-360 degree quadrant */

				/* active grid resolution */
extern int		grid_res;
				/* coverage/symmetry using INP_QUAD? flags */
extern int		inp_coverage;

				/* all incident angles in-plane so far? */
extern int		single_plane_incident;

				/* input/output orientations */
extern int		input_orient;
extern int		output_orient;

				/* processed incident DSF measurements */
extern RBFNODE		*dsf_list;

				/* RBF-linking matrices (edges) */
extern MIGRATION	*mig_list;

#define mtx_nrows(m)	(m)->rbfv[0]->nrbf
#define mtx_ncols(m)	(m)->rbfv[1]->nrbf
#define mtx_coef(m,i,j)	(m)->mtx[(i)*mtx_ncols(m) + (j)]
#define is_src(rbf,m)	((rbf) == (m)->rbfv[0])
#define is_dest(rbf,m)	((rbf) == (m)->rbfv[1])
#define nextedge(rbf,m)	(m)->enxt[is_dest(rbf,m)]
#define opp_rbf(rbf,m)	(m)->rbfv[is_src(rbf,m)]

#define	round(v)	(int)((v) + .5 - ((v) < -.5))

#define	BSDFREP_FMT	"BSDF_RBFmesh"

				/* global argv[0] */
extern char		*progname;

				/* get theta value in degrees [0,180) range */
#define get_theta180(v)	((180./M_PI)*acos((v)[2]))
				/* get phi value in degrees, [0,360) range */
#define	get_phi360(v)	((180./M_PI)*atan2((v)[1],(v)[0]) + 360.*((v)[1]<0))

				/* our loaded grid for this incident angle */
extern double		theta_in_deg, phi_in_deg;
extern GRIDVAL		dsf_grid[GRIDRES][GRIDRES];

/* Register new input direction */
extern int		new_input_direction(double new_theta, double new_phi);

#define	new_input_vector(v)\
			new_input_direction(get_theta180(v),get_phi360(v))

/* Apply symmetry to the given vector based on distribution */
extern int		use_symmetry(FVECT vec);

/* Reverse symmetry based on what was done before */
extern void		rev_symmetry(FVECT vec, int sym);

/* Reverse symmetry for an RBF distribution */
extern void		rev_rbf_symmetry(RBFNODE *rbf, int sym);

/* Rotate RBF to correspond to given incident vector */
extern void		rotate_rbf(RBFNODE *rbf, const FVECT invec);

/* Compute volume associated with Gaussian lobe */
extern double		rbf_volume(const RBFVAL *rbfp);

/* Compute outgoing vector from grid position */
extern void		ovec_from_pos(FVECT vec, int xpos, int ypos);

/* Compute grid position from normalized input/output vector */
extern void		pos_from_vec(int pos[2], const FVECT vec);

/* Evaluate RBF for DSF at the given normalized outgoing direction */
extern double		eval_rbfrep(const RBFNODE *rp, const FVECT outvec);

/* Insert a new directional scattering function in our global list */
extern int		insert_dsf(RBFNODE *newrbf);

/* Get the DSF indicated by its ordinal position */
extern RBFNODE *	get_dsf(int ord);

/* Get triangle surface orientation (unnormalized) */
extern void		tri_orient(FVECT vres, const FVECT v1,
					const FVECT v2, const FVECT v3);

/* Determine if vertex order is reversed (inward normal) */
extern int		is_rev_tri(const FVECT v1,
					const FVECT v2, const FVECT v3);

/* Find vertices completing triangles on either side of the given edge */
extern int		get_triangles(RBFNODE *rbfv[2], const MIGRATION *mig);

/* Clear our BSDF representation and free memory */
extern void		clear_bsdf_rep(void);

/* Write our BSDF mesh interpolant out to the given binary stream */
extern void		save_bsdf_rep(FILE *ofp);

/* Read a BSDF mesh interpolant from the given binary stream */
extern int		load_bsdf_rep(FILE *ifp);

/* Start new DSF input grid */
extern void		new_bsdf_data(double new_theta, double new_phi);

/* Add BSDF data point */
extern void		add_bsdf_data(double theta_out, double phi_out,
					double val, int isDSF);

/* Count up filled nodes and build RBF representation from current grid */ 
extern RBFNODE *	make_rbfrep(void);

/* Build our triangle mesh from recorded RBFs */
extern void		build_mesh(void);

/* Find edge(s) for interpolating the given vector, applying symmetry */
extern int		get_interp(MIGRATION *miga[3], FVECT invec);

/* Partially advect between recorded incident angles and allocate new RBF */
extern RBFNODE *	advect_rbf(const FVECT invec);
