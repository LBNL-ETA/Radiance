/* RCSid $Id: bsdf.h,v 2.1 2009/06/17 20:41:47 greg Exp $ */
/*
 * Header for BSDF i/o and access routines
 */
				/* up directions */
typedef enum {
	UDzneg=-3,
	UDyneg=-2,
	UDxneg=-1,
	UDunknown=0,
	UDxpos=1,
	UDypos=2,
	UDzpos=3
} UpDir;
				/* BSDF coordinate calculation routines */
				/* vectors always point away from surface */
typedef int	b_vecf(FVECT, int, void *);
typedef int	b_ndxf(FVECT, void *);
typedef double	b_radf(int, void *);

/* Bidirectional Scattering Distribution Function */
struct BSDF_data {
	int	ninc;			/* number of incoming directions */
	int	nout;			/* number of outgoing directions */
	void	*ib_priv;		/* input basis private data */
	b_vecf	*ib_vec;		/* get input vector from index */
	b_ndxf	*ib_ndx;		/* get input index from vector */
	b_radf	*ib_ohm;		/* get input radius for index */
	void	*ob_priv;		/* output basis private data */
	b_vecf	*ob_vec;		/* get output vector from index */
	b_ndxf	*ob_ndx;		/* get output index from vector */
	b_radf	*ob_ohm;		/* get output radius for index */
	float	*bsdf;			/* scattering distribution data */
};				/* bidirectional scattering distrib. func. */

#define getBSDF_incvec(v,b,i)	(*(b)->ib_vec)(v,i,(b)->ib_priv)
#define getBSDF_incndx(b,v)	(*(b)->ib_ndx)(v,(b)->ib_priv)
#define getBSDF_incohm(b,i)	(*(b)->ib_ohm)(i,(b)->ib_priv)
#define getBSDF_outvec(v,b,o)	(*(b)->ob_vec)(v,o,(b)->ob_priv)
#define getBSDF_outndx(b,v)	(*(b)->ob_ndx)(v,(b)->ob_priv)
#define getBSDF_outohm(b,o)	(*(b)->ob_ohm)(o,(b)->ob_priv)
#define BSDF_value(b,i,o)	(b)->bsdf[(o)*(b)->ninc + (i)]

extern struct BSDF_data *load_BSDF(char *fname);
extern void free_BSDF(struct BSDF_data *b);
extern int r_BSDF_incvec(FVECT v, struct BSDF_data *b, int i,
				double rv, MAT4 xm);
extern int r_BSDF_outvec(FVECT v, struct BSDF_data *b, int o,
				double rv, MAT4 xm);
extern int getBSDF_xfm(MAT4 xm, FVECT nrm, UpDir ud);
