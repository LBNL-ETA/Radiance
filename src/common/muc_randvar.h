/* RCSid $Id: muc_randvar.h,v 2.1 2016/07/14 17:32:12 greg Exp $ */
/*
**	Author: Christian Reetz (chr@riap.de)
*/
#ifndef __MUC_RANDVAR_H
#define	__MUC_RANDVAR_H

/**	\file muc_randvar.h
*	\author Ch. Reetz (chr@riap.de)
*	\brief representation of a random variables of arbitrary dimension
*/


#ifdef __cplusplus 
extern "C" {
#endif

#include "g3flist.h"

/** random variable structure*/
struct muc_rvar {
	int						n;			/*	number of samples */
	double					w;			/* sum of weights*/
	double*					sum;		/* sample values sum*/
	double*					sum_sqr;	/* squared sample values sum*/
	double*					min;		/* minimum component value*/
	double*					max;		/* maximum component value*/
	g3FList*				samples;	/* measured samples (if saved)*/
	int						save_samp;	/*	if true (default) save samples*/
};

/*	returns dimension of random variable*/
#define	muc_rvar_get_dim(rv)			(g3fl_get_comp_size(rv->samples))
/*	returns the number of samples*/
#define	muc_rvar_get_sample_size(rv)	(rv->n)

/* returns a reference to a sample list element*/
/* Take care to have save_samp flag set*/
#define	muc_rvar_get_sample(rv,id)		(g3fl_get(rv->samples, id))


int					muc_rvar_init(struct muc_rvar* rv);
struct muc_rvar*	muc_rvar_create();
void				muc_rvar_free(struct muc_rvar* rv);

/* set dimension of random variable (default 1). */
/*	random variable needs to be empty (\see muc_rvar_clear)*/
int		muc_rvar_set_dim(struct muc_rvar* rv,int dim);

/* set storing samples on (default) or off*/
/* random varianle must be cleared to succeed*/
int		muc_rvar_store_sample(struct muc_rvar* rv, int on);

/*	clear random variables samples. Don't change dimension.*/
void	muc_rvar_clear(struct muc_rvar* rv);

/*	Don't store measured sample values, just summed values.*/
/*	(Deletes samples if there are already stored some)*/
void	muc_rvar_no_samples(struct muc_rvar* rv);

/*	adds measured sample to sums (and stores them if \ref save_samp is true)*/
/*	The weight of the sample is assumed to be 1.0*/
int		muc_rvar_add_sample(struct muc_rvar* rv,double* sample);

/*	adds \param w weighted sample to sums */
/*	(and stores them if \ref save_samp is true)*/
int		muc_rvar_add_weighted_sample(struct muc_rvar* rv,double w,double* samp);

/*	returns current expected value in \param ex*/
int		muc_rvar_get_ex(struct muc_rvar* rv,double* ex);

/*	returns current variance in \param vx*/
int		muc_rvar_get_vx(struct muc_rvar* rv,double* vx);

/* returns current bounds*/
/*	\param bounds: array with len 2*component_size */
/* \return {min1, max1, min2, max2....}*/
int		muc_rvar_get_bounding_box(struct muc_rvar* rv, double* bounds);

/* return sample median (using averaging for even sample numbers)*/
/* \return false if samples are not stored*/
int		muc_rvar_get_median(struct muc_rvar* rv, double* median);

int		muc_rvar_get_percentile(struct muc_rvar* rv, double* median, double percentile);

#ifdef __cplusplus 
}
#endif
#endif
