#ifndef lint
static const char RCSid[] = "$Id: muc_randvar.c,v 2.1 2016/07/14 17:32:12 greg Exp $";
#endif
#include <stdio.h>
#include <stdlib.h>

#include "fvect.h"
#include "muc_randvar.h"

int		muc_rvar_init(struct muc_rvar* rv)
{
	rv->n = 0;
	rv->w = 0.0;
	if (!(rv->samples = g3fl_create(1)))
		return 0;
	rv->sum = NULL;
	rv->sum_sqr = NULL;
	rv->min = NULL;
	rv->max = NULL;
	rv->save_samp = 1;
	return muc_rvar_set_dim(rv,1);
}

struct  muc_rvar*	muc_rvar_create()
{
	struct muc_rvar* res;

	if (!(res = (struct  muc_rvar*) malloc(sizeof(struct  muc_rvar))))
		return NULL;
	if (!(muc_rvar_init(res))) {
		free(res);
		return NULL;
	}
	return res;
}

void	muc_rvar_free(struct muc_rvar* rv)
{
	free(rv->sum);
	free(rv->sum_sqr);
	free(rv->min);
	free(rv->max);
	g3fl_free(rv->samples);
	free(rv);
}

int		muc_rvar_set_dim(struct muc_rvar* rv,int dim)
{
	if (rv->n != 0) {
		fprintf(stderr,"muc_rvar_set_dim: variable not cleared\n");
		return 0;
	}
	if (!(rv->sum = realloc(rv->sum,dim*sizeof(double))))
		return 0;
	if (!(rv->sum_sqr = realloc(rv->sum_sqr,dim*sizeof(double))))
		return 0;
	if (!(rv->min = realloc(rv->min,dim*sizeof(double))))
		return 0;
	if (!(rv->max = realloc(rv->max,dim*sizeof(double))))
		return 0;
	rv->samples->comp_size = dim;
	muc_rvar_clear(rv);
	return 1;
}

void	muc_rvar_clear(struct muc_rvar* rv)
{
	int i;

	rv->samples->size=0;
	rv->w = 0.0;
	rv->n = 0;
	for(i=0;i<muc_rvar_get_dim(rv);i++) {
		rv->sum[i] = 0.0;
		rv->sum_sqr[i] = 0.0;
		rv->min[i] = FHUGE;
		rv->max[i] = -FHUGE;
	}
}

void	muc_rvar_no_samples(struct muc_rvar* rv)
{
	rv->samples->size = 0;
	rv->save_samp = 0;
}

int		muc_rvar_store_sample(struct muc_rvar* rv, int on)
{
	if (rv->n != 0) {
		fprintf(stderr,"muc_rvar_store_sample: clear samples first\n");
		return 0;
	}
	rv->save_samp = on;
	return 1;
}

int		muc_rvar_add_sample(struct muc_rvar* rv,double* s)
{
	return muc_rvar_add_weighted_sample(rv,1.0,s);
}

int	muc_rvar_add_weighted_sample(struct muc_rvar* rv,double w,double* s)
{
	int i;
	double val;

	if (rv->save_samp)
		if (!(g3fl_append(rv->samples,s)))
			return 0;
	for(i=0;i<muc_rvar_get_dim(rv);i++) {
		val = s[i]*w;
		rv->sum[i] += val;
		rv->sum_sqr[i] += s[i]*s[i]*w;
		if (rv->min[i] > val)
			rv->min[i] = val;
		if (rv->max[i] < val)
			rv->max[i] = val;
	}
	rv->w += w;
	rv->n++;
	return 1;
}
	
int		muc_rvar_get_vx(struct muc_rvar* rv,double* vx)
{
	int i;
	double ex,ex2;

	if (rv->w == 0.0) {
		return 0;
	}
	for(i=0;i<muc_rvar_get_dim(rv);i++) {
		ex = rv->sum[i]/rv->w;
		ex2 = rv->sum_sqr[i]/rv->w;
		vx[i] = ex2 - ex*ex;
	}
	return 1;
}

int		muc_rvar_get_ex(struct muc_rvar* rv,double* ex)
{
	int i;

	if (rv->w == 0.0) {
		return 0;
	}
	for(i=0;i<muc_rvar_get_dim(rv);i++) {
		ex[i] = rv->sum[i]/rv->w;
	}
	return 1;
}
int		muc_rvar_get_sum(struct muc_rvar* rv,double* sum2)
{
	int i;

	for(i=0;i<muc_rvar_get_dim(rv);i++) {
		sum2[i] = rv->sum[i];
	}
	return 1;
}

int		muc_rvar_get_median(struct muc_rvar* rv, double* median)
{
	double val;
	int i;

	if (!(rv->save_samp)) {
		fprintf(stderr,"muc_rvar_get_median: samples are not stored\n");
		return 0;
	}
	for(i=0;i<muc_rvar_get_dim(rv);i++) {
		g3fl_sort(rv->samples, i);
		val = g3fl_get(rv->samples, rv->n/2)[i];
		
		if (rv->n % 2 == 0) {
			val += g3fl_get(rv->samples,rv->n/2 - 1)[i];
			val /= 2.0;
		}
		median[i] = val;
	}
	return 1;
}
int		muc_rvar_get_percentile(struct muc_rvar* rv, double* median, double percentile)
{
	double val;
	int i;

	if (!(rv->save_samp)) {
		fprintf(stderr,"muc_rvar_get_median: samples are not stored\n");
		return 0;
	}
	for(i=0;i<muc_rvar_get_dim(rv);i++) {
		g3fl_sort(rv->samples, i);
		val = g3fl_get(rv->samples, rv->n*percentile)[i];
		
		if (rv->n % 1/percentile == 0) {
			val += g3fl_get(rv->samples,rv->n*percentile - 1)[i];
			val /= 2.0;
		}
		median[i] = val;
	}
	return 1;
}

int		muc_rvar_get_bounding_box(struct muc_rvar* rv, double* bounds)
{
	int i;

	if (rv->n == 0) 
		return 0;
	for(i=0;i<muc_rvar_get_dim(rv);i++) {
		bounds[2*i] = rv->min[i];
		bounds[2*i + 1] = rv->max[i];
	}
	return 1;
}

