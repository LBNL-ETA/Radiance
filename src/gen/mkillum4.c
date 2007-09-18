#ifndef lint
static const char RCSid[] = "$Id: mkillum4.c,v 2.1 2007/09/18 19:51:07 greg Exp $";
#endif
/*
 * Routines for handling BSDF data within mkillum
 */

#include "mkillum.h"
#include "paths.h"


struct BSDF_data *
load_BSDF(		/* load BSDF data from file */
	char *fname
)
{
	char			*path;
	FILE			*fp;
	struct BSDF_data	*dp;
	
	path = getpath(fname, getrlibpath(), R_OK);
	if (path == NULL) {
		sprintf(errmsg, "cannot find BSDF file \"%s\"", fname);
		error(WARNING, errmsg);
		return(NULL);
	}
	if ((fp = fopen(path, "r")) == NULL) {
		sprintf(errmsg, "cannot open BSDF \"%s\"", path);
		error(WARNING, errmsg);
		return(NULL);
	}
	dp = (struct BSDF_data *)malloc(sizeof(struct BSDF_data));
	if (dp == NULL)
		goto memerr;
	/* etc... */
	fclose(fp);
	return(dp);
memerr:
	error(SYSTEM, "out of memory in load_BSDF");
	return NULL;	/* pro forma return */
}


void
free_BSDF(		/* free BSDF data structure */
	struct BSDF_data *b
)
{
	if (b == NULL)
		return;
	free(b->inc_dir);
	free(b->inc_rad);
	free(b->out_dir);
	free(b->out_rad);
	free(b->bsdf);
	free(b);
}


void
r_BSDF_incvec(		/* compute random input vector at given location */
	FVECT v,
	struct BSDF_data *b,
	int i,
	double rv,
	MAT4 xm
)
{
	FVECT	pert;
	double	rad;
	int	j;
	
	getBSDF_incvec(v, b, i);
	rad = getBSDF_incrad(b, i);
	multisamp(pert, 3, rv);
	for (j = 0; j < 3; j++)
		v[j] += rad*(2.*pert[j] - 1.);
	if (xm != NULL)
		multv3(v, v, xm);
	normalize(v);
}


void
r_BSDF_outvec(		/* compute random output vector at given location */
	FVECT v,
	struct BSDF_data *b,
	int o,
	double rv,
	MAT4 xm
)
{
	FVECT	pert;
	double	rad;
	int	j;
	
	getBSDF_outvec(v, b, o);
	rad = getBSDF_outrad(b, o);
	multisamp(pert, 3, rv);
	for (j = 0; j < 3; j++)
		v[j] += rad*(2.*pert[j] - 1.);
	if (xm != NULL)
		multv3(v, v, xm);
	normalize(v);
}
