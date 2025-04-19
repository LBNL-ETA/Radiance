#ifndef lint
static const char RCSid[] = "$Id: rmconvert.c,v 2.1 2025/04/19 17:12:59 greg Exp $";
#endif
/*
 * Convert between RMATRIX and CMATRIX types
 *
 *	Segregated because it depends on both sublibraries
 */

#include "rmatrix.h"

/* Convert a color matrix to newly allocated RMATRIX buffer */
RMATRIX *
rmx_from_cmatrix(const CMATRIX *cm)
{
	RMATRIX	*dnew;

	if (!cm)
		return(NULL);
	dnew = rmx_alloc(cm->nrows, cm->ncols, 3);
	if (!dnew)
		return(NULL);

	dnew->dtype = sizeof(COLORV)==sizeof(float) ?
			DTfloat : DTdouble;

	if (sizeof(COLORV) == sizeof(rmx_dtype)) {
		memcpy(dnew->mtx, cm->cmem, rmx_array_size(dnew));
	} else {
		int	i, j;
		for (i = dnew->nrows; i--; )
		    for (j = dnew->ncols; j--; ) {
			const COLORV	*cv = cm_lval(cm,i,j);
			rmx_dtype	*dp = rmx_lval(dnew,i,j);
			dp[0] = cv[0];
			dp[1] = cv[1];
			dp[2] = cv[2];
		    }
	}
	return(dnew);
}

/* Convert general matrix to newly allocated CMATRIX buffer */
CMATRIX *
cm_from_rmatrix(const RMATRIX *rm)
{
	CMATRIX	*cnew;

	if (!rm || !rm->mtx | (rm->ncomp == 2) | (rm->ncomp > MAXCOMP))
		return(NULL);
	cnew = cm_alloc(rm->nrows, rm->ncols);
	if (!cnew)
		return(NULL);
	if ((sizeof(COLORV) == sizeof(rmx_dtype)) & (rm->ncomp == 3)) {
		memcpy(cnew->cmem, rm->mtx, rmx_array_size(rm));
	} else {
		int	i, j;
		for (i = cnew->nrows; i--; )
		    for (j = cnew->ncols; j--; ) {
			const rmx_dtype	*dp = rmx_val(rm,i,j);
			COLORV		*cv = cm_lval(cnew,i,j);
			switch (rm->ncomp) {
			case 3:
			    setcolor(cv, dp[0], dp[1], dp[2]);
			    break;
			case 1:
			    setcolor(cv, dp[0], dp[0], dp[0]);
			    break;
			default:
			    if (sizeof(COLORV) == sizeof(rmx_dtype)) {
				scolor2color(cv, (const COLORV *)dp,
						rm->ncomp, rm->wlpart);
			    } else {
				COLORV	scol[MAXCOMP];
				int	k = rm->ncomp;
				while (k--) scol[k] = dp[k];
				scolor2color(cv, scol, rm->ncomp, rm->wlpart);
			    }
			    break;
			}
		    }
	}
	return(cnew);
}
