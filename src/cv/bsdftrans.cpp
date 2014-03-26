#ifndef lint
static const char RCSid[] = "$Id: bsdftrans.cpp,v 2.2 2014/03/26 22:29:08 greg Exp $";
#endif
/*
 * Compute mass transport plan for RBF lobes
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <string.h>
#include "transportSimplex.h"
#include "bsdfrep.h"

using namespace t_simplex;
using namespace std;

/* Compute mass transport plan (internal call) */
extern "C" void
plan_transport(MIGRATION *mig)
{
	double			src[mtx_nrows(mig)];
	double			dst[mtx_ncols(mig)];
	TsSignature<RBFVAL>	srcSig(mtx_nrows(mig), mig->rbfv[0]->rbfa, src);
	TsSignature<RBFVAL>	dstSig(mtx_ncols(mig), mig->rbfv[1]->rbfa, dst);
	TsFlow			flow[mtx_nrows(mig)+mtx_ncols(mig)-1];
	int			n;
						/* clear flow matrix */
	memset(mig->mtx, 0, sizeof(float)*mtx_nrows(mig)*mtx_ncols(mig));
						/* set supplies & demands */
	for (n = mtx_nrows(mig); n--; )
		src[n] = rbf_volume(&mig->rbfv[0]->rbfa[n]) /
					mig->rbfv[0]->vtotal;
	for (n = mtx_ncols(mig); n--; )
		dst[n] = rbf_volume(&mig->rbfv[1]->rbfa[n]) /
					mig->rbfv[1]->vtotal;
					
	n = 0;					/* minimize EMD */
	try {
		transportSimplex(&srcSig, &dstSig, &lobe_distance, flow, &n);
	} catch (...) {
		fprintf(stderr, "%s: caught exception from transportSimplex()!\n",
				progname);
		exit(1);
	}
	while (n-- > 0)				/* assign sparse matrix */
		mtx_coef(mig, flow[n].from, flow[n].to) = flow[n].amount;
}
