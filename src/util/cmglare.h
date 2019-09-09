#ifndef _RAD_CMGLARE_H_
#define _RAD_CMGLARE_H_

#include "rtio.h"
#include "time.h"
#include "fvect.h"
#include "cmatrix.h"

#define DC_GLARE /* Perform glare autonomy calculation */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef DC_GLARE
#define TIMER(c, m) do {	\
	clock_t temp = clock(); \
	fprintf(stderr, "%s: %g s\n", m, 1.0 * (temp - c) / CLOCKS_PER_SEC); \
	c = temp; } while(0)
#else
#define TIMER(c, m)
#endif

float* cm_glare(const CMATRIX *dcmx, const CMATRIX *evmx, const CMATRIX *smx, const int *occupied, const double dgp_limit, const double dgp_threshold, const FVECT *views, const FVECT dir, const FVECT up);
int cm_load_schedule(const int count, int* schedule, FILE *fp);
FVECT* cm_load_views(const int nrows, const int inform, FILE *fp);
int cm_write_glare(const float *mp, const int nrows, const int ncols, const int dtype, FILE *fp);

#ifdef __cplusplus
}
#endif
#endif	/* _RAD_CMGLARE_H_ */
