/* RCSid: $Id: mkillum.h,v 2.3 2003/06/26 00:58:09 schorsch Exp $ */
/*
 * Common definitions for mkillum
 */
#ifndef _RAD_MKILLUM_H_
#define _RAD_MKILLUM_H_
#ifdef __cplusplus
extern "C" {
#endif


#include  "standard.h"
#include  "object.h"
#include  "otypes.h"
#include  "rtprocess.h"

				/* illum flags */
#define  IL_LIGHT	0x1		/* light rather than illum */
#define  IL_COLDST	0x2		/* use color distribution */
#define  IL_COLAVG	0x4		/* use average color */
#define  IL_DATCLB	0x8		/* OK to clobber data file */

struct illum_args {
	int	flags;			/* flags from list above */
	char	matname[MAXSTR];	/* illum material name */
	char	datafile[MAXSTR];	/* distribution data file name */
	int	dfnum;			/* data file number */
	char	altmat[MAXSTR];		/* alternate material name */
	int	sampdens;		/* point sample density */
	int	nsamps;			/* # of samples in each direction */
	float	minbrt;			/* minimum average brightness */
	float	col[3];			/* computed average color */
};				/* illum options */

struct rtproc {
	SUBPROC pd;			/* rtrace pipe descriptors */
	float	*buf;			/* rtrace i/o buffer */
	int	bsiz;			/* maximum rays for rtrace buffer */
	float	**dest;			/* destination for each ray result */
	int	nrays;			/* current length of rtrace buffer */
};				/* rtrace process */


#ifdef __cplusplus
}
#endif
#endif /* _RAD_MKILLUM_H_ */

