/* RCSid $Id: tmaptiff.h,v 3.3 2003/06/06 16:38:47 schorsch Exp $ */
/*
 *  tmaptiff.h
 *
 *  Header file for TIFF tone-mapping routines.
 *  Include after "tiffio.h" and "tonemap.h".
 *
 */
#ifndef _RAD_TMAPTIFF_H_
#define _RAD_TMAPTIFF_H_
#ifdef __cplusplus
extern "C" {
#endif


#include "copyright.h"


extern int	tmCvL16(TMbright *ls, uint16 *luvs, int len);
extern int	tmCvLuv24(TMbright *ls, BYTE *cs, uint32 *luvs, int len);
extern int	tmCvLuv32(TMbright *ls, BYTE *cs, uint32 *luvs, int len);
extern int	tmLoadTIFF(TMbright **lpp, BYTE **cpp, int *xp, int *yp,
				char *fname, TIFF *tp);
extern int	tmMapTIFF(BYTE **psp, int *xp, int *yp, int flags,
				RGBPRIMP monpri, double gamval,
				double Lddyn, double Ldmax,
				char *fname, TIFF *tp);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_TMAPTIFF_H_ */

