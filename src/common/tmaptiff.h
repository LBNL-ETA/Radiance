/* RCSid $Id: tmaptiff.h,v 3.2 2003/02/25 02:47:22 greg Exp $ */
/*
 *  tmaptiff.h
 *
 *  Header file for TIFF tone-mapping routines.
 *  Include after "tiffio.h" and "tonemap.h".
 *
 */

#include "copyright.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef NOPROTO

extern int	tmCvL16();
extern int	tmCvLuv24();
extern int	tmCvLuv32();
extern int	tmLoadTIFF();
extern int	tmMapTIFF();

#else

extern int	tmCvL16(TMbright *ls, uint16 *luvs, int len);
extern int	tmCvLuv24(TMbright *ls, BYTE *cs, uint32 *luvs, int len);
extern int	tmCvLuv32(TMbright *ls, BYTE *cs, uint32 *luvs, int len);
extern int	tmLoadTIFF(TMbright **lpp, BYTE **cpp, int *xp, int *yp,
				char *fname, TIFF *tp);
extern int	tmMapTIFF(BYTE **psp, int *xp, int *yp, int flags,
				RGBPRIMP monpri, double gamval,
				double Lddyn, double Ldmax,
				char *fname, TIFF *tp);

#endif

#ifdef __cplusplus
}
#endif
