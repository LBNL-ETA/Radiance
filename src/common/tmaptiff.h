/* RCSid $Id: tmaptiff.h,v 3.5 2005/01/07 20:33:02 greg Exp $ */
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

extern int	tmCvL16(TMstruct *tms, TMbright *ls, uint16 *luvs, int len);
extern int	tmCvLuv24(TMstruct *tms, TMbright *ls, BYTE *cs,
				uint32 *luvs, int len);
extern int	tmCvLuv32(TMstruct *tms, TMbright *ls, BYTE *cs,
				uint32 *luvs, int len);
extern int	tmLoadTIFF(TMstruct *tms, TMbright **lpp, BYTE **cpp,
				int *xp, int *yp, char *fname, TIFF *tp);
extern int	tmMapTIFF(BYTE **psp, int *xp, int *yp,
				int flags, RGBPRIMP monpri, double gamval,
				double Lddyn, double Ldmax,
				char *fname, TIFF *tp);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_TMAPTIFF_H_ */

