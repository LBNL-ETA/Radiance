/* RCSid: $Id: clrtab.h,v 2.2 2011/05/20 02:06:39 greg Exp $ */
/* color table routines
*/

#ifndef _RAD_CLRTAB_H_
#define _RAD_CLRTAB_H_

#ifdef __cplusplus
extern "C" {
#endif

	/* defined in clrtab.c */
extern int new_histo(int n);
extern void cnt_pixel(uby8 col[]);
extern void cnt_colrs(COLR *cs, int n);
extern int new_clrtab(int ncolors);
extern int map_pixel(uby8 col[]);
extern void map_colrs(uby8 *bs, COLR *cs, int n);
extern void dith_colrs(uby8 *bs, COLR *cs, int n);

	/* defined in neuclrtab.c */
extern int neu_init(long npixels);
extern void neu_pixel(uby8 col[]);
extern void neu_colrs(COLR *cs, int n);
extern int neu_clrtab(int ncolors);
extern int neu_map_pixel(uby8 col[]);
extern void neu_map_colrs(uby8 *bs, COLR *cs, int n);
extern void neu_dith_colrs(uby8 *bs, COLR *cs, int n);



#ifdef __cplusplus
}
#endif
#endif /* _RAD_CLRTAB_H_ */

