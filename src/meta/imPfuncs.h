/* RCSid: $Id: imPfuncs.h,v 1.4 2004/03/28 20:33:13 schorsch Exp $ */
/*
 *  imPfuncs - imPress functions
 *
 *  Written by William LeFebvre, LCSE, Rice University
 *
 *  This program can be freely redistributed to anyone, with the following
 *  provisions:  that this comment remain intact, and that no money is
 *  charged or collected for that redistribution.
 */

/*
 *  These functions are part of the imPress functional interface routines.
 */
#ifndef _RAD_IMPFUNCS_H_
#define _RAD_IMPFUNCS_H_

#include <varargs.h>

#ifdef __cplusplus
extern "C" {
#endif

#define im_putbyte(c)	putc((c) & 0377, imout)
#define im_putword(w)	(im_putbyte((w) >> 8), im_putbyte(w))

#define im_b(code, b)	(putc(code, imout), im_putbyte(b))

/* function definitions for imPress commands */

#define imBrule(wa, wb, wc)	im_www(imP_BRULE, wa, wb, wc)
#define imCircArc(wa, wb, wc)	im_www(imP_CIRC_ARC, wa, wb, wc)
#define imCircSegm(wa, wb, wc, wd, we) \
				im_wwwww(imP_CIRC_SEGM, wa, wb, wc, wd, we)
#define imCrLf()		putc(imP_CRLF, imout)
#define imDrawPath(b)		im_b(imP_DRAW_PATH, b)
#define imEllipseArc(wa, wb, wc, wd, we) \
				im_wwwww(imP_ELLIPSE_ARC, wa, wb, wc, wd, we)
#define imEndPage()		putc(imP_ENDPAGE, imout)
#define imEof()			putc(imP_EOF, imout)
#define imFillPath(b)		im_b(imP_FILL_PATH, b)
#define imMminus()		putc(imP_MMINUS, imout)
#define imMmove(w)		im_w(imP_MMOVE, w)
#define imMplus()		putc(imP_MPLUS, imout)
#define imMplus()		putc(imP_MPLUS, imout)
#define imNoOp()		putc(imP_NO_OP, imout)
#define imPage()		putc(imP_PAGE, imout)
#define imPop()			putc(imP_POP, imout)
#define imPush()		putc(imP_PUSH, imout)
#define imSetAbsH(w)		im_w(imP_SET_ABS_H, w)
#define imSetAbsV(w)		im_w(imP_SET_ABS_V, w)
#define imSetAdvDirs(b2, b1)	im_21(imP_SET_ADV_DIRS, b2, b1)
#define imSetBol(w)		im_w(imP_SET_BOL, w)
#define imSetFamily(b)		im_b(imP_SET_FAMILY, b)
#define imSetHVSystem(b2o, b2a, b3)   im_223(imP_SET_HV_SYSTEM, b2o, b2a, b3)
#define imSetIl(w)		im_w(imP_SET_IL, w)
#define imSetMagnification(b)	im_b(imP_SET_MAGN, b)
#define imSetPen(b)		im_b(imP_SET_PEN, b)
#define imSetPum(b1)		im_b(imP_SET_PUM, b1)
#define imSetPushMask(w)	im_w(imP_SET_PUSH_MASK, w)
#define imSetRelH(w)		im_w(imP_SET_REL_H, w)
#define imSetRelV(w)		im_w(imP_SET_REL_V, w)
#define imSetSp(w)		im_w(imP_SET_SP, w)
#define imSetTexture(b7a, b7b)	im_77(imP_SET_TEXTURE, b7a, b7b)
#define imSmove(w)		im_w(imP_SMOVE, w)
#define imSp()			putc(imP_SP, imout)
#define imSp1()			putc(imP_SP1, imout)

/* Push mask defines: */

#define	imPM_PEN_TEXTURE	0x100
#define	imPM_SP			0x080
#define	imPM_IL			0x040
#define	imPM_BOL		0x020
#define	imPM_FAMILY		0x010
#define	imPM_HV_POS		0x008
#define	imPM_ADV_DIRS		0x004
#define	imPM_ORIGIN		0x002
#define	imPM_ORIENTATION	0x001

extern FILE *imout;

extern void im_w(int code, unsigned int w);
extern void im_223(int code,unsigned int b2a,unsigned int b2b,unsigned int b3);
extern void im_77(int code, unsigned int b7a, unsigned int b7b);
/* extern void imCreateFamilyTable(char *s, ...); */ /* XXX */
/* extern void imCreatePath(char *s, ...); */ /* XXX */
extern void imCreatePathV(unsigned int count, unsigned int *vec);
extern void imInit(void);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_IMPFUNCS_H_ */

