/* RCSid $Id: font.h,v 2.7 2003/06/27 06:53:21 greg Exp $ */
/*
 * Header file for font handling routines
 */
#ifndef _RAD_FONT_H_
#define _RAD_FONT_H_
#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned char  GORD;

typedef struct {
	short  nverts;			/* number of vertices */
	GORD  left, right, top, bottom;	/* glyph extent */
					/* followed by vertex list */
}  GLYPH;

#define gvlist(g)	((GORD *)((g)+1))

typedef struct font {
	int  nref;			/* number of references */
	GLYPH  *fg[256];		/* font glyphs */
	short  mwidth, mheight;		/* mean glyph width and height */
	char  *name;			/* font file name */
	struct font  *next;		/* next font in list */
}  FONT;

extern int	retainfonts;		/* retain loaded fonts? */


extern FONT  *getfont(char *fname);
extern void  freefont(FONT *f);
extern int  uniftext(short *sp, char *tp, FONT *f);
extern int  squeeztext(short *sp, char *tp, FONT *f, int cis);
extern int  proptext(short *sp, char *tp, FONT *f, int cis, int nsi);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_FONT_H_ */

