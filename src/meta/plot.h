/* RCSid: $Id: plot.h,v 1.4 2003/11/15 02:13:37 schorsch Exp $ */
/*
 *   Definitions for plotting routines
 */
#ifndef _RAD_PLOT_H_
#define _RAD_PLOT_H_

#ifdef __cplusplus
extern "C" {
#endif

#define  NPATS  9		/* number of fill patterns */

#define  PATSIZE  16		/* pattern size (square) */


extern int  dxsize, dysize;	/* device size */

extern int  pati[];

extern unsigned char  pattern[][PATSIZE/8][PATSIZE];	/* fill patterns */

	/* macplot.c, plot.c, psplot.c */
extern void set(int attrib, char *value);
extern void reset(int attrib);
extern void unset(int attrib);

	/* xxxplot.c */
extern void plot(FILE *fp);
extern void fillpoly(PRIMITIVE *p);
extern void filltri(PRIMITIVE *p);
extern void printstr(PRIMITIVE *p);
extern void fillrect(PRIMITIVE *p);
extern void plotlseg(PRIMITIVE *p);

	/* meta2tga.c <-> rplot.c */
extern void nextblock(void);
extern void outputblock(void);
extern void printblock(void);

	/* primout.c */
extern void pglob(int co, int a0, char *s);
extern void pprim(int co,int a0,int xmin,int ymin,int xmax,int ymax,char *s);
extern void plseg(int a0, int xstart, int ystart, int xend, int yend);

	/* psplot, x11plot.c, xplot.c */
extern void thispage(void);
extern void nextpage(void);
extern void contpage(void);
extern void printspan(void);
extern void endpage(void);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_PLOT_H_ */

