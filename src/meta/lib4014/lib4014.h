/* RCSid: $Id: lib4014.h,v 1.1 2003/11/15 02:13:37 schorsch Exp $ */
/*
 *  random.h - header file for random(3) function.
 *
 *     10/1/85
 */
#ifndef _RAD_LIB4014_H_
#define _RAD_LIB4014_H_

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void openvt(void);
extern void openpl(void);
extern void closevt(void);
extern void closepl(void);

extern void step(int d);
extern void arc(int x, int y, int x0, int y0, int x1, int y1);
extern int quad(int x, int y, int xp, int yp);
extern int abs_(int a);
extern void box(int x0, int y0, int x1, int y1);
extern void circle(int x, int y, int r);
extern void dot(void);
extern void erase(void);
extern void label(char *s);
extern void line(int x0, int y0, int x1, int y1);
extern void linemod(char *s);
extern void move(int xi, int yi);
extern void point(int xi, int yi);
extern void scale(char i, float x, float y);
extern void space(int x0, int y0, int x1, int y1);
extern void cont(int x, int y);
extern void putch(int c);



#ifdef __cplusplus
}
#endif
#endif /* _RAD_LIB4014_H_ */

