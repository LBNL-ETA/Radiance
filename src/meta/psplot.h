/* RCSid: $Id: psplot.h,v 1.1 2003/11/15 02:13:37 schorsch Exp $ */
/*
 *   Structures for line segment output to dot matrix printers
 */
#ifndef _RAD_PSPLOT_H_
#define _RAD_PSPLOT_H_

#ifdef __cplusplus
extern "C" {
#endif

extern void init(char *s);
extern void done(void);
extern void segopen(char *s);
extern void segclose(void);
extern void doseg(register PRIMITIVE	*p);
extern void plotvstr(register PRIMITIVE	*p);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_PSPLOT_H_ */

