/* RCSid $Id: oconv.h,v 2.1 2004/03/27 12:41:45 schorsch Exp $ */
/*
 *  Header for oconv support routines
 *
 */
#ifndef _RAD_OCONV_H_
#define _RAD_OCONV_H_

#include "octree.h"
#include "object.h"

#ifdef __cplusplus
extern "C" {
#endif

	/* defined in initotypes.c */
void ot_initotypes(void);

	/* defined in writeoct.c */
extern void writeoct(int  store, CUBE  *scene, char  *ofn[]);

	/* defined in bbox.c */
extern void add2bbox(OBJREC  *o, FVECT  bbmin, FVECT  bbmax);

	/* defined in readobj2.c */
typedef void ro_cbfunc(OBJREC *o);
extern void readobj2(char  *input, ro_cbfunc callback);



#ifdef __cplusplus
}
#endif
#endif /* _RAD_OCONV_H_ */

