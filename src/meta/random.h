/* RCSid: $Id: random.h,v 1.2 2003/07/14 22:24:00 schorsch Exp $ */
/*
 *  random.h - header file for random(3) function.
 *
 *     10/1/85
 */
#ifndef _RAD_META_RANDOM_H_
#define _RAD_META_RANDOM_H_

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define  frandom()	(random()/2147483648.0)

#ifdef __cplusplus
}
#endif
#endif /* _RAD_META_RANDOM_H_ */

