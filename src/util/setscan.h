/* RCSid: $Id: setscan.h,v 2.3 2003/06/27 11:32:12 schorsch Exp $ */
/*
 * Defines for programs using setscan()
 */
#ifndef _RAD_SETSCAN_H_
#define _RAD_SETSCAN_H_
#ifdef __cplusplus
extern "C" {
#endif


#define  ANGLE		short
#define  AEND		(-1)

int setscan(register ANGLE *ang, register char *arg);


#ifdef __cplusplus
}
#endif
#endif /* _RAD_SETSCAN_H_ */

