/* RCSid $Id: tmerrmsg.h,v 3.5 2003/06/06 16:38:47 schorsch Exp $ */
/*
 * Error messages for tone mapping functions.
 * Included exclusively in "tonemap.c".
 * English version.
 */
#ifndef _RAD_TMERRMSG_H_
#define _RAD_TMERRMSG_H_
#ifdef __cplusplus
extern "C" {
#endif

#include "copyright.h"

char	*tmErrorMessage[] = {
		"no error",
		"out of memory",
		"illegal argument value",
		"invalid tone mapping",
		"cannot compute tone mapping",
		"file cannot be opened or is the wrong type",
		"code consistency error 1",
		"code consistency error 2",
	};


#ifdef __cplusplus
}
#endif
#endif /* _RAD_TMERRMSG_H_ */

