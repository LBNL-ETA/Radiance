/* RCSid $Id: tmerrmsg.h,v 3.6 2003/06/27 06:53:22 greg Exp $ */
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

