/* RCSid $Id: paths.h,v 2.12 2003/06/06 16:38:47 schorsch Exp $ */
/*
 * Definitions for paths on different machines
 */
#ifndef _RAD_PATHS_H_
#define _RAD_PATHS_H_
#ifdef __cplusplus
extern "C" {
#endif

#include "copyright.h"

#ifdef _WIN32

#define DIRSEP		'/'
#define ISDIRSEP(c)	((c)=='/' || (c)=='\\')
#define CASEDIRSEP	case '/': case '\\'
#define PATHSEP		';'
#define MAXPATH		128
#define TEMPLATE	"rtXXXXXX"
#define TEMPLEN		8
#define ULIBVAR		"RAYPATH"
#ifndef DEFPATH
#define DEFPATH		";c:/ray/lib"
#endif
extern char  *fixargv0();

#else
#ifdef AMIGA

#define DIRSEP		'/'
#define PATHSEP		';'
#define MAXPATH		128
#define TEMPLATE	"/tmp/rtXXXXXX"
#define TEMPLEN		13
#define ULIBVAR		"RAYPATH"
#ifndef DEFPATH
#define DEFPATH		";/ray/lib"
#endif
#define	 fixargv0(a0)	(a0)

#else

#define DIRSEP		'/'
#define PATHSEP		':'
#define MAXPATH		256
#define TEMPLATE	"/usr/tmp/rtXXXXXX"
#define TEMPLEN		17
#define ULIBVAR		"RAYPATH"
#ifndef DEFPATH
#define DEFPATH		":/usr/local/lib/ray"
#endif
#define	 fixargv0(a0)	(a0)

#endif
#endif

#ifndef ISDIRSEP
#define ISDIRSEP(c)	((c)==DIRSEP)
#endif
#ifndef CASEDIRSEP
#define CASEDIRSEP	case DIRSEP
#endif

extern char  *mktemp(), *getenv();

#ifdef BSD
extern char  *getwd();
#else
extern char  *getcwd();
#define  getwd(p)	getcwd(p, sizeof(p))
#endif


#ifdef __cplusplus
}
#endif
#endif /* _RAD_PATHS_H_ */

