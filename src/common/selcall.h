/* RCSid $Id: selcall.h,v 3.7 2003/06/26 00:58:09 schorsch Exp $ */
/*
 * header file for select call compatibility
 */
#ifndef _RAD_SELCALL_H_
#define _RAD_SELCALL_H_
#ifdef __cplusplus
extern "C" {
#endif

#include "copyright.h"

#include <sys/types.h>
#ifdef _WIN32
  /*#include <winsock2.h>*/
#else
  #include <sys/time.h>
#endif
#ifdef INCL_SEL_H
#include <sys/select.h>
#endif

#ifndef FD_SETSIZE
#ifdef _WIN32
#else
  #include <sys/param.h>
#endif
#define FD_SETSIZE	NOFILE		/* maximum # select file descriptors */
#endif
#ifndef FD_SET
#ifndef NFDBITS
#define NFDBITS		(8*sizeof(int))	/* number of bits per fd_mask */
#endif
#define	FD_SET(n, p)	((p)->fds_bits[(n)/NFDBITS] |= (1 << ((n) % NFDBITS)))
#define	FD_CLR(n, p)	((p)->fds_bits[(n)/NFDBITS] &= ~(1 << ((n) % NFDBITS)))
#define	FD_ISSET(n, p)	((p)->fds_bits[(n)/NFDBITS] & (1 << ((n) % NFDBITS)))
#ifdef BSD
#define FD_ZERO(p)	bzero((char *)(p), sizeof(*(p)))
#else
#define FD_ZERO(p)	memset((char *)(p), 0, sizeof(*(p)))
#endif
#endif


#ifdef __cplusplus
}
#endif
#endif /* _RAD_SELCALL_H_ */
