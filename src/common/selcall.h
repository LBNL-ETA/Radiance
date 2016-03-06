/* RCSid $Id: selcall.h,v 3.11 2016/03/06 01:13:17 schorsch Exp $ */
/*
 * header file for select call compatibility
 */
#ifndef _RAD_SELCALL_H_
#define _RAD_SELCALL_H_

#include <string.h>
#include <sys/types.h>
#if defined(_WIN32) || defined(_WIN64)
  /*#include <winsock2.h>*/
#else
  #include <sys/time.h>
#endif
#ifdef INCL_SEL_H
#include <sys/select.h>
#endif

#ifndef FD_SETSIZE
#if defined(_WIN32) || defined(_WIN64)
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
#define FD_ZERO(p)	memset((char *)(p), 0, sizeof(*(p)))
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* nothing yet */

#ifdef __cplusplus
}
#endif
#endif /* _RAD_SELCALL_H_ */
