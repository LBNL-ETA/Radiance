/* RCSid $Id: netproc.h,v 2.6 2003/06/27 11:32:12 schorsch Exp $ */
/*
 * Header file for network processing routines
 */
#ifndef _RAD_NETPROC_H_
#define _RAD_NETPROC_H_
#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
	char	*com;		/* command (stored on client side) */
	int	pid;		/* process ID */
	int	efd;		/* standard error file descriptor */
	char	*errs;		/* error output */
	int	elen;		/* error output length */
	int	(*cf)();	/* completion callback function */
} PROC;			/* process slot */

/* Callback function cf above passed process server, slot number and status */

typedef struct pserver {
	struct pserver	*next;	/* next process server in main list */
	char	hostname[64];	/* remote host ID */
	char	directory[128];	/* remote execution directory */
	char	username[32];	/* remote user ID */
	short	nprocs;		/* number of allocated process slots */
	PROC	proc[1];	/* process slot(s) (must be last in struct) */
} PSERVER;		/* process server */

extern PSERVER	*pslist;	/* global process server list */

extern PSERVER	*addpserver(char *host, char *dir, char *usr, int np);
extern PSERVER  *findjob(int *pnp);

extern char	*myhostname(void);

#define LHOSTNAME	"localhost"	/* accepted name for local host */

#define pserverOK(ps)	(wait4job(ps, startjob(ps, "true", NULL)) == 0)


#ifdef __cplusplus
}
#endif
#endif /* _RAD_NETPROC_H_ */

