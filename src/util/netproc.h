/* RCSid $Id: netproc.h,v 2.5 2003/06/27 06:53:23 greg Exp $ */
/*
 * Header file for network processing routines
 */

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

extern PSERVER	*addpserver(), *findjob();

extern char	*myhostname();

#define LHOSTNAME	"localhost"	/* accepted name for local host */

#define pserverOK(ps)	(wait4job(ps, startjob(ps, "true", NULL)) == 0)
