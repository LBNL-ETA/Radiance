#ifndef lint
static const char RCSid[] = "$Id: win_netproc.c,v 2.4 2004/10/23 18:55:53 schorsch Exp $";
#endif

/*
 * Skeleton netproc for WIN32
 */

/*
#include <windows.h>
*/
#include <direct.h>

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>

#include "rtprocess.h"
#include "netproc.h"
#include "paths.h"

static int maxfd;
static char	ourdir[PATH_MAX];
static char	ourhost[64];
static char	*ourshell;
static char	ouruser[32];
PSERVER	*pslist = NULL;

PSERVER*
addpserver(
	char* host,
	char* dir,
	char* usr,
	int np
)
{
	PSERVER* ps;
	if (np < 1)
		return(NULL);
	// For now
	if(strcmp(host, LHOSTNAME))
		return(NULL);
	if((ps = (PSERVER *)malloc(sizeof(PSERVER)+(np-1)*sizeof(NETPROC))) == NULL)
		return(NULL);
	if (!ourhost[0]) {
		char	dirtmp[PATH_MAX];
		register char	*cp;
		register int	len;

		strcpy(ourhost, myhostname());
		getcwd(dirtmp, sizeof(dirtmp));
		if ((cp = getenv("HOME")) != NULL) {
			if (!strcmp(cp, dirtmp))
				ourdir[0] = '\0';
			else if (!strncmp(cp, dirtmp, len=strlen(cp)) &&
					dirtmp[len] == '/')
				strcpy(ourdir, dirtmp+len+1);
			else
				strcpy(ourdir, dirtmp);
		} else
			strcpy(ourdir, dirtmp);
		if ((cp = getenv("USER")) != NULL)
			strcpy(ouruser, cp);
		if ((ourshell = getenv("SHELL")) == NULL)
			ourshell = "/bin/sh";
		//FD_ZERO(&errdesc);
		maxfd = -1;
	}

	if (host == NULL || !strcmp(host, ourhost) ||
			!strcmp(host, LHOSTNAME))
		ps->hostname[0] = '\0';
	else
		strcpy(ps->hostname, host);
	if (dir == NULL)
		strcpy(ps->directory, ourdir);
	else
		strcpy(ps->directory, dir);
	if (usr == NULL || !strcmp(usr, ouruser))
		ps->username[0] = '\0';
	else
		strcpy(ps->username, usr);
	ps->nprocs = np;
	while (np--) {
		ps->proc[np].com = NULL;
		ps->proc[np].pid = -1;
		ps->proc[np].efd = -1;
		ps->proc[np].errs = NULL;
		ps->proc[np].elen = 0;
		ps->proc[np].cf = NULL;
	}
	
	ps->next = pslist;
	pslist = ps;
	
	if (!pserverOK(ps)) {
		delpserver(ps);
		return(NULL);
	}
	return(ps);
}

void delpserver(PSERVER *ps)
{
	return;
}

PSERVER* findjob(int *pnp)
{
	static PSERVER ps;
	return &ps;
}


int startjob(
	PSERVER* pPS,
	char* command,
	pscompfunc *compf
)
{
	return 0;
}

int wait4job(PSERVER* pPS, int pID) {
	if(pID == -1)
		return -1;
	return 0;
}
