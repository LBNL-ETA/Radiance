#ifndef lint
static const char RCSid[] = "$Id: gethomedir.c,v 1.1 2003/10/27 10:20:58 schorsch Exp $";
#endif
/*
 *  gethomedir.c - search for a users home directory
 *
 */

#include  "copyright.h"

#include <stdlib.h>
#include <string.h>

#include "rtio.h"

#ifdef _WIN32

char *
gethomedir(char *uname, char *path, int plen)
{
	char *cd, *cp;

	if (uname == NULL || *uname == '\0') {	/* ours */
		/* pretend we're on unix first (eg. for Cygwin) */
		if ((cp = getenv("HOME")) != NULL) {
			strncpy(path, cp, plen);
			path[plen-1] = '\0';
			return path;
		}
		/* now let's see what Windows thinks */
		if ((cd = getenv("HOMEDRIVE")) != NULL
				&& (cp = getenv("HOMEPATH")) != NULL) {
			strncpy(path, cd, plen);
			strncat(path, cp, plen-2);
			path[plen-1] = '\0';
			return path;
		}
		return NULL;
	}
	/* No idea how to find the home directory of another user */
	return NULL;
}
	
#else /* _WIN32 */


#include <unistd.h>
#include <pwd.h>
#include <sys/types.h>

char *
gethomedir(char *uname, char *path, int plen)
{
	struct passwd *pwent;
	uid_t uid;
	char *cp;

	if (uname == NULL || *uname == '\0') {	/* ours */
		if ((cp = getenv("HOME")) != NULL) {
			strncpy(path, cp, plen);
			path[plen-1] = '\0';
			return path;
		}
		uid = getuid();
		if ((pwent = getpwuid(uid)) == NULL)
			return(NULL); /* we don't exist ?!? */
		strncpy(path, pwent->pw_dir, plen);
		path[plen-1] = '\0';
		return path;
	}
	/* someone else */
	if ((pwent = getpwnam(uname)) == NULL)
		return(NULL); /* no such user */

	strncpy(path, pwent->pw_dir, plen);
	path[plen-1] = '\0';
	return path;
}

#endif /* _WIN32 */

