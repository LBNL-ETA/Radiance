#ifndef lint
static const char RCSid[] = "$Id: gethomedir.c,v 1.4 2019/12/28 18:05:14 greg Exp $";
#endif
/*
 *  gethomedir.c - search for a users home directory
 *
 */

#include  "copyright.h"

#include <stdlib.h>

#include "rtio.h"

#if defined(_WIN32) || defined(_WIN64)

char *
gethomedir(char *uname, char *path, int plen)
{
	char *cd, *cp;

	if (uname == NULL || *uname == '\0') {	/* ours */
		/* pretend we're on unix first (eg. for Cygwin) */
		if ((cp = getenv("HOME")) != NULL) {
			strlcpy(path, cp, plen);
			return path;
		}
		/* now let's see what Windows thinks */
		if ((cd = getenv("HOMEDRIVE")) != NULL
				&& (cp = getenv("HOMEPATH")) != NULL) {
			strlcpy(path, cd, plen);
			strlcat(path, cp, plen);
			return path;
		}
		return NULL;
	}
	/* No idea how to find the home directory of another user */
	return NULL;
}
	
#else /* _WIN32 || _WIN64 */


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
			strlcpy(path, cp, plen);
			return path;
		}
		uid = getuid();
		if ((pwent = getpwuid(uid)) == NULL)
			return(NULL); /* we don't exist ?!? */
		strlcpy(path, pwent->pw_dir, plen);
		return path;
	}
	/* someone else */
	if ((pwent = getpwnam(uname)) == NULL)
		return(NULL); /* no such user */

	strlcpy(path, pwent->pw_dir, plen);
	return path;
}

#endif /* _WIN32 || _WIN64 */

