#ifndef lint
static const char	RCSid[] = "$Id: getpath.c,v 2.21 2019/12/28 18:05:14 greg Exp $";
#endif
/*
 *  getpath.c - function to search for file in a list of directories
 *
 *  External symbols declared in rtio.h
 */

#include "copyright.h"

#include  <ctype.h>

#include  "rtio.h"
#include  "paths.h"


#if defined(_WIN32) || defined(_WIN64)
static char *
core_getpath	/* wrapped below: expand fname, return full path */
#else
char *
getpath	/* expand fname, return full path */
#endif
(
	char  *fname,
	char  *searchpath,
	int  mode
)
{
	static char  pname[PATH_MAX];
	char uname[512];
	char  *cp;
	int i;

	if (fname == NULL) { return(NULL); }

	pname[0] = '\0';		/* check for full specification */

	if (ISABS(fname)) { /* absolute path */
		strlcpy(pname, fname, sizeof(pname));
	} else {
		switch (*fname) {
			case '.':		/* relative to cwd */
				strlcpy(pname, fname, sizeof(pname));
				break;
			case '~':		/* relative to home directory */
				fname++;
				cp = uname;
				for (i = 0; i < sizeof(uname) && *fname
						&& !ISDIRSEP(*fname); i++)
					*cp++ = *fname++;
				*cp = '\0';
				cp = gethomedir(uname, pname, sizeof(pname));
				if(cp == NULL) return NULL;
				strlcat(pname, fname, sizeof(pname));
				break;
		}
	}
	if (pname[0])		/* got it, check access if search requested */
		return(!searchpath || access(pname,mode)==0 ? pname : NULL);

	if (!searchpath) {			/* no search */
		strlcpy(pname, fname, sizeof(pname));
		return(pname);
	}
	/* check search path */
	do {
		cp = pname;
		while (*searchpath && (*cp = *searchpath++) != PATHSEP)
			cp += (cp-pname < sizeof(pname)-2);
		if (cp > pname && !ISDIRSEP(cp[-1]))
			*cp++ = DIRSEP;
		*cp = '\0';
		strlcat(pname, fname, sizeof(pname));
		if (access(pname, mode) == 0)		/* file accessable? */
			return(pname);
	} while (*searchpath);
	/* not found */
	return(NULL);
}


#if defined(_WIN32) || defined(_WIN64)
/* This is a wrapper around the above, "emulating" access mode X_OK,
   which is not supported on Windows.
   If we see X_OK and the filename has no extension, then we'll remove
   the X_OK from the mode, append ".exe" to the file name, and search
   with the resulting string. If that fails, we'll try again with ".bat".
   Theoretically, we might still not have execute rights on a file we find
   like this, but that's rare enough not to be worth checking the ACLs.
*/
char *
getpath(	/* expand fname, return full path */
	char  *ffname,
	char  *searchpath,
	int  mode
)
{
	char  *cp;
	char fname[PATH_MAX];

	if (ffname == NULL)
		return(NULL);

	/* if we have a dot in the string, we assume there is a file name
	   extension present */
	/* XXX We'd better test for .exe/.bat/.etc explicitly */
	if (!(mode & X_OK) || (strrchr(ffname, '.') != NULL)) {
		return core_getpath(ffname, searchpath, mode);
	}

	/* We're looking for an executable, Windows doesn't have X_OK. */
	mode &= ~X_OK;
	/* Append .exe */
	strncpy(fname, ffname, sizeof(fname)-5);
	strcat(fname, ".exe");
	cp = core_getpath(fname, searchpath, mode);
	if (cp != NULL) return cp;

	/* Try with .bat this time... */
	strncpy(fname, ffname, sizeof(fname)-5);
	strcat(fname, ".bat");
	cp = core_getpath(fname, searchpath, mode);
	return cp;
}
#endif /* _WIN32 */


#ifdef TEST_MODULE
int main()
{
	char * fp;
	char fmt[] = "%15s %-10s %s: %s\n";

	fp = getpath("rayinit.cal", getenv("RAYPATH"), R_OK);
	printf(fmt,  "rayinit.cal",        "RAYPATH", "R_OK", fp);
	fp = getpath("mkillum", getenv("PATH"), X_OK);
	printf(fmt,  "mkillum",        "PATH", "X_OK", fp);
	fp = getpath("/", getenv("PATH"), W_OK);
	printf(fmt,  "/",        "PATH", "W_OK", fp);
	fp = getpath("~", getenv("PATH"), F_OK);
	printf(fmt,  "~",        "PATH", "F_OK", fp);
	printf("Undefining HOME and HOMEPATH\n");
	unsetenv("HOME");
	unsetenv("HOMEPATH");
	fp = getpath("~", getenv("PATH"), F_OK);
	printf(fmt,  "~",        "PATH", "F_OK", fp);
	fp = getpath("~lp/blah", getenv("PATH"), F_OK);
	printf(fmt, "~lp/blah",         "PATH", "F_OK", fp);
}
#endif

