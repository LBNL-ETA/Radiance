#ifndef lint
static const char	RCSid[] = "$Id: getpath.c,v 2.17 2004/02/12 18:55:50 greg Exp $";
#endif
/*
 *  getpath.c - function to search for file in a list of directories
 *
 *  External symbols declared in standard.h
 */

#include "copyright.h"

#include  <string.h>
#include  <ctype.h>

#include  "rtio.h"
#include  "paths.h"


#ifdef _WIN32
static char *
core_getpath	/* wrapped below: expand fname, return full path */
#else
char *
getpath	/* expand fname, return full path */
#endif
(
register char  *fname,
register char  *searchpath,
int  mode
)
{
	static char  pname[PATH_MAX];
	char uname[512];
	register char  *cp;
	int i;

	if (fname == NULL) { return(NULL); }

	pname[0] = '\0';		/* check for full specification */

	if (ISABS(fname)) { /* absolute path */
		strncpy(pname, fname, sizeof(pname)-1);
	} else {
		switch (*fname) {
			case '.':				/* relative to cwd */
				strncpy(pname, fname, sizeof(pname)-1);
				break;
			case '~':				/* relative to home directory */
				fname++;
				cp = uname;
				for (i=0;i<sizeof(uname)&&*fname!='\0'&&!ISDIRSEP(*fname);i++)
					*cp++ = *fname++;
				*cp = '\0';
				cp = gethomedir(uname, pname, sizeof(pname));
				if(cp == NULL) return NULL;
				strncat(pname, fname, sizeof(pname)-strlen(pname)-1);
				break;
		}
	}
	if (pname[0])		/* got it, check access if search requested */
		return(searchpath==NULL||access(pname,mode)==0 ? pname : NULL);

	if (searchpath == NULL) {			/* don't search */
		strncpy(pname, fname, sizeof(pname)-1);
		return(pname);
	}
	/* check search path */
	do {
		cp = pname;
		while (*searchpath && (*cp = *searchpath++) != PATHSEP) {
			cp++;
		}
		if (cp > pname && !ISDIRSEP(cp[-1])) {
			*cp++ = DIRSEP;
		}
		strncpy(cp, fname, sizeof(pname)-strlen(pname)-1);
		if (access(pname, mode) == 0)		/* file accessable? */
			return(pname);
	} while (*searchpath);
	/* not found */
	return(NULL);
}


#ifdef _WIN32
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
register char  *ffname,
register char  *searchpath,
int  mode
)
{
	register char  *cp;
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

