#ifndef lint
static const char	RCSid[] = "$Id: getpath.c,v 2.13 2003/07/03 20:59:00 schorsch Exp $";
#endif
/*
 *  getpath.c - function to search for file in a list of directories
 *
 *  External symbols declared in standard.h
 */

#include "copyright.h"

#include  <string.h>
#include  <ctype.h>
#ifdef	 RHAS_GETPWNAM
  #include  <pwd.h>
  #include  <sys/types.h>
#endif

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
#ifdef	 RHAS_GETPWNAM
	struct passwd  *pwent;
#endif
	static char  pname[PATH_MAX];
	register char  *cp;

	if (fname == NULL) { return(NULL); }

	pname[0] = '\0';		/* check for full specification */

	if (ISABS(fname)) { /* Can't use CASEDIRSEP below on Windows */
		strncpy(pname, fname, sizeof(pname)-1);
	} else {
		switch (*fname) {
			/* XXX This doesn't work on Windows */
			/* CASEDIRSEP: */			/* relative to root */
			case '.':				/* relative to cwd */
				strncpy(pname, fname, sizeof(pname)-1);
				break;
			case '~':				/* relative to home directory */
				fname++;
				if (*fname == '\0' || ISDIRSEP(*fname)) {	/* ours */
					if ((cp = getenv("HOME")) == NULL)
#ifdef _WIN32  /* Windows sometimes uses a different var name */
						if ((cp = getenv("HOMEDIR")) == NULL)
#endif
							return(NULL);
					strncpy(pname, cp, sizeof(pname)-1);
					strncat(pname, fname, sizeof(pname)-strlen(pname)-1);
					break;
				}
#ifdef RHAS_GETPWNAM
				/* XXX Should we request our own home directory from the
				   XXX system as well if the above fails? */
				/* XXX In any case, we need do the same thing on Windows... */
				cp = pname;					/* user */
				do
					*cp++ = *fname++;
				while (*fname && !ISDIRSEP(*fname));
				*cp = '\0';
				if ((pwent = getpwnam(pname)) == NULL)
					return(NULL);
				strncpy(pname, pwent->pw_dir, sizeof(pname)-1);
				strncat(pname, fname, sizeof(pname)-strlen(pname)-1);
				break;
#endif
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
	printf("Undefining HOME and HOMEDIR\n");
	unsetenv("HOME");
	unsetenv("HOMEDIR");
	fp = getpath("~", getenv("PATH"), F_OK);
	printf(fmt,  "~",        "PATH", "F_OK", fp);
#ifndef RHAS_GETPWNAM
	fp = getpath("~lp", getenv("PATH"), F_OK);
	printf(fmt, "~lp",         "PATH", "F_OK", fp);
#endif
}
#endif

