#ifndef lint
static const char RCSid[] = "$Id: paths.c,v 2.1 2003/06/08 12:01:52 schorsch Exp $";
#endif
/*
 * Find a writeable tempfile directory.
 * Create unique filenames therein, and possibly open the file.
 *
 */

#include "copyright.h"

#include "paths.h"


#define TEMPFILE_TEMPLATE "rtXXXXXX"
static const char *defaultpaths[] = DEFAULT_TEMPDIRS;


/* Return a writeable directory for temporary files */
/* If s is NULL, we return a static string */
char *
temp_directory(char *s, size_t len)
{
	static char td[PATH_MAX] = "\0"; /* remember */
	char * ts = NULL;
	int i = 0;
	
	if (td[0] != '\0') { /* we already have one */
		if (s == NULL) { return td; }
		strncpy(s, td, len);
		s[len-1] = '\0';
		return s;
	}
	/* check where TMP and TEMP point to */
	ts = getenv("TMP");
	if (ts != NULL && access(ts, W_OK) == 0) {
		strncpy(td, ts, sizeof(td));
		td[sizeof(td)-1] = '\0';
	}
	if (ts == NULL) {
		ts = getenv("TEMP");
		if (ts != NULL && access(ts, W_OK) == 0) {
			strncpy(td, ts, sizeof(td));
			td[sizeof(td)-1] = '\0';
		}
	}
	/* check the platform specific default paths in order */
	while (ts == NULL) {
		if (defaultpaths[i] == NULL) {
			break;
		}
		if (access(defaultpaths[i], W_OK) == 0) {
			ts = strncpy(td, defaultpaths[i], sizeof(td));
			td[sizeof(td)-1] = '\0';
			break;
		}
		i++;
	}
	/* we found something */
	if (ts != NULL) {
		if (s == NULL) { return td; }
		strncpy(s, ts, len);
		s[len-1] = '\0';
		return s;
	}
	return NULL;
}


/* Concatenate two strings, leaving exactly one DIRSEP in between */
char *
append_filepath(char *s1, char *s2, size_t len)
{
	size_t siz;
	char *s;

	siz = strlen(s1);
	/* XXX siz > len is an error */
	while (siz > 1 && ISDIRSEP(s1[siz-1])) {
		s1[siz-1] = '\0';
		siz--;
	}
	if (siz+1 <= len) {
		s1[siz] = DIRSEP;
		siz++;
	}
	while (ISDIRSEP(s2[0])) {
		s2++;
	}
	s = strncat(s1, s2, len-siz);
	return s;
}


/* Do the actual work for tempfiles, except for the uniquification */
static char *
prepare_tmpname(char *s, size_t len, char *templ)
{
	static char lp[PATH_MAX] = "\0"; /* remember what we found last time */
	char *ts = NULL;

	if (s == NULL) { /* return our static string */
		s = lp;
		len = sizeof(lp);
	}

	ts = temp_directory(s, len);
	if (ts == NULL) { return NULL; }
	ts = append_filepath(ts, templ != NULL ? templ : TEMPFILE_TEMPLATE, len);
	return ts;
}


/* Compose a *currently* unique name within a temporary directory */
/* If s is NULL, we return a static string */
/* If templ is NULL, we take our default template */
/* WARNING: On Windows, there's a maximum of 27 unique names within
            one process for the same template. */
char *
temp_filename(char *s, size_t len, char *templ)
{
	char *ts = NULL;

	ts = prepare_tmpname(s, len, templ);
	if (ts == NULL) { return NULL; }
	return mktemp(ts);
}


/* Open a unique temp file in a safe way (not safe on Windows) */
/* If s is NULL, we use a static string the caller won't learn about */
/* If templ is NULL, we take our default template */
/* This one is supposed to protect against race conditions on unix */
/* WARNING: On Windows, there's no protection against race conditions */
/* WARNING: On Windows, there's a maximum of 27 unique names within
            one process for the same template. */
int
temp_file(char *s, size_t len, char *templ)
{
	char *ts = NULL;

	ts = prepare_tmpname(s, len, templ);
	if (ts == NULL) return -1;
#ifdef _WIN32
	ts = mktemp(ts);
	if (ts == NULL) return -1;
	return fopen(ts, "r+b");
#else
	return mkstemp(ts);
#endif
}


#ifdef MODULE_TEST
int main()
{
	static char p[PATH_MAX] = "\0";
	char * pp, *qq = NULL;
	pp = temp_directory(p, sizeof(p));

	printf("%s\n", pp);
	qq = temp_filename(pp, sizeof(p), "//something/else_XXXXXX");
	printf("%s\n", qq);
	qq = temp_filename(pp, sizeof(p), NULL);
	printf("%s\n", qq);
}
#endif

