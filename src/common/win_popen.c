#ifndef lint
static const char RCSid[] = "$Id: win_popen.c,v 1.3 2004/03/28 20:33:12 schorsch Exp $";
#endif
/*
Replacement for the posix popen() on Windows
 
We don't just let the shell do the work like the original.
Since per default there's no decent shell around, we implement
the most basic functionality right here.

currently supports mode "r" and | to chain several processes.
Substrings between matching quotes are considered one word,
ignoring any | within. Quotes don't nest.
 
*/

#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <io.h>     /* _open_osfhandle()  */
#include <fcntl.h>  /* _O_RDONLY          */

#include "rtio.h"
#include "rterror.h"


#define RAD_MAX_PIPES 32 /* maximum number of pipes */

#define R_MODE 1
#define W_MODE 2
#define A_MODE 3

static int parse_pipes(char *s, char *lines[], char **infn, char **outfn,
int *append, int maxl);
static BOOL createPipes(HANDLE*, HANDLE*, HANDLE*, HANDLE*);
static BOOL runChild(char*, char*, HANDLE, HANDLE, HANDLE);
static void resetStdHandles(HANDLE stdoutOrig, HANDLE stdinOrig);
HANDLE newFile(char *fn, int mode);


int
win_pclose(    /* posix pclose replacement */
FILE* p
)
{
	fclose(p);
	/* not sure if it's useful to wait for anything on Windows */
	return 0;
}


FILE *
win_popen(     /* posix popen replacement */
char* command,
char* type
)
{
	char *execfile, *args;
	char *cmdlines[RAD_MAX_PIPES];
	char *infn = NULL, *outfn = NULL;
	char executable[512], estr[512];
	int n, i, mode = 0, append = 0;
	int ncmds = 0;
	FILE *inf = NULL;
	HANDLE stdoutRd = NULL, stdoutWr = NULL;
	HANDLE stdinRd = NULL, stdinWr = NULL;
	HANDLE stderrWr = NULL;
	HANDLE stdoutOrig, stdinOrig;

	if (strchr(type, 'w')) {
		mode = W_MODE;
	} else if (strchr(type, 'r')) {
		mode = R_MODE;
	} else {
		_snprintf(estr, sizeof(estr),
				"Invalid mode \"%s\" for win_popen().", type);
		eputs(estr);
		return NULL;
	}

	stdoutOrig = GetStdHandle(STD_OUTPUT_HANDLE);
	stdinOrig = GetStdHandle(STD_INPUT_HANDLE);
	/* if we have a console, use it for error output */
	stderrWr = GetStdHandle(STD_ERROR_HANDLE);

	ncmds = parse_pipes(command,cmdlines,&infn,&outfn,&append,RAD_MAX_PIPES);
	if(ncmds <= 0) {
		eputs("Too many pipes or malformed command.");
		goto error;
	}

	if (infn != NULL) {
		stdoutRd = newFile(infn, mode);
	}

	for(n = 0; n < ncmds; ++n) {
		if(!createPipes(&stdoutRd, &stdoutWr,
			&stdinRd, &stdinWr)) {
			eputs("Error creating pipe");
			goto error;
		}
		if (outfn != NULL && n == ncmds - 1) {
			CloseHandle(stdoutWr);
			CloseHandle(stdoutRd);
			stdoutWr = newFile(outfn, mode);
		}
		if (n == 0 && mode == W_MODE) {
			/* create a standard C file pointer for writing to the input */
			inf = _fdopen(_open_osfhandle((long)stdinWr, _O_RDONLY), "w");
		}
		/* find the executable on the PATH */
		args = nextword(executable, sizeof(executable), cmdlines[n]);
		if (args == NULL) {
			eputs("Empty command.");
			goto error;
		}
		execfile = getpath(executable, getenv("PATH"), X_OK);
		if(execfile == NULL) {
			_snprintf(estr, sizeof(estr),
					"Can't find executable for \"%s\".", executable);
			eputs(estr);
			goto error;
		}
		if(!runChild(execfile, cmdlines[n], stdinRd, stdoutWr, stderrWr)) {
			_snprintf(estr, sizeof(estr),
					"Unable to execute executable \"%s\".", executable);
			eputs(estr);
			goto error;
		}
		/* close the stdout end just passed to the last process,
		   or the final read will block */
		CloseHandle(stdoutWr);
	}

	/* clean up */
	resetStdHandles(stdinOrig, stdoutOrig);
	for (i = 0; i < ncmds; i++) free(cmdlines[i]);
	if (infn != NULL) free(infn);
	if (outfn != NULL) free(outfn);

	if (mode == R_MODE) {
		/* return a standard C file pointer for reading the output */
		return _fdopen(_open_osfhandle((long)stdoutRd, _O_RDONLY), "r");
	} else if (mode == W_MODE) {
		/* return a standard C file pointer for writing to the input */
		return inf;
	}

error:
	resetStdHandles(stdinOrig, stdoutOrig);
	for (i = 0; i < ncmds; i++) free(cmdlines[i]);
	if (infn != NULL) free(infn);
	if (outfn != NULL) free(outfn);
	return NULL;
}


HANDLE
newFile(char *fn, int mode)
{
	SECURITY_ATTRIBUTES sAttr;
	HANDLE fh;

	sAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
	sAttr.bInheritHandle = TRUE;
	sAttr.lpSecurityDescriptor = NULL;

	if (mode == R_MODE) {
		fh = CreateFile(fn,
				GENERIC_READ,
				FILE_SHARE_READ,
				&sAttr,
				OPEN_EXISTING,
				FILE_ATTRIBUTE_NORMAL,
				NULL);
	} else if (mode == W_MODE || mode == A_MODE) {
		fh = CreateFile(fn,
				GENERIC_WRITE,
				FILE_SHARE_WRITE,
				&sAttr,
				mode==W_MODE?CREATE_ALWAYS:OPEN_ALWAYS,
				FILE_ATTRIBUTE_NORMAL,
				NULL);
	}
	if (fh == NULL) {
		int e = GetLastError();
	}
	return fh;
}


static BOOL
createPipes(       /* establish matching pipes for a subprocess */
HANDLE* stdoutRd,
HANDLE* stdoutWr,
HANDLE* stdinRd,
HANDLE* stdinWr
)
{
	HANDLE stdoutRdInh = NULL;
	HANDLE stdinWrInh = NULL;
	HANDLE curproc;
	SECURITY_ATTRIBUTES sAttr;
	
	curproc = GetCurrentProcess();
	
	/* The rules of inheritance for handles are a mess.
	   Just to be safe, make all handles we pass to the
	   child processes inheritable */
	sAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
	sAttr.bInheritHandle = TRUE;
	sAttr.lpSecurityDescriptor = NULL;
	
	if(*stdoutRd != NULL){
	/* if we have a previous stdout pipe,
	   assign the other end to stdin for the new process */
		CloseHandle(*stdinRd);
		if(!DuplicateHandle(curproc, *stdoutRd,
			curproc, stdinRd, 0,
			TRUE, DUPLICATE_SAME_ACCESS))
			return FALSE;
		CloseHandle(*stdoutRd);
		if(!SetStdHandle(STD_INPUT_HANDLE, *stdinRd))
			return FALSE;
	} else {
		/* there's no previous stdout, create a new stdin pipe */
		if(!CreatePipe(stdinRd, &stdinWrInh, &sAttr, 0))
			return FALSE;
		if(!SetStdHandle(STD_INPUT_HANDLE, *stdinRd))
			return FALSE;
		CloseHandle(stdinWrInh);
	}
	
	/* create the stdout pipe for the new process */
	if(!CreatePipe(&stdoutRdInh, stdoutWr, &sAttr, 0))
		return FALSE;
	if(!SetStdHandle(STD_OUTPUT_HANDLE, *stdoutWr))
		return FALSE;
	if(!DuplicateHandle(curproc, stdoutRdInh,
		curproc, stdoutRd, 0,
		FALSE, DUPLICATE_SAME_ACCESS))
		return FALSE;
	CloseHandle(stdoutRdInh);
	
	return TRUE;
}


static void
resetStdHandles(    /* clean up our std streams */
HANDLE stdoutOrig,
HANDLE stdinOrig
)
{
	SetStdHandle(STD_OUTPUT_HANDLE, stdoutOrig);
	SetStdHandle(STD_INPUT_HANDLE, stdinOrig);
}


static int
parse_pipes(     /* split a shell command pipe sequence */
char *s,
char *lines[],
char **infn,
char **outfn,
int *append,
int maxl
)
{
	int n = 0, i;
	char *se, *ws;
	char *curs;
	int llen = 0;
	int quote = 0;
	int last = 0;

	if (maxl<= 0) return 0;
	if (s == NULL) {
		return 0;
	}
	*infn = *outfn = NULL;
	while (isspace(*s)) s++; /* leading whitespace */
	se = s;
	while (n < maxl) {
		switch (*se) {
			case '"':
				if (quote == '"') quote = 0;
				else if (quote == 0) quote = '"';
				break;
			case '\'':
				if (quote == '\'') quote = 0;
				else if (quote == 0) quote = '\'';
				break;
			case '<':
			case '>':
			case '|':
			case '\0':
				if (*se != '\0' && quote)
					break;
				llen = se - s;
				curs = malloc(llen+1);
				strncpy(curs, s, llen);
				/* remove unix style line-end escapes */
				while((ws = strstr(curs, "\\\n")) != NULL)
					*ws = *(ws+1) = ' ';
				/* remove DOS style line-end escapes */
				while((ws = strstr(curs, "\\\r\n")) != NULL)
					*ws = *(ws+1) = *(ws+2) = ' ';
				while (isspace(*(curs + llen - 1)))
					llen--; /* trailing whitespace */
				curs[llen] = '\0';

				if (last == '|' || last == 0) { /* first or pipe */
					lines[n] = curs;
					n++;
					curs = NULL;
				} else if (last == '<') { /* input file */
					if (*infn != NULL) {
						eputs("win_popen(): ambiguous input redirection");
						goto error;
					}
					*infn = curs;
					curs = NULL;
				} else if (last == '>') { /* output file */
					if (*outfn != NULL) {
						eputs("win_popen(): ambiguous out redirection");
						goto error;
					}
					*outfn = curs;
					curs = NULL;
					if (*se != '\0' && *se+1 == '>') { /* >> */
						*append = 1;
						se++;
					}
				}
				last = *se;

				if (*se == '\0') return n;
				s = se + 1;
				while (isspace(*s)) s++; /* leading whitespace */
				se = s;
				break;
			default:
				break;
		}
		se++;
	}
	/* more jobs than slots */
error:
	for (i = 0; i < n; i++) free(lines[i]);
	if (*infn != NULL) free(*infn);
	if (*outfn != NULL) free(*outfn);
	if (curs != NULL) free(curs);
	return -1;
}


static BOOL
runChild(         /* start a child process with the right std streams */
char* executable,
char* cmdline,
HANDLE stdinRd,
HANDLE stdoutWr,
HANDLE stderrWr
)
{
	PROCESS_INFORMATION procInfo;
	STARTUPINFO startupInfo;

	/* use the given handles and don't display the console window */
	ZeroMemory(&startupInfo, sizeof(STARTUPINFO));
	startupInfo.cb = sizeof(STARTUPINFO);
	startupInfo.dwFlags = STARTF_USESTDHANDLES | STARTF_USESHOWWINDOW;
	startupInfo.wShowWindow = SW_HIDE;
	startupInfo.hStdInput = stdinRd;
	startupInfo.hStdOutput = stdoutWr;
	startupInfo.hStdError = stderrWr;

	return CreateProcess(executable, cmdline, NULL, NULL,
		TRUE, 0,
		NULL, NULL, &startupInfo, &procInfo);
}



