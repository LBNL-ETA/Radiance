#ifndef lint
static char RCSid[]="$Id: win_process.c,v 3.4 2003/10/21 19:19:28 schorsch Exp $";
#endif
/*
 * Routines to communicate with separate process via dual pipes.
 * Windows version.
 *
 * External symbols declared in standard.h
 */

#include "copyright.h"

#include <stdio.h>
#define STRICT
#include <windows.h> /* typedefs */
#include <io.h>      /* _open_osfhandle */
#include <fcntl.h>   /* _O_XXX */

#include "rterror.h"
#include "rtio.h"
#include "rtprocess.h"


int
win_nice(int inc) /* simple nice(2) replacement for Windows */
{
	/* We don't have much granularity available: IDLE_PRIORITY_CLASS
	   will run whenever no other higher priority process is running */
	if (inc > 0) {
		return (int)!SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS);
	}
	return 0;
}


/*
    Safely terminate a process by creating a remote thread
    in the process that calls ExitProcess.
	As presented by Andrew Tucker in Windows Developer Magazine.
*/
#ifndef OBSOLETE_WINDOWS  /* won't work on Win 9X/ME/CE. */
BOOL SafeTerminateProcess(HANDLE hProcess, UINT uExitCode)
{
	DWORD dwTID, dwCode, dwErr = 0;
	HANDLE hProcessDup = INVALID_HANDLE_VALUE;
	HANDLE hRT = NULL;
	HINSTANCE hKernel = GetModuleHandle("Kernel32");
	BOOL bSuccess = FALSE;

	BOOL bDup = DuplicateHandle(GetCurrentProcess(), 
			hProcess, 
			GetCurrentProcess(), 
			&hProcessDup, 
			PROCESS_ALL_ACCESS, 
			FALSE, 
			0);
	/* Detect the special case where the process is already dead... */
	if ( GetExitCodeProcess((bDup) ? hProcessDup : hProcess, &dwCode) && 
			(dwCode == STILL_ACTIVE) ) {
		FARPROC pfnExitProc;
		pfnExitProc = GetProcAddress(hKernel, "ExitProcess");
		hRT = CreateRemoteThread((bDup) ? hProcessDup : hProcess, 
				NULL, 
				0, 
				(LPTHREAD_START_ROUTINE)pfnExitProc,
				(PVOID)uExitCode, 0, &dwTID);
		if ( hRT == NULL ) dwErr = GetLastError();
	} else {
		dwErr = ERROR_PROCESS_ABORTED;
	}
	if ( hRT ) {
		/* Must wait process to terminate to guarantee that it has exited... */
		WaitForSingleObject((bDup) ? hProcessDup : hProcess, INFINITE);
		CloseHandle(hRT);
		bSuccess = TRUE;
	}
	if ( bDup ) CloseHandle(hProcessDup);
	if ( !bSuccess ) SetLastError(dwErr);
	return bSuccess;
}
#endif


static int
start_process(SUBPROC *proc, char *cmdstr)
{
	BOOL res;
	int Pflags = 0;
	SECURITY_ATTRIBUTES SAttrs;
	STARTUPINFO SInfo;
	PROCESS_INFORMATION PInfo;
	/* welcome to the world of resource handles */
	HANDLE hToChildRead = NULL;
	HANDLE hToChildWrite = NULL;
	HANDLE hFromChildRead = NULL;
	HANDLE hFromChildWrite = NULL;
	HANDLE hRead = NULL, hWrite = NULL;
	HANDLE hStdIn, hStdOut, hStdErr;
	HANDLE hCurProc;

	/* get process and standard stream handles */
	hCurProc = GetCurrentProcess();
	hStdIn = GetStdHandle(STD_INPUT_HANDLE);
	hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
	hStdErr = GetStdHandle(STD_ERROR_HANDLE);
	
	/* the remote pipe handles must be inheritable */
	SAttrs.bInheritHandle = 1;
	SAttrs.lpSecurityDescriptor = NULL;
	SAttrs.nLength = sizeof(SECURITY_ATTRIBUTES);

	/* make pipe, assign to stdout */
	/* we'll check for errors after CreateProcess()...*/
	res = CreatePipe(&hFromChildRead, &hFromChildWrite, &SAttrs, 0);
	res = SetStdHandle(STD_OUTPUT_HANDLE, hFromChildWrite);
	/* create non-inheritable dup of local end */
	res = DuplicateHandle(hCurProc, hFromChildRead, hCurProc, &hRead,
			0, FALSE, DUPLICATE_SAME_ACCESS);
	CloseHandle(hFromChildRead); hFromChildRead = NULL;

	res = CreatePipe(&hToChildRead, &hToChildWrite, &SAttrs, 0);
	res = SetStdHandle(STD_INPUT_HANDLE, hToChildRead);
	res = DuplicateHandle(hCurProc, hToChildWrite, hCurProc, &hWrite,
			0, FALSE, DUPLICATE_SAME_ACCESS);
	CloseHandle(hToChildWrite); hToChildWrite = NULL;

	CloseHandle(hCurProc); hCurProc = NULL;

	/* do some bookkeeping for Windows... */
	SInfo.cb = sizeof(STARTUPINFO);
	SInfo.lpReserved = NULL;
	SInfo.lpDesktop = NULL;
	SInfo.lpTitle = NULL;
	SInfo.cbReserved2 = 0;
	SInfo.lpReserved2 = NULL;
	/* don't open a console automatically, pass handles */
	SInfo.dwFlags = STARTF_USESHOWWINDOW | STARTF_USESTDHANDLES;
	SInfo.wShowWindow = SW_HIDE;
	SInfo.hStdInput = hToChildRead;
	SInfo.hStdOutput = hFromChildWrite;
	SInfo.hStdError = hStdErr; /* reuse original stderr */
	
	res = CreateProcess(NULL, /* command name */
		cmdstr, /* full command line */
		NULL,    /* default process attributes */
		NULL,    /* default security attributes */
		1,       /* inherit handles (pass doesn't work reliably) */
		Pflags,  /* process flags */
		NULL,    /* no new environment */
		NULL,    /* stay in current directory */
		&SInfo,
		&PInfo
		);
	/* reset stdin/stdout in any case */
	SetStdHandle(STD_OUTPUT_HANDLE, hStdOut);
	SetStdHandle(STD_INPUT_HANDLE, hStdIn);
	/* Oops... */
	if(res == 0) {
		char es[128];
		_snprintf(es, sizeof(es),
				"Error creating process (%d)\n", GetLastError());
		eputs(es);
		goto error;
	}
	/* close stuff we don't need */
	CloseHandle(PInfo.hThread);
	CloseHandle(hFromChildWrite); hFromChildWrite = NULL;
	CloseHandle(hToChildRead); hToChildRead = NULL;
	/* get the file descriptors */
	proc->r = _open_osfhandle((long)hRead, _O_RDONLY);
	proc->w = _open_osfhandle((long)hWrite, _O_APPEND);
	proc->pid = PInfo.dwProcessId;
	proc->running = 1;
	CloseHandle(hCurProc);
	/* Windows doesn't tell us the actual buffer size */
	return PIPE_BUF;

error: /* cleanup */
	if(PInfo.hThread) CloseHandle(PInfo.hThread);
	if(hToChildRead) CloseHandle(hToChildRead);
	if(hToChildWrite) CloseHandle(hToChildWrite);
	if(hFromChildRead) CloseHandle(hFromChildRead);
	if(hFromChildWrite) CloseHandle(hFromChildWrite);
	if(hRead) CloseHandle(hRead);
	if(hWrite) CloseHandle(hWrite);
	if(hCurProc) CloseHandle(hCurProc);
	proc->running = 0;
	return 0;
	/* There... Are we happy now? */
}


static int         /* copied size or -1 on error */
wordncopy(         /* copy (quoted) src to dest. */

char * dest,
char * src,
int dlen,
int insert_space,  /* prepend a space  */
int force_dq       /* turn 'src' into "dest" (for Win command line) */
)
{
	int slen;
	int pos = 0;

	slen = strlen(src);
	if (insert_space) {
		if (1 >= dlen) return -1;
		dest[pos++] = ' ';
	}
	if (strpbrk(src, " \f\n\r\t\v")) {
		if (force_dq && src[0] == '\'' && src[slen-1] == '\'') {
			if (slen + pos + 1 > dlen) return -1;
			dest[pos++] = '"';
			strncpy(dest + pos, src + 1, slen -2);
			pos += slen - 2;
			dest[pos++] = '"';
		} else if (src[0] == '"' && src[slen-1] == '"') {
			if (slen + pos + 1 > dlen) return -1;
			strncpy(dest + pos, src, slen);
			pos += slen;
		} else {
			if (slen + pos + 3 > dlen) return -1;
			dest[pos++] = '"';
			strncpy(dest + pos, src, slen);
			pos += slen;
			dest[pos++] = '"';
		}
	} else {
		if (slen + pos + 1 > dlen) return -1;
		strncpy(dest + pos, src, slen);
		pos += slen;
	}
	dest[pos] = '\0';
	return pos;
}	



static char *
quoted_cmdline(  /* compose command line for StartProcess() as static string */

char *cmdpath,  /* full path to executable */
char *sl[]       /* list of arguments */
)
{
	static char *cmdstr;
	static int clen;
	char *newcs;
	int newlen, pos, res, i;

	newlen = strlen(cmdpath) + 3; /* allow two quotes plus the final \0 */
	for (i = 0; sl[i] != NULL; i++) {
		newlen += strlen(sl[i]) + 3; /* allow two quotes and a space */
	}
	if (cmdstr == NULL) {
		cmdstr = (char *) malloc(newlen);
		if (cmdstr == NULL) return NULL;
	} else if (newlen > clen) {
		newcs = (char *) realloc(cmdstr, newlen);
		if (newcs == NULL) return NULL;
		cmdstr = newcs;
	}
	clen = newlen;
	pos = wordncopy(cmdstr, cmdpath, clen, 0, 1);
	if (pos < 0) return NULL;
	for (i = 0; sl[i] != NULL; i++) {
		res = wordncopy(cmdstr + pos, sl[i], clen - pos, 1, 1);
		if (res < 0) return NULL;
		pos += res;
	}
	return cmdstr;
}


int
open_process(SUBPROC *proc, char *av[])
{
	char *cmdpath;
	char *cmdstr;

	proc->running = 0;
	cmdpath = getpath(av[0], getenv("PATH"), X_OK);
	cmdstr = quoted_cmdline(cmdpath, av);
	if (cmdstr == NULL) { return 0; }
	return start_process(proc, cmdstr);
}


int
close_process(SUBPROC *proc) {
	int icres, ocres;
	DWORD pid;
	HANDLE hProc;

	ocres = close(proc->w);
	icres = close(proc->r);
	pid = proc->pid;
	if(ocres != 0 || icres != 0) {
		hProc = OpenProcess(SYNCHRONIZE|PROCESS_TERMINATE, FALSE, pid);
		/* something went wrong: enforce infanticide */
		/* other than that, it looks like we want to ignore errors here */
		if (proc->running) {
			if(hProc != NULL) {
#ifdef OBSOLETE_WINDOWS
#define KILL_TIMEOUT 10 * 1000 /* milliseconds */
				/* it might have some windows open... */
				EnumWindows((WNDENUMPROC)TerminateAppEnum, (LPARAM)pid);
				if(WaitForSingleObject(hProc, KILL_TIMEOUT)!=WAIT_OBJECT_0) {
					/* No way to avoid dangling DLLs here. */
					TerminateProcess(hProc, 0);
				}
#else
				SafeTerminateProcess(hProc, 0);
#endif
				/* WaitForSingleObject(hProc, 0); */
				/* not much use to wait on Windows */
				CloseHandle(hProc);
			}
		}
	}
	proc->running = 0;
	return 0; /* XXX we need to figure out more here... */
}


#ifdef TEST_MODULE
int
main( int argc, char **argv )
{
	SUBPROC proc;
	FILE *inf, *outf;
	int res;
	char ret[1024];
	char *command[]= {"word", "gappy word", "\"quoted words\"", "'squoted words'", NULL};

    res = open_process(&proc, command)
	if (res == 0) {
		printf("open_process() failed with return value 0\n");
		return -1;
	}
	printf("process opened with return value: %d, pid: %d,  r: %d,  w: %d\n",
			res, proc.pid, proc.r, proc.w);
	inf = fdopen(proc.r, "rb");
	outf = fdopen(proc.w, "wb");
	fprintf(outf,"0 0 0 0 1 0\n");
	fflush(outf);
	fgets(ret, sizeof(ret), inf);
	printf("%s\n",ret);
	close_process(&proc);
}
#endif

