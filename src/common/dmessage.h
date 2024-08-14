/*
 *  dmessage.h
 *  panlib
 *
 *  Depends on <stdio.h>
 *
 *  Debug and error message logging and recovery.
 *
 *  Created by gward on Thu May 10 2001.
 *  Modified by plonghurst on Thursday Feb 15 2007.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#ifndef _DMESSAGE_H_
#define _DMESSAGE_H_

#ifdef __cplusplus
extern "C" {
#endif

/********************
 * Type Definitions *
 ********************/

/***********************************************************************
 * Debug message classes
 *
 *  If you are unsure what class to report, guess on the high side, which
 *  is less serious.  For example, report an error like DMCdata even
 *  when you think it's DMCsystem, unless you really know.
 *  That way, more serious errors will be recorded from higher in the call
 *  tree.  The vast majority of errors fall into the DMCdata class.
 */
typedef enum {
		DMCassert,		/* assertion failure */
		DMCmemory,		/* out of memory */
		DMCsystem,		/* system error */
		DMCparameter,		/* program parameter error */
		DMCresource,		/* file/resource unavailable */
		DMCdata,		/* data error */
		DMCinput,		/* input error */
		DMCwarning,		/* warning message */
		DMCinfo,		/* informative output */
		DMCtrace,		/* trace point */
		DMCnerr			/* no error (terminator) */
} DMsgClass;

/***********************************************************************
 * Message flags:  record, write to logfile, send to stderr,
 *			alert user, exit program, abort
 */
typedef enum {
		DMFrecord=01, DMFlog=02, DMFstderr=04,
		DMFalert=010, DMFexit=020, DMFabort=040
} DMsgFlag;

/************************
 * Functions and Macros *
 ************************/

/***********************************************************************
 * int
 * dmessage(int cls, const char *msg, const char *file, int line);
 *
 *  Basic call interface to debug message handler.  The file name
 *  and line number are generally taken from the ANSI-C __FILE__
 *  and __LINE__ predefined macros.  The message class are taken
 *  from the list above, and the message should not contain newlines.
 *  All DMsgClass classes except DMCassert are recoverable and
 *  dmessage() will return.  The return value is an or'ing of
 *  DMsgFlags according to what was done.
 */
#ifdef __cplusplus
extern int	dmessage(DMsgClass cls, const char *msg,
			const char *file = NULL, int line = 0);
#else
extern int	dmessage(DMsgClass cls, const char *msg,
			const char *file, int line);
#endif

/***********************************************************************
 * int
 * DMESG(int cls, const char *msg)
 *
 *  Basic macro call for message reporting.  What actually happens to
 *  the message depends on the current settings, controlled by
 *  global variables described in the next section.
 */
#define DMESG(cls, msg)		dmessage(cls, msg, __FILE__, __LINE__)

/***********************************************************************
 * int
 * DTEST(int cnd, int cls, const char *msg)
 *
 *  Conditional call to DMESG() for convenience.  Zero is returned if
 *  the condition (cnd) evaluates to zero, and a collection of DMsgFlags
 *  saying what was done if the condition evaluates to non-zero.
 */
#define DTEST(cnd, cls, msg)	((cnd) ? DMESG(cls, msg) : 0)

/***********************************************************************
 * int
 * DMESGF(int cls, const char *fmt, va_list..)
 *
 *  Formatting version of DMESG() takes printf(3) format string
 *  and a variable argument list to format the message.
 */
#define DMESGF(cls, fmt, val)	(sprintf(dmessage_buf, fmt, val), \
					DMESG(cls, dmessage_buf))

/***********************************************************************
 * int
 * DTESTF(int cnd, DMsgClass cls, const char *fmt, va_list..)
 *
 *  Formatting version of DTEST() takes printf(3) format string
 *  and a variable argument list to format the message, but only
 *  if the condition (cnd) evaluates to non-zero.
 */
#define DTESTF(cnd, cls, fmt, val) \
				((cnd) ? DMESGF(cls, fmt, val) : 0)

/***********************************************************************
 * int
 * DRESET()
 *
 *  Reset last error message class.  Call before potential error-
 *  producing subroutine, then call DMESC() or DMESCF() afterwards.
 */
#define DRESET()	(dmessage_last_class = DMCnerr)

/***********************************************************************
 * int
 * DMESC(int cls, const char *msg)
 *
 *  Report only if this error is a higher priority than any since
 *  DRESET() was last called.
 */
#define DMESC(cls, msg)		DTEST((cls) < dmessage_last_class, \
					cls, msg)
 
/***********************************************************************
 * int
 * DMESCF(int cls, const char *fmt, va_list..)
 *
 *  Formatting version of DMESC() takes printf(3) format string, but
 *  only reports if this error is higher priority than any since
 *  DRESET() was last called.
 */
#define DMESCF(cls, fmt, val)	DTESTF((cls) < dmessage_last_class, \
					cls, fmt, val)

/***********************************************************************
 * void
 * DASSERT(int cnd)
 *
 *  Replacement for assert() macro; calls DTEST() with DMCassert class
 *  if asc yeilds 0.  If NDEBUG is defined, the macro becomes a no-op.
 */
#ifdef NDEBUG
#define DASSERT(asc)		((void)0)
#else
#if defined(__STDC__) || defined(__cplusplus)
#define DASSERT(asc)		((void)DTESTF(!(asc), DMCassert, \
					"assertion failed ! (%s)", #asc))
#else
#define DASSERT(asc)		((void)DTESTF(!(asc), DMCassert, \
					"assertion failed ! (%s)", "asc"))
#endif
#endif

/*********************
 * Global Variabless *
 *********************/

/***********************************************************************
 * char			dmessage_buf[DM_BUF_LEN]
 *
 *  Buffer to hold formatted message contents temporarily.
 */
#define DM_BUF_LEN	1024

extern char		dmessage_buf[DM_BUF_LEN];

/***********************************************************************
 * const char *		dmessage_class_name[]
 *
 *  The name of each of class in DMsgClass.
 */
extern const char *	dmessage_class_name[];

/***********************************************************************
 * int			dmessage_class_flags[]
 *
 *  Flags indicating what to do with each message class.  By default,
 *  the most recent message will be stored in the dmessage_last[] array
 *  (DMFrecord), but if other flags are set, the message will be sent
 *  other places as well.  These flags are available to be altered by the
 *  controlling application, and will not be changed from their initial
 *  defaults (DMFrecord at minimum) by any of the dmessage routines.
 *  Beware:  if the DMrecord flag is turned off for a given class and
 *  dmessage() fails to do anything, it will return zero and
 *  the DTEST() and DTESTF() macro calls may return zero even
 *  though their condition evaluated to non-zero.  This problem is
 *  avoided by keeping DMrecord set for all classes.
 */
extern int		dmessage_class_flags[];

/***********************************************************************
 * int			dmessage_class_done[]
 *
 *  Flags indicating what was last done with each message class.
 *  In other words, the last result of dmessage() for each class.
 */
extern int		dmessage_class_done[];

/***********************************************************************
 * const char *		dmessage_file[]
 *
 *  Most recent file reported for each error class.
 */
extern const char *	dmessage_file[];

/***********************************************************************
 * int			dmessage_line[]
 *
 *  Most recent line number reported for each error class.
 */
extern int		dmessage_line[];

/***********************************************************************
 * const char *		dmessage_record[]
 *
 *  Most recent message for each class with the DMFrecord flag set,
 *  exactly as passed to dmessage().  If no message of
 *  a given class has been reported, then the corresponding
 *  dmessage_record pointer will be NULL.
 */
extern const char *	dmessage_record[];

/***********************************************************************
 * int			dmessage_last_class
 *
 *  The lowest class message reported so far, set to DMCnerr initially.
 *  The DMESG_LAST macro provides convenient access to most recent error.
 */
extern DMsgClass	dmessage_last_class;

#define DMESG_LAST	(dmessage_record[dmessage_last_class])

/***********************************************************************
 * FILE *		dmessage_logfp
 *
 *  Pointer to open log file.  Each message is preceeded by the file
 *  name, line number, and class for log file and stderr output.  Data
 *  is flushed after each message.  If this pointer is NULL (the default),
 *  then logging is disabled even for classes with the DMFlog flag set.
 */
extern FILE *		dmessage_logfp;

/***********************************************************************
 * int			(*dmessage_call)(DMsgClass cls, const char *msg,
 *					const char *file, int line)
 *
 *  Pointer to function called by dmessage() to share error handling.
 *  If assigned, this function is called with the same arguments passed
 *  to dmessage(), independent of the class flag settings.  The returned
 *  value is then used as a mask on the class flags to turn off some or
 *  all of the default actions.  The flagged actions are considered
 *  "done" and set accordingly in the dmessage_class_done global.
 *  Defaults to NULL.
 */
extern int		(*dmessage_call)(DMsgClass cls, const char *msg,
					const char *file, int line);

/***********************************************************************
 * int			(*dmessage_alert)(const char *msg)
 *
 *  Pointer to function called by dmessage() to alert user.
 *  The function should return non-zero if the user was notified,
 *  or zero if the message could not be displayed for some reason.
 *  Defaults to NULL, which results in no DMFalert action.
 */
extern int		(*dmessage_alert)(const char *msg);

/***********************************************************************
 * void			(*dmessage_exit)(int status)
 *
 *  Pointer to function called by dmessage() for classes with the
 *  DMFexit flag set, which is usually reserved for the DMCmemory class.
 *  Defaults to system exit() function, but may be reassigned by the
 *  controlling application.  Assigning a value of NULL causes dmessage()
 *  to return from these calls, which might result in a memory fault if
 *  your program does not check for NULL pointer values.  In general,
 *  library authors should NOT assume the DMCmemory class calls exit().
 */
extern void		(*dmessage_exit)(int status);

/***********************************************************************
 * void			(*dmessage_abort)()
 *
 *  Pointer to function called by dmessage() for classes with the
 *  DMFabort flag set, which is usually reserved for the DMCassert class.
 *  Defaults to system abort() function, but may be reassigned by the
 *  controlling application.  Assigning a value of NULL causes dmessage()
 *  to return from these calls, which will probably prove disasterous.
 *  A better idea if the program wants to continue is to use setjmp() and
 *  longjmp() to recover control at a lower point in the call tree.  If
 *  assigned, this call should never return.
 */
extern void		(*dmessage_abort)();

#ifdef __cplusplus
}
#endif

#endif /* ! _DMESSAGE_H_ */
