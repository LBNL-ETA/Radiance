#ifndef lint
static const char RCSid[] = "$Id: dmessage.c,v 2.2 2025/06/01 03:25:17 greg Exp $";
#endif
/*
 *  dmessage.c
 *  panlib
 *
 *  Debug message reporting routines.
 *
 *  Created by gward on Thu May 10 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "dmessage.h"

#define DM_NCLASSES	DMCnerr		/* number of message classes */
#define DM_MOSTFLAGS	(DMFrecord|DMFlog|DMFstderr|DMFalert)

char			dmessage_buf[DM_BUF_LEN];
const char *		dmessage_class_name[DM_NCLASSES+1] = {
				"assert",
				"memory",
				"system",
				"parameter",
				"resource",
				"data",
				"input",
				"warning",
				"info",
				"trace",
				"none"
			};
int			dmessage_class_flags[DM_NCLASSES+1] = {
				DM_MOSTFLAGS|DMFabort,
				DM_MOSTFLAGS|DMFexit,
				DM_MOSTFLAGS,
				DM_MOSTFLAGS,
				DM_MOSTFLAGS,
				DM_MOSTFLAGS,
				DM_MOSTFLAGS,
				DMFrecord|DMFlog|DMFstderr,
				DMFrecord|DMFlog,
				DMFrecord,
				0
			};
int			dmessage_class_done[DM_NCLASSES];
const char *		dmessage_file[DM_NCLASSES];
int			dmessage_line[DM_NCLASSES];
const char *		dmessage_record[DM_NCLASSES+1];
DMsgClass		dmessage_last_class = DMCnerr;
FILE *			dmessage_logfp = NULL;
int			(*dmessage_call)(DMsgClass cls, const char *msg,
					const char *file, int line) = NULL;
int			(*dmessage_alert)(const char *msg) = NULL;
void			(*dmessage_exit)(int status) = exit;
void			(*dmessage_abort)() = abort;

#define	DM_MAX_RECORD	128

static char		dmessage_mine[DM_NCLASSES][DM_MAX_RECORD];

/* Report/record message, aborting program if cls is DMCassert */
int
dmessage(DMsgClass cls, const char *msg, const char *file, int line)
{
	int		todo;
        const char      *cp;
	if ((cls < 0) | (cls >= DMCnerr))
		return 0;
	todo = dmessage_class_flags[cls];
	if (msg == NULL || !*msg)
		msg = "(missing message)";
	if (file != NULL && ((cp = strrchr(file, '/')) != NULL ||
				(cp = strrchr(file, '\\')) != NULL))
		file = cp+1;				/* use tail */
	dmessage_file[cls] = file;			/* remember call */
	dmessage_line[cls] = line;
	dmessage_class_done[cls] = todo;
	if (cls < dmessage_last_class)			/* mark class */
		dmessage_last_class = cls;
	if (dmessage_call != NULL)			/* notify callback */
		todo &= ~(*dmessage_call)(cls, msg, file, line);
	if (todo & DMFrecord) {				/* save message */
		if (msg != dmessage_buf && strlen(msg) >= DM_MAX_RECORD) {
			dmessage_record[cls] = msg;
		} else {
			strncpy(dmessage_mine[cls], msg, DM_MAX_RECORD-1);
			dmessage_record[cls] = dmessage_mine[cls];
		}
		todo &= ~DMFrecord;
	}
	if (todo & DMFstderr) {				/* send to stderr */
		if (file != NULL)
			fprintf(stderr, "%s@%d>%s: %s", file, line,
					dmessage_class_name[cls], msg);
		else
			fprintf(stderr, "%s: %s",
					dmessage_class_name[cls], msg);
                if ((cls == DMCsystem) | (cls == DMCresource) && errno) {
			fputs(": ", stderr);
			perror(NULL);
		} else
			fputc('\n', stderr);
		if (fflush(stderr) != EOF) {
			todo &= ~DMFstderr;
			if (dmessage_logfp == stderr)
				todo &= ~DMFlog;
		}
	}
	if (todo & DMFlog && dmessage_logfp != NULL) {	/* log to file */
		if (file != NULL)
			fprintf(dmessage_logfp, "%s@%d>%s: %s\n", file, line,
					dmessage_class_name[cls], msg);
		else
			fprintf(dmessage_logfp, "%s: %s\n",
					dmessage_class_name[cls], msg);
		if (fflush(dmessage_logfp) != EOF)
			todo &= ~DMFlog;
	}
	if (todo & DMFalert && dmessage_alert != NULL)	/* alert user */
		if ((*dmessage_alert)(msg))
			todo &= ~DMFalert;

	if (todo & DMFabort && dmessage_abort != NULL) {/* abort */
		(*dmessage_abort)();
		todo &= ~DMFabort;
	}
	if (todo & DMFexit && dmessage_exit != NULL) {	/* exit */
		(*dmessage_exit)(1);
		todo &= ~DMFexit;
	}
							/* say what was done */
	return dmessage_class_done[cls] = (dmessage_class_flags[cls] & ~todo);
}
