#ifndef lint
static const char	RCSid[] = "$Id: setscan.c,v 2.4 2004/03/26 23:34:23 schorsch Exp $";
#endif
/*
 * Convert angle ranges of the form a-b:s,c to discrete values
 */

#include <stdlib.h>
#include <ctype.h>

#include "setscan.h"

int
setscan(			/* set up scan according to arg */
register ANGLE  *ang,
register char  *arg
)
{
	int  state = ',';
	int  start, finish, step;

	while (state) {
		switch (state) {
		case ',':
			start = atoi(arg);
			finish = start;
			step = 1;
			break;
		case '-':
			finish = atoi(arg);
			if (finish < start)
				return(-1);
			break;
		case ':':
			step = atoi(arg);
			break;
		default:
			return(-1);
		}
		if (!isdigit(*arg))
			return(-1);
		do
			arg++;
		while (isdigit(*arg));
		state = *arg++;
		if (!state || state == ',')
			while (start <= finish) {
				*ang++ = start;
				start += step;
			}
	}
	*ang = AEND;
	return(0);
}
