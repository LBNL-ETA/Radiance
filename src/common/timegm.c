/*
 *  timegm.c
 *
 *  Replacement for missing GNU library function.
 *
 */


#include <time.h>
#include <stdlib.h>

#ifdef _WIN32
char *
setGMT()
{
	static time_t	prevTZ;
	prevTZ = _timezone;
	_timezone = 0;
	return (char *)&prevTZ;
}
void
resetTZ(char *cp)
{
	_timezone = *(time_t *)cp;
}
#else
char *
setGMT()
{
	char  *tz = getenv("TZ");
	setenv("TZ", "", 1);
	tzset();
	return tz;
}
void
resetTZ(char *tz)
{
	if (tz)
		setenv("TZ", tz, 1);
	else
		unsetenv("TZ");
	tzset();
}
#endif

/* Convert GMT to UTC seconds from the epoch */
time_t
timegm(struct tm *tm)
{
	time_t  ret;
	char  *tz;

	tz = setGMT();
	ret = mktime(tm);
	resetTZ(tz);
	return(ret);
}
