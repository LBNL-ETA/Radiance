#ifndef lint
static const char RCSid[] = "$Id: win_usleep.c,v 3.1 2018/03/20 22:45:29 greg Exp $";
#endif
/*
 * Windows replacement usleep() function.
 */

#include "platform.h"

int
usleep(__int64 usec) 
{ 
    HANDLE timer; 
    LARGE_INTEGER ft; 

    // Convert to 100 nanosecond interval, negative value indicates relative time
    ft.QuadPart = -(10*usec); 

    timer = CreateWaitableTimer(NULL, TRUE, NULL);
    SetWaitableTimer(timer, &ft, 0, NULL, NULL, 0);
    WaitForSingleObject(timer, INFINITE);
    CloseHandle(timer);
    return(0);
}
