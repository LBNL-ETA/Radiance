#ifndef lint
static const char	RCSid[] = "$Id: tty.c,v 2.4 2004/03/30 16:13:01 schorsch Exp $";
#endif
/*
 *  tty.c - i/o for terminal drivers.
 */

#include "copyright.h"

#include  <stdio.h>
#include  <signal.h>
#include  <fcntl.h>
#include  <sys/ioctl.h>

#include  "driver.h"

/*
 *  Drivers using these routines must call getch()
 *  for input from stdin.
 */

struct sgttyb  ttymode;			/* original tty modes */
int  ttyfd;				/* tty file descriptor */
struct driver  *ttydev;			/* tty device */

static dr_cominf_t ttyin;
static dr_getchf_t;
static void ttydone(void);
static void newinp(void);


extern void /* XXX declare */
ttyset(				/* set up raw tty device */
	struct driver  *dev,
	int  fd
)
{
	struct sgttyb  flags;

	if (!isatty(fd))
		return(-1);

	ttyfd = fd;
	ioctl(ttyfd, TIOCGETP, &ttymode);
	ioctl(ttyfd, TIOCGETP, &flags);
	flags.sg_flags |= RAW;			/* also turns output */
	flags.sg_flags &= ~ECHO;		/* processing off */
	ioctl(ttyfd, TIOCSETP, &flags);
						/* install input handler */
	ttydev = dev;
	ttydev->comin = ttyin;
	ttydev->inpready = 0;
	signal(SIGIO, newinp);
	fcntl(ttyfd, F_SETFL, FASYNC);
	return(0);
}


static void
ttydone(void)				/* restore tty modes */
{
	fcntl(ttyfd, F_SETFL, 0);
	signal(SIGIO, SIG_DFL);
	ioctl(ttyfd, TIOCSETP, &ttymode);
	ttydev->comin = NULL;
	ttydev->inpready = 0;
	ttydev = NULL;
}


static int
getch(void)					/* get a character in raw mode */
{
	register int  c;

	fcntl(ttyfd, F_SETFL, 0);		/* allow read to block */
	if (ttydev->inpready > 0)
		ttydev->inpready--;
	c = getchar();
	fcntl(ttyfd, F_SETFL, FASYNC);
	return(c);
}


static void
ttyin(			/* read a line in raw mode */
	char  *buf,
	char  *prompt
)
{

	if (prompt != NULL)
		(*ttydev->comout)(prompt);
	editline(buf, getch, ttydev->comout);
}


static void
newinp(void)				/* new input */
{
	ttydev->inpready++;
}
