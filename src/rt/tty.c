#ifndef lint
static const char	RCSid[] = "$Id: tty.c,v 2.3 2003/02/25 02:47:23 greg Exp $";
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


ttyset(dev, fd)				/* set up raw tty device */
struct driver  *dev;
int  fd;
{
	int  ttyin(), newinp();
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


ttydone()				/* restore tty modes */
{
	fcntl(ttyfd, F_SETFL, 0);
	signal(SIGIO, SIG_DFL);
	ioctl(ttyfd, TIOCSETP, &ttymode);
	ttydev->comin = NULL;
	ttydev->inpready = 0;
	ttydev = NULL;
}


int
getch()					/* get a character in raw mode */
{
	register int  c;

	fcntl(ttyfd, F_SETFL, 0);		/* allow read to block */
	if (ttydev->inpready > 0)
		ttydev->inpready--;
	c = getchar();
	fcntl(ttyfd, F_SETFL, FASYNC);
	return(c);
}


static
ttyin(buf, prompt)			/* read a line in raw mode */
char  *buf, *prompt;
{
	int  getch();

	if (prompt != NULL)
		(*ttydev->comout)(prompt);
	editline(buf, getch, ttydev->comout);
}


static
newinp()				/* new input */
{
	ttydev->inpready++;
}
