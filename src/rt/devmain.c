#ifndef lint
static const char	RCSid[] = "$Id: devmain.c,v 2.4 2003/02/25 02:47:22 greg Exp $";
#endif
/*
 *  devmain.c - main for independent drivers.
 *
 *	Redefine your initialization routine to dinit.
 */

#include "copyright.h"

#include "standard.h"

#include "color.h"

#include "driver.h"

struct driver	*dev = NULL;			/* output device */

FILE	*devin, *devout;			/* communications channels */

char	*progname;				/* driver name */

int	r_clear(), r_paintr(), r_getcur(), r_comout(), r_comin(), r_flush();

int	(*dev_func[NREQUESTS])() = {		/* request handlers */
		r_clear, r_paintr,
		r_getcur, r_comout,
		r_comin, r_flush
	};


main(argc, argv)		/* set up communications and main loop */
int	argc;
char	*argv[];
{
	extern struct driver	*dinit();
	int  com;
						/* set up I/O */
	progname = argv[0];
	if (argc < 3) {
		eputs("arg count\n");
		quit(1);
	}
	devin = fdopen(atoi(argv[1]), "r");
	devout = fdopen(atoi(argv[2]), "w");
	if (devin == NULL || devout == NULL || getw(devin) != COM_SENDM) {
		eputs("connection failure\n");
		quit(1);
	}
						/* open device */
	if ((dev = dinit(argv[0], argv[3])) == NULL)
		quit(1);
	putw(COM_RECVM, devout);		/* verify initialization */
	sendstate();
	fflush(devout);
						/* loop on requests */
	while ((com = getc(devin)) != EOF) {
		if (com >= NREQUESTS || dev_func[com] == NULL) {
			eputs("invalid request\n");
			quit(1);
		}
		(*dev_func[com])();		/* process request */
	}
	quit(0);			/* all done, clean up and exit */
}


void
quit(code)				/* close device and exit */
int	code;
{
	if (dev != NULL)
		(*dev->close)();
	exit(code);
}


r_clear()				/* clear screen */
{
	int	xres, yres;

	xres = getw(devin);
	yres = getw(devin);
	(*dev->clear)(xres, yres);
}


r_paintr()				/* paint a rectangle */
{
	COLOR	col;
	int	xmin, ymin, xmax, ymax;

	fread((char *)col, sizeof(COLOR), 1, devin);
	xmin = getw(devin); ymin = getw(devin);
	xmax = getw(devin); ymax = getw(devin);
	(*dev->paintr)(col, xmin, ymin, xmax, ymax);
}


r_flush()				/* flush output */
{
	if (dev->flush != NULL)
		(*dev->flush)();
	putc(COM_FLUSH, devout);
	sendstate();
	fflush(devout);
}


r_getcur()				/* get and return cursor position */
{
	int	c;
	int	x, y;
					/* get it if we can */
	if (dev->getcur == NULL) {
		c = ABORT;
		x = y = 0;
	} else
		c = (*dev->getcur)(&x, &y);
					/* reply */
	putc(COM_GETCUR, devout);
	putc(c, devout);
	putw(x, devout);
	putw(y, devout);
	fflush(devout);
}


r_comout()				/* print string to command line */
{
	char	str[256];

	mygets(str, devin);
	(*dev->comout)(str);
}


r_comin()				/* read string from command line */
{
	char	buf[256], *prompt;
					/* get prompt */
	if (getc(devin)) {
		mygets(buf, devin);
		prompt = buf;
	} else
		prompt = NULL;
					/* get string */
	(*dev->comin)(buf, prompt);
					/* reply */
	putc(COM_COMIN, devout);
	myputs(buf, devout);
	sendstate();
	fflush(devout);
}


mygets(s, fp)				/* get string from file (with nul) */
register char	*s;
register FILE	*fp;
{
	register int	c;

	while ((c = getc(fp)) != EOF)
		if ((*s++ = c) == '\0')
			return;
	*s = '\0';
}


myputs(s, fp)				/* put string to file (with nul) */
register char	*s;
register FILE	*fp;
{
	do
		putc(*s, fp);
	while (*s++);
}


void
eputs(s)				/* put string to stderr */
register char  *s;
{
	static int  midline = 0;

	if (!*s)
		return;
	if (!midline++) {
		fputs(progname, stderr);
		fputs(": ", stderr);
	}
	fputs(s, stderr);
	if (s[strlen(s)-1] == '\n') {
		fflush(stderr);
		midline = 0;
	}
}


sendstate()				/* send driver state variables */
{
	fwrite((char *)&dev->pixaspect, sizeof(dev->pixaspect), 1, devout);
	putw(dev->xsiz, devout);
	putw(dev->ysiz, devout);
	putw(dev->inpready, devout);
}
