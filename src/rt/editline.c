#ifndef lint
static const char	RCSid[] = "$Id: editline.c,v 2.8 2004/03/30 16:13:01 schorsch Exp $";
#endif
/*
 *  editline.c - routine for editing raw input for rview.
 *
 *  External symbols declared in driver.h
 */

#include "copyright.h"

#include <stdio.h>
#include <string.h>

#include "color.h"

#include "driver.h"

#define iscntrl(c)	((c) < ' ')
#define isblank(c)	((c) == ' ')
#define iserase(c)	((c) == '\b' || (c) == 127)
#define iswerase(c)	((c) == 'W'-'@')
#define iskill(c)	((c) == 'U'-'@' || (c) == 'X'-'@')


extern void
editline(	/* edit input line */
	char  *buf,
	dr_getchf_t *c_get,
	dr_comoutf_t *s_put
)
{
	static char  erases[] = "\b \b";
	static char  obuf[4];
	register int  i;
	register int  c;
	
	i = 0;
	while ((c = (*c_get)()&0177) != '\n' && c != '\r')
		if (iserase(c)) {		/* single char erase */
			if (i > 0) {
				(*s_put)(erases);
				--i;
			}
		} else if (iswerase(c)) {	/* word erase */
			while (i > 0 && isblank(buf[i-1])) {
				(*s_put)(erases);
				--i;
			}
			while (i > 0 && !isblank(buf[i-1])) {
				(*s_put)(erases);
				--i;
			}
		} else if (iskill(c)) {		/* kill line */
			while (i > 0) {
				(*s_put)(erases);
				--i;
			}
		} else if (iscntrl(c)) {		/* control char */
			i = 0;
			buf[i++] = c;
			obuf[0] = '^'; obuf[1] = c|0100; obuf[2] = '\0';
			(*s_put)(obuf);
			break;
		} else {				/* regular char */
			buf[i++] = c;
			obuf[0] = c; obuf[1] = '\0';
			(*s_put)(obuf);
		}
	buf[i] = '\0';
	(*s_put)("\n");
}


static char  mybuf[512];


void
tocombuf(b, d)				/* add command(s) to my buffer */
register char  *b;
register struct driver  *d;
{
	register char  *cp;
	char  *comstart;

	for (cp = mybuf; *cp; cp++)
		;
	comstart = cp;
	while ( (*cp++ = *b) )
		if (cp >= mybuf+sizeof(mybuf)) {
			*comstart = '\0';
			return;		/* what should I do about this? */
		} else if (*b++ == '\n') {
			d->inpready++;
			comstart = cp;
		}
}


int
fromcombuf(b, d)			/* get command from my buffer */
char  *b;
struct driver  *d;
{
	register char	*cp;
						/* get next command */
	for (cp = mybuf; *cp != '\n'; cp++)
		if (!*cp)
			return(0);
	*cp++ = '\0';
#ifdef DEBUG
	(*d->comout)(mybuf);			/* echo my command */
	(*d->comout)("\n");
#endif
						/* send it as reply */
	strcpy(b, mybuf);
	d->inpready--;
						/* advance commands */
	strcpy(mybuf, cp);
	return(1);
}
