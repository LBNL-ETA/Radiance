#ifndef lint
static const char	RCSid[] = "$Id: linemod.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
#ifndef lint
static char sccsid[] = "@(#)linemod.c	4.1 (Berkeley) 6/27/83";
#endif

linemod(s)
char *s;
{
	char c;
	putch(033);
	switch(s[0]){
	case 'l':	
		c = 'd';
		break;
	case 'd':	
		if(s[3] != 'd')c='a';
		else c='b';
		break;
	case 's':
		if(s[5] != '\0')c='c';
		else c='`';
	}
	putch(c);
}
