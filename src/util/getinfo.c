#ifndef lint
static const char	RCSid[] = "$Id: getinfo.c,v 2.5 2003/02/22 02:07:30 greg Exp $";
#endif
/*
 *  getinfo.c - program to read info. header from file.
 *
 *     1/3/86
 */

#include  <stdio.h>

#ifdef MSDOS
#include <fcntl.h>
extern int  _fmode;
#endif

extern int  fputs();


int
tabstr(s)				/* put out line followed by tab */
register char  *s;
{
	while (*s) {
		putchar(*s);
		s++;
	}
	if (*--s == '\n')
		putchar('\t');
	return(0);
}


main(argc, argv)
int  argc;
char  **argv;
{
	int  dim = 0;
	FILE  *fp;
	int  i;

	if (argc > 1 && !strcmp(argv[1], "-d")) {
		argc--; argv++;
		dim = 1;
#ifdef MSDOS
		setmode(fileno(stdin), _fmode = O_BINARY);
#endif
	} else if (argc == 2 && !strcmp(argv[1], "-")) {
#ifdef MSDOS
		setmode(fileno(stdin), O_BINARY);
		setmode(fileno(stdout), O_BINARY);
#endif
		getheader(stdin, NULL, NULL);
		copycat();
		exit(0);
	}
	for (i = 1; i < argc; i++) {
		fputs(argv[i], stdout);
		if ((fp = fopen(argv[i], "r")) == NULL)
			fputs(": cannot open\n", stdout);
		else {
			if (dim) {
				fputs(": ", stdout);
				getdim(fp);
			} else {
				tabstr(":\n");
				getheader(fp, tabstr, NULL);
				putchar('\n');
			}
			fclose(fp);
		}
	}
	if (argc == 1)
		if (dim) {
			getdim(stdin);
		} else {
			getheader(stdin, fputs, stdout);
			putchar('\n');
		}
	exit(0);
}


getdim(fp)				/* get dimensions from file */
register FILE  *fp;
{
	int  j;
	register int  c;

	getheader(fp, NULL, NULL);	/* skip header */

	switch (c = getc(fp)) {
	case '+':		/* picture */
	case '-':
		do
			putchar(c);
		while (c != '\n' && (c = getc(fp)) != EOF);
		break;
	case 1:			/* octree */
		getc(fp);
		j = 0;
		while ((c = getc(fp)) != EOF)
			if (c == 0)
				if (++j >= 4)
					break;
				else
					putchar(' ');
			else
				putchar(c);
		putchar('\n');
		break;
	default:		/* ??? */
		fputs("unknown file type\n", stdout);
		break;
	}
}


copycat()			/* copy input to output */
{
	register int	c;

	while ((c = getchar()) != EOF)
		putchar(c);
}
