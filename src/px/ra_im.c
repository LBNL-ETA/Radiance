#ifndef lint
static const char	RCSid[] = "$Id: ra_im.c,v 2.2 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 *  ra_im.c - convert Radiance picture to imagetools raw format.
 *
 *	9/16/88
 */

#include <stdio.h>


#define PCOMM		"pvalue -h -b -db"

#define MINVAL		1
#define MAXVAL		252

extern FILE	*popen(), *freopen();


main(argc, argv)
int	argc;
char	*argv[];
{
	register int	c;
	register FILE	*fp;

	if (argc > 3) {
		fputs("Usage: ", stderr);
		fputs(argv[0], stderr);
		fputs(" [infile [outfile]]\n", stderr);
		exit(1);
	}
	if (argc > 1 && freopen(argv[1], "r", stdin) == NULL) {
		perror(argv[1]);
		exit(1);
	}
	if (argc > 2 && freopen(argv[2], "w", stdout) == NULL) {
		perror(argv[2]);
		exit(1);
	}
	if ((fp = popen(PCOMM, "r")) == NULL) {
		perror(argv[0]);
		exit(1);
	}
	while ((c = getc(fp)) != EOF) {
		if (c < MINVAL)
			putc(MINVAL, stdout);
		else if (c > MAXVAL)
			putc(MAXVAL, stdout);
		else
			putc(c, stdout);
	}
	exit(pclose(fp));
}
