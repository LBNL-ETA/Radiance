#ifndef lint
static const char	RCSid[] = "$Id: getinfo.c,v 2.12 2014/09/08 18:21:39 greg Exp $";
#endif
/*
 *  getinfo.c - program to read info. header from file.
 *
 *     1/3/86
 */

#include  <stdio.h>
#include  <string.h>

#include  "platform.h"
#include  "resolu.h"

#ifdef getc_unlocked		/* avoid nasty file-locking overhead */
#undef getchar
#undef putchar
#define getchar		getchar_unlocked
#define putchar		putchar_unlocked
#endif

#ifdef _WIN32
#include <process.h>
#define	execvp	_execvp
#endif

static gethfunc tabstr;
static void getdim(FILE *fp);
static void copycat(void);


static int
tabstr(				/* put out line followed by tab */
	char  *s,
	void *p
)
{
	while (*s) {
		putchar(*s);
		s++;
	}
	if (*--s == '\n')
		putchar('\t');
	return(0);
}


int
main(
	int  argc,
	char  **argv
)
{
	int  dim = 0;
	FILE  *fp;
	int  i;

	if (argc > 1 && !strcmp(argv[1], "-d")) {
		argc--; argv++;
		dim = 1;
		SET_FILE_BINARY(stdin);
	} else if (argc > 2 && !strcmp(argv[1], "-c")) {
		SET_FILE_BINARY(stdin);
		SET_FILE_BINARY(stdout);
		setvbuf(stdin, NULL, _IONBF, 2);
		getheader(stdin, (gethfunc *)fputs, stdout);
		printargs(argc-2, argv+2, stdout);
		fputc('\n', stdout);
		fflush(stdout);
		execvp(argv[2], argv+2);
		perror(argv[2]);
		return 1;
	} else if (argc == 2 && !strcmp(argv[1], "-")) {
		SET_FILE_BINARY(stdin);
		SET_FILE_BINARY(stdout);
		getheader(stdin, NULL, NULL);
		copycat();
		return 0;
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
				tabstr(":\n", NULL);
				getheader(fp, tabstr, NULL);
				fputc('\n', stdout);
			}
			fclose(fp);
		}
	}
	if (argc == 1) {
		if (dim) {
			getdim(stdin);
		} else {
			getheader(stdin, (gethfunc *)fputs, stdout);
			fputc('\n', stdout);
		}
	}
	return 0;
}


static void
getdim(				/* get dimensions from file */
	FILE  *fp
)
{
	int  j;
	int  c;

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
			if (c == 0) {
				if (++j >= 4)
					break;
				putchar(' ');
			} else {
				putchar(c);
			}
		putchar('\n');
		break;
	default:		/* ??? */
		fputs("unknown file type\n", stdout);
		break;
	}
}


static void
copycat(void)			/* copy input to output */
{
	char	buf[8192];
	int	n;

	fflush(stdout);
	while ((n = fread(buf, 1, sizeof(buf), stdin)) > 0)
		if (write(fileno(stdout), buf, n) != n)
			break;
}
