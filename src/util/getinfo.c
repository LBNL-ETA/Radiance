#ifndef lint
static const char	RCSid[] = "$Id: getinfo.c,v 2.21 2019/07/24 17:27:54 greg Exp $";
#endif
/*
 *  getinfo.c - program to read info. header from file.
 *
 *     1/3/86
 */

#include  "rtio.h"
#include  "platform.h"
#include  "rtprocess.h"
#include  "resolu.h"

#ifdef getc_unlocked		/* avoid nasty file-locking overhead */
#undef getc
#undef getchar
#undef putchar
#define getc		getc_unlocked
#define getchar		getchar_unlocked
#define putchar		putchar_unlocked
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

	if (argc > 1 && (argv[1][0] == '-') | (argv[1][0] == '+') &&
			argv[1][1] == 'd') {
		dim = 1 - 2*(argv[1][0] == '-');
		argc--; argv++;
	}
#ifdef getc_unlocked				/* avoid lock/unlock overhead */
	flockfile(stdin);
#endif
	SET_FILE_BINARY(stdin);
	if (argc > 2 && !strcmp(argv[1], "-c")) {
		SET_FILE_BINARY(stdout);
		setvbuf(stdin, NULL, _IONBF, 2);
		if (getheader(stdin, (gethfunc *)fputs, stdout) < 0) {
			fputs("Bad header!\n", stderr);
			return 1;
		}
		printargs(argc-2, argv+2, stdout);
		fputc('\n', stdout);
		if (dim) {			/* copy resolution string? */
			RESOLU	rs;
			if (!fgetsresolu(&rs, stdin)) {
				fputs("No resolution string!\n", stderr);
				return 1;
			}
			if (dim > 0)
				fputsresolu(&rs, stdout);
		}
		fflush(stdout);
		execvp(argv[2], argv+2);
		perror(argv[2]);
		return 1;
	} else if (argc > 2 && !strcmp(argv[1], "-a")) {
		SET_FILE_BINARY(stdout);
		if (getheader(stdin, (gethfunc *)fputs, stdout) < 0) {
			fputs("Bad header!\n", stderr);
			return 1;
		}
		for (i = 2; i < argc; i++) {
			int	len = strlen(argv[i]);
			if (!len) continue;
			fputs(argv[i], stdout);
			if (argv[i][len-1] != '\n')
				fputc('\n', stdout);
		}
		fputc('\n', stdout);
		copycat();
		return 0;
	} else if (argc == 2 && !strcmp(argv[1], "-")) {
		SET_FILE_BINARY(stdout);
		if (getheader(stdin, NULL, NULL) < 0) {
			fputs("Bad header!\n", stderr);
			return 1;
		}
		if (dim < 0) {			/* skip resolution string? */
			RESOLU	rs;
			if (!fgetsresolu(&rs, stdin)) {
				fputs("No resolution string!\n", stderr);
				return 1;
			}
		}
		copycat();
		return 0;
	}
	for (i = 1; i < argc; i++) {
		fputs(argv[i], stdout);
		if ((fp = fopen(argv[i], "r")) == NULL)
			fputs(": cannot open\n", stdout);
		else {
			if (dim < 0) {			/* dimensions only */
				if (getheader(fp, NULL, NULL) < 0) {
					fputs("bad header!\n", stdout);
					continue;	
				}
				fputs(": ", stdout);
				getdim(fp);
			} else {
				tabstr(":\n", NULL);
				if (getheader(fp, tabstr, NULL) < 0) {
					fputs(argv[i], stderr);
					fputs(": bad header!\n", stderr);
					return 1;
				}
				fputc('\n', stdout);
				if (dim > 0) {
					fputc('\t', stdout);
					getdim(fp);
				}
			}
			fclose(fp);
		}
	}
	if (argc == 1) {
		if (dim < 0) {
			if (getheader(stdin, NULL, NULL) < 0) {
				fputs("Bad header!\n", stderr);
				return 1;	
			}
			getdim(stdin);
		} else {
			if (getheader(stdin, (gethfunc *)fputs, stdout) < 0) {
				fputs("Bad header!\n", stderr);
				return 1;
			}
			fputc('\n', stdout);
			if (dim > 0)
				getdim(stdin);
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
		if (writebuf(fileno(stdout), buf, n) != n)
			break;
}
