#ifndef lint
static const char	RCSid[] = "$Id: getinfo.c,v 2.25 2022/04/13 15:43:06 greg Exp $";
#endif
/*
 *  getinfo.c - program to read info. header from file.
 *
 *     1/3/86
 */

#include  <ctype.h>
#include  "rtio.h"
#include  "platform.h"
#include  "rtprocess.h"
#include  "resolu.h"

static char	fmt[MAXFMTLEN] = "*";


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
		fputc('\t', stdout);
	return(0);
}


static int
adjheadline(			/* check for lines to remove */
	char  *s,
	void *p
)
{
	char	**av;

	if (formatval(fmt, s))
		return(0);	/* don't echo format */

				/* check for match to skip */
	for (av = (char **)p; *av; av++) {
		char	*s1 = s;
		char	*s2 = *av;
		if (isspace(*s2)) {
			if (!isspace(*s1))
				continue;
			while (isspace(*++s1)) ;
			while (isspace(*++s2)) ;
		}
		while (*s2 && *s1 == *s2) {
			if ((*s1 == '=') & (s1 > s))
				return(0);
			if (isspace(*s2))
				break;
			s1++; s2++;
		}
		if (s1 == s)
			continue;
		if ((*s1 == '=') & !*s2)
			return(0);
		if (isspace(*s1) & isspace(*s2))
			return(0);
	}
	fputs(s, stdout);	/* copy if no match */
	return(0);
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
				fputc(' ', stdout);
			} else {
				putchar(c);
			}
		fputc('\n', stdout);
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
	flockfile(stdout);
#endif
	SET_FILE_BINARY(stdin);
	if (argc > 2 && !strcmp(argv[1], "-c")) {
		SET_FILE_BINARY(stdout);
		setvbuf(stdin, NULL, _IONBF, 2);
		if (checkheader(stdin, fmt, stdout) < 0) {
			fputs("Bad header!\n", stderr);
			return 1;
		}
		printargs(argc-2, argv+2, stdout);
		if (fmt[0] != '*')		/* better be the same! */
			fputformat(fmt, stdout);
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
	}
	if (argc > 2 && !strcmp(argv[1], "-a")) {
		SET_FILE_BINARY(stdout);
		if (checkheader(stdin, fmt, stdout) < 0) {
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
		if (fmt[0] != '*')
			fputformat(fmt, stdout);
		fputc('\n', stdout);
		copycat();
		return 0;
	}
	if (argc > 2 && !strcmp(argv[1], "-r")) {
		SET_FILE_BINARY(stdout);
		if (getheader(stdin, adjheadline, argv+2) < 0) {
			fputs("Bad header!\n", stderr);
			return 1;
		}
		for (i = 2; i < argc; i++) {	/* add lines w/ var[= ]value */
			int	len = strlen(argv[i]);
			int	strt = 0, j;
			while(isspace(argv[i][strt])) strt++;
			if (strt == len) continue;
			while (isspace(argv[i][len-1])) len--;
			for (j = strt+1; j < len-1; j++)
				if (argv[i][j] == '=' || isspace(argv[i][j])) {
					for (j = 0; j < len; j++)
						putchar(argv[i][j]);
					fputc('\n', stdout);
					break;
				}
		}
		if (fmt[0] != '*')
			fputformat(fmt, stdout);
		fputc('\n', stdout);
		copycat();
		return 0;
	}
	if (argc == 2 && !strcmp(argv[1], "-")) {
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
#ifdef getc_unlocked				/* avoid lock/unlock overhead */
			flockfile(fp);
#endif
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
