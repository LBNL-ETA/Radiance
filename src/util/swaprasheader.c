#ifndef lint
static const char	RCSid[] = "$Id: swaprasheader.c,v 2.2 2003/02/22 02:07:30 greg Exp $";
#endif
#include <stdio.h>
#include "rasterfile.h"

main(argc, argv)
int	argc;
char	*argv[];
{
	int	i;

	for (i = 1; i < argc; i++)
		fixfile(argv[i]);

	if (argc == 1)
		fixstream(stdin, stdout);

	exit(0);
}


fixfile(name)
char	*name;
{
	struct rasterfile	head;
	FILE	*fp;

	if ((fp = fopen(name, "r+")) == NULL) {
		fprintf(stderr, "%s: cannon open for read/write\n", name);
		exit(1);
	}
	if (fread((char *)&head, sizeof(head), 1, fp) != 1) {
		fprintf(stderr, "%s: read error\n", name);
		exit(1);
	}
	swapbytes((char *)&head, sizeof(head), sizeof(int));
	if (fseek(fp, 0L, 0) == -1) {
		fprintf(stderr, "%s: seek error\n", name);
		exit(1);
	}
	if (fwrite((char *)&head, sizeof(head), 1, fp) != 1
			|| fclose(fp) == -1) {
		fprintf(stderr, "%s: write error\n", name);
		exit(1);
	}
}


fixstream(in, out)
register FILE	*in, *out;
{
	struct rasterfile	head;
	register int	c;

	if (fread((char *)&head, sizeof(head), 1, in) != 1) {
		fprintf(stderr, "read error\n");
		exit(1);
	}
	swapbytes((char *)&head, sizeof(head), sizeof(int));
	fwrite((char *)&head, sizeof(head), 1, out);
	while ((c = getc(in)) != EOF)
		putc(c, out);
	if (ferror(out)) {
		fprintf(stderr, "write error\n");
		exit(1);
	}
}


swapbytes(top, size, width)
char	*top;
int	size, width;
{
	register char	*b, *e;
	register int	c;

	for ( ; size > 0; size -= width, top += width) {
		b = top; e = top+width-1;
		while (b < e) {
			c = *e;
			*e-- = *b;
			*b++ = c;
		}
	}
}
