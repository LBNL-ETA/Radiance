#ifndef lint
static const char	RCSid[] = "$Id: colorscale.c,v 2.5 2004/03/28 20:33:13 schorsch Exp $";
#endif
/*
 *  colorscale.c - program to produce color pallets
 *
 *     9/9/87
 */

#include  <stdio.h>

#include  <math.h>

#include  "color.h"

int  primary = -1;
double  prival = 0.0;

static void colorscale(void);
static void printargs(int  ac, char  **av, FILE  *fp);


int
main(
	int  argc,
	char  *argv[]
)
{
	int  i;

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'r':
			primary = RED;
			prival = atof(argv[++i]);
			break;
		case 'g':
			primary = GRN;
			prival = atof(argv[++i]);
			break;
		case 'b':
			primary = BLU;
			prival = atof(argv[++i]);
			break;
		default:
			goto userr;
		}
	if (primary < 0 || i < argc)
		goto userr;
	printargs(argc, argv, stdout);
	printf("\n");
	printf("-Y 256 +X 256\n");
	colorscale();
	exit(0);
userr:
	fprintf(stderr, "Usage: %s -r|-g|-b value\n", argv[0]);
	exit(1);
}


static void
colorscale(void)			/* output our color scale */
{
	COLOR  scanline[256];
	int  j;
	register int  i;

	for (j = 256-1; j >= 0; j--) {
		for (i = 0; i < 256; i++)
			switch (primary) {
			case RED:
				setcolor(scanline[i],prival,i/256.0,j/256.0);
				break;
			case GRN:
				setcolor(scanline[i],i/256.0,prival,j/256.0);
				break;
			case BLU:
				setcolor(scanline[i],i/256.0,j/256.0,prival);
				break;
			}
		if (fwritescan(scanline, 256, stdout) < 0)
			goto writerr;
	}
	return;
writerr:
	fprintf(stderr, "write error in colorscale\n");
	exit(1);
}


static void
printargs(		/* print arguments to a file */
	int  ac,
	char  **av,
	FILE  *fp
)
{
	while (ac-- > 0) {
		fputs(*av++, fp);
		putc(' ', fp);
	}
	putc('\n', fp);
}
