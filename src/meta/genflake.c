#ifndef lint
static const char	RCSid[] = "$Id: genflake.c,v 1.2 2003/08/01 14:14:24 schorsch Exp $";
#endif
/*
 *  genflake.c - program to make snowflakes
 *
 *	11/12/86
 *
 *  cc genflake.c primout.o mfio.o syscalls.o misc.o
 */

#include  "meta.h"

#include  "random.h"

#define  FSIZE		6144
#define  FACT(a)	(int)((a)*5321L/FSIZE)
#define  CENTER		8192

char  *progname;

int  branch = 6;		/* number of branches */

main(argc, argv)
int  argc;
char  *argv[];
{
	int  i;
	int  seed = 0;
	int  nflakes = 1;
	
	progname = argv[0];
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 's':			/* seed */
			seed = atoi(argv[++i]);
			break;
		case 'b':			/* number of branches */
			branch = atoi(argv[++i]);
			break;
		case 'n':			/* number of flakes */
			nflakes = atoi(argv[++i]);
			break;
		default:
			fprintf(stderr, "%s: unknown option: %s\n",
					progname, argv[i]);
			exit(1);
		}
	while (seed-- > 0)
		random();
	while (nflakes-- > 0) {
		genflake(0, 0, 0, FSIZE, branch);
		pglob(PEOP, 0200, NULL);
	}
	writeof(stdout);
	exit(0);
}


genflake(x0, y0, x1, y1, nb)			/* make a flake */
int  x0, y0, x1, y1;
int  nb;
{
	int  bx0, by0, bx1, by1;
	int  i;
	
	symvect(x0, y0, x1, y1);
	for (i = 0; i < nb; i++) {
		bx0 = random()%FSIZE;
		by0 = y0 + (y1-y0)*(long)bx0/FSIZE;
		bx0 = x0 + (x1-x0)*(long)bx0/FSIZE;
		bx1 = FACT(y0-y1);
		by1 = FACT(x1-x0);
		if (random()&1) {
			bx1 += (x1-x0)/2;
			by1 += (y1-y0)/2;
		} else {
			bx1 -= (x1-x0)/2;
			by1 -= (y1-y0)/2;
		}
		bx1 = bx0 + bx1*(i+1)/(2*nb);
		by1 = by0 + by1*(i+1)/(2*nb);
		genflake(bx0, by0, bx1, by1, i);
	}
}


symvect(x0, y0, x1, y1)			/* output with symmetry */
int  x0, y0, x1, y1;
{
	symvect2(x0, y0, x1, y1);
	symvect2(x0, -y0, x1, -y1);
	symvect2(-x0, y0, -x1, y1);
	symvect2(-x0, -y0, -x1, -y1);
}


symvect2(x0, y0, x1, y1)
int  x0, y0, x1, y1;
{
	int  x0p, x1p, y0p, y1p;
	
	plseg(0, CENTER+x0, CENTER+y0, CENTER+x1, CENTER+y1);
	plseg(0, CENTER+x0/2+FACT(-y0), CENTER+y0/2+FACT(x0),
		 CENTER+x1/2+FACT(-y1), CENTER+y1/2+FACT(x1));
	plseg(0, CENTER+x0/2+FACT(y0), CENTER+y0/2+FACT(-x0),
		 CENTER+x1/2+FACT(y1), CENTER+y1/2+FACT(-x1));
}
