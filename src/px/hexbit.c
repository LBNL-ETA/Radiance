#ifndef lint
static const char	RCSid[] = "$Id: hexbit.c,v 3.1 2003/02/22 02:07:27 greg Exp $";
#endif
/*
 * Convert from Radiance picture to bitmap usable by hexbit4x1.fnt
 *
 * pipe input from "pvalue -h -b -di [picture]"
 */

#include <stdio.h>

main()
{
	int	xres, yres;
	int	i, j, i1,i2,i3,i4;

	if (scanf("-Y %d +X %d\n", &yres, &xres) != 2)
		exit(1);
	for (i = 0; i < yres; i++) {
		for (j = 0; j < xres/4; j++) {
			if (scanf("%d\n%d\n%d\n%d\n", &i1, &i2, &i3, &i4) != 4)
				exit(1);
			printf("%01X", (i1<128)<<3|(i2<128)<<2|(i3<128)<<1|(i4<128));
		}
		for (j = xres%4; j--; )		/* should get it, but too lazy */
			scanf("%*d\n");
		putchar('\n');
	}
	exit(0);
}
