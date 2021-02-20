#include <stdio.h>
main()
{
	int	c, n, x, y;
	char	com[256];
	FILE	*fp, *popen();

	while (scanf("%d", &c) == 1) {
		if (c >= '0' && c <= '9') {
	sprintf(com, "genprism alpha_mat %c - -c -l 0 0 .25 > %c.norm", c, c);
			fp = popen(com, "w");
		} else
			fp = NULL;
		n = 0;
		scanf("%d", &n);
		while (n-- > 0) {
			scanf("%d %d", &x, &y);
			if (fp != NULL)
				fprintf(fp, "%f %f\n", x/256.,y*(1.5/256.));
		}
		if (fp != NULL)
			pclose(fp);
	}
}
