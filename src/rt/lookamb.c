#ifndef lint
static const char	RCSid[] = "$Id: lookamb.c,v 2.10 2004/03/30 16:13:01 schorsch Exp $";
#endif
/*
 *  lookamb.c - program to examine ambient components.
 */

#include "copyright.h"

#include  "platform.h"
#include  "ray.h"
#include  "ambient.h"
#include  "resolu.h"


int  dataonly = 0;
int  header = 1;
int  reverse = 0;

AMBVAL  av;


int
main(		/* load ambient values from a file */
	int  argc,
	char  *argv[]
)
{
	FILE  *fp;
	int  i;

	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'd':
				dataonly = 1;
				break;
			case 'r':
				reverse = 1;
				break;
			case 'h':
				header = 0;
				break;
			default:
				fprintf(stderr, "%s: unknown option '%s'\n",
						argv[0], argv[i]);
				return(1);
			}
		else
			break;

	if (i >= argc)
		fp = stdin;
	else if ((fp = fopen(argv[i], "r")) == NULL) {
		fprintf(stderr, "%s: file not found\n", argv[i]);
		return(1);
	}
	if (reverse) {
		if (header) {
			if (checkheader(fp, "ascii", stdout) < 0)
				goto formaterr;
		} else {
			newheader("RADIANCE", stdout);
			printargs(argc, argv, stdout);
		}
		fputformat(AMBFMT, stdout);
		putchar('\n');
		SET_FILE_BINARY(stdout);
		putambmagic(stdout);
		writamb(fp);
	} else {
		SET_FILE_BINARY(fp);
		if (checkheader(fp, AMBFMT, header ? stdout : (FILE *)NULL) < 0)
			goto formaterr;
		if (!hasambmagic(fp))
			goto formaterr;
		if (header) {
			fputformat("ascii", stdout);
			putchar('\n');
		}
		lookamb(fp);
	}
	fclose(fp);
	return(0);
formaterr:
	fprintf(stderr, "%s: format error on input\n", argv[0]);
	exit(1);
}


extern void
lookamb(			/* get ambient values from a file */
	FILE  *fp
)
{
	while (readambval(&av, fp)) {
		if (dataonly) {
			printf("%f\t%f\t%f\t", av.pos[0], av.pos[1], av.pos[2]);
			printf("%f\t%f\t%f\t", av.dir[0], av.dir[1], av.dir[2]);
			printf("%d\t%f\t%f\t", av.lvl, av.weight, av.rad);
			printf("%e\t%e\t%e\t", colval(av.val,RED),
						colval(av.val,GRN),
						colval(av.val,BLU));
			printf("%f\t%f\t%f\t", av.gpos[0],
					av.gpos[1], av.gpos[2]);
			printf("%f\t%f\t%f\n", av.gdir[0],
					av.gdir[1], av.gdir[2]);
		} else {
			printf("\nPosition:\t%f\t%f\t%f\n", av.pos[0],
					av.pos[1], av.pos[2]);
			printf("Direction:\t%f\t%f\t%f\n", av.dir[0],
					av.dir[1], av.dir[2]);
			printf("Lvl,Wt,Rad:\t%d\t\t%f\t%f\n", av.lvl,
					av.weight, av.rad);
			printf("Value:\t\t%e\t%e\t%e\n", colval(av.val,RED),
					colval(av.val,GRN), colval(av.val,BLU));
			printf("Pos.Grad:\t%f\t%f\t%f\n", av.gpos[0],
					av.gpos[1], av.gpos[2]);
			printf("Dir.Grad:\t%f\t%f\t%f\n", av.gdir[0],
					av.gdir[1], av.gdir[2]);
		}
		if (ferror(stdout))
			exit(1);
	}
}


extern void
writamb(			/* write binary ambient values */
	FILE  *fp
)
{
	for ( ; ; ) {
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, "%f %f %f",
				&av.pos[0], &av.pos[1], &av.pos[2]) != 3)
			return;
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, "%f %f %f",
				&av.dir[0], &av.dir[1], &av.dir[2]) != 3)
			return;
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, "%d %f %f",
				&av.lvl, &av.weight, &av.rad) != 3)
			return;
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, "%f %f %f",
				&av.val[RED], &av.val[GRN], &av.val[BLU]) != 3)
			return;
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, "%f %f %f",
				&av.gpos[0], &av.gpos[1], &av.gpos[2]) != 3)
			return;
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, "%f %f %f",
				&av.gdir[0], &av.gdir[1], &av.gdir[2]) != 3)
			return;
		av.next = NULL;
		writambval(&av, stdout);
		if (ferror(stdout))
			exit(1);
	}
}
