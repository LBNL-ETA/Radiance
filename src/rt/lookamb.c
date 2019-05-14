#ifndef lint
static const char	RCSid[] = "$Id: lookamb.c,v 2.17 2019/05/14 17:39:10 greg Exp $";
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


static void
lookamb(			/* load & convert ambient values from a file */
	FILE  *fp
)
{
	FVECT	norm, uvec;

	while (readambval(&av, fp)) {
		decodedir(norm, av.ndir);
		decodedir(uvec, av.udir);
		if (dataonly) {
			printf("%f\t%f\t%f\t", av.pos[0], av.pos[1], av.pos[2]);
			printf("%f\t%f\t%f\t", norm[0], norm[1], norm[2]);
			printf("%f\t%f\t%f\t", uvec[0], uvec[1], uvec[2]);
			printf("%d\t%f\t%f\t%f\t", av.lvl, av.weight,
					av.rad[0], av.rad[1]);
			printf("%e\t%e\t%e\t", colval(av.val,RED),
						colval(av.val,GRN),
						colval(av.val,BLU));
			printf("%f\t%f\t", av.gpos[0], av.gpos[1]);
			printf("%f\t%f\t", av.gdir[0], av.gdir[1]);
			printf("%u\n", av.corral);
		} else {
			printf("Position:\t%f\t%f\t%f\n", av.pos[0],
					av.pos[1], av.pos[2]);
			printf("Normal:\t\t%f\t%f\t%f\n",
					norm[0], norm[1], norm[2]);
			printf("Uvector:\t%f\t%f\t%f\n",
					uvec[0], uvec[1], uvec[2]);
			printf("Lvl,Wt,UVrad:\t%d\t\t%f\t%f\t%f\n", av.lvl,
					av.weight, av.rad[0], av.rad[1]);
			printf("Value:\t\t%e\t%e\t%e\n", colval(av.val,RED),
					colval(av.val,GRN), colval(av.val,BLU));
			printf("Pos.Grad:\t%f\t%f\n", av.gpos[0], av.gpos[1]);
			printf("Dir.Grad:\t%f\t%f\n", av.gdir[0], av.gdir[1]);
			printf("Corral:\t\t%8X\n\n", av.corral);
		}
		if (ferror(stdout))
			exit(1);
	}
}


static void
writamb(			/* write binary ambient values to stdout */
	FILE  *fp
)
{
	FVECT	norm;

	for ( ; ; ) {
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, "%f %f %f",
				&av.pos[0], &av.pos[1], &av.pos[2]) != 3)
			return;
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, FVFORMAT, &norm[0], &norm[1], &norm[2]) != 3)
			return;
		av.ndir = encodedir(norm);
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, FVFORMAT, &norm[0], &norm[1], &norm[2]) != 3)
			return;
		av.udir = encodedir(norm);
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, "%hd %f %f %f", &av.lvl, &av.weight,
				&av.rad[0], &av.rad[1]) != 4)
			return;
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, "%f %f %f",
				&av.val[RED], &av.val[GRN], &av.val[BLU]) != 3)
			return;
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, "%f %f", &av.gpos[0], &av.gpos[1]) != 2)
			return;
		if (!dataonly)
			fscanf(fp, "%*s");
		if (fscanf(fp, "%f %f", &av.gdir[0], &av.gdir[1]) != 2)
			return;
		if (dataonly) {
			if (fscanf(fp, "%u", &av.corral) != 1)
				return;
		} else if (fscanf(fp, "%*s %X", &av.corral) != 1)
			return;
		av.next = NULL;
		writambval(&av, stdout);
		if (ferror(stdout))
			exit(1);
	}
}


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
