#ifndef lint
static const char	RCSid[] = "$Id: lookamb.c,v 2.7 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 *  lookamb.c - program to examine ambient components.
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

#include  "ray.h"

#include  "ambient.h"


int  dataonly = 0;
int  header = 1;
int  reverse = 0;

AMBVAL  av;


main(argc, argv)		/* load ambient values from a file */
int  argc;
char  *argv[];
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
#ifdef MSDOS
		setmode(fileno(stdout), O_BINARY);
#endif
		putambmagic(stdout);
		writamb(fp);
	} else {
#ifdef MSDOS
		setmode(fileno(fp), O_BINARY);
#endif
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


void
lookamb(fp)			/* get ambient values from a file */
FILE  *fp;
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


void
writamb(fp)			/* write binary ambient values */
FILE  *fp;
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
