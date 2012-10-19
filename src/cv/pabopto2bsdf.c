#ifndef lint
static const char RCSid[] = "$Id: pabopto2bsdf.c,v 2.1 2012/10/19 04:14:29 greg Exp $";
#endif
/*
 * Load measured BSDF data in PAB-Opto format.
 *
 *	G. Ward
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "platform.h"
#include "bsdfrep.h"
				/* global argv[0] */
char			*progname;

/* Load a set of measurements corresponding to a particular incident angle */
static int
load_pabopto_meas(const char *fname)
{
	FILE	*fp = fopen(fname, "r");
	int	inp_is_DSF = -1;
	double	new_theta, new_phi, theta_out, phi_out, val;
	char	buf[2048];
	int	n, c;
	
	if (fp == NULL) {
		fputs(fname, stderr);
		fputs(": cannot open\n", stderr);
		return(0);
	}
#ifdef DEBUG
	fprintf(stderr, "Loading measurement file '%s'...\n", fname);
#endif
				/* read header information */
	while ((c = getc(fp)) == '#' || c == EOF) {
		if (fgets(buf, sizeof(buf), fp) == NULL) {
			fputs(fname, stderr);
			fputs(": unexpected EOF\n", stderr);
			fclose(fp);
			return(0);
		}
		if (!strcmp(buf, "format: theta phi DSF\n")) {
			inp_is_DSF = 1;
			continue;
		}
		if (!strcmp(buf, "format: theta phi BSDF\n")) {
			inp_is_DSF = 0;
			continue;
		}
		if (sscanf(buf, "intheta %lf", &new_theta) == 1)
			continue;
		if (sscanf(buf, "inphi %lf", &new_phi) == 1)
			continue;
		if (sscanf(buf, "incident_angle %lf %lf",
				&new_theta, &new_phi) == 2)
			continue;
	}
	ungetc(c, fp);
	if (inp_is_DSF < 0) {
		fputs(fname, stderr);
		fputs(": unknown format\n", stderr);
		fclose(fp);
		return(0);
	}
					/* prepare input grid */
	new_bsdf_data(new_theta, new_phi);
					/* read actual data */
	while (fscanf(fp, "%lf %lf %lf\n", &theta_out, &phi_out, &val) == 3)
		add_bsdf_data(theta_out, phi_out, val, inp_is_DSF);
	n = 0;
	while ((c = getc(fp)) != EOF)
		n += !isspace(c);
	if (n) 
		fprintf(stderr,
			"%s: warning: %d unexpected characters past EOD\n",
				fname, n);
	fclose(fp);
	return(1);
}

/* Read in PAB-Opto BSDF files and output RBF interpolant */
int
main(int argc, char *argv[])
{
	extern int	nprocs;
	int		i;
						/* start header */
	SET_FILE_BINARY(stdout);
	newheader("RADIANCE", stdout);
	printargs(argc, argv, stdout);
	fputnow(stdout);
	progname = argv[0];			/* get options */
	while (argc > 2 && argv[1][0] == '-') {
		switch (argv[1][1]) {
		case 'n':
			nprocs = atoi(argv[2]);
			break;
		default:
			goto userr;
		}
		argv += 2; argc -= 2;
	}
	if (argc < 3)
		goto userr;
	for (i = 1; i < argc; i++) {		/* compile measurements */
		if (!load_pabopto_meas(argv[i]))
			return(1);
		make_rbfrep();
	}
	build_mesh();				/* create interpolation */
	save_bsdf_rep(stdout);			/* write it out */
	return(0);
userr:
	fprintf(stderr, "Usage: %s [-n nproc] meas1.dat meas2.dat .. > bsdf.sir\n",
					progname);
	return(1);
}
