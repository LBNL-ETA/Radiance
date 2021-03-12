#ifndef lint
static const char RCSid[] = "$Id: robjutil.c,v 2.4 2021/03/12 18:32:33 greg Exp $";
#endif
/*
 * Utility program for fixing up Wavefront .OBJ files.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rtio.h"
#include "objutil.h"

const char      *progname;
int		verbose = 0;

int
main(int argc, char *argv[])
{
	Scene		*myScene;
	int		radout = 0;
	const char      *grp[64];
	int		ngrps = 0;
	int		save_grps = 1;
	const char      *mat[64];
	int		nmats = 0;
	int		save_mats = 1;
	int		do_tex = 0;
	int		do_norm = 0;
	int		do_triangulate = 0;
	char		*xfm = NULL;
	char		cbuf[256];
	double		verteps = -1.;
	int     i, n;
	
	progname = argv[0];
						/* process options */
	for (i = 1; i < argc && (argv[i][0] == '-' || argv[i][0] == '+'); i++)
		switch (argv[i][1]) {
		case 'v':			/* verbose mode */
			verbose = (argv[i][0] == '+');
			break;
		case 'g':			/* save/delete group */
			if (save_grps != (argv[i][0] == '+')) {
				if (ngrps) {
					fprintf(stderr,
					"%s: -g and +g are mutually exclusive\n",
						argv[0]);
					exit(1);
				}
				save_grps = !save_grps;
			}
			if (ngrps >= 64) {
				fprintf(stderr, "%s: too many groups\n",
						argv[0]);
				exit(1);
			}
			grp[ngrps++] = argv[++i];
			break;
		case 'm':			/* save/delete material */
			if (save_mats != (argv[i][0] == '+')) {
				if (nmats) {
					fprintf(stderr,
					"%s: -m and +m are mutually exclusive\n",
						argv[0]);
					exit(1);
				}
				save_mats = !save_mats;
			}
			if (nmats >= 64) {
				fprintf(stderr, "%s: too many materials\n",
						argv[0]);
				exit(1);
			}
			mat[nmats++] = argv[++i];
			break;
		case 't':			/* do texture coord's */
			do_tex = (argv[i][0] == '+') ? 1 : -1;
			break;
		case 'n':			/* do normals */
			do_norm = (argv[i][0] == '+') ? 1 : -1;
			break;
		case 'c':			/* coalesce vertices */
			verteps = atof(argv[++i]);
			break;
		case 'r':			/* output to Radiance file? */
			radout = (argv[i][0] == '+');
			break;
		case 'T':			/* triangulate faces? */
			do_triangulate = (argv[i][0] == '+');
			break;
		case 'x':			/* apply a transform */
			if (xfm != NULL) {
				fprintf(stderr, "%s: only one '-x' option allowed\n",
						argv[0]);
				exit(1);
			}
			xfm = argv[++i];
			break;
		default:
			fprintf(stderr, "%s: unknown option: %s\n",
					argv[0], argv[i]);
			goto userr;
		}
	if (i == argc) {
		if (verbose)
			fputs("Loading scene from standard input...\n",
					stderr);
		myScene = loadOBJ(NULL, NULL);
	} else {
		myScene = newScene();
		for ( ; i < argc; i++) {
			if (verbose)
				fprintf(stderr,
					"Loading scene from \"%s\"...\n",
						argv[i]);
			if (loadOBJ(myScene, argv[i]) == NULL)
				exit(1);
		}
	}
	if (ngrps) {
		if (verbose)
			fputs("Deleting unwanted groups...\n", stderr);
		clearSelection(myScene, 0);
		sprintf(cbuf, "%s the following groups:",
				save_grps ? "Extracting" : "Deleting");
		addComment(myScene, cbuf);
		for (i = 0; i < ngrps; i++) {
			sprintf(cbuf, "\t%s", grp[i]);
			addComment(myScene, cbuf);
			selectGroup(myScene, grp[i], 0);
		}
		if (save_grps)
			n = deleteFaces(myScene, 0, FACE_SELECTED);
		else
			n = deleteFaces(myScene, FACE_SELECTED, 0);
		sprintf(cbuf, "\t(%d faces removed)", n);
		addComment(myScene, cbuf);
	}
	if (nmats) {
		if (verbose)
			fputs("Deleting unwanted materials...\n", stderr);
		clearSelection(myScene, 0);
		sprintf(cbuf, "%s the following materials:",
				save_mats ? "Extracting" : "Deleting");
		addComment(myScene, cbuf);
		for (i = 0; i < nmats; i++) {
			sprintf(cbuf, "\t%s", mat[i]);
			addComment(myScene, cbuf);
			selectMaterial(myScene, mat[i], 0);
		}
		if (save_mats)
			n = deleteFaces(myScene, 0, FACE_SELECTED);
		else
			n = deleteFaces(myScene, FACE_SELECTED, 0);
		sprintf(cbuf, "\t(%d faces removed)", n);
		addComment(myScene, cbuf);
	}
	if (verbose && do_tex < 0)
		fputs("Removing texture coordinates...\n", stderr);
	if (do_tex < 0 && (n = removeTexture(myScene, 0, 0))) {
		sprintf(cbuf, "Removed texture coordinates from %d faces", n);
		addComment(myScene, cbuf);
	}
	if (verbose && do_norm < 0)
		fputs("Removing surface normals...\n", stderr);
	if (do_norm < 0 && (n = removeNormals(myScene, 0, 0))) {
		sprintf(cbuf, "Removed surface normals from %d faces", n);
		addComment(myScene, cbuf);
	}
	if (verbose && verteps >= 0)
		fputs("Coalescing vertices...\n", stderr);
	if (verteps >= 0 && (n = coalesceVertices(myScene, verteps))) {
		sprintf(cbuf, "Coalesced %d vertices with epsilon of %g",
				n, verteps);
		addComment(myScene, cbuf);
	}
	if (verbose)
		fputs("Deleting degenerate faces...\n", stderr);
	n = deleteFaces(myScene, FACE_DEGENERATE, 0);
	if (n) {
		sprintf(cbuf, "Removed %d degenerate faces", n);
		addComment(myScene, cbuf);
	}
	if (verbose)
		fputs("Checking for duplicate faces...\n", stderr);
	if (findDuplicateFaces(myScene))
		n = deleteFaces(myScene, FACE_DUPLICATE, 0);
	if (n) {
		sprintf(cbuf, "Removed %d duplicate faces", n);
		addComment(myScene, cbuf);
	}
	if (do_triangulate) {
		if (verbose)
			fputs("Making sure all faces are triangles...\n", stderr);
		n = triangulateScene(myScene);
		if (n > 0) {
			sprintf(cbuf, "Added %d faces during triangulation", n);
			addComment(myScene, cbuf);
		}
	}
	if (xfm != NULL) {
		if (verbose)
			fputs("Applying transform...\n", stderr);
		if (!xfmScene(myScene, xfm)) {
			fprintf(stderr, "%s: transform error\n", argv[0]);
			exit(1);
		}
	}
	if (verbose)
		fputs("Writing out scene...\n", stderr);

	fputs("# File processed by: ", stdout);
	printargs(argc, argv, stdout);
	if (radout)
		toRadiance(myScene, stdout, 0, 0);
	else
		toOBJ(myScene, stdout);
	/* freeScene(myScene); */
	if (verbose)
		fputs("Done.                                    \n", stderr);
	return(0);
userr:
	fprintf(stderr, "Usage: %s [options] [input.obj ..]\n", argv[0]);
	fprintf(stderr, "Available options:\n");
	fprintf(stderr, "\t+/-r\t\t\t# Radiance scene output on/off\n");
	fprintf(stderr, "\t+/-v\t\t\t# on/off verbosity (progress reports)\n");
	fprintf(stderr, "\t+/-g name\t\t# save/delete group\n");
	fprintf(stderr, "\t+/-m name\t\t# save/delete faces with material\n");
	fprintf(stderr, "\t+/-t\t\t\t# keep/remove texture coordinates\n");
	fprintf(stderr, "\t+/-n\t\t\t# keep/remove vertex normals\n");
	fprintf(stderr, "\t-c epsilon\t\t# coalesce vertices within epsilon\n");
	fprintf(stderr, "\t+T\t\t\t# turn all faces into triangles\n");
	fprintf(stderr, "\t-x 'xf spec'\t# apply the quoted transform\n");
	return(1);
}

void
eputs(char *s)				/* put string to stderr */
{
	static int  midline = 0;

	if (!*s)
		return;
	if (!midline++) {
		fputs(progname, stderr);
		fputs(": ", stderr);
	}
	fputs(s, stderr);
	if (s[strlen(s)-1] == '\n')
		midline = 0;
}
