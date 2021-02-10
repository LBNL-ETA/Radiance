#ifndef lint
static const char RCSid[] = "$Id: readwfobj.c,v 2.8 2021/02/10 17:57:28 greg Exp $";
#endif
/*
 *  readobj.c
 *
 *  Routines for reading a Wavefront .OBJ file
 *
 *  Created by Greg Ward on Wed Feb 11 2004.
 */

#include "rtio.h"
#include "rterror.h"
#include "fvect.h"
#include "paths.h"
#include <ctype.h>
#include "objutil.h"

#define MAXARG		512	/* maximum # arguments in a statement */

static int	lineno;		/* current line number */

/* read the next statement from fp */
static int
get_stmt(char *av[MAXARG], FILE *fp)
{
	static char	sbuf[MAXARG*16];
	register char	*cp;
	register int	i;

	do {
		if (fgetline(cp=sbuf, sizeof(sbuf), fp) == NULL)
			return(0);
		i = 0;
		for ( ; ; ) {
			while (isspace(*cp) || *cp == '\\') {
				if (*cp == '\n')
					lineno++;
				*cp++ = '\0';
			}
			if (!*cp)
				break;
			if (i >= MAXARG-1) {
				sprintf(errmsg,
				"Too many arguments near line %d (limit %d)",
					lineno+1, MAXARG-1);
				error(WARNING, errmsg);
				break;
			}
			av[i++] = cp;
			while (*++cp && !isspace(*cp))
				;
		}
		av[i] = NULL;
		lineno++;
	} while (!i);

	return(i);
}

/* convert vertex string to index */
static int
cvtndx(VNDX vi, const Scene *sc, const VNDX ondx, const char *vs)
{
					/* get point */
	vi[0] = atoi(vs);
	if (vi[0] > 0) {
		if ((vi[0] += ondx[0] - 1) >= sc->nverts)
			return(0);
	} else if (vi[0] < 0) {
		vi[0] += sc->nverts;
		if (vi[0] < 0)
			return(0);
	} else
		return(0);
					/* get map coord. */
	while (*vs)
		if (*vs++ == '/')
			break;
	vi[1] = atoi(vs);
	if (vi[1] > 0) {
		if ((vi[1] += ondx[1] - 1) >= sc->ntex)
			return(0);
	} else if (vi[1] < 0) {
		vi[1] += sc->ntex;
		if (vi[1] < 0)
			return(0);
	} else
		vi[1] = -1;
					/* get normal */
	while (*vs)
		if (*vs++ == '/')
			break;
	vi[2] = atoi(vs);
	if (vi[2] > 0) {
		if ((vi[2] += ondx[2] - 1) >= sc->nnorms)
			return(0);
	} else if (vi[2] < 0) {
		vi[2] += sc->nnorms;
		if (vi[2] < 0)
			return(0);
	} else
		vi[2] = -1;
	return(1);
}

/* report syntax error */
static void
syntax(const char *fn, const char *er)
{
	sprintf(errmsg, "%s: Wavefront syntax error near line %d: %s",
			fn, lineno, er);
	error(USER, errmsg);
}

/* combine multi-group name into single identifier w/o spaces */
static char *
group_name(int ac, char **av)
{
	static char     nambuf[256];
	char		*cp;
	if (ac < 1)
		return("NO_NAME");
	if (ac == 1)
		return(av[0]);
	for (cp = nambuf; ac--; av++) {
		strcpy(cp, *av);
		while (*cp) {
			/* XXX white space disallowed by current get_stmt()
			if (isspace(*cp))
				*cp = '_';
			else
			*/
			if (*cp == '.')
				*cp = '_';
			++cp;
		}
		*cp++ = '.';
	}
	*--cp = '\0';
	return(nambuf);
}

/* add a new face to scene */
static int
add_face(Scene *sc, const VNDX ondx, int ac, char *av[])
{
	VNDX	vdef[4];
	VNDX	*varr = vdef;
	Face    *f = NULL;
	int     i;
	
	if (ac < 3)				/* legal polygon? */
		return(0);
	if (ac > 4)				/* need to allocate array? */
		varr = (VNDX *)emalloc(ac*sizeof(VNDX));
	for (i = ac; i--; )			/* index each vertex */
		if (!cvtndx(varr[i], sc, ondx, av[i]))
			break;
	if (i < 0)				/* create face if indices OK */
		f = addFace(sc, varr, ac);
	if (varr != vdef)
		efree((char *)varr);
	return(f != NULL);
}

/* Load a .OBJ file */
Scene *
loadOBJ(Scene *sc, const char *fspec)
{
	FILE    *fp;
	char	*argv[MAXARG];
	int	argc;
	char    buf[1024];
	int	nstats=0, nunknown=0;
	int     onfaces;
	VNDX    ondx;

	if (fspec == NULL) {
		fp = stdin;
		fspec = "<stdin>";
#if POPEN_SUPPORT
	} else if (fspec[0] == '!') {
		if ((fp = popen(fspec+1, "r")) == NULL) {
			sprintf(errmsg, "%s: cannot execute", fspec);
			error(SYSTEM, errmsg);
			return(NULL);
		}
#endif
	} else if ((fp = fopen(fspec, "r")) == NULL) {
		sprintf(errmsg, "%s: cannot open for reading", fspec);
		error(SYSTEM, errmsg);
		return(NULL);
	}
	if (sc == NULL)
		sc = newScene();
	lineno = 0;
	onfaces = sc->nfaces;
	ondx[0] = sc->nverts;
	ondx[1] = sc->ntex;
	ondx[2] = sc->nnorms;
	while ((argc = get_stmt(argv, fp)) > 0) {
		switch (argv[0][0]) {
		case 'v':		/* vertex */
			switch (argv[0][1]) {
			case '\0':			/* point */
				if (badarg(argc-1,argv+1,"fff")) {
					syntax(fspec, "bad vertex");
					goto failure;
				}
				addVertex(sc, atof(argv[1]),
						atof(argv[2]),
						atof(argv[3]));
				break;
			case 't':			/* texture coord. */
				if (argv[0][2])
					goto unknown;
				if (badarg(argc-1,argv+1,"ff"))
					goto unknown;
				addTexture(sc, atof(argv[1]), atof(argv[2]));
				break;
			case 'n':			/* normal */
				if (argv[0][2])
					goto unknown;
				if (badarg(argc-1,argv+1,"fff")) {
					syntax(fspec, "bad normal");
					goto failure;
				}
				if (addNormal(sc, atof(argv[1]),
						atof(argv[2]),
						atof(argv[3])) < 0) {
					syntax(fspec, "zero normal");
					goto failure;
				}
				break;
			default:
				goto unknown;
			}
			break;
		case 'f':				/* face */
			if (argv[0][1])
				goto unknown;
			if (!add_face(sc, ondx, argc-1, argv+1)) {
				syntax(fspec, "bad face");
				goto failure;
			}
			break;
		case 'u':				/* usemtl/usemap */
			if (!strcmp(argv[0], "usemap"))
				break;
			if (strcmp(argv[0], "usemtl"))
				goto unknown;
			if (argc != 2) {
				syntax(fspec, "bad # arguments");
				goto failure;
			}
			setMaterial(sc, argv[1]);
			break;
		case 'o':		/* object name */
		case 'g':		/* group name */
			if (argc < 2) {
				syntax(fspec, "missing argument");
				goto failure;
			}
			setGroup(sc, group_name(argc-1, argv+1));
			break;
		case '#':		/* comment */
			continue;
		default:;		/* something we don't deal with */
		unknown:
			nunknown++;
			break;
		}
		nstats++;
		if (verbose && !(nstats & 0x3fff))
			fprintf(stderr, " %8d statements\r", nstats);
	}
#if POPEN_SUPPORT
	if (fspec[0] == '!') {
		if (pclose(fp) != 0) {
			sprintf(errmsg, "Bad return status from: %s", fspec+1);
			error(USER, errmsg);
			freeScene(sc);
			return(NULL);
		}
	} else
#endif
	if (fp != stdin)
		fclose(fp);
	if (verbose)
		fprintf(stderr, "Read %d statements\n", nstats);
	if (strlen(fspec) < sizeof(buf)-32)
		sprintf(buf, "%d statements read from \"%s\"", nstats, fspec);
	else
		sprintf(buf, "%d statements read from (TOO LONG TO SHOW)", nstats);
	addComment(sc, buf);
	if (nunknown) {
		sprintf(buf, "\t%d unrecognized", nunknown);
		addComment(sc, buf);
	}
	sprintf(buf, "%d faces", sc->nfaces - onfaces);
	addComment(sc, buf);
	sprintf(buf, "\t%d vertices, %d texture coordinates, %d surface normals",
			sc->nverts - ondx[0],
			sc->ntex - ondx[1],
			sc->nnorms - ondx[2]);
	addComment(sc, buf);
	return(sc);
failure:
#if POPEN_SUPPORT
	if (fspec[0] == '!')
		pclose(fp);
	else
#endif
	if (fp != stdin)
		fclose(fp);
	freeScene(sc);
	return(NULL);
}
