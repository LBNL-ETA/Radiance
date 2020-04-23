#ifndef lint
static const char RCSid[] = "$Id: readwfobj.c,v 2.2 2020/04/23 03:19:48 greg Exp $";
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

typedef int	VNDX[3];	/* vertex index (point,map,normal) */

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

/* Add a vertex to our scene */
static void
add_vertex(Scene *sc, double x, double y, double z)
{
	sc->vert = chunk_alloc(Vertex, sc->vert, sc->nverts);
	sc->vert[sc->nverts].p[0] = x;
	sc->vert[sc->nverts].p[1] = y;
	sc->vert[sc->nverts].p[2] = z;
	sc->vert[sc->nverts++].vflist = NULL;
}

/* Add a texture coordinate to our scene */
static void
add_texture(Scene *sc, double u, double v)
{
	sc->tex = chunk_alloc(TexCoord, sc->tex, sc->ntex);
	sc->tex[sc->ntex].u = u;
	sc->tex[sc->ntex].v = v;
	sc->ntex++;
}

/* Add a surface normal to our scene */
static int
add_normal(Scene *sc, double xn, double yn, double zn)
{
	FVECT   nrm;

	nrm[0] = xn; nrm[1] = yn; nrm[2] = zn;
	if (normalize(nrm) == .0)
		return(0);
	sc->norm = chunk_alloc(Normal, sc->norm, sc->nnorms);
	VCOPY(sc->norm[sc->nnorms], nrm);
	sc->nnorms++;
	return(1);
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

/* set current group */
static void
set_group(Scene *sc, const char *nm)
{
	sc->lastgrp = findName(nm, (const char **)sc->grpname, sc->ngrps);
	if (sc->lastgrp >= 0)
		return;
	sc->grpname = chunk_alloc(char *, sc->grpname, sc->ngrps);
	sc->grpname[sc->lastgrp=sc->ngrps++] = savqstr((char *)nm);
}

/* set current material */
static void
set_material(Scene *sc, const char *nm)
{
	sc->lastmat = findName(nm, (const char **)sc->matname, sc->nmats);
	if (sc->lastmat >= 0)
		return;
	sc->matname = chunk_alloc(char *, sc->matname, sc->nmats);
	sc->matname[sc->lastmat=sc->nmats++] = savqstr((char *)nm);
}

/* Add a new face to scene */
static int
add_face(Scene *sc, const VNDX ondx, int ac, char *av[])
{
	Face    *f;
	int     i;
	
	if (ac < 3)
		return(0);
	f = (Face *)emalloc(sizeof(Face)+sizeof(VertEnt)*(ac-3));
	f->flags = 0;
	f->nv = ac;
	f->grp = sc->lastgrp;
	f->mat = sc->lastmat;
	for (i = 0; i < ac; i++) {		/* add each vertex */
		VNDX    vin;
		int     j;
		if (!cvtndx(vin, sc, ondx, av[i])) {
			efree((char *)f);
			return(0);
		}
		f->v[i].vid = vin[0];
		f->v[i].tid = vin[1];
		f->v[i].nid = vin[2];
		f->v[i].fnext = NULL;
		for (j = i; j-- > 0; )
			if (f->v[j].vid == vin[0])
				break;
		if (j < 0) {			/* first occurrence? */
			f->v[i].fnext = sc->vert[vin[0]].vflist;
			sc->vert[vin[0]].vflist = f;
		} else if (ac == 3)		/* degenerate triangle? */
			f->flags |= FACE_DEGENERATE;
	}
	f->next = sc->flist;			/* push onto face list */
	sc->flist = f;
	sc->nfaces++;
						/* check face area */
	if (!(f->flags & FACE_DEGENERATE) && faceArea(sc, f, NULL) <= FTINY)
		f->flags |= FACE_DEGENERATE;
	return(1);
}

/* Load a .OBJ file */
Scene *
loadOBJ(Scene *sc, const char *fspec)
{
	FILE    *fp;
	char	*argv[MAXARG];
	int	argc;
	char    buf[256];
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
				add_vertex(sc, atof(argv[1]),
						atof(argv[2]),
						atof(argv[3]));
				break;
			case 't':			/* texture coord. */
				if (argv[0][2])
					goto unknown;
				if (badarg(argc-1,argv+1,"ff"))
					goto unknown;
				add_texture(sc, atof(argv[1]), atof(argv[2]));
				break;
			case 'n':			/* normal */
				if (argv[0][2])
					goto unknown;
				if (badarg(argc-1,argv+1,"fff")) {
					syntax(fspec, "bad normal");
					goto failure;
				}
				if (!add_normal(sc, atof(argv[1]),
						atof(argv[2]),
						atof(argv[3]))) {
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
			set_material(sc, argv[1]);
			break;
		case 'o':		/* object name */
		case 'g':		/* group name */
			if (argc < 2) {
				syntax(fspec, "missing argument");
				goto failure;
			}
			set_group(sc, group_name(argc-1, argv+1));
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
	if (fspec[0] == '!')
		pclose(fp);
	else
#endif
	if (fp != stdin)
		fclose(fp);
	sprintf(buf, "%d statements read from \"%s\"", nstats, fspec);
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
