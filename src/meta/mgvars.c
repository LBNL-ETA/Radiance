#ifndef lint
static const char	RCSid[] = "$Id: mgvars.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  mgvars.c - routines dealing with graph variables.
 *
 *     6/23/86
 *
 *     Greg Ward Larson
 */

#include  <stdio.h>

#include  <stdlib.h>

#include  <math.h>

#include  <ctype.h>

#include  "mgvars.h"

#define  MAXLINE	512

#define  isid(c)  (isalnum(c) || (c) == '_' || (c) == '.')

#define  isnum(c) (isdigit(c)||(c)=='-'||(c)=='+'||(c)=='.'||(c)=='e'||(c)=='E')

char  *findfile(), *emalloc(), *ecalloc(), *erealloc(), *savestr(), *strcpy();
FILE  *fopen(), *popen();
extern char  *progname, *libpath[];

#ifdef  DCL_ATOF
double  atof();
#endif

IVAR  *ivhead = NULL;			/* intermediate variables */

VARIABLE  gparam[NVARS] = {		/* standard variables */
	{ "fthick", REAL, "frame thickness" },
	{ "grid", REAL, "grid on?" },
	{ "legend", STRING, "legend title" },
	{ "othick", REAL, "origin thickness" },
	{ "period", REAL, "period of polar plot" },
	{ "subtitle", STRING },
	{ "symfile", STRING, "symbol file" },
	{ "tstyle", REAL, "tick mark style" },
	{ "title", STRING },
	{ "xlabel", STRING },
	{ "xmap", FUNCTION, "x axis mapping function" },
	{ "xmax", REAL },
	{ "xmin", REAL },
	{ "xstep", REAL },
	{ "ylabel", STRING },
	{ "ymap", FUNCTION, "y axis mapping function" },
	{ "ymax", REAL },
	{ "ymin", REAL },
	{ "ystep", REAL },
};

VARIABLE  cparam[MAXCUR][NCVARS] = {	/* curve variables */
	{
		{ "A", FUNCTION, "function for curve A" },
		{ "Acolor", REAL, "color for A" },
		{ "Adata", DATA, "point data for A" },
		{ "Alabel", STRING },
		{ "Alintype", REAL, "line type for A" },
		{ "Anpoints", REAL, "number of points for A" },
		{ "Asymsize", REAL, "symbol size for A" },
		{ "Asymtype", STRING, "symbol type for A" },
		{ "Athick", REAL, "line thickness for A" },
	}, {
		{ "B", FUNCTION, "function for curve B" },
		{ "Bcolor", REAL, "color for B" },
		{ "Bdata", DATA, "point data for B" },
		{ "Blabel", STRING },
		{ "Blintype", REAL, "line type for B" },
		{ "Bnpoints", REAL, "number of points for B" },
		{ "Bsymsize", REAL, "symbol size for B" },
		{ "Bsymtype", STRING, "symbol type for B" },
		{ "Bthick", REAL, "line thickness for B" },
	}, {
		{ "C", FUNCTION, "function for curve C" },
		{ "Ccolor", REAL, "color for C" },
		{ "Cdata", DATA, "point data for C" },
		{ "Clabel", STRING },
		{ "Clintype", REAL, "line type for C" },
		{ "Cnpoints", REAL, "number of points for C" },
		{ "Csymsize", REAL, "symbol size for C" },
		{ "Csymtype", STRING, "symbol type for C" },
		{ "Cthick", REAL, "line thickness for C" },
	}, {
		{ "D", FUNCTION, "function for curve D" },
		{ "Dcolor", REAL, "color for D" },
		{ "Ddata", DATA, "point data for D" },
		{ "Dlabel", STRING },
		{ "Dlintype", REAL, "line type for D" },
		{ "Dnpoints", REAL, "number of points for D" },
		{ "Dsymsize", REAL, "symbol size for D" },
		{ "Dsymtype", STRING, "symbol type for D" },
		{ "Dthick", REAL, "line thickness for D" },
	}, {
		{ "E", FUNCTION, "function for curve E" },
		{ "Ecolor", REAL, "color for E" },
		{ "Edata", DATA, "point data for E" },
		{ "Elabel", STRING },
		{ "Elintype", REAL, "line type for E" },
		{ "Enpoints", REAL, "number of points for E" },
		{ "Esymsize", REAL, "symbol size for E" },
		{ "Esymtype", STRING, "symbol type for E" },
		{ "Ethick", REAL, "line thickness for E" },
	}, {
		{ "F", FUNCTION, "function for curve F" },
		{ "Fcolor", REAL, "color for F" },
		{ "Fdata", DATA, "point data for F" },
		{ "Flabel", STRING },
		{ "Flintype", REAL, "line type for F" },
		{ "Fnpoints", REAL, "number of points for F" },
		{ "Fsymsize", REAL, "symbol size for F" },
		{ "Fsymtype", STRING, "symbol type for F" },
		{ "Fthick", REAL, "line thickness for F" },
	}, {
		{ "G", FUNCTION, "function for curve G" },
		{ "Gcolor", REAL, "color for G" },
		{ "Gdata", DATA, "point data for G" },
		{ "Glabel", STRING },
		{ "Glintype", REAL, "line type for G" },
		{ "Gnpoints", REAL, "number of points for G" },
		{ "Gsymsize", REAL, "symbol size for G" },
		{ "Gsymtype", STRING, "symbol type for G" },
		{ "Gthick", REAL, "line thickness for G" },
	}, {
		{ "H", FUNCTION, "function for curve H" },
		{ "Hcolor", REAL, "color for H" },
		{ "Hdata", DATA, "point data for H" },
		{ "Hlabel", STRING },
		{ "Hlintype", REAL, "line type for H" },
		{ "Hnpoints", REAL, "number of points for H" },
		{ "Hsymsize", REAL, "symbol size for H" },
		{ "Hsymtype", STRING, "symbol type for H" },
		{ "Hthick", REAL, "line thickness for H" },
	},
};


mgclearall()			/* clear all variable settings */
{
	int  j;
	register IVAR  *iv;
	register int  i;

	for (iv = ivhead; iv != NULL; iv = iv->next) {
		dremove(iv->name);
		freestr(iv->name);
		freestr(iv->dfn);
		efree((char *)iv);
	}
	ivhead = NULL;

	for (i = 0; i < NVARS; i++)
		if (gparam[i].flags & DEFINED)
			undefine(&gparam[i]);
	
	for (j = 0; j < MAXCUR; j++)
		for (i = 0; i < NCVARS; i++)
			if (cparam[j][i].flags & DEFINED)
				undefine(&cparam[j][i]);
}


mgload(file)			/* load a file */
char  *file;
{
	FILE  *fp;
	char  sbuf[MAXLINE], *fgets();
	int  inquote;
	register char  *cp, *cp2;

	if (file == NULL) {
		fp = stdin;
		file = "<stdin>";
	} else if ((fp = fopen(file, "r")) == NULL) {
		fprintf(stderr, "%s: Cannot open: %s\n", progname, file);
		quit(1);
	}
	while (fgets(sbuf+1, sizeof(sbuf)-1, fp) != NULL) {
		inquote = 0;
		cp2 = sbuf;
		for (cp = sbuf+1; *cp; cp++)	/* condition the input line */
			switch (*cp) {
			case '#':
				if (!inquote) {
					cp[0] = '\n';
					cp[1] = '\0';
					break;
				}
				*cp2++ = *cp;
				break;
			case '"':
				inquote = !inquote;
				break;
			case '\\':
				if (!cp[1])
					break;
				if (cp[1] == '\n') {
					cp[0] = '\0';
					fgets(cp, sizeof(sbuf)-(cp-sbuf), fp);
					cp--;
					break;
				}
				*cp2++ = *++cp;
				break;
			case ' ':
			case '\t':
			case '\n':
				if (!inquote)
					break;
				*cp2++ = *cp;
				break;
			default:
				*cp2++ = *cp;
				break;
			}
		*cp2 = '\0';
		if (inquote) {
			fputs(sbuf, stderr);
			fprintf(stderr, "%s: %s: Missing quote\n",
					progname, file);
			quit(1);
		}
		if (sbuf[0])
			setmgvar(file, fp, sbuf);
	}
	if (fp != stdin)
		fclose(fp);
}


mgsave(file)				/* save our variables */
char  *file;
{
	FILE  *fp;
	int  j;
	register IVAR  *iv;
	register int  i;

	if (file == NULL)
		fp = stdout;
	else if ((fp = fopen(file, "w")) == NULL) {
		fprintf(stderr, "%s: Cannot write: %s\n", progname, file);
		quit(1);
	}
	for (iv = ivhead; iv != NULL; iv = iv->next)
		fprintf(fp, "%s\n", iv->dfn);

	for (i = 0; i < NVARS; i++)
		if (gparam[i].flags & DEFINED)
			mgprint(&gparam[i], fp);
	
	for (j = 0; j < MAXCUR; j++)
		for (i = 0; i < NCVARS; i++)
			if (cparam[j][i].flags & DEFINED)
				mgprint(&cparam[j][i], fp);

	if (fp != stdout)
		fclose(fp);
}


setmgvar(fname, fp, string)		/* set a variable */
char  *fname;
FILE  *fp;
char  *string;
{
	char  name[128];
	FILE  *fp2;
	register int  i;
	register char  *s;
	register VARIABLE  *vp;

	if (!strncmp(string, "include=", 8)) {	/* include file */
		if ((s = findfile(string+8, libpath)) == NULL) {
			fprintf(stderr, "%s\n", string);
			fprintf(stderr, "%s: %s: File not found: %s\n",
					progname, fname, string+8);
			quit(1);
		}
		strcpy(name, s);
		mgload(name);
		return;
	}
	s = string;
	i = 0;
	while (i < sizeof(name)-1 && isid(*s))
		name[i++] = *s++;
	name[i] = '\0';
	vp = vlookup(name);
	if (vp != NULL) {
		undefine(vp);
		switch (vp->type) {
		case REAL:
		case FUNCTION:
			if ((*s == '(') != (vp->type == FUNCTION)) {
				fprintf(stderr, "%s\n", string);
				fprintf(stderr,
					"%s: %s: Bad %s declaration: %s\n",
					progname, fname,
					vp->type == FUNCTION ?
					"function" : "variable",
					name);
				quit(1);
			}
			scompile(string, fname, 0);
			vp->v.dfn = savestr(string);
			break;
		case STRING:
			if (*s++ != '=') {
				fprintf(stderr, "%s\n", string);
				fprintf(stderr, "%s: %s: Missing '='\n",
						progname, fname);
				quit(1);
			}
			vp->v.s = savestr(s);
			break;
		case DATA:
			if (*s++ != '=') {
				fprintf(stderr, "%s\n", string);
				fprintf(stderr, "%s: %s: Missing '='\n",
						progname, fname);
				quit(1);
			}
			if (!*s) {
				loaddata(fname, fp, &vp->v.d);
#ifdef  UNIX
			} else if (*s == '!') {
				if ((fp2 = popen(s+1, "r")) == NULL) {
					fprintf(stderr, "%s\n", string);
					fprintf(stderr,
					"%s: %s: Cannot execute: %s\n",
							progname, fname, s+1);
					quit(1);
				}
				loaddata(s, fp2, &vp->v.d);
				pclose(fp2);
#endif
			} else {
				if ((fp2 = fopen(s, "r")) == NULL) {
				    fprintf(stderr, "%s\n", string);
				    fprintf(stderr,
					    "%s: %s: Data file not found: %s\n",
					    progname, fname, s);
				    quit(1);
				}
				loaddata(s, fp2, &vp->v.d);
				fclose(fp2);
			}
			break;
		}
		vp->flags |= DEFINED;
	} else
		setivar(name, fname, string);		/* intermediate */
}


setivar(vname, fname, definition)	/* set an intermediate variable */
char  *vname;
char  *fname;
char  *definition;
{
	IVAR  ivbeg;
	register IVAR  *iv;

	scompile(definition, fname, 0);		/* compile the string */

	ivbeg.next = ivhead;
	for (iv = &ivbeg; iv->next != NULL; iv = iv->next)
		if (!strcmp(vname, iv->next->name)) {
			iv = iv->next;
			freestr(iv->dfn);
			iv->dfn = savestr(definition);
			return;
		}

	iv->next = (IVAR *)emalloc(sizeof(IVAR));
	iv = iv->next;
	iv->name = savestr(vname);
	iv->dfn = savestr(definition);
	iv->next = NULL;
	ivhead = ivbeg.next;
}


mgtoa(s, vp)			/* get a variable's value in ascii form */
register char  *s;
VARIABLE  *vp;
{
	register  char  *sv;

	if (!(vp->flags & DEFINED)) {
		strcpy(s, "UNDEFINED");
		return;
	}
	switch (vp->type) {
	case REAL:
	case FUNCTION:
		sv = vp->v.dfn;
		while (*sv != '=' && *sv != ':')
			sv++;
		while (*++sv && *sv != ';')
			*s++ = *sv;
		*s = '\0';
		break;
	case STRING:
		strcpy(s, vp->v.s);
		break;
	case DATA:
		strcpy(s, "DATA");
		break;
	}
}


mgprint(vp, fp)			/* print a variable definition */
register VARIABLE  *vp;
FILE  *fp;
{
	register int  i;
	
	switch (vp->type) {
	case REAL:
	case FUNCTION:
		fprintf(fp, "%s\n", vp->v.dfn);
		break;
	case STRING:
		fprintf(fp, "%s=\"", vp->name);
		for (i = 0; vp->v.s[i]; i++)
			switch (vp->v.s[i]) {
			case '"':
			case '\\':
				putc('\\', fp);
			/* fall through */
			default:
				putc(vp->v.s[i], fp);
				break;
			}
		fprintf(fp, "\"\n");
		break;
	case DATA:
		fprintf(fp, "%s=", vp->name);
		for (i = 0; i < vp->v.d.size; i++) {
			if (i % 4 == 0)
				fprintf(fp, "\n");
			fprintf(fp, "\t%10e", vp->v.d.data[i]);
		}
		fprintf(fp, "\n;\n");
		break;
	}
}


undefine(vp)			/* undefine a variable */
register VARIABLE  *vp;
{
	if (vp == NULL || !(vp->flags & DEFINED))
		return;

	switch (vp->type) {
	case REAL:
	case FUNCTION:
		dremove(vp->name);
		freestr(vp->v.dfn);
		break;
	case STRING:
		freestr(vp->v.s);
		break;
	case DATA:
		efree((char *)vp->v.d.data);
		break;
	}
	vp->flags &= ~DEFINED;
}


VARIABLE *
vlookup(vname)			/* look up a variable by its name */
char  *vname;
{
	register int  i;
	register VARIABLE  *vp;

	i = vname[0] - 'A';
	if (i >= 0 && i < MAXCUR)		/* curve variables */
		for (vp = cparam[i], i = 0; i < NCVARS; vp++, i++)
			if (!strcmp(vp->name, vname))
				return(vp);
						/* standard variables */
	for (vp = gparam; vp < &gparam[NVARS]; vp++)
		if (!strcmp(vp->name, vname))
			return(vp);
	return(NULL);				/* not found */
}


loaddata(fname, fp, dp)			/* load data from a stream */
char  *fname;
FILE  *fp;
register DARRAY  *dp;
{
	char  sbuf[MAXLINE], *fgets();
	register char  *cp;

	dp->size = 0;
	dp->data = NULL;
	while (fgets(sbuf, sizeof(sbuf), fp) != NULL) {
		cp = sbuf;
		while (*cp) {
			while (isspace(*cp) || *cp == ',')
				cp++;
			if (isnum(*cp)) {
				dp->data = (float *)erealloc((char *)dp->data,
						(dp->size+1)*sizeof(float));
				dp->data[dp->size++] = atof(cp);
				do
					cp++;
				while (isnum(*cp));
			} else if (*cp == ';') {
				return;
			} else if (*cp) {
				fputs(sbuf, stderr);
				fprintf(stderr, "%s: %s: Bad data\n",
						progname, fname);
				quit(1);
			}
		}
	}
}


char *
findfile(fname, pathlist)	/* find the file fname, return full path */
char  *fname;
register char  **pathlist;
{
	static char  ffname[128];
	char  *strcpy(), *strcat();
	register int  fd;

	if (fname[0] == '/')
		return(fname);

	while (*pathlist != NULL) {
		strcpy(ffname, *pathlist);
		strcat(ffname, fname);
		if ((fd = open(ffname, 0)) != -1) {
			close(fd);
			return(ffname);
		}
		pathlist++;
	}
	return(NULL);
}


int
mgcurve(c, f)			/* get a curve's (unmapped) values */
int  c;
int  (*f)();
{
	int  nargs;
	double  x[2], step;
	register VARIABLE  *cv;
	register float  *p;
	register int  npts;

	if (c < 0 || c >= MAXCUR)
		return(-1);
	cv = cparam[c];

	if (cv[C].flags & DEFINED) {		/* function or map */

		nargs = fundefined(cv[C].name);
		if (nargs < 1 || nargs > 2) {
			fprintf(stderr, "%s: bad # of arguments for '%c'\n",
					progname, c+'A');
			quit(1);
		}

		if (cv[CDATA].flags & DEFINED) {		/* map */
			npts = cv[CDATA].v.d.size / nargs;
			p = cv[CDATA].v.d.data;
			while (npts--) {
				x[0] = *p++;
				if (nargs == 2)
					x[1] = *p++;
				(*f)(c, x[0],
					funvalue(cv[C].name, nargs, x));
			}
			npts = cv[CDATA].v.d.size / nargs;
		} else if ( nargs == 1 &&			/* function */
				gparam[XMIN].flags & DEFINED &&
				gparam[XMAX].flags & DEFINED &&
				cv[CNPOINTS].flags & DEFINED ) {
			npts = varvalue(cv[CNPOINTS].name);
			if (npts > 1)
				step = (varvalue(gparam[XMAX].name) -
					varvalue(gparam[XMIN].name)) /
					(npts - 1);
			else
				step = 0.0;
			for (x[0] = varvalue(gparam[XMIN].name);
					npts--; x[0] += step)
				(*f)(c, x[0],
					funvalue(cv[C].name, 1, x));
			npts = varvalue(cv[CNPOINTS].name);
		} else {
			fprintf(stderr,
			"%s: function '%c' needs %cdata or xmin, xmax, %cnpoints\n",
					progname, c+'A', c+'A', c+'A');
			quit(1);
		}

	} else if (cv[CDATA].flags & DEFINED) {	/* data */

		npts = cv[CDATA].v.d.size / 2;
		p = cv[CDATA].v.d.data;
		while (npts--) {
			(*f)(c, p[0], p[1]);
			p += 2;
		}
		npts = cv[CDATA].v.d.size / 2;

	} else

		npts = 0;

	return(npts);
}
