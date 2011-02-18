#ifndef lint
static const char	RCSid[] = "$Id: mgf_xf.c,v 3.1 2011/02/18 00:40:25 greg Exp $";
#endif
/*
 * Routines for 4x4 homogeneous, rigid-body transformations
 */

#include <stdio.h>
#include <stdlib.h>
#include "mgf_parser.h"

#define  d2r(a)		((PI/180.)*(a))

#define  checkarg(a,l)	if (av[i][a] || badarg(ac-i-1,av+i+1,l)) goto done

XF_SPEC	*xf_context = NULL;	/* current context */
char	**xf_argend;		/* end of transform argument list */
static char	**xf_argbeg;	/* beginning of transform argument list */


int
xf_handler(int ac, char **av)		/* handle xf entity */
{
	register XF_SPEC	*spec;
	register int	n;
	int	rv;

	if (ac == 1) {			/* something with existing transform */
		if ((spec = xf_context) == NULL)
			return(MG_ECNTXT);
		n = -1;
		if (spec->xarr != NULL) {	/* check for iteration */
			register struct xf_array	*ap = spec->xarr;

			(void)xf_aname((struct xf_array *)NULL);
			n = ap->ndim;
			while (n--) {
				if (++ap->aarg[n].i < ap->aarg[n].n)
					break;
				(void)strcpy(ap->aarg[n].arg, "0");
				ap->aarg[n].i = 0;
			}
			if (n >= 0) {
				if ((rv = mg_fgoto(&ap->spos)) != MG_OK)
					return(rv);
				sprintf(ap->aarg[n].arg, "%d", ap->aarg[n].i);
				(void)xf_aname(ap);
			}
		}
		if (n < 0) {			/* pop transform */
			xf_context = spec->prev;
			free_xf(spec);
			return(MG_OK);
		}
	} else {			/* else allocate transform */
		if ((spec = new_xf(ac-1, av+1)) == NULL)
			return(MG_EMEM);
		if (spec->xarr != NULL)
			(void)xf_aname(spec->xarr);
		spec->prev = xf_context;	/* push onto stack */
		xf_context = spec;
	}
					/* translate new specification */
	n = xf_ac(spec);
	n -= xf_ac(spec->prev);		/* incremental comp. is more eff. */
	if (xf(&spec->xf, n, xf_av(spec)) != n)
		return(MG_ETYPE);
					/* check for vertex reversal */
	if ((spec->rev = (spec->xf.sca < 0.)))
		spec->xf.sca = -spec->xf.sca;
					/* compute total transformation */
	if (spec->prev != NULL) {
		multmat4(spec->xf.xfm, spec->xf.xfm, spec->prev->xf.xfm);
		spec->xf.sca *= spec->prev->xf.sca;
		spec->rev ^= spec->prev->rev;
	}
	spec->xid = comp_xfid(spec->xf.xfm);	/* compute unique ID */
	return(MG_OK);
}


XF_SPEC *
new_xf(int ac, char **av)			/* allocate new transform structure */
{
	register XF_SPEC	*spec;
	register int	i;
	char	*cp;
	int	n, ndim;

	ndim = 0;
	n = 0;				/* compute space req'd by arguments */
	for (i = 0; i < ac; i++)
		if (!strcmp(av[i], "-a")) {
			ndim++;
			i++;
		} else
			n += strlen(av[i]) + 1;
	if (ndim > XF_MAXDIM)
		return(NULL);
	spec = (XF_SPEC *)malloc(sizeof(XF_SPEC) + n);
	if (spec == NULL)
		return(NULL);
	if (ndim) {
		spec->xarr = (struct xf_array *)malloc(sizeof(struct xf_array));
		if (spec->xarr == NULL)
			return(NULL);
		mg_fgetpos(&spec->xarr->spos);
		spec->xarr->ndim = 0;		/* incremented below */
	} else
		spec->xarr = NULL;
	spec->xac = ac + xf_argc;
					/* and store new xf arguments */
	if (xf_argbeg == NULL || xf_av(spec) < xf_argbeg) {
		register char	**newav =
				(char **)malloc((spec->xac+1)*sizeof(char *));
		if (newav == NULL)
			return(NULL);
		for (i = xf_argc; i-- > 0; )
			newav[ac+i] = xf_argend[i-xf_context->xac];
		*(xf_argend = newav + spec->xac) = NULL;
		if (xf_argbeg != NULL)
			free(xf_argbeg);
		xf_argbeg = newav;
	}
	cp = (char *)(spec + 1);	/* use memory allocated above */
	for (i = 0; i < ac; i++)
		if (!strcmp(av[i], "-a")) {
			xf_av(spec)[i++] = "-i";
			xf_av(spec)[i] = strcpy(
					spec->xarr->aarg[spec->xarr->ndim].arg,
					"0");
			spec->xarr->aarg[spec->xarr->ndim].i = 0;
			spec->xarr->aarg[spec->xarr->ndim++].n = atoi(av[i]);
		} else {
			xf_av(spec)[i] = strcpy(cp, av[i]);
			cp += strlen(av[i]) + 1;
		}
	return(spec);
}


void
free_xf(XF_SPEC *spec)			/* free a transform */
{
	if (spec->xarr != NULL)
		free(spec->xarr);
	free(spec);
}


int
xf_aname(struct xf_array *ap)			/* put out name for this instance */
{
	static char	oname[10*XF_MAXDIM];
	static char	*oav[3] = {mg_ename[MG_E_OBJECT], oname};
	register int	i;
	register char	*cp1, *cp2;

	if (ap == NULL)
		return(mg_handle(MG_E_OBJECT, 1, oav));
	cp1 = oname;
	*cp1 = 'a';
	for (i = 0; i < ap->ndim; i++) {
		for (cp2 = ap->aarg[i].arg; *cp2; )
			*++cp1 = *cp2++;
		*++cp1 = '.';
	}
	*cp1 = '\0';
	return(mg_handle(MG_E_OBJECT, 2, oav));
}


long
comp_xfid(MAT4 xfm)			/* compute unique ID from matrix */
{
	static char	shifttab[64] = { 15, 5, 11, 5, 6, 3,
				9, 15, 13, 2, 13, 5, 2, 12, 14, 11,
				11, 12, 12, 3, 2, 11, 8, 12, 1, 12,
				5, 4, 15, 9, 14, 5, 13, 14, 2, 10,
				10, 14, 12, 3, 5, 5, 14, 6, 12, 11,
				13, 9, 12, 8, 1, 6, 5, 12, 7, 13,
				15, 8, 9, 2, 6, 11, 9, 11 };
	register int	i;
	register long	xid;

	xid = 0;			/* compute unique transform id */
	for (i = 0; i < sizeof(MAT4)/sizeof(unsigned short); i++)
		xid ^= (long)(((unsigned short *)xfm)[i]) << shifttab[i&63];
	return(xid);
}


void
xf_clear(void)			/* clear transform stack */
{
	register XF_SPEC	*spec;

	if (xf_argbeg != NULL) {
		free(xf_argbeg);
		xf_argbeg = xf_argend = NULL;
	}
	while ((spec = xf_context) != NULL) {
		xf_context = spec->prev;
		free_xf(spec);
	}
}


void
xf_xfmpoint(FVECT v1, FVECT v2)	/* transform a point by the current matrix */
{
	if (xf_context == NULL) {
		VCOPY(v1, v2);
		return;
	}
	multp3(v1, v2, xf_context->xf.xfm);
}


void
xf_xfmvect(FVECT v1, FVECT v2)	/* transform a vector using current matrix */
{
	if (xf_context == NULL) {
		VCOPY(v1, v2);
		return;
	}
	multv3(v1, v2, xf_context->xf.xfm);
}


void
xf_rotvect(FVECT v1, FVECT v2)	/* rotate a vector using current matrix */
{
	xf_xfmvect(v1, v2);
	if (xf_context == NULL)
		return;
	v1[0] /= xf_context->xf.sca;
	v1[1] /= xf_context->xf.sca;
	v1[2] /= xf_context->xf.sca;
}


double
xf_scale(double d)		/* scale a number by the current transform */
{
	if (xf_context == NULL)
		return(d);
	return(d*xf_context->xf.sca);
}
