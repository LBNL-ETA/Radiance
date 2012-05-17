#ifndef lint
static const char	RCSid[] = "$Id: mgf_context.c,v 3.2 2012/05/17 17:10:23 greg Exp $";
#endif
/*
 * Context handlers
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mgf_parser.h"
#include "lookup.h"

				/* default context values */
static C_MATERIAL	c_dfmaterial = C_DEFMATERIAL;
static C_VERTEX		c_dfvertex = C_DEFVERTEX;

				/* the unnamed contexts */
static C_COLOR		c_uncolor = C_DEFCOLOR;
static C_MATERIAL	c_unmaterial = C_DEFMATERIAL;
static C_VERTEX		c_unvertex = C_DEFVERTEX;

				/* the current contexts */
C_COLOR		*c_ccolor = &c_uncolor;
char		*c_ccname = NULL;
C_MATERIAL	*c_cmaterial = &c_unmaterial;
char		*c_cmname = NULL;
C_VERTEX	*c_cvertex = &c_unvertex;
char		*c_cvname = NULL;

static LUTAB	clr_tab = LU_SINIT(free,free);	/* color lookup table */
static LUTAB	mat_tab = LU_SINIT(free,free);	/* material lookup table */
static LUTAB	vtx_tab = LU_SINIT(free,free);	/* vertex lookup table */

static int	setspectrum(C_COLOR *, int ac, char **av);


int
c_hcolor(int ac, char **av)		/* handle color entity */
{
	double	w, wsum;
	int	i;
	LUENT	*lp;

	switch (mg_entity(av[0])) {
	case MG_E_COLOR:	/* get/set color context */
		if (ac > 4)
			return(MG_EARGC);
		if (ac == 1) {		/* set unnamed color context */
			c_uncolor = c_dfcolor;
			c_ccolor = &c_uncolor;
			c_ccname = NULL;
			return(MG_OK);
		}
		if (!isname(av[1]))
			return(MG_EILL);
		lp = lu_find(&clr_tab, av[1]);	/* lookup context */
		if (lp == NULL)
			return(MG_EMEM);
		c_ccname = lp->key;
		c_ccolor = (C_COLOR *)lp->data;
		if (ac == 2) {		/* reestablish previous context */
			if (c_ccolor == NULL)
				return(MG_EUNDEF);
			return(MG_OK);
		}
		if (av[2][0] != '=' || av[2][1])
			return(MG_ETYPE);
		if (c_ccolor == NULL) {	/* create new color context */
			lp->key = (char *)malloc(strlen(av[1])+1);
			if (lp->key == NULL)
				return(MG_EMEM);
			strcpy(lp->key, av[1]);
			lp->data = (char *)malloc(sizeof(C_COLOR));
			if (lp->data == NULL)
				return(MG_EMEM);
			c_ccname = lp->key;
			c_ccolor = (C_COLOR *)lp->data;
			c_ccolor->clock = 0;
			c_ccolor->client_data = NULL;
		}
		i = c_ccolor->clock;
		if (ac == 3) {		/* use default template */
			*c_ccolor = c_dfcolor;
			c_ccolor->clock = i + 1;
			return(MG_OK);
		}
		lp = lu_find(&clr_tab, av[3]);	/* lookup template */
		if (lp == NULL)
			return(MG_EMEM);
		if (lp->data == NULL)
			return(MG_EUNDEF);
		*c_ccolor = *(C_COLOR *)lp->data;
		c_ccolor->clock = i + 1;
		return(MG_OK);
	case MG_E_CXY:		/* assign CIE XY value */
		if (ac != 3)
			return(MG_EARGC);
		if (!isflt(av[1]) | !isflt(av[2]))
			return(MG_ETYPE);
		c_ccolor->cx = atof(av[1]);
		c_ccolor->cy = atof(av[2]);
		c_ccolor->flags = C_CDXY|C_CSXY;
		if ((c_ccolor->cx < 0.) | (c_ccolor->cy < 0.) |
				(c_ccolor->cx + c_ccolor->cy > 1.))
			return(MG_EILL);
		c_ccolor->clock++;
		return(MG_OK);
	case MG_E_CSPEC:	/* assign spectral values */
		return(setspectrum(c_ccolor, ac-1, av+1));
	case MG_E_CCT:		/* assign black body spectrum */
		if (ac != 2)
			return(MG_EARGC);
		if (!isflt(av[1]))
			return(MG_ETYPE);
		if (!c_bbtemp(c_ccolor, atof(av[1])))
			return(MG_EILL);
		c_ccolor->clock++;
		return(MG_OK);
	case MG_E_CMIX:		/* mix colors */
		if (ac < 5 || (ac-1)%2)
			return(MG_EARGC);
		if (!isflt(av[1]))
			return(MG_ETYPE);
		wsum = atof(av[1]);
		if ((lp = lu_find(&clr_tab, av[2])) == NULL)
			return(MG_EMEM);
		if (lp->data == NULL)
			return(MG_EUNDEF);
		*c_ccolor = *(C_COLOR *)lp->data;
		for (i = 3; i < ac; i += 2) {
			if (!isflt(av[i]))
				return(MG_ETYPE);
			w = atof(av[i]);
			if ((lp = lu_find(&clr_tab, av[i+1])) == NULL)
				return(MG_EMEM);
			if (lp->data == NULL)
				return(MG_EUNDEF);
			c_cmix(c_ccolor, wsum, c_ccolor,
					w, (C_COLOR *)lp->data);
			wsum += w;
		}
		if (wsum <= 0.)
			return(MG_EILL);
		c_ccolor->clock++;
		return(MG_OK);
	}
	return(MG_EUNK);
}


int
c_hmaterial(int ac, char **av)		/* handle material entity */
{
	int	i;
	LUENT	*lp;

	switch (mg_entity(av[0])) {
	case MG_E_MATERIAL:	/* get/set material context */
		if (ac > 4)
			return(MG_EARGC);
		if (ac == 1) {		/* set unnamed material context */
			c_unmaterial = c_dfmaterial;
			c_cmaterial = &c_unmaterial;
			c_cmname = NULL;
			return(MG_OK);
		}
		if (!isname(av[1]))
			return(MG_EILL);
		lp = lu_find(&mat_tab, av[1]);	/* lookup context */
		if (lp == NULL)
			return(MG_EMEM);
		c_cmname = lp->key;
		c_cmaterial = (C_MATERIAL *)lp->data;
		if (ac == 2) {		/* reestablish previous context */
			if (c_cmaterial == NULL)
				return(MG_EUNDEF);
			return(MG_OK);
		}
		if (av[2][0] != '=' || av[2][1])
			return(MG_ETYPE);
		if (c_cmaterial == NULL) {	/* create new material */
			lp->key = (char *)malloc(strlen(av[1])+1);
			if (lp->key == NULL)
				return(MG_EMEM);
			strcpy(lp->key, av[1]);
			lp->data = (char *)malloc(sizeof(C_MATERIAL));
			if (lp->data == NULL)
				return(MG_EMEM);
			c_cmname = lp->key;
			c_cmaterial = (C_MATERIAL *)lp->data;
			c_cmaterial->clock = 0;
			c_cmaterial->client_data = NULL;
		}
		i = c_cmaterial->clock;
		if (ac == 3) {		/* use default template */
			*c_cmaterial = c_dfmaterial;
			c_cmaterial->clock = i + 1;
			return(MG_OK);
		}
		lp = lu_find(&mat_tab, av[3]);	/* lookup template */
		if (lp == NULL)
			return(MG_EMEM);
		if (lp->data == NULL)
			return(MG_EUNDEF);
		*c_cmaterial = *(C_MATERIAL *)lp->data;
		c_cmaterial->clock = i + 1;
		return(MG_OK);
	case MG_E_IR:		/* set index of refraction */
		if (ac != 3)
			return(MG_EARGC);
		if (!isflt(av[1]) | !isflt(av[2]))
			return(MG_ETYPE);
		c_cmaterial->nr = atof(av[1]);
		c_cmaterial->ni = atof(av[2]);
		if (c_cmaterial->nr <= FTINY)
			return(MG_EILL);
		c_cmaterial->clock++;
		return(MG_OK);
	case MG_E_RD:		/* set diffuse reflectance */
		if (ac != 2)
			return(MG_EARGC);
		if (!isflt(av[1]))
			return(MG_ETYPE);
		c_cmaterial->rd = atof(av[1]);
		if ((c_cmaterial->rd < 0.) | (c_cmaterial->rd > 1.))
			return(MG_EILL);
		c_cmaterial->rd_c = *c_ccolor;
		c_cmaterial->clock++;
		return(MG_OK);
	case MG_E_ED:		/* set diffuse emittance */
		if (ac != 2)
			return(MG_EARGC);
		if (!isflt(av[1]))
			return(MG_ETYPE);
		c_cmaterial->ed = atof(av[1]);
		if (c_cmaterial->ed < 0.)
			return(MG_EILL);
		c_cmaterial->ed_c = *c_ccolor;
		c_cmaterial->clock++;
		return(MG_OK);
	case MG_E_TD:		/* set diffuse transmittance */
		if (ac != 2)
			return(MG_EARGC);
		if (!isflt(av[1]))
			return(MG_ETYPE);
		c_cmaterial->td = atof(av[1]);
		if ((c_cmaterial->td < 0.) | (c_cmaterial->td > 1.))
			return(MG_EILL);
		c_cmaterial->td_c = *c_ccolor;
		c_cmaterial->clock++;
		return(MG_OK);
	case MG_E_RS:		/* set specular reflectance */
		if (ac != 3)
			return(MG_EARGC);
		if (!isflt(av[1]) | !isflt(av[2]))
			return(MG_ETYPE);
		c_cmaterial->rs = atof(av[1]);
		c_cmaterial->rs_a = atof(av[2]);
		if ((c_cmaterial->rs < 0.) | (c_cmaterial->rs > 1.) |
				(c_cmaterial->rs_a < 0.))
			return(MG_EILL);
		c_cmaterial->rs_c = *c_ccolor;
		c_cmaterial->clock++;
		return(MG_OK);
	case MG_E_TS:		/* set specular transmittance */
		if (ac != 3)
			return(MG_EARGC);
		if (!isflt(av[1]) | !isflt(av[2]))
			return(MG_ETYPE);
		c_cmaterial->ts = atof(av[1]);
		c_cmaterial->ts_a = atof(av[2]);
		if ((c_cmaterial->ts < 0.) | (c_cmaterial->ts > 1.) |
				(c_cmaterial->ts_a < 0.))
			return(MG_EILL);
		c_cmaterial->ts_c = *c_ccolor;
		c_cmaterial->clock++;
		return(MG_OK);
	case MG_E_SIDES:	/* set number of sides */
		if (ac != 2)
			return(MG_EARGC);
		if (!isint(av[1]))
			return(MG_ETYPE);
		i = atoi(av[1]);
		if (i == 1)
			c_cmaterial->sided = 1;
		else if (i == 2)
			c_cmaterial->sided = 0;
		else
			return(MG_EILL);
		c_cmaterial->clock++;
		return(MG_OK);
	}
	return(MG_EUNK);
}


int
c_hvertex(int ac, char **av)		/* handle a vertex entity */
{
	int	i;
	LUENT	*lp;

	switch (mg_entity(av[0])) {
	case MG_E_VERTEX:	/* get/set vertex context */
		if (ac > 4)
			return(MG_EARGC);
		if (ac == 1) {		/* set unnamed vertex context */
			c_unvertex = c_dfvertex;
			c_cvertex = &c_unvertex;
			c_cvname = NULL;
			return(MG_OK);
		}
		if (!isname(av[1]))
			return(MG_EILL);
		lp = lu_find(&vtx_tab, av[1]);	/* lookup context */
		if (lp == NULL)
			return(MG_EMEM);
		c_cvname = lp->key;
		c_cvertex = (C_VERTEX *)lp->data;
		if (ac == 2) {		/* reestablish previous context */
			if (c_cvertex == NULL)
				return(MG_EUNDEF);
			return(MG_OK);
		}
		if (av[2][0] != '=' || av[2][1])
			return(MG_ETYPE);
		if (c_cvertex == NULL) {	/* create new vertex context */
			lp->key = (char *)malloc(strlen(av[1])+1);
			if (lp->key == NULL)
				return(MG_EMEM);
			strcpy(lp->key, av[1]);
			lp->data = (char *)malloc(sizeof(C_VERTEX));
			if (lp->data == NULL)
				return(MG_EMEM);
			c_cvname = lp->key;
			c_cvertex = (C_VERTEX *)lp->data;
			c_cvertex->clock = 0;
			c_cvertex->client_data = NULL;
		}
		i = c_cvertex->clock;
		if (ac == 3) {		/* use default template */
			*c_cvertex = c_dfvertex;
			c_cvertex->clock = i + 1;
			return(MG_OK);
		}
		lp = lu_find(&vtx_tab, av[3]);	/* lookup template */
		if (lp == NULL)
			return(MG_EMEM);
		if (lp->data == NULL)
			return(MG_EUNDEF);
		*c_cvertex = *(C_VERTEX *)lp->data;
		c_cvertex->clock = i + 1;
		return(MG_OK);
	case MG_E_POINT:	/* set point */
		if (ac != 4)
			return(MG_EARGC);
		if (!isflt(av[1]) | !isflt(av[2]) | !isflt(av[3]))
			return(MG_ETYPE);
		c_cvertex->p[0] = atof(av[1]);
		c_cvertex->p[1] = atof(av[2]);
		c_cvertex->p[2] = atof(av[3]);
		c_cvertex->clock++;
		return(MG_OK);
	case MG_E_NORMAL:	/* set normal */
		if (ac != 4)
			return(MG_EARGC);
		if (!isflt(av[1]) | !isflt(av[2]) | !isflt(av[3]))
			return(MG_ETYPE);
		c_cvertex->n[0] = atof(av[1]);
		c_cvertex->n[1] = atof(av[2]);
		c_cvertex->n[2] = atof(av[3]);
		(void)normalize(c_cvertex->n);
		c_cvertex->clock++;
		return(MG_OK);
	}
	return(MG_EUNK);
}


void
c_clearall(void)		/* empty context tables */
{
	c_uncolor = c_dfcolor;
	c_ccolor = &c_uncolor;
	c_ccname = NULL;
	lu_done(&clr_tab);
	c_unmaterial = c_dfmaterial;
	c_cmaterial = &c_unmaterial;
	c_cmname = NULL;
	lu_done(&mat_tab);
	c_unvertex = c_dfvertex;
	c_cvertex = &c_unvertex;
	c_cvname = NULL;
	lu_done(&vtx_tab);
}


C_MATERIAL *
c_getmaterial(char *name)	/* get a named material */
{
	LUENT	*lp;

	if ((lp = lu_find(&mat_tab, name)) == NULL)
		return(NULL);
	return((C_MATERIAL *)lp->data);
}


C_VERTEX *
c_getvert(char *name)			/* get a named vertex */
{
	LUENT	*lp;

	if ((lp = lu_find(&vtx_tab, name)) == NULL)
		return(NULL);
	return((C_VERTEX *)lp->data);
}


C_COLOR *
c_getcolor(char *name)		/* get a named color */
{
	LUENT	*lp;

	if ((lp = lu_find(&clr_tab, name)) == NULL)
		return(NULL);
	return((C_COLOR *)lp->data);
}


static int
setspectrum(			/* convert a spectrum */
	C_COLOR	*clr,
	int	ac,
	char	**av
)
{
	double	wlmin, wlmax;
	int	i;
	float	*va;
					/* get bounds */
	if (ac < 4)
		return(MG_EARGC);
	if (!isflt(av[0]) || !isflt(av[1]))
		return(MG_ETYPE);
	wlmin = atof(av[0]);
	wlmax = atof(av[1]);
	ac -= 2; av += 2;
	va = (float *)malloc(sizeof(float)*ac);
	if (va == NULL)
		return(MG_EMEM);
	for (i = ac; i--; ) {
		if (!isflt(av[i]))
			return(MG_ETYPE);
		va[i] = atof(av[i]);
	}
	if (c_sset(clr, wlmin, wlmax, va, ac) <= FTINY) {
		free(va);
		return(MG_EILL);
	}
	free(va);
	clr->clock++;
	return(MG_OK);
}
