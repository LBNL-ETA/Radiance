#ifndef lint
static const char	RCSid[] = "$Id: sm_sets.c,v 3.7 2003/04/23 00:52:34 greg Exp $";
#endif
/*
 * Quadtree-specific set operations.
 */

#include "standard.h"

#include "sm_qtree.h"


#define QTSETIBLK	2048     	/* set index allocation block size */
#define QTELEMBLK	8		/* set element allocation block */
#define QTELEMMOD(ne)	((ne)&7)	/* (ne) % QTELEMBLK */

#define QTONTHRESH(ne)	(!QTELEMMOD((ne)+1))
#define QTNODESIZ(ne)	((ne) + QTELEMBLK - QTELEMMOD(ne))

OBJECT	**qtsettab= NULL;		/* quadtree leaf node table */
QUADTREE  qtnumsets;			/* number of used set indices */
static int  qtfreesets = EMPTY;		/* free set index list */


QUADTREE
qtnewleaf(oset)			/* return new leaf node for object set */
OBJECT  *oset;
{
  register QUADTREE  osi;

  if (*oset <= 0)
    return(EMPTY);		/* should be error? */
  if (qtfreesets != EMPTY) {
    osi = qtfreesets;
    qtfreesets = (int)qtsettab[osi];
  } else if ((osi = qtnumsets++) % QTSETIBLK == 0) {
    qtsettab = (OBJECT **)realloc((void *)qtsettab,
				  (unsigned)(osi+QTSETIBLK)*sizeof(OBJECT *));
    if (qtsettab == NULL)
      goto memerr;
  }
  qtsettab[osi] = (OBJECT *)malloc(QTNODESIZ(*oset)*sizeof(OBJECT));
  if (qtsettab[osi] == NULL)
    goto memerr;
  setcopy(qtsettab[osi], oset);
  return(QT_INDEX(osi));
memerr:
  error(SYSTEM, "out of memory in qtnewleaf");
}


QUADTREE
qtdelelem(qt, id)		/* delete element from leaf node */
QUADTREE  qt;
OBJECT  id;
{
	register QUADTREE  lf;

#ifdef DEBUG
	if(id < 0)
		eputs("qtdelelem():Invalid triangle id\n");
	if (!QT_IS_LEAF(qt) || (lf = QT_SET_INDEX(qt)) >= qtnumsets)
		error(CONSISTENCY, "empty/bad node in qtdelelem");
	if (!inset(qtsettab[lf], id))
		error(CONSISTENCY, "id not in leaf in qtdelelem");
#else
	lf = QT_SET_INDEX(qt);
#endif
	if (qtsettab[lf][0] <= 1) {		/* blow leaf away */
		free((void *)qtsettab[lf]);
		qtsettab[lf] = (OBJECT *)qtfreesets;
		qtfreesets = lf;
		return(EMPTY);
	}
	deletelem(qtsettab[lf], id);
	if (QTONTHRESH(qtsettab[lf][0]))
		qtsettab[lf] = (OBJECT *)realloc((void *)qtsettab[lf],
				QTNODESIZ(qtsettab[lf][0])*sizeof(OBJECT));
	return(qt);
}


QUADTREE
qtaddelem(qt, id)		/* add element to leaf node */
QUADTREE  qt;
OBJECT  id;
{
	OBJECT	newset[2];
	register QUADTREE  lf;

#ifdef DEBUG
	if(id < 0)
		eputs("qtaddelem():Invalid triangle id\n");
#endif
	if (QT_IS_EMPTY(qt)) {		/* create new leaf */
		newset[0] = 1; newset[1] = id;
		return(qtnewleaf(newset));
	}				/* else add element */
#ifdef DEBUG
	if (!QT_IS_LEAF(qt) || (lf = QT_SET_INDEX(qt)) >= qtnumsets)
		error(CONSISTENCY, "bad node in qtaddelem");
	if (inset(qtsettab[lf], id))
		error(CONSISTENCY, "id already in leaf in qtaddelem");
#else
	lf = QT_SET_INDEX(qt);
#endif
	if (QTONTHRESH(qtsettab[lf][0])) {
		qtsettab[lf] = (OBJECT *)realloc((void *)qtsettab[lf],
				QTNODESIZ(qtsettab[lf][0]+1)*sizeof(OBJECT));
		if (qtsettab[lf] == NULL)
			error(SYSTEM, "out of memory in qtaddelem");
	}
	insertelem(qtsettab[lf], id);
	return(qt);
}


#ifdef DEBUG
OBJECT *
qtqueryset(qt)			/* return object set for leaf node */
QUADTREE  qt;
{
	register QUADTREE  lf;

	if (!QT_IS_LEAF(qt) || (lf = QT_SET_INDEX(qt)) >= qtnumsets)
		error(CONSISTENCY, "bad node in qtqueryset");
	return(qtsettab[lf]);
}
#endif


qtfreeleaf(qt)			/* free set and leaf node */
QUADTREE  qt;
{
	register QUADTREE	osi;

	if (!QT_IS_LEAF(qt))
		return;
	osi = QT_SET_INDEX(qt);
	if (osi >= qtnumsets)
		return;
	free((void *)qtsettab[osi]);
	qtsettab[osi] = (OBJECT *)qtfreesets;
	qtfreesets = osi;
}



qtfreeleaves()			/* free ALL sets and leaf nodes */
{
	register int	i;

	while ((i = qtfreesets) != EMPTY) {
		qtfreesets = (int)qtsettab[i];
		qtsettab[i] = NULL;
	}
	for (i = qtnumsets; i--; )
		if (qtsettab[i] != NULL)
			free((void *)qtsettab[i]);
	free((void *)qtsettab);
	qtsettab = NULL;
	qtnumsets = 0;
}


check_set(os, cs)                /* modify checked set and set to check */
register OBJECT  *os;                   /* os' = os - cs */
register OBJECT  *cs;                   /* cs' = cs + os */
{
    OBJECT  cset[QT_MAXCSET+QT_MAXSET+1];
    register int  i, j;
    int  k;
                                        /* copy os in place, cset <- cs */
    cset[0] = 0;
    k = 0;
    for (i = j = 1; i <= os[0]; i++) {
	while (j <= cs[0] && cs[j] < os[i])
	   cset[++cset[0]] = cs[j++];
	if (j > cs[0] || os[i] != cs[j]) {      /* object to check */
	    os[++k] = os[i];
	    cset[++cset[0]] = os[i];
	}
    }
    if (!(os[0] = k))               /* new "to check" set size */
       return;                 /* special case */
    while (j <= cs[0])              /* get the rest of cs */
       cset[++cset[0]] = cs[j++];
    if (cset[0] > QT_MAXCSET)          /* truncate "checked" set if nec. */
       cset[0] = QT_MAXCSET;
    /* setcopy(cs, cset); */        /* copy cset back to cs */
    os = cset;
    for (i = os[0]; i-- >= 0; )
       *cs++ = *os++;
}









