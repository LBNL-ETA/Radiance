#ifndef lint
static const char RCSid[] = "$Id: m_alias.c,v 2.6 2006/07/12 05:47:05 greg Exp $";
#endif
/*
 * Handler for modifier alias
 */

#include "copyright.h"

#include "ray.h"
#include "otypes.h"
#include "rtotypes.h"
#include "otspecial.h"

/*
 *	If the alias has a single string argument, it's the
 *  name of the target modifier, and we must substitute the
 *  target and its arguments in place of this alias.  The
 *  only difference is that we use our modifier rather than
 *  theirs.
 *	If the alias has no string arguments, then we simply
 *  pass through to our modifier as if we weren't in the
 *  chain at all.
 */

extern int
m_alias(			/* transfer shading to alias target */
	OBJREC	*m,
	RAY	*r
)
{
	OBJECT	aobj;
	OBJREC	*aop;
	OBJREC	arec;
	int	rval;
					/* straight replacement? */
	if (!m->oargs.nsargs)
		return(rayshade(r, m->omod));
					/* else replace alias */
	if (m->oargs.nsargs != 1)
		objerror(m, INTERNAL, "bad # string arguments");
	aop = m;
	aobj = objndx(aop);
	do {				/* follow alias trail */
		if (aop->oargs.nsargs == 1)
			aobj = lastmod(aobj, aop->oargs.sarg[0]);
		else
			aobj = aop->omod;
		if (aobj < 0)
			objerror(aop, USER, "bad reference");
		aop = objptr(aobj);
	} while (aop->otype == MOD_ALIAS);
					/* copy struct */
	arec = *aop;
					/* irradiance hack */
	if (do_irrad && !(r->crtype & ~(PRIMARY|TRANS)) &&
			m->otype != MAT_CLIP &&
			(ofun[arec.otype].flags & (T_M|T_X))) {
		if (irr_ignore(arec.otype)) {
			raytrans(r);
			return(1);
		}
		if (!islight(arec.otype))
			return((*ofun[Lamb.otype].funp)(&Lamb, r));
	}
					/* substitute modifier */
	arec.omod = m->omod;
					/* replacement shader */
	rval = (*ofun[arec.otype].funp)(&arec, r);
					/* save allocated struct */
	if (arec.os != aop->os) {
		if (aop->os != NULL)	/* should never happen */
			free_os(aop);
		aop->os = arec.os;
	}
	return(rval);
}
