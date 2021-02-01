#ifndef lint
static const char RCSid[] = "$Id: m_alias.c,v 2.11 2021/02/01 16:19:49 greg Exp $";
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

int
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

	aop = m;			/* else look it up */
	aobj = objndx(aop);
	do {				/* follow entire alias trail */
		if (!aop->oargs.nsargs)
			aobj = aop->omod;
		else if (aop->oargs.nsargs == 1)
			aobj = lastmod(aobj, aop->oargs.sarg[0]);
		else
			objerror(aop, INTERNAL, "bad # string arguments");
		if (aobj == OVOID)
			objerror(aop, USER, "bad reference");
		aop = objptr(aobj);
	} while (aop->otype == MOD_ALIAS);
					/* copy struct */
	arec = *aop;
					/* substitute modifier */
	arec.omod = m->omod;
					/* irradiance hack */
	if (do_irrad && !(r->crtype & ~(PRIMARY|TRANS)) && raytirrad(&arec, r))
		return(1);
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
