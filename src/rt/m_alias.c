#ifndef lint
static const char RCSid[] = "$Id: m_alias.c,v 2.4 2004/01/20 22:16:53 greg Exp $";
#endif
/*
 * Handler for modifier alias
 */

#include "copyright.h"
#include "ray.h"
#include "otypes.h"
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
m_alias(m, r)			/* transfer shading to alias target */
OBJREC	*m;
RAY	*r;
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
	aobj = lastmod(objndx(m), m->oargs.sarg[0]);
	if (aobj < 0)
		objerror(m, USER, "bad reference");
	aop = objptr(aobj);
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
