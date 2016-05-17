#ifndef lint
static const char RCSid[] = "$Id: pmapamb.c,v 2.8 2016/05/17 17:39:47 rschregle Exp $";
#endif

/* 
   ==================================================================
   Photon map interface to RADIANCE ambient calculation

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================   
   
   $Id: pmapamb.c,v 2.8 2016/05/17 17:39:47 rschregle Exp $
*/


#include "pmapamb.h"
#include "pmap.h"


int ambPmap (COLOR aval, RAY *r, int rdepth)
/* Factor irradiance from photon map into ambient coefficient aval;
 * return 1 on success, else 0 (with aval unmodified) */
{
   COLOR rcoef, photonIrrad;
   /* Handle precedence in case of multiple photon maps:
    * contrib > precomp > global */
   PhotonMap *pmap = contribPhotonMapping ? contribPmap 
                                          : preCompPmap ? preCompPmap 
                                                        : globalPmap;
      
   /* Get photon irradiance either via 1 ambient bounce (final
    * gather) if ambounce > 0, or directly if ambounce < 0. */
   if (pmap && (rdepth || ambounce < 0)) {
      /* Temporarily factor ambient value into ray coefficient
       * (required for contribution photon map) */
      copycolor(rcoef, r -> rcoef);
      multcolor(r -> rcoef, aval);
      
      /* Get photon irradiance via callback */
      pmap -> lookupFlags = 0;
      (pmap -> lookup)(pmap, r, photonIrrad);
         
      /* Factor irradiance into ambient value and restore ray coeficient */
      multcolor(aval, photonIrrad);
      copycolor(r -> rcoef, rcoef);
      
      return 1;
   }
      
   return 0;
}


int ambPmapCaustic (COLOR aval, RAY *r, int rdepth)
/* Factor specular-diffuse (caustic) irradiance from photon map into ambient
 * coeff aval; return 1 if successful, else 0 (with aval set to zero) */
{
   COLOR rcoef, photonIrrad;
   /* Handle precedence in case of multiple photon maps: contrib > caustic */
   PhotonMap *pmap = contribPhotonMapping ? contribPmap : causticPmap;

   /* Get caustic photon density estimate at primary rays or when
    * filling in ambient rays that have no global photon map to use */
   if (pmap && (!rdepth || !globalPmap & !contribPmap & !preCompPmap)) {
      /* Temporarily factor ambient value into ray coefficient
       * (required for contribution photon map) */
      copycolor(rcoef, r -> rcoef);      
      multcolor(r -> rcoef, aval);
      
      /* Set caustic flag and get photon irradiance via callback */
      pmap -> lookupCaustic = 1;
      (pmap -> lookup)(pmap, r, photonIrrad);
      
      /* Factor irradiance into ambient value and restore ray coeficient */
      multcolor(aval, photonIrrad);
      copycolor(r -> rcoef, rcoef);
      
      return 1;
   }

   setcolor(aval, 0, 0, 0);
   
   return 0;
}
