/* RCSid $Id: pmapcontrib.h,v 2.6 2024/09/17 16:36:05 greg Exp $ */

/* 
   ==================================================================
   Photon map support for light source contributions

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmapcontrib.h,v 2.6 2024/09/17 16:36:05 greg Exp $
*/

#ifndef PMAPCONTRIB_H
   #define PMAPCONTRIB_H

   #include "pmapdata.h"   

#ifdef __cplusplus
extern "C" {
#endif

   void initPmapContrib (LUTAB *srcContrib, unsigned numSrcContrib);
   /* Set up photon map contributions (interface to rcmain.c) */

   void distribPhotonContrib (PhotonMap *pmap, unsigned numProc);
   /* Emit photons from light sources with tagged contributions, and
    * build photon map */

   void photonContrib (PhotonMap *pmap, RAY *ray, COLOR irrad);
   /* Accumulate light source contributions in pmap -> srcMods from
    * photons, and return cumulative irradiance from density esimate */

#ifdef __cplusplus
}
#endif

#endif
