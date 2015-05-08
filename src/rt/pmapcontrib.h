/* 
   ==================================================================
   Photon map support for light source contributions

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmapcontrib.h,v 2.2 2015/05/08 13:20:22 rschregle Exp $
*/

#ifndef PMAP_CONTRIB_H
   #define PMAP_CONTRIB_H

   #include "pmapdata.h"   
   
   void initPmapContrib (LUTAB *srcContrib, unsigned numSrcContrib);
   /* Set up photon map contributions (interface to rcmain.c) */

   void distribPhotonContrib (PhotonMap *pmap);
   /* Emit photons from light sources with tagged contributions, and
    * build photon map */

   void photonContrib (PhotonMap *pmap, RAY *ray, COLOR irrad);
   /* Accumulate light source contributions in pmap -> srcMods from
    * photons, and return cumulative irradiance from density esimate */

#endif
