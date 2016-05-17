/* RCSid $Id: pmapcontrib.h,v 2.5 2016/05/17 17:39:47 rschregle Exp $ */

/* 
   ==================================================================
   Photon map support for light source contributions

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmapcontrib.h,v 2.5 2016/05/17 17:39:47 rschregle Exp $
*/

#ifndef PMAPCONTRIB_H
   #define PMAPCONTRIB_H

   #include "pmapdata.h"   
   
   void initPmapContrib (LUTAB *srcContrib, unsigned numSrcContrib);
   /* Set up photon map contributions (interface to rcmain.c) */

   void distribPhotonContrib (PhotonMap *pmap, unsigned numProc);
   /* Emit photons from light sources with tagged contributions, and
    * build photon map */

   void photonContrib (PhotonMap *pmap, RAY *ray, COLOR irrad);
   /* Accumulate light source contributions in pmap -> srcMods from
    * photons, and return cumulative irradiance from density esimate */

#endif
