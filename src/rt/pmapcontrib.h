/* 
   ==================================================================
   Photon map support for light source contributions

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts   
   ==================================================================
   
   $Id: pmapcontrib.h,v 2.1 2015/02/24 19:39:26 greg Exp $
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
