/* RCSid $Id: pmap.h,v 2.9 2018/01/24 19:39:05 rschregle Exp $ */

/* 
   ======================================================================
   Photon map main header

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: pmap.h,v 2.9 2018/01/24 19:39:05 rschregle Exp $
*/


#ifndef PMAP_H
   #define PMAP_H

   #ifndef NIX
      #if defined(_WIN32) || defined(_WIN64)
         #define NIX 0
      #else
         #define NIX 1
      #endif
   #endif

   #include "pmapparm.h"
   #include "pmapdata.h"


   #ifndef min
      #define min(a, b)          ((a) < (b) ? (a) : (b))
   #endif
   
   #ifndef max
      #define max(a, b)          ((a) > (b) ? (a) : (b))
   #endif
   
   #define sqr(a)                ((a) * (a))

   /* Average over colour channels */
   #define colorAvg(col) ((col [0] + col [1] + col [2]) / 3)

   /* Macros to test for enabled photon maps */
   #define photonMapping         (globalPmap || preCompPmap || \
                                  causticPmap || contribPmap)
   #define causticPhotonMapping  (causticPmap != NULL)
   #define directPhotonMapping   (directPmap != NULL)
   #define volumePhotonMapping   (volumePmap != NULL)
   /* #define contribPhotonMapping  (contribPmap && contribPmap -> srcContrib)
    */
   #define contribPhotonMapping  (contribPmap)
      

   extern void (*pmapLookup [])(PhotonMap*, RAY*, COLOR);
   /* Photon map lookup functions per type */
   
   void loadPmaps (PhotonMap **pmaps, const PhotonMapParams *params);
   /* Load photon and set their respective parameters, checking timestamps
    * relative to octree for possible staleness */

   void savePmaps (const PhotonMap **pmaps, int argc, char **argv);
   /* Save all defined photon maps with specified command line */

   void cleanUpPmaps (PhotonMap **pmaps);
   /* Trash all photon maps after processing is complete */

   void distribPhotons (PhotonMap **pmaps, unsigned numProc);
   /* Emit photons from light sources and build photon maps for non-NULL
    * entries in photon map array */

   void tracePhoton (RAY*);
   /* Follow photon as it bounces around the scene. Analogon to
    * raytrace(). */

   void photonDensity (PhotonMap*, RAY*, COLOR irrad);
   /* Perform surface density estimate from incoming photon flux at
      ray's intersection point. Returns irradiance from found photons. */

   void photonPreCompDensity (PhotonMap *pmap, RAY *r, COLOR irrad);
   /* Returns precomputed photon density estimate at ray -> rop. */
      
   void volumePhotonDensity (PhotonMap*, RAY*, COLOR);
   /* Perform volume density estimate from incoming photon flux at 
      ray's intersection point. Returns irradiance. */

   void colorNorm (COLOR);
   /* Normalise colour channels to average of 1 */

#endif
