/* RCSid $Id: pmapparm.h,v 2.10 2021/04/14 11:26:25 rschregle Exp $ */

/* 
   ======================================================================
   Parameters for photon map generation and rendering; used by mkpmap
   and rpict/rvu/rtrace.   
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com}
   (c) Fraunhofer Institute for Solar Energy Systems,
       supported by the German Research Foundation 
       (DFG LU-204/10-2, "Fassadenintegrierte Regelsysteme" (FARESYS)) 
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation 
       (SNSF #147053, "Daylight Redirecting Components")
   ======================================================================
   
   $Id: pmapparm.h,v 2.10 2021/04/14 11:26:25 rschregle Exp $
*/


#ifndef PMAPPARAMS_H
   #define PMAPPARAMS_H

   #include "pmaptype.h"
   
      
   /* Struct for passing params per photon map from rpict/rtrace/rvu */
   typedef struct {
      char *fileName;                /* Photon map file */
      unsigned minGather, maxGather; /* Num photons to gather */
      unsigned long distribTarget;   /* Num photons to store */
   } PhotonMapParams;
   
   /* Region of interest */
   typedef struct {
      /* siz [1], siz [2] <= 0 --> sphere, else rectangle */
      float pos [3], siz [3];
   } PhotonMapROI;
   
   #define PMAP_ROI_ISSPHERE(roi) ((roi)->siz[1] <= 0 && (roi)->siz[2] <= 0)
   #define PMAP_ROI_SETSPHERE(roi) ((roi)->siz[1] = (roi)->siz[2] = -1)


   extern PhotonMapParams pmapParams [NUM_PMAP_TYPES];
   
   /* Macros for type specific photon map parameters */
   #define globalPmapParams	(pmapParams [PMAP_TYPE_GLOBAL])
   #define preCompPmapParams  (pmapParams [PMAP_TYPE_PRECOMP])
   #define causticPmapParams	(pmapParams [PMAP_TYPE_CAUSTIC])
   #define volumePmapParams	(pmapParams [PMAP_TYPE_VOLUME])
   #define directPmapParams	(pmapParams [PMAP_TYPE_DIRECT])
   #define contribPmapParams  (pmapParams [PMAP_TYPE_CONTRIB])
   
   
   extern float         pdfSamples, preDistrib, finalGather,
                        gatherTolerance, maxDistFix, pmapMaxDist,
                        photonMaxDist;
   extern unsigned long photonHeapSizeInc, photonMaxBounce; 
   extern unsigned      photonRepTime, maxPreDistrib, defaultGather,
                        verbose;

   extern unsigned      pmapNumROI;
   extern PhotonMapROI  *pmapROI;

#ifdef PMAP_OOC
   extern float         pmapCachePageSize;
   extern unsigned long pmapCacheSize;
#endif   


   struct PhotonMap;
   
   int setPmapParam (struct PhotonMap **pm, const PhotonMapParams *parm);
   /* Allocate photon map and set its parameters from parm */
   
   unsigned long parseMultiplier (const char *num);
   /* Evaluate numeric parameter string with optional multiplier suffix
      (G = 10^9, M = 10^6, K = 10^3). Returns 0 if parsing fails. */
#endif
