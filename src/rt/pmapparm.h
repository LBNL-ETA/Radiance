/* 
   ==================================================================
   Parameters common to all photon map modules. 
   For inclusion in main program.
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com}
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts
   ==================================================================
   
   $Id: pmapparm.h,v 2.1 2015/02/24 19:39:27 greg Exp $
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


   extern PhotonMapParams pmapParams [NUM_PMAP_TYPES];
   
   /* Macros for type specific photon map parameters */
   #define globalPmapParams	(pmapParams [PMAP_TYPE_GLOBAL])
   #define preCompPmapParams  (pmapParams [PMAP_TYPE_PRECOMP])
   #define causticPmapParams	(pmapParams [PMAP_TYPE_CAUSTIC])
   #define volumePmapParams	(pmapParams [PMAP_TYPE_VOLUME])
   #define directPmapParams	(pmapParams [PMAP_TYPE_DIRECT])
   #define contribPmapParams  (pmapParams [PMAP_TYPE_CONTRIB])
   
   
   extern float pdfSamples, preDistrib, finalGather, gatherTolerance,
                maxDistCoeff;
   extern unsigned long photonHeapSizeInc, photonMaxBounce;
   extern unsigned photonRepTime, maxPreDistrib, defaultGather;
#ifdef PMAP_ROI                
   extern float pmapROI [6];
#endif   

   struct PhotonMap;
   
   int setPmapParam (struct PhotonMap **pm, const PhotonMapParams *parm);
   /* Allocate photon map and set its parameters from parm */
   
   unsigned long parseMultiplier (const char *num);
   /* Evaluate numeric parameter string with optional multiplier suffix
      (G = 10^9, M = 10^6, K = 10^3). Returns 0 if parsing fails. */
#endif
