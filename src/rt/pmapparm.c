/* 
   ==================================================================
   Parameters for photon map generation; used by MKPMAP
   For inclusion in mkpmap.c
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com}
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmapparm.c,v 2.4 2015/05/26 13:31:19 rschregle Exp $
*/


#include "pmapparm.h"
#include "pmapdata.h"
#include "standard.h"
#include <ctype.h>


float pdfSamples = 1000,                /* PDF samples per steradian */
      finalGather = 0.25,               /* fraction of global photons for
                                           irradiance precomputation */
      preDistrib = 0.25,                /* fraction of num photons for
                                           distribution prepass */
      gatherTolerance = 0.5,            /* Photon map lookup tolerance;
                                           lookups returning fewer than this
                                           fraction of minGather/maxGather
                                           are restarted with a larger
                                           search radius */
      maxDistFix = 0;                   /* Static maximum photon search
                                           radius (radius is adaptive if
                                           this is zero) */

#ifdef PMAP_ROI
   /* Region of interest bbox: {xmin, xmax, ymin, ymax, zmin, zmax} */
   float pmapROI [6] = {-FHUGE, FHUGE, -FHUGE, FHUGE, -FHUGE, FHUGE};                                        
#endif                                        

unsigned long photonHeapSizeInc = 1000, /* Photon heap size increment */
              photonMaxBounce = 5000;   /* Runaway photon bounce limit */
              
unsigned photonRepTime = 0,             /* Seconds between reports */
         maxPreDistrib = 4,             /* Max predistrib passes */
         defaultGather = 40;            /* Default num photons for lookup */

/* Per photon map params */
PhotonMapParams pmapParams [NUM_PMAP_TYPES] = {
   {NULL, 0, 0, 0}, {NULL, 0, 0, 0}, {NULL, 0, 0, 0},  {NULL, 0, 0, 0}, 
   {NULL, 0, 0, 0}, {NULL, 0, 0, 0}
};


int setPmapParam (PhotonMap** pm, const PhotonMapParams* parm)
{
   if (parm && parm -> fileName) {
      if (!(*pm = (PhotonMap*)malloc(sizeof(PhotonMap)))) 
         error(INTERNAL, "failed to allocate photon map");
      
      (*pm) -> fileName = parm -> fileName;
      (*pm) -> minGather = parm -> minGather;
      (*pm) -> maxGather = parm -> maxGather;
      (*pm) -> distribTarget = parm -> distribTarget;
      (*pm) -> heapSizeInc = photonHeapSizeInc;      
      (*pm) -> maxDist0 = FHUGE;
      (*pm) -> srcContrib = NULL;

      return 1;
   }
   
   return 0;
}


unsigned long parseMultiplier (const char *num)
{
   unsigned long mult = 1;
   const int strEnd = strlen(num) - 1;
   
   if (strEnd <= 0)
      return 0;
   
   if (!isdigit(num [strEnd]))
      switch (toupper(num [strEnd])) {
         case 'G': mult *= 1000;
         case 'M': mult *= 1000;
         case 'K': mult *= 1000;
                   break;
         default : error(USER, "unknown multiplier");
      }
      
   return (unsigned long)(mult * atof(num));
}
