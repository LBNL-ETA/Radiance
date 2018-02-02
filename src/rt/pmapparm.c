#ifndef lint
static const char RCSid[] = "$Id: pmapparm.c,v 2.9 2018/02/02 19:47:55 rschregle Exp $";
#endif

/* 
   ======================================================================
   Parameters for photon map generation and rendering; used by mkpmap
   and rpict/rvu/rtrace.
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com}
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: pmapparm.c,v 2.9 2018/02/02 19:47:55 rschregle Exp $
*/


#include "pmapparm.h"
#include "pmapdata.h"
#include <ctype.h>


float pdfSamples        = 1000,        /* PDF samples per steradian */

      finalGather       = 0.25,        /* fraction of global photons for
                                          irradiance precomputation */
                                          
      preDistrib        = 0.25,        /* fraction of num photons for
                                          distribution prepass */
                                          
      gatherTolerance   = 0.5,         /* Photon map lookup tolerance;
                                          lookups returning fewer than this
                                          fraction of minGather/maxGather
                                          are restarted with a larger
                                          search radius */
                                          
      maxDistFix        = 0,           /* Static maximum photon search
                                          radius (radius is adaptive if
                                          this is zero) */
      
      photonMaxDist     = 0;           /* Maximum cumulative distance of
                                          photon path */
                                          
#ifdef PMAP_OOC
float          pmapCachePageSize = 8;     /* OOC cache pagesize as multiple
                                           * of maxGather */
unsigned long  pmapCacheSize     = 1e6;   /* OOC cache size in photons */
#endif


/* Regions of interest */
unsigned pmapNumROI = 0;
PhotonMapROI *pmapROI = NULL;


unsigned verbose = 0;                  /* Verbose console output */
unsigned long photonMaxBounce = 5000;  /* Runaway photon bounce limit */
unsigned photonRepTime        = 0,     /* Seconds between reports */
         maxPreDistrib        = 4,     /* Max predistrib passes */
         defaultGather        = 40;    /* Default num photons for lookup */


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
