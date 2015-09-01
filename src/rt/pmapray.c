#ifndef lint
static const char RCSid[] = "$Id: pmapray.c,v 2.5 2015/09/01 16:27:53 greg Exp $";
#endif
/* 
   ==================================================================
   Photon map interface to RADIANCE raycalls

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================   
   
*/


#include "pmapray.h"
#include "pmap.h"


void ray_init_pmap ()
/* Interface to ray_init(); init & load photon maps */
{
   loadPmaps(photonMaps, pmapParams);
}


void ray_done_pmap ()
/* Interface to ray_done(); free photon maps */
{
   cleanUpPmaps(photonMaps);
}


void ray_save_pmap (RAYPARAMS *rp)
/* Interface to ray_save(); save photon map params */
{
   unsigned t;
   
   for (t = 0; t < NUM_PMAP_TYPES; t++) {
      if (pmapParams [t].fileName)
         rp -> pmapParams [t].fileName = savqstr(pmapParams [t].fileName);
        
      else rp -> pmapParams [t].fileName = NULL;
   
      rp -> pmapParams [t].minGather = pmapParams [t].minGather;
      rp -> pmapParams [t].maxGather = pmapParams [t].maxGather;
   }
}


void ray_restore_pmap (RAYPARAMS *rp)
/* Interface to ray_restore(); restore photon mapping params */
{
   unsigned t;
   
   for (t = 0; t < NUM_PMAP_TYPES; t++) {
      pmapParams [t].fileName = rp -> pmapParams [t].fileName;
      pmapParams [t].minGather = rp -> pmapParams [t].minGather;
      pmapParams [t].maxGather = rp -> pmapParams [t].maxGather;
   }
}
