#ifndef lint
static const char RCSid[] = "$Id: pmapray.c,v 2.9 2024/03/23 16:40:52 greg Exp $";
#endif

/* 
   ==================================================================
   Photon map interface to RADIANCE raycalls

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================   
   
   $Id: pmapray.c,v 2.9 2024/03/23 16:40:52 greg Exp $
*/


#include "ray.h"
#include "pmapray.h"
#include "pmap.h"


void ray_init_pmap ()
/* Interface to ray_init(); init & load photon maps */
{
   loadPmaps(photonMaps, pmapParams);

   if (NCSAMP > 3 && photonMapping)
	error(USER, "hyperspectral photon-mapping currently unsupported");
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
      rp -> pmapParams [t].distribTarget = pmapParams [t].distribTarget;
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
      pmapParams [t].distribTarget = rp -> pmapParams [t].distribTarget;
   }
}

void ray_defaults_pmap (RAYPARAMS *rp)
/* Interface to ray_defaults(); set photon mapping defaults */
{
   unsigned t;
   
   for (t = 0; t < NUM_PMAP_TYPES; t++) {
      rp -> pmapParams [t].fileName = NULL;
      rp -> pmapParams [t].minGather = 0;
      rp -> pmapParams [t].maxGather = 0;
      rp -> pmapParams [t].distribTarget = 0;
   }
}
