/* 
   ==================================================================
   In-core kd-tree for photon map

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmapkdt.h,v 1.2 2020/04/08 15:14:21 rschregle Exp $
*/


#ifndef PMAPKDT_H
   #define PMAPKDT_H

/*   #include <stdio.h>
   #include "fvect.h" */
   #include "pmapdata.h"
   

   /* Forward declarations to break dependency loop with pmapdata.h */
   struct PhotonMap;


   /* k-d tree as linear array of photons */
   typedef struct {
      Photon   *nodes;  /* k-d tree as linear array */
   } PhotonKdTree;

   typedef  PhotonKdTree   PhotonStorage;
   typedef  Photon*        PhotonIdx;

   /* Priority queue node for kd-tree lookups */
   typedef struct {                  
      PhotonIdx   idx;        /* Pointer to photon */
      float       dist2;      /* Photon's SQUARED distance to query point */
   } kdT_SearchQueueNode;
   
   typedef struct {
      kdT_SearchQueueNode  *node;
      unsigned             len, tail;
   } kdT_SearchQueue;

   typedef  kdT_SearchQueueNode  PhotonSearchQueueNode;
   typedef  kdT_SearchQueue      PhotonSearchQueue;



   void kdT_Null (PhotonKdTree *kdt);
   /* Initialise kd-tree prior to storing photons */
   
   void kdT_BuildPhotonMap (struct PhotonMap *pmap);
   /* Build a balanced kd-tree pmap -> store from photons in unsorted
    * heapfile pmap -> heap to guarantee logarithmic search times.  The heap
    * is destroyed on return.  */

   int kdT_SavePhotons (const struct PhotonMap *pmap, FILE *out);
   /* Save photons in kd-tree to file. Return -1 on error, else 0 */
   
   int kdT_LoadPhotons (struct PhotonMap *pmap, FILE *in);
   /* Load photons in kd-tree from file. Return -1 on error, else 0 */

   void kdT_InitFindPhotons (struct PhotonMap *pmap);
   /* Initialise NN search queue prior to calling kdT_FindPhotons() */
   
   int kdT_FindPhotons (struct PhotonMap* pmap, const FVECT pos, 
                        const FVECT norm);   
   /* Locate pmap -> squeue.len nearest photons to pos with similar normal
    * (NULL for volume photons) and return in search queue pmap -> squeue,
    * starting with the further photon at pmap -> squeue.node. Return -1
    * if none found, else 0. */

   int kdT_Find1Photon (struct PhotonMap* pmap, const FVECT pos, 
                        const FVECT norm, Photon *photon);
   /* Locate single nearest photon to pos with similar normal. Return -1
    * if none found, else 0. */
   
   int kdT_GetPhoton (const struct PhotonMap *pmap, PhotonIdx idx,
                      Photon *photon);
   /* Retrieve photon referenced by idx from kd-tree and return -1 on
    * error, else 0. */

   Photon *kdT_GetNearestPhoton (const PhotonSearchQueue *squeue, 
                                 PhotonIdx idx);
   /* Retrieve photon from NN search queue after OOC_FindPhotons() */
   
   PhotonIdx kdT_FirstPhoton (const struct PhotonMap* pmap);
   /* Return pointer to first photon in kd-Tree */
   
   void kdT_Delete (PhotonKdTree *kdt);
   /* Self-destruct */

#endif
