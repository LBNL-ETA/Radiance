/* 
   ==================================================================
   Photon map interface to out-of-core octree

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmapooc.h,v 1.2 2020/04/08 15:14:21 rschregle Exp $
*/



#ifndef PMAPOOC_H
   #define PMAPOOC_H
   
   #include "oocnn.h"

   
   
   /* Suffixes for octree filenames */
   /* #define PMAP_OOC_NODESUFFIX   ".node"
      #define PMAP_OOC_HEAPSUFFIX   ".heap" */
   #define PMAP_OOC_LEAFSUFFIX   ".leaf"
   
   /* Out-of-core octree constants */
   #define PMAP_OOC_NUMBLK       32    /* Num blocks for external sort */
   #define PMAP_OOC_BLKSIZE      1e8   /* Block size for external sort */
   #define PMAP_OOC_LEAFMAX      (OOC_OCTCNT_MAX)  /* Max photons per leaf */
   #define PMAP_OOC_MAXDEPTH     (OOC_MORTON_BITS) /* Max octree depth */



   typedef  OOC_SearchQueueNode  PhotonSearchQueueNode;
   typedef  OOC_SearchQueue      PhotonSearchQueue;
   typedef  OOC_Octree           PhotonStorage;   
   typedef  unsigned             PhotonIdx;

   /* Forward declarations to break dependency loop with pmapdata.h */
   struct PhotonMap;


      
   void OOC_BuildPhotonMap (struct PhotonMap *pmap, unsigned numProc);
   /* Build out-of-core octree pmap -> store from photons in unsorted
    * heapfile pmap -> heap and generate nodes and leaf file with prefix
    * pmap -> fileName.  Photon map construction may be parallelised if
    * numProc > 1, if supported.  The heap is destroyed on return.  */
   
   int OOC_SavePhotons (const struct PhotonMap *pmap, FILE *out);
   /* Save photons in out-of-core octree to file. Return -1 on error, 
    * else 0 */
   
   int OOC_LoadPhotons (struct PhotonMap *pmap, FILE *in);
   /* Load photons in out-of-core octree from file. Return -1 on error, 
    * else 0 */

   void OOC_InitFindPhotons (struct PhotonMap *pmap);
   /* Initialise NN search queue prior to calling kdT_FindPhotons() */
   
   int OOC_FindPhotons (struct PhotonMap* pmap, const FVECT pos, 
                        const FVECT norm);
   /* Locate pmap -> squeue.len nearest photons to pos with similar normal
    * (NULL for volume photons) and return in search queue pmap -> squeue,
    * starting with the further photon at pmap -> squeue.node. Return -1 
    * if none found, else 0. */

   int OOC_Find1Photon (struct PhotonMap* pmap, const FVECT pos, 
                        const FVECT norm, Photon *photon);
   /* Locate single nearest photon to pos with similar normal. Return -1 
    * if none found, else 0. */   
    
   int OOC_GetPhoton (struct PhotonMap *pmap, PhotonIdx idx, 
                      Photon *photon);
   /* Retrieve photon referenced by idx from leaf file and return -1 on
    * error, else 0. */                      

   Photon *OOC_GetNearestPhoton (const PhotonSearchQueue *squeue, 
                                 PhotonIdx idx);
   /* Retrieve photon from NN search queue after OOC_FindPhotons() */

   PhotonIdx OOC_FirstPhoton (const struct PhotonMap* pmap);
   /* Return index to first photon in octree */    
#endif
