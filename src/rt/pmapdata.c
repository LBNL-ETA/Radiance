/* 
   ==================================================================
   Photon map data structures and kd-tree handling

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================   
   
   $Id: pmapdata.c,v 2.5 2015/05/20 14:44:12 greg Exp $
*/



#include "pmap.h"
#include "pmaprand.h"
#include "pmapmat.h"
#include "otypes.h"
#include "source.h"
#include "rcontrib.h"



PhotonMap *photonMaps [NUM_PMAP_TYPES] = {
   NULL, NULL, NULL, NULL, NULL, NULL
};



void initPhotonMap (PhotonMap *pmap, PhotonMapType t)
/* Init photon map 'n' stuff... */
{
   if (!pmap) 
      return;
      
   pmap -> heapSize = pmap -> heapEnd = 0;
   pmap -> heap = NULL;
   pmap -> squeue = NULL;
   pmap -> biasCompHist = NULL;
   pmap -> maxPos [0] = pmap -> maxPos [1] = pmap -> maxPos [2] = -FHUGE;
   pmap -> minPos [0] = pmap -> minPos [1] = pmap -> minPos [2] = FHUGE;
   pmap -> minGathered = pmap -> maxGathered = pmap -> totalGathered = 0;
   pmap -> gatherTolerance = gatherTolerance;
   pmap -> minError = pmap -> maxError = pmap -> rmsError = 0;
   pmap -> numDensity = 0;
   pmap -> distribRatio = 1;
   pmap -> type = t;

   /* Init local RNG state */
   pmap -> randState [0] = 10243;
   pmap -> randState [1] = 39829;
   pmap -> randState [2] = 9433;
   /* pmapSeed(25999, pmap -> randState); */
   pmapSeed(randSeed, pmap -> randState);
   
   /* Set up type-specific photon lookup callback */
   pmap -> lookup = pmapLookup [t];

   pmap -> primary = NULL;
   pmap -> primarySize = pmap -> primaryEnd = 0;
}



const PhotonPrimary* addPhotonPrimary (PhotonMap *pmap, const RAY *ray)
{
   PhotonPrimary *prim = NULL;
   FVECT dvec;
   
   if (!pmap || !ray)
      return NULL;
      
   /* Check if last primary ray has spawned photons (srcIdx >= 0, see
    * addPhoton()), in which case we keep it and allocate a new one;
    * otherwise we overwrite the unused entry */
   if (pmap -> primary && pmap -> primary [pmap -> primaryEnd].srcIdx >= 0)
      pmap -> primaryEnd++;
      
   if (!pmap -> primarySize || pmap -> primaryEnd >= pmap -> primarySize) {
      /* Allocate/enlarge array */
      pmap -> primarySize += pmap -> heapSizeInc;
      
      /* Counter wraparound? */
      if (pmap -> primarySize < pmap -> heapSizeInc)
         error(INTERNAL, "photon primary overflow");
         
      pmap -> primary = (PhotonPrimary *)realloc(pmap -> primary,
                                                 sizeof(PhotonPrimary) * 
                                                 pmap -> primarySize);
      if (!pmap -> primary)
         error(USER, "can't allocate photon primaries");
   }
   
   prim = pmap -> primary + pmap -> primaryEnd;
   
   /* Mark unused with negative source index until path spawns a photon (see
    * addPhoton()) */
   prim -> srcIdx = -1;
      
   /* Reverse incident direction to point to light source */
   dvec [0] = -ray -> rdir [0];
   dvec [1] = -ray -> rdir [1];
   dvec [2] = -ray -> rdir [2];
   prim -> dir = encodedir(dvec);

   VCOPY(prim -> pos, ray -> rop);
   
   return prim;
}



const Photon* addPhoton (PhotonMap* pmap, const RAY* ray)
{
   unsigned i;
   Photon* photon = NULL;
   COLOR photonFlux;
   
   /* Account for distribution ratio */
   if (!pmap || pmapRandom(pmap -> randState) > pmap -> distribRatio) 
      return NULL;
      
   /* Don't store on sources */
   if (ray -> robj > -1 && islight(objptr(ray -> ro -> omod) -> otype)) 
      return NULL;

#if 0
   if (inContribPmap(pmap)) {
      /* Adding contribution photon */
      if (ray -> parent && ray -> parent -> rtype & PRIMARY)
         /* Add primary photon ray if parent is primary; by putting this
          * here and checking the ray's immediate parent, we only add
          * primaries that actually contribute photons, and we only add them
          * once.  */
         addPhotonPrimary(pmap, ray -> parent);
      
      /* Save index to primary ray (remains unchanged if primary already in
       * array) */
      primary = pmap -> primaryEnd;
   }
#endif
   
#ifdef PMAP_ROI
   /* Store photon if within region of interest -- for egg-spurtz only! */
   if (ray -> rop [0] >= pmapROI [0] && ray -> rop [0] <= pmapROI [1] &&
       ray -> rop [1] >= pmapROI [2] && ray -> rop [1] <= pmapROI [3] &&
       ray -> rop [2] >= pmapROI [4] && ray -> rop [2] <= pmapROI [5])
#endif
   {       
      if (pmap -> heapEnd >= pmap -> heapSize) {
         /* Enlarge heap */
         pmap -> heapSize += pmap -> heapSizeInc;
         
         /* Counter wraparound? */
         if (pmap -> heapSize < pmap -> heapSizeInc)
            error(INTERNAL, "photon heap overflow");
         
         pmap -> heap = (Photon *)realloc(pmap -> heap, 
                                          sizeof(Photon) * pmap -> heapSize);
         if (!pmap -> heap) 
            error(USER, "can't allocate photon heap");
      } 

      photon = pmap -> heap + pmap -> heapEnd++;

      /* Adjust flux according to distribution ratio and ray weight */
      copycolor(photonFlux, ray -> rcol);   
      scalecolor(photonFlux, 
                 ray -> rweight / (pmap -> distribRatio ? pmap -> distribRatio
                                                        : 1));
      setPhotonFlux(photon, photonFlux);
               
      /* Set photon position and flags */
      VCOPY(photon -> pos, ray -> rop);
      photon -> flags = PMAP_CAUSTICRAY(ray) ? PMAP_CAUST_BIT : 0;

      /* Set primary ray index and mark as used for contrib photons */
      if (isContribPmap(pmap)) {
         photon -> primary = pmap -> primaryEnd;
         pmap -> primary [pmap -> primaryEnd].srcIdx = ray -> rsrc;
      }
      else photon -> primary = 0;
      
      /* Update min and max positions & set normal */
      for (i = 0; i <= 2; i++) {
         if (photon -> pos [i] < pmap -> minPos [i]) 
            pmap -> minPos [i] = photon -> pos [i];
         if (photon -> pos [i] > pmap -> maxPos [i]) 
            pmap -> maxPos [i] = photon -> pos [i];
         photon -> norm [i] = 127.0 * (isVolumePmap(pmap) ? ray -> rdir [i] 
                                                          : ray -> ron [i]);
      }
   }
            
   return photon;
}



static void nearestNeighbours (PhotonMap* pmap, const float pos [3], 
                               const float norm [3], unsigned long node)
/* Recursive part of findPhotons(..).
   Note that all heap and priority queue index handling is 1-based, but 
   accesses to the arrays are 0-based! */
{
   Photon* p = &pmap -> heap [node - 1];
   unsigned i, j;
   /* Signed distance to current photon's splitting plane */
   float d = pos [photonDiscr(*p)] - p -> pos [photonDiscr(*p)], 
         d2 = d * d;
   PhotonSQNode* sq = pmap -> squeue;
   const unsigned sqSize = pmap -> squeueSize;
   float dv [3];
   
   /* Search subtree closer to pos first; exclude other subtree if the 
      distance to the splitting plane is greater than maxDist */
   if (d < 0) {
      if (node << 1 <= pmap -> heapSize) 
         nearestNeighbours(pmap, pos, norm, node << 1);
      if (d2 < pmap -> maxDist && node << 1 < pmap -> heapSize) 
         nearestNeighbours(pmap, pos, norm, (node << 1) + 1);
   }
   else {
      if (node << 1 < pmap -> heapSize) 
         nearestNeighbours(pmap, pos, norm, (node << 1) + 1);
      if (d2 < pmap -> maxDist && node << 1 <= pmap -> heapSize) 
         nearestNeighbours(pmap, pos, norm, node << 1);
   }

   /* Reject photon if normal faces away (ignored for volume photons) */
   if (norm && DOT(norm, p -> norm) <= 0)
      return;
      
   if (isContribPmap(pmap) && pmap -> srcContrib) {
      /* Lookup in contribution photon map */
      OBJREC *srcMod; 
      const int srcIdx = photonSrcIdx(pmap, p);
      
      if (srcIdx < 0 || srcIdx >= nsources)
         error(INTERNAL, "invalid light source index in photon map");
      
      srcMod = objptr(source [srcIdx].so -> omod);

      /* Reject photon if contributions from light source which emitted it
       * are not sought */
      if (!lu_find(pmap -> srcContrib, srcMod -> oname) -> data)
         return;

      /* Reject non-caustic photon if lookup for caustic contribs */
      if (pmap -> lookupFlags & PMAP_CAUST_BIT & ~p -> flags)
         return;
   }
   
   /* Squared distance to current photon */
   dv [0] = pos [0] - p -> pos [0];
   dv [1] = pos [1] - p -> pos [1];
   dv [2] = pos [2] - p -> pos [2];
   d2 = DOT(dv, dv);

   /* Accept photon if closer than current max dist & add to priority queue */
   if (d2 < pmap -> maxDist) {
      if (pmap -> squeueEnd < sqSize) {
         /* Priority queue not full; append photon and restore heap */
         i = ++pmap -> squeueEnd;
         
         while (i > 1 && sq [(i >> 1) - 1].dist <= d2) {
            sq [i - 1].photon = sq [(i >> 1) - 1].photon;
            sq [i - 1].dist = sq [(i >> 1) - 1].dist;
            i >>= 1;
         }
         
         sq [--i].photon = p;
         sq [i].dist = d2;
         /* Update maxDist if we've just filled the queue */
         if (pmap -> squeueEnd >= pmap -> squeueSize)
            pmap -> maxDist = sq [0].dist;
      }
      else {
         /* Priority queue full; replace maximum, restore heap, and 
            update maxDist */
         i = 1;
         
         while (i <= sqSize >> 1) {
            j = i << 1;
            if (j < sqSize && sq [j - 1].dist < sq [j].dist) 
               j++;
            if (d2 >= sq [j - 1].dist) 
               break;
            sq [i - 1].photon = sq [j - 1].photon;
            sq [i - 1].dist = sq [j - 1].dist;
            i = j;
         }
         
         sq [--i].photon = p;
         sq [i].dist = d2;
         pmap -> maxDist = sq [0].dist;
      }
   }
}



/* Dynamic max photon search radius increase and reduction factors */
#define PMAP_MAXDIST_INC 4
#define PMAP_MAXDIST_DEC 0.9

/* Num successful lookups before reducing in max search radius */
#define PMAP_MAXDIST_CNT 1000

/* Threshold below which we assume increasing max radius won't help */
#define PMAP_SHORT_LOOKUP_THRESH 1

void findPhotons (PhotonMap* pmap, const RAY* ray)
{
   float pos [3], norm [3];
   int redo = 0;
   
   if (!pmap -> squeue) {
      /* Lazy init priority queue */
      pmap -> squeueSize = pmap -> maxGather + 1;
      pmap -> squeue = (PhotonSQNode*)malloc(pmap -> squeueSize * 
                                             sizeof(PhotonSQNode));
      if (!pmap -> squeue) 
         error(USER, "can't allocate photon priority queue");
         
      pmap -> minGathered = pmap -> maxGather;
      pmap -> maxGathered = pmap -> minGather;
      pmap -> totalGathered = 0;
      pmap -> numLookups = pmap -> numShortLookups = 0;
      pmap -> shortLookupPct = 0;
      pmap -> minError = FHUGE;
      pmap -> maxError = -FHUGE;
      pmap -> rmsError = 0;
#ifdef PMAP_MAXDIST_ABS
      /* Treat maxDistCoeff as an *absolute* and *fixed* max search radius.
         Primarily intended for debugging; FOR ZE ECKSPURTZ ONLY! */
      pmap -> maxDist0 = pmap -> maxDistLimit = maxDistCoeff;
#else
      /* Maximum search radius limit based on avg photon distance to
       * centre of gravity */
      pmap -> maxDist0 = pmap -> maxDistLimit = 
         maxDistCoeff * pmap -> squeueSize * pmap -> CoGdist / 
         pmap -> heapSize;
#endif      
   }

   do {
      pmap -> squeueEnd = 0;
      pmap -> maxDist = pmap -> maxDist0;
         
      /* Search position is ray -> rorg for volume photons, since we have no 
         intersection point. Normals are ignored -- these are incident 
         directions). */   
      if (isVolumePmap(pmap)) {
         VCOPY(pos, ray -> rorg);
         nearestNeighbours(pmap, pos, NULL, 1);
      }
      else {
         VCOPY(pos, ray -> rop);
         VCOPY(norm, ray -> ron);
         nearestNeighbours(pmap, pos, norm, 1);
      }

#ifndef PMAP_MAXDIST_ABS      
      if (pmap -> squeueEnd < pmap -> squeueSize * pmap -> gatherTolerance) {
         /* Short lookup; too few photons found */
         if (pmap -> squeueEnd > PMAP_SHORT_LOOKUP_THRESH) {
            /* Ignore short lookups which return fewer than
             * PMAP_SHORT_LOOKUP_THRESH photons under the assumption there
             * really are no photons in the vicinity, and increasing the max
             * search radius therefore won't help */
   #ifdef PMAP_LOOKUP_WARN
            sprintf(errmsg, 
                    "%d/%d %s photons found at (%.2f,%.2f,%.2f) on %s",
                    pmap -> squeueEnd, pmap -> squeueSize, 
                    pmapName [pmap -> type], pos [0], pos [1], pos [2], 
                    ray -> ro ? ray -> ro -> oname : "<null>");
            error(WARNING, errmsg);         
   #endif             

            if (pmap -> maxDist0 < pmap -> maxDistLimit) {
               /* Increase max search radius if below limit & redo search */
               pmap -> maxDist0 *= PMAP_MAXDIST_INC;
   #ifdef PMAP_LOOKUP_REDO
               redo = 1;
   #endif               

   #ifdef PMAP_LOOKUP_WARN
               sprintf(errmsg, 
                       redo ? "restarting photon lookup with max radius %.1e"
                            : "max photon lookup radius adjusted to %.1e",
                       pmap -> maxDist0);
               error(WARNING, errmsg);
   #endif
            }
   #ifdef PMAP_LOOKUP_REDO
            else {
               sprintf(errmsg, "max photon lookup radius clamped to %.1e",
                       pmap -> maxDist0);
               error(WARNING, errmsg);
            }
   #endif
         }
         
         /* Reset successful lookup counter */
         pmap -> numLookups = 0;
      }
      else {
         /* Increment successful lookup counter and reduce max search radius if
          * wraparound */
         pmap -> numLookups = (pmap -> numLookups + 1) % PMAP_MAXDIST_CNT;
         if (!pmap -> numLookups)
            pmap -> maxDist0 *= PMAP_MAXDIST_DEC;
            
         redo = 0;
      }
#endif
   } while (redo);
}



static void nearest1Neighbour (PhotonMap *pmap, const float pos [3], 
                               const float norm [3], Photon **photon, 
                               unsigned long node)
/* Recursive part of find1Photon(..).
   Note that all heap index handling is 1-based, but accesses to the 
   arrays are 0-based! */
{
   Photon *p = pmap -> heap + node - 1;
   /* Signed distance to current photon's splitting plane */
   float d = pos [photonDiscr(*p)] - p -> pos [photonDiscr(*p)], 
         d2 = d * d;
   float dv [3];
   
   /* Search subtree closer to pos first; exclude other subtree if the 
      distance to the splitting plane is greater than maxDist */
   if (d < 0) {
      if (node << 1 <= pmap -> heapSize) 
         nearest1Neighbour(pmap, pos, norm, photon, node << 1);
      if (d2 < pmap -> maxDist && node << 1 < pmap -> heapSize) 
         nearest1Neighbour(pmap, pos, norm, photon, (node << 1) + 1);
   }
   else {
      if (node << 1 < pmap -> heapSize) 
         nearest1Neighbour(pmap, pos, norm, photon, (node << 1) + 1);
      if (d2 < pmap -> maxDist && node << 1 <= pmap -> heapSize) 
         nearest1Neighbour(pmap, pos, norm, photon, node << 1);
   }
   
   /* Squared distance to current photon */
   dv [0] = pos [0] - p -> pos [0];
   dv [1] = pos [1] - p -> pos [1];
   dv [2] = pos [2] - p -> pos [2];
   d2 = DOT(dv, dv);
   
   if (d2 < pmap -> maxDist && DOT(norm, p -> norm) > 0) {
      /* Closest photon so far with similar normal */
      pmap -> maxDist = d2;
      *photon = p;
   }
}



Photon* find1Photon (PhotonMap *pmap, const RAY* ray)
{
   float fpos [3], norm [3];
   Photon* photon = NULL;

   VCOPY(fpos, ray -> rop);
   VCOPY(norm, ray -> ron);
   pmap -> maxDist = thescene.cusize;
   nearest1Neighbour(pmap, fpos, norm, &photon, 1);
   
   return photon;
}



static unsigned long medianPartition (const Photon* heap, 
                                      unsigned long* heapIdx,
                                      unsigned long* heapXdi, 
                                      unsigned long left, 
                                      unsigned long right, unsigned dim)
/* Returns index to median in heap from indices left to right 
   (inclusive) in dimension dim. The heap is partitioned relative to 
   median using a quicksort algorithm. The heap indices in heapIdx are
   sorted rather than the heap itself. */
{
   register const float* p;
   const unsigned long n = right - left + 1;
   register unsigned long l, r, lg2, n2, m;
   register unsigned d;
   
   /* Round down n to nearest power of 2 */
   for (lg2 = 0, n2 = n; n2 > 1; n2 >>= 1, ++lg2);
   n2 = 1 << lg2;
   
   /* Determine median position; this takes into account the fact that 
      only the last level in the heap can be partially empty, and that 
      it fills from left to right */
   m = left + ((n - n2) > (n2 >> 1) - 1 ? n2 - 1 : n - (n2 >> 1));
   
   while (right > left) {
      /* Pivot node */
      p = heap [heapIdx [right]].pos;
      l = left;
      r = right - 1;
      
      /* l & r converge, swapping elements out of order with respect to 
         pivot node. Identical keys are resolved by cycling through
         dim. The convergence point is then the pivot's position. */
      do {
         while (l <= r) {
            d = dim;
            
            while (heap [heapIdx [l]].pos [d] == p [d]) {
               d = (d + 1) % 3;
               
               if (d == dim) {
                  /* Ignore dupes? */
                  error(WARNING, "duplicate keys in photon heap");
                  l++;
                  break;
               }
            }
            
            if (heap [heapIdx [l]].pos [d] < p [d]) 
               l++;
            else break;
         }
         
         while (r > l) {
            d = dim;
            
            while (heap [heapIdx [r]].pos [d] == p [d]) {
               d = (d + 1) % 3;
               
               if (d == dim) {
                  /* Ignore dupes? */
                  error(WARNING, "duplicate keys in photon heap");
                  r--;
                  break;
               }
            }
            
            if (heap [heapIdx [r]].pos [d] > p [d]) 
               r--;
            else break;
         }
         
         /* Swap indices (not the nodes they point to) */
         n2 = heapIdx [l];
         heapIdx [l] = heapIdx [r];
         heapIdx [r] = n2;
         /* Update reverse indices */
         heapXdi [heapIdx [l]] = l;
         heapXdi [n2] = r;
      } while (l < r);
      
      /* Swap indices of convergence and pivot nodes */
      heapIdx [r] = heapIdx [l];
      heapIdx [l] = heapIdx [right];
      heapIdx [right] = n2;
      /* Update reverse indices */
      heapXdi [heapIdx [r]] = r;
      heapXdi [heapIdx [l]] = l;
      heapXdi [n2] = right;
      if (l >= m) right = l - 1;
      if (l <= m) left = l + 1;
   }
   
   /* Once left & right have converged at m, we have found the median */
   return m;
}



void buildHeap (Photon* heap, unsigned long* heapIdx,
                unsigned long* heapXdi, const float min [3], 
                const float max [3], unsigned long left, 
                unsigned long right, unsigned long root)
/* Recursive part of balancePhotons(..). Builds heap from subarray
   defined by indices left and right. min and max are the minimum resp. 
   maximum photon positions in the array. root is the index of the
   current subtree's root, which corresponds to the median's 1-based 
   index in the heap. heapIdx are the balanced heap indices. The heap
   is accessed indirectly through these. heapXdi are the reverse indices 
   from the heap to heapIdx so that heapXdi [heapIdx [i]] = i. */
{
   float maxLeft [3], minRight [3];
   Photon rootNode;
   unsigned d;
   
   /* Choose median for dimension with largest spread and partition 
      accordingly */
   const float d0 = max [0] - min [0], 
               d1 = max [1] - min [1], 
               d2 = max [2] - min [2];
   const unsigned char dim = d0 > d1 ? d0 > d2 ? 0 : 2 
                                     : d1 > d2 ? 1 : 2;
   const unsigned long median = 
      left == right ? left 
                    : medianPartition(heap, heapIdx, heapXdi, left, right, dim);
   
   /* Place median at root of current subtree. This consists of swapping 
      the median and the root nodes and updating the heap indices */
   memcpy(&rootNode, heap + heapIdx [median], sizeof(Photon));
   memcpy(heap + heapIdx [median], heap + root - 1, sizeof(Photon));
   setPhotonDiscr(rootNode, dim);
   memcpy(heap + root - 1, &rootNode, sizeof(Photon));
   heapIdx [heapXdi [root - 1]] = heapIdx [median];
   heapXdi [heapIdx [median]] = heapXdi [root - 1];
   heapIdx [median] = root - 1;
   heapXdi [root - 1] = median;
   
   /* Update bounds for left and right subtrees and recurse on them */
   for (d = 0; d <= 2; d++)
      if (d == dim) 
         maxLeft [d] = minRight [d] = rootNode.pos [d];
      else {
         maxLeft [d] = max [d];
         minRight [d] = min [d];
      }
      
   if (left < median) 
      buildHeap(heap, heapIdx, heapXdi, min, maxLeft, 
                left, median - 1, root << 1);
                
   if (right > median) 
      buildHeap(heap, heapIdx, heapXdi, minRight, max, 
                median + 1, right, (root << 1) + 1);
}



void balancePhotons (PhotonMap* pmap, double *photonFlux)
{
   Photon *heap = pmap -> heap;
   unsigned long i;
   unsigned long *heapIdx;        /* Photon index array */
   unsigned long *heapXdi;        /* Reverse index to heapIdx */
   unsigned j;
   COLOR flux;
   /* Need doubles here to reduce errors from increment */
   double avgFlux [3] = {0, 0, 0}, CoG [3] = {0, 0, 0}, CoGdist = 0;
   FVECT d;
   
   if (pmap -> heapEnd) {
      pmap -> heapSize = pmap -> heapEnd;
      heapIdx = (unsigned long*)malloc(pmap -> heapSize * 
                                       sizeof(unsigned long));
      heapXdi = (unsigned long*)malloc(pmap -> heapSize * 
                                       sizeof(unsigned long));
      if (!heapIdx || !heapXdi)
         error(USER, "can't allocate heap index");
         
      for (i = 0; i < pmap -> heapSize; i++) {
         /* Initialize index arrays */
         heapXdi [i] = heapIdx [i] = i;
         getPhotonFlux(heap + i, flux);

         /* Scale photon's flux (hitherto normalised to 1 over RGB); in case
          * of a contrib photon map, this is done per light source, and
          * photonFlux is assumed to be an array */
         if (photonFlux) {
            scalecolor(flux, photonFlux [isContribPmap(pmap) ? 
                                         photonSrcIdx(pmap, heap + i) : 0]);
            setPhotonFlux(heap + i, flux);
         }

         /* Need a double here */
         addcolor(avgFlux, flux);
         
         /* Add photon position to centre of gravity */
         for (j = 0; j < 3; j++)
            CoG [j] += heap [i].pos [j];
      }
      
      /* Average photon positions to get centre of gravity */
      for (j = 0; j < 3; j++)
         pmap -> CoG [j] = CoG [j] /= pmap -> heapSize;
      
      /* Compute average photon distance to CoG */
      for (i = 0; i < pmap -> heapSize; i++) {
         VSUB(d, heap [i].pos, CoG);
         CoGdist += DOT(d, d);
      }
      
      pmap -> CoGdist = CoGdist /= pmap -> heapSize;
      
      /* Average photon flux based on RGBE representation */
      scalecolor(avgFlux, 1.0 / pmap -> heapSize);
      copycolor(pmap -> photonFlux, avgFlux);
      
      /* Build kd-tree */
      buildHeap(pmap -> heap, heapIdx, heapXdi, pmap -> minPos, 
                pmap -> maxPos, 0, pmap -> heapSize - 1, 1);
                
      free(heapIdx);
      free(heapXdi);
   }
}



void deletePhotons (PhotonMap* pmap)
{
   free(pmap -> heap);
   free(pmap -> squeue);
   free(pmap -> biasCompHist);
   
   pmap -> heapSize = 0;
   pmap -> minGather = pmap -> maxGather =
      pmap -> squeueSize = pmap -> squeueEnd = 0;
}
