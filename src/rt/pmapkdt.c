/* 
   ======================================================================
   In-core kd-tree for photon map

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: pmapkdt.c,v 1.6 2020/04/08 15:14:21 rschregle Exp $
*/



#include "pmapdata.h"   /* Includes pmapkdt.h */
#include "source.h"
#include "otspecial.h"
#include "random.h"




void kdT_Null (PhotonKdTree *kdt)
{
   kdt -> nodes = NULL;
}



static unsigned long kdT_MedianPartition (const Photon *heap, 
                                          unsigned long *heapIdx,
                                          unsigned long *heapXdi, 
                                          unsigned long left, 
                                          unsigned long right, unsigned dim)
/* Returns index to median in heap from indices left to right 
   (inclusive) in dimension dim. The heap is partitioned relative to 
   median using a quicksort algorithm. The heap indices in heapIdx are
   sorted rather than the heap itself. */
{
   const float    *p;
   unsigned long  l, r, lg2, n2, m, n = right - left + 1;
   unsigned       d;
   
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
      
      if (l >= m) 
         right = l - 1;
      if (l <= m) 
         left = l + 1;
   }
   
   /* Once left & right have converged at m, we have found the median */
   return m;
}



static void kdT_Build (Photon *heap, unsigned long *heapIdx,
                       unsigned long *heapXdi, const float min [3], 
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
   float                maxLeft [3], minRight [3];
   Photon               rootNode;
   unsigned             d;
   
   /* Choose median for dimension with largest spread and partition 
      accordingly */
   const float          d0 = max [0] - min [0], 
                        d1 = max [1] - min [1], 
                        d2 = max [2] - min [2];
   const unsigned char  dim = d0 > d1 ? d0 > d2 ? 0 : 2 
                                      : d1 > d2 ? 1 : 2;
   const unsigned long  median = left == right 
                                 ? left 
                                 : kdT_MedianPartition(heap, heapIdx, heapXdi,
                                                       left, right, dim);
   
   /* Place median at root of current subtree. This consists of swapping 
      the median and the root nodes and updating the heap indices */
   memcpy(&rootNode, heap + heapIdx [median], sizeof(Photon));
   memcpy(heap + heapIdx [median], heap + root - 1, sizeof(Photon));
   rootNode.discr = dim;
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
      kdT_Build(heap, heapIdx, heapXdi, min, maxLeft, left, median - 1, 
                root << 1);
                
   if (right > median) 
      kdT_Build(heap, heapIdx, heapXdi, minRight, max, median + 1, right, 
                (root << 1) + 1);
}



void kdT_BuildPhotonMap (struct PhotonMap *pmap)
{
   Photon         *nodes;
   unsigned long  i;
   unsigned long  *heapIdx,        /* Photon index array */
                  *heapXdi;        /* Reverse index to heapIdx */
   
   /* Allocate kd-tree nodes and load photons from heap file */
   if (!(nodes = calloc(pmap -> numPhotons, sizeof(Photon))))
      error(SYSTEM, "failed in-core heap allocation in kdT_BuildPhotonMap");
     
   rewind(pmap -> heap);
   i = fread(nodes, sizeof(Photon), pmap -> numPhotons, pmap -> heap);
   if (i !=
       pmap -> numPhotons)
      error(SYSTEM, "failed loading photon heap in kdT_BuildPhotonMap");
      
   pmap -> store.nodes = nodes;
   heapIdx = calloc(pmap -> numPhotons, sizeof(unsigned long));
   heapXdi = calloc(pmap -> numPhotons, sizeof(unsigned long));
   if (!heapIdx || !heapXdi)
      error(SYSTEM, "failed heap index allocation in kdT_BuildPhotonMap");
         
   /* Initialize index arrays */
   for (i = 0; i < pmap -> numPhotons; i++)
      heapXdi [i] = heapIdx [i] = i;
      
   /* Build kd-tree */
   kdT_Build(nodes, heapIdx, heapXdi, pmap -> minPos, pmap -> maxPos, 0, 
             pmap -> numPhotons - 1, 1);
                
   /* Cleanup */
   free(heapIdx);
   free(heapXdi);
}



int kdT_SavePhotons (const struct PhotonMap *pmap, FILE *out)
{
   unsigned long        i, j;
   Photon               *p = (Photon*)pmap -> store.nodes;
   
   for (i = 0; i < pmap -> numPhotons; i++, p++) {
      /* Write photon attributes */
      for (j = 0; j < 3; j++)
         putflt(p -> pos [j], out);
         
      /* Bytewise dump otherwise we have portability probs */
      for (j = 0; j < 3; j++) 
         putint(p -> norm [j], 1, out);
         
#ifdef PMAP_FLOAT_FLUX
      for (j = 0; j < 3; j++) 
         putflt(p -> flux [j], out);
#else
      for (j = 0; j < 4; j++) 
         putint(p -> flux [j], 1, out);
#endif

      putint(p -> primary, sizeof(p -> primary), out);
      putint(p -> flags, 1, out);
            
      if (ferror(out))
         return -1;
   }
   
   return 0;
}



int kdT_LoadPhotons (struct PhotonMap *pmap, FILE *in)
{
   unsigned long  i, j;
   Photon         *p;
   
   /* Allocate kd-tree based on initialised pmap -> numPhotons */
   pmap -> store.nodes = calloc(sizeof(Photon), pmap -> numPhotons);
   if (!pmap -> store.nodes) 
      error(SYSTEM, "failed kd-tree allocation in kdT_LoadPhotons");
      
   /* Get photon attributes */   
   for (i = 0, p = pmap -> store.nodes; i < pmap -> numPhotons; i++, p++) {
      for (j = 0; j < 3; j++) 
         p -> pos [j] = getflt(in);
         
      /* Bytewise grab otherwise we have portability probs */
      for (j = 0; j < 3; j++) 
         p -> norm [j] = getint(1, in);
         
#ifdef PMAP_FLOAT_FLUX
      for (j = 0; j < 3; j++) 
         p -> flux [j] = getflt(in);
#else      
      for (j = 0; j < 4; j++) 
         p -> flux [j] = getint(1, in);
#endif

      p -> primary = getint(sizeof(p -> primary), in);
      p -> flags = getint(1, in);
      
      if (feof(in))
         return -1;
   }
   
   return 0;
}



void kdT_InitFindPhotons (struct PhotonMap *pmap)
{
   pmap -> squeue.len = pmap -> maxGather + 1;
   pmap -> squeue.node = calloc(pmap -> squeue.len,
                                sizeof(PhotonSearchQueueNode));
   if (!pmap -> squeue.node) 
      error(SYSTEM, "can't allocate photon search queue");
}



static void kdT_FindNearest (PhotonMap *pmap, const float pos [3], 
                             const float norm [3], unsigned long node)
/* Recursive part of kdT_FindPhotons(). Locate pmap -> squeue.len nearest
 * neighbours to pos with similar normal and return in search queue starting
 * at pmap -> squeue.node.  Note that all heap and queue indices are
 * 1-based, but accesses to the arrays are 0-based!  */
{
   Photon                  *p = (Photon*)pmap -> store.nodes + node - 1;
   unsigned                i, j;
   /* Signed distance to current photon's splitting plane */
   float                   d = pos [p -> discr] - p -> pos [p -> discr], 
                           d2 = d * d, dv [3];
   PhotonSearchQueueNode*  sq = pmap -> squeue.node;
   const unsigned          sqSize = pmap -> squeue.len;
   
   /* Search subtree closer to pos first; exclude other subtree if the 
      distance to the splitting plane is greater than maxDist */
   if (d < 0) {
      if (node << 1 <= pmap -> numPhotons) 
         kdT_FindNearest(pmap, pos, norm, node << 1);
         
      if (d2 < pmap -> maxDist2 && node << 1 < pmap -> numPhotons) 
         kdT_FindNearest(pmap, pos, norm, (node << 1) + 1);
   }
   else {
      if (node << 1 < pmap -> numPhotons) 
         kdT_FindNearest(pmap, pos, norm, (node << 1) + 1);
         
      if (d2 < pmap -> maxDist2 && node << 1 <= pmap -> numPhotons) 
         kdT_FindNearest(pmap, pos, norm, node << 1);
   }

   /* Reject photon if normal faces away (ignored for volume photons) with
    * tolerance to account for perturbation; note photon normal is coded
    * in range [-127,127], hence we factor this in */
   if (norm && DOT(norm, p -> norm) <= PMAP_NORM_TOL * 127 * frandom())
      return;
      
   if (isContribPmap(pmap)) {
      /* Lookup in contribution photon map; filter according to emitting
       * light source if contrib list set, else accept all */
       
      if (pmap -> srcContrib) {
         OBJREC *srcMod; 
         const int srcIdx = photonSrcIdx(pmap, p);
         
         if (srcIdx < 0 || srcIdx >= nsources)
            error(INTERNAL, "invalid light source index in photon map");
         
         srcMod = findmaterial(source [srcIdx].so);

         /* Reject photon if contributions from light source which emitted it
          * are not sought */
         if (!lu_find(pmap -> srcContrib, srcMod -> oname) -> data)
            return;
      }

      /* Reject non-caustic photon if lookup for caustic contribs */
      if (pmap -> lookupCaustic & !p -> caustic)
         return;
   }
   
   /* Squared distance to current photon (note dist2() requires doubles) */
   VSUB(dv, pos, p -> pos);
   d2 = DOT(dv, dv);

   /* Accept photon if closer than current max dist & add to priority queue */
   if (d2 < pmap -> maxDist2) {
      if (pmap -> squeue.tail < sqSize) {
         /* Priority queue not full; append photon and restore heap */
         i = ++pmap -> squeue.tail;
         
         while (i > 1 && sq [(i >> 1) - 1].dist2 <= d2) {
            sq [i - 1].idx    = sq [(i >> 1) - 1].idx;
            sq [i - 1].dist2  = sq [(i >> 1) - 1].dist2;
            i >>= 1;
         }
         
         sq [--i].idx = (PhotonIdx)p;
         sq [i].dist2 = d2;
         /* Update maxDist if we've just filled the queue */
         if (pmap -> squeue.tail >= pmap -> squeue.len)
            pmap -> maxDist2 = sq [0].dist2;
      }
      else {
         /* Priority queue full; replace maximum, restore heap, and 
            update maxDist */
         i = 1;
         
         while (i <= sqSize >> 1) {
            j = i << 1;
            if (j < sqSize && sq [j - 1].dist2 < sq [j].dist2) 
               j++;
            if (d2 >= sq [j - 1].dist2) 
               break;
            sq [i - 1].idx    = sq [j - 1].idx;
            sq [i - 1].dist2  = sq [j - 1].dist2;
            i = j;
         }
         
         sq [--i].idx = (PhotonIdx)p;
         sq [i].dist2 = d2;
         pmap -> maxDist2 = sq [0].dist2;
      }
   }
}



int kdT_FindPhotons (struct PhotonMap *pmap, const FVECT pos, 
                     const FVECT norm)
{
   float p [3], n [3];
   
   /* Photon pos & normal stored at lower precision */
   VCOPY(p, pos);
   if (norm)
      VCOPY(n, norm);
   kdT_FindNearest(pmap, p, norm ? n : NULL, 1);
   
   /* Return success or failure (empty queue => none found) */
   return pmap -> squeue.tail ? 0 : -1;
}



static void kdT_Find1Nearest (PhotonMap *pmap, const float pos [3], 
                              const float norm [3], Photon **photon, 
                              unsigned long node)
/* Recursive part of kdT_Find1Photon().  Locate single nearest neighbour to
 * pos with similar normal.  Note that all heap and queue indices are
 * 1-based, but accesses to the arrays are 0-based!  */
{
   Photon   *p = (Photon*)pmap -> store.nodes + node - 1;
   /* Signed distance to current photon's splitting plane */
   float    d  = pos [p -> discr] - p -> pos [p -> discr], d2 = d * d, 
            dv [3];
   
   /* Search subtree closer to pos first; exclude other subtree if the 
      distance to the splitting plane is greater than maxDist */
   if (d < 0) {
      if (node << 1 <= pmap -> numPhotons) 
         kdT_Find1Nearest(pmap, pos, norm, photon, node << 1);
         
      if (d2 < pmap -> maxDist2 && node << 1 < pmap -> numPhotons) 
         kdT_Find1Nearest(pmap, pos, norm, photon, (node << 1) + 1);
   }
   else {
      if (node << 1 < pmap -> numPhotons) 
         kdT_Find1Nearest(pmap, pos, norm, photon, (node << 1) + 1);
         
      if (d2 < pmap -> maxDist2 && node << 1 <= pmap -> numPhotons) 
         kdT_Find1Nearest(pmap, pos, norm, photon, node << 1);
   }
   
   /* Squared distance to current photon */
   VSUB(dv, pos, p -> pos);
   d2 = DOT(dv, dv);
   
   if (d2 < pmap -> maxDist2 && 
       (!norm || DOT(norm, p -> norm) > PMAP_NORM_TOL * 127 * frandom())) {
      /* Closest photon so far with similar normal. We allow for tolerance
       * to account for perturbation in the latter; note the photon normal
       * is coded in the range [-127,127], hence we factor this in  */
      pmap -> maxDist2 = d2;
      *photon = p;
   }
}



int kdT_Find1Photon (struct PhotonMap *pmap, const FVECT pos, 
                     const FVECT norm, Photon *photon)
{
   float    p [3], n [3];
   Photon   *pnn = NULL;
   
   /* Photon pos & normal stored at lower precision */
   VCOPY(p, pos);
   if (norm)
      VCOPY(n, norm);   
   kdT_Find1Nearest(pmap, p, norm ? n : NULL, &pnn, 1);   
   if (!pnn)
      /* No photon found => failed */
      return -1;
   else {
      /* Copy found photon => successs */
      memcpy(photon, pnn, sizeof(Photon));
      return 0;
   }
}



int kdT_GetPhoton (const struct PhotonMap *pmap, PhotonIdx idx,
                   Photon *photon)
{
   memcpy(photon, idx, sizeof(Photon));   
   return 0;
}



Photon *kdT_GetNearestPhoton (const PhotonSearchQueue *squeue, PhotonIdx idx)
{
   return idx;
}



PhotonIdx kdT_FirstPhoton (const struct PhotonMap* pmap)
{
   return pmap -> store.nodes;
}



void kdT_Delete (PhotonKdTree *kdt)
{
   free(kdt -> nodes);
   kdt -> nodes = NULL;
}
