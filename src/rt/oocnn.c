#ifndef lint
static const char RCSid[] = "$Id: oocnn.c,v 2.3 2021/03/23 00:16:38 rschregle Exp $";
#endif


/* 
   =========================================================================
   k-nearest neighbour lookup routines for out-of-core octree data structure
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   =========================================================================
   
   $Id: oocnn.c,v 2.3 2021/03/23 00:16:38 rschregle Exp $
*/


#if !defined(_WIN32) && !defined(_WIN64) || defined(PMAP_OOC)
/* No Windoze support for now */

#include "oocnn.h"
#include "oocsort.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>



#ifdef DEBUG_OOC_NN
static int OOC_SearchQueueCheck (OOC_SearchQueue *queue, const FVECT key,
                                 RREAL *(*keyFunc)(const void*), 
                                 unsigned root)
/* Priority queue sanity check */
{
   unsigned                   kid;
   const OOC_SearchQueueNode  *qn = queue -> node;
   void                       *rec;
   float                      d2;
   
   if (root < queue -> tail) {
      rec = OOC_GetNearest(queue, qn [root].idx);
      d2 = dist2(keyFunc(rec), key);
      
      if (qn [root].dist2 != d2)
         return -1;
      
      if ((kid = (root << 1) + 1) < queue -> tail) {
         if (qn [kid].dist2 > qn [root].dist2)
            return -1;
         else return OOC_SearchQueueCheck(queue, key, keyFunc, kid);
      }
         
      if (++kid < queue -> tail) {
         if (qn [kid].dist2 > qn [root].dist2)
            return -1;
         else return OOC_SearchQueueCheck(queue, key, keyFunc, kid);
      }
   }
   
   return 0;
}
#endif



static float OOC_PutNearest (OOC_SearchQueue *queue, float d2, void *rec)
/* Insert data record with SQUARED distance to query point into search
 * priority queue, maintaining the most distant record at the queue head. 
 * If the queue is full, the new record is only inserted if it is closer
 * than the queue head.
 * Returns the new maximum SQUARED distance at the head if the queue is
 * full.  Otherwise returns -1, indicating a maximum for the entire queue is
 * as yet undefined
 * The data record is copied into the queue's local record buffa for
 * post-search retrieval to minimise redundant disk access. Note that it
 * suffices to only rearrange the corresponding indices in the queue nodes
 * when restoring the priority queue after every insertion, rather than
 * moving the actual records. */
{
   OOC_SearchQueueNode  *qn = queue -> node;
   unsigned             root, kid, kid2, rootIdx;
   float                d2max = -1;    /* Undefined max distance ^2 */
   
   /* The queue is represented as a linearised binary tree with the root
    * corresponding to the queue head, and the tail corresponding to the
    * last leaf */
   if (queue -> tail < queue -> len) {
      /* Enlarge queue if not full, insert at tail and resort towards head */
      kid = queue -> tail++;
   
      while (kid) {
         root = (kid - 1) >> 1;
      
         /* Compare with parent and swap if necessary, else terminate */
         if (d2 > qn [root].dist2) {
            qn [kid].dist2 = qn [root].dist2;
            qn [kid].idx  = qn [root].idx;
            kid = root;
         }
         else break;
      }
      
      /* Assign tail position as linear index into record buffa 
       * queue -> nnRec and append record */
      qn [kid].dist2 = d2;
      qn [kid].idx   = queue -> tail - 1;
      memcpy(OOC_GetNearest(queue, qn [kid].idx), rec, queue -> recSize);
   }
   else if (d2 < qn [0].dist2) {
      /* Queue full and new record closer than maximum at head; replace head
       * and resort towards tail */
      root = 0;
      rootIdx = qn [root].idx;
      
      while ((kid = (root << 1) + 1) < queue -> tail) {
         /* Compare with larger kid & swap if necessary, else terminate */
         if ((kid2 = (kid + 1)) < queue -> tail && 
             qn [kid2].dist2 > qn [kid].dist2)
            kid = kid2;
         
         if (d2 < qn [kid].dist2) {
            qn [root].dist2 = qn [kid].dist2;
            qn [root].idx  = qn [kid].idx;
         }
         else break;
         
         root = kid;
      }
      
      /* Reassign head's previous buffa index and overwrite corresponding
       * record */
      qn [root].dist2   = d2;
      qn [root].idx     = rootIdx;
      memcpy(OOC_GetNearest(queue, qn [root].idx), rec, queue -> recSize);
      
      /* Update SQUARED maximum distance from head node */
      d2max = qn [0].dist2;
   }

   return d2max;
}



int OOC_InitNearest (OOC_SearchQueue *squeue, 
                     unsigned len, unsigned recSize)
{
   squeue -> len = len;
   squeue -> recSize = recSize;
   squeue -> node = calloc(len, sizeof(OOC_SearchQueueNode));
   squeue -> nnRec = calloc(len, recSize); 
   if (!squeue -> node || !squeue -> nnRec) {
      perror("OOC_InitNearest: failed search queue allocation");
      return -1;
   }
   
   return 0;
}



void *OOC_GetNearest (const OOC_SearchQueue *squeue, unsigned idx)
{
   return squeue -> nnRec + idx * squeue -> recSize;
}



static float OOC_BBoxDist2 (const FVECT bbOrg, RREAL bbSiz, const FVECT key)
/* Return minimum *SQUARED* distance between key and bounding box defined by
 * bbOrg and bbSiz; a distance of 0 implies the key lies INSIDE the bbox */
{
   /* Explicit comparison with bbox corners */
   int   i;
   FVECT d;
   
   for (i = 0; i < 3; i++) {
      d [i] = key [i] - bbOrg [i];
      d [i] = d [i] < 0 ? -d [i] : d [i] - bbSiz;
      
      /* Set to 0 if inside bbox */
      if (d [i] < 0)
         d [i] = 0;
   }
   
   return DOT(d, d);
}



float OOC_FindNearest (OOC_Octree *oct, OOC_Node *node, OOC_DataIdx dataIdx,
                       const FVECT org, float size, const FVECT key, 
                       const OOC_SearchFilter *filter, 
                       OOC_SearchQueue *queue, float maxDist2)
{
   const float kidSize = size * 0.5;
   unsigned    i, kid, kid0;
   float       d2;
   char        rec [oct -> recSize];
   FVECT       kidOrg;
   OOC_DataIdx kidDataIdx, recIdx;
   OOC_Node    *kidNode;
   
   /* Start with suboctant closest to key */
   for (kid0 = 0, i = 0; i < 3; i++)
      kid0 |= (key [i] > org [i] + kidSize) << i;
      
   for (i = 0; i < 8; i++) {
      kid = kid0 ^ i;
      kidNode = node;
      kidDataIdx = dataIdx + OOC_GetKid(oct, &kidNode, kid);
      
      /* Prune empty suboctant */
      if ((!kidNode && !OOC_ISLEAF(node)) ||
          (OOC_ISLEAF(node) && !node -> leaf.num [kid]))
         continue;
         
      /* Set up suboctant */
      VCOPY(kidOrg, org);
      OOC_OCTORIGIN(kidOrg, kid, kidSize);
    
      /* Prune suboctant if not overlapped by maxDist2 */
      if (OOC_BBoxDist2(kidOrg, kidSize, key) > maxDist2)
         continue;
         
      if (kidNode) {
         /* Internal node; recurse into non-empty suboctant */
         maxDist2 = OOC_FindNearest(oct, kidNode, kidDataIdx, kidOrg, 
                                    kidSize, key, filter, queue, maxDist2);
         if (maxDist2 < 0) 
            /* Bail out on error */
            break;
      }                
      else if (OOC_ISLEAF(node))
         /* Leaf node; do linear check of all records in suboctant */
         for (recIdx = kidDataIdx; 
              recIdx < kidDataIdx + node -> leaf.num [kid]; recIdx++) {
            if (OOC_GetData(oct, recIdx, rec))
               return -1;
               
            if (!filter || filter -> func(rec, filter -> data))
               /* Insert record in search queue SQUARED dist to key below
                * maxDist2 and passes filter */
               if ((d2 = dist2(key, oct -> key(rec))) < maxDist2) {
                  if ((d2 = OOC_PutNearest(queue, d2, rec)) >= 0)
                     /* Update maxDist2 if queue is full */
                     maxDist2 = d2;
#ifdef DEBUG_OOC_NN
                     if (OOC_SearchQueueCheck(queue, key, oct -> key, 0)) {
                        fprintf(stderr, "OOC_SearchPush: priority queue "
                                "inconsistency\n");
                        return -1;
                     }
#endif
               }
         }
   }
   
   return maxDist2;
}



float OOC_Find1Nearest (OOC_Octree *oct, OOC_Node *node, OOC_DataIdx dataIdx,
                        const FVECT org, float size, const FVECT key,  
                        const OOC_SearchFilter *filter, void *nnRec, 
                        float maxDist2)
{
   const float kidSize = size * 0.5;
   unsigned    i, kid, kid0;
   float       d2;
   char        rec [oct -> recSize];
   FVECT       kidOrg;
   OOC_DataIdx kidDataIdx, recIdx;
   OOC_Node    *kidNode;
   
   /* Start with suboctant closest to key */
   for (kid0 = 0, i = 0; i < 3; i++)
      kid0 |= (key [i] > org [i] + kidSize) << i;
      
   for (i = 0; i < 8; i++) {
      kid = kid0 ^ i;
      kidNode = node;
      kidDataIdx = dataIdx + OOC_GetKid(oct, &kidNode, kid);
      
      /* Prune empty suboctant */
      if ((!kidNode && !OOC_ISLEAF(node)) ||
          (OOC_ISLEAF(node) && !node -> leaf.num [kid]))
         continue;
         
      /* Set up suboctant */
      VCOPY(kidOrg, org);
      OOC_OCTORIGIN(kidOrg, kid, kidSize);
    
      /* Prune suboctant if not overlapped by maxDist2 */
      if (OOC_BBoxDist2(kidOrg, kidSize, key) > maxDist2)
         continue;
         
      if (kidNode) {
         /* Internal node; recurse into non-empty suboctant */
         maxDist2 = OOC_Find1Nearest(oct, kidNode, kidDataIdx, kidOrg, 
                                     kidSize, key, filter, nnRec, maxDist2);
         if (maxDist2 < 0) 
            /* Bail out on error */
            break;                                     
      }
      else if (OOC_ISLEAF(node))
         /* Leaf node; do linear check of all records in suboctant */
         for (recIdx = kidDataIdx; 
              recIdx < kidDataIdx + node -> leaf.num [kid]; recIdx++) {
            if (OOC_GetData(oct, recIdx, rec))
               return -1;
               
            if (!filter || filter -> func(rec, filter -> data))
               /* Update closest record and max SQUARED dist to key if it
                * passes filter */
               if ((d2 = dist2(key, oct -> key(rec))) < maxDist2) {
                  memcpy(nnRec, rec, oct -> recSize);
                  maxDist2 = d2;
               }
         }
   }
   
   return maxDist2;
}

#endif /* NIX / PMAP_OOC */
