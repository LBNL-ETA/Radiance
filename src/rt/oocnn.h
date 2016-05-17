/* 
   =========================================================================
   k-nearest neighbour lookup routines for out-of-core octree data structure
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   =========================================================================
   
   $Id: oocnn.h,v 2.1 2016/05/17 17:39:47 rschregle Exp $
*/



#ifndef OOC_SEARCH_H
   #define OOC_SEARCH_H
   
   #include "oococt.h"
   
   /* Priority queue node for octree lookups */
   typedef struct {
      float       dist2;   /* Record's *SQUARED* distance to query point */
      unsigned    idx;     /* Record's index in queue buffer */
   } OOC_SearchQueueNode;


   /* Priority queue for octree lookups. Note that we explicitly store the
    * NN candidates in a local in-core buffer nnRec for later retrieval
    * without incurring additional disk I/O. */
   typedef struct {
      OOC_SearchQueueNode  *node;
      unsigned             len, tail, recSize;
      void                 *nnRec;
   } OOC_SearchQueue;
   
   
   /* Filter for k-NN search; records are accepted for which func(rec, data)
    * returns nonzero */
   typedef struct {
      int   (*func)(void *rec, void *data);
      void  *data;
   } OOC_SearchFilter;
   

   
   float OOC_FindNearest (OOC_Octree *oct, OOC_Node *node, 
                          OOC_DataIdx dataIdx, const FVECT org, float size, 
                          const FVECT key, const OOC_SearchFilter *filter, 
                          OOC_SearchQueue *queue, float maxDist2);
   /* Find k nearest neighbours (where k = queue -> len) within a maximum
    * SQUARED distance of maxDist2 around key. Returns the corresponding
    * data records with their SQUARED distances in the search queue 
    * (queue -> node [0] .. queue -> node [queue -> tail - 1]), with the
    * furthest neighbour at the queue head and queue->tail <= queue->len.  
    *
    * This is a recursive function requiring params for the current node:
    * DataIdx is the data offset for records in the current node, which is
    * necessary for implied addressing into the leaf file. Org and size are
    * the origin and size of the current node's bounding box.
    *
    * An optional filter may be specified; if filter != NULL, filter->func()
    * is called for each potential nearest neighbour along with additional
    * data provided by the caller; see OOC_SearchFilter typedef above.  
    *
    * Return value is the SQUARED distance to furthest neighbour on success,
    * else -1 on failure. */          


   float OOC_Find1Nearest (OOC_Octree *oct, OOC_Node *node, 
                           OOC_DataIdx dataIdx, const FVECT org, float size, 
                           const FVECT key, const OOC_SearchFilter *filter, 
                           void *nnRec, float maxDist2);
    /* Find single nearest neighbour within max SQUARED distance maxDist2
     * around key and copy corresponding record in nnRec. This is an
     * optimised version of OOC_FindNearest() without a search queue */
    
    
    int OOC_InitNearest (OOC_SearchQueue *squeue, 
                         unsigned len, unsigned recSize);
    /* Initialise NN search queue of length len and local buffa for records
     * of size recSize; precondition to calling OOC_FindNearest() */
    
    
    void *OOC_GetNearest (const OOC_SearchQueue *queue, unsigned idx);
    /* Return pointer to record at pos idx in search queue after calling
     * OOC_FindNearest() */
#endif
