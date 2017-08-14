#ifndef lint
static const char RCSid[] = "$Id: oocbuild.c,v 2.4 2017/08/14 21:12:10 rschregle Exp $";
#endif


/* 
   =======================================================================
   Routines for building out-of-core octree data structure
   
   Adapted from:  Kontkanen J., Tabellion E. and Overbeck R.S.,
                  "Coherent Out-of-Core Point-Based Global Illumination",
                  EGSR 2011.   
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   =======================================================================
   
   $Id: oocbuild.c,v 2.4 2017/08/14 21:12:10 rschregle Exp $
*/


#if !defined(_WIN32) && !defined(_WIN64) || defined(PMAP_OOC)
/* No Windoze support for now */

#include "oococt.h"
#include "oocsort.h"
#include <stdlib.h>
#include <string.h>


/* Test for empty/full input queue, return pointer to head/tail */
#define QueueFull(q)    ((q) -> len == (q) -> cap)
#define QueueEmpty(q)   (!(q) -> len)
#define QueueHead(q)    ((q) -> data + (q) -> head * (q) -> recSize)
#define QueueTail(q)    ((q) -> data + \
                           ((q)->head + (q)->len-1) % (q)->cap * (q)->recSize)



/* Input queue for bottom-up octree construction */
typedef struct {
   void     *data;
   unsigned head, len, cap, recSize;   /* Queue head, length (from head),
                                          capacity and record size */
   FILE     *in;                       /* Input file for data records */
} OOC_BuildQueue;


static OOC_BuildQueue *QueueInit (OOC_BuildQueue *q, unsigned recSize,
                                  unsigned capacity)
/* Initialise queue of #capacity records of size recSize each; returns queue
 * pointer or NULL if failed. */
{
   if (!(q && (q -> data = calloc(capacity, recSize))))
      return NULL;
   
   q -> cap = capacity;
   q -> recSize = recSize;
   q -> head = q -> len = 0;
   
   return q;
}
   
   
static int QueuePush (OOC_BuildQueue *q, const void *rec)
/* Append record to queue tail; return new queue length or -1 on failure */
{
   if (!q || !rec || QueueFull(q))
      return -1;
   
   ++q->len;
   memcpy(QueueTail(q), rec, q -> recSize);

   return q -> len;
}


static int QueuePop (OOC_BuildQueue *q, void *rec)
/* Remove record from queue head and return in rec if not NULL; return new
 * queue length or -1 on failure */
{
   if (!q || QueueEmpty(q))
      return -1;
      
   /* Return head if rec != NULL */
   if (rec)
      memcpy(rec, QueueHead(q), q -> recSize);      
      
   q -> head = (q -> head + 1) % q -> cap;
   
   return --q -> len;
}   


static int QueueFill (OOC_BuildQueue *q)
/* Read records from q -> in until the queue is full; return queue
 * length or -1 on failure */
{
   static void *rec = NULL;
   
   if (!rec && !(rec = malloc(q -> recSize)))
      return -1;

   while (!QueueFull(q) && !feof(q -> in) && 
          fread(rec, q -> recSize, 1, q -> in))
      QueuePush(q, rec);
      
   return q -> len;
}



static OOC_DataIdx OOC_BuildRecurse (OOC_Octree *oct, OOC_Node* node, 
                                     FVECT org, RREAL size, unsigned depth,
                                     OOC_BuildQueue *queue)
/* Recursive part of OOC_Build(); insert records from input queue into
 * octree node and subdivide into kids if necessary. Returns number of
 * records in subtree or OOC_DATAIDX_ERR if failed. */
{
   int            k;
   const RREAL    kidSize = size * 0.5;
   
   if (!oct || !node)
      return OOC_DATAIDX_ERR;

   if (QueueEmpty(queue) || 
       !OOC_InBBox(oct, org, size, oct -> key(QueueHead(queue))))
      /* Input exhausted or queue head outside node */
      return 0;
   
   if (QueueFull(queue) && depth < oct -> maxDepth && 
       OOC_InBBox(oct, org, size, oct -> key(QueueTail(queue)))) {
      /*************************** SUBDIVIDE NODE *************************
       * At least leafMax + 1 records (since the queue is full) lie inside
       * the current node's bounding box, and maxDepth has not been reached
       * ==> subdivide this node.
       * (Note it suffices to only test the queue tail against the bounding
       * box, as the records are in Z-order)
       ********************************************************************/
      OOC_Node    kid [8];
      OOC_DataIdx dataCnt;
      FVECT       kidOrg;

#ifdef DEBUG_OOC_BUILD
FVECT       key;
unsigned    k2 = 0;
#endif      
      
      /* We recurse on the nonempty kids first, then finalise their nodes so
       * they are ordered consecutively, since the parent only indexes the 1st kid */
      for (k = 0; k < 8; k++) {
         /* Clear kid node and get its octant origin */
         OOC_CLEARNODE(kid + k);

         VCOPY(kidOrg, org);
         OOC_OCTORIGIN(kidOrg, k, kidSize);
         
         /* Recurse on kid and check for early bailout */
         if (OOC_BuildRecurse(oct, kid + k, kidOrg, kidSize, 
                              depth + 1, queue) == OOC_DATAIDX_ERR)
            return OOC_DATAIDX_ERR;

#ifdef DEBUG_OOC_BUILD
if (!QueueEmpty(queue)) {
VCOPY(key, oct -> key(QueueHead(queue)));
k2 = OOC_Branch(oct, org, kidSize, key, NULL);
if (OOC_InBBox(oct, org, size, key) && k2 < k) {
   fprintf(stderr, 
           "OOC_BuildRecurse, node subdiv: unsorted key [%f, %f, %f] with "
           "octant %d (last %d with bbox [%f-%f, %f-%f, %f-%f])\n",
           key [0], key [1], key [2], k2, k, 
           kidOrg [0], kidOrg [0] + kidSize, kidOrg [1], kidOrg [1] + kidSize,
           kidOrg [2], kidOrg [2] + kidSize);
}}
#endif
      }

      /* Now finalise consecutive kid nodes, skipping empty ones */
      for (k = 0; k < 8; k++)
         if ((dataCnt = OOC_DATACNT(kid + k))) {         
            /* Nonzero kid ==> allocate and set node */
            if (!NewNode(oct)) {
               fputs("OOC_BuildRecurse: failed to allocate new node\n",
                     stderr);
               return OOC_DATAIDX_ERR;
            }
            OOC_SETROOT(oct, kid + k);
         
            /* Sum kid's data count to parent's and check for overflow */
            if ((dataCnt += node -> node.num) <= OOC_DATAIDX_MAX)
               node -> node.num = dataCnt;
            else {
               fputs("OOC_BuildRecurse: data count overflow in node\n",
                     stderr);
               return OOC_DATAIDX_ERR;
            }            
         
            /* Set kid index in parent (if first kid) and corresponding
             * branch bit. The kid is the most recent node and thus at the
             * end of the node array, which coincides with the current
             * subtree root */
            if (!node -> node.branch)
               node -> node.kid = OOC_ROOTIDX(oct);            
               
            OOC_SETBRANCH(node, k);
         }
   }
   else {
      /****************************** MAKE LEAF ****************************
       * Queue contains no more than leafMax records, queue tail lies
       * outside node's bounding box, or maxDepth reached 
       * ==> make this node a leaf.
       *********************************************************************/
      RREAL          *key;

#ifdef DEBUG_OOC_BUILD
OOC_MortonIdx      zIdx, lastzIdx = 0;
FVECT             /* key, */
                  lastKey, kidOrg;
unsigned          lastk = 0;
#endif
      
      /* Mark as leaf (note it's been cleared by the parent call) */
      OOC_SETLEAF(node);

      while (!QueueEmpty(queue) && 
             OOC_InBBox(oct, org, size, (key = oct->key(QueueHead(queue))))) {
         /* Record lies inside leaf; increment data counter for octant 
          * containing record. */

         if ((k = OOC_Branch(oct, org, kidSize, key, NULL)) < 0) {
            /* Shouldn't happen, as key tested within bbox above */
            fprintf(stderr, "OOC_BuildRecurse: buggered Morton code, "
                    "disruption in space-time continuum?\n");
            return OOC_DATAIDX_ERR;
         }
         
         if (node -> leaf.num [k] == OOC_OCTCNT_MAX) {
            /* Currently we're buggered here; merge records instead? */
            fprintf(stderr, "OOC_BuildRecurse: data count overflow in "
                    "leaf: depth = %d, count = %d\n", 
                    depth, node -> leaf.num [k]);
            return OOC_DATAIDX_ERR;
         }
         
         ++node -> leaf.num [k];

#ifdef DEBUG_OOC_BUILD
/* VCOPY(key, oct -> key(QueueHead(queue))); */
if ((zIdx = OOC_KEY2MORTON(key, oct)) < lastzIdx) {
   fprintf(stderr, "OOC_BuildRecurse, make leaf: unsorted zIdx %020ld for "
           "key [%f, %f, %f] (previous zIdx %020ld for "
           "key [%f, %f, %f]\n", zIdx, key [0], key [1], key [2],
           lastzIdx, lastKey [0], lastKey [1], lastKey [2]);
}

VCOPY(kidOrg, org);
OOC_OCTORIGIN(kidOrg, k, kidSize);
if (k < lastk || zIdx < lastzIdx) {
   fprintf(stderr, 
           "OOC_BuildRecurse, make leaf: unsorted octant %d (last %d) with "
           "bbox [%f-%f, %f-%f, %f-%f] for key [%f, %f, %f] with zIdx %020ld "
           "(last [%f, %f, %f], zIdx %020ld)\n",
           k, lastk, kidOrg [0], kidOrg [0] + kidSize, 
           kidOrg [1], kidOrg [1] + kidSize, kidOrg [2], kidOrg [2] + kidSize, 
           key [0], key [1], key [2], zIdx, 
           lastKey [0], lastKey [1], lastKey [2], lastzIdx);
}           
lastk = k;
lastzIdx = zIdx;
VCOPY(lastKey, key);
#endif

         /* Remove record from queue */
         QueuePop(queue, NULL);          
      }

      /* Refill queue for next node(s) */
      if (QueueFill(queue) < 0) {
         fputs("OOC_Build: failed input queue fill\n", stderr);
         return OOC_DATAIDX_ERR;
      }
   }
   
   return OOC_DATACNT(node);
}                                   

   
   
OOC_Octree *OOC_Build (OOC_Octree *oct, unsigned leafMax, unsigned maxDepth)
/* Bottom-up construction of out-of-core octree in postorder traversal. The
 * octree oct is assumed to be initialised with its origin (oct -> org),
 * size (oct -> size), key callback (oct -> key), and its associated leaf
 * file (oct -> leafFile).
 
 * Records are read from the leafFile and assumed to be sorted in Z-order,
 * which defines an octree leaf ordering.  Leaves (terminal nodes) are
 * constructed such that they contain <= leafMax records and have a maximum
 * depth of maxDepth.
 
 * Note that the following limits apply:
 *    leafMax  <= OOC_OCTCNT_MAX     (see oococt.h)
 *    maxDepth <= OOC_MORTON_BITS    (see oocsort.h)
 
 * The maxDepth limit arises from the fact that the Z-ordering has a limited
 * resolution and will map node coordinates beyond a depth of
 * OOC_MORTON_BITS to the same Z-index, causing nodes to be potentially read
 * out of sequence and winding up in the wrong octree nodes.
  
 * On success, the octree pointer oct is returned, with the constructed
 * nodes in oct -> nodes, and the node count in oct -> numNodes.  On
 * failure, NULL is returned.  */
{
   OOC_BuildQueue queue;
   OOC_Node       root;
   
   if (!oct || !oct -> size) {
      fputs("OOC_Build: octree not initialised", stderr);
      return NULL;
   }
   
   if (!oct -> leafFile) {
      fputs("OOC_Build: empty leaf file", stderr);
      return NULL;
   }
      
   oct -> leafMax = leafMax;
   oct -> maxDepth = maxDepth;
   queue.in = oct -> leafFile;
   rewind(queue.in);
   
   /* Init queue and fill from leaf file */
   if (!QueueInit(&queue, oct -> recSize, leafMax + 1) || 
       QueueFill(&queue) < 0) {
      fputs("OOC_Build: failed input queue init\n", stderr);
      return NULL;
   }

   /* Clear octree root and recurse */
   OOC_CLEARNODE(&root);
   if (OOC_BuildRecurse(oct, &root, oct -> org, oct -> size, 0, &queue) ==
       OOC_DATAIDX_ERR)
      return NULL;
    
   /* Finalise octree root */
   if (!NewNode(oct)) {
      fputs("OOC_Build: failed to allocate octree root\n", stderr);
      return NULL;
   }      
   OOC_SETROOT(oct, &root);
   /* Passing OOC_ROOT(oct) avoids annoying compiler warnings about (&root)
    * always evaluating to true when calling OOC_DATAIDX() */
   oct -> numData = OOC_DATACNT(OOC_ROOT(oct));
   
   return oct;
}

#endif /* NIX / PMAP_OOC */
