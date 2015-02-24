/* 
   =======================================================================
   Routines for building out-of-core octree data structure
   
   Adapted from:  Kontkanen J., Tabellion E. and Overbeck R.S.,
                  "Coherent Out-of-Core Point-Based Global Illumination",
                  EGSR 2011.   
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts   
   =======================================================================
   
   $Id: oocbuild.c,v 2.1 2015/02/24 19:39:26 greg Exp $
*/


#include "oococt.h"
#include "oocsort.h"


/* Test for empty/full input queue, return pointer to head/tail */
#define QueueFull(q)    ((q) -> len == (q) -> cap)
#define QueueEmpty(q)   (!(q) -> len)
#define QueueHead(q)    ((q) -> data + (q) -> head * (q) -> recSize)
#define QueueTail(q)    ((q) -> data + \
                         (q) -> head + (q) -> len - 1) * (q) -> recSize)


/* Input queue from heap */
type struct {
   void *data = NULL;
   /* Queue head, length (from head), capacity and record size */
   unsigned head, len, cap, recSize;
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
/* Append record to queue tail; returns new queue length or -1 on failure */
{
   int tail;
   
   if (!q || !rec)
      return -1;
      
   if (q -> len >= q -> cap)
      /* Queue full */
      return -1;
   
   tail = (q -> head + q-> len++) % q -> cap;
   memcpy(q -> data + tail * q -> recSize, rec, q -> recSize);
   
   return q -> len;
}


static int QueuePop (OOC_BuildQueue *q, void *rec)
/* Remove record from queue head; returns new queue length or -1 on failure */
{
   if (!q || !rec)
      return -1;
      
   if (!q -> len)
      /* Queue empty */
      return -1;
      
   memcpy(rev, q -> data + q -> head * q -> recSize, q -> recSize);
   q -> head = (q -> head + 1) % q -> cap;
   
   return --q -> len;
}   


static OOC_Idx *OOC_BuildRecurse (OOC_Octree *oct, OOC_Node *node, 
                                  FVECT org, RREAL size, unsigned depth,
                                  OOC_Queue *queue, FILE *heapFile, 
                                  FILE *nodeFile, FILE *leafFile)
/* Recursive part of OOC_Build(); insert records from input queue into
 * octree node and subdivide if necessary. Returns number of records in
 * subtree. */
{
   if (QueueEmpty(queue) || !InNode(QueueHead(queue)))
      /* Input exhausted or queue head outside node */
      return 0;
      
   if (InNode(QueueTail(queue)) && QueueFull(queue) && 
       depth < oct -> maxDepth) {
       
   }
   else {
   
   }      
}                                   
   
   
OOC_Octree *OOC_Build (OOC_Octree *oct, unsigned recSize, unsigned leafMax,
                       unsigned maxDepth, const char *heapFileName, 
                       const char *nodeFileName, const char *leafFileName)
/* Build out-of-core octree from unsorted records of size recSize in
 * heapFileName, such that leaves contain <= leafMax records, except those
 * at maxDepth. Nodes and leaves are output separately to nodeFileName and leafFileName,
 * respectively. Returns octree pointer on success, else NULL. */
{
   OOC_BuildQueue queue;
   void *rec = malloc(recSize);
   FILE *heapFile = fopen(heapFileName, "rb"), 
        *nodeFile = fopen(nodeFileName, "wb"), 
        *leafFile = fopen(leafFileName, "wb");
   
   if (OOC_Init(oct, recSize, leafMax, maxDepth) < 0) {
      perror("OOC_Build: failed octree init");
      return NULL;
   }
      
   if (!heapFile) {
      perror("OOC_Build: failed opening heap file");
      return NULL;
   }
   
   if (!nodeFile || !leafFile) {
      perror("OOC_Build: failed opening node/leaf output file");
      return NULL;
   }
   
   /* Init and fill queue from heap file */
   if (!rec || !QueueInit(&queue, recSize, leafMax + 1)) {
      perror("OOC_Build: failed input queue init");
      return NULL;
   }
   
   do {
      if (feof(heapFile) || !fread(rec, recSize, 1, heapFile)) {
         perror("OOC_Build: failed reading from heap file");
         
      QueuePush(&queue, rec);
   } while (!QueueFull(&queue));
   
   if (!OOC_BuildRecurse(oct, oct -> root, queue, 
                         heapFile, nodeFile, leafFile))
      return NULL;
      
   fclose(heapFile);
   fclose(nodeFile);
   fclose(leafFile);
   
   return oct;
}
