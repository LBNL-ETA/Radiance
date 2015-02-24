/* 
   ==================================================================
   N-way hybrid out-of-core merge sort
   Sorts N blocks internally using quicksort, followed by N-way
   merge using priority queue.
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts   
   ==================================================================
   
   $Id $
*/


#include "oocsort.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* Priority queue node */
typedef struct {
   OOC_Sort_Key pri;   /* Record's priority (sort key) */
   unsigned blk;           /* Block containing record */
} OOC_Sort_PQNode;

typedef struct {
   OOC_Sort_PQNode *node;
   unsigned len, tail;
} OOC_Sort_PQueue;

/* Block descriptor */
typedef struct {
   FILE *file;                /* Temporary file */
   char *buf, *head, *tail;   /* Associated I/O buffa, 
                                 head/tail pointers */
} OOC_Sort_Block;


/* Record priority evaluation function */
static OOC_Sort_Key (*pri)(const void *);


static int OOC_Sort_Compare (const void *rec1, const void *rec2)
/* Comparison function for internal sorting */
{
   const OOC_Sort_Key pri1 = pri(rec1), pri2 = pri(rec2);
   
   if (!rec1 || !rec2)
      return 0;
   
   if (pri1 < pri2)
      return -1;
   else if (pri1 > pri2)
      return 1;
   else 
      return 0;
}


static char *OOC_Sort_Read (OOC_Sort_Block *blk, unsigned recSize)
/* Read next record from specified block; return pointer to buffa entry
   or NULL on error */
{
   char *rec = NULL;
         
   if (blk -> head >= blk -> tail) {
      /* Input buffa empty; read next block (potentially truncated if last),
       * where tail marks end of buffa */ 
      const unsigned long bufSize = blk -> tail - blk -> buf;
    
      if (feof(blk -> file))
         return NULL;

      blk -> head = blk -> buf;
      blk -> tail = blk -> buf + fread(blk -> buf, 1, bufSize, blk -> file);
                    
      if (blk -> tail == blk -> head)
         return NULL;
   }
   
   rec = blk -> head;
   blk -> head += recSize;
   
   return rec;
}


static char *OOC_Sort_Peek (OOC_Sort_Block *blk, unsigned recSize)
/* Return next record from specified block WITHOUT removing from buffa */
{
   char *rec = NULL;
   
   if (rec = OOC_Sort_Read(blk, recSize)) {
      /* Restore buffa head */
      blk -> head -= recSize;
   }
   
   return rec;
}


static char *OOC_Sort_Write (OOC_Sort_Block *blk, unsigned recSize,
                             const char *rec)
/* Output record to specified block and return pointer to buffa entry
   or NULL on error */
{
   char *res = NULL;
   
   if (blk -> head >= blk -> tail) {
      /* Flush output buffa (tail marks end of buffa) */
      const unsigned long bufSize = blk -> tail - blk -> buf;

      blk -> head = blk -> buf;
      if (fwrite(blk -> buf, 1, bufSize, blk -> file) != bufSize)
         return NULL;
   }
   
   if (!rec)
      return NULL;
      
   memcpy(blk -> head, rec, recSize);
   res = blk -> head;
   blk -> head += recSize;
   
   return res;
}


#ifdef DEBUG
   static int OOC_Sort_PQCheck (OOC_Sort_PQueue *pq, unsigned root)
   /* Priority queue sanity check */
   {
      unsigned kid;
      
      if (root < pq -> tail) {
         if ((kid = (root << 1) + 1) < pq -> tail)
            if (pq -> node [kid].pri < pq -> node [root].pri)
               return -1;
            else return OOC_Sort_PQCheck(pq, kid);
            
         if ((kid = kid + 1) < pq -> tail)
            if (pq -> node [kid].pri < pq -> node [root].pri)
               return -1;
            else return OOC_Sort_PQCheck(pq, kid);
         
      }
      
      return 0;
   }
#endif


static int OOC_Sort_Push (OOC_Sort_PQueue *pq, OOC_Sort_Key pri,
                          unsigned blk)
/* Insert specified block index into priority queue; return block index 
 * or -1 if queue is full */
{
   OOC_Sort_PQNode *pqn = pq -> node;
   unsigned kid, root;
   
   if (pq -> tail >= pq -> len)
      /* Queue full */
      return -1;
   
   /* Append node at tail and re-sort */
   kid = pq -> tail++;
   
   while (kid) {
      root = kid - 1 >> 1;
      
      /* Compare with parent node and swap if necessary, else terminate */
      if (pri < pqn [root].pri) {
         pqn [kid].pri = pqn [root].pri;
         pqn [kid].blk = pqn [root].blk;
         kid = root;
      }
      else break;
   }
      
   pqn [kid].pri = pri;
   pqn [kid].blk = blk;   

#ifdef DEBUG
   if (OOC_Sort_PQCheck(pq, 0) < 0) {
      fprintf(stderr, "OOC_Sort: priority queue inconsistency\n");
      return -1;
   }
#endif
   
   return blk;
}


static int OOC_Sort_Pop (OOC_Sort_PQueue *pq)
/* Remove head of priority queue and return its block index 
 * or -1 if queue is empty */
{
   OOC_Sort_PQNode *pqn = pq -> node;
   OOC_Sort_Key pri;
   unsigned kid, kid2, root = 0, blk, res;
   
   if (!pq -> tail)
      /* Queue empty */
      return -1;
   
   res = pqn -> blk;
   pri = pqn [--pq -> tail].pri;
   blk = pqn [pq -> tail].blk;

   /* Replace head with tail node and re-sort */
   while ((kid = (root << 1) + 1) < pq -> tail) {
      /* Compare with with smaller kid and swap if necessary, else
       * terminate */
      if ((kid2 = (kid + 1)) < pq -> tail && pqn [kid2].pri < pqn [kid].pri)
         kid = kid2;
            
      if (pri > pqn [kid].pri) {
         pqn [root].pri = pqn [kid].pri;
         pqn [root].blk = pqn [kid].blk;
      }
      else break;
         
      root = kid;
   }
   
   pqn [root].pri = pri;
   pqn [root].blk = blk;

#ifdef DEBUG
   if (OOC_Sort_PQCheck(pq, 0) < 0) {
      fprintf(stderr, "OOC_Sort: priority queue inconsistency\n");
      return -1;
   }
#endif
   
   return res;
}
 
#if 0 
int OOC_Sort (const char *inFile, const char *outFile,
              unsigned long blkSize, unsigned recSize,
              OOC_Sort_Key (*priority)(const void *))
/* Sort records in inFile and output to outFile by (a) internally
 * quicksorting block of blkSize bytes at a time, then (b) merging them
 * via a priority queue. RecSize specifies the size in bytes of each data
 * record. The priority() callback evaluates a record's priority and must be
 * supplied by the caller. */
{
   FILE *in = NULL;
   OOC_Sort_PQueue pqueue;                         /* Priority queue */
   OOC_Sort_Block *iBlk = NULL, oBlk;              /* Block descriptors */ 
   char *rec, *sortBuf = NULL;                     /* Internal sort buffa */
   unsigned bufSize, numBlk;
   int b;
   
   /* Set record priority evaluation callback */                
   pri = priority;
   
   /* Round block and buffa size down to nearest multiple of record size */
   blkSize = (blkSize / recSize) * recSize;
   bufSize = (OOC_SORT_BUFSIZE / recSize) * recSize;

   if (!(in = fopen(inFile, "rb"))) {
      fprintf(stderr, "OOC_Sort: failure opening input file %s\n", inFile);
      return -1;
   }
   
   /* Get input file size and number of blocks (rounded up to nearest
    * multiple of block size) */
   fseek(in, 0, SEEK_END);
   numBlk = (ftell(in) + blkSize - 1) / blkSize;
   rewind(in);
#else

int OOC_Sort (const char *inFile, const char *outFile,
              unsigned numBlk, unsigned recSize,
              OOC_Sort_Key (*priority)(const void *))
/* Sort records in inFile and output to outFile by (a) internally
 * quicksorting numBlk blocks at a time, then (b) merging them via a priority
 * queue.  RecSize specifies the size in bytes of each data record.  The
 * priority() callback evaluates a record's priority (ordinal index) and
 * must be supplied by the caller.  */
{
   FILE *in = NULL;
   OOC_Sort_PQueue pqueue;                         /* Priority queue */
   OOC_Sort_Block *iBlk = NULL, oBlk;              /* Block descriptors */ 
   char *rec, *sortBuf = NULL;                     /* Internal sort buffa */
   unsigned long blkSize;
   unsigned bufSize;
   int b;
   
   /* Set record priority evaluation callback */                
   pri = priority;
   
   /* Round buffa size down to nearest multiple of record size */
   bufSize = (OOC_SORT_BUFSIZE / recSize) * recSize;

   if (!(in = fopen(inFile, "rb"))) {
      fprintf(stderr, "OOC_Sort: failure opening input file %s\n", inFile);
      return -1;
   }
   
   /* Get input file size and block size (in number of records) rounded up
    * to nearest multiple of number of blocks */
   fseek(in, 0, SEEK_END);
   blkSize = (ftell(in) / recSize + numBlk - 1) / numBlk;   
   rewind(in);
#endif

   /* Allocate buffa for internal sorting */
   if (!(sortBuf = malloc(blkSize))) {
      fprintf(stderr, "OOC_Sort: failure allocating internal sort buffer\n");
      return -1;
   }
   
   /* Allocate input blocks */
   if (!(iBlk = calloc(numBlk, sizeof(OOC_Sort_Block)))) {
      fprintf(stderr, "OOC_Sort: failure allocating input blocks\n");
      return -1;
   }
   
   /* ===================================================================
    * Pass 1: Internal quicksort of each input block in sortBuf
    * =================================================================== */

   for (b = 0; b < numBlk; b++) {
      unsigned long numRec;
      
      /* Open temporary file associated with block */
      if (!(iBlk [b].file = tmpfile())) {
         /* fprintf(stderr, "OOC_Sort: failure opening block file\n"); */
         perror("OOC_Sort: failure opening block file");
         return -1;
      }
      
      if (feof(in)) {
         fprintf(stderr, "OOC_Sort: unexpected end of input file %s\n",
                 inFile);
         return -1;
      }
      
      /* Read next block (potentially truncated if last) */
      if (!(numRec = fread(sortBuf, 1, blkSize, in) / recSize)) {
         fprintf(stderr, "OOC_Sort: error reading from input file %s\n",
                 inFile);
         return -1;
      }
         
      /* Quicksort block internally and write out to temp file */
      qsort(sortBuf, numRec, recSize, OOC_Sort_Compare);
      
      if (fwrite(sortBuf, recSize, numRec, iBlk [b].file) != numRec) {
         fprintf(stderr, "OOC_Sort: error writing to block file\n");
         return -1;
      }
      
      rewind(iBlk [b].file);
   }

   /* Close input file and free sort buffa */
   fclose(in);
   free(sortBuf);

   /* Allocate priority queue with numBlk nodes */   
   pqueue.tail = 0;
   pqueue.len = numBlk;
   if (!(pqueue.node = calloc(numBlk, sizeof(OOC_Sort_PQNode)))) {
      fprintf(stderr, "OOC_Sort: failure allocating priority queue\n");
      return -1;
   }
   
   /* Prepare for pass 2: allocate buffa for each input block and initialise
    * fill priority queue with single record from each block */
   for (b = 0; b < numBlk; b++) {
      if (!(iBlk [b].buf = malloc(bufSize))) {
         fprintf(stderr, "OOC_Sort: failure allocating input buffer\n");
         return -1;
      }
      
      /* Peek first record in block file without modifying buffa */
      iBlk [b].head = iBlk [b].tail = iBlk [b].buf + bufSize;
      if (!(rec = OOC_Sort_Peek(iBlk + b, recSize))) {
         fprintf(stderr, "OOC_Sort: error reading from block file\n");
         return -1;
      }
      
      /* Insert record into priority queue */
      if (OOC_Sort_Push(&pqueue, priority(rec), b) < 0) {
         fprintf(stderr, "OOC_Sort: failed priority queue insertion\n");
         return -1;
      }
   }
   
   /* Allocate output buffa and open output file */
   if (!(oBlk.file = fopen(outFile, "wb"))) {
      fprintf(stderr, "OOC_Sort: failure opening output file %s\n", outFile);
      return -1;
   }
   
   if (!(oBlk.head = oBlk.buf = malloc(bufSize))) {
      fprintf(stderr, "OOC_Sort: failure allocating output buffer\n");
      return -1;
   }
   
   /* tail marks end of output buffa */
   oBlk.tail = oBlk.buf + bufSize;
   
   /* ===================================================================
    * Pass 2: External merge of all blocks using priority queue
    * =================================================================== */
      
   do {
      /* Get next node in priority queue, read next record in corresponding
       * block, and send to output */
      b = OOC_Sort_Pop(&pqueue);
            
      if (b >= 0) {
         /* Priority queue non-empty */
         if (!(rec = OOC_Sort_Read(iBlk + b, recSize))) {
            /* Corresponding record should still be in the buffa, so EOF
             * should not happen */
            fprintf(stderr, "OOC_Sort: unexpected EOF in block file\n");
            return -1;
         }
         
         if (!OOC_Sort_Write(&oBlk, recSize, rec)) {
            fprintf(stderr, "OOC_Sort; error writing to output file %s\n",
                    outFile);
            return -1;
         }
         
         /* Peek next record from same block and insert in priority queue */
         if (rec = OOC_Sort_Peek(iBlk + b, recSize))
            /* Buffa not exhausted */
            if (OOC_Sort_Push(&pqueue, priority(rec), b) < 0) {
               fprintf(stderr, "OOC_Sort: failed priority queue insertion\n");
               return -1;
            }
      }
   } while (b >= 0);
   
   /* Priority queue now empty; flush output buffa & close file */
   oBlk.tail = oBlk.head;
   OOC_Sort_Write(&oBlk, 0, NULL);
   fclose(oBlk.file);
   
   for (b = 0; b < numBlk; b++) {
      fclose(iBlk [b].file);
      free(iBlk [b].buf);
   }

   free(iBlk);   
   free(pqueue.node);
   
   return 0;           
}   