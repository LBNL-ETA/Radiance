#ifndef lint
static const char RCSid[] = "$Id: oocsort.c,v 2.5 2021/02/10 21:48:50 rschregle Exp $";
#endif


/* 
   =========================================================================
   N-way out-of-core merge sort for records with 3D keys.  Recursively
   subdivides input into N blocks until these are of sufficiently small size
   to be sorted in-core according to Z-order (Morton code), then N-way
   merging them out-of-core using a priority queue as the stack unwinds.
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com}) 
   (c) Lucerne University of Applied Sciences and Arts, 
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ==========================================================================
   
   $Id: oocsort.c,v 2.5 2021/02/10 21:48:50 rschregle Exp $
*/


#if !defined(_WIN32) && !defined(_WIN64) || defined(PMAP_OOC)
/* No Windoze support for now */

#include "oocsort.h"
#include "oocmorton.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>


/* Priority queue node */
typedef struct {
   OOC_MortonIdx     pri;        /* Record's priority (sort key) */
   unsigned          blk;        /* Block containing record */
} OOC_SortQueueNode;


/* Priority queue */
typedef struct {
   OOC_SortQueueNode *node;
   unsigned          len, tail;
} OOC_SortQueue;


/* Additional data for qsort() compare function. We resort to instancing
 * this as a global variable instead of passing it to the compare func via
 * qsort_r(), since the latter is a non-portable GNU extension. */
typedef struct {
   RREAL    *(*key)(const void*);   /* Callback to access 3D key in record */
   FVECT    bbOrg;                  /* Origin of bbox containing keys */
   RREAL    mortonScale;            /* Scale for generating Morton code */
} OOC_KeyData;

static OOC_KeyData keyData;



static int OOC_KeyCompare (const void *rec1, const void *rec2)
/* Comparison function for in-core quicksort */
{
   OOC_MortonIdx  pri1, pri2;
   
   if (!rec1 || !rec2)
      return 0;
      
   pri1 = OOC_Key2Morton(keyData.key(rec1), keyData.bbOrg,
                         keyData.mortonScale);
   pri2 = OOC_Key2Morton(keyData.key(rec2), keyData.bbOrg, 
                         keyData.mortonScale);
   
   if (pri1 < pri2)
      return -1;
   else if (pri1 > pri2)
      return 1;
   else 
      return 0;
}



static int OOC_SortRead (FILE *file, unsigned recSize, char *rec)
/* Read next record from file; return 0 and record in rec on success, 
 * else -1 */
{
   if (!rec || feof(file) || !fread(rec, recSize, 1, file))
      return -1;

   return 0;
}



static int OOC_SortPeek (FILE *file, unsigned recSize, char *rec)
/* Return next record from file WITHOUT advancing file position;
 * return 0 and record in rec on success, else -1 */
{
   const unsigned long filePos = ftell(file);

   if (OOC_SortRead(file, recSize, rec))
      return -1;
      
   /* Workaround; fseek(file, -recSize, SEEK_CUR) causes subsequent
    * fread()'s to fail until rewind() */
   rewind(file);    
   if (fseek(file, filePos, SEEK_SET) < 0)
      return -1;

   return 0;
}



static int OOC_SortWrite (FILE *file, unsigned recSize, const char *rec)
/* Output record to file; return 0 on success, else -1 */
{
   if (!rec || !fwrite(rec, recSize, 1, file))
      return -1;
   
   return 0;
}



#ifdef DEBUG_OOC_SORT
   static int OOC_SortQueueCheck (OOC_SortQueue *pq, unsigned root)
   /* Priority queue sanity check */
   {
      unsigned kid;
      
      if (root < pq -> tail) {
         if ((kid = (root << 1) + 1) < pq -> tail) {
            if (pq -> node [kid].pri < pq -> node [root].pri)
               return -1;
            else return OOC_SortQueueCheck(pq, kid);
         }
            
         if ((kid = kid + 1) < pq -> tail) {
            if (pq -> node [kid].pri < pq -> node [root].pri)
               return -1;
            else return OOC_SortQueueCheck(pq, kid);
         }
      }
      
      return 0;
   }
#endif



static int OOC_SortPush (OOC_SortQueue *pq, OOC_MortonIdx pri, unsigned blk)
/* Insert specified block index into priority queue; return block index 
 * or -1 if queue is full */
{
   OOC_SortQueueNode *pqn = pq -> node;
   unsigned          kid, root;
   
   if (pq -> tail >= pq -> len)
      /* Queue full */
      return -1;
   
   /* Append node at tail and re-sort */
   kid = pq -> tail++;
   
   while (kid) {
      root = (kid - 1) >> 1;
      
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

#ifdef DEBUG_OOC_SORT
   if (OOC_SortQueueCheck(pq, 0) < 0) {
      fprintf(stderr, "OOC_Sort: priority queue inconsistency\n");
      return -1;
   }
#endif
   
   return blk;
}



static int OOC_SortPop (OOC_SortQueue *pq)
/* Remove head of priority queue and return its block index 
   or -1 if queue empty */
{
   OOC_SortQueueNode *pqn = pq -> node;
   OOC_MortonIdx     pri;
   unsigned          kid, kid2, root = 0, blk, res;
   
   if (!pq -> tail)
      /* Queue empty */
      return -1;
   
   res = pqn -> blk;
   pri = pqn [--pq -> tail].pri;
   blk = pqn [pq -> tail].blk;

   /* Replace head with tail node and re-sort */
   while ((kid = (root << 1) + 1) < pq -> tail) {
      /* Compare with smaller kid and swap if necessary, else terminate */
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

#ifdef DEBUG_OOC_SORT
   if (OOC_SortQueueCheck(pq, 0) < 0) {
      fprintf(stderr, "OOC_Sort: priority queue inconsistency\n");
      return -1;
   }
#endif
   
   return res;
}



/* Active subprocess counter updated by parent process; must persist when
 * recursing into OOC_SortRecurse(), hence global */
static unsigned procCnt = 0;

static int OOC_SortRecurse (FILE *in, unsigned long blkLo, 
                            unsigned long blkHi, FILE *out, 
                            unsigned numBlk, unsigned long maxBlkSize, 
                            unsigned numProc, unsigned recSize, 
                            char *sortBuf, OOC_SortQueue *pqueue, 
                            const OOC_KeyData *keyData)
/* Recursive part of OOC_Sort(). Reads block of records from input file
 * within the interval [blkLo, blkHi] and divides it into numBlk blocks
 * until the size (in bytes) does not exceed maxBlkSize, then quicksorts
 * each block into a temporary file.  These files are then mergesorted via a
 * priority queue to the output file as the stack unwinds.  NOTE: Critical
 * sections are prepended with '!!!'
 *
 * Parameters are as follows:
 * in          Input file containing unsorted records (assumed to be open)
 * blkLo       Start of current block in input file, in number of records
 * blkHi       End of current block in input file, in number of records
 * out         Output file containing sorted records (assumed to be open)
 * numBlk      Number of blocks to divide into / merge from
 * maxBlkSize  Max block size and size of in-core sort buffer, in bytes
 * numProc     Number of parallel processes for in-core sort
 * recSize     Size of input records in bytes
 * sortBuf     Preallocated in-core sort buffer of size maxBlkSize
 * pqueue      Preallocated priority queue of length numBlk for block merge
 * keyData     Aggregate data for Morton code generation and comparison
 */
{
   const unsigned long  blkLen = blkHi - blkLo + 1, 
                        blkSize = blkLen * recSize;
   int                  stat;
   
   if (!blkLen || blkHi < blkLo)
      return 0;
      
   if (blkSize > maxBlkSize) {
      unsigned long     blkLo2 = blkLo, blkHi2 = blkLo, blkSize2;
      const double      blkLen2 = (double)blkLen / numBlk;
      FILE              *blkFile [numBlk];         /* Violates C89!  */
      char              rec [recSize];             /* Ditto          */
      OOC_MortonIdx     pri;
      int               b, pid;
#ifdef DEBUG_OOC_SORT      
      unsigned long     pqCnt = 0;
#endif      
      
      /* ======================================================
       * Block too large for in-core sort -> divide into numBlk
       * subblocks and recurse
       * ====================================================== */

#ifdef DEBUG_OOC_SORT
      fprintf(stderr, "OOC_Sort: splitting block [%lu - %lu]\n", 
              blkLo, blkHi);
#endif

      for (b = 0; b < numBlk; b++) {
         /* Open temporary file as output for subblock */
         if (!(blkFile [b] = tmpfile())) {
            perror("OOC_Sort: failed opening temporary block file");
            return -1;
         }
         
         /* Subblock interval [blkLo2, blkHi2] of size blkSize2 */
         blkHi2 = blkLo + (b + 1) * blkLen2 - 1;
         blkSize2 = (blkHi2 - blkLo2 + 1) * recSize;

         if (blkSize2 <= maxBlkSize) {
            /* !!! Will be in-core sorted on recursion -> fork kid process,
             * !!! but don't fork more than numProc kids; wait if necessary */
            while (procCnt >= numProc && wait(&stat) >= 0) {
               if (!WIFEXITED(stat) || WEXITSTATUS(stat))
                  return -1;

               procCnt--;
            }

            /* Now fork kid process */
            if (!(pid = fork())) {
               /* Recurse on subblocks with new input filehandle; */
               if (OOC_SortRecurse(in, blkLo2, blkHi2, blkFile [b], numBlk,
                                   maxBlkSize, numProc, recSize, sortBuf,
                                   pqueue, keyData))
                  exit(-1);
                                                                         
               /* !!! Apparently the parent's tmpfile isn't deleted when the
                * !!! child exits (which is what we want), though some
                * !!! sources suggest using _Exit() instead; is this
                * !!! implementation specific?  */
               exit(0);
            }
            else if (pid < 0) {
               perror("OOC_Sort: failed to fork subprocess");
               return -1;
            }

#ifdef DEBUG_OOC_FORK
            fprintf(stderr, "OOC_Sort: forking kid %d for block %d\n",
                    procCnt, b);
#endif            
      
            /* Parent continues here */
            procCnt++;
         }
         else {
            /* Recurse on subblock; without forking */
            if (OOC_SortRecurse(in, blkLo2, blkHi2, blkFile [b], numBlk,
                                maxBlkSize, numProc, recSize, sortBuf,
                                pqueue, keyData))
               return -1;
         }

         /* Prepare for next block */ 
         blkLo2 = blkHi2 + 1;
      }

      /* !!! Wait for any forked processes to terminate */
      while (procCnt && wait(&stat) >= 0) {
         if (!WIFEXITED(stat) || WEXITSTATUS(stat))
            return -1;

         procCnt--;
      }

      /* ===============================================================
       * Subblocks are now sorted; prepare priority queue by peeking and
       * enqueueing first record from corresponding temporary file
       * =============================================================== */

#ifdef DEBUG_OOC_SORT
      fprintf(stderr, "OOC_Sort: merging block [%lu - %lu]\n", blkLo, blkHi);
#endif

      for (b = 0; b < numBlk; b++) {
#ifdef DEBUG_OOC_SORT
         fseek(blkFile [b], 0, SEEK_END);
         if (ftell(blkFile [b]) % recSize) {
            fprintf(stderr, "OOC_Sort: truncated record in tmp block "
                    "file %d\n", b);
            return -1;
         }
          
         fprintf(stderr, "OOC_Sort: tmp block file %d contains %ld rec\n",
                 b, ftell(blkFile [b]) / recSize);
#endif

         rewind(blkFile [b]);

         if (OOC_SortPeek(blkFile [b], recSize, rec)) {
            fprintf(stderr, "OOC_Sort: error reading from block file\n");
            return -1;
         }      
         
         /* Enqueue record along with its Morton code as priority */
         pri = OOC_Key2Morton(keyData -> key(rec), keyData -> bbOrg, 
                              keyData -> mortonScale);
                              
         if (OOC_SortPush(pqueue, pri, b) < 0) {
            fprintf(stderr, "OOC_Sort: failed priority queue insertion\n");
            return -1;
         }
      }
      
      /* ==========================================================
       * Subblocks now sorted and priority queue filled; merge from
       * temporary files 
       * ========================================================== */

      do {
         /* Get next node in priority queue, read next record in corresponding
          * block, and send to output */
         b = OOC_SortPop(pqueue);
               
         if (b >= 0) {
            /* Priority queue non-empty */
            if (OOC_SortRead(blkFile [b], recSize, rec)) {
               /* Corresponding record should still be in the input block */
               fprintf(stderr, "OOC_Sort: unexpected end reading block file\n");
               return -1;
            }

            if (OOC_SortWrite(out, recSize, rec)) {
               fprintf(stderr, "OOC_Sort; error writing output file\n");
               return -1;
            }

#ifdef DEBUG_OOC_SORT
            pqCnt++;
#endif            
            
            /* Peek next record from same block and insert in priority queue */
            if (!OOC_SortPeek(blkFile [b], recSize, rec)) {
               /* Block not exhausted */
               pri = OOC_Key2Morton(keyData -> key(rec), keyData -> bbOrg, 
                                    keyData -> mortonScale);
               
               if (OOC_SortPush(pqueue, pri, b) < 0) {
                  fprintf(stderr, "OOC_Sort: failed priority queue insert\n");
                  return -1;
               }
            }
         }
      } while (b >= 0);
      
#ifdef DEBUG_OOC_SORT
      fprintf(stderr, "OOC_Sort: dequeued %lu rec\n", pqCnt);
      fprintf(stderr, "OOC_Sort: merged file contains %lu rec\n",
              ftell(out) / recSize);
#endif
      
      /* Priority queue now empty -> done; close temporary subblock files,
       * causing them to be automatically deleted. */
      for (b = 0; b < numBlk; b++) {
         fclose(blkFile [b]);
      }
         
      /* !!! Commit output file to disk before caller reads it; omitting
       * !!! this apparently leads to corrupted files (race condition?) when
       * !!! the caller tries to read them! */
      fflush(out);
      fsync(fileno(out));
   }

   else {   
      /* ======================================
       * Block is small enough for in-core sort
       * ====================================== */   
      int   ifd = fileno(in), ofd = fileno(out);

#ifdef DEBUG_OOC_SORT
      fprintf(stderr, "OOC_Sort: Proc %d (%d/%d) sorting block [%lu - %lu]\n", 
              getpid(), procCnt, numProc - 1, blkLo, blkHi);

#endif

      /* Atomically seek and read block into in-core sort buffer */
      /* !!! Unbuffered I/O via pread() avoids potential race conditions
       * !!! and buffer corruption which can occur with lseek()/fseek()
       * !!! followed by read()/fread(). */
      if (pread(ifd, sortBuf, blkSize, blkLo * recSize) != blkSize) {
         perror("OOC_Sort: error reading from input file");
         return -1;
      }            

      /* Quicksort block in-core and write to output file */
      qsort(sortBuf, blkLen, recSize, OOC_KeyCompare);

      if (write(ofd, sortBuf, blkSize) != blkSize) {
         perror("OOC_Sort: error writing to block file");
         return -1;
      }
      
      fsync(ofd);

#ifdef DEBUG_OOC_SORT
      fprintf(stderr, "OOC_Sort: proc %d wrote %ld records\n",
              getpid(), lseek(ofd, 0, SEEK_END) / recSize);
#endif      
   }      

   return 0;
}  

 

int OOC_Sort (FILE *in, FILE *out, unsigned numBlk, 
              unsigned long blkSize, unsigned numProc, unsigned recSize, 
              FVECT bbOrg, RREAL bbSize, RREAL *(*key)(const void*))
/* Sort records in inFile and append to outFile by subdividing inFile into
 * small blocks, sorting these in-core, and merging them out-of-core via a
 * priority queue.
 *
 * This is implemented as a recursive (numBlk)-way sort; the input is
 * successively split into numBlk smaller blocks until these are of size <=
 * blkSize, i.e. small enough for in-core sorting, then merging the sorted
 * blocks as the stack unwinds. The in-core sort is parallelised over
 * numProx processes. 
 *
 * Parameters are as follows:
 * inFile      Opened input file containing unsorted records
 * outFile     Opened output file containing sorted records
 * numBlk      Number of blocks to divide into / merge from 
 * blkSize     Max block size and size of in-core sort buffer, in bytes 
 * numProc     Number of parallel processes for in-core sort 
 * recSize     Size of input records in bytes
 * bbOrg       Origin of bounding box containing record keys for Morton code
 * bbSize      Extent of bounding box containing record keys for Morton code
 * key         Callback to access 3D coords from records for Morton code
 */
{
   unsigned long     numRec;
   OOC_SortQueue     pqueue;
   char              *sortBuf = NULL;
   int               stat;
   
   if (numBlk < 1)
      numBlk = 1;

   /* Open input file & get size in number of records */
   if (fseek(in, 0, SEEK_END) < 0) {
      fputs("OOC_Sort: failed seek in input file\n", stderr);
      return -1;
   }
   
   if (!(numRec = ftell(in) / recSize)) {
      fputs("OOC_Sort: empty input file\n", stderr);
      return -1;
   }

   /* Allocate & init priority queue */
   if (!(pqueue.node = calloc(numBlk, sizeof(OOC_SortQueueNode)))) {
      fputs("OOC_Sort: failure allocating priority queue\n", stderr);
      return -1;
   }
   pqueue.tail = 0;
   pqueue.len = numBlk;

   /* Allocate in-core sort buffa */
   if (!(sortBuf = malloc(blkSize))) {
      fprintf(stderr, "OOC_Sort: failure allocating in-core sort buffer");
      return -1;
   }
   
   /* Set up key data to pass to qsort()'s comparison func */
   keyData.key = key;
   keyData.mortonScale = OOC_MORTON_MAX / bbSize;
   VCOPY(keyData.bbOrg, bbOrg);
     
   stat = OOC_SortRecurse(in, 0, numRec - 1, out, numBlk, blkSize, numProc,
                          recSize, sortBuf, &pqueue, &keyData);

   /* Cleanup */
   free(pqueue.node);
   free(sortBuf);
   
   return stat;           
}

#endif /* NIX / PMAP_OOC */
