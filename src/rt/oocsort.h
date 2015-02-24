/* 
   ==================================================================
   Header for N-way hybrid out-of-core merge sort
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts   
   ==================================================================
   
   $Id $
*/

#ifndef OOC_SORT_H
#define OOC_SORT_H
   #include <stdint.h>
   
   #define OOC_SORT_BUFSIZE (1 << 13)     /* 8kb */
   
   typedef uint64_t OOC_Sort_Key;

#if 0   
   int OOC_Sort (const char *inFile, const char *outFile, 
                 unsigned long blkSize, unsigned recSize, 
                 OOC_Sort_Key (*priority)(const void *));                 
   /* Sort records in inFile and output to outFile by (a) internally
    * quicksorting block of blkSize bytes at a time, then (b) merging them
    * via a priority queue.  RecSize specifies the size in bytes of each
    * data record.  The priority() callback evaluates a record's priority
    * and must be supplied by the caller.  */
#else    
   int OOC_Sort (const char *inFile, const char *outFile, 
                 unsigned numBlocks, unsigned recSize, 
                 OOC_Sort_Key (*priority)(const void *));                 
   /* Sort records in inFile and output to outFile by (a) internally
    * quicksorting numBlocks blocks at a time, then (b) merging them
    * via a priority queue.  RecSize specifies the size in bytes of each
    * data record.  The priority() callback evaluates a record's priority
    * (ordinal index) and must be supplied by the caller.  */
#endif

#endif
