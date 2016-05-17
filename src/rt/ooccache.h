/* 
   =======================================================================
   Cache for out-of-core octree.
   
   Uses an LRU page replacement scheme. Page table nodes are stored in
   order of reference in a doubly linked list, with the most recently used
   (MRU) page at the head, and the least recently used (LRU) at the tail.
   
   Nodes are referenced by the page's base address via a hashtable for
   constant time access (in the absence of hash collisions).  Hash
   collisions can be limited by adjusting the hashtable load factor
   OOC_CACHE_LOAD; this is the fraction of the number of pages that will
   actually be cached.
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   =======================================================================
   
   $Id: ooccache.h,v 2.1 2016/05/17 17:39:47 rschregle Exp $
*/


#ifndef OOC_CACHE_H
   #define OOC_CACHE_H
   
   #include <stdio.h>
   #include <stdint.h>


   
   /* Hashtable load factor to limit hash collisions */
   #define  OOC_CACHE_LOAD       0.75
   
   /* Minimum number of cache pages */
   #define  OOC_CACHE_MINPAGES   8

   /* Hashing function for key -> pagetable mapping (assumes n is prime) */
   #define  OOC_CACHE_PRIME      19603
   #define  OOC_CACHE_HASH(k,n)  ((k) * OOC_CACHE_PRIME % (n))
   
   /* Next index in pagetable probe sequence for collision resolution */
   #define  OOC_CACHE_NEXT(i,n)  (((i) + 1) % (n))

   /* Cache key and index types */
   typedef  uint32_t OOC_CacheKey;
   typedef  uint32_t OOC_CacheIdx;   
   
   /* Undefined value for cache index */
   #define  OOC_CACHEIDX_NULL    UINT32_MAX

   
      
   typedef struct {
      OOC_CacheKey   key;        /* Key to resolve hashing collisions */
      OOC_CacheIdx   prev, next; /* Previous/next page in history */
      void           *data;      /* Start of page data */
   } OOC_CacheNode;
      
   typedef struct {
      unsigned       recSize,    /* Data record size in bytes */
                     recPerPage; /* Page size in number of records */
      unsigned long  pageSize;   /* Page size in bytes (precomp) */
      OOC_CacheIdx   pageCnt,    /* Num pages used */
                     numPages;   /* Pagetable size */

      /* NOTE: The hashtable load factor limits the number of pages actually
       * used, such that pageCnt <= OOC_CACHE_LOAD * numPages */
                     
      OOC_CacheNode  *pageTable; /* Pagetable with numPages nodes */
      OOC_CacheIdx   mru, lru;   /* Most/least recently used page
                                    == head/tail of page history */
      OOC_CacheKey   lastKey;    /* Previous key to detect repeat lookups */
      OOC_CacheNode  *lastPage;  /* Previous page for repeat lookups */
      
      unsigned long  numReads,   /* Statistics counters */
                     numHits,
                     numColl,
                     numRept;
   } OOC_Cache;
   
   
   
   /* Initialise cache containing up to numPages pages; this is rounded up
    * the next prime number to reduce hash clustering and collisions.  The
    * page size recPerPage is specified in units of records of size recSize
    * bytes.  Returns 0 on success, else -1.  */
   int OOC_CacheInit (OOC_Cache *cache, unsigned numPages, 
                      unsigned recPerPage, unsigned recSize);
                      
   /* Return pointer to cached record with index recIdx, loading from file
    * and possibly evicting the LRU page */
   void *OOC_CacheData (OOC_Cache *cache, FILE *file, unsigned long recIdx);
   
   /* Delete cache and free allocated pages */
   void OOC_DeleteCache (OOC_Cache *cache);
#endif
