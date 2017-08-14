#ifndef lint
static const char RCSid[] = "$Id: ooccache.c,v 2.2 2017/08/14 21:12:10 rschregle Exp $";
#endif


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
   
   $Id: ooccache.c,v 2.2 2017/08/14 21:12:10 rschregle Exp $
*/


#if !defined(_WIN32) && !defined(_WIN64) || defined(PMAP_OOC)
/* No Windoze support for now */

#include <stdlib.h>
#include <unistd.h>
#include "ooccache.h"



static OOC_CacheIdx OOC_PrevPrime (OOC_CacheIdx n)
/* Return largest prime number <= n */
{
   OOC_CacheIdx n2 = n + 1, i;
   
   do {
      n2--;
      for (i = 2; n2 % i && i <= n2 >> 1; i++);
   } while (!(n2 % i));
   
   return n2;
}



int OOC_CacheInit (OOC_Cache *cache, unsigned numPages, 
                   unsigned recPerPage, unsigned recSize)
{
   OOC_CacheIdx i;

   if (!cache)
      return -1;

   /* Clamp number of pages if necessary */
   if (numPages < OOC_CACHE_MINPAGES) {
      numPages = OOC_CACHE_MINPAGES;
      fputs("OOC_CacheInit: clamping num pages to min\n", stderr);
   }
      
   if (numPages >= OOC_CACHEIDX_NULL) {
      numPages = OOC_CACHEIDX_NULL - 1;
      fputs("OOC_CacheInit: clamping num pages to max\n", stderr); 
   }

   /* Round numPages to next prime to reduce hash clustering & collisions */      
   numPages = OOC_PrevPrime(numPages);
      
   cache -> recPerPage = recPerPage;
   cache -> numPages = numPages;
   cache -> recSize = recSize;
   cache -> pageSize = (unsigned long)cache -> recPerPage * 
                       cache -> recSize;
   cache -> pageCnt = 0;
   cache -> numHits = cache -> numReads = cache -> numColl = 
      cache -> numRept = 0;
   cache -> mru = cache -> lru = OOC_CACHEIDX_NULL;
   cache -> lastKey = 0;
   cache -> lastPage = NULL;
   
   /* Allocate & init hashtable */
   if (!(cache -> pageTable = calloc(numPages, sizeof(OOC_CacheNode)))) {
      perror("OOC_InitCache: failed pagetable allocation");
      return -1;
   }
   
   for (i = 0; i < numPages; i++) {
      cache -> pageTable [i].data = NULL;
      cache -> pageTable [i].key = 0;
      cache -> pageTable [i].prev = cache -> pageTable [i].next = 
         OOC_CACHEIDX_NULL;
   }
   
#ifdef OOC_CACHE_STATS
   fprintf(stderr, 
           "OOC_InitCache: caching %d pages @ %d records (%.2f Mb)\n",
           cache -> numPages, cache -> recPerPage, 
           (float)cache -> numPages * cache -> pageSize / (1L << 20));
#endif 
  
   return 0;
}       



static OOC_CacheIdx OOC_CacheFind (OOC_Cache *cache, OOC_CacheKey key)
/* Return index to pagetable node for given key; collisions are resolved
 * via the probe sequence defined by OOC_CACHE_NEXT() */ 
{
   OOC_CacheIdx   idx = OOC_CACHE_HASH(key, cache -> numPages);
   OOC_CacheNode  *page = cache -> pageTable + idx;

   /* Find matching pagetable node starting from idx and stop at the first
    * empty node; a full table would result in an infinite loop! */
   while (page -> data && page -> key != key) {
      idx = OOC_CACHE_NEXT(idx, cache -> numPages);
      page = cache -> pageTable + idx;
      cache -> numColl++;
   }

   return idx;
}



static void OOC_CacheUnlink (OOC_Cache *cache, OOC_CacheIdx pageIdx)
/* Unlink page from history chain, updating LRU in cache if necessary */
{
   OOC_CacheNode  *page = cache -> pageTable + pageIdx;
   
   /* Unlink prev/next index */
   if (page -> prev != OOC_CACHEIDX_NULL)
      cache -> pageTable [page -> prev].next = page -> next;
      
   if (page -> next != OOC_CACHEIDX_NULL)
      cache -> pageTable [page -> next].prev = page -> prev;
     
   /* Update LRU index */ 
   if (pageIdx == cache -> lru)
      cache -> lru = page -> prev;
}



static void OOC_CacheRelink (OOC_Cache *cache, OOC_CacheIdx pageIdx)
/* Update history chain with new page index */
{
   OOC_CacheNode  *page = cache -> pageTable + pageIdx;
   
   if (page -> prev != OOC_CACHEIDX_NULL)
      cache -> pageTable [page -> prev].next = pageIdx;
      
   if (page -> next != OOC_CACHEIDX_NULL)
      cache -> pageTable [page -> next].prev = pageIdx;
}      



static void OOC_CacheSetMRU (OOC_Cache *cache, OOC_CacheIdx pageIdx)
/* Set page as most recently used and move to head of history chain */
{
   OOC_CacheNode  *page = cache -> pageTable + pageIdx;
   
   OOC_CacheUnlink(cache, pageIdx);
   page -> prev = OOC_CACHEIDX_NULL;
   page -> next = cache -> mru;
   cache -> pageTable [cache -> mru].prev = pageIdx;
   cache -> mru = pageIdx;
}



static void *OOC_CacheDelete (OOC_Cache *cache, OOC_CacheIdx pageIdx)
/* Delete cache node at pageIdx and tidy up fragmentation in pagetable to
 * minimise collisions.  Returns pointer to deleted node's page data
 *
 * !!! Note cache -> pageCnt is NOT decremented as this is immediately
 * !!! followed by a call to OOC_CacheNew() in OOC_CacheData() */
{
   OOC_CacheNode  *page = cache -> pageTable + pageIdx, *nextPage;
   OOC_CacheIdx   nextIdx;
   void           *delData = page -> data, *nextData;
   
   /* Unlink page from history chain */
   OOC_CacheUnlink(cache, pageIdx);
   
   /* Free pagetable node */
   page -> data = NULL;
   
   nextIdx = OOC_CACHE_NEXT(pageIdx, cache -> numPages);
   nextPage = cache -> pageTable + nextIdx;
   
   /* Close gaps in probe sequence until empty node found */
   while ((nextData = nextPage -> data)) {
      /* Mark next node as free and reprobe for its key */
      nextPage -> data = NULL;
      pageIdx = OOC_CacheFind(cache, nextPage -> key);
      page = cache -> pageTable + pageIdx;
      
      /* Set page data pointer & reactivate page */
      page -> data = nextData;
      
      if (pageIdx != nextIdx) {
         /* Next page maps to freed node; relocate to close gap and update
          * history chain and LRU/MRU if necessary */
         page -> key = nextPage -> key;         /* memcpy() faster ? */
         page -> prev = nextPage -> prev;
         page -> next = nextPage -> next;
         OOC_CacheRelink(cache, pageIdx);
         
         if (cache -> lru == nextIdx)
            cache -> lru = pageIdx;
            
         if (cache -> mru == nextIdx)
            cache -> mru = pageIdx;
      }
            
      nextIdx = OOC_CACHE_NEXT(nextIdx, cache -> numPages);
      nextPage = cache -> pageTable + nextIdx;      
   }
   
   return delData;
}



static OOC_CacheNode *OOC_CacheNew (OOC_Cache *cache, OOC_CacheIdx pageIdx, 
                                    OOC_CacheKey pageKey, void *pageData)
/* Init new pagetable node with given index, key and data. Returns pointer
 * to new node */
{
   OOC_CacheNode  *page = cache -> pageTable + pageIdx;
   
   page -> prev = page -> next = OOC_CACHEIDX_NULL;
   page -> key = pageKey;
   page -> data = pageData;
   
   return page;
}
 
   

#ifdef OOC_CACHE_CHECK
static int OOC_CacheCheck (OOC_Cache *cache)
/* Run sanity checks on hashtable and page list */
{
   OOC_CacheIdx   i, pageCnt = 0;
   OOC_CacheNode  *page;
   
   /* Check pagetable mapping */
   for (i = 0, page = cache->pageTable; i < cache->numPages; i++, page++)
      if (page -> data) {
         /* Check hash func maps page key to this node */
         if (OOC_CacheFind(cache, page -> key) != i) {
            fputs("OOC_CacheCheck: hashing inconsistency\n", stderr);
            return -1;
         }
         
         pageCnt++;
      }

   if (pageCnt != cache -> pageCnt) {
      fputs("OOC_CacheCheck: pagetable count inconsistency\n", stderr);
      return -1;
   }
      
   if (cache -> mru == OOC_CACHEIDX_NULL || 
       cache -> lru == OOC_CACHEIDX_NULL) {
      fputs("OOC_CacheCheck: undefined MRU/LRU\n", stderr);
      return -1;
   }

   pageCnt = 0;
   i = cache -> mru;
   
   /* Check page history */
   while (i != OOC_CACHEIDX_NULL) {
      /* Check link to previous page (should be none for MRU) */
      page = cache -> pageTable + i;
      
      if (i == cache -> mru)
         if (page -> prev != OOC_CACHEIDX_NULL) {
            fputs("OOC_CacheCheck: MRU inconsistency\n", stderr);
            return -1;
         }
         else;
      else if (page -> prev == OOC_CACHEIDX_NULL || 
               cache -> pageTable [page -> prev].next != i) {
         fputs("OOC_CacheCheck: prev page inconsistency\n", stderr);
         return -1;
      }
      
      /* Check link to next page node (should be none for LRU) */
      if (i == cache -> lru)
         if (page -> next != OOC_CACHEIDX_NULL) {
            fputs("OOC_CacheCheck: LRU inconsistency\n", stderr);
            return -1;
         }
         else;
      else if (page -> next != OOC_CACHEIDX_NULL &&
               cache -> pageTable [page -> next].prev != i) {
         fputs("OOC_CacheCheck: next page inconsistency\n", stderr);
         return -1;
      }
      
      /* Check LRU is at tail of page history */
      if (page -> next == OOC_CACHEIDX_NULL && i != cache -> lru) {
         fputs("OOC_CacheCheck: page history tail not LRU\n", stderr);
         return -1;
      }
      
      i = page -> next;
      pageCnt++;
   }
   
   if (pageCnt != cache -> pageCnt) {
      fputs("OOC_CacheCheck: page history count inconsistency\n", stderr);
      return -1;
   }
   
   return 0;
}
#endif



void *OOC_CacheData (OOC_Cache *cache, FILE *file, unsigned long recIdx)
{
   const  OOC_CacheKey  pageKey = recIdx / cache -> recPerPage;
   OOC_CacheNode        *page = NULL;

   /* We assume locality of reference, and that it's likely the same page
    * will therefore be accessed sequentially; in this case we just reuse
    * the pagetable index */
   if (pageKey != cache -> lastKey || !cache -> pageCnt) {
      OOC_CacheIdx   pageIdx = OOC_CacheFind(cache, pageKey);
      void           *pageData;
      
      page = cache -> pageTable + pageIdx;
      
      if (!page -> data) {
         /* Page not cached; insert in pagetable */
         const unsigned long  filePos = cache -> pageSize * pageKey;
            
         if (cache -> pageCnt < OOC_CACHE_LOAD * cache -> numPages) {
            /* Cache not fully loaded; allocate new page data */
            if (!(pageData = calloc(cache -> recPerPage, cache -> recSize))) {
               perror("OOC_CacheData: failed page allocation");
               return NULL;
            }
            
            /* Init LRU and MRU if necessary */
            if (cache -> mru == OOC_CACHEIDX_NULL ||  
                cache -> lru == OOC_CACHEIDX_NULL)
               cache -> mru = cache -> lru = pageIdx;
               
            cache -> pageCnt++;
         }
         else {
            /* Cache fully loaded; adding more entries would aggravate
             * hash collisions, so evict LRU page. We reuse its
             * allocated page data and then overwrite it. */
            pageData = OOC_CacheDelete(cache, cache -> lru);
            
            /* Update page index as it may have moved forward in probe
             * sequence after deletion */
            pageIdx = OOC_CacheFind(cache, pageKey);
         }
         
         /* Init new pagetable node, will be linked in history below */
         page = OOC_CacheNew(cache, pageIdx, pageKey, pageData);
         
         /* Load page data from file; the last page may be truncated */
         if (!pread(fileno(file), page -> data, cache -> pageSize, filePos)) {
            fputs("OOC_CacheData: failed seek/read from data stream", stderr);
            return NULL;
         }
      }
      else 
         /* Page in cache */
         cache -> numHits++;

      if (cache -> mru != pageIdx)
         /* Cached page is now most recently used; unlink and move to front
          * of history chain and update MRU */
         OOC_CacheSetMRU(cache, pageIdx);
      
      /* Update last page to detect repeat accesses */
      cache -> lastKey = pageKey;
      cache -> lastPage = page;
   }
   else {
      /* Repeated page, skip hashing and just reuse last page */
      page = cache -> lastPage;
      cache -> numHits++;
      cache -> numRept++;
   }
   
   cache -> numReads++;

#ifdef OOC_CACHE_CHECK
   if (OOC_CacheCheck(cache))
      return NULL;
#endif
           
   /* Return pointer to record within cached page */
   return page -> data + (recIdx % cache -> recPerPage) * cache -> recSize;
}

   
                      
void OOC_DeleteCache (OOC_Cache *cache)
{
   OOC_CacheIdx   i;

   for (i = 0; i < cache -> numPages; i++)
      if (cache -> pageTable [i].data)
         free(cache -> pageTable [i].data);
         
   free(cache -> pageTable);
   cache -> pageTable = NULL;
   cache -> pageCnt = 0;
   cache -> mru = cache -> lru = OOC_CACHEIDX_NULL;
}

#endif /* NIX / PMAP_OOC */
