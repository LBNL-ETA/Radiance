#ifndef lint
static const char RCSid[] = "$Id: oococt.c,v 2.4 2017/08/14 21:12:10 rschregle Exp $";
#endif


/* 
   ======================================================================
   Out-of-core octree data structure
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: oococt.c,v 2.4 2017/08/14 21:12:10 rschregle Exp $
*/


#if !defined(_WIN32) && !defined(_WIN64) || defined(PMAP_OOC)
/* No Windoze support for now */

#include "oococt.h"
#include "rtio.h"
#include <stdlib.h>
#include <unistd.h>



void OOC_Null (OOC_Octree *oct)
/* Init empty octree */
{
   oct -> maxNodes = oct -> numNodes = oct -> leafMax = oct -> maxDepth =
      oct -> numData = oct -> size = oct -> recSize = oct -> mortonScale = 0;
   oct -> org [0] = oct -> org [1] = oct -> org [2] = 
      oct -> bound [0] = oct -> bound [1] = oct -> bound [2] = 0;
   oct -> key = NULL;
   oct -> nodes = NULL;
   oct -> leafFile = NULL;
   oct -> cache = NULL;
}



void OOC_Init (OOC_Octree *oct, unsigned recSize, FVECT org, RREAL size, 
               RREAL *(*key)(const void*), FILE *leafFile)
{
   oct -> maxNodes = oct -> numNodes = oct -> leafMax = 
      oct -> maxDepth = oct -> numData = 0;
   VCOPY(oct -> org, org);
   oct -> size = oct -> bound[0] = oct -> bound[1] = oct -> bound[2] = size;
   VADD(oct -> bound, oct -> bound, org);
   oct -> mortonScale = size > 0 ? OOC_MORTON_MAX / size : 0;
   oct -> recSize = recSize;
   oct -> key = key;
   oct -> nodes = NULL;
   oct -> leafFile = leafFile;
   oct -> cache = NULL;    /* Cache currently initialised externally */
}



OOC_Node *NewNode (OOC_Octree *oct)
/* Allocate new octree node, enlarge array if necessary.
   Return pointer to new node or NULL if failed. */
{
   OOC_Node *n = NULL;

   if (oct -> numNodes >= OOC_NODEIDX_MAX) {
      /* Node index overflow */
      fprintf(stderr, "OOC_NewNode: node index overflow (numNodes = %d)\n",
              oct -> numNodes);
      return NULL;
   }
      
   if (++oct -> numNodes > oct -> maxNodes) {
      /* Need to allocate root or enlarge array */
      oct -> maxNodes += OOC_BLKSIZE / sizeof(OOC_Node);
      oct -> nodes = realloc(oct -> nodes, 
                             oct -> maxNodes * sizeof(OOC_Node));
      if (!oct -> nodes) {
         perror("OOC_NewNode: couldn't realloc() nodes");
         return NULL;
      }
   }
   
   n = oct -> nodes + oct -> numNodes - 1;
   n -> node.num = n -> node.kid = n -> node.branch = n -> node.type = 0;
   
   return n;
}



int OOC_GetData (OOC_Octree *oct, OOC_DataIdx idx, void *rec)
/* Reads indexed data record from leaf file and copies it to rec, else
 * returns -1 on failure */
{
   if (!oct || !rec)
      return -1;
   
   if (idx >= oct -> numData) {
      fprintf(stderr, "OOC_GetData: invalid data record index %u\n", idx);
      return -1;
   }   

   if (oct -> cache) {
      /* Retrieve record from leaf file via I/O cache */
      void  *cachedRec;

      if (!(cachedRec = OOC_CacheData(oct -> cache, oct -> leafFile, idx)))
         return -1;
            
      /* It's safer to copy the record to the caller's local buffer as a
       * returned pointer may be invalidated by subsequent calls if the
       * page it points to is swapped out */
      memcpy(rec, cachedRec, oct -> recSize);
   }
   else {
      /* No I/O caching; do (SLOW!) per-record I/O from leaf file */ 
      const unsigned long  pos = (unsigned long)oct -> recSize * idx;

      if (pread(fileno(oct -> leafFile), rec, oct -> recSize, pos) !=
          oct -> recSize) {
         perror("OOC_GetData: failed seek/read in leaf file");
         return -1;
      }
   }

   return 0;
}



int OOC_InBBox (const OOC_Octree *oct, const FVECT org, RREAL size, 
                const FVECT key)
/* Check whether key lies within bounding box defined by org and size.
 * Compares Morton codes rather than the coordinates directly to account
 * for dicretisation error introduced by the former. */
{
   const OOC_MortonIdx  keyZ = OOC_KEY2MORTON(key, oct);
   FVECT                bbox = {size, size, size};
   
   VADD(bbox, org, bbox);
   return keyZ > OOC_KEY2MORTON(org, oct) && 
          keyZ < OOC_KEY2MORTON(bbox, oct);
}



int OOC_Branch (const OOC_Octree *oct, const FVECT org, RREAL size, 
                const FVECT key, FVECT nuOrg)
/* Return index of octant containing key and corresponding origin if nuOrg
 * != NULL, or -1 if key lies outside all octants.  Size is already assumed
 * to be halved, i.e.  corresponding to the length of an octant axis.
 * Compares Morton codes rather than the coordinates directly to account
 * for dicretisation error introduced by the former. */
{
   int   o;
   FVECT octOrg;
   
   for (o = 0; o < 8; o++) {
      VCOPY(octOrg, org);
      OOC_OCTORIGIN(octOrg, o, size);
      
      if (OOC_InBBox(oct, octOrg, size, key)) {
         if (nuOrg)
            VCOPY(nuOrg, octOrg);
            
         return o;
      }
   }
   
   return -1;   
}



int OOC_BranchZ (const OOC_Octree *oct, const FVECT org, RREAL size, 
                 OOC_MortonIdx keyZ, FVECT nuOrg)
/* Optimised version of OOC_Branch() with precomputed key Morton code */
{
   int         o;
   const FVECT cent = {size, size, size};
   FVECT       octOrg, bbox;
   
   for (o = 0; o < 8; o++) {
      VCOPY(octOrg, org);
      OOC_OCTORIGIN(octOrg, o, size);
      VADD(bbox, octOrg, cent);
      
      if (keyZ > OOC_KEY2MORTON(octOrg, oct) && 
          keyZ < OOC_KEY2MORTON(bbox, oct)) {
         if (nuOrg)
            VCOPY(nuOrg, octOrg);
            
         return o;
      }
   }
   
   return -1;   
}



OOC_DataIdx OOC_GetKid (const OOC_Octree *oct, OOC_Node **node, unsigned kid)
/* For a given octree node, return the sum of the data counters for kids
 * [0..k-1].  On return, the node either points to the k'th kid, or
 * NULL if the kid is nonexistent or the node is a leaf */                        
{
   unsigned       k;
   OOC_Node       *kidNode = OOC_KID1(oct, *node); /* NULL if leaf */
   OOC_DataIdx    dataIdx = 0;
   
   for (k = 0; k < kid; k++) {
      if (OOC_ISLEAF(*node))
         /* Sum data counters of leaf octants */
         dataIdx += (*node) -> leaf.num [k];
      else 
         /* Sum data counter of inner node's nonempty kids and advance
          * pointer to sibling */
         if (OOC_ISBRANCH(*node, k)) {
            dataIdx += OOC_DATACNT(kidNode);
            kidNode++;
         }
   }
   
   /* Update node pointer to selected kid (or NULL if nonexistent or node is
    * a leaf) */
   *node = OOC_ISBRANCH(*node, kid) ? kidNode : NULL;
   return dataIdx;
}



#ifdef DEBUG_OOC
int OOC_Check (OOC_Octree *oct, const OOC_Node *node, 
               FVECT org, RREAL size, OOC_DataIdx dataIdx)
/* Traverse tree & check for valid nodes; oct -> leafFile must be open;
 * return 0 if ok, otherwise -1 */
{
   unsigned    k;
   FVECT       kidOrg;
   const RREAL kidSize = size * 0.5;

   if (!oct || !node)
      return -1;
      
   if (OOC_ISLEAF(node)) {
      /* Node is a leaf; check records in each octant */
      char        rec [oct -> recSize];   /* Violates C89! */
      OOC_OctCnt  d;
      RREAL       *key;
      
      for (k = 0; k < 8; k++) {
         VCOPY(kidOrg, org);
         OOC_OCTORIGIN(kidOrg, k, kidSize);
         
         for (d = node -> leaf.num [k]; d; d--, dataIdx++) {
#ifdef DEBUG_OOC_CHECK
            fprintf(stderr, "OOC_Check: checking record %lu\n", 
                    (unsigned long)dataIdx);
#endif            
            if (OOC_GetData(oct, dataIdx, rec)) {
               fprintf(stderr, "OOC_Check: failed getting record at %lu\n",
                       (unsigned long)dataIdx);
               return -1;
            }
            
            key = oct -> key(rec);
            if (!OOC_InBBox(oct, kidOrg, kidSize, key)) {
               fprintf(stderr, "OOC_Check: key [%f, %f, %f] (zIdx %020lu) "
                       "outside bbox [[%f-%f], [%f-%f], [%f-%f]] "
                       "in octant %d (should be %d)\n",
                       key [0], key [1], key [2], 
                       OOC_KEY2MORTON(key, oct),
                       kidOrg [0], kidOrg [0] + kidSize, 
                       kidOrg [1], kidOrg [1] + kidSize, 
                       kidOrg [2], kidOrg [2] + kidSize,
                       k, OOC_Branch(oct, org, kidSize, key, NULL));
               /* return -1; */               
            }
         }
      }
   }
   else {
      /* Internal node; recurse on all kids */
      const OOC_Node *kid = OOC_KID1(oct, node);
      OOC_DataIdx    numData = 0;
      
      if (node -> node.kid >= oct -> numNodes) {
         fprintf(stderr, "OOC_Check: invalid node index %u\n", 
                 node -> node.kid);
         return -1;
      }
         
      if (!node -> node.num) {
         fputs("OOC_Check: empty octree node\n", stderr);
         return -1;
      }
      
      for (k = 0; k < 8; k++)
         if (OOC_ISBRANCH(node, k)) {
            VCOPY(kidOrg, org);
            OOC_OCTORIGIN(kidOrg, k, kidSize);
            
            if (OOC_Check(oct, kid, kidOrg, kidSize, dataIdx))
               return -1;
            
            dataIdx += OOC_DATACNT(kid);
            numData += OOC_DATACNT(kid);
            kid++;
         }
      
      if (node -> node.num != numData) {
         fprintf(stderr, 
                 "Parent/kid data count mismatch; expected %u, found %u\n",
                 node -> node.num, numData);
         return -1;
      }
   }
   
   return 0;
}
#endif
         


int OOC_SaveOctree (const OOC_Octree *oct, FILE *out)
/* Appends octree nodes to specified file along with metadata. Uses
 * RADIANCE's portable I/O routines.  Returns 0 on success, else -1.  Note
 * the internal representation of the nodes is platform specific as they
 * contain unions, therefore we use the portable I/O routines from the
 * RADIANCE library */
{
   int         i;
   OOC_NodeIdx nc;
   OOC_Node    *np = NULL;
   
   if (!oct || !oct->nodes || !oct->numData || !oct->numNodes || !out) {
      fputs("OOC_SaveOctree: empty octree\n", stderr);
      return -1;
   }
      
   /* Write octree origin, size, data record size, number of records, and
    * number of nodes */
   for (i = 0; i < 3; i++)
      putflt(oct -> org [i], out);
      
   putflt(oct -> size, out);
   putint(oct -> recSize,     sizeof(oct -> recSize),       out);
   putint(oct -> numData,     sizeof(oct -> numData),       out);
   putint(oct -> numNodes,    sizeof(oct -> numNodes),      out);
   
   /* Write nodes to file */
   for (nc = oct -> numNodes, np = oct -> nodes; nc; nc--, np++) {
      if (OOC_ISLEAF(np)) {
         /* Write leaf marker */
         putint(-1, 1, out);
         
         /* Write leaf octant counters */
         for (i = 0; i < 8; i++)
            putint(np -> leaf.num [i], sizeof(np -> leaf.num [i]), out);
      }
      else {
         /* Write node marker */
         putint(0, 1, out);
         
         /* Write node, rounding field size up to nearest number of bytes */
         putint(np -> node.kid,     (OOC_NODEIDX_BITS + 7) >> 3, out);
         putint(np -> node.num,     (OOC_DATAIDX_BITS + 7) >> 3, out);
         putint(np -> node.branch,  1,                           out);
      }
         
      if (ferror(out)) {
         fputs("OOC_SaveOctree: failed writing to node file", stderr);
         return -1;
      }
   }
   
   fflush(out);
   return 0;
}



int OOC_LoadOctree (OOC_Octree *oct, FILE *nodeFile, 
                    RREAL *(*key)(const void*), FILE *leafFile)
/* Load octree nodes and metadata from nodeFile and associate with records
 * in leafFile.  Uses RADIANCE's portable I/O routines.  Returns 0 on
 * success, else -1.  */
{
   int         i;
   OOC_NodeIdx nc;
   OOC_Node    *np = NULL;
   
   if (!oct || !nodeFile)
      return -1;

   OOC_Null(oct);
   
   /* Read octree origin, size, record size, #records and #nodes */
   for (i = 0; i < 3; i++)
      oct -> org [i] = getflt(nodeFile);

   oct -> size = getflt(nodeFile);
   oct -> bound [0] = oct -> bound [1] = oct -> bound [2] = oct -> size;
   VADD(oct -> bound, oct -> bound, oct -> org);   
   oct -> mortonScale   = OOC_MORTON_MAX / oct -> size;
   oct -> recSize       = getint(sizeof(oct -> recSize),       nodeFile);
   oct -> numData       = getint(sizeof(oct -> numData),       nodeFile);
   oct -> numNodes      = getint(sizeof(oct -> numNodes),      nodeFile);
   oct -> key           = key;
   oct -> leafFile      = leafFile;
   
   if (feof(nodeFile) || !oct -> numData || !oct -> numNodes) {
      fputs("OOC_LoadOctree: empty octree\n", stderr);
      return -1;
   }
   
   /* Allocate octree nodes */
   if (!(oct -> nodes = calloc(oct -> numNodes, sizeof(OOC_Node)))) {
      fputs("OOC_LoadOctree: failed octree allocation\n", stderr);
      return -1;
   }

   /* Read nodes from file */   
   for (nc = oct -> numNodes, np = oct -> nodes; nc; nc--, np++) {
      OOC_CLEARNODE(np);
      
      /* Read node type marker */
      if (getint(1, nodeFile)) {
         /* Read leaf octant counters and set node type */
         for (i = 0; i < 8; i++)
            np -> leaf.num [i] = getint(sizeof(np -> leaf.num [i]), 
                                        nodeFile);
         
         OOC_SETLEAF(np);
      }
      else {
         /* Read node, rounding field size up to nearest number of bytes */
         np -> node.kid    = getint((OOC_NODEIDX_BITS + 7) >> 3,  nodeFile);
         np -> node.num    = getint((OOC_DATAIDX_BITS + 7) >> 3,  nodeFile);
         np -> node.branch = getint(1,                            nodeFile);
      }
         
      if (ferror(nodeFile)) {
         fputs("OOC_LoadOctree: failed reading from node file\n", stderr);
         return -1;
      }
   }
   
   return 0;
}   



void OOC_Delete (OOC_Octree *oct)
/* Self-destruct */
{
   free(oct -> nodes);
   fclose(oct -> leafFile);
   
   if (oct -> cache)
      OOC_DeleteCache(oct -> cache);
}

#endif /* NIX / PMAP_OOC */
