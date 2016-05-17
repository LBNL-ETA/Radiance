/* 
   ======================================================================
   Out-of-core octree data structure
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: oococt.h,v 2.3 2016/05/17 17:39:47 rschregle Exp $
*/


#ifndef OOC_OCT_H
#define OOC_OCT_H

   #include "oocmorton.h"
   #include "ooccache.h"
   #include <stdint.h>
   #include <stdio.h>



   /* ===================================================================
    * CONSTS
    * =================================================================== */
    
   /* Octree index & data counter types */
   typedef uint32_t           OOC_NodeIdx;
   typedef uint32_t           OOC_DataIdx;  
   typedef uint8_t            OOC_OctCnt;
   
   /* Octree node field sizes for above */   
   #define OOC_NODEIDX_BITS   31
   #define OOC_DATAIDX_BITS   32
   
   /* Limits for above */
   #define OOC_NODEIDX_MAX    ((1L << OOC_NODEIDX_BITS) - 1)
   #define OOC_DATAIDX_MAX    ((1L << OOC_DATAIDX_BITS) - 2)
   #define OOC_OCTCNT_MAX     UINT8_MAX
   
   /* Reserved values for above */
   #define OOC_DATAIDX_ERR    (OOC_DATAIDX_MAX + 1)
   
   /* Block size for node allocation (8k) */
   #define OOC_BLKSIZE        (1L << 13)



   /* ===================================================================
    * UTILZ
    * =================================================================== */   

   /* Test/set leaf flag for node n */
   #define OOC_ISLEAF(n)      ((n) -> node.type)
   #define OOC_SETLEAF(n)     ((n) -> node.type = 1)
   
   /* Test/set branch bit for octant k in node n */
   #define OOC_ISBRANCH(n,k)  (!OOC_ISLEAF(n) && (n) -> node.branch & 1 << (k))
   #define OOC_SETBRANCH(n,k) ((n) -> node.branch |=  1 << (k))
   
   /* Return index to node's 1st kid */
   #define OOC_KID1(o,n)      (!OOC_ISLEAF(n) ? (o)->nodes + (n)->node.kid \
                                              : NULL)
                                 
   /* Get index to last node in octree; -1 indicates empty octree */
   #define OOC_ROOTIDX(o)     ((o) -> nodes && (o) -> numNodes > 0 \
                                 ? (o) -> numNodes - 1 : -1)
                                 
   /* Get octree root; this is at the *end* of the node array! */
   #define OOC_ROOT(o)        (OOC_ROOTIDX(o) >= 0 \
                                 ? (o) -> nodes + OOC_ROOTIDX(o) : NULL)
   
   /* Copy node to octree root */
   #define OOC_SETROOT(o,n)   (memcpy(OOC_ROOT(o), (n), sizeof(OOC_Node)))

   /* Set all node fields to 0 */   
   #define OOC_CLEARNODE(n)   (memset((n), 0, sizeof(OOC_Node)))

   /* Modify origin o for octant k of size s */
   #define OOC_OCTORIGIN(o,k,s)  ((o) [0] += (s) * ((k)      & 1), \
                                  (o) [1] += (s) * ((k) >> 1 & 1), \
                                  (o) [2] += (s) * ((k) >> 2 & 1))

   /* Get num records stored in node (also handles leaves and NULL pointer) */
   #define OOC_DATACNT(n) \
      ((n) ? OOC_ISLEAF(n) ? (n) -> leaf.num [0] + (n) -> leaf.num [1] + \
                             (n) -> leaf.num [2] + (n) -> leaf.num [3] + \
                             (n) -> leaf.num [4] + (n) -> leaf.num [5] + \
                             (n) -> leaf.num [6] + (n) -> leaf.num [7]   \
                           : (n) -> node.num \
           : 0)

   /* Shortcuts for morton code from key and data record */
   #define OOC_KEY2MORTON(k, o)  OOC_Key2Morton((k), (o) -> org, \
                                                (o) -> mortonScale)
   #define OOC_REC2MORTON(r, o)  OOC_Key2Morton((o) -> key(r), (o) -> org, \
                                                (o) -> mortonScale)



   /* ===================================================================
    * DATA TYPEZ
    * =================================================================== */   
                           
   /* Octree node */
   typedef struct {
      union {
         struct {
            /* Interior node (node.type = 0) with:
               node.kid    = index to 1st child node in octree array (a.k.a
                             "Number One Son"), immediately followed by its
                             nonzero siblings as indicated by branch
               node.num    = num records stored in this subtree,
               node.branch = branch bitmask indicating nonzero kid nodes */
            char           type  : 1;
            OOC_NodeIdx    kid   : OOC_NODEIDX_BITS;             
            OOC_DataIdx    num   : OOC_DATAIDX_BITS;
            uint8_t        branch;
            
            /* NOTE that we can't make the kid's index relative to parent
             * (which would be more compact), since we don't know the
             * parent's index a priori during the bottom-up construction! */
         } node;

         struct {
            /* Leaf node (leaf.type = node.type = 1 with:
               leaf.num [k] = num records stored in octant k */
            char           type  : 1;
            OOC_OctCnt     num   [8];
         } leaf;
      };
   } OOC_Node;



   /* Top level octree container */
   typedef struct {
      FVECT          org, bound;    /* Cube origin (min. XYZ), size, and
                                       resulting bounds (max. XYZ) */
      RREAL          size,
                     *(*key)(const void*); /* Callback for data rec coords */
      RREAL          mortonScale;   /* Scale factor for generating Morton
                                       codes from coords */
      OOC_Node       *nodes;        /* Pointer to node array */
                                    /* *************************************
                                     * **** NODES ARE STORED IN POSTFIX ****
                                     * **** ORDER, I.E. ROOT IS LAST!   ****
                                     * *************************************/
      OOC_DataIdx    numData;       /* Num records in leaves */
      OOC_NodeIdx    maxNodes,      /* Num currently allocated nodes */
                     numNodes;      /* Num used nodes (<= maxNodes) */
      unsigned       recSize,       /* Size of data record in bytes */
                     leafMax,       /* Max #records per leaf (for build) */
                     maxDepth;      /* Max subdivision depth (for build) */
      FILE           *leafFile;     /* File containing data in leaves */
      OOC_Cache      *cache;        /* I/O cache for leafFile */
   } OOC_Octree;

   /* DOAN' FORGET: NODES IN ***POSTFIX*** ORDAH!!! */



   /* ===================================================================
    * FUNKZ
    * =================================================================== */   

   /* Init empty octree */
   void OOC_Null (OOC_Octree *oct);

    
   /* Initialises octree for records of size recSize with 3D coordinates,
    * accessed via key() callback.  The octree's bounding box is defined by
    * its origin org and size per axis, and must contain all keys. The
    * octree is associated with the specified leaf file containing the
    * records in Morton order */
   void OOC_Init (OOC_Octree *oct, unsigned recSize, FVECT org, RREAL size, 
                  RREAL *(*key)(const void*), FILE *leafFile);


   /* Allocate new octree node, enlarge array if necessary.  Returns pointer
      to new node or NULL if failed.  */
   OOC_Node *NewNode (OOC_Octree *oct);


   /* Reads indexed data record from leaf file and copies it to rec, else
    * returns -1 on failure */
   int OOC_GetData (OOC_Octree *oct, OOC_DataIdx idx, void *rec);
   

   /* Appends octree nodes to specified file along with metadata. Uses
    * RADIANCE's portable I/O routines. Returns 0 on success, else -1. */
   int OOC_SaveOctree (const OOC_Octree *oct, FILE *nodeFile);


   /* Load octree nodes and metadata from nodeFile and associate with
    * records in leafFile.  Uses RADIANCE's portable I/O routines.  Returns
    * 0 on success, else -1.  */
   int OOC_LoadOctree (OOC_Octree *oct, FILE *nodeFile, 
                       RREAL *(*key)(const void*), FILE *leafFile);

#ifdef DEBUG_OOC
   /* Traverse tree & check for valid nodes; oct -> leafFile must be open;
    * return 0 if ok, otherwise -1 */
   int OOC_Check (OOC_Octree *oct, const OOC_Node *node, 
                  FVECT org, RREAL size, OOC_DataIdx dataIdx);   
#endif

   /* Check whether key lies within bounding box defined by org and size */
   int OOC_InBBox (const OOC_Octree *oct, const FVECT org, RREAL size, 
                   const FVECT key); 

   /* Return index of octant containing key and corresponding origin if
    * nuOrg = NULL, or -1 if key lies outside all octants. Size is already
    * assumed to be halved, i.e. corresponding to the length of an octant
    * axis. */
   int OOC_Branch (const OOC_Octree *oct, const FVECT org, RREAL size, 
                   const FVECT key, FVECT nuOrg);                  

   /* Optimised version of OOC_Branch() with precomputed key Morton code */
   int OOC_BranchZ (const OOC_Octree *oct, const FVECT org, RREAL size, 
                    OOC_MortonIdx keyZ, FVECT nuOrg);                  

   /* For a given octree node, return the sum of the data counters for kids
    * [0..k-1].  On return, the node either points to the k'th kid, or
    * NULL if the kid is nonexistent or the node is a leaf */
   OOC_DataIdx OOC_GetKid (const OOC_Octree *oct, OOC_Node **node, 
                           unsigned k);

   /* Self-destruct */
   void OOC_Delete (OOC_Octree *oct);

   /* DID WE MENTION NODES ARE IN POSTFIX ORDAH? */
#endif
