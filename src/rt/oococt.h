#ifndef OOC_OCT_H
#define OOC_OCT_H

/* 
   ==================================================================
   Out-of-core octree data structure
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts   
   ==================================================================
   
   $Id: oococt.h,v 2.1 2015/02/24 19:39:26 greg Exp $
*/


#include "fvect.h"


/* Pointer to i-th record */
#define OOC_RECPTR(p,i) ((p) + 

/* Default data block size (8k) */
#define OOC_BLKBITS     13
#define OOC_BLKSIZE     (1 << OOC_BLKBITS)

/* Index to data block and index within the block */
#define OOC_IDXTOBLK(i) ((i) >> OOC_BLKBITS)
#define OOC_IDXINBLK(i) ((i) & (1 << OOC_BLKBITS) - 1)

/* Check for block overflow when adding node
#define OOC_BLKFULL(n)  ((!(n) * sizeof(OOC_Node) % OOC_BLKSIZE) */

/* Check for internal node, empty node, interior leaf, or terminal node 
   (8 leaves) */
#define OOC_INTNODE(n)     ((n).num > 0 && (n).kid > 0)
#define OOC_EMPTYNODE(n)   (!(n).num)
#define OOC_INTLEAF(n)     ((n).kid < 0)
#define OOC_TERMNODE(n)    ((n).num < 0)

/* Return record count and kid/leaf index for node */
#define OOC_NUMREC(n)      abs((n).num)
#define OOC_KID(n)         abs((n).kid)

/* Check for loaded data block */
#define OOC_LOADED(b)      ((b).dataBlk != NULL)

/* (Mutually exclusive) operations on octree to be specified as flags */
#define OOC_INSERT      1
#define OOC_SEARCH      2

/* Counter for records at leaves (2 bytes) */
typedef unsigned short OOC_Leaf;

/* Node and record index type */
typedef int OOC_Idx; 

/* Octree node */
typedef struct {
   /* 
      For interior nodes (num > 0): 
         num = num records stored in this (sub)tree,
         kid = index to 1st child node (immediately followed by its
               siblings) in octree array.
      
      For interior leaves (kid < 0):
         num = num records stored in leaf,
         abs(kid) = index to single leaf, no siblings.
         
      For terminal nodes (num < 0):
         abs(num) = num records stored in leaves,
         kid = index to 1st of 8 leaves (immediately followed by its
               siblings) in leaf array.
   */
   OOC_Idx        num, kid;
} OOC_Node;


/* Load map entry; points to start of data block of size OOC_BLKSIZE, or
 * NULL if not loaded. LRU counter counts accesses for paging from disk */
typedef struct {
   void           *dataBlk;
   unsigned long  lruCnt;
} OOC_Load;


/* Top level octree container */
typedef struct {
   FVECT          org;        /* Cube origin (min. XYZ) and size */
   RREAL          size;
   OOC_Node       *root;      /* Pointer to node array */
   OOC_Idx        maxNodes,   /* Num currently allocated nodes */
                  numNodes;   /* Num used nodes (< maxNodes) */
   OOC_Load       *loadMap;   /* Pointer to load map */
   OOC_Leaf       *leaves;    /* Pointer to leaf array */
   unsigned       recSize,    /* Size of data record in bytes */
                  leafMax,    /* Max num records per leaf */
                  maxDepth;   /* Max subdivision depth */
} OOC_Octree;


/* Initialise octree for records of size recSize, return -1 on failure */
int OOC_Init (OOC_Octree *oct, unsigned recSize, unsigned leafMax);

/* Return index of kid containing key, update origin and size accordingly */
int OOC_Branch (FVECT key, FVECT org, RREAL *size);

/* Return data index to leaf containing key, or -1 if not found */
OOC_Idx OOC_Find (OOC_Octree *oct, FVECT key);

#if 0
/* Return pointer to terminal node containing key */
OOC_Node *FindNode (OOC_Octree *oct, FVECT key);

/* Return index into data block for leaf containing key */
OOC_Idx FindLeaf (OOC_Leaf *leaf, FVECT key);
#endif

#endif
