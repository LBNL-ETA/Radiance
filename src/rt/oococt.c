/* 
   ==================================================================
   Out-of-core octree data structure
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts   
   ==================================================================
   
   $Id: oococt.c,v 2.1 2015/02/24 19:39:26 greg Exp $
*/


#include <stdlib.h>
#include "oococt.h"
#include "oocsort.h"


static OOC_Sort_Key zInterleave (OOC_Sort_Key k)
/* Interleave lower 21 bits of k with 00, resulting in 63 bits. 
   This code taken from 
   http://www.forceflow.be/2013/10/07/
      morton-encodingdecoding-through-bit-interleaving-implementations/ */
{
   /* Mask out lower 21 bits, the rest is cryptic */
   OOC_Sort_Key i = k & 0x1fffff;   
   
   i = (i | i << 32) & 0x001f00000000ffff;
   i = (i | i << 16) & 0x001f0000ff0000ff;  
   i = (i | i << 8)  & 0x100f00f00f00f00f;
   i = (i | i << 4)  & 0x10c30c30c30c30c3;
   i = (i | i << 2)  & 0x1249249249249249;
   
   return i;
}


static OOC_Sort_Key Key2zOrder (OOC_Octree *oct, FVECT key)
/* Compute 63-bit Morton (Z-order) index for 3D key */
{
   int i;
   OOC_Sort_Key k [3];
      
   /* Normalise key and map to 21-bit int */
   for (i = 0; i < 3; i++)
      k [i] = (((OOC_Sort_Key)1 << 21) - 1) / oct -> size * 
              (key [i] - oct -> org [i]);
   
   /* Interleave each dim with zeros and merge */
   return zInterleave(k [0]) | zInterleave(k [1]) << 1 | 
          zInterleave(k [2]) << 2;
}

      
static OOC_Node *NewNode (OOC_Octree *oct)
/* Allocate new octree node, enlarge array if necessary.
   Return pointer to new node or NULL if failed. */
{
   OOC_Node *n = NULL;
   
   if (++oct -> numNodes > oct -> maxNodes) {
      /* Need to allocate root or enlarge array */
      oct -> maxNodes += OOC_BLKSIZE / sizeof(OOC_Node);
      realloc(oct -> root, oct -> maxNodes * sizeof(OOC_Node));
      if (!oct -> root)
         return NULL;
   }
   
   n = oct -> root + oct -> numNodes - 1;
   n -> num = n -> kid = 0;
}
      

int OOC_Init (OOC_Octree *oct, unsigned recSize, unsigned leafMax, 
              unsigned maxDepth)
{
   ooc -> org [0] = oct -> org [1] = oct -> org [2] = FHUGE;
   ooc -> size = oct -> maxNodes = oct -> lastNode = 0;
   ooc -> recSize = recSize;
   ooc -> leafMax = leafMax;
   ooc -> maxDepth = maxDepth;
   ooc -> root = NULL;

   return NewNode(oct) ? 0 : -1;
}


int OOC_Branch (FVECT key, FVECT org, RREAL *size);
/* Return index of kid containing key, update origin and size accordingly */
{
   size *= 0.5;
   int kidIdx = 0;
   
   for (int i = 0; i < 3; i++) {
      RREAL cent = org [i] + size;
      if (key [i] > cent) {
         kidIdx |= 1 << i;
         org [i] = cent;
      }
   }
   
   return kidIdx;
}


OOC_Idx OOC_Find (OOC_Octree *oct, FVECT key)
/* Iterative search for given key */
{
   FVECT    org;
   RREAL    size;
   OOC_Idx  dataIdx = 0;
   OOC_Node *node = NULL;
   
   if (!oct || !oct -> root)
      return -1;
   
   node = oct -> root;
   VCOPY(org, oct -> org);
   size = oct -> size;
   
   /* Quit at terminal/empty node or internal leaf */
   while (OOC_INTNODE(*node)) {
      /* Determine kid to branch into and update data index with skipped
       * kids' data counts */
      int kidIdx = OOC_Branch(key, org, size);
      node = oct -> root + OOC_KID(*node);
      while (kidIdx--)
         dataIdx += OOC_NUMDATA(node++);
   }
   
   /* Reached empty node -> not found */
   if (OOC_EMPTYNODE(*node))
      return -1;
   
   /* Reached terminal node -> find leaf to branch into and update data
    * index with skipped leaves' data counts */
   if (OOC_TERMNODE(*node)) {
      OOC_Leaf *leaf = ooc -> leaves [OOC_KID(*node)];
      int leafIdx = OOC_Branch(key, org, size);
      while (leafIdx--)
         dataIdx += *leaf++;
   }

   /* (Reached internal leaf -> already have data index) */
   return dataIdx;
}


#if 0 
OOC_Node *FindNode (OOC_Octree *oct, FVECT key, OOC_Idx *dataIdx, int flags)
/* ITERATIVE node search for given key */
{
   FVECT    org;
   double   size;
   OOC_Idx  dataIdx = 0;
   OOC_Node *node = NULL;
   
   if (!oct || !oct -> root)
      return NULL;
   
   node = oct -> root;
   VCOPY(org, oct -> org);
   size = oct -> size;
   
   /* Quit at terminal or empty node */
   while (!OOC_TERM(*node) && !OOC_EMPTY(*node)) {
      size *= 0.5;
      int kidIdx = 0;
      
      /* Determine kid to branch into and update origin accordingly */
      for (int i = 0; i < 3; i++) {
         RREAL cent = org [i] + size;
         if (key [i] > cent) {
            kidIdx |= 1 << i;
            org [i] = cent;
         }
      }
       
      /* Increment node's data counter if inserting */
      if (flags & OOC_INSERT) 
         node -> numData++;
           
      /* Update data index with skipped kids' data counts */
      node = oct -> root + node -> kid;
      while (kidIdx--)
         dataIdx += node++ -> numData;
   }
   
   return node;
}

 
OOC_Idx FindLeaf (OOC_Leaf *leaf, FVECT key, OOC_Idx *dataIdx);
/* Return index into data block for leaf containing key */
{
   int kidIdx = 0;
   OOC_Idx = 
   
   for (int i = 0; i < 3; i++) {
      if (key [i] 
#endif
