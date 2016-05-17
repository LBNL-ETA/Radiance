/* 
   =======================================================================
   Header for building out-of-core octree data structure
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   =======================================================================
   
   $Id: oocbuild.h,v 2.1 2016/05/17 17:39:47 rschregle Exp $
*/


#ifndef OOC_BUILD_H
   #define OOC_BUILD_H

   /* Bottom-up construction of out-of-core octree in postorder traversal.
    * The octree oct is assumed to be initialised with its origin (oct ->
    * org), size (oct -> size), key callback (oct -> key), and its
    * associated leaf file (oct -> leafFile).  
    
    * Records are read from the leafFile and assumed to be sorted in
    * Z-order, which defines an octree leaf ordering.  Leaves (terminal
    * nodes) are constructed such that they contain <= leafMax records and
    * have a maximum depth of maxDepth.
    
    * Note that the following limits apply:
    *    leafMax  <= OOC_OCTCNT_MAX     (see oococt.h)
    *    maxDepth <= OOC_MORTON_BITS    (see oocsort.h)
    
    * The maxDepth limit arises from the fact that the Z-ordering has a
    * limited resolution and will map node coordinates beyond a depth of
    * OOC_MORTON_BITS to the same Z-index, causing nodes to be potentially
    * read out of sequence and winding up in the wrong octree nodes.
     
    * On success, the octree pointer oct is returned, with the constructed
    * nodes in oct -> nodes, and the node count in oct -> numNodes.  On
    * failure, NULL is returned.  */
    OOC_Octree *OOC_Build (OOC_Octree *oct, unsigned leafMax, 
                           unsigned maxDepth);
#endif
    