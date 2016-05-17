/* 
   ======================================================================
   Routines to generate and compare Morton Codes, i.e. indices on space
   filling Z-curves.
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com}) 
   (c) Lucerne University of Applied Sciences and Arts, 
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================   
   
   $Id: oocmorton.h,v 2.1 2016/05/17 17:39:47 rschregle Exp $
*/

#ifndef OOC_MORTON_H
   #define OOC_MORTON_H
   
   #include "fvect.h"
   #include <stdint.h>

   
   /* Type to represent Morton codes (Z-curve indices) as sort keys */
   typedef uint64_t OOC_MortonIdx;   
   
   
   /* Number of bits per dimension for Morton code (Z-curve index) used to
    * sort records.  This corresponds to the maximum number of bounding box
    * subdivisions which can be resolved per axis; further subdivisions will
    * map to the same Morton code, thus invalidating the sort! */
   #define OOC_MORTON_BITS 21
   #define OOC_MORTON_MAX  (((OOC_MortonIdx)1 << OOC_MORTON_BITS) - 1)


   /* Interleave lower OOC_MORTON_BITS bits of k with 00, resulting in
    * 3*OOC_MORTON_BITS bits. Optimised "bitmask hack" version. This code
    * taken from: 
    * http://www.forceflow.be/2013/10/07/
    *   morton-encodingdecoding-through-bit-interleaving-implementations/ */
   #define OOC_BitInterleave(k) \
      ((k) &= OOC_MORTON_MAX, \
       (k) = ((k) | (k) << 32) & 0x001f00000000ffff, \
       (k) = ((k) | (k) << 16) & 0x001f0000ff0000ff, \
       (k) = ((k) | (k) << 8)  & 0x100f00f00f00f00f, \
       (k) = ((k) | (k) << 4)  & 0x10c30c30c30c30c3, \
       (k) = ((k) | (k) << 2)  & 0x1249249249249249)

   
   /* Compute Morton code (Z-curve index) of length 3 * OOC_MORTON_BITS bits
    * for 3D key within bounding box defined by org and scaled to maximum
    * index with scale */
   OOC_MortonIdx OOC_Key2Morton (const FVECT key, const FVECT org, 
                                 RREAL scale);
#endif
 