#ifndef lint
static const char RCSid[] = "$Id: oocmorton.c,v 2.2 2017/08/14 21:12:10 rschregle Exp $";
#endif


/* 
   =========================================================================
   Routines to generate and compare Morton Codes, i.e. indices on space
   filling Z-curves.
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com}) 
   (c) Lucerne University of Applied Sciences and Arts, 
       supported by the Swiss National Science Foundation (SNSF, #147053)
   =========================================================================
   
   $Id: oocmorton.c,v 2.2 2017/08/14 21:12:10 rschregle Exp $
*/


#if !defined(_WIN32) && !defined(_WIN64) || defined(PMAP_OOC)
/* No Windoze support for now */

#include "oocmorton.h"



OOC_MortonIdx OOC_Key2Morton (const FVECT key, const FVECT org, RREAL scale)
/* Compute Morton code (Z-curve index) of length 3 * OOC_MORTON_BITS bits
 * for 3D key within bounding box defined by org and scaled to maximum index
 * with scale */
{
   unsigned       i;
   OOC_MortonIdx  k [3];

   /* Normalise key and map each dim to int of OOC_MORTON_BITS */
   for (i = 0; i < 3; i++)
      k [i] = scale * (key [i] - org [i]);
                  
   /* Interleave each dim with zeros and merge */
   return OOC_BitInterleave(k [0]) | OOC_BitInterleave(k [1]) << 1 |
          OOC_BitInterleave(k [2]) << 2;
}

#endif /* NIX / PMAP_OOC */
