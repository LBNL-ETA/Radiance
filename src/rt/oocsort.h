/* 
   ======================================================================
   Header for N-way out-of-core merge sort for records with 3D keys.
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================   
   
   $Id: oocsort.h,v 2.3 2016/05/17 17:39:47 rschregle Exp $
*/

#ifndef OOC_SORT_H
   #define OOC_SORT_H
   
   #include "fvect.h"
   #include <stdio.h>   
   

   /* Sort records in inFile and append to outFile by subdividing inFile
    * into small blocks, sorting these in-core, and merging them out-of-core
    * via a priority queue.
    *
    * This is implemented as a recursive (numBlk)-way sort; the input is
    * successively split into numBlk smaller blocks until these are of size
    * <= blkSize, i.e.  small enough for in-core sorting, then merging the
    * sorted blocks as the stack unwinds.  The in-core sort is parallelised
    * over numProc processes.
    *
    * Parameters are as follows:
    * inFile      Opened input file containing unsorted records
    * outFile     Opened output file containing sorted records
    * numBlk      Number of blocks to divide into / merge from 
    * blkSize     Max block size and size of in-core sort buffer, in bytes 
    * numProc     Number of parallel processes for in-core sort 
    * recSize     Size of input records in bytes
    * bbOrg       Origin of bounding box containing record keys for Morton code
    * bbSize      Extent of bounding box containing record keys for Morton code
    * key         Callback to access 3D coords from records for Morton code
    */                                  
   int OOC_Sort (FILE *inFile, FILE *outFile, unsigned numBlk, 
                 unsigned long blkSize, unsigned numProc, unsigned recSize,
                 FVECT bbOrg, RREAL bbSize, RREAL *(*key)(const void*));
#endif
