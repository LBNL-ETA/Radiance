#ifndef lint
static const char	RCSid[] = "$Id: sm_samp.c,v 3.9 2003/06/30 14:59:12 schorsch Exp $";
#endif
/*
 * sm_samp.c
 */

#include <string.h>

#include "standard.h"
#include "sm_flag.h"
#include "rhd_sample.h"

SAMP rsL;
int32 *samp_flag=NULL;

/* Each sample has a world coord point, and direction, brightness,chrominance,
   and RGB triples
 */

#define TMSIZE sizeof(TMbright)
#define SAMPSIZ (3*sizeof(SFLOAT)+sizeof(int32)+ 6*sizeof(BYTE) + TMSIZE + 2*sizeof(int))

 /* Extra points world space point, vert flag and qt flag */
#define POINTSIZ (3*sizeof(SFLOAT) + 2*sizeof(int))

/* Clear the flags for all samples */ 
sClear_all_flags(s)
SAMP *s;
{
  if(samp_flag)
    memset((char *)samp_flag, '\0', FLAG_BYTES(S_MAX_BASE_PT(s)));
}

sInit(s)
   SAMP *s;
{
  S_REPLACE_SAMP(s) = 0;
  S_NUM_SAMP(s) = 0;
  S_TONE_MAP(s) = 0;
  S_FREE_SAMP(s) = -1;
  sClear_base_points(s);

}

sFree(s)
SAMP *s;
{
    if(S_BASE(s))
    {
	free(S_BASE(s));
	S_BASE(s)=NULL;
    }
    if(samp_flag)
    {
	free(samp_flag);
	samp_flag = NULL;
    }
}
	/* Initialize/clear global smL sample list for at least n samples */
SAMP
*sAlloc(nptr,extra_points)
int *nptr,extra_points;
{
  unsigned nbytes;
  register unsigned	i;
  SAMP *s;
  int n;
  
  s = &rsL;
  /* round space up to nearest power of 2 */
  nbytes = (*nptr)*SAMPSIZ + extra_points*POINTSIZ + 8;
  for (i = 1024; nbytes > i; i <<= 1)
		;
  /* get number of samples that fit in block */
  n = (i - 8 - extra_points*POINTSIZ) / SAMPSIZ;	

  /* Must make sure n + extra_points can fit in a S_ID identifier */
  if ( n > (S_ID_MAX - extra_points))
    n = (S_ID_MAX - extra_points);

  S_BASE(s) = (char *)malloc(n*SAMPSIZ + extra_points*POINTSIZ);
  if (!S_BASE(s))
    error(SYSTEM,"sAlloc(): Unable to allocate memory");

  /* assign larger alignment types earlier */
  S_W_PT(s) = (SFLOAT(*)[3])S_BASE(s);
  S_W_DIR(s) = (int32 *)(S_W_PT(s) + n + extra_points);
  S_BRT(s) = (TMbright *)(S_W_DIR(s) + n);
  S_CHR(s) = (BYTE (*)[3])(S_BRT(s) + n);
  S_RGB(s) = (BYTE (*)[3])(S_CHR(s) + n);
  S_INFO(s) = (int *)(S_RGB(s)+n);
  S_INFO1(s) = (int *)(S_INFO(s)+n+extra_points);
  S_MAX_SAMP(s) = n;
  S_MAX_BASE_PT(s) = n + extra_points;

  /* Allocate memory for a per/sample bit flag */
  if(!(samp_flag = (int32 *)malloc(FLAG_BYTES(n+extra_points))))
    error(SYSTEM,"sAlloc(): Unable to allocate flag memory");
  sInit(s);
  sClear_all_flags(s);

  *nptr = n;
  
  return(s);
}


/* Add a base point to the sample structure: only the world point
   is added: These points are not displayed-they are used to form the
   initial mesh
 */
S_ID
sAdd_base_point(s,v)
     SAMP *s;
     FVECT v;
{
    S_ID id;

    /* Get pointer to next available point */
    id = S_NEXT_BASE_PT(s);
    if(id >= S_MAX_BASE_PT(s))
       return(-1);

    /* Update free aux pointer */
    S_NEXT_BASE_PT(s)++;

    /* Copy vector into point */
    VCOPY(S_NTH_W_PT(s,id),v);
    return(id);
}

/* Copy the values of sample n_id into sample id: Note must not call with
   Base point n_id
*/
int
sCopy_samp(s,n_id,id)
   SAMP *s;
   S_ID n_id,id;
{

#ifdef DEBUG
   if(id <0||id >= S_MAX_SAMP(s)||n_id <0||n_id >= S_MAX_SAMP(s))
   {
     eputs("smReset_sample():invalid sample id\n");
     return(0);
   }
#endif
   if(n_id == id)
     return(1);
 /* Copy vector into point */
   VCOPY(S_NTH_W_PT(s,n_id),S_NTH_W_PT(s,id));
   S_NTH_W_DIR(s,n_id) = S_NTH_W_DIR(s,id);

   S_NTH_BRT(s,n_id) = S_NTH_BRT(s,id);
   S_NTH_CHR(s,n_id)[0] = S_NTH_CHR(s,id)[0];
   S_NTH_CHR(s,n_id)[1] = S_NTH_CHR(s,id)[1];
   S_NTH_CHR(s,n_id)[2] = S_NTH_CHR(s,id)[2];

   return(1);
}


/* Initialize the sample 'id' to contain world point 'p', direction 'd' and
   color 'c'. If 'p' is NULL, a direction only is represented in 'd'. In
   this case, the world point is set to the projection of the direction on
   the view sphere, and the direction value is set to -1
 */
void
sInit_samp(s,id,c,d,p,o_id)
   SAMP *s;
   S_ID id;
   COLR c;
   FVECT d,p;
   S_ID o_id;
{

  if(o_id != INVALID)
    sCopy_samp(s,id,o_id);
  else
  {
    /* Copy vector into point */
    VCOPY(S_NTH_W_PT(s,id),p);
    if(d)
      S_NTH_W_DIR(s,id) = encodedir(d);
    else
      S_NTH_W_DIR(s,id) = -1;

    /* calculate the brightness and chrominance,RGB will be set by
       tonemapping
       */
    tmCvColrs(&S_NTH_BRT(s,id),S_NTH_CHR(s,id),c,1);
  }
    /* Set ACTIVE bit upon creation */
  S_SET_FLAG(id);
  if(id < S_TONE_MAP(s))
    tmMapPixels(S_NTH_RGB(s,id),&S_NTH_BRT(s,id),S_NTH_CHR(s,id),1);
}


/* Allocate a new sample. If an existing sample was replaced: set flag */
S_ID
sAlloc_samp(s,replaced)
   SAMP *s;
   int *replaced;
   
{
    S_ID id;

    /* First check if there are any freed sample available */
    if((id = S_FREE_SAMP(s)) != INVALID)
    {
      *replaced = 0;
#if 0
      fprintf(stderr,"allocating previously freed sample %d\n",id);
#endif
      S_FREE_SAMP(s) = S_NTH_NEXT(s,id);
      return(id);
    }
    /* If havn't reached end of sample array:return next sample */
    if(S_NUM_SAMP(s) < S_MAX_SAMP(s))
    {
      id = S_NUM_SAMP(s)++;
#if 0
      fprintf(stderr,"allocating sample %d\n",id);
#endif
      *replaced = 0;
      return(id);
    }

#ifdef DEBUG 
{
  static int replace = 0;
  if(replace == 0)
  {
    eputs("Out of samples: using replacement strategy\n");
    replace =1 ;
  }
}
#endif
    
    /* CLOCKED LRU replacement policy: Search through samples starting 
       with S_REPLACE_SAMP for a sample that does not have its active 
       bit set. Clear bit along the way 
     */
    id = S_REPLACE_SAMP(s);
    while(S_IS_FLAG(id))
    {
      S_CLR_FLAG(id);
      id = (id +1)% S_MAX_SAMP(s);
    }
    S_REPLACE_SAMP(s) = (id +1)% S_MAX_SAMP(s);
    *replaced = 1;
    return(id);
}

















