/////////////////////////////////////////////////////////////////////////////
//  									   //
// Copyright (C) 1992-1995 by Alex Keller (keller@informatik.uni-kl.de)	   //
// All rights reserved							   //
//									   //
// This software may be freely copied, modified, and redistributed	   //
// provided that this copyright notice is preserved on all copies.	   //
// 									   //
// You may not distribute this software, in whole or in part, as part of   //
// any commercial product without the express consent of the authors.	   //
// 									   //
// There is no warranty or other guarantee of fitness of this software	   //
// for any purpose.  It is provided solely "as is".			   //
//									   //
/////////////////////////////////////////////////////////////////////////////

int Face(int ac, char ** av)
{
   Point Center;
   Vector E1, E2, E3, nn, ConvexNormal;
   int i, max, Corners, p0, p1, p2, n, Start;
   double maxd, d, x, y;
   Boolean * Out, Empty;
   C_VERTEX * Pnt;
   FVECT P;
   Vertex * Points;
   Surface * CurrentSurface;
   Texture * Luminaire;
   Color Tone;
   
   if(ac--)
      av++;
   
   if(Zero(c_cmaterial->ed))
      Luminaire = NULL;
   else
   {
      x = c_cmaterial->ed_c.cx;
      y = c_cmaterial->ed_c.cy;
      Tone.fromCIE(x, y, 1.0 - x - y);
      Tone *= c_cmaterial->ed / (2 * M_PI);
      Luminaire = new Monoton(Tone);
   }
   
   x = c_cmaterial->rd_c.cx;
   y = c_cmaterial->rd_c.cy;
   Tone.fromCIE(x, y, 1.0 - x - y);
   Tone *= c_cmaterial->rd;
   CurrentSurface = new Ward(new Monoton(Tone), c_cmaterial->rd,
			     NULL, c_cmaterial->rs, c_cmaterial->rs_a, c_cmaterial->rs_a,
			     NULL, c_cmaterial->td,
			     NULL, c_cmaterial->ts, c_cmaterial->ts_a, c_cmaterial->ts_a,
			     Luminaire);
   
   Points = new Vertex[ac];
   n = 0;   
   
   while(ac-- > 0)
   {
      Pnt = c_getvert(*av++);
      
      //if(Pnt == NULL)
      //	 cout << "Error" << endl;
      
      xf_xfmpoint(P, Pnt->p);
      
      Points[n].Origin.x = P[0];
      Points[n].Origin.y = P[1];
      Points[n].Origin.z = P[2];
      
      n++;
   }
   
   if(n == 3)
      *World << new Triangle(*World, &Points[0], &Points[1], &Points[2], CurrentSurface);
   else if(n > 3)
   {
      Center = Points[0].Origin;
      
      for(i = 1; i < n; i++)
      {
	 Center.x += Points[i].Origin.x;
	 Center.y += Points[i].Origin.y;
	 Center.z += Points[i].Origin.z;
      }
      
      Center.x /= (double) n;
      Center.y /= (double) n;
      Center.z /= (double) n;
      
      maxd = Length(Center - Points[0].Origin);
      max = 0;
      
      for(i = 1; i < n; i++)
      {
	 d = Length(Center - Points[i].Origin);
	 
	 if(d > maxd)
	 {
	    max = i;
	    maxd = d;
	 }
      }
      
      Out = new Boolean[n];
      
      for(i = 0; i < n; i++)
	 Out[i] = False;
      
      p1 = max;
      p0 = p1 - 1;
      
      if(p0 < 0)
	 p0 = n - 1;
      
      p2 = p1 + 1;
      
      if(p2 == n)
	 p2 = 0;
      
      ConvexNormal = (Points[p2].Origin - Points[p0].Origin) | (Points[p1].Origin - Points[p0].Origin);
      Normalize(ConvexNormal);
      
      Corners = n;
      p0 = - 1;
      
      while(Corners >= 3)
      {
	 Start = p0;
	 
	 do
	 {
	    p0 = (p0 + 1) % n;
	    
	    while(Out[p0])
	       p0 = (p0 + 1) % n;
	    
	    p1 = (p0 + 1) % n;
	    
	    while(Out[p1])
	       p1 = (p1 + 1) % n;
	    
	    p2 = (p1 + 1) % n;
	    
	    while(Out[p2])
	       p2 = (p2 + 1) % n;
	    
	    if(p0 == Start)
	       break;
	    
	    nn = (Points[p2].Origin - Points[p0].Origin) | (Points[p1].Origin - Points[p0].Origin);
	    Normalize(nn);
	    d = Length(nn - ConvexNormal);
	    
	    E1 = nn | (Points[p1].Origin - Points[p0].Origin);
	    E2 = nn | (Points[p2].Origin - Points[p1].Origin);
	    E3 = nn | (Points[p0].Origin - Points[p2].Origin);
	    
	    Empty = True;
	    
	    for(i = 0; i < n; i++)
	       if(! Out[i])
		  if((i != p0) && (i != p1) && (i != p2))
		  {
		     Empty = Empty && ! ((E1 * (Points[i].Origin - Points[p0].Origin) <= - Epsilon)
					 && (E2 * (Points[i].Origin - Points[p1].Origin) <= - Epsilon)
					 && (E3 * (Points[i].Origin - Points[p2].Origin) <= - Epsilon));
		  }
	 }
	 while((d > 1.0) || (! Empty));
	 
	 if(p0 == Start)
	 {
	    cout << "misbuilt polygonal face..." << endl;
	    break;
	 }
	 
	 *World << new Triangle(*World, &Points[p0], &Points[p1], &Points[p2], CurrentSurface);
	 
	 Out[p1] = True;
	 Corners--;
      }
      
      delete [] Out;
   }
   
   return(MG_OK);
}

/////////////////////////////////////////////////////////////////////////////
