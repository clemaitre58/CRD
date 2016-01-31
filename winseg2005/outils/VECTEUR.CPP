#include "vecteur.h"
#include <math.h>
   

 double   produit_scalaire(point P1,point P2,point P3)
{
 point p1,p2;
 p1.X = P2.X-P1.X; p1.Y = P2.Y-P1.Y ; p1.Z = P2.Z-P1.Z;
 p2.X = P3.X-P1.X; p2.Y = P3.Y-P1.Y ; p2.Z = P3.Z-P1.Z;
 return  p1.X*p2.X+p1.Y*p2.Y+p1.Z*p2.Z;
}


 point    produit_vectoriel(point P1,point P2,point P3)
{
 point p1,p2,p;

 p1.X = P2.X-P1.X; p1.Y = P2.Y-P1.Y ; p1.Z = P2.Z-P1.Z;
 p2.X = P3.X-P1.X; p2.Y = P3.Y-P1.Y ; p2.Z = P3.Z-P1.Z;
 p.X = p1.Y*p2.Z-p1.Z*p2.Y;
 p.Y = p1.Z*p2.X-p1.X*p2.Z;
 p.Z = p1.X*p2.Y-p1.Y*p2.X;
 return p;
}


 double   norme(point P)
{
  return sqrt( pow(P.X,(double)2.0) + pow(P.Y,(double)2.0) + pow(P.Z,(double)2.0) );
}
      


 double   distance_eucli(point P1,point P2)
{
 double d1,d2,d3;
 d1 = pow((P1.X-P2.X),(double)2.0);
 d2 = pow((P1.Y-P2.Y),(double)2.0);
 d3 = pow((P1.Z-P2.Z),(double)2.0);
 return sqrt( d1 + d2 + d3 );
}

 double  surface(point P1,point P2,point P3)
{
 return norme(produit_vectoriel(P1,P2,P3));
}

