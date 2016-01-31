#include <windows.h>
                                                      
struct point
{
	double X;
	double Y;
	double Z;
};
       
 point    produit_vectoriel(point P1,point P2,point P3) ;
 double   norme(point P);
 double   produit_scalaire(point P1,point P2,point P3);
 double   distance_eucli(point P1,point P2);
 double  surface(point P1,point P2,point P3);