#include "inverse.h"


 void  Inversion(SImage im)
{                                                      
   long l;
   for (long j=im.Hmin;j<=im.Hmax;j++)
   { 
    	l = j*im.Larg+im.Lmin;
    	for (long i=im.Lmin;i<=im.Lmax;i++) im.Image[l++] =  255-im.Image[l];  
   } 

}
