#include "rotation.h"

#include <stdlib.h> 

 void  RotationD(SImage im, SImage imRotD)
{
	long Fin_j= im.Hmax;
	div_t res = div((int)(im.Hmax-im.Hmin+1),4);   
	if (res.rem !=0) Fin_j -= res.rem;
	
	for (long i=im.Lmin;i<=im.Lmax;i++) 
	{	
		for(long j=im.Hmin;j<=Fin_j;j++)
			imRotD.Image[imRotD.Larg*imRotD.Haut-imRotD.Larg*(i-im.Lmin+1)+j-im.Hmin] = im.Image[j*im.Larg+i];
	}
}                                                                     
 void  RotationD(SImageFloat im, SImageFloat imRotD)
{
	long Fin_j= im.Hmax;
	div_t res = div((int)(im.Hmax-im.Hmin+1),4);   
	if (res.rem !=0) Fin_j -= res.rem;
	
	for (long i=im.Lmin;i<=im.Lmax;i++) 
	{	
		for(long j=im.Hmin;j<=Fin_j;j++)
			imRotD.Image[imRotD.Larg*imRotD.Haut-imRotD.Larg*(i-im.Lmin+1)+j-im.Hmin] = im.Image[j*im.Larg+i];
	}
}                                                                     

  void  RotationG(SImage im, SImage imRotG)
{  
	long Fin_j= im.Hmax;
	div_t res = div((int)(im.Hmax-im.Hmin+1),4);   
	if (res.rem !=0) Fin_j -= res.rem;
	
	for (long i=im.Lmin;i<=im.Lmax;i++) 
	{	
		for(long j=im.Hmin;j<=Fin_j;j++)
			imRotG.Image[imRotG.Larg-1+imRotG.Larg*(i-im.Lmin)-(j-im.Hmin)] = im.Image[j*im.Larg+i];
	}
}
  void  RotationG(SImageFloat im, SImageFloat imRotG)
{  
	long Fin_j= im.Hmax;
	div_t res = div((int)(im.Hmax-im.Hmin+1),4);   
	if (res.rem !=0) Fin_j -= res.rem;
	
	for (long i=im.Lmin;i<=im.Lmax;i++) 
	{	
		for(long j=im.Hmin;j<=Fin_j;j++)
			imRotG.Image[imRotG.Larg-1+imRotG.Larg*(i-im.Lmin)-(j-im.Hmin)] = im.Image[j*im.Larg+i];
	}
}