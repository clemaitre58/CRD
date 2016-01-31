#include "miroir.h"



 void  MiroirH(SImage im)
{                                                                                    
	DWORD T = im.Larg*im.Haut;
	HGLOBAL ic;
	ic = GlobalAlloc(GMEM_MOVEABLE,T);
	unsigned char  * tab;
	tab =  (unsigned char  *) GlobalLock(ic) ;  
   	for (long l=0;l<im.Larg*im.Haut;l++) tab[l]=im.Image[l]; 
	for (long j=im.Hmin;j<=im.Hmax;j++)
	{
			l = j*im.Larg+im.Lmin;  
			for (long i=im.Lmin;i<=im.Lmax;i++) im.Image[l++] = tab[im.Larg*(im.Hmax-j+im.Hmin)+i] ;
	}
	GlobalUnlock(ic);
	GlobalFree(ic); 	
}


 void  MiroirV(SImage im)
{
	DWORD T = im.Larg*im.Haut;
	HGLOBAL ic;
	ic = GlobalAlloc(GMEM_MOVEABLE,T);
	unsigned char  * tab;
	tab =  (unsigned char  *) GlobalLock(ic) ;  
   	for (long l=0;l<im.Larg*im.Haut;l++) tab[l]=im.Image[l]; 
	for (long j=im.Hmin;j<=im.Hmax;j++)
	{
			l = j*im.Larg+im.Lmin;   	
    		for (long i=im.Lmin;i<=im.Lmax;i++) im.Image[l++] = tab[im.Larg*j+(im.Lmax+im.Lmin-i)];
     }
	GlobalUnlock(ic);
	GlobalFree(ic); 
}		
