#include "color.h"
#include <windows.h>  
#include <windowsx.h>  


void SeparationRGB(SImageRGB ImageRGB, SImage r, SImage v, SImage b)	
{
	rgb RGB;
	long l = 0;
    unsigned char   * ImR = (unsigned char   *) r.Image;
    unsigned char   * ImV = (unsigned char   *) v.Image;
    unsigned char   * ImB = (unsigned char   *) b.Image;
	for (long j= ImageRGB.Hmin;j<=ImageRGB.Hmax;j++) 
	{                              
	   	for (long i=ImageRGB.Lmin;i<=ImageRGB.Lmax;i++) 
	    {       
			RGB = ImageRGB.Image[j*ImageRGB.Larg+i];
			((unsigned char   * )ImR)[l] = RGB.r; 
			((unsigned char   * )ImV)[l] = RGB.g; 
			((unsigned char   * )ImB)[l] = RGB.b; 
			l++;
		} 
	}	

}
void LuminanceChrominance(SImageRGB ImageRGB, SImage Y, SImage C0, SImage C1)
{   
	rgb RGB;                                                        
	long l = 0;
    int lum;
    long i,j;

    unsigned char   * ImY = (unsigned char   *) Y.Image;   
    unsigned char   * ImC0 = (unsigned char   *) C0.Image;   
    unsigned char   * ImC1 = (unsigned char   *) C1.Image;   
    
    int  * Chromi0  = (int  * ) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,(ImageRGB.Lmax-ImageRGB.Lmin+1)*(ImageRGB.Hmax-ImageRGB.Hmin+1)*sizeof(int));
    int  * Chromi1  = (int  *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,(ImageRGB.Lmax-ImageRGB.Lmin+1)*(ImageRGB.Hmax-ImageRGB.Hmin+1)*sizeof(int));
      
    int MaxRouge=-255,MaxBleu=-255;
    int MinRouge=255,MinBleu=255;
    for (j= ImageRGB.Hmin;j<=ImageRGB.Hmax;j++) 
	{                              
		for (i=ImageRGB.Lmin;i<=ImageRGB.Lmax;i++) 
		{                              
	    		RGB = ImageRGB.Image[j*ImageRGB.Larg+i]; 
	    		lum =  (int)(0.3*RGB.r + 0.59*RGB.g + 0.11*RGB.b);
				((unsigned char   * )ImY)[l] = (unsigned char) lum; 
				Chromi0[l] = (int)RGB.r-lum;
				Chromi1[l] = (int)RGB.b-lum;
				if (Chromi0[l] 	> MaxRouge) MaxRouge =Chromi0[l];
				if (Chromi1[l]  > MaxBleu) MaxBleu = Chromi1[l] ;
				if (Chromi0[l]  < MinRouge) MinRouge = Chromi0[l] ;
				if (Chromi1[l]  < MinBleu) MinBleu = Chromi1[l] ;
				l++;
		}
	} 
	l = 0;   
	// normalisation des chrominances 
	int normR,normB,norm; 
	if (MaxRouge==MinRouge) normR= 1; else normR = (MaxRouge-MinRouge)/255;
	if (MaxBleu==MinBleu) normB= 1; else normB = (MaxBleu-MinBleu)/255;      
	if (normR>normB) norm = normR; else norm = normB;   
	for (j= ImageRGB.Hmin;j<=ImageRGB.Hmax;j++) 
	{                              
		for (i=ImageRGB.Lmin;i<=ImageRGB.Lmax;i++) 
		{                              
	    		RGB = ImageRGB.Image[j*ImageRGB.Larg+i]; 
				((unsigned char   * )ImC0)[l] = (unsigned char) (Chromi0[l]-MinRouge)/norm;
				((unsigned char   * )ImC1)[l] = (unsigned char) (Chromi1[l] -MinBleu)/norm;
				l++;
		}
	}

	GlobalFreePtr(Chromi0);	
	GlobalFreePtr(Chromi1);	
}