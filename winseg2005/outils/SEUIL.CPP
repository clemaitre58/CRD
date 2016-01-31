#include "seuil.h"  

 void  Seuillage1(SImage im, unsigned char s1,BOOL flag)
{   
  /* le flag a TRUE effectue une binarisation 0 - 255
   sinon une binarisation 0 - idem
  */ 
 long l;           
 unsigned char tabseuil[256];
 for (int k=0;k<=s1;k++) tabseuil[k] = 0x00;
 for (k=s1+1;k<=255;k++) tabseuil[k] =  0xFF;
 for (long j=im.Hmin;j<=im.Hmax;j++)
 { 
   	l = j*im.Larg+im.Lmin;
   	for (long i=im.Lmin;i<=im.Lmax;i++) 
   	{ 	
   		if (flag) im.Image[l] = tabseuil[im.Image[l]];  
   		else if (im.Image[l]<s1) im.Image[l] =  tabseuil[im.Image[l]]; 
   		l++;
   	}
 } 
}                                  

 void  Seuillage2(SImage im, unsigned char s1, unsigned char s2)
{
 long l;           
 unsigned char tabseuil[256];
 for (int k=0;k<s1;k++) tabseuil[k] = 0x00;
 for (k=s2;k<=255;k++) tabseuil[k] =  0x00;
 for (k=s1+1;k<=s2;k++) tabseuil[k] =  0xFF;
 for (long j=im.Hmin;j<=im.Hmax;j++)
 { 
   	l = j*im.Larg+im.Lmin;
   	for (long i=im.Lmin;i<=im.Lmax;i++)  im.Image[l++] =  tabseuil[im.Image[l]];  
 } 
}


