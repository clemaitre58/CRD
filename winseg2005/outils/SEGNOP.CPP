#include "segnop.h"
#include <math.h>

 void  Roberts(SImage im,unsigned char sb,unsigned char sh)
{ 
// fabrice M. 			
 	DWORD g2;
	long l;
	for (long j=im.Hmin+1;j<im.Hmax;j++)
	{             
	 	l = j*im.Larg+im.Lmin+1;
		for (long i=im.Lmin+1;i<im.Lmax;i++)
		{
			g2 = abs(im.Image[l]-im.Image[l+im.Larg+1]) + abs(im.Image[l+1]-im.Image[l+im.Larg]);
 			if ((g2>=sb)&&(g2<=sh)) im.Image[l]=255;	else	im.Image[l]=0;
 			l++;
		}
   	}
}	    	

 void  Sobel(SImage im,unsigned char  sb,unsigned char  sh)
{ 
// fabrice M.

		HGLOBAL ic;
		unsigned char   * ImageBuf;
		DWORD TailleBuf = im.Larg*im.Haut;
		ic = GlobalAlloc(GMEM_MOVEABLE,TailleBuf);            
		ImageBuf =  (unsigned char   *) GlobalLock(ic);
		long i,j,l,niv; 
		for (j=im.Hmin+1;j<im.Hmax;j++)
		{   
	 		l = j*im.Larg+im.Lmin+1;
			for (i=im.Lmin+1;i<im.Lmax;i++)
			{
  				niv = (long)(labs(2*(long)im.Image[l-1]+(long)im.Image[l+im.Larg-1]+(long)im.Image[l-im.Larg-1]-2*(long)im.Image[l+1]-(long)im.Image[l+im.Larg+1]-(long)im.Image[l-im.Larg+1])+
  					 labs((long)im.Image[l+im.Larg-1]+2*(long)im.Image[l+im.Larg]+(long)im.Image[l+im.Larg+1]-(long)im.Image[l-im.Larg-1]-2*(long)im.Image[l-im.Larg]-(long)im.Image[l-im.Larg+1]));
  				if ((niv>=sb)&&(niv<=sh)) 	ImageBuf[l]  =255;	else ImageBuf[l] =0;
  				l ++;
  			}
	    }
		for ( j=im.Hmin+1;j<im.Hmax;j++)
		{   
			l = j*im.Larg+im.Lmin+1;
			for (i=im.Lmin+1;i<im.Lmax;i++) im.Image[l++] = ImageBuf[l];
	    }
		GlobalUnlock(ic);
  		GlobalFree(ic);
}

  void  Prewitt(SImage im,unsigned char   sb,unsigned char  sh)
{           
		HGLOBAL ic;
		unsigned char   * ImageBuf;
		DWORD T = im.Larg*im.Haut;
		ic = GlobalAlloc(GMEM_MOVEABLE,T);            
		ImageBuf =  (unsigned char   *) GlobalLock(ic);
		long i,j,l,niv;

		for ( j=im.Hmin+1;j<im.Hmax;j++)
		{    
			l = j*im.Larg+im.Lmin+1;
			for ( i=im.Lmin+1;i<im.Lmax;i++) 
			{
  				niv = labs(im.Image[l-1]+im.Image[l+im.Larg-1]+im.Image[l-im.Larg-1]-im.Image[l+1]-im.Image[l+im.Larg+1]-im.Image[l-im.Larg+1])+
  					 labs(im.Image[l+im.Larg-1]+im.Image[l+im.Larg]+im.Image[l+im.Larg+1]-im.Image[l-im.Larg-1]-im.Image[l-im.Larg]-im.Image[l-im.Larg+1]);
  				if ((niv>=sb)&&(niv<=sh)) ImageBuf[l]=255;	else ImageBuf[l]=0;
  				l++;
  			}
	    }
                               
		for (j=im.Hmin+1;j<im.Hmax;j++)
		{   
			l = j*im.Larg+im.Lmin+1;
			for ( i=im.Lmin+1;i<im.Lmax;i++) im.Image[l++] = ImageBuf[l];
	    }

		GlobalUnlock(ic);
  		GlobalFree(ic);
}

