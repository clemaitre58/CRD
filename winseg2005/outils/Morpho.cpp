#include "morpho.h"


unsigned char Mini(unsigned char *b,int n)
{
	unsigned char min=255;
	for (int i=0;i<n;i++) if (b[i]<min) min=b[i];
	return min;
}

unsigned char Maxi(unsigned char *b,int n)
{
	unsigned char max=0;
	for (int i=0;i<n;i++) if (b[i]>max) max=b[i];
	return max;
}


void  Erode(SImage im,int n)
{                                                                                    
	DWORD T = im.Larg*im.Haut;
	HGLOBAL ic;
	ic = GlobalAlloc(GMEM_MOVEABLE,T);
	unsigned char  * tab;
	tab =  (unsigned char  *) GlobalLock(ic) ;  
	unsigned char b[9];
	long last,l,x,y;
	int ind;
	for (ind=0;ind<n;ind++)
	{
   		for (l=0;l<im.Larg*im.Haut;l++) tab[l]=im.Image[l]; 
		for (y=im.Hmin+1;y<im.Hmax;y++)
		{
			for (x=im.Lmin+1;x<im.Lmax;x++)
			{
				last=y*im.Larg+x;
				b[0]=tab[last-im.Larg-1];
				b[1]=tab[last-im.Larg];
				b[2]=tab[last-im.Larg+1];
				b[3]=tab[last-1];
				b[4]=tab[last];
				b[5]=tab[last+1];
				b[6]=tab[last+im.Larg-1];
				b[7]=tab[last+im.Larg];
				b[8]=tab[last+im.Larg+1];
				im.Image[last] = Mini(b,9);
			}
		}
		for (x=im.Lmin+1;x<im.Lmax;x++)
		{
			b[0]=tab[x-1];
			b[1]=tab[x];
			b[2]=tab[x+1];
			b[3]=tab[x-1+im.Larg];
			b[4]=tab[x+im.Larg];
			b[5]=tab[x+1+im.Larg];
			im.Image[x] = Mini(b,6);
		}
		for (x=(im.Lmin+1)+(im.Hmax*im.Larg);x<im.Lmax+(im.Hmax*im.Larg);x++)
		{
			b[0]=tab[x-1-im.Larg];
			b[1]=tab[x-im.Larg];
			b[2]=tab[x+1-im.Larg];
			b[3]=tab[x-1];
			b[4]=tab[x];
			b[5]=tab[x+1];
			im.Image[x] = Mini(b,6);
		}
		for (y=im.Hmin+1;y<im.Hmax;y++)
		{
			b[0]=tab[(y-1)*im.Larg];
			b[1]=tab[(y-1)*im.Larg+1];
			b[2]=tab[y*im.Larg];
			b[3]=tab[y*im.Larg+1];
			b[4]=tab[(y+1)*im.Larg];
			b[5]=tab[(y+1)*im.Larg+1];
			im.Image[x] = Mini(b,6);
		}
		last=im.Lmax-1;
		for (y=im.Hmin+1;y<im.Hmax;y++)
		{
			b[0]=tab[(y-1)*im.Larg+last];
			b[1]=tab[(y-1)*im.Larg+1+last];
			b[2]=tab[y*im.Larg+last];
			b[3]=tab[y*im.Larg+1+last];
			b[4]=tab[(y+1)*im.Larg+last];
			b[5]=tab[(y+1)*im.Larg+1+last];
			im.Image[x] = Mini(b,6);
		}
		{
			b[0]=tab[0];
			b[1]=tab[1];
			b[2]=tab[im.Larg];
			b[3]=tab[im.Larg+1];
			im.Image[0]= Mini(b,4);
		}
		{
			b[0]=tab[im.Lmax-1];
			b[1]=tab[im.Lmax];
			b[2]=tab[im.Lmax+im.Larg-1];
			b[3]=tab[im.Lmax+im.Larg];
			im.Image[im.Lmax]= Mini(b,4);
		}
		{
			b[0]=tab[im.Hmax*(im.Larg-1)];
			b[1]=tab[im.Hmax*(im.Larg-1)+1];
			b[2]=tab[im.Hmax*im.Larg];
			b[3]=tab[im.Hmax*im.Larg+1];
			im.Image[im.Hmax*im.Larg]= Mini(b,4);
		}
		{
			b[0]=tab[im.Hmax*(im.Larg-1)+im.Lmax-1];
			b[1]=tab[im.Hmax*(im.Larg-1)+im.Lmax];
			b[2]=tab[im.Hmax*im.Larg+im.Lmax-1];
			b[3]=tab[im.Hmax*im.Larg+im.Lmax];
			im.Image[im.Hmax*im.Larg+im.Lmax]= Mini(b,4);
		}
	}
	GlobalUnlock(ic);
	GlobalFree(ic); 	
}


void  Dilate(SImage im,int n)
{                                                                                    
	DWORD T = im.Larg*im.Haut;
	HGLOBAL ic;
	ic = GlobalAlloc(GMEM_MOVEABLE,T);
	unsigned char  * tab;
	tab =  (unsigned char  *) GlobalLock(ic) ;  
	unsigned char b[9];
	long last,l,x,y;
	int ind;
	for (ind=0;ind<n;ind++)
	{
   		for (l=0;l<im.Larg*im.Haut;l++) tab[l]=im.Image[l]; 
		for (y=im.Hmin+1;y<im.Hmax;y++)
		{
			for (x=im.Lmin+1;x<im.Lmax;x++)
			{
				last=y*im.Larg+x;
				b[0]=tab[last-im.Larg-1];
				b[1]=tab[last-im.Larg];
				b[2]=tab[last-im.Larg+1];
				b[3]=tab[last-1];
				b[4]=tab[last];
				b[5]=tab[last+1];
				b[6]=tab[last+im.Larg-1];
				b[7]=tab[last+im.Larg];
				b[8]=tab[last+im.Larg+1];
				im.Image[last] = Maxi(b,9);
			}
		}
		for (x=im.Lmin+1;x<im.Lmax;x++)
		{
			b[0]=tab[x-1];
			b[1]=tab[x];
			b[2]=tab[x+1];
			b[3]=tab[x-1+im.Larg];
			b[4]=tab[x+im.Larg];
			b[5]=tab[x+1+im.Larg];
			im.Image[x] = Maxi(b,6);
		}
		for (x=(im.Lmin+1)+(im.Hmax*im.Larg);x<im.Lmax+(im.Hmax*im.Larg);x++)
		{
			b[0]=tab[x-1-im.Larg];
			b[1]=tab[x-im.Larg];
			b[2]=tab[x+1-im.Larg];
			b[3]=tab[x-1];
			b[4]=tab[x];
			b[5]=tab[x+1];
			im.Image[x] = Maxi(b,6);
		}
		for (y=im.Hmin+1;y<im.Hmax;y++)
		{
			b[0]=tab[(y-1)*im.Larg];
			b[1]=tab[(y-1)*im.Larg+1];
			b[2]=tab[y*im.Larg];
			b[3]=tab[y*im.Larg+1];
			b[4]=tab[(y+1)*im.Larg];
			b[5]=tab[(y+1)*im.Larg+1];
			im.Image[x] = Maxi(b,6);
		}
		last=im.Lmax-1;
		for (y=im.Hmin+1;y<im.Hmax;y++)
		{
			b[0]=tab[(y-1)*im.Larg+last];
			b[1]=tab[(y-1)*im.Larg+1+last];
			b[2]=tab[y*im.Larg+last];
			b[3]=tab[y*im.Larg+1+last];
			b[4]=tab[(y+1)*im.Larg+last];
			b[5]=tab[(y+1)*im.Larg+1+last];
			im.Image[x] = Maxi(b,6);
		}
		{
			b[0]=tab[0];
			b[1]=tab[1];
			b[2]=tab[im.Larg];
			b[3]=tab[im.Larg+1];
			im.Image[0]= Maxi(b,4);
		}
		{
			b[0]=tab[im.Lmax-1];
			b[1]=tab[im.Lmax];
			b[2]=tab[im.Lmax+im.Larg-1];
			b[3]=tab[im.Lmax+im.Larg];
			im.Image[im.Lmax]= Maxi(b,4);
		}
		{
			b[0]=tab[im.Hmax*(im.Larg-1)];
			b[1]=tab[im.Hmax*(im.Larg-1)+1];
			b[2]=tab[im.Hmax*im.Larg];
			b[3]=tab[im.Hmax*im.Larg+1];
			im.Image[im.Hmax*im.Larg]= Maxi(b,4);
		}
		{
			b[0]=tab[im.Hmax*(im.Larg-1)+im.Lmax-1];
			b[1]=tab[im.Hmax*(im.Larg-1)+im.Lmax];
			b[2]=tab[im.Hmax*im.Larg+im.Lmax-1];
			b[3]=tab[im.Hmax*im.Larg+im.Lmax];
			im.Image[im.Hmax*im.Larg+im.Lmax]= Maxi(b,4);
		}
	}
	GlobalUnlock(ic);
	GlobalFree(ic); 	
}
