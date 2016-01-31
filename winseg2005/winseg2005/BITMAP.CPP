#include "stdafx.h"
#include "bitmap.h"

	
void SetInfoBitmap(BITMAPINFO * info,long L,long H,NBCOLOR NBRGB)
{   
	info->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	info->bmiHeader.biWidth = L;
	info->bmiHeader.biHeight = H;
	info->bmiHeader.biPlanes = 1;   
	if (NBRGB == NB) 
	{  
	    info->bmiHeader.biSizeImage = L*H;
		info->bmiHeader.biBitCount = 8;  
	    info->bmiHeader.biClrUsed = 256;   
	    for (int i = 0 ; i<256; i ++)
	    {
    		info->bmiColors[i].rgbBlue = i;
    		info->bmiColors[i].rgbGreen = i;
    		info->bmiColors[i].rgbRed = i;
    		info->bmiColors[i].rgbReserved = 0;
    	}	
	}
	else
	{
	    info->bmiHeader.biSizeImage = sizeof(rgb)*L*H;
	  	info->bmiHeader.biBitCount = 24;
	    info->bmiHeader.biClrUsed = 0;
	}
	
    info->bmiHeader.biCompression =  BI_RGB;
    info->bmiHeader.biXPelsPerMeter = 0;
    info->bmiHeader.biYPelsPerMeter = 0;
    info->bmiHeader.biClrImportant = 0;
}

void SetInfoBitmap(BITMAPINFO * infoDest,long L,long H,BITMAPINFO * infoSource)
{                    

	infoDest->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	infoDest->bmiHeader.biWidth = L;
	infoDest->bmiHeader.biHeight = H;
	infoDest->bmiHeader.biPlanes = 1;   
	infoDest->bmiHeader.biBitCount = 8;  
    infoDest->bmiHeader.biClrUsed = 256;   
    for (int i = 0 ; i<256; i ++)
    {
   		infoDest->bmiColors[i].rgbBlue = infoSource->bmiColors[i].rgbBlue;
   		infoDest->bmiColors[i].rgbGreen = infoSource->bmiColors[i].rgbGreen;
   		infoDest->bmiColors[i].rgbRed = infoSource->bmiColors[i].rgbRed ;
   		infoDest->bmiColors[i].rgbReserved = infoSource->bmiColors[i].rgbReserved;
   	}	
	
    infoDest->bmiHeader.biCompression =  BI_RGB;
    infoDest->bmiHeader.biSizeImage = L*H;
    infoDest->bmiHeader.biXPelsPerMeter = 0;
    infoDest->bmiHeader.biYPelsPerMeter = 0;
    infoDest->bmiHeader.biClrImportant = 0;
}

void SetTeteBitmap(BITMAPFILEHEADER * tete,long L,long H,NBCOLOR NBRGB)
{
	if (NBRGB == COLOR) 
	{
		tete->bfOffBits = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER); 
		tete->bfSize = tete->bfOffBits+L*H*sizeof(rgb);
	}	
	else 
	{
		tete->bfOffBits = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER)+256*sizeof(RGBQUAD); 
		tete->bfSize = tete->bfOffBits+L*H;
	}	
	tete->bfType = 'MB';
	tete->bfReserved1 = 0;
	tete->bfReserved2 = 0;  
	tete->bfOffBits;
}

