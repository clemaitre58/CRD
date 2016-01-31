#include "fourier.h"
#include <math.h>
#include <stdlib.h>       
//#include <windows.h>
#include <stdio.h>

float Pi = 3.141592654F; 


 long  puissance2(long p)
{
  	return (long) pow((double)2,(double)p);
}

 long  TestPuissance2(long valeur)
{
 long v=1,i=0;
 
 do 
 {  v = v*2;   
   	i++;
 }
 while (v<=valeur); 
 return i-1;
}


void FFT1D(long q, SComplexe  *TFFT1D_S)
{
    // bloc FFT 
    
    long Npoint = puissance2(q); 
   	DWORD Taille = Npoint;
	HGLOBAL Hdle = GlobalAlloc(GMEM_MOVEABLE,Taille*sizeof(SComplexe));            
	SComplexe  * TFFT1D_D =  (SComplexe  *) GlobalLock(Hdle);    

    long i,j;
    long nx,ni;
   	div_t r;
	for (i=0;i<Npoint;i++) 
	{                 
		nx = 0;
		ni = i;    
		for (j=1;j<=q;j++) 
		{   
			r = div((int)ni,2);
			nx = nx*2 + r.rem;
			ni = r.quot;
		}   
		TFFT1D_D[nx] =  TFFT1D_S[i];           
	}
	                               
	for (i=0;i<Npoint;i++) TFFT1D_S[i] = TFFT1D_D[i];
	
    long jn,nn = 1,p,ind,indice;     
    float zreel,zimaj,wreel,wimaj,theta;
    for (long l=1; l<=q;l++)
    {                       
    	jn = (Npoint/2)/nn -1;
	    for (j=0;j<=jn;j++)
	    {
	    	p = 2*nn*j;
	    	for (i=0;i<nn;i++)
	    	{             
				ind = p + i;
				indice = ind + nn;
				theta = Pi*(float)i/(float)nn;
				wreel = (float)cos((double)theta);
				wimaj = -(float)sin((double)theta);
				zreel = wreel*TFFT1D_S[indice].Re - wimaj*TFFT1D_S[indice].Im;
				zimaj = wimaj*TFFT1D_S[indice].Re + wreel*TFFT1D_S[indice].Im;   
				
				TFFT1D_S[indice].Re = TFFT1D_S[ind].Re - zreel;
				TFFT1D_S[indice].Im = TFFT1D_S[ind].Im - zimaj;
				
				TFFT1D_S[ind].Re =	TFFT1D_S[ind].Re + zreel;
				TFFT1D_S[ind].Im = 	TFFT1D_S[ind].Im + zimaj;
	    	}
	    }  
	    nn = nn*2;
    }
    
    // delete de l'espace mémoire
  	GlobalUnlock(Hdle);  
  	GlobalFree(Hdle);     
}

       
 void  FFT2D(SImage imS,SImage imModule,SImage imPhase)
{
    // Largeur et hauteur de l'image en puissance de 2 !!  
    char c[50];      
	long i,j,l,k;
    long Larg = imS.Lmax-imS.Lmin+1;
    long Haut = imS.Hmax-imS.Hmin+1;
     
    long qL = TestPuissance2(Larg);
    long qH = TestPuissance2(Haut);
    
    long L = puissance2(qL);
    long H = puissance2(qH);    
    
   if ((L>0) && (H>0))
   {
	imModule.Larg = L;
	imModule.Haut = H;
    imPhase.Larg = L;
    imPhase.Haut = H;
	
    // bloc memoire d'une image TFFT 
	DWORD TailleImage = L*H;
	HGLOBAL HdleImage = GlobalAlloc(GMEM_MOVEABLE,TailleImage*sizeof(SComplexe));            
	SComplexe  * TFFT2D =  (SComplexe  *) GlobalLock(HdleImage);    
    
    // bloc Ligne FFT
	DWORD TailleLigne = L;
	HGLOBAL HdleLigne = GlobalAlloc(GMEM_MOVEABLE,TailleLigne*sizeof(SComplexe));            
	SComplexe  * TFFT1DLig =  (SComplexe  *) GlobalLock(HdleLigne);    
    
    // bloc Colonne FFT 
	DWORD TailleColonne = H;
	HGLOBAL HdleColonne = GlobalAlloc(GMEM_MOVEABLE,TailleColonne*sizeof(SComplexe));            
	SComplexe  * TFFT1DCol =  (SComplexe  *) GlobalLock(HdleColonne);    
 
	// parcourt dans le sens horizontal    
    for (j = 0;j<H;j++)
    {
	    k = (j+imS.Hmin)*imS.Larg+imS.Lmin;
    	for (i=0;i<L;i++)
		{   
        	TFFT1DLig[i].Re = imS.Image[k++];                  // 1
    		TFFT1DLig[i].Im = 0.F;	
    	}   
    	FFT1D(qL,TFFT1DLig); 
    	k = j*L;
    	for (i=0;i<L;i++) TFFT2D[k++] = TFFT1DLig[i];
    }			
      
	// parcourt dans le sens vertical  
	long jj = H/2;
    float Rmin,Rmax=-1e12F;
    float R1;        
    float PhaseMin = 1e12F, PhaseMax = -1e12F,Phase1,Phase2;
    
    for (i = 0;i<L;i++)
    {
	    l = 0;
    	for (j=0;j<H;j++) TFFT1DCol[l++] = TFFT2D[j*L+i];      // 2
    	FFT1D(qH,TFFT1DCol);     
		for (j=0;j<jj;j++)    	
		{   
			R1 = (float) sqrt((double) (TFFT1DCol[j].Re*TFFT1DCol[j].Re + TFFT1DCol[j].Im*TFFT1DCol[j].Im));  
			if (R1 > Rmax) Rmax = R1; 

			if (TFFT1DCol[j].Re==0.0)  
			{                                          
				if (TFFT1DCol[j].Im==0.F) Phase1 = 0.F;
				else if (TFFT1DCol[j].Im>0) Phase1 = Pi/(float)2.0; 
				else Phase1 = -Pi/(float)2.0; 
			}			
			else
			{
				Phase1 = (float) atan((double)(TFFT1DCol[j].Im/TFFT1DCol[j].Re));
				if (TFFT1DCol[j].Re<0.0)
				{	if (Phase1>0) Phase1+=Pi/(float)2.0; 
					else Phase1-=Pi/(float)2.0;
				}
			}	
			k = j*L+i;
			TFFT2D[k].Re = (float) sqrt((double)(TFFT1DCol[j+jj].Re*TFFT1DCol[j+jj].Re  + TFFT1DCol[j+jj].Im*TFFT1DCol[j+jj].Im)); 
			TFFT2D[k+jj*L].Re = R1; 
			
			
			if (TFFT1DCol[j+jj].Re==0.0)  
			{                                          
				if (TFFT1DCol[j+jj].Im==0.0) Phase2 = 0.F;
				else if (TFFT1DCol[j+jj].Im>0) Phase2 = Pi/(float)2.0; 
				else Phase2 = -Pi/(float)2.0; 
			}			
			else
			{
				Phase2 = (float) atan((double)(TFFT1DCol[j+jj].Im/TFFT1DCol[j+jj].Re));  
				if (TFFT1DCol[j+jj].Re<0.0)
				{	if (Phase2>0.0) Phase2+= Pi/(float)2.0; 
					else Phase2-=Pi/(float)2.0;
				}
			}	

			TFFT2D[k].Im 		= Phase2;
			TFFT2D[k+jj*L].Im 	= Phase1;
			
		}
	}			 
    
    Rmin = 0.F;
 
	// Normalisation
	float normph,norm;                           
	Rmax = (float) log((double)(Rmax+1.0)); 
    
    if (Rmax!=0) norm = (float)255/Rmax; else norm = 1.F;                         

	PhaseMin = -Pi;
	PhaseMax = Pi;
    normph = (float)log((double)(PhaseMax-PhaseMin));
    long ii = L/2;
    for (j=0;j<H;j++)
    {   
    	l = j*L;
    	for (i=0;i<ii;i++)
    	{   
 			imModule.Image[l+i+ii] 	= (unsigned char) ((log((double)(1.0+TFFT2D[l+i].Re)))*norm);
 			imModule.Image[l+i] 	= (unsigned char) ((log((double)(1.0+TFFT2D[l+i+ii].Re)))*norm);  
 			imPhase.Image[l+i+ii] 	= (unsigned char) (255.0*(float)log((double)(TFFT2D[l+i].Im-PhaseMin+1.0))/	normph);
 			imPhase.Image[l+i] 		= (unsigned char) (255.0*(float)log((double)(TFFT2D[l+i+ii].Im-PhaseMin+1.0)) /	normph); 
    	}
    }
    
    
    // delete des blocs mémoires Image et Ligne           // 4
  	GlobalUnlock(HdleImage); 	GlobalFree(HdleImage);
  	GlobalUnlock(HdleLigne);  	GlobalFree(HdleLigne);    
  	GlobalUnlock(HdleColonne);  GlobalFree(HdleColonne); 
           
	sprintf(c,"min = %e max = %e",Rmin,Rmax);                 
  	MessageBox(NULL,c,"information",MB_OK); 
  }	
  else MessageBox(NULL,"dimension de la fenêtre trop petite","information",MB_OK);
}
       
