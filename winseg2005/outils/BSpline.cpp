#include "bspline.h"


void OnConvBspline(SImageFloat im,SImageFloat dest)
{   
	// Conversion d'une image float en image SPLINE
   float cp,al;
   long j,i,l,lar;  
   
   al = float(2-sqrt(3)); 
   // traitement sur l'image de dest 
   lar=im.Lmax-im.Lmin+1;
   l=0;
   j=im.Hmin;i=im.Lmin;
   cp=im.Image[j*im.Larg+i];
   dest.Image[l++]=cp;
   for (i=im.Lmin+1;i<=im.Lmax;i++)
   {
   		cp=im.Image[j*im.Larg+i]+al*(im.Image[j*im.Larg+i]-cp);
   		dest.Image[l++]=cp;
   }
   for (j=im.Hmin+1;j<=im.Hmax;j++) 
   {	                   
   		i=im.Lmin;
   		cp=im.Image[j*im.Larg+i]; 
   		dest.Image[l]=cp+al*(cp-dest.Image[l-lar]);  
   		l++;
   		for (i=im.Lmin+1;i<=im.Lmax;i++) 
   		{	            
   			cp=im.Image[j*im.Larg+i]+al*(im.Image[j*im.Larg+i]-cp);
   			dest.Image[l]=cp+al*(cp-dest.Image[l-lar]);  
   			l++;
   		}	
   }	    
   l--;
   j=im.Hmax;i=im.Lmax;
   cp=dest.Image[l--];
   for (i=im.Lmax-1;i>=im.Lmin;i--)
   {
   		cp=dest.Image[l]+al*(dest.Image[l]-cp);
   		dest.Image[l--]=cp;
   }
   for (j=im.Hmax-1;j>=im.Hmin;j--) 
   {	                   
   		i=im.Lmax;
   		cp=dest.Image[l]; 
   		dest.Image[l]=cp+al*(cp-dest.Image[l+lar]);  
   		l--;
   		for (i=im.Lmax-1;i>=im.Lmin;i--) 
   		{	            
   			cp=dest.Image[l]+al*(dest.Image[l]-cp);
   			dest.Image[l]=cp+al*(cp-dest.Image[l+lar]);  
   			l--;
   		}	
   }
   // fin du traitement
}




void OnConvBsplineFloat(SImage im,SImageFloat dest)
{  
   	//  Conversion d'une image byte en image Spline

   float cp,al;
   long j,i,l,lar;  
   
   al = float(2-sqrt(3)); 
   // traitement sur l'image de dest 
   lar=im.Lmax-im.Lmin+1;
   l=0;
   j=im.Hmin;i=im.Lmin;
   cp=im.Image[j*im.Larg+i];
   dest.Image[l++]=cp;
   for (i=im.Lmin+1;i<=im.Lmax;i++)
   {
   		cp=im.Image[j*im.Larg+i]+al*(im.Image[j*im.Larg+i]-cp);
   		dest.Image[l++]=cp;
   }
   for (j=im.Hmin+1;j<=im.Hmax;j++) 
   {	                   
   		i=im.Lmin;
   		cp=im.Image[j*im.Larg+i]; 
   		dest.Image[l]=cp+al*(cp-dest.Image[l-lar]);  
   		l++;
   		for (i=im.Lmin+1;i<=im.Lmax;i++) 
   		{	            
   			cp=im.Image[j*im.Larg+i]+al*(im.Image[j*im.Larg+i]-cp);
   			dest.Image[l]=cp+al*(cp-dest.Image[l-lar]);  
   			l++;
   		}	
   }	    
   l--;
   j=im.Hmax;i=im.Lmax;
   cp=dest.Image[l--];
   for (i=im.Lmax-1;i>=im.Lmin;i--)
   {
   		cp=dest.Image[l]+al*(dest.Image[l]-cp);
   		dest.Image[l--]=cp;
   }
   for (j=im.Hmax-1;j>=im.Hmin;j--) 
   {	                   
   		i=im.Lmax;
   		cp=dest.Image[l]; 
   		dest.Image[l]=cp+al*(cp-dest.Image[l+lar]);  
   		l--;
   		for (i=im.Lmax-1;i>=im.Lmin;i--) 
   		{	            
   			cp=dest.Image[l]+al*(dest.Image[l]-cp);
   			dest.Image[l]=cp+al*(cp-dest.Image[l+lar]);  
   			l--;
   		}	
   }
   // fin du traitement
}        





void OnGradientBspline(SImageFloat im,SImageFloat dest)
{   
	// 	Calcul du gradient d'une image a partir de son image Spline
   	long j,i,l,ll;  
    float p,q;
    
   	// traitement sur l'image de dest 
   	ll=0;
   	for (i=im.Lmin;i<=im.Lmax;i++)
   	{
   		dest.Image[ll++]=0;
   	}
   	for (j=im.Hmin+1;j<im.Hmax;j++) 
   	{	                   
   		dest.Image[ll++]=0;  
   		for (i=im.Lmin+1;i<im.Lmax;i++) 
   		{	                     
   			l=j*im.Larg+i;
   			p=((im.Image[l+1+im.Larg]+4*im.Image[l+1]+im.Image[l+1-im.Larg])
   				-(im.Image[l-1+im.Larg]+4*im.Image[l-1]+im.Image[l-1-im.Larg]))/12;
   			q=((im.Image[l-1+im.Larg]+4*im.Image[l+im.Larg]+im.Image[l+1+im.Larg])
   				-(im.Image[l-1-im.Larg]+4*im.Image[l-im.Larg]+im.Image[l+1-im.Larg]))/12;
   			dest.Image[ll++]=float(sqrt(p*p+q*q));  
   		}
   		dest.Image[ll++]=0;	
   	}
   	for (i=im.Lmin;i<=im.Lmax;i++)
   	{
   		dest.Image[ll++]=0;
   	}	    
   	// fin du traitement
}



void OnDerDirSec(SImageFloat im,SImageFloat dest)
{                                                     
    
    //Calcul de la derivee directionnelle seconde d'une image
    //dans la direction du gradient et a partir de son image Spline 
   	long j,i,l,ll;  
    float p,q,r,s,t,result1,result2,result3;
    
   	// traitement sur l'image de dest 
   	ll=0;
   	for (i=im.Lmin;i<=im.Lmax;i++)
   	{
   		dest.Image[ll++]=0;
   	}
   	for (j=im.Hmin+1;j<im.Hmax;j++) 
   	{	                   
   		dest.Image[ll++]=0;  
   		for (i=im.Lmin+1;i<im.Lmax;i++) 
   		{	                         
   			l=j*im.Larg+i;
   			p=((im.Image[l+1+im.Larg]+4*im.Image[l+1]+im.Image[l+1-im.Larg])
   				-(im.Image[l-1+im.Larg]+4*im.Image[l-1]+im.Image[l-1-im.Larg]))/12;
   			q=((im.Image[l-1+im.Larg]+4*im.Image[l+im.Larg]+im.Image[l+1+im.Larg])
   				-(im.Image[l-1-im.Larg]+4*im.Image[l-im.Larg]+im.Image[l+1-im.Larg]))/12;   
   			r=((im.Image[l-1+im.Larg]+4*im.Image[l-1]+im.Image[l-1-im.Larg])
   				-2*(im.Image[l+im.Larg]+4*im.Image[l]+im.Image[l-im.Larg])
   				+(im.Image[l+1+im.Larg]+4*im.Image[l+1]+im.Image[l+1-im.Larg]))/6;   
   			s=((-im.Image[l-1+im.Larg]+im.Image[l+1+im.Larg])
   				+(im.Image[l-1-im.Larg]-im.Image[l+1-im.Larg]))/4;
   			t=((im.Image[l-1+im.Larg]+4*im.Image[l+im.Larg]+im.Image[l+1+im.Larg])
   				-2*(im.Image[l-1]+4*im.Image[l]+im.Image[l+1])
   				+(im.Image[l-1-im.Larg]+4*im.Image[l-im.Larg]+im.Image[l+1-im.Larg]))/6;

			result1=r*p*p;
			result2=2*s*p*q;
			result3=t*q*q;
			dest.Image[ll++]=result1+result2+result3;
   		}
   		dest.Image[ll++]=0;	
   	}
   	for (i=im.Lmin;i<=im.Lmax;i++)
   	{
   		dest.Image[ll++]=0;
   	}	    
   	// fin du traitement
}





void OnLigneZero(SImageFloat im,SImage dest)
{
	//Recherche du passage par zero d'une image float 
   	long j,i,l,ll;
   	double p,q,moy;  
    int s;
   
   	// traitement sur l'image de dest
   	ll=0;q=0;
   	for (j=im.Hmin;j<=im.Hmax;j++) 
   	{	                   
   		for (i=im.Lmin;i<=im.Lmax;i++) 
   		{	 
			 q=q+fabs(im.Image[j*im.Larg+i]);
			 ll++;  			      
   		}
   	}                   
   	moy=8*q/ll;
   	ll=0;
   	for (j=im.Hmin;j<im.Hmax;j++) 
   	{	                   
   		for (i=im.Lmin;i<im.Lmax;i++) 
   		{	                         
   			l=j*im.Larg+i;
   			s=0;q=0;
   			p=im.Image[l+1+im.Larg];
   			q=q+fabs(p); 
   			if(p<0)s++;
   			p=im.Image[l+im.Larg];
   			q=q+fabs(p); 
   			if(p<0)s++;
   			p=im.Image[l+1];
   			q=q+fabs(p); 
   			if(p<0)s++;
   			p=im.Image[l];
   			q=q+fabs(p); 
   			if(p<0)s++;
  			if((s>0) && (s<4) && (q>moy))
   				dest.Image[ll++]=255; 
   				else  dest.Image[ll++]=0;
   		}
   		i=im.Lmax;
   		l=j*im.Larg+i;
   		s=0;q=0;
   		p=im.Image[l+im.Larg];
   		q=q+fabs(p); 
   		if(p<0)s++;
   		p=im.Image[l];
   		q=q+fabs(p); 
   		if(p<0)s++;
  		if((s>0) && (s<2) && (q>(moy/2)))
   			dest.Image[ll++]=255; 
   			else  dest.Image[ll++]=0;
   	}   
   	j=im.Hmax;
   	for (i=im.Lmin;i<im.Lmax;i++)
   	{
   		l=j*im.Larg+i;
   		s=0;q=0;
   		p=im.Image[l+1];
   		q=q+fabs(p); 
   		if(p<0)s++;
   		p=im.Image[l];
   		q=q+fabs(p); 
   		if(p<0)s++;
  		if((s>0) && (s<2) && (q>(moy/2)))
   			dest.Image[ll++]=255; 
   			else  dest.Image[ll++]=0;
   	}	                
   	dest.Image[ll++]=0;
   	// fin du traitement
}



////////////////////////////////////////
void OnDecompo(SImageFloat im,SImageFloat dest,int prof)
{
	// Decomposition en ondelettes d'une image float 	
   float cp,ca;
   long j,i,l,lar,larg,hau; 
   float * buffer;
   float * buf2;
   float al = 0.535280431F; 
   float be = 0.1225546152F;
   float ga = 0.00914869485F;
   int pr;   
   
   // traitement sur l'image de dest 
   larg=im.Lmax-im.Lmin+1;    
   lar = larg;
   hau=im.Hmax-im.Hmin+1;
	l=0;                      
   	for (j=im.Hmin;j<=im.Hmax;j++)
   	{            
    	for (i=im.Lmin;i<=im.Lmax;i++)
   		{
   			dest.Image[l++]=im.Image[j*im.Larg+i];
   		}
   	}      
	for (pr=0;pr<prof;pr++)
	
	{	
	// Traitement horizontal
	// traitement image * h    	
   	buffer 	= new float[lar];
   	buf2	= new float[lar];
   	for (j=0;j<hau;j++)
   	{
   	    for (i=0;i<lar;i++)
   		{
   			buf2[i]= dest.Image[j*larg+i];
   			buffer[i]=  buf2[i];
   		}
   		cp = buffer[0];
   		for (i=0;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=(ca+cp)/2;
   			cp=ca;
   		}  
   		cp = buffer[lar-1];
   		for (i=lar-1;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=(ca+cp)/2;
   			cp=ca;
   		}
    	cp = buffer[0];
   		for (i=0;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=(ca+cp)/2;
   			cp=ca;
   		}  
   		cp = buffer[lar-1];
   		for (i=lar-1;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=(ca+cp)/2;
   			cp=ca;
   		}                      
   		for (i=0;i<(lar/2);i++)
   		{
   			dest.Image[j*larg+i] =buffer[2*i];
   		}
  		// traitement image * g
    	for (i=0;i<lar;i++)
   		{
   			buffer[i]=buf2[i];
   		}
   		cp = buffer[0];
   		buffer[0]=0;
   		for (i=1;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=(-ca+cp)/2;
   			cp=ca;
   		}  
   		cp = buffer[lar-1];
   		buffer[lar-1]=-cp/2;
   		for (i=lar-2;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=(-ca+cp)/2;
   			cp=ca;
   		}
    	cp = buffer[0];
    	buffer[0]=-cp/2;
   		for (i=1;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=(-ca+cp)/2;
   			cp=ca;
   		}  
   		cp = buffer[lar-1];
   		buffer[lar-1]=-cp/2;
   		for (i=lar-2;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=(-ca+cp)/2;
   			cp=ca;
   		}        
   		cp = buffer[0];
   		buffer[0]=-cp;       
   		ca = buffer[1];
   		buffer[1]=-ca+al*cp;
   		cp=ca;
   		for (i=2;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=-ca+al*(cp-buffer[i-2]);
   			cp=ca;
   		}  
   		cp = buffer[lar-1];
   		buffer[lar-1]=-cp;
   		ca=buffer[lar-2];
   		buffer[lar-2]=-ca+al*cp;
   		cp=ca;
   		for (i=lar-3;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=-ca+al*(cp-buffer[i+2]);
   			cp=ca;
   		}   
   		cp = buffer[0];
   		buffer[0]=-cp;       
   		ca = buffer[1];
   		buffer[1]=-ca+be*cp;
   		cp=ca;
   		for (i=2;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=-ca+be*(cp-buffer[i-2]);
   			cp=ca;
   		}  
   		cp = buffer[lar-1];
   		buffer[lar-1]=-cp;
   		ca=buffer[lar-2];
   		buffer[lar-2]=-ca+be*cp;
   		cp=ca;
   		for (i=lar-3;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=-ca+be*(cp-buffer[i+2]);
   			cp=ca;
   		}                   
   		cp = buffer[0];
   		buffer[0]=-cp;       
   		ca = buffer[1];
   		buffer[1]=-ca+ga*cp;
   		cp=ca;
   		for (i=2;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=-ca+ga*(cp-buffer[i-2]);
   			cp=ca;
   		}  
   		cp = buffer[lar-1];
   		buffer[lar-1]=-cp;
   		ca=buffer[lar-2];
   		buffer[lar-2]=-ca+ga*cp;
   		cp=ca;
   		for (i=lar-3;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=-ca+ga*(cp-buffer[i+2]);
   			cp=ca;
   		} 
   		for (i=0;i<(lar/2);i++)
   		{
   			dest.Image[j*larg+i+lar/2] =buffer[2*i+1];
   		}
    }
  	delete buffer;
  	delete buf2;         

	// traitement vertical
   buffer = new float[hau];                     
   buf2 =	new float[hau];
   	for (i=0;i<lar;i++)
   	{            
   	   	for (j=0;j<hau;j++)
   		{
   			buf2[j]= dest.Image[j*larg+i];
   		}
   		
		// traitement image * h 
   	    for (j=0;j<hau;j++)
   		{
   			buffer[j]=  buf2[j];
   		}
   		cp = buffer[0];
   		for (j=0;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=(ca+cp)/2;
   			cp=ca;
   		}  
   		cp = buffer[hau-1];
   		for (j=hau-1;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=(ca+cp)/2;
   			cp=ca;
   		}
    	cp = buffer[0];
   		for (j=0;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=(ca+cp)/2;
   			cp=ca;
   		}  
   		cp = buffer[hau-1];
   		for (j=hau-1;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=(ca+cp)/2;
   			cp=ca;
   		}                      
   		for (j=0;j<(hau/2);j++)
   		{
   			dest.Image[j*larg+i] =buffer[2*j];
   		}
    
  		// traitement image * g 
    	for (j=0;j<hau;j++)
   		{
   			buffer[j]=buf2[j];
   		}
   		cp = buffer[0];
   		buffer[0]=0;
   		for (j=1;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=(-ca+cp)/2;
   			cp=ca;
   		}  
   		cp = buffer[hau-1];
   		buffer[hau-1]=-cp/2;
   		for (j=hau-2;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=(-ca+cp)/2;
   			cp=ca;
   		}
    	cp = buffer[0];
    	buffer[0]=-cp/2;
   		for (j=1;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=(-ca+cp)/2;
   			cp=ca;
   		}  
   		cp = buffer[hau-1];
   		buffer[hau-1]=-cp/2;
   		for (j=hau-2;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=(-ca+cp)/2;
   			cp=ca;
   		}        
   		cp = buffer[0];
   		buffer[0]=-cp;       
   		ca = buffer[1];
   		buffer[1]=-ca+al*cp;
   		cp=ca;
   		for (j=2;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=-ca+al*(cp-buffer[j-2]);
   			cp=ca;
   		}  
   		cp = buffer[hau-1];
   		buffer[hau-1]=-cp;
   		ca=buffer[hau-2];
   		buffer[hau-2]=-ca+al*cp;
   		cp=ca;
   		for (j=hau-3;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=-ca+al*(cp-buffer[j+2]);
   			cp=ca;
   		}   
   		cp = buffer[0];
   		buffer[0]=-cp;       
   		ca = buffer[1];
   		buffer[1]=-ca+be*cp;
   		cp=ca;
   		for (j=2;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=-ca+be*(cp-buffer[j-2]);
   			cp=ca;
   		}  
   		cp = buffer[hau-1];
   		buffer[hau-1]=-cp;
   		ca=buffer[hau-2];
   		buffer[hau-2]=-ca+be*cp;
   		cp=ca;
   		for (j=hau-3;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=-ca+be*(cp-buffer[j+2]);
   			cp=ca;
   		}                   
   		cp = buffer[0];
   		buffer[0]=-cp;       
   		ca = buffer[1];
   		buffer[1]=-ca+ga*cp;
   		cp=ca;
   		for (j=2;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=-ca+ga*(cp-buffer[j-2]);
   			cp=ca;
   		}  
   		cp = buffer[hau-1];
   		buffer[hau-1]=-cp;
   		ca=buffer[hau-2];
   		buffer[hau-2]=-ca+ga*cp;
   		cp=ca;
   		for (j=hau-3;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=-ca+ga*(cp-buffer[j+2]);
   			cp=ca;
   		} 
   		for (j=0;j<(hau/2);j++)
   		{
   			dest.Image[(j+hau/2)*larg+i] =buffer[2*j+1];
   		}
    }
   	delete buffer; 
   	delete buf2;
    hau=hau/2;
    lar=lar/2;
    }  
   // fin du traitement
}
	



void OnRecomp(SImageFloat im,SImageFloat dest, int prof)
{
	// Reconstruction en image float d'une image ondelette 	
   float cp,ca;
   long j,i,l,lar,larg,hau;  
   float * buffer;
   float * buf2;
   float al = 0.535280431F; 
   float be = 0.1225546152F;
   float ga = 0.00914869485F;
   int pr;
   
   	// traitement sur l'image de dest 
   	larg=im.Lmax-im.Lmin+1;    
   	lar = larg;
   	hau=im.Hmax-im.Hmin+1;
	l=0;                      
   	for (j=im.Hmin;j<=im.Hmax;j++)
   	{            
    	for (i=im.Lmin;i<=im.Lmax;i++)
   		{
   			dest.Image[l++]=im.Image[j*im.Larg+i];
   		}
   	}
   	for (pr=0; pr<prof-1;pr++)
	{
		hau=hau/2;
		lar=lar/2;		
	}
    for(pr=0;pr<prof;pr++)
    {
    
  	// traitement vertical
   	buffer = new float[hau];                     
   	buf2 =	new float[hau];
   	for (i=0;i<lar;i++)
   	{   
   	   	for (j=0;j<hau;j++)
   		{
   			buf2[j]= dest.Image[j*larg+i];
   		}
   		
   		// traitement image * h    		
   	    for (j=0;j<(hau/2);j++)
   	    {
   	    	buffer[2*j]=buf2[j];
   	    	buffer[2*j+1]=0;
   	    }  
   	    cp = 0;
   		for (j=0;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=(ca+cp)/2;
   			cp=ca;
   		}  
   		cp = buffer[hau-1];
   		for (j=hau-1;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=(ca+cp)/2;
   			cp=ca;
   		}
    	cp = buffer[0];
   		for (j=0;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=(ca+cp)/2;
   			cp=ca;
   		}  
   		cp = buffer[hau-1];
   		for (j=hau-1;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=(ca+cp)/2;
   			cp=ca;
   		}        
   		
   		cp= buffer[1];
   		for (j=2;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=ca+al*(cp-buffer[j-2]);
   			cp=ca;
   		}  
   		cp=buffer[hau-2];
   		for (j=hau-3;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=ca+al*(cp-buffer[j+2]);
   			cp=ca;
   		}   
   		cp = buffer[1];
   		for (j=2;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=ca+be*(cp-buffer[j-2]);
   			cp=ca;
   		}  
   		cp=buffer[hau-2];
   		for (j=hau-3;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=ca+be*(cp-buffer[j+2]);
   			cp=ca;
   		}                   
   		cp = buffer[1];
   		for (j=2;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=ca+ga*(cp-buffer[j-2]);
   			cp=ca;
   		}  
   		cp=buffer[hau-2];
   		for (j=hau-3;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=ca+ga*(cp-buffer[j+2]);
   			cp=ca;
   		}             
   		
   		for (j=0;j<hau;j++)
   		{
   			dest.Image[j*larg+i] =2*buffer[j];
   		}
    
  		// traitement image * g 
    	for (j=0;j<(hau/2);j++)
   		{
   			buffer[2*j]=buf2[j+(hau/2)];
   			buffer[2*j+1]=0;
   		}
   	    cp = 0;
   		for (j=0;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=(-ca+cp)/2;
   			cp=ca;
   		}  
   		cp = 0;
   		for (j=hau-1;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=(-ca+cp)/2;
   			cp=ca;
   		}
    	cp = 0;
   		for (j=0;j<hau;j++)
   		{          
   			ca=buffer[j];
   			buffer[j]=(-ca+cp)/2;
   			cp=ca;
   		}  
   		cp = 0;
   		for (j=hau-1;j>=0;j--)
   		{          
   			ca=buffer[j];
   			buffer[j]=(-ca+cp)/2;
   			cp=ca;
   		}                      
   		for (j=1;j<hau;j++)
   		{
   			dest.Image[j*larg+i] = dest.Image[j*larg+i]+2*buffer[j-1];
   		}
    }
   delete buffer; 
   delete buf2;
 
 	// Traitement horizontal 
   	buffer 	= new float[lar];
   	buf2	= new float[lar];

   	for (j=0;j<hau;j++)
   	{
   	    for (i=0;i<lar;i++)
   		{
   			buf2[i]= dest.Image[j*larg+i];
   		}
   	    // traitement image * h 
   	    
   	    for (i=0;i<(lar/2);i++)
   		{
   	    	buffer[2*i]=buf2[i];
   	    	buffer[2*i+1]=0;
   	    }  
   		cp = 0;
   		for (i=0;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=(ca+cp)/2;
   			cp=ca;
   		}  
   		cp = buffer[lar-1];
   		for (i=lar-1;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=(ca+cp)/2;
   			cp=ca;
   		}
    	cp = buffer[0];
   		for (i=0;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=(ca+cp)/2;
   			cp=ca;
   		}  
   		cp = buffer[lar-1];
   		for (i=lar-1;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=(ca+cp)/2;
   			cp=ca;
   		}        
   		cp= buffer[1];
   		for (i=2;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=ca+al*(cp-buffer[i-2]);
   			cp=ca;
   		}  
   		cp=buffer[lar-2];
   		for (i=lar-3;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=ca+al*(cp-buffer[i+2]);
   			cp=ca;
   		}   
   		cp = buffer[1];
   		for (i=2;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=ca+be*(cp-buffer[i-2]);
   			cp=ca;
   		}  
   		cp=buffer[lar-2];
   		for (i=lar-3;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=ca+be*(cp-buffer[i+2]);
   			cp=ca;
   		}                   
   		cp = buffer[1];
   		for (i=2;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=ca+ga*(cp-buffer[i-2]);
   			cp=ca;
   		}  
   		cp=buffer[lar-2];
   		for (i=lar-3;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=ca+ga*(cp-buffer[i+2]);
   			cp=ca;
   		}             
   		for (i=0;i<lar;i++)
   		{
   			dest.Image[j*larg+i] =2*buffer[i];
   		}
  	
  		// traitement image * g 
    	for (i=0;i<(lar/2);i++)
   		{
   			buffer[2*i]=buf2[i+(lar/2)];
   			buffer[2*i+1]=0;
   		}
   	    cp = 0;
   		for (i=0;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=(-ca+cp)/2;
   			cp=ca;
   		}  
   		cp = 0;
   		for (i=lar-1;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=(-ca+cp)/2;
   			cp=ca;
   		}
    	cp = 0;
   		for (i=0;i<lar;i++)
   		{          
   			ca=buffer[i];
   			buffer[i]=(-ca+cp)/2;
   			cp=ca;
   		}  
   		cp = 0;
   		for (i=lar-1;i>=0;i--)
   		{          
   			ca=buffer[i];
   			buffer[i]=(-ca+cp)/2;
   			cp=ca;
   		}                      
   		for (i=1;i<lar;i++)
   		{
   			dest.Image[j*larg+i] = dest.Image[j*larg+i]+2*buffer[i-1];
   		}
    }
  	delete buffer;
  	delete buf2;         
  	
  	hau=hau*2;
  	lar=lar*2;
  	}
   // fin du traitement
}

void OnTraitOnd(SImageFloat im,SImageFloat dest,int prof,float *m_a)
{
	// Traitement d'une image ondelette
	long j,i,l,lar,larg,hau,haut; 
   	float a[10][2][2];
   	int pr,li,co,prr;    
   	
   	for ( pr = prof;pr>0;pr--)
   	{             
   		a[pr][0][0]=m_a[0+(pr*4)];
		a[pr][0][1]=m_a[1+(pr*4)];
		a[pr][1][0]=m_a[2+(pr*4)];
		a[pr][1][1]=m_a[3+(pr*4)];
	}   
   	/////////////////////////
	// a mettre dans image //
	/////////////////////////
	//CString texte;
	//CDialCoef boite;
	//for (int pr=prof;pr>0;pr--)
	//{
	//	texte.Format("%d",pr);
	//	boite.m_ProfAct=texte;
	//	boite.DoModal();
	//	m_a[0+(pr*4)]=boite.m_a00;
	//	m_a[1+(pr*4)]=boite.m_a01;
	//	m_a[2+(pr*4)]=boite.m_a10;
	//	m_a[3+(pr*4)]=boite.m_a11;
	//}


   	// traitement sur l'image de dest 
   	larg=im.Lmax-im.Lmin+1;    
   	lar = larg;
   	haut=im.Hmax-im.Hmin+1;
   	hau=haut;
	l=0;                      
   	for (j=im.Hmin;j<=im.Hmax;j++)
   	{            
    	for (i=im.Lmin;i<=im.Lmax;i++)
   		{
   			dest.Image[l++]=im.Image[j*im.Larg+i];
   		}
   	}      
	float F;
	for (prr=1;prr<=prof;prr++)
	{
		hau=hau/2;
		lar=lar/2;
	}	
	for (pr=prof;pr>0;pr--)
	{              
		for (li=0;li<2;li++)
		{
			for (co=0;co<2;co++)
			{
				for (j=li*hau;j<(li+1)*hau;j++)
				{
					for (i=co*lar;i<(co+1)*lar;i++)
					{                    
						F= dest.Image[j*larg+i];
						dest.Image[j*larg+i]=F*a[pr][li][co];
					}
				}
			}
			
		}
		hau=hau*2;
		lar=lar*2;
	}                             
//	dest->prof=prof;
}



