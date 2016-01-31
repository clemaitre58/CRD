#include "wavelet.h"

void ss_echant_lig(float  *im,float  *f, long haut,long larg,int echelle)
{
	// ---------------- sous - échantillonnage (colonne) ------------ 
	// a c e 0 0
	// f h j 0 0
	// k m o 0 0
	// p r t 0 0

	long cpt,i,j;
	cpt=0;
	for (j=0;j<haut;j++)
	{        
		cpt=0;
		for(i=0;i<larg;i=i+2)
			{
			f[cpt+j*larg*echelle]=im[i+j*larg*echelle];
			cpt++;	
			}
	}               

}
 
void sur_echant_lig(float  *im,float  *f, long haut,long larg,int echelle)
{


	long cpt,i,j;
	cpt=0;
	for (j=0;j<haut;j++)
	{        
		cpt=0;
		for(i=0;i<larg;i=i+2)
			{
			f[i+j*larg*echelle]=im[cpt+j*larg*echelle];
			cpt++;	
			}
	}               

} 

void ss_echant_col(float  *im,float  *f, long haut, long larg,int echelle)
{            
	long cpt,i,j;
	cpt=0;    
	for (j=0;j<haut;j=j+2)
		{    
		for(i=0;i<larg;i++)
			f[i+cpt*larg*echelle]=im[i+j*larg*echelle];	
		cpt++;	
		}       
}
              
void sur_echant_col(float  *im,float  *f, long haut, long larg, int echelle)
{            
	
	long cpt,i,j;
	cpt=0;    
	for (j=0;j<haut;j=j+2)
		{    
		for(i=0;i<larg;i++)
			f[i+j*larg*echelle]=im[i+cpt*larg*echelle];	
		cpt++;	
		}       
}
              
void reset(float  *f,long taille)
{ 
	long i;
	for(i=0;i<(long)taille;i++)
		f[i]=0.F;	
	
}                              


void convolution_ligne(float  *f, float  *f_conv,float  *filtre,int fin, int centre,long haut, long larg,int echelle)
{                 
long i,j,k,indice;
float tampon;              
for (j=0;j<haut;j++)
	{	
	for(i=0;i<larg;i++)
		{
			tampon=0.F;		
			for(k=fin;k>centre;k--)
			{
	//			if (i-(k-centre)<0) indice=0;    // traitement des bords
	//				else indice=i-(k-centre);
				if (i-(k-centre)<0) indice=-(i-(k-centre));        // symétrie
					else indice=i-(k-centre);
					
			tampon+=f[indice+j*larg*echelle]*filtre[k];			
			}   
			                        
			for(k=0;k<=centre;k++)
			{
			//	if (i+k<larg) indice=i+k;   // traitement des bords
			//		else indice=larg-1;      
				if (i+k<larg) indice=i+k;   
					else indice=2*(larg-1)-i-k;      
					
			tampon+=f[indice+j*larg*echelle]*filtre[centre-k];		
			}                                  
			f_conv[i+j*larg*echelle]=tampon;	
		}
	}
              
}	        

void convolution_colonne(float  *f, float  *f_conv,float  *filtre,int fin, int centre, long haut, long larg,int echelle)
{             
long i,j,k,indice;
float tampon;
for (i=0;i<larg;i++)
	{	
	for(j=0;j<haut;j++)
		{
			tampon=0.F;		
			for(k=fin;k>centre;k--)
			{
			//	if (j-(k-centre)<0) indice=0;	// Traitement des bords 
			//		else indice=j-(k-centre);
				if (j-(k-centre)<0) indice=-(j-(k-centre));       // symétrie
					else indice=j-(k-centre);
					
			tampon+=f[i+indice*larg*echelle]*filtre[k];			
			}   
			                        
			for(k=0;k<=centre;k++)
			{
			//	if (j+k<haut) indice=j+k;	//Traitement des bords
			//		else indice=haut-1;         
				if (j+k<haut) indice=j+k;                     // symétrie
					else indice=2*(haut-1)-j-k;         
					
			tampon+=f[i+indice*larg*echelle]*filtre[centre-k];		
			}                                  
			f_conv[i+j*larg*echelle]=tampon;	
		}
	}	
	
}                        

void min_max(float  *f,float  * min,float * max,long taille)
{          
	long i;
	float temp;
	*min=*max=f[0];		
	for(i =1 ;i<(long)taille;i++)
	{
	 	temp=f[i];
	 	if (temp>*max) *max=temp;
	 	if (temp<*min) *min=temp;
	}  
	if (min==max) max+=1;
}

                       
void affichage(float  *f,unsigned char  *im,long offset,float maxi,float mini,long haut,long larg,int echelle)
{   
	long i,j;
	long max;
	if (haut%2) max=haut/2+1;
	     else max=haut/2;
	
	for (j=0;j<(int)max;j++)	
		{    
			for(i=0;i<(int)(larg/2);i++)		
				im[i+j*larg*echelle+offset]=(unsigned char)(((f[i+j*larg*echelle]-mini)*255)/(maxi-mini));	
		}                       
	
}
void affichage(float  *f,float  *im,long offset,long haut,long larg,int echelle)
{   
	long i,j;            
	long max;
	if (haut%2) max=haut/2+1;
	     else max=haut/2;
	for (j=0;j<(int)max;j++)	
		{    
			for(i=0;i<(int)(larg/2);i++)		
				im[i+j*larg*echelle+offset]=f[i+j*larg*echelle];	
		}                       
	
}

void r_affichage(float  *f,unsigned char  *im,float maxi,float mini,long haut,long larg,int echelle)
{   
	long i,j;   
	float norm = 255/(maxi-mini);
	for (j=0;j<haut;j++)	
		{    
			for(i=0;i<larg;i++)		
			{   
				im[i+j*larg*(long)echelle]=(unsigned char) ((f[i+j*larg*(long)echelle]-mini)*norm);	
				
				//im[i+j*larg*echelle]=(unsigned char)f[i+j*larg*echelle];	
			}					
		}                       
	
}   
void r_affichage(float  *f,float  *im,long haut,long larg,int echelle)
{   
	long i,j;   
	for (j=0;j<haut;j++)	
		{    
			for(i=0;i<larg;i++)		
			{   
				im[i+j*larg*(long)echelle]=f[i+j*larg*(long)echelle];	
			}					
		}                       
	
}   

void recop_image(unsigned char  *im,float  *dest, long taille)	
{                                                       
long i;
for (i=0;i<(long)taille;i++)
	dest[i]=(float)im[i];
}
void recop_image(float  *im,float  *dest, long taille)	
{                                                       
long i;
for (i=0;i<(long)taille;i++)
	dest[i]=im[i];
}

void reconstruction(float  *h, int indice_h_centre,int indice_h_fin,float  *g,int indice_g_centre,int indice_g_fin,float  *v_j,float  *d_h,
 float  *d_v,float  *d_o, float  *f_conv,long larg, long haut,int echelle)
{   
	
	long taille;
    taille=larg*haut;
	float  * f;
	f = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
	float  * poub_v_j;
	poub_v_j = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
    float  * poub_d_h;
	poub_d_h = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
    float  * poub_d_v;
	poub_d_v = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
    float  * poub_d_o;
	poub_d_o = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
     
    float *poub_g=new float[indice_g_fin*sizeof(float)];
    float *poub_h=new float[indice_h_fin*sizeof(float)];
    
	long i;
	float tampi;
	// --- inversion des filtres pour la reconstruction -------
    for (i=0;i<=indice_h_fin;i++)
    	poub_h[i]=2*h[i];
 	for (i=0;i<=indice_g_fin;i++)
    	poub_g[i]=2*g[i];
       	
	indice_g_centre=indice_g_fin-indice_g_centre;
	indice_h_centre=indice_h_fin-indice_h_centre;


	for (i=indice_g_fin;i>=indice_g_fin/2+1;i--) 
		{   
			
			tampi=poub_g[i];
			poub_g[i]=poub_g[indice_g_fin-i];
			poub_g[indice_g_fin-i]=tampi; 
		}
	for (i=indice_h_fin;i>=indice_h_fin/2+1;i--) 
		{
			tampi=poub_h[i];
			poub_h[i]=poub_h[indice_h_fin-i];
			poub_h[indice_h_fin-i]=tampi; 
		}
		                            
	// --------- sur - echantillonnage et convolution ------
	reset(f_conv,taille);
    sur_echant_col(v_j,f,haut/echelle,larg/echelle,echelle);
    convolution_colonne(f,f_conv,poub_h,indice_h_fin,indice_h_centre,haut/echelle,larg/echelle,echelle);  	
    //  --- recopier le résultat dans le tableau d'origine ----
    for (i=0;i<(long)taille;i++)
    	 poub_v_j[i]=f_conv[i];
    	 
    reset(f,taille);
    reset(f_conv,taille);
    sur_echant_col(d_h,f,haut/echelle,larg/echelle,echelle);
    convolution_colonne(f,f_conv,poub_g,indice_g_fin,indice_g_centre,haut/echelle,larg/echelle,echelle);  	
	//  --- recopier le résultat dans le tableau d'origine ----
    for (i=0;i<(long)taille;i++)
    	 poub_d_h[i]=f_conv[i]+poub_v_j[i];
    reset(f,taille);    
    sur_echant_lig(poub_d_h,f,haut/echelle,larg/echelle,echelle);
    reset(f_conv,taille);    
	convolution_ligne(f,f_conv,poub_h,indice_h_fin,indice_h_centre,haut/echelle,larg/echelle,echelle);
    for (i=0;i<(long)taille;i++)
    	 poub_d_h[i]=f_conv[i];
   
     
	// --------------------- 2 partie ---
	reset(f,taille);
	reset(f_conv,taille);	
    sur_echant_col(d_v,f,haut/echelle,larg/echelle,echelle);
   
    convolution_colonne(f,f_conv,poub_h,indice_h_fin,indice_h_centre,haut/echelle,larg/echelle,echelle);  	
    //  --- recopier le résultat dans le tableau d'origine ----
    for (i=0;i<(long)taille;i++)
    	 poub_d_v[i]=f_conv[i];
    reset(f,taille);
    reset(f_conv,taille);
    sur_echant_col(d_o,f,haut/echelle,larg/echelle,echelle);
    convolution_colonne(f,f_conv,poub_g,indice_g_fin,indice_g_centre,haut/echelle,larg/echelle,echelle);  	
	//  --- recopier le résultat dans le tableau d'origine ----
    for (i=0;i<(long)taille;i++)
    	 poub_d_o[i]=f_conv[i]+poub_d_v[i];
    reset(f,taille);
    reset(f_conv,taille);    
    sur_echant_lig(poub_d_o,f,haut/echelle,larg/echelle,echelle);
	convolution_ligne(f,f_conv,poub_g,indice_g_fin,indice_g_centre,haut/echelle,larg/echelle,echelle);
	
    for (i=0;i<(long)taille;i++)
    	 poub_d_o[i]=f_conv[i];
	
	
	reset(f_conv,taille);    
    	 
    // --------------------------- On Somme ----
    
    for (i=0;i<(long)taille;i++)
    	//f_conv[i]=abs((d_o[i]+d_h[i])-Image[i]);
    	f_conv[i]=(poub_d_o[i]+poub_d_h[i]);    
                   
                   
    GlobalFreePtr(f);
    GlobalFreePtr(poub_v_j);
    GlobalFreePtr( poub_d_h);
	GlobalFreePtr( poub_d_v);
	GlobalFreePtr( poub_d_o);
    delete poub_g;
    delete poub_h;
}
	
void TOD (SFiltre h,SFiltre g,SImage source,SImage result,SImage reconst,SImage cont,BOOL ok_reconst,BOOL ok_cont,BOOL ok_basse,BOOL ok_extra,int echelle)
	{
    

// ---- divers ---
    int poub_echelle;
	long i;
	float mini,maxi;                  
	long haut,larg,taille;
	//unsigned char  * image;
// ---------- informations relatives a l'image -----
	haut=source.Haut;
	larg=source.Larg;
	taille=larg*haut;
	
// --- tampons pour la reconstruction -- 

	float  * d_h;
	d_h = (float  *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
	float  * d_v;
	d_v = (float  *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
    float  * d_o;
	d_o = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
    float  * v_j ;
	v_j = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));	
	
// --- tampons pour les calculs -----	

	float  * f;
	f = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
    
    float  * f_conv;
	f_conv = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
// ----------------------------------------------------------------------------------------
// 							Espace j+1
// ----------------------------------------------------------------------------------------	                
recop_image(source.Image,f,taille);	
reset(f_conv,taille);	
convolution_ligne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
reset(f_conv,taille);
convolution_colonne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);
for (i=0;i<(long)taille;i++)
	v_j[i]=f[i];

// --------------- Image destination (analyse multirésolution -----
    
min_max(f,&mini,&maxi,taille);
    
// On recopie l'image initiale afin de conserver lés détails précédents ---
for (i=0;i<(long)taille;i++)
	result.Image[i]=source.Image[i];  
	
affichage(f,result.Image,0,maxi,mini,haut/echelle,larg/echelle,echelle);	
	
// ----------------------------------------------------------------------------------------
// 							Ondelette -
// ----------------------------------------------------------------------------------------
            
recop_image(source.Image,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);
reset(f_conv,taille);                  
convolution_colonne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
for (i=0;i<(long)taille;i++)
	d_h[i]=f[i];

min_max(f,&mini,&maxi,taille );	 
affichage(f,result.Image,(long)((larg/echelle)/2),maxi,mini,haut/echelle,larg/echelle,echelle);	

// ----------------------------------------------------------------------------------------
// 							Ondelette |
// ----------------------------------------------------------------------------------------
recop_image(source.Image,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);
reset(f_conv,taille);                  
convolution_colonne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
for (i=0;i<(long)taille;i++)
	d_v[i]=f[i];

min_max(f,&mini,&maxi,taille );	 
affichage(f,result.Image,(long)((haut/echelle)/2)*larg,maxi,mini,haut/echelle,larg/echelle,echelle);	

// ----------------------------------------------------------------------------------------
// 							Ondelette /
// ----------------------------------------------------------------------------------------
recop_image(source.Image,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);
reset(f_conv,taille);                  
convolution_colonne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
for (i=0;i<(long)taille;i++)
	d_o[i]=f[i];

min_max(f,&mini,&maxi,taille );	 
affichage(f,result.Image,(long)((haut/echelle)/2)*larg+(larg/echelle)/2,maxi,mini,haut/echelle,larg/echelle,echelle);	

//------------------------------------------------------------------
//--------------------------------- reconstruction -----------------
//------------------------------------------------------------------


    
if (ok_reconst)
{   
	reset(f_conv,taille);
	reconstruction(h.Filtre,h.Centre,h.Fin,g.Filtre,g.Centre,g.Fin,v_j,d_h,d_v,d_o,f_conv,source.Larg,source.Haut,echelle);	
	// le résultat de la reconstruction est dans f-conv
	r_affichage(f_conv,reconst.Image,(float)255,(float)0,haut/echelle,larg/echelle,echelle);	                
	
} 

// -------- At At tention cette fonction va modifier echelle !!!!!!!
// --------           
//                   et v_j,d_h....

poub_echelle=echelle; 
if (ok_cont)
{   
	//--- On reconstruit jusqu'a l'échèlle de départ ---			
	reset(v_j,taille); // correspond aux a_j => on affecte 0
	while (poub_echelle>=1)
	{
	reset(f_conv,taille);
    reconstruction(h.Filtre,h.Centre,h.Fin,g.Filtre,g.Centre,g.Fin,v_j,d_h,d_v,d_o,f_conv,source.Larg,source.Haut,poub_echelle);    
    //-- il est necessaire de recopier v_j pour reconstruire au niveau sup
    for (i=0;i<taille;i++)
     v_j[i]=f_conv[i];    
     // --- par contre il faut maintenant annuler les details ---
     reset(d_o,taille);
     reset(d_v,taille);
     reset(d_h,taille);
    if (poub_echelle/2>=1) poub_echelle=poub_echelle/2;
		else poub_echelle=-2; // on sort du while	
     }
      
	min_max(f_conv,&mini,&maxi,taille );
	r_affichage(f_conv,cont.Image,maxi,mini,haut,larg,1);	                    
	
} 


poub_echelle=echelle; 
if (ok_basse)
{   
	//--- On reconstruit jusqu'a l'échèlle de départ ---			
	//reset(v_j,taille); // correspond aux a_j => on affecte 0
	 reset(d_o,taille);
     reset(d_v,taille);
     reset(d_h,taille);
    
	while (poub_echelle>=1)
	{
	reset(f_conv,taille);
    reconstruction(h.Filtre,h.Centre,h.Fin,g.Filtre,g.Centre,g.Fin,v_j,d_h,d_v,d_o,f_conv,source.Larg,source.Haut,poub_echelle);    
    //-- il est necessaire de recopier v_j pour reconstruire au niveau sup
    for (i=0;i<taille;i++)
     v_j[i]=f_conv[i];    
     // --- par contre il faut maintenant annuler les details ---
     reset(d_o,taille);
     reset(d_v,taille);
     reset(d_h,taille);
    if (poub_echelle/2>=1) poub_echelle=poub_echelle/2;
		else poub_echelle=-2; // on sort du while	
     }
      
	min_max(f_conv,&mini,&maxi,taille );
	r_affichage(f_conv,cont.Image,maxi,mini,haut,larg,1);	                    
	
} 
       


poub_echelle=echelle; 
if (ok_extra)
{   
	//--- On reconstruit jusqu'a l'échèlle de départ ---			
	//reset(v_j,taille); // correspond aux a_j => on affecte 0
echelle=echelle*2;	
// ----------------------------------------------------------------------------------------
// 							Détails -
// ----------------------------------------------------------------------------------------
            
recop_image(v_j,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
convolution_colonne(f_conv,f,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);

for (i=0;i<(long)taille;i++)
	d_h[i]=f[i];
// ----------------------------------------------------------------------------------------
// 							Détails |
// ----------------------------------------------------------------------------------------
recop_image(v_j,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
convolution_colonne(f_conv,f,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);


for (i=0;i<(long)taille;i++)
	d_v[i]=f[i];
// ----------------------------------------------------------------------------------------
// 							Détails /
// ----------------------------------------------------------------------------------------
recop_image(v_j,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);                  
convolution_colonne(f_conv,f,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
for (i=0;i<(long)taille;i++)
	d_o[i]=f[i];
	
echelle=echelle/2;	                                                          

// ------------- boucle de reconstruction  -----------------------------------    

	while (poub_echelle>=1)
	{
	reset(f_conv,taille);
    reconstruction(h.Filtre,h.Centre,h.Fin,g.Filtre,g.Centre,g.Fin,v_j,d_h,d_v,d_o,f_conv,source.Larg,source.Haut,poub_echelle);    
    //-- il est necessaire de recopier v_j pour reconstruire au niveau sup
    for (i=0;i<taille;i++)
     v_j[i]=f_conv[i];    
     // --- par contre il faut maintenant recalculer les details ---
	

    if (poub_echelle/2>=1) 
    	{
    		poub_echelle=poub_echelle/2;
    		// ----------------------------------------------------------------------------------------
// 							Détails -
// ----------------------------------------------------------------------------------------
            
recop_image(v_j,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/poub_echelle,larg/poub_echelle,poub_echelle);
reset(f,taille);
convolution_colonne(f_conv,f,g.Filtre,g.Fin,g.Centre,haut/poub_echelle,larg/poub_echelle,poub_echelle);
for (i=0;i<(long)taille;i++)
	d_h[i]=f[i];
// ----------------------------------------------------------------------------------------
// 							Détails |
// ----------------------------------------------------------------------------------------
recop_image(v_j,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/poub_echelle,larg/poub_echelle,poub_echelle);
reset(f,taille);
convolution_colonne(f_conv,f,h.Filtre,h.Fin,h.Centre,haut/poub_echelle,larg/poub_echelle,poub_echelle);
for (i=0;i<(long)taille;i++)
	d_v[i]=f[i];
// ----------------------------------------------------------------------------------------
// 							Détails /
// ----------------------------------------------------------------------------------------
recop_image(v_j,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/poub_echelle,larg/poub_echelle,poub_echelle);
reset(f,taille);                  
convolution_colonne(f_conv,f,g.Filtre,g.Fin,g.Centre,haut/poub_echelle,larg/poub_echelle,poub_echelle);
for (i=0;i<(long)taille;i++)
	d_o[i]=f[i];

    	}
		else poub_echelle=-2; // on sort du while	
     }
      
	min_max(f_conv,&mini,&maxi,taille );
	r_affichage(f_conv,cont.Image,maxi,mini,haut,larg,1);	                    
	
} 
    
    
    
// ---- On efface les variables nécessaire à la reconstruction ----
GlobalFreePtr(v_j);
GlobalFreePtr(d_o);
GlobalFreePtr(d_v);
GlobalFreePtr(d_h);
// ---- On efface les variables tampons ----
GlobalFreePtr(f);
GlobalFreePtr(f_conv);	

} 

//
//
//----------------------- Surcharge de TOD pour les Float ----------------------------
//
//

void TOD (SFiltre h,SFiltre g,SImageFloat source,SImageFloat result,SImageFloat reconst,SImageFloat cont,BOOL ok_reconst,BOOL ok_cont,BOOL ok_basse,int echelle)
	{
    

// ---- divers ---
    int poub_echelle;
	long i;
	long haut,larg,taille;
	//unsigned char  * image;
// ---------- informations relatives a l'image -----
	haut=source.Haut;
	larg=source.Larg;
	taille=larg*haut;
	
// --- tampons pour la reconstruction -- 

	float  * d_h;
	d_h = (float  *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
	float  * d_v;
	d_v = (float  *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
    float  * d_o;
	d_o = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
    float  * v_j ;
	v_j = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));	
	
// --- tampons pour les calculs -----	

	float  * f;
	f = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
    
    float  * f_conv;
	f_conv = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
// ----------------------------------------------------------------------------------------
// 							Espace j+1
// ----------------------------------------------------------------------------------------	                

recop_image(source.Image,f,taille);	
reset(f_conv,taille);	
convolution_ligne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
reset(f_conv,taille);
convolution_colonne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);
for (i=0;i<(long)taille;i++)
	v_j[i]=f[i];

// --------------- Image destination (analyse multirésolution -----
    
    
// On recopie l'image initiale afin de conserver lés détails précédents ---
for (i=0;i<(long)taille;i++)
	result.Image[i]=source.Image[i];  
	
affichage(f,result.Image,0,haut/echelle,larg/echelle,echelle);	
	
// ----------------------------------------------------------------------------------------
// 							Ondelette -
// ----------------------------------------------------------------------------------------
            
recop_image(source.Image,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);
reset(f_conv,taille);                  
convolution_colonne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
for (i=0;i<(long)taille;i++)
	d_h[i]=f[i];

affichage(f,result.Image,(long)((larg/echelle)/2),haut/echelle,larg/echelle,echelle);
	

// ----------------------------------------------------------------------------------------
// 							Ondelette |
// ----------------------------------------------------------------------------------------
recop_image(source.Image,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);
reset(f_conv,taille);                  
convolution_colonne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
for (i=0;i<(long)taille;i++)
	d_v[i]=f[i];

affichage(f,result.Image,(long)((haut/echelle)/2)*larg,haut/echelle,larg/echelle,echelle);	

// ----------------------------------------------------------------------------------------
// 							Ondelette /
// ----------------------------------------------------------------------------------------
recop_image(source.Image,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);
reset(f_conv,taille);                  
convolution_colonne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
for (i=0;i<(long)taille;i++)
	d_o[i]=f[i];

affichage(f,result.Image,(long)((haut/echelle)/2)*larg+(larg/echelle)/2,haut/echelle,larg/echelle,echelle);	

//------------------------------------------------------------------
//--------------------------------- reconstruction -----------------
//------------------------------------------------------------------


    
if (ok_reconst)
{   
	reset(f_conv,taille);
	reconstruction(h.Filtre,h.Centre,h.Fin,g.Filtre,g.Centre,g.Fin,v_j,d_h,d_v,d_o,f_conv,source.Larg,source.Haut,echelle);	
	// le résultat de la reconstruction est dans f-conv
	r_affichage(f_conv,reconst.Image,haut/echelle,larg/echelle,echelle);	                
	
} 

// -------- At At tention cette fonction va modifier echelle !!!!!!!
// --------           
//                   et v_j,d_h....
poub_echelle=echelle; 
if (ok_cont)
{   
	//--- On reconstruit jusqu'a l'échèlle de départ ---			
	reset(v_j,taille); // correspond aux a_j => on affecte 0
	while (poub_echelle>=1)
	{
	reset(f_conv,taille);
    reconstruction(h.Filtre,h.Centre,h.Fin,g.Filtre,g.Centre,g.Fin,v_j,d_h,d_v,d_o,f_conv,source.Larg,source.Haut,poub_echelle);    
    //-- il est necessaire de recopier v_j pour reconstruire au niveau sup
    for (i=0;i<taille;i++)
     v_j[i]=f_conv[i];    
     // --- par contre il faut maintenant annuler les details ---
     reset(d_o,taille);
     reset(d_v,taille);
     reset(d_h,taille);
    if (poub_echelle/2>=1) poub_echelle=poub_echelle/2;
		else poub_echelle=-2; // on sort du while	
     }
      
	r_affichage(f_conv,cont.Image,haut,larg,1);	                    
	
} 

poub_echelle=echelle; 
if (ok_basse)
{   
	//--- On reconstruit jusqu'a l'échèlle de départ ---			
	//reset(v_j,taille); // correspond aux a_j => on affecte 0
	 reset(d_o,taille);
     reset(d_v,taille);
     reset(d_h,taille);
   
	while (poub_echelle>=1)
	{
	reset(f_conv,taille);
    reconstruction(h.Filtre,h.Centre,h.Fin,g.Filtre,g.Centre,g.Fin,v_j,d_h,d_v,d_o,f_conv,source.Larg,source.Haut,poub_echelle);    
    //-- il est necessaire de recopier v_j pour reconstruire au niveau sup
    for (i=0;i<taille;i++)
     v_j[i]=f_conv[i];    
     // --- par contre il faut maintenant annuler les details ---
     reset(d_o,taille);
     reset(d_v,taille);
     reset(d_h,taille);
    if (poub_echelle/2>=1) poub_echelle=poub_echelle/2;
		else poub_echelle=-2; // on sort du while	
     }
      
	r_affichage(f_conv,cont.Image,haut,larg,1);	                    
	
} 

       


    
    
    
// ---- On efface les variables nécessaire à la reconstruction ----
GlobalFreePtr(v_j);
GlobalFreePtr(d_o);
GlobalFreePtr(d_v);
GlobalFreePtr(d_h);
// ---- On efface les variables tampons ----
GlobalFreePtr(f);
GlobalFreePtr(f_conv);	

} 

// ------------------------------------------------------------------------------------------
// -
// -
// -
// --------------------- Procédure Permettant de reconstruire avec les détails précédents ---
// -
// -
// -
// ------------------------------------------------------------------------------------------
void TOD_EX (SFiltre h,SFiltre g,SImage source,SImage dest,SImage extrac,float * *tab_v_j,float * *tab_d_h,float * *tab_d_v,float * *tab_d_o,int echelle,int cpt)
{
    

// ---- divers ---
    int poub_echelle;
	long i;
	float mini,maxi;                  
	long haut,larg,taille;
	//unsigned char  * image;
// ---------- informations relatives a l'image -----
	haut=source.Haut;
	larg=source.Larg;
	taille=larg*haut;

// --- tampons pour les calculs -----	
    
	float  * f;
	f = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
    
    float  * f_conv;
	f_conv = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
// ----------------------------------------------------------------------------------------
// 							Espace j+1
// ----------------------------------------------------------------------------------------	                

recop_image(source.Image,f,taille);	
reset(f_conv,taille);	
convolution_ligne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
reset(f_conv,taille);
convolution_colonne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);
for (i=0;i<(long)taille;i++)
	tab_v_j[cpt][i]=f[i];
min_max(f,&mini,&maxi,taille);
    
// On recopie l'image initiale afin de conserver lés détails précédents ---
for (i=0;i<(long)taille;i++)
	dest.Image[i]=source.Image[i];  
	
affichage(f,dest.Image,0,maxi,mini,haut/echelle,larg/echelle,echelle);	

// ----------------------------------------------------------------------------------------
// 							Ondelette -
// ----------------------------------------------------------------------------------------
            
recop_image(source.Image,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);
reset(f_conv,taille);                  
convolution_colonne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
for (i=0;i<(long)taille;i++)
	tab_d_h[cpt][i]=f[i];
min_max(f,&mini,&maxi,taille );	 
affichage(f,dest.Image,(long)((larg/echelle)/2),maxi,mini,haut/echelle,larg/echelle,echelle);	

// ----------------------------------------------------------------------------------------
// 							Ondelette |
// ----------------------------------------------------------------------------------------
recop_image(source.Image,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);
reset(f_conv,taille);                  
convolution_colonne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
for (i=0;i<(long)taille;i++)
	tab_d_v[cpt][i]=f[i];
	
min_max(f,&mini,&maxi,taille );	 
affichage(f,dest.Image,(long)((haut/echelle)/2)*larg,maxi,mini,haut/echelle,larg/echelle,echelle);	

// ----------------------------------------------------------------------------------------
// 							Ondelette /
// ----------------------------------------------------------------------------------------
recop_image(source.Image,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);
reset(f_conv,taille);                  
convolution_colonne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
for (i=0;i<(long)taille;i++)
	tab_d_o[cpt][i]=f[i];
	
min_max(f,&mini,&maxi,taille );	 
affichage(f,dest.Image,(long)((haut/echelle)/2)*larg+(larg/echelle)/2,maxi,mini,haut/echelle,larg/echelle,echelle);	

//------------------------------------------------------------------
//--------------------------------- reconstruction -----------------
//------------------------------------------------------------------


    
poub_echelle=echelle; 
//--- On reconstruit jusqu'a l'échèlle de départ ---			
reset(tab_v_j[cpt],taille); // correspond aux a_j => on affecte 0

// --- Il faut recopier les détails dans le tableau  tab_vj[ ---
		
while (poub_echelle>=1)
	{
	reset(f_conv,taille);
    reconstruction(h.Filtre,h.Centre,h.Fin,g.Filtre,g.Centre,g.Fin,tab_v_j[cpt],tab_d_h[cpt],tab_d_v[cpt],tab_d_o[cpt],f_conv,source.Larg,source.Haut,poub_echelle);    
    //-- il est necessaire de recopier v_j pour reconstruire au niveau sup
    if (cpt>0)
    	{
    		cpt=cpt-1;
    		for (i=0;i<taille;i++)
     		tab_v_j[cpt][i]=f_conv[i];    
        }
     // --- Modifier l'indice de Cpt afin de récuperer les détails précédents  ---
   
         
    if (poub_echelle/2>=1) poub_echelle=poub_echelle/2;
		else poub_echelle=-2; // on sort du while	
     }
      
min_max(f_conv,&mini,&maxi,taille );
r_affichage(f_conv,extrac.Image,maxi,mini,haut,larg,1);	                    

// ---- On efface les variables tampons ----
GlobalFreePtr(f);
GlobalFreePtr(f_conv);	

} 
//
//
// ---- surcharge de l'extraction 
//
//
 
void TOD_EX (SFiltre h,SFiltre g,SImageFloat source,SImageFloat dest,SImageFloat extrac,float * *tab_v_j,float * *tab_d_h,float * *tab_d_v,float * *tab_d_o,int echelle,int cpt)
{
    

// ---- divers ---
    int poub_echelle;
	long i;
	long haut,larg,taille;
	//unsigned char  * image;
// ---------- informations relatives a l'image -----
	haut=source.Haut;
	larg=source.Larg;
	taille=larg*haut;

// --- tampons pour les calculs -----	
    
	float  * f;
	f = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
    
    float  * f_conv;
	f_conv = (float   *) GlobalAllocPtr(GMEM_FIXED | GMEM_ZEROINIT,larg*haut*sizeof(float));   
	
// ----------------------------------------------------------------------------------------
// 							Espace j+1
// ----------------------------------------------------------------------------------------	                

recop_image(source.Image,f,taille);	
reset(f_conv,taille);	
convolution_ligne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
reset(f_conv,taille);
convolution_colonne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);
for (i=0;i<(long)taille;i++)
	tab_v_j[cpt][i]=f[i];
    
// On recopie l'image initiale afin de conserver lés détails précédents ---
for (i=0;i<(long)taille;i++)
	dest.Image[i]=source.Image[i];  
	
affichage(f,dest.Image,0,haut/echelle,larg/echelle,echelle);	

// ----------------------------------------------------------------------------------------
// 							Ondelette -
// ----------------------------------------------------------------------------------------
            
recop_image(source.Image,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);
reset(f_conv,taille);                  
convolution_colonne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
for (i=0;i<(long)taille;i++)
	tab_d_h[cpt][i]=f[i];
affichage(f,dest.Image,(long)((larg/echelle)/2),haut/echelle,larg/echelle,echelle);	

// ----------------------------------------------------------------------------------------
// 							Ondelette |
// ----------------------------------------------------------------------------------------
recop_image(source.Image,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);
reset(f_conv,taille);                  
convolution_colonne(f,f_conv,h.Filtre,h.Fin,h.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
for (i=0;i<(long)taille;i++)
	tab_d_v[cpt][i]=f[i];
	
affichage(f,dest.Image,(long)((haut/echelle)/2)*larg,haut/echelle,larg/echelle,echelle);	

// ----------------------------------------------------------------------------------------
// 							Ondelette /
// ----------------------------------------------------------------------------------------
recop_image(source.Image,f,taille);
reset(f_conv,taille);
convolution_ligne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_lig(f_conv,f,haut/echelle,larg/echelle,echelle);
reset(f_conv,taille);                  
convolution_colonne(f,f_conv,g.Filtre,g.Fin,g.Centre,haut/echelle,larg/echelle,echelle);
reset(f,taille);
ss_echant_col(f_conv,f,haut/echelle,larg/echelle,echelle);                                     
for (i=0;i<(long)taille;i++)
	tab_d_o[cpt][i]=f[i];
	
affichage(f,dest.Image,(long)((haut/echelle)/2)*larg+(larg/echelle)/2,haut/echelle,larg/echelle,echelle);	

//------------------------------------------------------------------
//--------------------------------- reconstruction -----------------
//------------------------------------------------------------------


    
poub_echelle=echelle; 
//--- On reconstruit jusqu'a l'échèlle de départ ---			
reset(tab_v_j[cpt],taille); // correspond aux a_j => on affecte 0

// --- Il faut recopier les détails dans le tableau  tab_vj[ ---
		
while (poub_echelle>=1)
	{
	reset(f_conv,taille);
    reconstruction(h.Filtre,h.Centre,h.Fin,g.Filtre,g.Centre,g.Fin,tab_v_j[cpt],tab_d_h[cpt],tab_d_v[cpt],tab_d_o[cpt],f_conv,source.Larg,source.Haut,poub_echelle);    
    //-- il est necessaire de recopier v_j pour reconstruire au niveau sup
    if (cpt>0)
    	{
    		cpt=cpt-1;
    		for (i=0;i<taille;i++)
     		tab_v_j[cpt][i]=f_conv[i];    
        }
     // --- Modifier l'indice de Cpt afin de récuperer les détails précédents  ---
   
         
    if (poub_echelle/2>=1) poub_echelle=poub_echelle/2;
		else poub_echelle=-2; // on sort du while	
     }
      
r_affichage(f_conv,extrac.Image,haut,larg,1);	                    

// ---- On efface les variables tampons ----
GlobalFreePtr(f);
GlobalFreePtr(f_conv);	

} 

void masque_convolution_2D(float  *source,float  *dest,SImage info,SMasque masq)
{                                                       
// ---------------- attention il faudra retourner le masque !!!!! ---------


long i,j,haut,larg;
long indice_x,indice_y;
int l,c;               
int centre_c,centre_l,f_larg,f_haut;          
float tampon=0.F,valeur=0.F; 
//------- recopie des structures image et masque
haut=info.Haut;
larg=info.Larg;                  

centre_c=masq.Centre_col;
centre_l=masq.Centre_lig;
f_larg=masq.Larg;
f_haut=masq.Haut;


// --- les bords sont symétrisés ---
for (j=0;j<haut;j++)
	{	
	for (i=0;i<larg;i++)
		{       
		tampon=0.F;
		for(l=0;l<centre_l;l++)
			{
			for (c=0;c<centre_c;c++)
			{                       
			       	indice_x=(i-(centre_c-c));
		           	indice_y=j-(centre_l-l);             
		           	// -- cas des bords --             
				   	if (indice_x<0) indice_x=indice_x*-1;
				   	if (indice_y<0) indice_y=indice_y*-1;
					valeur=source[indice_x+indice_y*larg]*masq.Filtre[c+l*f_larg];
				 tampon+=valeur;
				}
			for (c=centre_c;c<f_larg;c++)
				{ 
				  	indice_x=(i+(c-centre_c));
				  	indice_y=j-(centre_l-l);
					// -- cas des bords --  				  
					if (indice_x>larg-1) indice_x=larg-(indice_x-(larg-1));
					if (indice_y<0)  indice_y=indice_y*-1;
					valeur=source[indice_x+indice_y*larg]*masq.Filtre[c+l*f_larg];
			        tampon+=valeur;
				}
		    }    
		 for(l=centre_l;l<f_haut;l++)
			{
			for (c=0;c<centre_c;c++)
				{                   
					indice_x=(i-(centre_c-c));
					indice_y=j+(l-centre_l);
					// -- cas des bords --					  				  
					if (indice_y>haut-1) indice_y=haut-(indice_y-(haut-1));
					if (indice_x<0) indice_x=indice_x*-1;				   	
					valeur=source[indice_x+indice_y*larg]*masq.Filtre[c+l*f_larg];				
					tampon+=valeur;
				}
			for (c=centre_c;c<f_larg;c++)
				{                                   
					indice_x=(i+(c-centre_c));
					indice_y=j+(l-centre_l);
					// -- cas des bords --					  				  
					if (indice_y>haut-1) indice_y=haut-(indice_y-(haut-1));				
					if (indice_x>larg-1) indice_x=larg-(indice_x-(larg-1));
					valeur=source[indice_x+indice_y*larg]*masq.Filtre[c+l*f_larg];
			        tampon+=valeur;
				}        
				
		    }    
		 dest[i+j*larg]=tampon;  
		 }
	}
}

