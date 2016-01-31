//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!			fichiers include personnels													!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "etiquett.h"
#include "stdio.h"
#include "seuil.h"

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!			Défintion des fonctions Etiquetage et Centre de Gravité						!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

long Etiquett(SImage imS,SImage imEtiqu,unsigned char seuil)
{
  // dimension de l'Image Larg, Haut
  // dimension de la fenetre (x)imS.Lmin,imS.Lmax,(y)imS.Hmin,imS.Hmax
  long f;  // Compteur dans le tableau d'étiquettes
  long i;
  long j;
  long k;    
  long l;    // Compteur dans l'image
  long m;
  long mu;
  long Compteur=1;
  long beta1;
  long beta2;
  long beta3;
  float dilat;
  unsigned int fond=0;   // fond blanc : 255  fond noir : 0 
  long Taille=imS.Larg*imS.Haut;
  long larg=imS.Lmax-imS.Lmin+1;   
  long taille = (imS.Hmax-imS.Hmin+1)*(imS.Lmax-imS.Lmin+1);   

	HGLOBAL idTab_equ;
	HGLOBAL idEtiquette;
	   
	long  * Etiquette;
	long  * Tab_equ;  
	
			 
	idEtiquette = GlobalAlloc(GMEM_MOVEABLE,(Taille*sizeof(long)));            
	Etiquette =  (long  *) GlobalLock(idEtiquette);
	
	idTab_equ = GlobalAlloc(GMEM_MOVEABLE,(long)(Taille/4)*sizeof(long));            
	Tab_equ =  (long  *) GlobalLock(idTab_equ); 
	
			
	Seuillage1(imS,seuil,TRUE);
	f=0;
	    
	for (j=1;j<(long)(Taille/4);j++) Tab_equ[j]=j;   // initialisation de la table d'équivalence
  
          
  for (j=imS.Hmin;j<=imS.Hmax;j++)           //
  {                                      //   Parcours de
  	l=j*imS.Larg+imS.Lmin;                    //   l'image
  	for (i=imS.Lmin;i<=imS.Lmax;i++)         //
  	{	   
  		if ((j==imS.Hmin) || (i==imS.Lmin) || (i==imS.Lmax))   // 
  		{                                                 // Bord gauche, droit et haut
  		 	imS.Image[l]=fond;                                // considerés comme du fond
  		 	Etiquette[f]=0;                               //
  		}
		else
		{                                         
		 	if (imS.Image[l]==fond) Etiquette[f]=0;   
		 	else                                  
		 		if ((imS.Image[l-imS.Larg-1]==fond) && (imS.Image[l-imS.Larg]==fond) &&                  
		 		(imS.Image[l-imS.Larg+1]==fond) && (imS.Image[l-1]==fond)) Etiquette[f]=Compteur++;  
		 		else
		 		{
		 		 	if ((imS.Image[l-imS.Larg-1]!=fond) && ((Tab_equ[Etiquette[f-larg-1]]<=Tab_equ[Etiquette[f-larg]])
		 		 		|| (Etiquette[f-larg]==0)) && ((Tab_equ[Etiquette[f-larg-1]]<=Tab_equ[Etiquette[f-larg+1]])
		 		 		|| (Etiquette[f-larg+1]==0)) && ((Tab_equ[Etiquette[f-larg-1]]<=Tab_equ[Etiquette[f-1]])
		 		 		|| (Etiquette[f-1]==0)))
		 		 	{
		 		 		 Etiquette[f]=Tab_equ[Etiquette[f-larg-1]];
		 		 		 beta1=Etiquette[f-larg];
		 		 		 beta2=Etiquette[f-larg+1];
		 		 		 beta3=Etiquette[f-1]; 
		 		 	}
		 		 	else if ((imS.Image[l-imS.Larg]!=fond) && ((Tab_equ[Etiquette[f-larg]]<=Tab_equ[Etiquette[f-larg+1]])
		 		 				|| (Etiquette[f-larg+1]==0)) && ((Tab_equ[Etiquette[f-larg]]<=Tab_equ[Etiquette[f-1]])
		 		 				|| (Etiquette[f-1]==0)))
		 		 			{
		 		 		 		Etiquette[f]=Tab_equ[Etiquette[f-larg]];
		 		 		 		beta1=Etiquette[f-larg-1];
		 		 		 		beta2=Etiquette[f-larg+1];
		 		 		 		beta3=Etiquette[f-1]; 
		 		 			}
		 		 			else if ((imS.Image[l-imS.Larg+1]!=fond)  && ((Tab_equ[Etiquette[f-larg+1]]<=Tab_equ[Etiquette[f-1]])
		 		 						|| (Etiquette[f-1]==0)))
		 		 					{
		 		 		 				Etiquette[f]=Tab_equ[Etiquette[f-larg+1]];
		 		 		 				beta1=Etiquette[f-larg-1];
		 		 		 				beta2=Etiquette[f-larg];
		 		 		 				beta3=Etiquette[f-1]; 
		 		 					}
		 		 					else
		 		 					{
		 		 						Etiquette[f]=Tab_equ[Etiquette[f-1]];
		 		 		 				beta1=Etiquette[f-larg-1];
		 		 		 				beta2=Etiquette[f-larg];
		 		 		 				beta3=Etiquette[f-larg+1]; 
		 		 					}
		 		 	 
		 		 	 if (beta1!=0) 
		 		 		while (Tab_equ[beta1]!=Etiquette[f])
  						{
  							mu=Tab_equ[beta1];
  							Tab_equ[beta1]=Etiquette[f];
  							beta1=mu;
  						}         
  					
		 		 	 if (beta2!=0)
		 		 		while (Tab_equ[beta2]!=Etiquette[f])
  						{
  							mu=Tab_equ[beta2];
  							Tab_equ[beta2]=Etiquette[f];
  							beta2=mu;
  						}   
  					 if (beta3!=0) 
		 		 	 	while (Tab_equ[beta3]!=Etiquette[f])
  						{
  							mu=Tab_equ[beta3];
  							Tab_equ[beta3]=Etiquette[f];
  							beta3=mu;
  						} 
		 		 	 						
		 		}      
		}  		
  		l++;
  		f++;	                   
  	}            
  }	     		
   	  
   	  
  // Actualisation de la table 
  
   for (k=1;k<Compteur;k++)
   {
		m=k;
  		while (Tab_equ[m]!=m) m=Tab_equ[m];
  		Tab_equ[k]=m;
  	} 
  	
  	
    // Simplification de la table
    
    m=1;
	for(k=1;k<Compteur;k++)
	{
		if (Tab_equ[k]==k) Tab_equ[k]=m++;
		else Tab_equ[k]=Tab_equ[Tab_equ[k]];
	}
	 	
	// Remplacement des valeurs, taille des défauts
	
	if (m!=1)
	{ 
			// Répartition sur les valeurs de gris          
		dilat= (float) 255/(m-1);
		
//		for (k=1;k<Compteur;k++) Tab_equ[k]*=dilat;
				
		for (j=imS.Hmin;j<=imS.Hmax;j++)           //
		{                                      //   Parcours de
			l=j*imS.Larg+imS.Lmin;
			k=(j-imS.Hmin)*larg;                    //   l'image
			for (i=imS.Lmin;i<=imS.Lmax;i++)         //
			{
			 	if (Etiquette[k]==0) imEtiqu.Image[k] =fond;
				else 
				{
					Etiquette[k]=Tab_equ[Etiquette[k]];   
					imEtiqu.Image[k] = (unsigned char) (Etiquette[k]*dilat);  
				}	
				k++;
				l++;
			}
		} 	
    
		  
	}

	GlobalUnlock(idTab_equ);
  	GlobalFree(idTab_equ);
			
	GlobalUnlock(idEtiquette);
	GlobalFree(idEtiquette);  
	
	return m-1;	
}  


void CentresGravites(SImage imS,SImage imGravit, long NbEtiqu)
/*===============================================================================================
Paramètres:
    imS = structure SImage correspond a l image SOURCE
    imGravit = structure SImage correspond a l image Resultat                
	pNbEtiqu = Nb de taches apres etiquettage
=================================================================================================
Description:
===============================================================================================*/
{
	long  *Xc;  
	long  *Yc;
	long cpt;
	HGLOBAL id_Xc;
	HGLOBAL id_Yc;
	id_Xc = GlobalAlloc(GMEM_FIXED,(NbEtiqu * sizeof(long)));
	Xc = (LONG  *) GlobalLock(id_Xc);
	id_Yc = GlobalAlloc(GMEM_FIXED,(NbEtiqu * sizeof(long)));
	Yc = (LONG  *) GlobalLock(id_Yc);
    for (long l=0 ; l<(imS.Larg*imS.Haut) ; l++)
    	imGravit.Image[l] = 0;//imS.Image[l];
    
	for (long i=1 ; i<=NbEtiqu ; i++)
	{                    
		Xc[i]=0;
		Yc[i]=0;
		cpt=0;
		for(long k=0 ; k<imS.Haut ; k++)
		{
			for(long j=0 ; j<imS.Larg ; j++)
			{
				if (imS.Image[k*imS.Larg+j] == (unsigned char) (i*255/NbEtiqu))
				{       
					Xc[i]+=j;
					Yc[i]+=k;
					cpt++;
			    }           
			}
		}                     
		if (cpt!=0)
		{
			Xc[i]=(long)(Xc[i]/cpt);
			Yc[i]=(long)(Yc[i]/cpt);
			if((Yc[i]>0)&&(Yc[i]<imS.Haut)&&(Xc[i]>0)&&(Xc[i]<imS.Larg))
			{
//				imGravit.Image[(Yc[i]*imS.Larg+Xc[i]+1)]=255;
				imGravit.Image[(Yc[i]*imS.Larg+Xc[i])]=255;
//				imGravit.Image[((Yc[i]+1)*imS.Larg+Xc[i]+1)]=255;
//				imGravit.Image[((Yc[i]+1)*imS.Larg+Xc[i])]=255;
			}
		}
	}
 	GlobalUnlock(id_Xc);
  	GlobalFree(id_Xc);  
 	GlobalUnlock(id_Yc);
  	GlobalFree(id_Yc);  
}

