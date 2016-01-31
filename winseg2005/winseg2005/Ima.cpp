// Ima.cpp: implementation of the CIma class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
//#include "demovid.h"
#include "Ima.h"
#include "math.h"
#include "float.h"
//#include "FFTReal.h"
#include "svm.h"
#include "petrou.h"

#include <afxtempl.h>
#include <direct.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
#ifndef M_PI
# define M_PI    3.14159265358979323846
#endif

BYTE Lettre[][8]=

{

            {0x3C,0x24,0x24,0x24,0x24,0x24,0x24,0x3C},    // définition de 0
            {0x8, 0x18,0x28,0x08,0x08,0x08,0x08,0x1C},     // définition de 1
            {0x30,0x48,0x48,0x8,0x10,0x20,0x40,0x78},      // définition de 2
            {0x38,0x44,0x04,0x08,0x04,0x04,0x04,0x38},     // définition de 3 
            {0x08,0x18,0x28,0x48,0x88,0xFE,0x08,0x08},    // définition de 4 
            {0x3E,0x20,0x20,0x20,0x3C,0x02,0x02,0x3C},   // définition de 5
            {0x1E,0x10,0x20,0x20,0x3C,0x22,0x22,0x1C},   // définition de 6
            {0x3E,0x02,0x02,0x3E,0x08,0x10,0x10,0x20},    // définition de 7
            {0x3C,0x24,0x24,0x3C,0x24,0x24,0x24,0x3C},   // définition de 8
            {0x3C,0x24,0x24,0x3C,0x04,0x04,0x04,0x3C},   // définition de 9
            {0x18,0x24,0x42,0x42,0x7e,0x41,0x41,0x41},   // définition de A
            {0x7c,0x42,0x42,0x44,0x4c,0x42,0x42,0x7c},   // définition de B
            {0x7e,0x40,0x40,0x40,0x40,0x40,0x40,0x7e},   // définition de C
            {0x78,0x44,0x42,0x42,0x42,0x42,0x42,0x7c},   // définition de D
            {0x7e,0x40,0x40,0x40,0x78,0x40,0x40,0x7f},   // définition de E
            {0x7e,0x40,0x40,0x40,0x78,0x40,0x40,0x40},   // définition de F
            {0x3e,0x40,0x40,0x40,0x5c,0x42,0x42,0x3c},   // définition de E
            {0x42,0x42,0x42,0x42,0x7e,0x42,0x42,0x42},   // définition de H
            {0x10,0x10,0x10,0x10,0x10,0x10,0x10,0x10},   // définition de I
            {0x02,0x02,0x02,0x02,0x22,0x22,0x22,0x13},   // définition de J
            {0x24,0x28,0x50,0x60,0x50,0x48,0x44,0x42},   // définition de K
            {0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x7e},   // définition de L
            {0x63,0x55,0x49,0x41,0x41,0x41,0x41,0x41},   // définition de M
            {0x41,0x61,0x51,0x49,0x45,0x43,0x41,0x41},   // définition de N
            {0x3C,0x24,0x24,0x24,0x24,0x24,0x24,0x3C},    // définition de 0
            {0x7C,0x42,0x42,0x42,0x7c,0x40,0x40,0x40},    // définition de P
            {0x3C,0x42,0x42,0x42,0x42,0x4A,0x43,0x4E},    // définition de Q
            {0x7E,0x42,0x42,0x42,0x7c,0x48,0x44,0x42},    // définition de R
            {0x1E,0x20,0x20,0x10,0x80,0x40,0x20,0x3c},    // définition de S
            {0x7E,0x18,0x18,0x18,0x18,0x18,0x18,0x18},    // définition de T
            {0x42,0x42,0x42,0x42,0x42,0x42,0x42,0x3C},    // définition de U
            {0x41,0x41,0x22,0x22,0x14,0x14,0x08,0x08},    // définition de V
            {0x41,0x41,0x22,0x22,0x2A,0x14,0x14,0x14},    // définition de W
            {0x41,0x22,0x14,0x01,0x14,0x22,0x41,0x41},    // définition de X
            {0x41,0x22,0x14,0x08,0x08,0x08,0x08,0x08},    // définition de Y
            {0x7e,0x02,0x04,0x08,0x10,0x20,0x40,0x7e},    // définition de Z


} ;




//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CIma::CIma()
{
	m_minMax=0.1;
	m_TailleFFT=2048;


	TailleMasque=128;
	m_Npar=((TailleMasque/2)-1)*3;
	LFinal=TailleMasque;
	CFinal = TailleMasque;


	LCarre=240;
	CCarre=240;
	m_TailleFFT2D=CFinal;

	m_Echelle=TRUE;
	m_Tab=NULL;
	m_TabDistance=NULL;
	m_i_Nb_Modeles=0;
	m_TailleModele=NULL;
	m_Tab_Modele=NULL;
	//pBoostMc=NULL;
	Tmin=NULL;
	Tmax=NULL;
	model=NULL;

	m_PImage=NULL;
	m_Tete=NULL;
	m_Info=NULL;

	bContour=FALSE;


	_getcwd( CheminDepart.GetBuffer(_MAX_PATH), _MAX_PATH );
	CheminDepart.ReleaseBuffer();
	m_TabModule=(double*)malloc(sizeof(double)*(m_TailleFFT/2+1));
	
	//pfft=new FFTReal (m_TailleFFT);

//	pfft2D=new  FFTReal (m_TailleFFT2D);

	PImageCarre = (rgb*) malloc(CFinal*LFinal*sizeof(rgb)); // mémoire pour les données


//	_mkdir("Modeles");
//	m_CheminDesModeles		= "\\Modeles\\";


	InitPourFFT();

}

CIma::CIma(BYTE *pImage, int lar, int hau,int s) // reçoit le pointeur vers l'image couleur
{
	m_Tab_Modele=NULL;
	//pBoostMc=NULL;
	Tmin=NULL;
	Tmax=NULL;
	model=NULL;

	m_PImage=NULL;
	m_Tete=NULL;
	m_Info=NULL;
	
	TailleMasque=128;
	m_Npar=((TailleMasque/2)-1)*3;

	m_TailleFFT=2048;
	LFinal=TailleMasque;
	CFinal = TailleMasque;
	LCarre=240;
	CCarre=240;

		m_Tab=NULL;
	m_TabDistance=NULL;
	m_i_Nb_Modeles=0;
	m_TailleModele=NULL;
	m_Tab_Modele=NULL;

	bContour=FALSE;
	
	m_TailleFFT2D=CFinal;
	m_minMax=0.1;
	m_TabModule=(double*)malloc(sizeof(double)*(m_TailleFFT/2+1));
	PImageCarre = (rgb*) malloc(CFinal*LFinal*sizeof(rgb)); // mémoire pour les données
	m_Echelle=TRUE;
	_getcwd( CheminDepart.GetBuffer(_MAX_PATH), _MAX_PATH );
	CheminDepart.ReleaseBuffer();
	Init(pImage,lar,hau,s);
//	_mkdir("Modeles");
//	m_CheminDesModeles		= "\\Modeles\\";
	//pfft=new FFTReal (m_TailleFFT);
	//pfft2D=new  FFTReal (m_TailleFFT2D);
	
	InitPourFFT();

}

CIma::~CIma()
{
	free(PImageCarre);
	LibereMemoire();
	LibereModeles();
	//delete pfft;
	//delete pfft2D;
	if (m_TabModule!=NULL) free (m_TabModule);
	//if (pBoostMc!=NULL) delete pBoostMc;
	if (Tmin!=NULL) free(Tmin);
	if (Tmax!=NULL) free(Tmax);
	if (model!=NULL)svm_destroy_model(model);
	if (m_PImage!=NULL) free(m_PImage);
	if (m_Info!=NULL) delete m_Info;
	if (m_Tete!=NULL) delete m_Tete;

	// pour FFT
	for(int s=0;s<(N+1);s++)
		free (fct[s]);
	free(fct);



}

void CIma::Fin()
{
	LibereMemoire();
	LibereModeles();
	//delete pfft;
	//delete pfft2D;
	if (m_TabModule!=NULL) free (m_TabModule);
	//if (pBoostMc!=NULL) delete pBoostMc;

}

void CIma::LibereModeles()
{
	if (m_Tab_Modele!=NULL)
	{
		for (int i=0;i<m_i_Nb_Modeles;i++)
		{
			free(m_Tab_Modele[i]);
			m_Tab_Modele[i]=NULL;
		}
	}
	free(m_Tab_Modele);
	m_Tab_Modele=NULL;
	TabRefNom.RemoveAll();
	if (m_TailleModele!=NULL) free(m_TailleModele);
	m_TailleModele=NULL;


}


void CIma::LibereMemoire()
{

	if (m_Tab!=NULL) free(m_Tab);
	if (m_TabDistance!=NULL) free(m_TabDistance);
	m_Tab=NULL;
	m_TabDistance=NULL;


}



void CIma::Init(BYTE *pImage, int lar, int hau,int Seuil)
{
	m_TabDistance=NULL;
	m_Taille=lar*hau;
	m_Largeur=lar;
	m_Hauteur=hau;

	m_Tete = (BITMAPFILEHEADER*) new BYTE[sizeof(BITMAPFILEHEADER)];    
	m_Info = (BITMAPINFO*) new BYTE[ sizeof(BITMAPINFOHEADER)];
	m_Info->bmiHeader.biWidth=m_Largeur;
	m_Info->bmiHeader.biHeight=m_Hauteur;
	m_Info->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	m_Info->bmiHeader.biPlanes = 1;   
	m_Info->bmiHeader.biSizeImage = sizeof(rgb)*m_Largeur*m_Hauteur;
	m_Info->bmiHeader.biBitCount = 24;
	m_Info->bmiHeader.biClrUsed = 0;
	m_Info->bmiHeader.biCompression =  BI_RGB;
	m_Info->bmiHeader.biXPelsPerMeter = 0;
	m_Info->bmiHeader.biYPelsPerMeter = 0;
	m_Info->bmiHeader.biClrImportant = 0;

	/* pointeur membre pour image en niveau de gris */
	m_Tab=(BYTE * ) malloc(m_Taille);
	int j,	i=0;
	for (j=0;j<m_Taille;j++) {m_Tab[j]=((int)pImage[i]+(int)pImage[i+1]+(int)pImage[i+2])/3.0;i+=3;}
 
	/* seuillage fixe pour suivi de contour fond blanc pièce noire */
	// a modifier
	for (j=0;j<m_Taille;j++) {if (m_Tab[j]>Seuil) m_Tab[j]=127; else m_Tab[j]=0; }

}


void CIma::Contour(BYTE *pImage, int Larg, int Haut)
{
	int coul,l,c;

	for (c=0;c<Larg-1;c++)
	{
		for (l=0;l<Haut;l++)
		{
			for (coul=0;coul<3;coul++)
			{
				pImage[(Larg*l+c)*3+coul]=abs((int)pImage[(Larg*l+c)*3+coul]-(int)pImage[(Larg*l+c+1)*3+coul]);
			} // for coul
		} // for l
	} // for c

}


long CIma::Etiquette(BYTE *pImage, int Larg, int Haut, int Seuil)
{
	long f;  // Compteur dans le tableau d'étiquettes
	long i;
	long j;
	long k;    
	long l;    // Compteur dans l'image
	long mu;
	long Compteur=1;
	long beta1;
	long beta2;
	long beta3;
	long temp;
	float dilat;
	BYTE fond=0;   // fond blanc : 255  fond noir : 0
//	BYTE pixel1, pixel2, pixel3, pixel4;
	BYTE m; // on supppose qu'il n'y a pas plus de  255 formes dans l'image
	
	long x_debut, x_fin, y_debut, y_fin;

	long *Etiquette;
	long *Tab_equ;

	Etiquette = NULL;
	Tab_equ = NULL;
	long Taille = Larg * Haut;



	// passage en niveau de gris
	BYTE *pImNb;
	BYTE *pTemp;
	pTemp=(BYTE *) malloc(Taille);
	pImNb=(BYTE *) malloc(Taille);
	i=0;
	for (j=0;j<Taille;j++) {pImNb[j]=(pImage[i]+pImage[i+1]+pImage[i+2])/3.0;i+=3;pTemp[j]=0;}

	//seuillage à modifier pour le rendre variable
	Seuillage(pImNb,Taille,Seuil);

	for (j=0;j<256;j++) HistoSurf[j]=0;

	// diverses initialisations
	x_debut = 0;
	x_fin = Larg-1;
	y_debut = 0;
	y_fin = Haut-1;
	long larg = x_fin-x_debut+1;	   
	temp = long(Taille);
	Etiquette =  (long*)malloc(Taille*sizeof(long));           
	Tab_equ =(long*)malloc(temp*sizeof(long));

			
	f=0;
	    
	for (j=1;j<(long)(Taille/4);j++) Tab_equ[j]=j;   // initialisation de la table d'équivalence
  
          
  for (j=y_debut;j<=y_fin;j++)           //
  {                                      //   Parcours de
  	l=j*Larg+x_debut;                    //   l'image
  	for (i=x_debut;i<=x_fin;i++)         //
  	{	   
  		if ((j==y_debut) || (i==x_debut) || (i==x_fin))   // 
  		{                                                 // Bord gauche, droit et haut
  		 	pImNb[l]=fond;                                // considerés comme du fond
  		 	Etiquette[f]=0;                               //
  		}
		else
		{                                         
		 	if (pImNb[l]==fond) Etiquette[f]=0;   
		 	else                                  
		 		if ((pImNb[l-Larg-1]==fond) && (pImNb[l-Larg]==fond) &&                  
		 		(pImNb[l-Larg+1]==fond) && (pImNb[l-1]==fond)) Etiquette[f]=Compteur++;  
		 		else
		 		{
		 		 	if ((pImNb[l-Larg-1]!=fond) && ((Tab_equ[Etiquette[f-larg-1]]<=Tab_equ[Etiquette[f-larg]])
		 		 		|| (Etiquette[f-larg]==0)) && ((Tab_equ[Etiquette[f-larg-1]]<=Tab_equ[Etiquette[f-larg+1]])
		 		 		|| (Etiquette[f-larg+1]==0)) && ((Tab_equ[Etiquette[f-larg-1]]<=Tab_equ[Etiquette[f-1]])
		 		 		|| (Etiquette[f-1]==0)))
		 		 	{
		 		 		 Etiquette[f]=Tab_equ[Etiquette[f-larg-1]];
		 		 		 beta1=Etiquette[f-larg];
		 		 		 beta2=Etiquette[f-larg+1];
		 		 		 beta3=Etiquette[f-1]; 
		 		 	}
		 		 	else if ((pImNb[l-Larg]!=fond) && ((Tab_equ[Etiquette[f-larg]]<=Tab_equ[Etiquette[f-larg+1]])
		 		 				|| (Etiquette[f-larg+1]==0)) && ((Tab_equ[Etiquette[f-larg]]<=Tab_equ[Etiquette[f-1]])
		 		 				|| (Etiquette[f-1]==0)))
		 		 	{
		 		 		Etiquette[f]=Tab_equ[Etiquette[f-larg]];
		 		 		beta1=Etiquette[f-larg-1];
		 		 		beta2=Etiquette[f-larg+1];
		 		 		beta3=Etiquette[f-1]; 
		 		 	}
		 		 	else if ((pImNb[l-Larg+1]!=fond)  && ((Tab_equ[Etiquette[f-larg+1]]<=Tab_equ[Etiquette[f-1]])
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
	//	dilat= (float) 255/(m-1);
		dilat=1; // pour analyse
//		for (k=1;k<Compteur;k++) Tab_equ[k]*=dilat;
				
		for (j=y_debut;j<=y_fin;j++)           //
		{                                      //   Parcours de
			l=j*Larg+x_debut;
			k=(j-y_debut)*larg;                    //   l'image
			for (i=x_debut;i<=x_fin;i++)         //
			{
			 	if (Etiquette[k]==0) pTemp[k] =fond;
				else 
				{
					Etiquette[k]=Tab_equ[Etiquette[k]];   
					pTemp[k] = (unsigned char) (Etiquette[k]*dilat); 
				}	
				HistoSurf[pTemp[k]]++;
				k++;
				l++;
			}
		} 	
    
		  
	}


	i=0;
	for (j=0;j<Taille;j++)
	{
		pImage[i]=pTemp[j];
		pImage[i+1]=pTemp[j];
		pImage[i+2]=pTemp[j];
		i+=3;
	}


	if(Tab_equ)
		free(Tab_equ);
	if(Etiquette)
		free(Etiquette);
	free(pImNb);
	free(pTemp);		
	
	NbEtiq=m-1;
	return m-1;
}




void CIma::Seuillage(BYTE *pIm, long Taille, BYTE Seuil)
{
	for (int j=0;j<Taille;j++) {if (pIm[j]<Seuil) pIm[j]=0; else pIm[j]=255;}

}

void CIma::Erosion(BYTE *pImage, int Larg, int Haut)
{
	long p,i,j,c,l,cc,ll;
	long Taille = Larg * Haut;

	int tc,tl;
	tc=1;
	tl=1;
	BYTE min;

	// passage en niveau de gris
	BYTE *pImNb;
	pImNb=(BYTE *) malloc(Taille);
	BYTE *pTemp;
	pTemp=(BYTE *) malloc(Taille);
	i=0;
	for (j=0;j<Taille;j++) {pImNb[j]=(pImage[i]+pImage[i+1]+pImage[i+2])/3.0;i+=3;pTemp[j]=0;}


	for (l=tl;l<Haut-tl-1;l++)
	{
		for (c=tc;c<Larg-tc-1;c++)
		{
			min=255;
			for (ll=l-tl;ll<=l+tl;ll++)
			{
				for (cc=c-tc;cc<=c+tc;cc++)
				{
					p=Larg*(ll)+cc;
					if (pImNb[p]<min) min = pImNb[p];
				}
			}
			pTemp[Larg*l+c]=min;
		}
	}

	
	i=0;
	for (j=0;j<Taille;j++)
	{
		pImage[i]=pTemp[j];
		pImage[i+1]=pTemp[j];
		pImage[i+2]=pTemp[j];
		i+=3;
	}
	free(pImNb);
	free(pTemp);
}

void CIma::Masquage(BYTE *pImage, int Larg, int Haut, int c1, int l1, int c2, int l2)
{

int coul,l,c,p;
	for (l=0;l<Haut;l++)
	{
		for (c=0;c<Larg;c++)
		{
			for (coul=0;coul<3;coul++) 
			{
				p=(Larg*l+c)*3+coul;
				if (!(l>l1 && l<l2 && c>c1 && c<c2))
					pImage[p]=0;

			}
		}
	}

}

void CIma::Masquage(CPoint p1, CPoint p2)
{
	if (m_Tab==NULL) return;
	int l,c,p;
	int l1=m_Hauteur - p1.y;
	int l2=m_Hauteur - p2.y;
	int c1=p1.x;
	int c2=p2.x;
	for (l=0;l<m_Hauteur;l++)
	{
		for (c=0;c<m_Largeur;c++)
		{
				p=(m_Largeur*l+c);
				if (!(l>l1 && l<l2 && c>c1 && c<c2)) m_Tab[p]=0;
		}
	}

}


BOOL CIma::ChargeModeles()
{

	CString Fichier_Modele_courant;
	CString NomModele;

	m_i_Nb_Modeles			= 0;

	CString Fich_Modele	= "";
	CString buffer_Modele	= "";

	CFileFind finder;
	BOOL bWorking = finder.FindFile(CheminDepart+m_CheminDesModeles + "*.txtmod",0);


	while (bWorking) 
	{ 
		bWorking = finder.FindNextFile(); 

		if(!finder.IsDots() || !finder.IsDirectory())
		{
			NomModele=finder.GetFileName();
			NomModele=NomModele.Mid( 0, NomModele.GetLength() -7 );
			TabRefNom.Add(NomModele); 
			m_i_Nb_Modeles++;
		}
	}

	if (m_i_Nb_Modeles!=0)
	{
		m_Tab_Modele   = (double **) malloc(m_i_Nb_Modeles * sizeof(double *));
		if (m_TailleModele!=NULL) free(m_TailleModele);
		m_TailleModele=(int *) malloc(m_i_Nb_Modeles * sizeof(int));

		int i = 0;
		bWorking = finder.FindFile(CheminDepart+m_CheminDesModeles + "*.txtmod",0);
		while (bWorking) 
		{ 
		
			bWorking = finder.FindNextFile(); 

			if(!finder.IsDots() || !finder.IsDirectory())

			Fichier_Modele_courant = finder.GetFilePath();

			// charger fichier modele dans tableau m_TabModele
			CStdioFile FichierModele(Fichier_Modele_courant, CFile::modeRead);


			if (FichierModele == NULL ) 
				FichierModele.Close();


			FichierModele.ReadString(buffer_Modele);

			m_TailleModele[i] = atoi(buffer_Modele);

			m_Tab_Modele[i] = (double *)malloc(sizeof(double)*m_TailleModele[i]);	
																				

			// --- empile le flux du fichier modele dans un tableau.
			for(int j =0;j<m_TailleModele[i];j++)
			{
					FichierModele.ReadString(buffer_Modele);
					m_Tab_Modele[i][j] = (double) atof(buffer_Modele);

			} // endfor


			FichierModele.Close();
			if (m_Echelle) Normalise(	m_Tab_Modele[i],m_TailleModele[i]); // test avec normaliseation		

			i++;
		}
	} // si on a au moins un modèle à charger

if (m_i_Nb_Modeles==0) return FALSE; else return TRUE;

}

void CIma::CalculDistance()
{

	Barycentre.x = 0;
	Barycentre.y = 0;
	int i=0;
	
	pos = listePointContour.GetHeadPosition();		
	while(pos)
	{				
		currentPoint = listePointContour.GetAt(pos);
				
		Barycentre.x += currentPoint.x; 
		Barycentre.y += currentPoint.y; 

		listePointContour.GetNext(pos);		
	}

	Barycentre.x = Barycentre.x / nbelement;
	Barycentre.y = Barycentre.y / nbelement;

	pos = listePointContour.GetHeadPosition();
	while(pos)
	{
		currentPoint = listePointContour.GetAt(pos);
		Centre.x = (currentPoint.x - Barycentre.x)*(currentPoint.x - Barycentre.x);
		Centre.y = (currentPoint.y - Barycentre.y)*(currentPoint.y - Barycentre.y);
		m_TabDistance[i++]=sqrt((double)(Centre.x+Centre.y)); // derni
		
		listePointContour.GetNext(pos);
	}
//	CalculFFTProfile();
}

// ==========================================================================
//    Fonction de type structure qui renvoie la valeur la plus 
//    petite et son indice dans un tableau a 1 dimension.
// ==========================================================================


indicemin CIma::db_min(double *tab, unsigned long t)
{
 indicemin	R;
			R.min		= tab[0];
			R.indice	= 0;


unsigned long	i		= 0;

	for(i=1; i< t ;i++)
	{
			if(tab[i] < R.min) 
					{ R.min = tab[i]; R.indice = i; }			
	}
	
   return R; 
}

indicemin CIma::Identification()
{
/////////////////////////////////////////////////////////////////////////////
//         ANALYSE PAR MODELE
/////////////////////////////////////////////////////////////////////////////

// reçois le resultat final de chaque fichier modele.
m_db_final  = (double *)malloc(sizeof(double)*m_i_Nb_Modeles);	


for(int i=0;i<m_i_Nb_Modeles;i++)
 {

//		m_TailleAnalyse	=__min(nbelement,m_TailleModele[i]);
		m_TailleAnalyse	=m_TailleModele[i];
		//int M 	=__max(nbelement,m_TailleModele[i]);

		m_di_min  = (double *)malloc(sizeof(double)*m_TailleAnalyse); 
//		TabC	  = (double *)malloc(sizeof(double)*m_TailleAnalyse); //der
		TabC	  = (double *)malloc(sizeof(double)*m_TailleAnalyse*2); 
		int j,k;


		/////////////////////////////////////////////////////////////////////////////
		// discrimination de la taille des tableaux de +25%
		/////////////////////////////////////////////////////////////////////////////
		int pourcent = ( (abs(m_TailleModele[i] - nbelement) * 100) / m_TailleAnalyse );


		if(pourcent > pourcent/*20*/)  // pas de test
		{
			m_db_final[i] = DBL_MAX; 
		} 
		else
		{

			/////////////////////////////////////////////////////////////////////////////
			// rééchantillonnage avec règle de 3    ( M*i) / N
			/////////////////////////////////////////////////////////////////////////////
			for(j=0;j<m_TailleAnalyse;j++)
			{
				k=(j*nbelement)/m_TailleAnalyse;
				TabC[j] = m_TabDistance[k] ;

			}

			if (m_Echelle) Normalise(TabC,m_TailleAnalyse); // test avec normaliseation		

			// prépare rotation en doublant le tableau
			for(j=0;j<m_TailleAnalyse;j++)
			{
				TabC[j+m_TailleAnalyse]=TabC[j];
			}
			/////////////////////////////////////////////////////////////////////////////
			//         OUVRIRE UNE BOUCLE TAND QU'IL Y A DES ROTATIONS EFFECTUER
			/////////////////////////////////////////////////////////////////////////////			
			for(int k=0;k<m_TailleAnalyse;k++)
			{

			// faire la somme et stocker la valeur dans un tableau di_min.
			
			//m_di_min[k] = sum_fabs(m_TabDistance, m_Tab_Modele[i], m_TailleAnalyse, m_TailleModele[i]);

				m_di_min[k] = sum_fabs(TabC, m_Tab_Modele[i], m_TailleAnalyse, m_TailleModele[i]);
				TabC++; // on incrémente le pointeur = rotation
				// der	TabRotation(TabC, m_TailleAnalyse); // der
                
			} // EndFor k 'rotation'

			
			TabC-=m_TailleAnalyse; // on restitue le pointeur d'origine
			
			
			/////////////////////////////////////////////////////////////////////////////
			//              EFFECTUER LA VALEUR MIN DU TABLEAU ET LA STOCKER 
			/////////////////////////////////////////////////////////////////////////////

			// faire le db_min de l'ensemble du tableau d0i.

			indicemin	R;
			R = db_min(m_di_min,m_TailleAnalyse);

			m_db_final[i] = R.min/(double)m_TailleAnalyse; // moyenne du min
	
		} // fin else

		free(m_di_min); 
		free(TabC);
} // EndFor i 

/////////////////////////////////////////////////////////////////////////////
//              faire le db_final de l'ensemble m_db_final.
/////////////////////////////////////////////////////////////////////////////


	indicemin	Rf;
	Rf = db_min(m_db_final,m_i_Nb_Modeles);

	free(m_db_final);

	return Rf;
}

void CIma::AjouterRef(CString nomfichier)
{
	int i;
	CString ch;
	CStdioFile temp;
	CalculDistance();
	temp.Open(nomfichier,CFile::modeCreate|CFile::modeWrite);

	ch.Format("%d\n",nbelement);
	temp.WriteString(ch);
	for (i=0;i<nbelement;i++)
	{
		ch.Format("%f\n",m_TabDistance[i]);
		temp.WriteString(ch);
	}
	temp.Close();

	CalculFFTProfile();
	AjouterDonneeDansRef(nomfichier+".dat");

}

BOOL CIma::SuiviContour()
{
long LargIm = m_Largeur; 
long HautIm = m_Hauteur;

BOOL PasContour=TRUE;
int essai=0;

BordAZero();	// par sécurité

do
{
	bool FirstPoint = false;
	long i = 0,j = 0;

	nbelement = 0;

	listePointContour.RemoveAll();


	do
	{
		if (m_Tab[(j*LargIm)+i] !=0) 
		{
			firstPoint.x = i;
			firstPoint.y = j;
			FirstPoint = true;
		}		
		i++;
		if (i == LargIm)
		{
			i=0;
			j++;
		}

	}while(FirstPoint==false && j<HautIm);

	if (FirstPoint==FALSE) return FALSE; // si image toute noire, on quitte la fonction

	currentPoint = firstPoint;
	int r,choix,spp=0;
	r = 1;
	do
	{
		choix=spp+r;
		switch(choix)
		{
			case 0 : currentPoint.y++; spp=6; break;
			case 1 : currentPoint.y--; spp=2; break;
			case 2 : currentPoint.x++; spp=0; break;
			case 3 : currentPoint.x--; spp=4; break;
			case 4 : currentPoint.y--; spp=2; break;
			case 5 : currentPoint.y++; spp=6; break;
			case 6 : currentPoint.x--; spp=4; break;
			case 7 : currentPoint.x++; spp=0; break;
		}

		if(m_Tab[currentPoint.x + LargIm*currentPoint.y] !=0)
		{
			listePointContour.AddTail(currentPoint);
			r=1;
			nbelement++; // <- Ajout de code ALEX
		}
		else
		{
			r=0;
		}

	}while (firstPoint != currentPoint); 

	if (nbelement>1) PasContour=FALSE;
	else
	{
		m_Tab[(currentPoint.y*LargIm)+currentPoint.x] =0;
		essai++;
		if (essai>100) PasContour=FALSE;
	}

} while (PasContour);


if (!PasContour)
{
	if (m_TabDistance!=NULL) free(m_TabDistance);
	m_TabDistance=(double *) malloc(nbelement*sizeof(double));

	//////////////////////////////////////////////////////////////////////////////////////////////
	//    pour matérialiser le contour 
	//////////////////////////////////////////////////////////////////////////////////////////////

	//POSITION pos;
	pos = listePointContour.GetHeadPosition();
	long p=0;
	Cmin=m_Largeur;
	Cmax=0;
	Lmin=m_Hauteur;
	Lmax=0;
	while(pos)
	{
		currentPoint = listePointContour.GetAt(pos);
		p=currentPoint.x +LargIm*currentPoint.y;
		if(Cmin>currentPoint.x) Cmin=currentPoint.x;
		if(Cmax<currentPoint.x) Cmax=currentPoint.x;
		if(Lmin>currentPoint.y) Lmin=currentPoint.y;
		if(Lmax<currentPoint.y) Lmax=currentPoint.y;
		if (p<0 || p>m_Taille) p=0;
		m_Tab[p] = 255; 
		listePointContour.GetNext(pos);
	}
}

/*
Rayon=__max(abs(Cmin-Cmax),abs(Lmin-Lmax))*0.5;
CPoint Centre;
Centre.x=(Cmin+Cmax)/2;
Centre.y=(Lmin+Lmax)/2;
// visu cercle
for (int l=0;l<m_Hauteur;l++)
{ 
	for (int c=0;c<m_Largeur;c++)
	{
		int p=m_Largeur*l+c;
		double ll=l-Centre.y;
		double cc=c-Centre.x;
		double r=sqrt((ll*ll)+(cc*cc));
		if (r<=Rayon) m_Tab[p]=255;
	}
}*/
return TRUE;
} // fin suivi contour

inline double CIma::sum_fabs(double *tab1, double *tab2, unsigned long t1, unsigned long t2)
{
	unsigned long	i		= 0;
	unsigned long	t_min	= 0;
	double			db_sum	= 0;

	t_min = __min( t1, t2 );

	for(i=0; i < t_min ;i++)
	{
		db_sum += fabs( tab1[i] - tab2[i] )/*/t_min*/; // si on veut normaliser
	}

	return db_sum;
}

void CIma::TabRotation(double *tab, unsigned long t)
{
	unsigned	long	i	= 0;
	double  tmp = tab[t-1];

	for(i=t-1;i>0;i--)
	{
		tab[i] = tab[i-1];
	}

	tab[0] = tmp;
}

void CIma::Analyse()
{

	if (SuiviContour())
	{
		CalculDistance();

		indicemin R;
		R = Identification();	
		m_identite	= R.indice ;	
		m_min		= (float)R.min;		


		//Chaine.Format("i=%d    Taille = %d    min = %f Nom = %s",m_identite,nbelement,m_min,NomModele); // modif code ALEX
		if (m_Echelle==FALSE) m_minMax=6.0; // a modifier
		if (m_min<m_minMax)
		m_Chaine.Format("%s ",TabRefNom[m_identite]); 
		else
		m_Chaine="Inconnu";

	}
}



void CIma::BordAZero()
{
	// --- périmètre de sécurité --------------------------------------------- //

	int k;
	for (k=0;k<m_Largeur;k++) m_Tab[k]=0;
	for (k=m_Taille-m_Largeur+1;k<m_Taille;k++) m_Tab[k]=0;
	for (k=0;k<m_Hauteur;k++) m_Tab[k*m_Largeur]=0;
	for (k=0;k<m_Hauteur;k++) m_Tab[k*m_Largeur+1]=0;
	for (k=0;k<m_Hauteur;k++) m_Tab[k*m_Largeur+2]=0;
	for (k=0;k<m_Hauteur;k++) m_Tab[k*m_Largeur+3]=0;
	for (k=0;k<m_Hauteur;k++) m_Tab[k*m_Largeur+(m_Largeur-1)]=0;


}

void CIma::VisuNB(BYTE *pImage)
{
	int j,i=0;
	if (m_Tab==NULL) return;
	for (j=0;j<m_Taille;j++)
	{
		pImage[i]=m_Tab[j];
		pImage[i+1]=m_Tab[j];
		pImage[i+2]=m_Tab[j];
		i+=3;
	}

}





void CIma::InitAvecEtiquette(BYTE *pImage, int lar, int hau, int num)
{
	m_TabDistance=NULL;
	m_Taille=lar*hau;
	m_Largeur=lar;
	m_Hauteur=hau;
	// init boite enblobante
	Cmin=m_Largeur;
	Cmax=0;
	Lmin=m_Hauteur;
	Lmax=0;
	/* pointeur membre pour image en niveau de gris */
	if (m_Tab!=NULL) free(m_Tab);
	m_Tab=(BYTE * ) malloc(m_Taille);
	int j,	i=0;
	for (j=0;j<m_Taille;j++) 
	{
		if (pImage[i]==num) m_Tab[j]=127;
		else m_Tab[j]=0;
		i+=3;
	}
 

}




/*! 
 *  \fn void CIma::Negatif(BYTE *pIm, long Taille)
 *  \brief fonction membre permettant de calculer le negatif.
 *  \param *pIm un tableau d'octet(image)
 *  \param Taille un long contenant la taille du tableau(image)
 */

void CIma::Negatif(BYTE *pIm, long Taille)
{
	DWORD *Ptr;
	Ptr=(DWORD *) pIm;
	Taille=Taille>>2;
	for (int j=0;j<Taille;j++) {Ptr[j]=0xFFFFFFFF-Ptr[j];}
//	for (int j=0;j<Taille;j++) {pIm[j]=255-pIm[j];}


}

void CIma::Expand(BYTE *pIm, long Taille)
{

if(NbEtiq!=0)	for (int j=0;j<Taille;j++) pIm[j]=((long)pIm[j]*255)/NbEtiq;


}

void CIma::IncrusteBoite(BYTE *pIm,int num)
{
	int l,c;
	for (c=Cmin;c<=Cmax;c++) 
	{
		pIm[(m_Largeur*Lmin+c)*3+2]=num;
//		pIm[(m_Largeur*Lmin+c)*3]=255;
		pIm[(m_Largeur*Lmax+c)*3+2]=num;
//		pIm[(m_Largeur*Lmax+c)*3]=255;
	}
	for (l=Lmin;l<=Lmax;l++) 
	{
		pIm[(m_Largeur*l+Cmin)*3+2]=num;
		pIm[(m_Largeur*l+Cmax)*3+2]=num;
//		pIm[(m_Largeur*l+Cmin)*3]=255;
//		pIm[(m_Largeur*l+Cmax)*3]=255;
	}

}

void CIma::Derive(float *tab1, unsigned long t1)
{
	unsigned long i;
	
	CString ch;
	CStdioFile temp;
	temp.Open("prof.txt",CFile::modeCreate|CFile::modeWrite);
	for(i=0; i < t1 ;i++)
	{
		ch.Format("%f,\n",tab1[i]);
		temp.WriteString(ch);
	}
	temp.Close();
/*
	for(i=0; i < t1-2 ;i++)
	{
		tab1[i]=  -tab1[i] + 2.0*tab1[i+1] -tab1[i+2];
	}
	tab1[t1-2]=tab1[t1-3];
	tab1[t1-1]=tab1[t1-2];
*/
}

void CIma::Normalise(double *tab1, unsigned long t1)
{
	unsigned long i;
	double max=tab1[0];
	for(i=1; i < t1 ;i++)
	{
		if (tab1[i]>max) max=tab1[i];
	}
	for(i=0; i < t1 ;i++)
	{
		tab1[i]=tab1[i]/max;
	}


}


void CIma::CalculFFTProfile()
{

	int nb2=m_TailleFFT>>1;
	int i,j,k;
	double	real;
	double img;

	if (nbelement==0) return;
	float *Tab=(float*)malloc(sizeof(float)*m_TailleFFT);
	float *fft=(float*)malloc(sizeof(float)*m_TailleFFT);
	
	Normalise(m_TabDistance,nbelement);
	
	for(j=0;j<m_TailleFFT;j++)
	{
		k=(j*nbelement)/m_TailleFFT;
		Tab[j] = m_TabDistance[k] ;

	}
//	Derive(Tab,m_TailleFFT);


	//pfft->do_fft(fft,Tab); // attention fft supprimee A REMETTRE avec l'autre fft
	double Accu=0;
	for (i = 0; i <= nb2; i ++)
	{
		real = fft [i];
		if (i > 0 && i < nb2)
		{
			img = fft [i + nb2];
		}
		else
		{
			img = 0;
		}
		m_TabModule[i]=sqrt (real * real + img * img);
		m_TabModule[i]=(real * real) + (img * img);
		Accu+=m_TabModule[i];
	}

//	Normalise(m_TabModule,nb2);

	for(i=0; i <nb2 ;i++)
	{
		m_TabModule[i]=(m_TabModule[i+1]*255)/Accu; // i+1 permet de ne pas prendre le premier
		if (m_TabModule[i]>0.98) m_TabModule[i]=0.98;
	}

	free (Tab);
	free (fft);


}

void CIma::AjouterDonneeDansRef(CString nomfichier)
{
	int par,nbd;
	CStdioFile Fic;
	CString ChaineVal,ChainePar,ChaineDonnee,ChaineTemp;
	ChainePar.Format("%d\n",m_Npar);
	ChaineDonnee.Format("%010d\n",1);
	CString Temp=nomfichier.Right(4);
	if (Temp!=".dat") nomfichier+=".dat";

	bool bfic=Fic.Open(nomfichier,CFile::modeCreate|CFile::modeReadWrite|CFile::modeNoTruncate,NULL);
	if (bfic)
	{
		if (Fic.GetLength()!=0) // le fichier existait
		{
			Fic.ReadString(ChaineTemp);
			Fic.ReadString(ChaineDonnee);
			nbd=atoi(ChaineDonnee);
			nbd+=1;
			ChaineDonnee.Format("%010d\n",nbd);
			Fic.SeekToBegin();
			Fic.WriteString(ChainePar);
			Fic.WriteString(ChaineDonnee);
			Fic.SeekToEnd();
		}
		else
		{
			nbd=1;
			ChaineDonnee.Format("%010d\n",nbd);
			Fic.WriteString(ChainePar);
			Fic.WriteString(ChaineDonnee);
		}
			
		for(par=0;par<m_Npar;par++) 
 		{ 
			ChaineVal.Format("%f ",m_TabModule[par]);
			Fic.WriteString(ChaineVal);
		}//for par
		Fic.WriteString("\n");
		Fic.Flush();
		Fic.Close();
	}
}

BOOL CIma::ChargeBoost(CString NomBoost)
{
/*
	if (pBoostMc!=NULL) delete pBoostMc;
	pBoostMc=new CBoostMcBis();
	
	pBoostMc->Lit(NomBoost);
		if (pBoostMc->Npar!=m_Npar) return FALSE;
		else	return TRUE;
*/
	return FALSE;
}

void CIma::AnalyseAvecBoosting()
{
	/*
	int par,cf;
//	BYTE b;
	SBoostOutput Fin;
	if (pBoostMc==NULL) return;
	float *x_test=new float[m_Npar];
	if (SuiviContour())
	{
// distance et FFT !
		CalculDistance();
		CalculFFTProfile();
		for (par=0;par<m_Npar;par++) x_test[par]=m_TabModule[par];
		Fin=pBoostMc->FinalClass(x_test);
		cf=Fin.cf;
		m_Chaine.Format("%d",cf); 
	} // si suivi contour ok
	delete [] x_test;*/
}



void CIma::IncrusteChaine(BYTE *pImage, int x, int y, CString ch,int num)
{
	long p;
	int l,c,cc;
	int Cini,Lini;
	int i;
	i=-1;
	// pour un caractère
	int lon=ch.GetLength();
	if (lon>3) lon=3;
	for (cc=0;cc<lon;cc++)
	{
		char cara=ch[cc]; 
		if (cara>='a' && cara <='z')
			i=cara-'a'+10;
		if (cara>='A' && cara <='Z')
			i=cara-'A'+10;
		if (cara>='0' && cara <='9')
			i=cara-'0';
		if (i==-1) return;// pas prévu pour l'instant

		int bin;
		Cini=x;
		Lini=y;
		for (l=0;l<16;l++)
		{
			for (c=0;c<16;c++)
			{
					p=(m_Largeur*(m_Hauteur-(l+Lini))+(c+Cini+cc*16))*3 ;
					if (p<0 || p>m_Taille*3) break;
					bin=(Lettre[i][l>>1]<<(c>>1))&0x80;
					if (bin) 
					{
						pImage[p]=0;pImage[p+1]=0;pImage[p+2]=num;
					}
					else 
					{
				//		pImage[p]=0;pImage[p+1]=0;pImage[p+2]=0;
					}

			} // for l
		} // for c

	}
}

void CIma::ReSampleCarre(BYTE *Image, int Larg, int Haut, CClientDC *pdc,CString Nom,BOOL dec,BOOL Pos,BOOL Sauve,CString NomFich)
{
 // on suppose que l'image fait 320x240
 // on prend 240x240
 // et on reechantillonne en 128x128
	int p,pp,cp,lp,l,c;
	float Ratio;
	int Taille=Larg*Haut*3;
	Ratio=(float)LCarre/(float)LFinal;

	if (Pos)
	{


		for (c=0;c<CFinal;c++)
		{
			cp=c*Ratio;
			for (l=0;l<LFinal;l++)
			{
				lp=l*Ratio;
				p=LFinal*l+c;
				pp=(Larg*lp+(cp+40))*3;
				if (pp>Taille) pp=0;
				PImageCarre[p].r=Image[pp];
				PImageCarre[p].g=Image[pp+1];
				PImageCarre[p].b=Image[pp+2];
			} // for l
		} // for c

		// pour visu carre

		// pour visu
		
		BITMAPINFO	 Info; // entête de l'image

		Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
		Info.bmiHeader.biWidth = CFinal;
		Info.bmiHeader.biHeight = LFinal;
		Info.bmiHeader.biPlanes = 1;   
		Info.bmiHeader.biSizeImage = sizeof(rgb)*CFinal*LFinal;
		Info.bmiHeader.biBitCount = 24;
		Info.bmiHeader.biClrUsed = 0;
		Info.bmiHeader.biCompression =  BI_RGB;
		Info.bmiHeader.biXPelsPerMeter = 0;
		Info.bmiHeader.biYPelsPerMeter = 0;
		Info.bmiHeader.biClrImportant = 0;
//		SetDIBitsToDevice(pdc->GetSafeHdc(),330,0, CFinal, LFinal , 0,0,0,LFinal,PImageCarre,&Info,DIB_RGB_COLORS);
		if (Sauve && !dec)
		{
			CFile Fichier; 
			if (Fichier.Open(NomFich,CFile::modeWrite|CFile::modeCreate))
			  {   
				BITMAPFILEHEADER Tete;
				Tete.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER); 
				Tete.bfSize = Tete.bfOffBits+CFinal*LFinal*sizeof(rgb);
				Tete.bfType = 'MB'; // obligatoire
				Tete.bfReserved1 = 0;
				Tete.bfReserved2 = 0;  
				Fichier.Write(&Tete,(UINT) sizeof(BITMAPFILEHEADER)); // entête du fichier
				Fichier.Write(&Info,sizeof(BITMAPINFOHEADER));   // entête de l'image
				Fichier.Write(PImageCarre,CFinal*LFinal*sizeof(rgb));	 // écriture des données
				Fichier.Close();

			}// si ouverture du fichier ok
	
		}
// fin visu
		
		// pretraitement eventuel
		if (bContour) Sobel(PImageCarre,CFinal,LFinal);
		// on affiche ce qui est traité
		SetDIBitsToDevice(pdc->GetSafeHdc(),330,0, CFinal, LFinal , 0,0,0,LFinal,PImageCarre,&Info,DIB_RGB_COLORS);

		double *DesFour;
		DesFour=DescripteursFourier1(PImageCarre);
		

		
		int i;
		for (i=0;i<m_Npar;i++)
		{
			m_TabModule[i]=DesFour[i];
		}
		if (!dec) AjouterDonneeDansRef(Nom); // si apprentissage
		else // si decision temps reel
		{ /*
			int par,cf;
			BYTE b;
			SBoostOutput Fin;
			if (pBoostMc==NULL) return;
			float *x_test=new float[m_Npar];
			for (par=0;par<m_Npar;par++) x_test[par]=m_TabModule[par+1];
			
			NormaliseAvecTable(x_test,m_Npar);

			Fin=pBoostMc->FinalClass(x_test);
			cf=Fin.cf;
			m_Chaine.Format("%d",cf); 
			delete [] x_test;*/

		
			struct svm_node *x;
			int par;
			double v;
			x = (struct svm_node *) malloc(m_Npar*sizeof(struct svm_node)*2);
			for(par=0;par<m_Npar;par++)
			{ 
				x[par].index=par+1;
				x[par].value=(double)m_TabModule[par];

			}//for par
			x[m_Npar].index = -1;
			v = svm_predict(model,x);
		//	v=1;
			m_Chaine=TabRefNom[(int)v-1];
			//m_Chaine.Format("%d",(int)v); 
			free(x);
		
		
		
		}
/*
		// pour visu
		
		pdc->Rectangle(0,500,400,280);
		pdc->MoveTo(0,500-DesFour[0]*200);
		
		for (i=0;i<m_Npar*4;i+=4)
		{
				pdc->LineTo(i,500-DesFour[i/4]*200);
		}
		
		// fin visu
*/
		if (DesFour!=NULL) free(DesFour);
	} // si position ok

	for (c=0;c<(Larg-CCarre)/2;c++)
	{
		for (l=0;l<Haut;l++)
		{
			pp=(Larg*l+c)*3;
			Image[pp]=127;
			Image[pp+1]=127;
			Image[pp+2]=127;
		} // for l
	} // for c
	for (c=Larg-(Larg-CCarre)/2;c<Larg;c++)
	{
		for (l=0;l<Haut;l++)
		{
			pp=(Larg*l+c)*3;
			Image[pp]=127;
			Image[pp+1]=127;
			Image[pp+2]=127;
		} // for l
	} // for c
	
}

double * CIma::DescripteursFourier2(rgb *pImage)
{

	
	
	
	double *coef=NULL;
	int p;
	long LargIm,HautIm;
	
	CStdioFile Fic;
	CString ChaineTemp;

	double pi=3.1415926535897932384626433832795;
	
	
	//N=1024; // initialisé dans InitPourFFT
	
	LargIm=LFinal;
	HautIm=CFinal;
	LargIm=N;
	HautIm=N;
	int N2=N/2;



	if (LargIm!=HautIm) return 0; // si l'image n'est pas carree c'est pas bon

	long l,c,i,j,m,r;



	// calcul de coefficient de correction
/*
	int s,r1;
	_complex **fct;
	_complex corr,cors;
	double dN=(double)N;
	double dN2=N/2;
	fct = (_complex **) malloc ((N+1)*sizeof(_complex *));

	for(s=0;s<(N+1);s++)
		fct[s]=	(_complex*)malloc((N+1)*sizeof(_complex));

	for( r1=1; r1<= N; r1++)
	{
		for (s=1;s<=N;s++)
		{
			corr=cor((r1-dN2)/dN);
			cors=cor((s-dN2)/dN);

			fct[r1][s].x= (corr.x*cors.x-corr.y*cors.y)/dN ;
			fct[r1][s].y= (corr.x*cors.y+corr.y*cors.x)/dN ;
		}
	}
*/

/*
for (s=2;s<=128;s++)
		{
			mu=(double)(s-1)/128;
			fct[1][s].x= (sin(2*pi*mu))/(2*pi*128*mu);
			fct[1][s].y= (1-cos(2*pi*mu))/(2*pi*128*mu);
}

for (r1=2;r1<=128;r1++)
		{
			lamda=(double)(r1-1)/128;
			fct[r1][1].x= (sin(2*pi*lamda))/(2*pi*128*lamda);
			fct[r1][1].y= (1-cos(2*pi*lamda))/(2*pi*128*lamda);
}
fct[1][1].x=1/128;
fct[1][1].y=0;


// recntrage de facteur correctif

long irc,jrc;
	_complex tft;
  for(j=1;j<=HautIm/2;j++)
		{
			jrc=(j+HautIm/2);//%HautIm;
			for(i=1;i<=LargIm;i++)
			{
				tft=fct[j][i];
				fct[j][i]=fct[jrc][i];
				fct[jrc][i]=tft;
			}
		}
		for(i=1;i<=LargIm/2;i++)
		{
			irc=(i+LargIm/2);//%LargIm;
			for(j=1;j<=HautIm;j++)
			{
				tft=fct[j][i];
				fct[j][i]=fct[j][irc];
				fct[j][irc]=tft;
			}
		}

*/
/*
Fic.Open("fctx.dat",CFile::modeWrite|CFile::modeCreate);
	for (l=1;l<=HautIm;l++)
	{
		for (c=1;c<=LargIm;c++)
		{
			ChaineTemp.Format("%f;",fct[l][c].x);
			Fic.WriteString(ChaineTemp);
		}
		Fic.WriteString("\n");
	}
	Fic.Close();

Fic.Open("fcty.dat",CFile::modeWrite|CFile::modeCreate);
	for (l=1;l<=HautIm;l++)
	{
		for (c=1;c<=LargIm;c++)
		{
			ChaineTemp.Format("%f;",fct[l][c].y);
			Fic.WriteString(ChaineTemp);
		}
		Fic.WriteString("\n");
	}
	Fic.Close();
*/	

//	float *TabData=(float*)malloc(m_TailleFFT2D*sizeof(float));
//	float *fft=(float*)malloc(sizeof(float)*m_TailleFFT2D);

	double* real;
	double* imag;

	_complex ***grid;
	grid = (_complex ***) malloc (3*sizeof(_complex **)); //new _complex**[3];
	for(m=0;m<3;m++)
	{
		grid[m] =(_complex**) malloc(LargIm*sizeof(_complex*));// new _complex*[LargIm];
		for(j=0;j<LargIm;j++) grid[m][j] =(_complex*)malloc(HautIm*sizeof(_complex));// new _complex[HautIm];
	}

	
/*	_complex ***grid1;
	grid1 = (_complex ***) malloc (3*sizeof(_complex **)); //new _complex**[3];
	for(m=0;m<3;m++)
	{
		grid1[m] =(_complex**) malloc(LargIm*sizeof(_complex*));// new _complex*[LargIm];
		for(j=0;j<LargIm;j++) grid1[m][j] =(_complex*)malloc(HautIm*sizeof(_complex));// new _complex[HautIm];
	}
*/
//	RGBQUAD colorReal,colorImag;

	//collecte des données
	
	for (l=0;l<HautIm;l++)  
	{
		for (c=0;c<LargIm;c++) 
		{
				
			//p=LargIm*l+c;
				grid[0][c][l].x = 0.0;
				grid[0][c][l].y = 0.0;
				grid[1][c][l].x = 0.0;
				grid[1][c][l].y =0.0;
				grid[2][c][l].x = 0.0;
				grid[2][c][l].y =0.0; 
			//complex z= complex (grid[0][c][l].x,grid[0][c][l].y);
		
		} //c
	} // l
	
	for (l=0;l<m_Hauteur;l++)  
	{
		for (c=0;c<m_Largeur;c++) 
		{
				
				p=m_Largeur*l+c;
				grid[0][c][l].x = pImage[p].r;
				grid[0][c][l].y = 0.0;//pImage[p].r;
				grid[1][c][l].x = pImage[p].g;
				grid[1][c][l].y =0.0;//pImage[p].g;
				grid[2][c][l].x = pImage[p].b;
				grid[2][c][l].y =0.0;// pImage[p].b; 
		
		} //c
	} // l

	/* Transform the rows */
	real = (double *)malloc(LargIm * sizeof(double));
	imag = (double *)malloc(LargIm * sizeof(double));


	m=0;
	while((1<<m)<HautIm) m++;

	for (r=0;r<3;r++)
	{
		for (l=0;l<HautIm;l++) 
		{
			for (c=0;c<LargIm;c++) 
			{
				real[c] = grid[r][c][l].x;
				imag[c] = grid[r][c][l].y;
			}

		//	pfft2D->do_fft(fft,TabData);

			FFT(1,m,real,imag);

			for (c = 0; c <LargIm; c ++)
			{
				grid[r][c][l].x = real[c];
				grid[r][c][l].y = imag[c];
			}

		}
	}

	/* Transform the columns */
	for (r=0;r<3;r++)
	{	
		for (c=0;c<LargIm;c++)
		{
			for (l=0;l<HautIm;l++) 
			{
	//			TabData[j] = grid[r][j][i].x;
				real[l] = grid[r][c][l].x;
				imag[l] = grid[r][c][l].y;
			}


			FFT(1,m,real,imag);

		for (l = 0; l <HautIm; l ++)
			{
				grid[r][c][l].x = real[l];
				grid[r][c][l].y = imag[l];
			} // ligne

		} // colonne
	} // couleur


/*
		Fic.Open("fftxlc.dat",CFile::modeWrite|CFile::modeCreate);
	for (l=0;l<HautIm;l++)
	{
		for (c=0;c<LargIm;c++)
		{
			ChaineTemp.Format("%f;",grid[0][l][c].x);
			Fic.WriteString(ChaineTemp);
		}
		Fic.WriteString("\n");
	}
	Fic.Close();

	Fic.Open("ffty.dat",CFile::modeWrite|CFile::modeCreate);
	for (c=0;c<LargIm;c++)
	{
		ChaineTemp.Format("%f\n",grid[0][0][c].y);
		Fic.WriteString(ChaineTemp);
	}
	Fic.Close();
*/

/*	Fic.Open("fft.dat",CFile::modeWrite|CFile::modeCreate);
	for (c=0;c<LargIm;c++)
	{
		ChaineTemp.Format("%f\n",_cabs(grid[0][65][c]));
		Fic.WriteString(ChaineTemp);
	}
	Fic.Close();
	*/

	RecentrageGrilleCouleur(grid, HautIm, LargIm);

/*
	Fic.Open("fftxlc.dat",CFile::modeWrite|CFile::modeCreate);
	for (l=0;l<HautIm;l++)
	{
		for (c=0;c<LargIm;c++)
		{
			ChaineTemp.Format("%f;",grid[0][l][c].x);
			Fic.WriteString(ChaineTemp);
		}
		Fic.WriteString("\n");
	}
	Fic.Close();

	Fic.Open("fftx1.dat",CFile::modeWrite|CFile::modeCreate);
	for (c=0;c<LargIm;c++)
	{
		ChaineTemp.Format("%f\n",grid[0][64][c].x);
		Fic.WriteString(ChaineTemp);
	}
	Fic.Close();
*/
	for(r=0;r<3;r++)
	{
		for(c=1;c<N;c++)
		{
			for(l=1;l<N;l++)
			{ 
				grid[r][c][l].x=((fct[c][l].x * grid[r][c][l].x) - (fct[c][l].y * grid[r][c][l].y));

				grid[r][c][l].y=((fct[c][l].x*grid[r][c][l].y)+(fct[c][l].y*grid[r][c][l].x));
			}
		}
	}

/*
	Fic.Open("fftxlcCor.dat",CFile::modeWrite|CFile::modeCreate);
	for (l=0;l<HautIm;l++)
	{
		for (c=0;c<LargIm;c++)
		{
			ChaineTemp.Format("%f;",grid[0][l][c].x);
			Fic.WriteString(ChaineTemp);
		}
		Fic.WriteString("\n");
	}
	Fic.Close();
	Fic.Open("fftylcCorIma.dat",CFile::modeWrite|CFile::modeCreate);
	for (l=0;l<HautIm;l++)
	{
		for (c=0;c<LargIm;c++)
		{
			ChaineTemp.Format("%f;",grid[0][l][c].y);
			Fic.WriteString(ChaineTemp);
		}
		Fic.WriteString("\n");
	}
	Fic.Close();
*/
	//RecentrageGrilleCouleur(grid, HautIm, LargIm);

//	CStdioFile Fic;
//	CString ChaineTemp;
/*	Fic.Open("ffty.dat",CFile::modeWrite|CFile::modeCreate);
	for (c=0;c<LargIm;c++)
	{
		ChaineTemp.Format("%f\n",grid[0][65][c].y);
		Fic.WriteString(ChaineTemp);
	}
	Fic.Close();

	Fic.Open("fftx.dat",CFile::modeWrite|CFile::modeCreate);
	for (c=0;c<LargIm;c++)
	{
		ChaineTemp.Format("%f\n",grid[0][65][c].x);
		Fic.WriteString(ChaineTemp);
	}
	Fic.Close();
*/
/*	double positionpix;
	double positionpix_x;
	double positionpix_y;
	double module;
	double moduleReel;
	int Ray;*/
	double MaxR,MaxV,MaxB;
	int rmax,lmax;	
	rmax = LargIm/2;
	lmax = HautIm/2;


	// pour visu imaginaire FFT

/*	MaxR=(grid[0][rmax][rmax].y);

	for(l = 0;l<HautIm ;l++)
	{
		for( c = 0;c<LargIm ;c++)
		{
			//module=_cabs(grid[0][c][l]);
			module=(grid[0][c][l].y);
			if (MaxR<module) MaxR=module;
		} // colonne
	} // ligne
	p=0;
	for(l = 0;l<HautIm ;l++)
	{
		for( c = 0;c<LargIm ;c++)
		{
			//module=(_cabs(grid[0][c][l])/MaxR)*10240;
			module=((grid[0][c][l].y)/MaxR)*10240;

			if (module>255.0) module=255;
			pImage[p].r=module;
			pImage[p].g=module;
			pImage[p].b=module;
			p++;
		} // colonne
	} // ligne

*/



	// modification de la calcul pour la nouvelle famille de Dscripteurs

	
	int t1=0;
	int t2=0,l3=0;

	
	_complex *coefR;
	_complex *coefV;
	_complex *coefB;

	
	_complex *coef1R;
	_complex *coef1V;
	_complex *coef1B;

	_complex *coef2R;
	_complex *coef2V;
	_complex *coef2B;

	_complex *coef3R;
	_complex *coef3V;
	_complex *coef3B;

	_complex *coefRf;
	_complex *coefVf;
	_complex *coefBf;

/*	_complex **tabR;
	_complex **tabV;
	_complex **tabB;
*/	
	
	coef = (double*)malloc(6*m_Npar*sizeof(double));

	//pour le résultat
	coefR = (_complex *)malloc((N2*lmax)*sizeof(_complex));
	coefV = (_complex*)malloc((N2*lmax)*sizeof(_complex));
	coefB = (_complex*)malloc((N2*lmax)*sizeof(_complex));

	//pour le résultat finale aprés d'éliminer les descripteurs qui ont une circonf < seuil
	coefRf = (_complex*)malloc((N2*lmax)*sizeof(_complex));
	coefVf = (_complex*)malloc((N2*lmax)*sizeof(_complex));
	coefBf = (_complex*)malloc((N2*lmax)*sizeof(_complex));

	//pour le vecteur V1
	coef1R = (_complex*)malloc((N2*lmax)*sizeof(_complex));
	coef1V = (_complex*)malloc((N2*lmax)*sizeof(_complex));
	coef1B = (_complex*)malloc((N2*lmax)*sizeof(_complex));

	//pour le vecteur V2
	coef2R = (_complex*)malloc((N2*lmax)*sizeof(_complex));
	coef2V = (_complex*)malloc((N2*lmax)*sizeof(_complex));
	coef2B = (_complex*)malloc((N2*lmax)*sizeof(_complex));

	//pour le vecteur V3
	coef3R = (_complex*)malloc((N2*lmax)*sizeof(_complex));
	coef3V = (_complex*)malloc((N2*lmax)*sizeof(_complex));
	coef3B = (_complex*)malloc((N2*lmax)*sizeof(_complex));

/*	tabR = (_complex **) malloc (N2*sizeof(_complex *));
	for(s=0;s<N2;s++)
		tabR[s]=	(_complex*)malloc(N2*sizeof(_complex));

	tabV = (_complex **) malloc (N2*sizeof(_complex *));

	for(s=0;s<N2;s++)
		tabV[s]=	(_complex*)malloc(N2*sizeof(_complex));

	tabB = (_complex **) malloc (N2*sizeof(_complex *));

	for(s=0;s<N2;s++)
		tabB[s]=	(_complex*)malloc(N2*sizeof(_complex));
	for(i=0;i<N2;i++)

	{
		for(j=0;j<N2;j++)
		{
			tabR[i][j].x=99.0;
			tabR[i][j].y=99.0;
			tabV[i][j].x=99.0;
			tabV[i][j].y=99.0;
			tabB[i][j].x=99.0;
			tabB[i][j].x=99.0;
		}
	}

	*/
	long *Circonf=(long * )malloc(N2*lmax*sizeof(long));


	for(i=0;i<(N2*lmax);i++)
	{
		coefR[i].x=0.0;
		coefR[i].y=0.0;
		coefV[i].x=0.0;
		coefV[i].y=0.0;
		coefB[i].x=0.0;
		coefB[i].y=0.0;

		coef1R[i].x=0.0;
		coef1R[i].y=0.0;
		coef1V[i].x=0.0;
		coef1V[i].y=0.0;
		coef1B[i].x=0.0;
		coef1B[i].x=0.0;

		coef2R[i].x=0.0;
		coef2R[i].y=0.0;
		coef2V[i].x=0.0;
		coef2V[i].y=0.0;
		coef2B[i].x=0.0;
		coef2B[i].y=0.0;

		coef3R[i].x=0.0;
		coef3R[i].y=0.0;
		coef3V[i].x=0.0;
		coef3V[i].y=0.0;
		coef3B[i].x=0.0;
		coef3B[i].y=0.0;

		coefRf[i].x=0.0;
		coefRf[i].y=0.0;
		coefVf[i].x=0.0;
		coefVf[i].y=0.0;
		coefBf[i].x=0.0;
		coefBf[i].y=0.0;
		
		
		Circonf[i]=0;
	}

	MaxR=0.0;
	MaxV=0.0;
	MaxB=0.0;

		
			
	


int l1,l2,z=0,teta,w,NL=0,NC=0;
int lamda1,lamda2;
lamda1=1;
lamda2=32;
BOOL Test;
Test=true;
	//CStdioFile Fic;
	//CString ChaineTemp;


//Fic.Open("dl3.dat",CFile::modeWrite|CFile::modeCreate);
	
	//for(w=0;w<360;w=w+10)
//	{
//	if(Test) z++;
//	Test=false;
w=10;
l1=1;
//l2=13;

	double wRd=(double)(w*pi)/180.0;
	
	//for(l1=lamda1;l1<32;l1=l1+4)

	//{	
		//NC=0;
		//for(l2=lamda2;l2<64;l2=l2++)

		//{ //NC++;
//for (w=0;w<=40;w=w+20)
//{
	for (l2=1;l2<64;l2++)
{

		for(teta=0;teta<360;teta=teta++)

	{ 
		double tetaRd=(double)(teta*pi)/180.0; 
		double wtRd=(double)((teta+w)*pi)/180.0;


		/*
		double dV1x= (double)(-l1*sin(wtRd));
		double dV1y= (double)l1*cos(wtRd);//wtRd=tRd+wRd

		double dV2x= (double)(-l2*sin(tetaRd)); 
		double dV2y= (double)l2*cos(tetaRd);
		*/
		

		
	/*	double dV1x= (double)l1*cos(wtRd);
		double dV1y= (double)l1*sin(wtRd);//wtRd=tRd+wRd

		double dV2x= (double)l2*cos(tetaRd); 
		double dV2y= (double)l2*sin(tetaRd);
*/
		double dV2x= (double)l2*cos(wtRd);
		double dV2y= (double)l2*sin(wtRd);//wtRd=tRd+wRd

		double dV1x= (double)l1*cos(tetaRd); 
		double dV1y= (double)l1*sin(tetaRd);
		
		// Vecteur V1

			// Vecteur V3

		double dV3x=dV1x+dV2x;
		double dV3y=dV1y+dV2y;

		double dl3=sqrt((dV3x*dV3x)+(dV3y*dV3y));

		//ChaineTemp.Format("%f\n",dl3);
		//Fic.WriteString(ChaineTemp);
		
	/*
		int V1x=ceil(dV1x);
		int V1y=ceil(dV1y);

		int V2x=ceil(dV2x);
		int V2y=ceil(dV2y);
		
		int V3x=ceil(dV3x);
		int V3y=ceil(dV3y);
	*/
		int l3=ceil(dl3);


		double V1x=dV1x+N2;
		double V1y=dV1y+N2;
		double V2x=dV2x+N2;
		double V2y=dV2y+N2;
		double V3x=dV3x+N2;
		double V3y=dV3y+N2;

			int V1xp;
			if(V1x>N2)
				V1xp=floor(V1x);
			else 
				 V1xp=ceil(V1x);

			int V1yp;

			if(V1y>N2)
				V1yp=floor(V1y);
			else 
				 V1yp=ceil(V1y);

			int V2xp;
			if(V2x>N2)
				V2xp=floor(V2x);
			else 
				 V2xp=ceil(V2x);

			int V2yp;

			if(V2y>N2)
				V2yp=floor(V2y);
			else 
				 V2yp=ceil(V2y);
			
			int V3xp;
			if(V3x>N2)
				V3xp=floor(V3x);
			else 
				 V3xp=ceil(V3x);

			int V3yp;

			if(V3y>N2)
				V3yp=floor(V3y);
			else 
				 V3yp=ceil(V3y);

			

		if(l3<N2)
		{
			
			//Test=true;
		/*
			int V1xp=V1x+64;
			int V1yp=V1y+64;

			int V2xp=V2x+64;
			int V2yp=V2y+64;

			int V3xp=V3x+64;
			int V3yp=V3y+64;
			*/

			//transforméée de fourier au vecteur V1

					coef1R[z].x= grid[0][V1xp][V1yp].x;
					coef1R[z].y= grid[0][V1xp][V1yp].y;
					coef1V[z].x= grid[1][V1xp][V1yp].x;
					coef1V[z].y= grid[1][V1xp][V1yp].y;
					coef1B[z].x= grid[2][V1xp][V1yp].x;
					coef1B[z].y= grid[2][V1xp][V1yp].y;

			//transforméée de fourier au vecteur V2

					coef2R[z].x= grid[0][V2xp][V2yp].x;
					coef2R[z].y= grid[0][V2xp][V2yp].y;
					coef2V[z].x= grid[1][V2xp][V2yp].x; 
					coef2V[z].y= grid[1][V2xp][V2yp].y;
					coef2B[z].x= grid[2][V2xp][V2yp].x;
					coef2B[z].y= grid[2][V2xp][V2yp].y;

			// transformée de Fourier au vecteur V3

					coef3R[z].x=	grid[0][V3xp][V3yp].x;
					coef3R[z].y=	(-grid[0][V3xp][V3yp].y);
					coef3V[z].x=	grid[1][V3xp][V3yp].x;
					coef3V[z].y=	(-grid[1][V3xp][V3yp].y);
					coef3B[z].x=	grid[2][V3xp][V3yp].x;
					coef3B[z].y=	(-grid[2][V3xp][V3yp].y);
					
			// transforméée de Fourier au produit V1*V2*(conjugué V3)
					
			
					coefR[z].x+=((coef1R[z].x*coef2R[z].x*coef3R[z].x)-(coef1R[z].x*coef2R[z].y*coef3R[z].y)-(coef1R[z].y*coef2R[z].x*coef3R[z].y)-(coef1R[z].y*coef2R[z].y*coef3R[z].x));
					coefR[z].y+=((coef1R[z].x*coef2R[z].x*coef3R[z].y)+(coef1R[z].x*coef2R[z].y*coef3R[z].x)+(coef1R[z].y*coef2R[z].x*coef3R[z].x)-(coef1R[z].y*coef2R[z].y*coef3R[z].y));

					coefV[z].x+=((coef1V[z].x*coef2V[z].x*coef3V[z].x)-(coef1V[z].x*coef2V[z].y*coef3V[z].y)-(coef1V[z].y*coef2V[z].x*coef3V[z].y)-(coef1V[z].y*coef2V[z].y*coef3V[z].x));
					coefV[z].y+=((coef1V[z].x*coef2V[z].x*coef3V[z].y)+(coef1V[z].x*coef2V[z].y*coef3V[z].x)+(coef1V[z].y*coef2V[z].x*coef3V[z].x)-(coef1V[z].y*coef2V[z].y*coef3V[z].y));
					
					coefB[z].x+=((coef1B[z].x*coef2B[z].x*coef3B[z].x)-(coef1B[z].x*coef2B[z].y*coef3B[z].y)-(coef1B[z].y*coef2B[z].x*coef3B[z].y)-(coef1B[z].y*coef2B[z].y*coef3B[z].x));
					coefB[z].y+=((coef1B[z].x*coef2B[z].x*coef3B[z].y)+(coef1B[z].x*coef2B[z].y*coef3B[z].x)+(coef1B[z].y*coef2B[z].x*coef3B[z].x)-(coef1B[z].y*coef2B[z].y*coef3B[z].y));
						
					Circonf[z]++;
				
			;
		}//if
	
	
}//teta
		
	
			if (fabs(coefR[z].x )> MaxR) MaxR=fabs(coefR[z].x);
			if (fabs(coefV[z].x )> MaxV) MaxV=fabs(coefV[z].x);	
			if (fabs(coefB[z].x )> MaxB) MaxB=fabs(coefB[z].x);
	
	z++;
}//fin l2
		
//}//fin w

//	tabR[NL][NC].x=coefR[z].x;


	//NC++;
//	z++;

	// recherche de Maximum pour lapartie réelle 
			//if (fabs(coefR[z].x )> MaxR) MaxR=fabs(coefR[z].x);
			//if (fabs(coefV[z].x )> MaxV) MaxV=fabs(coefV[z].x);	
			//if (fabs(coefB[z].x )> MaxB) MaxB=fabs(coefB[z].x);

		// recherche de Maximum pour lapartie imaginaire
			//if (fabs(coefR[z].y )> MaxR) MaxR=fabs(coefR[z].y);
			//if (fabs(coefV[z].y) > MaxV) MaxV=fabs(coefV[z].y);	
			//if (fabs(coefB[z].y) > MaxB) MaxB=fabs(coefB[z].y);

	//z++;

	//}//l2
	//NL++;
	
	//}//l1
	
	//}//w
/*	
for (l=0;l<NL;l++)
	{
		for (c=0;c<NC;c++)
		{
			tabR[l][c].x=tabR[l][c].x/MaxR;
		}
}
	Fic.Open("coefLC1.dat",CFile::modeWrite|CFile::modeCreate);
	for (l=0;l<NL;l++)
	{
		for (c=0;c<NC;c++)
		{
			ChaineTemp.Format("%f;",tabR[l][c].x);
			Fic.WriteString(ChaineTemp);
		}
		Fic.WriteString("\n");
	}
	Fic.Close();
 */

 
//	Fic.Close();
/*
	Fic.Open("ree1.dat",CFile::modeWrite|CFile::modeCreate);
	for (c=0;c<=z;c++)
	{
		ChaineTemp.Format("%f\n",coefR[c].x);
		Fic.WriteString(ChaineTemp);
	}
	Fic.Close();

	Fic.Open("ime1.dat",CFile::modeWrite|CFile::modeCreate);
	for (c=0;c<=z;c++)
	{
		ChaineTemp.Format("%f\n",coefR[c].y);
		Fic.WriteString(ChaineTemp);
	}
	Fic.Close();
*/


int t=0;
	
	// éliminer les descripteurs qui ont un Circonf faibles
	for(i=0;i<z;i++)
	{
		if(Circonf[i]>15)

		{ 
			coefRf[t].x=coefR[i].x;
			coefVf[t].x=coefV[i].x;
			coefBf[t].x=coefB[i].x;
			
			t++;
		 }
	}
  
	Nb_desc = t;
	// Normalisation
	
	
	for( i=0;  i<t; i++)
	{
		if(MaxR!=0)
		coefRf[i].x = coefRf[i].x/MaxR;//(MaxR*Circonf[i]);
		if(MaxV!=0)
		coefVf[i].x = coefVf[i].x/MaxV;//(MaxV*Circonf[i]);
		if(MaxB!=0)
		coefBf[i].x = coefBf[i].x/MaxB;//(MaxB*Circonf[i]);
	}
	
	// de meme recopie sans prendre le premier de chaque composante



	for( i=0; i<t; i++)
	{
		coef[i]=(coefRf[i+1].x ); 
		coef[i+t]=(coefVf[i+1].x );
		coef[i+2*t]=(coefBf[i+1].x );
	}

	for( i=0; i<3*t; i++)
	{
		coef[i]=fabs(coef[i]);
	}
	/*
	Fic.Open("Robjet1.dat",CFile::modeWrite|CFile::modeCreate);
	for (c=0;c<3*t;c++)
	{
		ChaineTemp.Format("%f\n",coef[c]);
		Fic.WriteString(ChaineTemp);
	}
	Fic.Close();
*/
/*	// pour la partie imaginaire
	for(i=0;i<z;i++)
	{
		if(Circonf[i]>15)

		{ 
			coefRf[t].y=coefR[i].y;
			coefVf[t].y=coefV[i].y;
			coefBf[t].y=coefB[i].y;
			t++;
		 }
	}
  
	Nb_desc = t;
	// Normalisation
	
	for( i=0;  i<t; i++)
	{
		if(MaxR!=0)
		coefRf[i].y = coefRf[i].y/MaxR;//(MaxR*Circonf[i]);
		if(MaxV!=0)
		coefVf[i].y = coefVf[i].y/MaxV;//(MaxV*Circonf[i]);
		if(MaxB!=0)
		coefBf[i].y = coefBf[i].y/MaxB;//(MaxB*Circonf[i]);
	}
	
	
  // de meme recopie sans prendre le premier de chaque composante

	for( i=0; i<t; i++)
	{
		coef[i]=(coefRf[i+1].y+1)/2; //i+1
		coef[i+t]=(coefVf[i+1].y+1)/2;
		coef[i+2*t]=(coefBf[i+1].y+1)/2;
	}

	Fic.Open("im.dat",CFile::modeWrite|CFile::modeCreate);
	for (c=0;c<t;c++)
	{
		ChaineTemp.Format("%f\n",coef[c]);
		Fic.WriteString(ChaineTemp);
	}
	Fic.Close();
*/
	free(coefR);
	free(coefV);
	free(coefB);

	free(coef1R);
	free(coef1V);
	free(coef1B);

	free(coef2R);
	free(coef2V);
	free(coef2B);

	free(coef3R);
	free(coef3V);
	free(coef3B);

	free(coefRf);
	free(coefVf);
	free(coefBf);



	for(m=0;m<3;m++)
	{
		for(j=0;j<LargIm;j++) 
		{
			free(grid[m][j]);
			//free(grid[m][j]);
		}
		free(grid[m]);
		//free(grid1[m]);
	}
	free(grid);
	//free(grid1);

	free(real);
	free(imag); 
	free(Circonf); 
	
		return coef;

	// fin de modification
	

/*
//	coef = new double[HautIm/2];
	coef = (double*)malloc(m_Npar*sizeof(double));
	coefR = (double*)malloc((lmax)*sizeof(double));
	coefV = (double*)malloc((lmax)*sizeof(double));
	coefB = (double*)malloc((lmax)*sizeof(double));
	long *Circonf=(long * )malloc(lmax*sizeof(long));

	for( i=0;  i<lmax; i++)
	{
		coefR[i] = 0.0;
		coefV[i] = 0.0;
		coefB[i] = 0.0;
		Circonf[i]=0;
	}
	for(l = 0;l<HautIm ;l++)
	{
		for( c = 0;c<LargIm ;c++)
		{
			positionpix_x = c - rmax;
			positionpix_y = l - rmax;

			positionpix_x*=positionpix_x;
			positionpix_y*=positionpix_y;
			positionpix=(int)sqrt(positionpix_x+positionpix_y);
			if(positionpix<rmax && positionpix>0)
			{
				Ray=(int)positionpix;
				coefR[Ray] += _cabs(grid[0][c][l]);
				coefV[Ray] += _cabs(grid[1][c][l]);
				coefB[Ray] += _cabs(grid[2][c][l]);
				Circonf[Ray]++;
			}
		
		} // colonne
	} // ligne
	coefR[0]=_cabs(grid[0][rmax][rmax]);
	coefV[0]=_cabs(grid[1][rmax][rmax]);
	coefB[0]=_cabs(grid[2][rmax][rmax]);
	MaxR=coefR[0];

	MaxV=coefV[0];
	MaxB=coefB[0];
	for( i=1;  i<lmax; i++)
	{
		coefR[i] = coefR[i]/(MaxR*Circonf[i]);
		coefV[i] = coefV[i]/(MaxV*Circonf[i]);
		coefB[i] = coefB[i]/(MaxB*Circonf[i]);
	}


   // recopie sans prendre le premier de chaque composante
	int np=m_Npar/3;
	for( i=0; i<np; i++)
	{
		coef[i]=coefR[i+1]; 
		coef[i+np]=coefV[i+1];
		coef[i+2*np]=coefB[i+1];
	}

	free(coefR);
	free(coefV);
	free(coefB);



	for(m=0;m<3;m++)
	{
		for(j=0;j<LargIm;j++) 
			free(grid[m][j]);
		free(grid[m]);
	}
	free(grid);

	free(real);
	free(imag);
	free(Circonf);

//	free(TabData);
//	free(fft);
	return coef;
*/
}

void CIma::RecentrageGrilleCouleur(_complex ***grid, long HautIm, long LargIm)
{
	long i,j,irc,jrc;
	_complex tft;

	for(int r=0; r<3; r++)
	{
		for(j=0;j<HautIm/2;j++)
		{
			jrc=(j+HautIm/2);//%HautIm;
			for(i=0;i<LargIm;i++)
			{
				tft=grid[r][j][i];
				grid[r][j][i]=grid[r][jrc][i];
				grid[r][jrc][i]=tft;
			}
		}
		for(i=0;i<LargIm/2;i++)
		{
			irc=(i+LargIm/2);//%LargIm;
			for(j=0;j<HautIm;j++)
			{
				tft=grid[r][j][i];
				grid[r][j][i]=grid[r][j][irc];
				grid[r][j][irc]=tft;
			} 
		}
	}
}

void CIma::LitMinMax(CString NomFichier)
{
	CFile fi;
	if(fi.Open(NomFichier,CFile::modeRead))
	{
		fi.Read(&m_Npar,sizeof(m_Npar));
		if (Tmin!=NULL) free(Tmin);
		if (Tmax!=NULL) free(Tmax);
		Tmin=(float *) malloc(m_Npar*sizeof(float));
		Tmax=(float *) malloc(m_Npar*sizeof(float));
		fi.Read(Tmin,m_Npar*sizeof(float));
		fi.Read(Tmax,m_Npar*sizeof(float));
		fi.Close();
	}

}

void CIma::NormaliseAvecTable(float *tab, int Taille)
{
	int i;
	float coef;
	float val;
	float binf=0.1F;
	float bsup=0.9F;
	if (bsup<=binf) {return;}
	for(i=0;i<Taille;i++) 
		{
			coef=(Tmax[i]-Tmin[i])/(bsup-binf);
			val=tab[i];
			if (_isnan(val) || _isnan(-val) ) 
			{val=0; MessageBox(NULL,"attention","infini detect",MB_OK);}
	//		if (val==0.0) 
	//		{ ch.Format("%d",i);MessageBox(NULL,"attention",ch,MB_OK);}
			if (val<Tmin[i]) val=Tmin[i];
			else if (val>Tmax[i]) val=Tmax[i];
			else 
			{
				if (coef!=0) val=binf+((val-Tmin[i])/coef);
				else val=bsup;
			}
			tab[i]=val;
		}


} // fin normalise

void CIma::ChargeModelSVM(CString ChaineModel)
{
	if (model!=NULL)svm_destroy_model(model);
	model=svm_load_model(ChaineModel);



}
/**	\brief	The CIma::LireCmc fonction membre permettant de relire un fichier CMC existant afin de retrouver les différent fichier "*.dat"
	\author J Miteran
	\param	NomCmc	 CString contenant le nom du fichier CMC à relire

	\return	void

	
*/

void CIma::LireCmc(CString NomCmc)
{

	int i,NBclass,*NBFichParclass=NULL;
	CStdioFile Fichier; // variable fichier
	CString text,Temp;
	if (Fichier.Open(NomCmc, CFile::modeRead))
	{ 

		TabRefNom.RemoveAll();
		Fichier.ReadString(text);
		NBclass=atoi(LPCSTR(text));
		NBFichParclass=new int[NBclass];
		for( i=0;i<NBclass;i++)
		{
			Fichier.ReadString(text);
			NBFichParclass[i]=atoi(LPCSTR(text));
		}
		
		for(i=0;i<NBclass;i++)
		{
			for(int j=0;j<NBFichParclass[i];j++)
			{
				Fichier.ReadString(text);		
			}
			Temp=text.Right(text.GetLength()-text.ReverseFind('\\')-1);
			TabRefNom.Add(Temp.Left(Temp.GetLength()-4));
		}
		if (NBFichParclass!=NULL) delete NBFichParclass;
	} // ouverture ok

}

void CIma::LitBMP(CString Nom)
{
	CFile Fichier;
	long Taille;
	if (Fichier.Open(Nom, CFile::modeRead))
	{   
		if (m_PImage!=NULL) free(m_PImage);
		if (m_Info!=NULL) delete m_Info;
		if (m_Tete!=NULL) delete m_Tete;
	    m_Tete = (BITMAPFILEHEADER*) new BYTE[sizeof(BITMAPFILEHEADER)];    
		m_Info = (BITMAPINFO*) new BYTE[ sizeof(BITMAPINFOHEADER) /*+ sizeof(RGBQUAD)*256*/];
		Fichier.Read(m_Tete,sizeof(BITMAPFILEHEADER));
		if (m_Tete->bfType==19778)
		{
			Fichier.Read(m_Info,sizeof(BITMAPINFOHEADER) /*+ sizeof(RGBQUAD)*256*/);
			m_Largeur = m_Info->bmiHeader.biWidth;    
			m_Hauteur =m_Info->bmiHeader.biHeight;
			Taille=m_Largeur*m_Hauteur*sizeof(rgb);
			m_PImage = (rgb*) malloc(Taille);
			Fichier.Read(m_PImage,Taille);
		}	 
		Fichier.Close();
	}	

}

void CIma::VisuBMP(CClientDC *pdc)
{
	
	if (m_PImage!=NULL) 
		SetDIBitsToDevice(pdc->GetSafeHdc(),330,0, m_Largeur, m_Hauteur , 
	0, 0 ,0,m_Hauteur,m_PImage,m_Info,DIB_RGB_COLORS);
}

void CIma::SauveBMP(CString NomFich)
{

	CFile Fichier; // variable fichier BITMAP
	if (Fichier.Open(NomFich,CFile::modeWrite|CFile::modeCreate))
	  {   
		m_Tete->bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER); 
		m_Tete->bfSize = m_Tete->bfOffBits+m_Largeur*m_Hauteur*sizeof(rgb);
		m_Tete->bfType = 'MB'; // obligatoire
		m_Tete->bfReserved1 = 0;
		m_Tete->bfReserved2 = 0;  
		Fichier.Write(m_Tete,(UINT) sizeof(BITMAPFILEHEADER)); // entête du fichier
		Fichier.Write(m_Info,sizeof(BITMAPINFOHEADER));   // entête de l'image
		Fichier.Write(m_PImage,m_Largeur*m_Hauteur*sizeof(rgb));	 // écriture des données
		Fichier.Close();

	}// si ouverture du fichier ok
}

void CIma::ReSampleCarreBMP()
{
	int p,pp,cp,lp,l,c;
	float RatioL,RatioC;
	int Taille=m_Largeur*m_Hauteur;
	RatioL=(float)m_Hauteur/(float)LFinal;
	RatioC=(float)m_Largeur/(float)CFinal;
//	rgb* PImageCarre;
	
	if (m_PImage==NULL) return;

//	PImageCarre = (rgb*) malloc(CFinal*LFinal*sizeof(rgb)); // mémoire pour les données

	for (c=0;c<CFinal;c++)
	{
		cp=c*RatioC;
		for (l=0;l<LFinal;l++)
		{
			lp=l*RatioL;
			p=LFinal*l+c;
			pp=m_Largeur*lp+cp;
			if (pp>Taille) pp=0;
			PImageCarre[p].r=m_PImage[pp].r;
			PImageCarre[p].g=m_PImage[pp].g;
			PImageCarre[p].b=m_PImage[pp].b;
		} // for l
	} // for c
	
	m_Largeur=CFinal;
	m_Hauteur=LFinal;
	free(m_PImage);
	m_PImage=(rgb*)malloc(m_Largeur*m_Hauteur*sizeof(rgb));
	m_Info->bmiHeader.biWidth=m_Largeur;    
	m_Info->bmiHeader.biHeight=m_Hauteur;
	memcpy(m_PImage,PImageCarre,CFinal*LFinal*sizeof(rgb));
//	free(PImageCarre);

}

bool CIma::FFT(int dir, int m, double *x, double *y)
{
	long nn,i,i1,j,k,i2,l,l1,l2;
	double c1,c2,tx,ty,t1,t2,u1,u2,z;

	/* Calculate the number of points */
	nn = 1<<m;

	/* Do the bit reversal */
	i2 = nn >> 1;
	j = 0;
	for (i=0;i<nn-1;i++) {
		if (i < j) {
			tx = x[i];
			ty = y[i];
			x[i] = x[j];
			y[i] = y[j];
			x[j] = tx;
			y[j] = ty;
		}
		k = i2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}

	/* Compute the FFT */
	c1 = -1.0;
	c2 = 0.0;
	l2 = 1;
	for (l=0;l<m;l++) {
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0;
		u2 = 0.0;
		for (j=0;j<l1;j++) {
			for (i=j;i<nn;i+=l2) {
				i1 = i + l1;
				t1 = u1 * x[i1] - u2 * y[i1];
				t2 = u1 * y[i1] + u2 * x[i1];
				x[i1] = x[i] - t1;
				y[i1] = y[i] - t2;
				x[i] += t1;
				y[i] += t2;
			}
			z =  u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		c2 = sqrt((1.0 - c1) / 2.0);
		if (dir == 1)
			c2 = -c2;
		c1 = sqrt((1.0 + c1) / 2.0);
	}

	/* Scaling for forward transform */
	//if (dir == 1) {
	//	for (i=0;i<nn;i++) {
	//		x[i] /= (double)nn;
	//		y[i] /= (double)nn;
	//	}
	//}

   return true;
}

void CIma::ReSampleCarreBMPOffset()
{
	int p,pp,cp,lp,l,c;
	float RatioL,RatioC;
	int Taille=m_Largeur*m_Hauteur;
	BYTE MonSeuil=50;
	if (m_PImage==NULL) return;

	int CInitial=0;
	int LInitial=0;

	for (l=0;l<m_Hauteur;l++)
	{
		for (c=0;c<m_Largeur;c++)
		{
			pp=m_Largeur*l+c;
			if (m_PImage[pp].b>MonSeuil) 
			{
				LInitial=l-10;
				break;
			}
		} // for l
		if (LInitial!=0) break;
	} // for c
	for (c=0;c<m_Largeur;c++)
	{
		for (l=0;l<m_Hauteur;l++)
		{
			pp=m_Largeur*l+c;
			if (m_PImage[pp].b>MonSeuil) 
			{
				CInitial=c-10;
				break;
			}
		} // for l
		if (CInitial!=0) break;
	} // for c

/*	for (l=0;l<m_Hauteur;l++)
	{
			p=m_Largeur*l+CInitial;
			m_PImage[p].r=255;
			m_PImage[p].g=255;
			m_PImage[p].b=255;
	} // for l



	for (l=0;l<280;l++)
	{
		for (c=0;c<280;c++)
		{
			pp=m_Largeur*(l+LInitial)+c+CInitial;
			p=m_Largeur*l+c;
			m_PImage[p].r=m_PImage[pp].r;
			m_PImage[p].g=m_PImage[pp].g;
			m_PImage[p].b=m_PImage[pp].b;
		} // for l
	} // for c

*/

	RatioL=(float)280/(float)LFinal;
	RatioC=(float)280/(float)CFinal;

	for (c=0;c<CFinal;c++)
	{
		cp=c*RatioC;
		for (l=0;l<LFinal;l++)
		{
			lp=l*RatioL;
			p=LFinal*l+c;
			pp=m_Largeur*(lp+LInitial)+(cp+CInitial);
			if (pp>Taille) pp=0;
			PImageCarre[p].r=m_PImage[pp].r;
			PImageCarre[p].g=m_PImage[pp].g;
			PImageCarre[p].b=m_PImage[pp].b;
		} // for l
	} // for c
	
	m_Largeur=CFinal;
	m_Hauteur=LFinal;
	free(m_PImage);
	m_PImage=(rgb*)malloc(m_Largeur*m_Hauteur*sizeof(rgb));
	m_Info->bmiHeader.biWidth=m_Largeur;    
	m_Info->bmiHeader.biHeight=m_Hauteur;
	m_Info->bmiHeader.biSizeImage = sizeof(rgb)*LFinal*CFinal;
	memcpy(m_PImage,PImageCarre,CFinal*LFinal*sizeof(rgb));
	


}


void CIma::ReCherche(BYTE *Image, int Larg, int Haut)
{
	
	int p,p1,p0,co,li,l,c,lidef2,coldef2;
//	int i;
	float Ratio;
	int Taille=Larg*Haut*3;
	Ratio=(float)LCarre/(float)LFinal;
	int Ech,EchMin=4,EchMax=5;
	struct svm_node *x;
	double *DesFour;
	int par;
	double v;
	int cptmask=0;
	BYTE *ImageC;
	ImageC=(BYTE*) malloc(Larg*Haut*3);
	memcpy(ImageC,Image,Larg*Haut*3);


	LFinal=CFinal=TailleMasque;
	m_Npar=((TailleMasque/2)-1)*3;

	BITMAPINFO	 Info; // entête de l'image
	Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	Info.bmiHeader.biWidth = CFinal;
	Info.bmiHeader.biHeight = LFinal;
	Info.bmiHeader.biPlanes = 1;   
	Info.bmiHeader.biSizeImage = sizeof(rgb)*CFinal*LFinal;
	Info.bmiHeader.biBitCount = 24;
	Info.bmiHeader.biClrUsed = 0;
	Info.bmiHeader.biCompression =  BI_RGB;
	Info.bmiHeader.biXPelsPerMeter = 0;
	Info.bmiHeader.biYPelsPerMeter = 0;
	Info.bmiHeader.biClrImportant = 0;
	BITMAPFILEHEADER Tete;
	Tete.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER); 
	Tete.bfSize = Tete.bfOffBits+CFinal*LFinal*sizeof(rgb);
	Tete.bfType = 'MB'; // obligatoire
	Tete.bfReserved1 = 0;
	Tete.bfReserved2 = 0;  

	
	
	
	x = (struct svm_node *) malloc(m_Npar*sizeof(struct svm_node)*2);
	for(Ech=EchMin;Ech<EchMax;Ech++) // 3 a 6
	{
		for (l=2;l<(Haut-TailleMasque*Ech);l+=8)
		{
			for (c=2;c<(Larg-TailleMasque*Ech);c+=8)
			{

				for (li=0;li<TailleMasque;li++)
				{
					for (co=0;co<TailleMasque;co++)
					{
					
						p0=(TailleMasque*li+co);
						p1=(Larg*(l+(li*Ech))+c+(co*Ech))*3;
						PImageCarre[p0].r=ImageC[p1];
						PImageCarre[p0].g=ImageC[p1+1];
						PImageCarre[p0].b=ImageC[p1+2];
					} // co
				} // li

/*
			CFile Fichier; 
			CString NomFich;
			NomFich.Format("temp%d.bmp",cptmask++);
			if (Fichier.Open(NomFich,CFile::modeWrite|CFile::modeCreate))
			  {   
				Fichier.Write(&Tete,(UINT) sizeof(BITMAPFILEHEADER)); // entête du fichier
				Fichier.Write(&Info,sizeof(BITMAPINFOHEADER));   // entête de l'image
				Fichier.Write(PImageCarre,CFinal*LFinal*sizeof(rgb));	 // écriture des données
				Fichier.Close();

			}// si ouverture du fichier ok
*/	




				DesFour=DescripteursFourier1(PImageCarre);
			
			
				for(par=0;par<m_Npar;par++)
				{ 
					x[par].index=par+1;
					x[par].value=DesFour[par];

				}//for par
				x[m_Npar].index = -1;
				v = svm_predict(model,x);
				if (v==1.0) 
				{

					lidef2=l+(TailleMasque*Ech);
					coldef2=c+(TailleMasque*Ech);
					if (lidef2>Haut) lidef2 =Haut;
					if (coldef2>Larg) coldef2=Larg;
	
					for (li=l;li<lidef2;li++) 	
					{
						p=(li*Larg+c)*3;
						Image[p]=0;		Image[p+1]=0;		Image[p+2]=255;
						p=(li*Larg+c+(TailleMasque*Ech))*3;	
						Image[p]=0;		Image[p+1]=0;		Image[p+2]=255;

					}
					for (co=c;co<coldef2;co++) 
					{
						p=	(l*Larg+co)*3;
						Image[p]=0;		Image[p+1]=0;		Image[p+2]=255;
						p=	((l+TailleMasque*Ech)*Larg+co)*3;
						Image[p]=0;		Image[p+1]=0;		Image[p+2]=255;
					}

				}
				if (DesFour!=NULL) free(DesFour);

			} //c
		} // l
	} // Ech


	free(x);
	free(ImageC);



	
}

double* CIma::MomentsZernike(rgb *pImage)
{


	double *coef;


	int p,i,j;
	double i_0, j_0;
	double i_scale, j_scale;
	double x,y;
	double coef0R,coef0V,coef0B;
	_complex resR;
	_complex resV;
	_complex resB;
	_complex v;
	int isize=CFinal;
	int jsize=LFinal;
	int k,n,m;
	int OrdreMax=5;

	/*
	if (Ipixtyp(inband) != Iu_byte_typ)
	Error(1, "zer_mom: Input pixel type must be unsigned byte.\n");

	if ((n<0) || (abs(m) > n) || ((n-abs(m))%2!=0))
	Error(2, "zer_mom: n=%i, m=%i, n-|m|=%i\n", n, m, (n-abs(m)));
	*/

	i_0 = (isize+1)/2.0;
	j_0 = (jsize+1)/2.0;
	i_scale = isize/2.0;
	j_scale = jsize/2.0;

	m_Npar=0;
	for (n=0;n<5;n++)
	{
		for (m=-n; m<=n; m++)
		{
			if ((n-abs(m)) %2 == 0)
			{
				m_Npar++;
			}
		}
	}
	m_Npar--;
	int np=m_Npar; // pour une couleur
	m_Npar*=3; // pour les 3 couleurs
	coef = (double*)malloc(m_Npar*sizeof(double));

	// rouge
	k=0;
	for (n=0;n<OrdreMax;n++)
	{
		for (m=-n; m<=n; m++)
		{
			if ((n-abs(m)) %2 == 0)
			{

				resR.x = 0.0;
				resR.y = 0.0;
				resV.x = 0.0;
				resV.y = 0.0;
				resB.x = 0.0;
				resB.y = 0.0;

				for (i=0; i<isize; i++)
				{
					for (j=0; j< jsize; j++)
					  {
						x = (i-i_0)/i_scale;
						y = (j-j_0)/j_scale;
						/*	printf ("zer_mom: x=%6.3f y=%6.3f i=%i j=%i\n", x,y,i,j);
						*/
						if ((x*x + y*y) <= 1.0)
						{
							zer_pol(n,m,x,y, &v);
							p=CFinal*j+i;
							resR.x += pImage[p].r*v.x;
							resR.y += pImage[p].r*(-v.y);
							resV.x += pImage[p].g*v.x;
							resV.y += pImage[p].g*(-v.y);
							resB.x += pImage[p].b*v.x;
							resB.y += pImage[p].b*(-v.y);
						} // if 
					} // j
				} // i
				resR.x = resR.x*(n+1)/M_PI;
				resR.y = resR.y*(n+1)/M_PI;
				resV.x = resV.x*(n+1)/M_PI;
				resV.y = resV.y*(n+1)/M_PI;
				resB.x = resB.x*(n+1)/M_PI;
				resB.y = resB.y*(n+1)/M_PI;
				if (k!=0) 
				{
					coef[k-1]=_cabs(resR)/coef0R;
					coef[k-1+np]=_cabs(resV)/coef0V;
					coef[k-1+(2*np)]=_cabs(resB)/coef0B;
				}
				else 
				{
					coef0R=_cabs(resR);
					coef0V=_cabs(resV);
					coef0B=_cabs(resB);
				}
				k++;
			}// if
		}// for m
	} // for n
	

	return coef;

}

void CIma::zer_pol(int n, int m, double x, double y, _complex *res)
{
  double R;
  double arg;

  if ((x*x + y*y) > 1.0)
    {
      res->x = 0.0;
      res->y = 0.0;
    }
  else
    {
      zer_pol_R(n,m,x,y, &R);
      arg = m*atan2(y,x);
      res->x = R*cos(arg);
      res->y = R*sin(arg);
    }

}

void CIma::zer_pol_R(int n, int m_in, double x, double y, double *res)
{
  int i;
  int m;
  int s;
  int a; /* (n-s)! */
  int b; /*   s!   */
  int c; /* [(n+|m|)/2-s]! */
  int d; /* [(n-|m|)/2-s]! */
  int sign;

  m = abs(m_in);
/*
  if ((n-m)%2!=0)
    Warning(1,
	    "zer_pol_R: Warning. R(%i,%i,%5.2f,%5.2f). n-|m| = %i is odd!\n", 
	    n, m, x, y, m);
*/
  /* The code is optimized with respect to the faculty operations */

  res[0] = 0.0;
  if ((x*x + y*y) <= 1.0)
    {
      sign = 1;
      a = 1;
      for (i=2; i<=n; i++)
	a*=i;
      b=1;
      c = 1;
      for (i=2; i <= (n+m)/2; i++)
	c*=i;
      d = 1;
      for (i=2; i <= (n-m)/2; i++)
	d*=i;
      
      /* Before the loop is entered, all the integer variables (sign, a, */
      /* b, c, d) have their correct values for the s=0 case. */
      for (s=0; s<= (n-m)/2; s++)
	{
	  /*printf("zer_pol_R: s=%i, n=%i, m=%i, x=%6.3f, y=%6.3f, a=%i, */
	  /*b=%i, c=%i, d=%i, sign=%i\n", s,n,m,x,y,a,b,c,d,sign); */
	  res[0] += sign * (a*1.0/(b*c*d)) * pow((x*x + y*y),(n/2.0)-s);
	  
	  /* Now update the integer variables before the next iteration of */
	  /* the loop. */
	  
	  if (s < (n-m)/2)
	    {
	      sign = -sign;
	      a = a/(n-s);
	      b = b*(s+1);
	      c = c / ((n+m)/2 - s);
	      d = d / ((n-m)/2 - s);
	    }
	}      
    }

}

void CIma::ReSampleZernike(BYTE *Image, int Larg, int Haut, CClientDC *pdc, CString Nom, BOOL dec, BOOL Pos, BOOL Sauve, CString NomFich)
{
	int p,pp,cp,lp,l,c;
	float Ratio;
	int Taille=Larg*Haut*3;
	Ratio=(float)LCarre/(float)LFinal;

	if (Pos)
	{


		for (c=0;c<CFinal;c++)
		{
			cp=c*Ratio;
			for (l=0;l<LFinal;l++)
			{
				lp=l*Ratio;
				p=LFinal*l+c;
				pp=(Larg*lp+(cp+40))*3;
				if (pp>Taille) pp=0;
				PImageCarre[p].r=Image[pp];
				PImageCarre[p].g=Image[pp+1];
				PImageCarre[p].b=Image[pp+2];
			} // for l
		} // for c

		// pour visu carre


		// pour visu
		
		BITMAPINFO	 Info; // entête de l'image

		Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
		Info.bmiHeader.biWidth = CFinal;
		Info.bmiHeader.biHeight = LFinal;
		Info.bmiHeader.biPlanes = 1;   
		Info.bmiHeader.biSizeImage = sizeof(rgb)*CFinal*LFinal;
		Info.bmiHeader.biBitCount = 24;
		Info.bmiHeader.biClrUsed = 0;
		Info.bmiHeader.biCompression =  BI_RGB;
		Info.bmiHeader.biXPelsPerMeter = 0;
		Info.bmiHeader.biYPelsPerMeter = 0;
		Info.bmiHeader.biClrImportant = 0;
		if (Sauve)
		{
			CFile Fichier; 
			if (Fichier.Open(NomFich,CFile::modeWrite|CFile::modeCreate))
			  {   
				BITMAPFILEHEADER Tete;
				Tete.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER); 
				Tete.bfSize = Tete.bfOffBits+CFinal*LFinal*sizeof(rgb);
				Tete.bfType = 'MB'; // obligatoire
				Tete.bfReserved1 = 0;
				Tete.bfReserved2 = 0;  
				Fichier.Write(&Tete,(UINT) sizeof(BITMAPFILEHEADER)); // entête du fichier
				Fichier.Write(&Info,sizeof(BITMAPINFOHEADER));   // entête de l'image
				Fichier.Write(PImageCarre,CFinal*LFinal*sizeof(rgb));	 // écriture des données
				Fichier.Close();

			}// si ouverture du fichier ok
	
		}
// fin visu
		if (bContour) Sobel(PImageCarre,CFinal,LFinal);
		SetDIBitsToDevice(pdc->GetSafeHdc(),330,0, CFinal, LFinal , 0,0,0,LFinal,PImageCarre,&Info,DIB_RGB_COLORS);

		double *DesFour;
		
		DesFour=MomentsZernike(PImageCarre);
		int i;
		for (i=0;i<m_Npar;i++)
		{
			m_TabModule[i]=DesFour[i];
		}
		if (!dec) AjouterDonneeDansRef(Nom);
		else
		{ /*
			int par,cf;
			BYTE b;
			SBoostOutput Fin;
			if (pBoostMc==NULL) return;
			float *x_test=new float[m_Npar];
			for (par=0;par<m_Npar;par++) x_test[par]=m_TabModule[par+1];
			
			NormaliseAvecTable(x_test,m_Npar);

			Fin=pBoostMc->FinalClass(x_test);
			cf=Fin.cf;
			m_Chaine.Format("%d",cf); 
			delete [] x_test;*/

		
			struct svm_node *x;
			int par;
			double v;
			x = (struct svm_node *) malloc(m_Npar*sizeof(struct svm_node)*2);
			for(par=0;par<m_Npar;par++)
			{ 
				x[par].index=par+1;
				x[par].value=(double)m_TabModule[par];

			}//for par
			x[m_Npar].index = -1;
			v = svm_predict(model,x);
			//	v=1;
			m_Chaine=TabRefNom[(int)v-1];
			//m_Chaine.Format("%d",(int)v); 
			free(x);
		
		
		
		}

		// pour visu
	/*	
		pdc->Rectangle(0,500,400,280);
		pdc->MoveTo(0,500-DesFour[0]*200);
		
		for (i=0;i<m_Npar*4;i+=4)
		{
				pdc->LineTo(i,500-DesFour[i/4]*200);
		}
		
		// fin visu
*/
		if (DesFour!=NULL) free(DesFour);
	} // si position ok

	for (c=0;c<(Larg-CCarre)/2;c++)
	{
		for (l=0;l<Haut;l++)
		{
			pp=(Larg*l+c)*3;
			Image[pp]=127;
			Image[pp+1]=127;
			Image[pp+2]=127;
		} // for l
	} // for c
	for (c=Larg-(Larg-CCarre)/2;c<Larg;c++)
	{
		for (l=0;l<Haut;l++)
		{
			pp=(Larg*l+c)*3;
			Image[pp]=127;
			Image[pp+1]=127;
			Image[pp+2]=127;
		} // for l
	} // for c
}

/**	\brief	The CIma::ReSampleMixteZernikeFourier Fonction membre permettant de stocker les valeurs Des invariants de Zernike et des Descripteurs de Fouriers.

	\param	Image	un tableau BYTE conteant l'image
	\param	Larg	 int contenant la largeur de l'image
	\param	Haut	 int contenant la hauteur de l'image
	\param	pdc	CClientDC  pointeur sur le client Dc pour permettre l'affichage
	\param	Nom	CString contient le nom du fichier *.dat
	\param	dec	 BOOL bit de validation pour choix decision ou creation fichier d'apprentissage
	\param	Pos	 BOOL bit de validation pour choix création d'image pour l'affichage
	\param	Sauve BOOL bit de validation pour choix de sauvegarde de l'image
	\param	NomFich	CString contenant le nom du fichier bmp

	\return	void

	
*/


void CIma::ReSampleMixteZernikeFourier(BYTE *Image, int Larg, int Haut, CClientDC *pdc, CString Nom, BOOL dec, BOOL Pos, BOOL Sauve, CString NomFich)
{


	int p,pp,cp,lp,l,c;
	float Ratio;
	int Taille=Larg*Haut*3;
	Ratio=(float)LCarre/(float)LFinal;

	if (Pos)
	{


		for (c=0;c<CFinal;c++)
		{
			cp=c*Ratio;
			for (l=0;l<LFinal;l++)
			{
				lp=l*Ratio;
				p=LFinal*l+c;
				pp=(Larg*lp+(cp+40))*3;
				if (pp>Taille) pp=0;
				PImageCarre[p].r=Image[pp];
				PImageCarre[p].g=Image[pp+1];
				PImageCarre[p].b=Image[pp+2];
			} // for l
		} // for c

		// pour visu carre


		// pour visu
		
		BITMAPINFO	 Info; // entête de l'image

		Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
		Info.bmiHeader.biWidth = CFinal;
		Info.bmiHeader.biHeight = LFinal;
		Info.bmiHeader.biPlanes = 1;   
		Info.bmiHeader.biSizeImage = sizeof(rgb)*CFinal*LFinal;
		Info.bmiHeader.biBitCount = 24;
		Info.bmiHeader.biClrUsed = 0;
		Info.bmiHeader.biCompression =  BI_RGB;
		Info.bmiHeader.biXPelsPerMeter = 0;
		Info.bmiHeader.biYPelsPerMeter = 0;
		Info.bmiHeader.biClrImportant = 0;
		if (Sauve)
		{
			CFile Fichier; 
			if (Fichier.Open(NomFich,CFile::modeWrite|CFile::modeCreate))
			  {   
				BITMAPFILEHEADER Tete;
				Tete.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER); 
				Tete.bfSize = Tete.bfOffBits+CFinal*LFinal*sizeof(rgb);
				Tete.bfType = 'MB'; // obligatoire
				Tete.bfReserved1 = 0;
				Tete.bfReserved2 = 0;  
				Fichier.Write(&Tete,(UINT) sizeof(BITMAPFILEHEADER)); // entête du fichier
				Fichier.Write(&Info,sizeof(BITMAPINFOHEADER));   // entête de l'image
				Fichier.Write(PImageCarre,CFinal*LFinal*sizeof(rgb));	 // écriture des données
				Fichier.Close();

			}// si ouverture du fichier ok
	
		}
// fin visu
		// pretraitement eventuel
		if (bContour) Sobel(PImageCarre,CFinal,LFinal);
		// on affiche ce qui est traité

		SetDIBitsToDevice(pdc->GetSafeHdc(),330,0, CFinal, LFinal , 0,0,0,LFinal,PImageCarre,&Info,DIB_RGB_COLORS);

		/*double *DesFour;
		DesFour=MomentsZernike(PImageCarre);
		int i;
		for (i=0;i<m_Npar;i++)
		{
			m_TabModule[i]=DesFour[i];
		}
		if (!dec) AjouterDonneeDansRef(Nom);*/

		///////////////rajout debut

		

		
		long m_NparZer = 42;
		long m_NparFour = 189;

		  m_Npar=m_NparZer;

		  double *DesFour=MomentsZernike(PImageCarre);
		  int i;
		  for (i=0;i<m_NparZer;i++)
		  {
		   m_TabModule[i]=DesFour[i];
		  }
		  if (DesFour!=NULL) free(DesFour);
		  m_Npar=m_NparFour;
 
		  DesFour=DescripteursFourier1(PImageCarre);

		  for (i=0;i<m_NparFour;i++)
		  {
		   m_TabModule[i+m_NparZer]=DesFour[i];
		  }
		  m_Npar=m_NparZer+m_NparFour;
		  if (!dec) AjouterDonneeDansRef(Nom);


		////////////////rajout fin
		else
		{ /*
			int par,cf;
			BYTE b;
			SBoostOutput Fin;
			if (pBoostMc==NULL) return;
			float *x_test=new float[m_Npar];
			for (par=0;par<m_Npar;par++) x_test[par]=m_TabModule[par+1];
			
			NormaliseAvecTable(x_test,m_Npar);

			Fin=pBoostMc->FinalClass(x_test);
			cf=Fin.cf;
			m_Chaine.Format("%d",cf); 
			delete [] x_test;*/

		
			struct svm_node *x;
			int par;
			double v;
			x = (struct svm_node *) malloc(m_Npar*sizeof(struct svm_node)*2);
			for(par=0;par<m_Npar;par++)
			{ 
				x[par].index=par+1;
				x[par].value=(double)m_TabModule[par];

			}//for par
			x[m_Npar].index = -1;
			v = svm_predict(model,x);
			//	v=1;
			m_Chaine=TabRefNom[(int)v-1];
			//m_Chaine.Format("%d",(int)v); 
			free(x);
		
		
		
		}

		// pour visu
		
		pdc->Rectangle(0,500,400,280);
		pdc->MoveTo(0,500-DesFour[0]*200);
		
		for (i=0;i<m_Npar*4;i+=4)
		{
				pdc->LineTo(i,500-DesFour[i/4]*200);
		}
		
		// fin visu

		if (DesFour!=NULL) free(DesFour);
	} // si position ok

	for (c=0;c<(Larg-CCarre)/2;c++)
	{
		for (l=0;l<Haut;l++)
		{
			pp=(Larg*l+c)*3;
			Image[pp]=127;
			Image[pp+1]=127;
			Image[pp+2]=127;
		} // for l
	} // for c
	for (c=Larg-(Larg-CCarre)/2;c<Larg;c++)
	{
		for (l=0;l<Haut;l++)
		{
			pp=(Larg*l+c)*3;
			Image[pp]=127;
			Image[pp+1]=127;
			Image[pp+2]=127;
		} // for l
	} // for c

}


/**	\brief	CIma::DitherFloydSteinberg fonction membre permettant de de transformer une image couleur en une image binaire via un dither.
C'est ici l'algorithme de Floyd Steinberg (propagation d'erreur)

	\return	BYTE* tableau contenant les valeurs binaires.

	
*/

BYTE* CIma::DitherFloydSteinberg()
{
	// Temporary veriable
	int i, j,width,height;

	height= m_Hauteur;
	width = m_Largeur;

	long Taille = m_Hauteur*m_Largeur;
	BYTE * gray_data = (BYTE*)malloc(sizeof(BYTE)*m_Hauteur*m_Largeur);

	// Current error
	double error;
	// the current offset(x, y) and the offset of next row(x, y+1)
	int    now, next_row;
	// the starting point
//	int offset;

	// Transform the image to gray scale
	//To_Grayscale();
	j=0;
	for(i=0;i<Taille;i=i+3)
	{
	long temp = m_PImage[i].r+ m_PImage[i].g+m_PImage[i].b;
	temp = temp /3;
	gray_data[j] = temp;
	j++;
	}



	// Go through entire image with the zigzag form
	for(i = 0; i < height; i++){

		// For the odd number row 1, 3, 5, 7 and it should 
		// go from end to start
		if(i % 2){
			for(j = width - 1; j >= 0; j--){
				// Calculate the current offset and next row offest
				now = j + i * width;
				next_row = j + (i + 1) * width;

				// Determine whether the current pixel is black or white with
				// the threshold 0.5
				// black the error is current - 0;
				if(gray_data[now] < 0.5){
					error = gray_data[now];
					gray_data[now] = 0;
				}
				// white the error is current - 1.0
				else{
					error = gray_data[now] - 1.0;
					gray_data[now] = 1.0;
				}

				// Relay the error down to the following elements
				// The point is the end of the row
				if(j == 0){
					// The last row means the last element not need to relay
					// Not the last row needs to down, down right
					if(i != height - 1){
						gray_data[next_row]     += (error *  0.3125);
						gray_data[next_row + 1] += (error *  0.1875);
					}
				}// end of the last element of the row
				else if(j == width - 1)	{
					// The last row and needs to relay to the next element
					if(i == height - 1)
						gray_data[now - 1]      += (error *  0.4375);
					// Else need to relay to the next element, down, down_left
					else{
						gray_data[now - 1]      += (error *  0.4375);
						gray_data[next_row]     += (error *  0.3125);
						gray_data[next_row - 1] += (error *  0.0625);
					}
				}// end of the start element of the row
				// In the center
				else{
					// at the last row and need to relay to next element
					if(i == height - 1)
						gray_data[now - 1]      += (error *  0.4375);
					// Needs to relay to next, down, down_left, down_right
					else{
						gray_data[now - 1]      += (error *  0.4375);
						gray_data[next_row + 1] += (error *  0.1875);
						gray_data[next_row]     += (error *  0.3125);
						gray_data[next_row - 1] += (error *  0.0625);
					
					}
				}// end of the if(j...) series
			}//end of the for
		}// end of the odd number row
		// The even row 0, 2, 4 go from left to right
		else{
			for(j = 0; j < width; j++){
				// Calculate the offset of current and next row
				now = j + i * width;
				next_row = j + (i + 1) * width;

				// Determine the current dot
				if(gray_data[now] < 0.5){
					error = gray_data[now];
					gray_data[now] = 0;
				}
				else{
					error = gray_data[now] - 1.0;
					gray_data[now] = 1.0;
				}

				// Relay the error down
				// The start element of the row
				if(j == 0){
					// The last row and needs to relay to next
					if(i == height - 1)
						gray_data[now + 1]      += (error *  0.4375);
					// Else needs to relay to next, down, down_right
					else{
						gray_data[now + 1]      += (error *  0.4375);
						gray_data[next_row]     += (error *  0.3125);
						gray_data[next_row + 1] += (error *  0.0625);
					}
				}// end of the start element
				// The last element
				else if(j == width - 1){
					// The last row means the last element of the image 
					// and no need to relay to next
					// else needs to relay to down, down_left
					if(i != height - 1){
						gray_data[next_row]     += (error *  0.3125);
						gray_data[next_row - 1] += (error *  0.1875);
					}
				}// end of the last element
				else{
					// at the last row needs to relay to next
					if(i == height - 1)
						gray_data[now + 1]      += (error *  0.4375);
					// needs to relay to next, down, down_left, down_right
					else{
						gray_data[now + 1]      += (error *  0.4375);
						gray_data[next_row - 1] += (error *  0.1875);
						gray_data[next_row]     += (error *  0.3125);
						gray_data[next_row + 1] += (error *  0.0625);
					}
				}// end of the if(j..) series
			}// end of for j loop
		}//end of the even number operation
	}//end of for i loop

	/*for(i = 0; i < image_size; i++){
	// Check the Alpha
 		offset = 4 * i;
		if(data[offset + 3] < 128){
			data[offset + 3] = 0;
			for(j = 0; j < 3; j++)
				data[offset + j] = 0;
		}
		else{
			data[offset + 3] = 255;
			for(j = 0; j < 3; j++)
				data[offset + j] = gray_data[i] * 255.0;
		}
	}*/

	return gray_data;

}

_complex CIma::cor(double lambda)
{

	double pi=3.1415926535897932384626433832795;
	_complex c;
	double deuxpil=2*pi*lambda;
	if (fabs(lambda)<1E-10)
	{
		c.x=1.0;
		c.y=0.0;
	}
	else
	{
		c.x=sin(deuxpil)/(deuxpil);
		c.y=-(cos(deuxpil)-1)/(deuxpil);
	}
	return c;
}

void CIma::CalculCoefCorr()
{
	// calcul de coefficient de correction

	int s,r1;
	int N2=N/2;
	_complex corr,cors;
	double dN=(double)N;
	double dN2=N/2;
	for( r1=1; r1<= N; r1++)
	{
		for (s=1;s<=N;s++)
		{
			corr=cor((r1-dN2)/dN);
			cors=cor((s-dN2)/dN);

			fct[r1][s].x= (corr.x*cors.x-corr.y*cors.y)/dN ;
			fct[r1][s].y= (corr.x*cors.y+corr.y*cors.x)/dN ;
		}
	}
}

void CIma::InitPourFFT()
{

	N=1024;
	fct = (_complex **) malloc ((N+1)*sizeof(_complex *));
	for(int s=0;s<(N+1);s++)
		fct[s]=	(_complex*)malloc((N+1)*sizeof(_complex));

	CalculCoefCorr();

}



double * CIma::DescripteursFourier1(rgb *pImage)
{
	double *coef=NULL;
	double *coefR=NULL;
	double *coefV=NULL;
	double *coefB=NULL;

	double* real;
	double* imag;

	double pi=3.1415926535897932384626433832795;
	
	int p;
	long LargIm,HautIm;
	long l,c,i,j,m,r;
	LargIm=LFinal;
	HautIm=CFinal;
	//LargIm=N;
	//HautIm=N;
	if (LargIm!=HautIm) return 0; // si l'image n'est pas carree c'est pas bon

	
//	if (bContour) Sobel(pImage,CFinal,LFinal);


//	float *TabData=(float*)malloc(m_TailleFFT2D*sizeof(float));
//	float *fft=(float*)malloc(sizeof(float)*m_TailleFFT2D);


	_complex ***grid;
	grid = (_complex ***) malloc (3*sizeof(_complex **)); //new _complex**[3];
	for(m=0;m<3;m++)
	{
		grid[m] =(_complex**) malloc(LargIm*sizeof(_complex*));// new _complex*[LargIm];
		for(j=0;j<LargIm;j++) grid[m][j] =(_complex*)malloc(HautIm*sizeof(_complex));// new _complex[HautIm];
	}



//	RGBQUAD colorReal,colorImag;

	//collecte des données
	for (l=0;l<HautIm;l++) 
	{
		for (c=0;c<LargIm;c++) 
		{
				p=LargIm*l+c;
				grid[0][c][l].x = pImage[p].r;
				grid[0][c][l].y = pImage[p].r;
				grid[1][c][l].x = pImage[p].g;
				grid[1][c][l].y = pImage[p].g;
				grid[2][c][l].x = pImage[p].b;
				grid[2][c][l].y = pImage[p].b;
		} //c
	} // l

	/* Transform the rows */
	real = (double *)malloc(LargIm * sizeof(double));
	imag = (double *)malloc(LargIm * sizeof(double));


	m=0;
	while((1<<m)<HautIm) m++;

	for (r=0;r<3;r++)
	{
		for (l=0;l<HautIm;l++) 
		{
			for (c=0;c<LargIm;c++) 
			{
				real[c] = grid[r][c][l].x;
				imag[c] = grid[r][c][l].y;
			}

		//	pfft2D->do_fft(fft,TabData);

			FFT(1,m,real,imag);

			for (c = 0; c <LargIm; c ++)
			{
				grid[r][c][l].x = real[c];
				grid[r][c][l].y = imag[c];
			}

		}
	}

	/* Transform the columns */
	for (r=0;r<3;r++)
	{	
		for (c=0;c<LargIm;c++)
		{
			for (l=0;l<HautIm;l++) 
			{
	//			TabData[j] = grid[r][j][i].x;
				real[l] = grid[r][c][l].x;
				imag[l] = grid[r][c][l].y;
			}

			FFT(1,m,real,imag);

			for (l = 0; l <HautIm; l ++)
			{
				grid[r][c][l].x = real[l];
				grid[r][c][l].y = imag[l];
			} // ligne

		} // colonne
	} // couleur

	RecentrageGrilleCouleur(grid, HautIm, LargIm);

	double positionpix;
	double positionpix_x;
	double positionpix_y;
	double MaxR,MaxV,MaxB;
	double module;
	int Ray;
	int rmax,lmax;	
	rmax = LargIm/2;
	lmax = HautIm/2;


	// pour visu FFT

	/*
	MaxR=_cabs(grid[0][rmax][rmax]);
	for(l = 0;l<HautIm ;l++)
	{
		for( c = 0;c<LargIm ;c++)
		{
			module=_cabs(grid[0][c][l]);
			if (MaxR<module) MaxR=module;
		} // colonne
	} // ligne
	p=0;
	for(l = 0;l<HautIm ;l++)
	{
		for( c = 0;c<LargIm ;c++)
		{
			module=(_cabs(grid[0][c][l])/MaxR)*1024;
			if (module>255.0) module=255;
			pImage[p].r=module;
			pImage[p].g=module;
			pImage[p].b=module;
			p++;
		} // colonne
	} // ligne

   */


//	coef = new double[HautIm/2];
	coef = (double*)malloc(m_Npar*sizeof(double));
	coefR = (double*)malloc((lmax)*sizeof(double));
	coefV = (double*)malloc((lmax)*sizeof(double));
	coefB = (double*)malloc((lmax)*sizeof(double));
	long *Circonf=(long * )malloc(lmax*sizeof(long));

	for( i=0;  i<lmax; i++)
	{
		coefR[i] = 0.0;
		coefV[i] = 0.0;
		coefB[i] = 0.0;
		Circonf[i]=0;
	}
	for(l = 0;l<HautIm ;l++)
	{
		for( c = 0;c<LargIm ;c++)
		{
			positionpix_x = c - rmax;
			positionpix_y = l - rmax;

			positionpix_x*=positionpix_x;
			positionpix_y*=positionpix_y;
			positionpix=(int)sqrt(positionpix_x+positionpix_y);
			if(positionpix<rmax && positionpix>0)
			{
				Ray=(int)positionpix;
				coefR[Ray] += _cabs(grid[0][c][l]);
				coefV[Ray] += _cabs(grid[1][c][l]);
				coefB[Ray] += _cabs(grid[2][c][l]);
				Circonf[Ray]++;
			}
		
		} // colonne
	} // ligne
	coefR[0]=_cabs(grid[0][rmax][rmax]);
	coefV[0]=_cabs(grid[1][rmax][rmax]);
	coefB[0]=_cabs(grid[2][rmax][rmax]);
	MaxR=coefR[0];

	MaxV=coefV[0];
	MaxB=coefB[0];
	for( i=1;  i<lmax; i++)
	{
		coefR[i] = coefR[i]/(MaxR*Circonf[i]);
		coefV[i] = coefV[i]/(MaxV*Circonf[i]);
		coefB[i] = coefB[i]/(MaxB*Circonf[i]);
	}


   // recopie sans prendre le premier de chaque composante
	int np=m_Npar/3;
	for( i=0; i<np; i++)
	{
		coef[i]=coefR[i+1];
		coef[i+np]=coefV[i+1];
		coef[i+2*np]=coefB[i+1];
	}

	free(coefR);
	free(coefV);
	free(coefB);



	for(m=0;m<3;m++)
	{
		for(j=0;j<LargIm;j++) 
			free(grid[m][j]);
		free(grid[m]);
	}
	free(grid);

	free(real);
	free(imag);
	free(Circonf);

//	free(TabData);
//	free(fft);
	return coef;


}


void CIma::Sobel(rgb *PImagel,long Largl,long Hautl)
{
	int v1,v2,v3,h1,h2,h3,k2,k3,g1,g2,g;
	long li,co,p;
	int Taillel=Largl*Hautl*sizeof(rgb);
	rgb *Temp=(rgb*)calloc(Largl*Hautl,3);
	for (li=1;li<Hautl-1;li++) 	
	{
		for (co=1;co<Largl-1;co++) 
		{
			p=li*Largl+co;
			v1=(int)PImagel[p-Largl-1].r;
			v2=(int)PImagel[p-1].r;
			v3=(int)PImagel[p+Largl-1].r;
			g1=v1+2*v2+v3;
			h1=(int)PImagel[p-Largl+1].r;
			h2=(int)PImagel[p+1].r;
			h3=(int)PImagel[p+Largl+1].r;
			g2=h1+2*h2+h3;
			k2=(int)PImagel[p-Largl].r;
			k3=(int)PImagel[p+Largl].r;
			g=abs(g1-g2);
			g1=v1+2*k2+h1;
			g2=v3+2*k3+h3;
			g=g+abs(g1-g2);
			g=g/8;
			Temp[p].r=g;


			v1=(int)PImagel[p-Largl-1].g;
			v2=(int)PImagel[p-1].g;
			v3=(int)PImagel[p+Largl-1].g;
			g1=v1+2*v2+v3;
			h1=(int)PImagel[p-Largl+1].g;
			h2=(int)PImagel[p+1].g;
			h3=(int)PImagel[p+Largl+1].g;
			g2=h1+2*h2+h3;
			k2=(int)PImagel[p-Largl].g;
			k3=(int)PImagel[p+Largl].g;
			g=abs(g1-g2);
			g1=v1+2*k2+h1;
			g2=v3+2*k3+h3;
			g=g+abs(g1-g2);
			g=g/8;
			Temp[p].g=g;

			v1=(int)PImagel[p-Largl-1].b;
			v2=(int)PImagel[p-1].b;
			v3=(int)PImagel[p+Largl-1].b;
			g1=v1+2*v2+v3;
			h1=(int)PImagel[p-Largl+1].b;
			h2=(int)PImagel[p+1].b;
			h3=(int)PImagel[p+Largl+1].b;
			g2=h1+2*h2+h3;
			k2=(int)PImagel[p-Largl].b;
			k3=(int)PImagel[p+Largl].b;
			g=abs(g1-g2);
			g1=v1+2*k2+h1;
			g2=v3+2*k3+h3;
			g=g+abs(g1-g2);
			g=g/8;
			Temp[p].b=g;
		} // c
	}//l
	memcpy(PImagel,Temp,Taillel);
	free(Temp);
}

int CIma::LireDats(CString Chemin)
{
	BOOL bFound = TRUE;
	CFileFind ff;
	CString Temp;
	int n=0;

	TabRefNom.RemoveAll();

	if( ff.FindFile(Chemin) )
	{
		while(bFound)
		{
			bFound = ff.FindNextFile();
			if((!ff.IsDots()) && (!ff.IsDirectory()))
			{
				Temp=ff.GetFileName();
				TabRefNom.Add(Temp.Left(Temp.GetLength()-4));
				n++;
			}
		}
	}
	ff.Close();
	return n;

}

void CIma::FourierSurPointsContour(int SeuilG,int Apprend)
{

	// filtre de petrou
	// extrema
	int NbPtRnd; // nombre de points aléatoires sur le contour
	int LargFen=32;
	
	int SeuilPts=5; // si plus de 5 points d'un objet on le reconnait

	int i,li,co,p,t;
	long Taille=m_Largeur*m_Hauteur;

	int NbClassesMax=1000;

	int Histo[1000]; // 1000 classes max, à modifier

	BYTE *Cop=(BYTE *) malloc(Taille*sizeof(rgb));
	BYTE *Mir=(BYTE *) malloc(Taille*sizeof(rgb));

	for (i=0;i<NbClassesMax;i++) Histo[i]=0;
	
	if (Apprend) NbPtRnd=100;
	else NbPtRnd=200;

	SeuilPts=(NbPtRnd*5)/100;  // 5 %

	for (i=0;i<Taille;i++) 
	{ 
		int lum=(m_PImage[i].r/3)+(m_PImage[i].g/3)+(m_PImage[i].b/3);
/*		if (PImage[j]>10) PImage[j]-=10; else PImage[j]=0;
		if (PImage[j+1]>10) PImage[j+1]-=10; else PImage[j+1]=0;
		if (PImage[j+2]>10) PImage[j+2]-=10; else PImage[j+2]=0;*/
		Cop[i]=(BYTE) lum;
		Mir[i]=Cop[i];
	}
#ifdef MP	
	CPetrou Petrou;

	Petrou.m_TE = 4;
	Petrou.m_MaxLength = 8000;
	Petrou.m_MaxWidth = 15;
	Petrou.m_MaxDw = 3;
	Petrou.m_CheckDw = TRUE;
	Petrou.m_Check_Reg=TRUE;
	Petrou.m_MinLength = 10;
	Petrou.m_MaxReg=float(SeuilG)/255.0;
//	Petrou.m_PasSuivi=2;

	Petrou.Go(Cop,Mir,m_Hauteur,m_Largeur,Petrou.m_TE,0);
#endif   

	// Echantillonnage aleatoire sur les points de contour
	int ialeat,NbPts=0;
	for (i=0;i<Taille;i++) 
	{ 
		if (Mir[i]) {NbPts++;}
	}
	
	if (NbPts==0) return;

	CPoint *pPtCtr=new CPoint[NbPts];

	int k=0;
	for (li=LargFen;li<m_Hauteur-LargFen-1;li++) 	
	{
		for (co=LargFen;co<m_Largeur-LargFen-1;co++) 
		{
			if (Mir[m_Largeur*li+co])
			{
				pPtCtr[k].x=co;
				pPtCtr[k++].y=li;
			}
		}
	}

	NbPts=k; // on actualise car on a pas parcouru les bords 
	CPoint *pPtCtrAlea=new CPoint[NbPts];
	if (NbPts<NbPtRnd) NbPtRnd=NbPts;
	if (NbPts==0) return;

	m_Npar=(LargFen/2-1)*3;
	int LargFen2=LargFen/2;
	LFinal=LargFen;
	CFinal=LargFen;
	rgb *pFen=(rgb *)malloc(LargFen*LargFen*3);
	double *coef1;

	
	struct svm_node *x;
	int par;
	double v;
	int Vert;

	x = (struct svm_node *) malloc(m_Npar*sizeof(struct svm_node)*2);

	for (i=0;i<NbPtRnd;i++) 
	{ 
		  ialeat=rand();
		  ialeat*=(NbPts-1);
		  ialeat/=RAND_MAX;
		  pPtCtrAlea[i]=pPtCtr[ialeat];

	// descripteur f1 en chaque point (à l'interieur de l'objet)

	
		k=0;
		Vert=0;
		for (li=pPtCtrAlea[i].y-LargFen2;li<pPtCtrAlea[i].y+LargFen2;li++) 	
		{
			for (co=pPtCtrAlea[i].x-LargFen2;co<pPtCtrAlea[i].x+LargFen2;co++) 
			{
					p=m_Largeur*li+co; // indice dans image originale
					pFen[k].b=m_PImage[p].b;
					pFen[k].g=m_PImage[p].g;
					pFen[k++].r=m_PImage[p].r;
					if (m_PImage[p].g==255 && m_PImage[p].r==0 && m_PImage[p].b==0) Vert++;

			}
		}
/*
		for (li=pPtCtrAlea[i].y;li<pPtCtrAlea[i].y+LargFen;li++) 	
		{
			for (co=pPtCtrAlea[i].x;co<pPtCtrAlea[i].x+LargFen;co++) 
			{
					p=m_Largeur*li+co; // indice dans image originale
					pFen[k].b=m_PImage[p].b;
					pFen[k].g=m_PImage[p].g;
					pFen[k++].r=m_PImage[p].r;
					if (m_PImage[p].g==255 && m_PImage[p].r==0 && m_PImage[p].b==0) Vert++;

			}
		}*/
	


		if (Vert<(LargFen*LargFen)/4)
		{
		BITMAPINFO	 Info; // entête de l'image
		Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
		Info.bmiHeader.biWidth = CFinal;
		Info.bmiHeader.biHeight = LFinal;
		Info.bmiHeader.biPlanes = 1;   
		Info.bmiHeader.biSizeImage = sizeof(rgb)*CFinal*LFinal;
		Info.bmiHeader.biBitCount = 24;
		Info.bmiHeader.biClrUsed = 0;
		Info.bmiHeader.biCompression =  BI_RGB;
		Info.bmiHeader.biXPelsPerMeter = 0;
		Info.bmiHeader.biYPelsPerMeter = 0;
		Info.bmiHeader.biClrImportant = 0;
		CFile Fichier; 
		CString NomFichFenBmp;
		NomFichFenBmp.Format("fen%03d.bmp",i);
		/*if (Fichier.Open(NomFichFenBmp,CFile::modeWrite|CFile::modeCreate))
		  {   
			BITMAPFILEHEADER Tete;
			Tete.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER); 
			Tete.bfSize = Tete.bfOffBits+CFinal*LFinal*sizeof(rgb);
			Tete.bfType = 'MB'; // obligatoire
			Tete.bfReserved1 = 0;
			Tete.bfReserved2 = 0;  
			Fichier.Write(&Tete,(UINT) sizeof(BITMAPFILEHEADER)); // entête du fichier
			Fichier.Write(&Info,sizeof(BITMAPINFOHEADER));   // entête de l'image
			Fichier.Write(pFen,CFinal*LFinal*sizeof(rgb));	 // écriture des données
			Fichier.Close();

		}*/// si ouverture du fichier ok
			coef1=DescripteursFourier1(pFen);
			/*if (Apprend) // si apprentissage
			{
				// sauvegarde du .dat
				for (t=0;t<m_Npar;t++)
			}
			else // si mode décision
				{
					m_TabModule[t]=coef1[t];
				}*/
				AjouterDonneeDansRef("essailocal.dat");
			{
				for(par=0;par<m_Npar;par++)
				{ 
					x[par].index=par+1;
					x[par].value=(double)coef1[par];

				}//for par
				x[m_Npar].index = -1;
				v = svm_predict(model,x);
				Histo[(int)v-1]++;
			//	v=1;
				//m_Chaine=TabRefNom[(int)v-1];
				//m_Chaine.Format("%d",(int)v); 

			} // si app ou decide
			free(coef1);
		} // si pas fenetre verte
		else i--; // on prend un autre point


	} // fin pour tous les points

	free(x);
	if (!Apprend)
	{
		CString Temp;
		m_Chaine="";
		for (i=0;i<NbClassesMax;i++)
		{
			if (Histo[i]!=0)//>SeuilPts) 
			{
				Temp.Format("classe %d:%d  ",i,Histo[i]);
				m_Chaine+=Temp;
			}
		}
	}
/*
	for (i=0;i<Taille;i++) 
	{ 
		if (Mir[i]) {m_PImage[i].r=Mir[i]; m_PImage[i].g=m_PImage[i].b=0;}
	}
*/
	for (i=0;i<NbPtRnd;i++) 
	{
		p=m_Largeur*pPtCtrAlea[i].y+pPtCtrAlea[i].x; // indice dans image originale
		m_PImage[p].r=255; m_PImage[p].g=m_PImage[p].b=0;
	}

	free (Cop);
	free (Mir);

	delete [] pPtCtr;
	delete [] pPtCtrAlea;

	free (pFen);
//	SauveBMP("testpetrou.bmp");
}


//-----------------------------------------------
//Calcul d'un HyperGradient sur bloc descripteurs
//-----------------------------------------------
void CIma::HyperGradientDescripteurBloc(rgb *RGBrtn,rgb *pImage,int LargeurBloc, int HauteurBloc,int LargeurFenetre,int HauteurFenetre)
{
//------------------
//	Nombre de bloc
//------------------
int NombreBlocLargeur = LargeurFenetre/LargeurBloc;
int NombreBlocHauteur = HauteurFenetre/HauteurBloc;
//--------------------------------
// Variables de travail
//--------------------------------
int i=0;
int j=0;
int l=0;
int c=0;
int coul=0;
long p=0;
long p1=0;
long p00=0;
long p01=0;
long p10=0;
long p11=0;
long indicebloc=0;
long compteurresult=0;
double temp=0;
double minib=0,minir=0,minig=0;
double maxib=0,maxig=0,maxir=0;
int Nsur2 = LargeurBloc/2;
int m_NParsur3=0;
int deuxm_NParsur3=0;
//int HauteurMoinsNsur2 = HauteurFenetre - Nsur2;
//int LargeurMoinsNsur2 = LargeurFenetre - Nsur2;
//int idebut=0,ifin=0,jdebut=0,jfin=0;
CString NomFichier = "", ChaineTemp="";
int TailleImageBloc = NombreBlocLargeur*NombreBlocHauteur;
int NombrePixelBloc = LargeurBloc*HauteurBloc;
int TailleImage = LargeurFenetre*HauteurFenetre*3;

rgb * BlocPartieRecopie = (rgb *)malloc(sizeof(rgb)*NombrePixelBloc);
double * ResultatGradientb = (double *)malloc(sizeof(double)*TailleImageBloc);
double * ResultatGradientg = (double *)malloc(sizeof(double)*TailleImageBloc);
double * ResultatGradientr = (double *)malloc(sizeof(double)*TailleImageBloc);
//RGBrtn = (rgb *)malloc(sizeof(rgb)*TailleImageBloc);

for(i=0;i<TailleImageBloc;i++)
{
	ResultatGradientb[i]=0.0f;
	ResultatGradientg[i]=0.0f;
	ResultatGradientr[i]=0.0f;
	RGBrtn[i].b=0;
	RGBrtn[i].g=0;
	RGBrtn[i].r=0;
}



//------------------------------------
// Récupération du pointeur d'image
//------------------------------------

//m_PImage = pImage;

//------------------------------------
//	Récupération des caractéristiques
//------------------------------------
/*m_Largeur = LargeurFenetre;
m_Hauteur = HauteurFenetre;*/
LFinal = LargeurBloc;
CFinal = HauteurBloc;

//-1 si pas premier

m_Npar=(LargeurBloc/2)*3;

m_NParsur3=LargeurBloc/2;
deuxm_NParsur3=2*m_NParsur3;

//----------------------
//	HyperImage resultat
//----------------------

double * BlocPartieResultat;// = (double *)malloc(sizeof(double)*m_Npar);

double ** drg = (double **) malloc(sizeof(double*)*TailleImageBloc);
double ** drb = (double **) malloc(sizeof(double*)*TailleImageBloc);
double ** drr = (double **) malloc(sizeof(double*)*TailleImageBloc);

for(i=0;i<TailleImageBloc;i++)
{
		drg[i]= (double*) malloc(sizeof(double)*m_NParsur3);
		drb[i]= (double*) malloc(sizeof(double)*m_NParsur3);
		drr[i]= (double*) malloc(sizeof(double)*m_NParsur3);
}



//------------------------------------
//	Calcul des descripteurs
//------------------------------------
indicebloc=0;
for (l=0;l<HauteurFenetre;l=l+HauteurBloc) // hauteur
	{
	for (c=0;c<LargeurFenetre;c=c+LargeurBloc) //largeur
		{
		p1=0;
		for(i=0;i<HauteurBloc;i++)//balayage hauteur sur la petite fenetre
		{
			for(j=0;j<LargeurBloc;j++) //balayage largeur sur la petite fentre
			{
			p=(LargeurFenetre*(l+i)+(c+j));

			BlocPartieRecopie[p1].r=pImage[p].r;
			BlocPartieRecopie[p1].g=pImage[p].g;
			BlocPartieRecopie[p1].b=pImage[p].b;

			p1++;
			}
		}
		//calcul de la FFT pour le bloc
		BlocPartieResultat=DescripteursFourier1CompletGaussien(BlocPartieRecopie,2);
		//recopie des descripteurs dans le gros tableau de résultats
		for(i=0;i<m_NParsur3;i++)
		{
			drb[indicebloc][i]=BlocPartieResultat[i];
			drg[indicebloc][i]=BlocPartieResultat[i+m_NParsur3];
			drr[indicebloc][i]=BlocPartieResultat[i+deuxm_NParsur3];
		}
		indicebloc++;
		free(BlocPartieResultat);
		} // c
	}//l


//--------------------------------
//		Normalisation des blocs
//--------------------------------

// recherche du mini


minib= maxib=drb[0][0];
minig= maxig=drg[0][0];
minir= maxir=drr[0][0];


for(i=0;i<m_NParsur3;i++)
	{
	minib= maxib=drb[0][i];
	minig= maxig=drg[0][i];
	minir= maxir=drr[0][i];

	for(j=0;j<TailleImageBloc;j++)
	{
		if(drb[j][i]<minib)
		{
			minib = drb[j][i];
		}
		if(drb[j][i]>maxib)
		{
			maxib = drb[j][i];
		}
		if(drg[j][i]<minig)
		{
			minig = drg[j][i];
		}
		if(drg[j][i]>maxig)
		{
			maxig = drg[j][i];
		}
		if(drr[j][i]<minir)
		{
			minir = drr[j][i];
		}
		if(drr[j][i]>maxir)
		{
			maxir = drr[j][i];
		}
	}
	for(j=0;j<TailleImageBloc;j++)
	{
		//drb[j][i]-=minib;
		drb[j][i]*=255.0f;
		drb[j][i]/=(maxib);//-minib);

		//drg[j][i]-=minig;
		drg[j][i]*=255.0f;
		drg[j][i]/=(maxig);//-minig);

		//drr[j][i]-=minir;
		drr[j][i]*=255.0f;
		drr[j][i]/=(maxir);//-minir);
	}

}


//--------------------------------
//		Calcul de l'hypergradient
//--------------------------------

for(l=0;l<NombreBlocHauteur-1;l++)
{
	for(c=0;c<NombreBlocLargeur-1;c++)
	{
		p00=NombreBlocLargeur*l+c;
		p01=p00+1;
		p10=p00+NombreBlocLargeur;
		p11=p10+1;
		compteurresult = p00;
		temp = HyperDistance(drb[p00],drb[p11],m_Npar/3)+HyperDistance(drb[p10],drb[p01],m_Npar/3);
		temp /=2;
		ResultatGradientb[compteurresult]=temp;
		temp = HyperDistance(drg[p00],drg[p11],m_Npar/3)+HyperDistance(drg[p10],drg[p01],m_Npar/3);
		temp /=2;
		ResultatGradientg[compteurresult]=temp;
		temp = HyperDistance(drr[p00],drr[p11],m_Npar/3)+HyperDistance(drr[p10],drr[p01],m_Npar/3);
		temp /=2;
		ResultatGradientr[compteurresult]=temp;
		//compteurresult++;
	}
}



//-------------------------------
//	Normalisation et affichage 
//-------------------------------
//on recherche le mini
minib= maxib=ResultatGradientb[0];
minig= maxig=ResultatGradientg[0];
minir= maxir=ResultatGradientr[0];


	for(i=0;i<m_NParsur3;i++)
	{
		if(ResultatGradientb[i]<minib)
		{
			minib = ResultatGradientb[i];
		}
		if(ResultatGradientb[i]>maxib)
		{
			maxib = ResultatGradientb[i];
		}
		if(ResultatGradientg[i]<minig)
		{
			minig = ResultatGradientg[i];
		}
		if(ResultatGradientg[i]>maxig)
		{
			maxig = ResultatGradientg[i];
		}
		if(ResultatGradientr[i]<minir)
		{
			minir = ResultatGradientr[i];
		}
		if(ResultatGradientr[i]>maxir)
		{
			maxir = ResultatGradientr[i];
		}
	}


//on retranche le mini

	for(i=0;i<m_NParsur3;i++)
	{
	//ResultatGradientb[i]-=minib;
	ResultatGradientb[i]*=255.0f;
	ResultatGradientb[i]/=(maxib);//-minib);

	//ResultatGradientg[i]-=minig;
	ResultatGradientg[i]*=255.0f;
	ResultatGradientg[i]/=(maxig);//-minig);

	//ResultatGradientr[i]-=minir;
	ResultatGradientr[i]*=255.0f;
	ResultatGradientr[i]/=(maxir);//-minir);


	}


//on recopie pour l'affichage

if(m_PImage !=NULL) free(m_PImage);
m_PImage = (rgb*)malloc(sizeof(rgb)*TailleImageBloc);


for(j=0;j<m_NParsur3;j++)
{
		for(i=0;i<TailleImageBloc;i++)
		{
			RGBrtn[i].b=ResultatGradientb[i];
			RGBrtn[i].g=ResultatGradientg[i];
			RGBrtn[i].r=ResultatGradientr[i];

					

			m_PImage[i].b=drr[i][j];
			m_PImage[i].g=drr[i][j];
			m_PImage[i].r=drr[i][j];
			

			

		}
	NomFichier = "Fichier_num_r";
	ChaineTemp.Format("%d",j);
	NomFichier = NomFichier + ChaineTemp + ".bmp";
	SauveBMP(NomFichier);
	NomFichier = "";


		for(i=0;i<TailleImageBloc;i++)
		{
				

			m_PImage[i].b=drg[i][j];
			m_PImage[i].g=drg[i][j];
			m_PImage[i].r=drg[i][j];
			

		}
	NomFichier = "Fichier_num_g";
	ChaineTemp.Format("%d",j);
	NomFichier = NomFichier + ChaineTemp + ".bmp";
	SauveBMP(NomFichier);
	NomFichier = "";


		for(i=0;i<TailleImageBloc;i++)
		{


			m_PImage[i].b=drb[i][j];
			m_PImage[i].g=drb[i][j];
			m_PImage[i].r=drb[i][j];
		}

	NomFichier = "Fichier_num_b";
	ChaineTemp.Format("%d",j);
	NomFichier = NomFichier + ChaineTemp + ".bmp";
	SauveBMP(NomFichier);
	NomFichier = "";
}



//----------------------------
//		Liberation memoire
//----------------------------

free(BlocPartieRecopie);
free(ResultatGradientr);
free(ResultatGradientg);
free(ResultatGradientb);
//free(BlocPartieResultat);
for(i=0;i<TailleImageBloc;i++)
{
		free(drb[i]);
			free(drg[i]);
				free(drr[i]);
}
free(drb);
free(drg);
free(drr);

}

//--------------------------------------------
//		Calcul de la distance euclidienne 
//		entre 2 vecteurs de double
//--------------------------------------------

double CIma::HyperDistance(double *Vecteur1, double *Vecteur2, long TailleVecteur)
{
//-----------------------
//	Initialisation
//-----------------------
double Rtn=0;
long i=0;

double RacineTaille = sqrt((double)TailleVecteur);


//------------------------
// calcul de la distance
//------------------------

for(i=0;i<TailleVecteur;i++)
{
Rtn += (Vecteur1[i]-Vecteur2[i]) * (Vecteur1[i]-Vecteur2[i]) ;
}

Rtn = sqrt(Rtn)/RacineTaille;

//--------------------
// retour du resultat
//--------------------
return Rtn;
}



double * CIma::DescripteursFourier1Complet(rgb *pImage)
{
	double *coef=NULL;
	double *coefR=NULL;
	double *coefV=NULL;
	double *coefB=NULL;

	double* real;
	double* imag;

	double pi=3.1415926535897932384626433832795;
	
	int p;
	long LargIm,HautIm;
	long l,c,i,j,m,r;
	LargIm=LFinal;
	HautIm=CFinal;
	//LargIm=N;
	//HautIm=N;
	if (LargIm!=HautIm) return 0; // si l'image n'est pas carree c'est pas bon

	
//	if (bContour) Sobel(pImage,CFinal,LFinal);


//	float *TabData=(float*)malloc(m_TailleFFT2D*sizeof(float));
//	float *fft=(float*)malloc(sizeof(float)*m_TailleFFT2D);


	_complex ***grid;
	grid = (_complex ***) malloc (3*sizeof(_complex **)); //new _complex**[3];
	for(m=0;m<3;m++)
	{
		grid[m] =(_complex**) malloc(LargIm*sizeof(_complex*));// new _complex*[LargIm];
		for(j=0;j<LargIm;j++) grid[m][j] =(_complex*)malloc(HautIm*sizeof(_complex));// new _complex[HautIm];
	}



//	RGBQUAD colorReal,colorImag;

	//collecte des données
	for (l=0;l<HautIm;l++) 
	{
		for (c=0;c<LargIm;c++) 
		{
				p=LargIm*l+c;
				grid[0][c][l].x = pImage[p].r;
				grid[0][c][l].y = pImage[p].r;
				grid[1][c][l].x = pImage[p].g;
				grid[1][c][l].y = pImage[p].g;
				grid[2][c][l].x = pImage[p].b;
				grid[2][c][l].y = pImage[p].b;
		} //c
	} // l

	/* Transform the rows */
	real = (double *)malloc(LargIm * sizeof(double));
	imag = (double *)malloc(LargIm * sizeof(double));


	m=0;
	while((1<<m)<HautIm) m++;

	for (r=0;r<3;r++)
	{
		for (l=0;l<HautIm;l++) 
		{
			for (c=0;c<LargIm;c++) 
			{
				real[c] = grid[r][c][l].x;
				imag[c] = grid[r][c][l].y;
			}

		//	pfft2D->do_fft(fft,TabData);

			FFT(1,m,real,imag);

			for (c = 0; c <LargIm; c ++)
			{
				grid[r][c][l].x = real[c];
				grid[r][c][l].y = imag[c];
			}

		}
	}

	/* Transform the columns */
	for (r=0;r<3;r++)
	{	
		for (c=0;c<LargIm;c++)
		{
			for (l=0;l<HautIm;l++) 
			{
	//			TabData[j] = grid[r][j][i].x;
				real[l] = grid[r][c][l].x;
				imag[l] = grid[r][c][l].y;
			}

			FFT(1,m,real,imag);

			for (l = 0; l <HautIm; l ++)
			{
				grid[r][c][l].x = real[l];
				grid[r][c][l].y = imag[l];
			} // ligne

		} // colonne
	} // couleur

	RecentrageGrilleCouleur(grid, HautIm, LargIm);

	double positionpix;
	double positionpix_x;
	double positionpix_y;
	double MaxR,MaxV,MaxB;
	double module;
	int Ray;
	int rmax,lmax;	
	rmax = LargIm/2;
	lmax = HautIm/2;


	// pour visu FFT

	/*
	MaxR=_cabs(grid[0][rmax][rmax]);
	for(l = 0;l<HautIm ;l++)
	{
		for( c = 0;c<LargIm ;c++)
		{
			module=_cabs(grid[0][c][l]);
			if (MaxR<module) MaxR=module;
		} // colonne
	} // ligne
	p=0;
	for(l = 0;l<HautIm ;l++)
	{
		for( c = 0;c<LargIm ;c++)
		{
			module=(_cabs(grid[0][c][l])/MaxR)*1024;
			if (module>255.0) module=255;
			pImage[p].r=module;
			pImage[p].g=module;
			pImage[p].b=module;
			p++;
		} // colonne
	} // ligne

   */


//	coef = new double[HautIm/2];
	coef = (double*)malloc(m_Npar*sizeof(double));
	coefR = (double*)malloc((lmax)*sizeof(double));
	coefV = (double*)malloc((lmax)*sizeof(double));
	coefB = (double*)malloc((lmax)*sizeof(double));
	long *Circonf=(long * )malloc(lmax*sizeof(long));

	for( i=0;  i<lmax; i++)
	{
		coefR[i] = 0.0;
		coefV[i] = 0.0;
		coefB[i] = 0.0;
		Circonf[i]=0;
	}
	for(l = 0;l<HautIm ;l++)
	{
		for( c = 0;c<LargIm ;c++)
		{
			positionpix_x = c - rmax;
			positionpix_y = l - rmax;

			positionpix_x*=positionpix_x;
			positionpix_y*=positionpix_y;
			positionpix=(int)sqrt(positionpix_x+positionpix_y);
			if(positionpix<rmax && positionpix>0)
			{
				Ray=(int)positionpix;
				coefR[Ray] += _cabs(grid[0][c][l]);
				coefV[Ray] += _cabs(grid[1][c][l]);
				coefB[Ray] += _cabs(grid[2][c][l]);
				Circonf[Ray]++;
			}
		
		} // colonne
	} // ligne
	coefR[0]=_cabs(grid[0][rmax][rmax]);
	coefV[0]=_cabs(grid[1][rmax][rmax]);
	coefB[0]=_cabs(grid[2][rmax][rmax]);
	MaxR=1;//coefR[0];

	MaxV=1;//coefV[0];
	MaxB=1;//coefB[0];
	for( i=1;  i<lmax; i++)
	{
		coefR[i] = coefR[i]/(MaxR*Circonf[i]);
		coefV[i] = coefV[i]/(MaxV*Circonf[i]);
		coefB[i] = coefB[i]/(MaxB*Circonf[i]);
	}


   // recopie sans prendre le premier de chaque composante
	int np=m_Npar/3;
	for( i=0; i<np; i++)
	{
		coef[i]=coefR[i];
		coef[i+np]=coefV[i];
		coef[i+2*np]=coefB[i];
	}

	free(coefR);
	free(coefV);
	free(coefB);



	for(m=0;m<3;m++)
	{
		for(j=0;j<LargIm;j++) 
			free(grid[m][j]);
		free(grid[m]);
	}
	free(grid);

	free(real);
	free(imag);
	free(Circonf);

//	free(TabData);
//	free(fft);
	return coef;


}

//------------------------------------------------------------
//	Calcul l'hypergradient sur les descripteurs de fourier
//	en utilisant une fenêtre glissante pour la ROI de la FFT
//------------------------------------------------------------


void CIma::HyperGradientDescripteurGlissant(rgb *RGBrtn, rgb *pImage, int LargeurBloc, int HauteurBloc, int LargeurFenetre, int HauteurFenetre)
{
//--------------------------------
// Variables de travail
//--------------------------------
int i=0;
int j=0;
int l=0;
int c=0;
int coul=0;
long p=0;
long p1=0;
long p00=0;
long p01=0;
long p10=0;
long p11=0;
long indicebloc=0;
long compteurresult=0;
double temp=0;
double minib=0,minir=0,minig=0;
double maxib=0,maxig=0,maxir=0;
CString NomFichier = "", ChaineTemp="";
int NombrePixelBloc = LargeurBloc*HauteurBloc;
long TailleImage = LargeurFenetre*HauteurFenetre;
int Nsur2 = LargeurBloc/2;
int m_NParsur3=0;
int deuxm_NParsur3=0;
int HauteurMoinsNsur2 = HauteurFenetre - Nsur2;
int LargeurMoinsNsur2 = LargeurFenetre - Nsur2;
int idebut=0,ifin=0,jdebut=0,jfin=0;
rgb * BlocPartieRecopie = (rgb *)malloc(sizeof(rgb)*NombrePixelBloc);
double * ResultatGradientb = (double *)malloc(sizeof(double)*TailleImage);
double * ResultatGradientg = (double *)malloc(sizeof(double)*TailleImage);
double * ResultatGradientr = (double *)malloc(sizeof(double)*TailleImage);

for(i=0;i<TailleImage;i++)
{
	ResultatGradientb[i]=0.0f;
	ResultatGradientg[i]=0.0f;
	ResultatGradientr[i]=0.0f;
	RGBrtn[i].b=0;
	RGBrtn[i].g=0;
	RGBrtn[i].r=0;
}

//------------------------------------
//	Récupération des caractéristiques
//------------------------------------

LFinal = LargeurBloc;
CFinal = HauteurBloc;

//-1 si pas premier

m_Npar=(LargeurBloc/2)*3;
m_NParsur3 = LargeurBloc/2;
deuxm_NParsur3 = 2*m_NParsur3;
//----------------------
//	HyperImage resultat
//----------------------

double * BlocPartieResultat;

double ** drg = (double **) malloc(sizeof(double*)*TailleImage);
double ** drb = (double **) malloc(sizeof(double*)*TailleImage);
double ** drr = (double **) malloc(sizeof(double*)*TailleImage);

for(i=0;i<TailleImage;i++)
{
		drg[i]= (double*) malloc(sizeof(double)*m_NParsur3);
		drb[i]= (double*) malloc(sizeof(double)*m_NParsur3);
		drr[i]= (double*) malloc(sizeof(double)*m_NParsur3);
}



//------------------------------------
//	Calcul des descripteurs
//------------------------------------
indicebloc=0;
for (l=Nsur2;l<HauteurMoinsNsur2;l++) // hauteur
	{
	for (c=Nsur2;c<LargeurMoinsNsur2;c++) //largeur
		{
		p1=0;
		idebut=-Nsur2;
		ifin=Nsur2;
		jdebut=-Nsur2;
		jfin=+Nsur2;
		for(i=idebut;i<ifin;i++)//balayage hauteur sur la petite fenetre
		{
			for(j=jdebut;j<jfin;j++) //balayage largeur sur la petite fentre
			{
			p=(LargeurFenetre*(l+i)+(c+j));

			BlocPartieRecopie[p1].r=pImage[p].r;
			BlocPartieRecopie[p1].g=pImage[p].g;
			BlocPartieRecopie[p1].b=pImage[p].b;

			
			p1++;
			}
		}
		//calcul de la FFT pour le bloc
		BlocPartieResultat=DescripteursFourier1Complet(BlocPartieRecopie);
		//recopie des descripteurs dans le gros tableau de résultats
		indicebloc=LargeurFenetre*l+c;
		for(i=0;i<m_Npar/3;i++)
		{
			drb[indicebloc][i]=BlocPartieResultat[i];
			drg[indicebloc][i]=BlocPartieResultat[i+m_NParsur3];
			drr[indicebloc][i]=BlocPartieResultat[i+deuxm_NParsur3]; //attention ici
		}
		
		free(BlocPartieResultat);//allocation dans la fonction descripteur
		} // c
	}//l


//--------------------------------
//		Normalisation des blocs
//--------------------------------

// recherche du mini


minib= maxib=drb[0][0];
minig= maxig=drg[0][0];
minir= maxir=drr[0][0];


for(i=0;i<m_NParsur3;i++)
	{
	minib= maxib=drb[0][i];
	minig= maxig=drg[0][i];
	minir= maxir=drr[0][i];

	for(j=0;j<TailleImage;j++)
	{
		if(drb[j][i]<minib)
		{
			minib = drb[j][i];
		}
		if(drb[j][i]>maxib)
		{
			maxib = drb[j][i];
		}
		if(drg[j][i]<minig)
		{
			minig = drg[j][i];
		}
		if(drg[j][i]>maxig)
		{
			maxig = drg[j][i];
		}
		if(drr[j][i]<minir)
		{
			minir = drr[j][i];
		}
		if(drr[j][i]>maxir)
		{
			maxir = drr[j][i];
		}
	}
	for(j=0;j<TailleImage;j++)
	{
		//drb[j][i]-=minib;
		drb[j][i]*=255.0f;
		drb[j][i]/=(maxib);//-minib);

		//drg[j][i]-=minig;
		drg[j][i]*=255.0f;
		drg[j][i]/=(maxig);//-minig);

		//drr[j][i]-=minir;
		drr[j][i]*=255.0f;
		drr[j][i]/=(maxir);//-minir);
	}

}


//--------------------------------
//		Calcul de l'hypergradient
//--------------------------------

for(l=0;l<HauteurFenetre-1;l++)
{
	for(c=0;c<LargeurFenetre-1;c++)
	{
		p00=LargeurFenetre*l+c;
		p01=p00+1;
		p10=p00+LargeurFenetre;
		p11=p10+1;
		compteurresult = p00;
		temp = HyperDistance(drb[p00],drb[p11],m_NParsur3)+HyperDistance(drb[p10],drb[p01],m_NParsur3);
		temp /=2;
		ResultatGradientb[compteurresult]=temp;
		temp = HyperDistance(drg[p00],drg[p11],m_NParsur3)+HyperDistance(drg[p10],drg[p01],m_NParsur3);
		temp /=2;
		ResultatGradientg[compteurresult]=temp;
		temp = HyperDistance(drr[p00],drr[p11],m_NParsur3)+HyperDistance(drr[p10],drr[p01],m_NParsur3);
		temp /=2;
		ResultatGradientr[compteurresult]=temp;
	}
}

//-------------------------------
//	Normalisation et affichage 
//-------------------------------
//on recgerche le mini
minib= maxib=ResultatGradientb[0];
minig= maxig=ResultatGradientg[0];
minir= maxir=ResultatGradientr[0];


	for(i=0;i<m_NParsur3;i++)
	{
		if(ResultatGradientb[i]<minib)
		{
			minib = ResultatGradientb[i];
		}
		if(ResultatGradientb[i]>maxib)
		{
			maxib = ResultatGradientb[i];
		}
		if(ResultatGradientg[i]<minig)
		{
			minig = ResultatGradientg[i];
		}
		if(ResultatGradientg[i]>maxig)
		{
			maxig = ResultatGradientg[i];
		}
		if(ResultatGradientr[i]<minir)
		{
			minir = ResultatGradientr[i];
		}
		if(ResultatGradientr[i]>maxir)
		{
			maxir = ResultatGradientr[i];
		}
	}


//on retranche le mini

	for(i=0;i<m_NParsur3;i++)
	{
	//ResultatGradientb[i]-=minib;
	ResultatGradientb[i]*=255.0f;
	ResultatGradientb[i]/=(maxib);//-minib);

	//ResultatGradientg[i]-=minig;
	ResultatGradientg[i]*=255.0f;
	ResultatGradientg[i]/=(maxig);//-minig);

	//ResultatGradientr[i]-=minir;
	ResultatGradientr[i]*=255.0f;
	ResultatGradientr[i]/=(maxir);//-minir);


	}


//on recopie pour l'affichage

if(m_PImage !=NULL) free(m_PImage);
m_PImage = (rgb*)malloc(sizeof(rgb)*TailleImage);


for(j=0;j<m_NParsur3;j++)
{
	for(i=0;i<TailleImage;i++)
	{
		RGBrtn[i].b=ResultatGradientb[i];
		RGBrtn[i].g=ResultatGradientg[i];
		RGBrtn[i].r=ResultatGradientr[i];

		

		m_PImage[i].b=drr[i][j];
		m_PImage[i].g=drr[i][j];
		m_PImage[i].r=drr[i][j];
		

	}
	NomFichier = "Fichier_num_r";
	ChaineTemp.Format("%d",j);
	NomFichier = NomFichier + ChaineTemp + ".bmp";
	SauveBMP(NomFichier);
	NomFichier = "";


		for(i=0;i<TailleImage;i++)
	{

		
		m_PImage[i].b=drg[i][j];
		m_PImage[i].g=drg[i][j];
		m_PImage[i].r=drg[i][j];
		

	}
	NomFichier = "Fichier_num_g";
	ChaineTemp.Format("%d",j);
	NomFichier = NomFichier + ChaineTemp + ".bmp";
	SauveBMP(NomFichier);
	NomFichier = "";


		for(i=0;i<TailleImage;i++)
	{	

		m_PImage[i].b=drb[i][j];
		m_PImage[i].g=drb[i][j];
		m_PImage[i].r=drb[i][j];
		

		
		}
	NomFichier = "Fichier_num_b";
	ChaineTemp.Format("%d",j);
	NomFichier = NomFichier + ChaineTemp + ".bmp";
	SauveBMP(NomFichier);
	NomFichier = "";
}



//----------------------------
//		Liberation memoire
//----------------------------

free(BlocPartieRecopie);
free(ResultatGradientr);
free(ResultatGradientg);
free(ResultatGradientb);


for(i=0;i<TailleImage;i++)
{
		free(drb[i]);
			free(drg[i]);
				free(drr[i]);
}
free(drb);
free(drg);
free(drr);
}

void CIma::convolution_2D(float *source, float *dest, SImage info, SMasque masq)
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



//--------------------------------------------
//		Calcul descripteur de fourier
//		avec multiplication par une gaussien
//		dans le cadre de l'hypergradient à
//		fenetre glissante
//--------------------------------------------
double * CIma::DescripteursFourier1CompletGaussien(rgb *pImage, double sigma)
{
	double *coef=NULL;
	double *coefR=NULL;
	double *coefV=NULL;
	double *coefB=NULL;

	double* real;
	double* imag;

	double pi=3.1415926535897932384626433832795;
	
	int p;
	long LargIm,HautIm;
	long l,c,i,j,m,r;
	LargIm=LFinal;
	HautIm=CFinal;
	//LargIm=N;
	//HautIm=N;
	if (LargIm!=HautIm) return 0; // si l'image n'est pas carree c'est pas bon

	
//	if (bContour) Sobel(pImage,CFinal,LFinal);


//	float *TabData=(float*)malloc(m_TailleFFT2D*sizeof(float));
//	float *fft=(float*)malloc(sizeof(float)*m_TailleFFT2D);


//------------------------------------
//		Préparation de la gaussienne
//------------------------------------
	//-------------------------
	// mettre une valeur à E-T
	//-------------------------
	//double sigma=2;

	SMasque masq;
	//int l,c,p;
	double deuxpi=6.28318;
	double sigma2=sigma*sigma;
	double unsdeuxsigc=1.0/(2.0*sigma2);
	double norm=1.0/(deuxpi*sigma2);
	double norm2=1.0/(deuxpi*sigma2*sigma2);
	double norm3=1.0/(deuxpi*sigma2*sigma2*sigma2);
	double x,y;
	double arg;
	long TailleMasque=0;

	// --- parametres du masque ----
	masq.Haut=16;
	masq.Larg=16;
	masq.Centre_col=masq.Larg/2;
	masq.Centre_lig=masq.Haut/2;
	masq.Filtre=new float[masq.Haut*masq.Larg];

	for (l=0;l<masq.Haut;l++)
	{
		for (c=0;c<masq.Larg;c++)
		{
			p=masq.Larg*l+c;
			x=double(c-masq.Centre_col);
			y=double(l-masq.Centre_lig);
			arg=(x*x+y*y)*unsdeuxsigc;
			// gaussien
			masq.Filtre[p]=(float)(norm*exp(-arg));
			// dérivée de gaussien suivant c
		//	masq.Filtre[p]=(float)(-x*norm2*exp(-arg));
		//	derivee seconde
		//	masq.Filtre[p]=(float)(-x*norm2*exp(-arg));

		}
	}
	TailleMasque = masq.Haut*masq.Larg;
	drgb * ApresGauss = (drgb*) malloc(sizeof(drgb)*masq.Haut*masq.Larg);
	//Multiplication masque x imageRGB
	for(l=0;l<TailleMasque;l++)
	{
		ApresGauss[l].b = pImage[l].b*masq.Filtre[l];
		ApresGauss[l].g = pImage[l].g*masq.Filtre[l];
		ApresGauss[l].r = pImage[l].r*masq.Filtre[l];

	}



	_complex ***grid;
	grid = (_complex ***) malloc (3*sizeof(_complex **)); //new _complex**[3];
	for(m=0;m<3;m++)
	{
		grid[m] =(_complex**) malloc(LargIm*sizeof(_complex*));// new _complex*[LargIm];
		for(j=0;j<LargIm;j++) grid[m][j] =(_complex*)malloc(HautIm*sizeof(_complex));// new _complex[HautIm];
	}



//	RGBQUAD colorReal,colorImag;

	//collecte des données
	for (l=0;l<HautIm;l++) 
	{
		for (c=0;c<LargIm;c++) 
		{
				p=LargIm*l+c;
				grid[0][c][l].x = ApresGauss[p].r;
				grid[0][c][l].y = ApresGauss[p].r;
				grid[1][c][l].x = ApresGauss[p].g;
				grid[1][c][l].y = ApresGauss[p].g;
				grid[2][c][l].x = ApresGauss[p].b;
				grid[2][c][l].y = ApresGauss[p].b;
		} //c
	} // l

	/* Transform the rows */
	real = (double *)malloc(LargIm * sizeof(double));
	imag = (double *)malloc(LargIm * sizeof(double));


	m=0;
	while((1<<m)<HautIm) m++;

	for (r=0;r<3;r++)
	{
		for (l=0;l<HautIm;l++) 
		{
			for (c=0;c<LargIm;c++) 
			{
				real[c] = grid[r][c][l].x;
				imag[c] = grid[r][c][l].y;
			}

		//	pfft2D->do_fft(fft,TabData);

			FFT(1,m,real,imag);

			for (c = 0; c <LargIm; c ++)
			{
				grid[r][c][l].x = real[c];
				grid[r][c][l].y = imag[c];
			}

		}
	}

	/* Transform the columns */
	for (r=0;r<3;r++)
	{	
		for (c=0;c<LargIm;c++)
		{
			for (l=0;l<HautIm;l++) 
			{
	//			TabData[j] = grid[r][j][i].x;
				real[l] = grid[r][c][l].x;
				imag[l] = grid[r][c][l].y;
			}

			FFT(1,m,real,imag);

			for (l = 0; l <HautIm; l ++)
			{
				grid[r][c][l].x = real[l];
				grid[r][c][l].y = imag[l];
			} // ligne

		} // colonne
	} // couleur

	RecentrageGrilleCouleur(grid, HautIm, LargIm);

	double positionpix;
	double positionpix_x;
	double positionpix_y;
	double MaxR,MaxV,MaxB;
	double module;
	int Ray;
	int rmax,lmax;	
	rmax = LargIm/2;
	lmax = HautIm/2;


	// pour visu FFT

	/*
	MaxR=_cabs(grid[0][rmax][rmax]);
	for(l = 0;l<HautIm ;l++)
	{
		for( c = 0;c<LargIm ;c++)
		{
			module=_cabs(grid[0][c][l]);
			if (MaxR<module) MaxR=module;
		} // colonne
	} // ligne
	p=0;
	for(l = 0;l<HautIm ;l++)
	{
		for( c = 0;c<LargIm ;c++)
		{
			module=(_cabs(grid[0][c][l])/MaxR)*1024;
			if (module>255.0) module=255;
			pImage[p].r=module;
			pImage[p].g=module;
			pImage[p].b=module;
			p++;
		} // colonne
	} // ligne

   */


//	coef = new double[HautIm/2];
	coef = (double*)malloc(m_Npar*sizeof(double));
	coefR = (double*)malloc((lmax)*sizeof(double));
	coefV = (double*)malloc((lmax)*sizeof(double));
	coefB = (double*)malloc((lmax)*sizeof(double));
	long *Circonf=(long * )malloc(lmax*sizeof(long));

	for( i=0;  i<lmax; i++)
	{
		coefR[i] = 0.0;
		coefV[i] = 0.0;
		coefB[i] = 0.0;
		Circonf[i]=0;
	}
	for(l = 0;l<HautIm ;l++)
	{
		for( c = 0;c<LargIm ;c++)
		{
			positionpix_x = c - rmax;
			positionpix_y = l - rmax;

			positionpix_x*=positionpix_x;
			positionpix_y*=positionpix_y;
			positionpix=(int)sqrt(positionpix_x+positionpix_y);
			if(positionpix<rmax && positionpix>0)
			{
				Ray=(int)positionpix;
				coefR[Ray] += _cabs(grid[0][c][l]);
				coefV[Ray] += _cabs(grid[1][c][l]);
				coefB[Ray] += _cabs(grid[2][c][l]);
				Circonf[Ray]++;
			}
		
		} // colonne
	} // ligne
	coefR[0]=_cabs(grid[0][rmax][rmax]);
	coefV[0]=_cabs(grid[1][rmax][rmax]);
	coefB[0]=_cabs(grid[2][rmax][rmax]);
	MaxR=1;//coefR[0];

	MaxV=1;//coefV[0];
	MaxB=1;//coefB[0];
	for( i=1;  i<lmax; i++)
	{
		coefR[i] = coefR[i]/(MaxR*Circonf[i]);
		coefV[i] = coefV[i]/(MaxV*Circonf[i]);
		coefB[i] = coefB[i]/(MaxB*Circonf[i]);
	}


   // recopie sans prendre le premier de chaque composante
	int np=m_Npar/3;
	for( i=0; i<np; i++)
	{
		coef[i]=coefR[i];
		coef[i+np]=coefV[i];
		coef[i+2*np]=coefB[i];
	}

	free(coefR);
	free(coefV);
	free(coefB);



	for(m=0;m<3;m++)
	{
		for(j=0;j<LargIm;j++) 
			free(grid[m][j]);
		free(grid[m]);
	}
	free(grid);

	free(real);
	free(imag);
	free(Circonf);

//	free(TabData);
//	free(fft);
	return coef;
}

void CIma::HyperGradientGlissantZeroPadding(rgb *RGBrtn, rgb *pImage, int LargeurBloc, int HauteurBloc, int LargeurFenetre, int HauteurFenetre)
{
//--------------------------------
// Variables de travail
//--------------------------------
int i=0;
int j=0;
int l=0;
int c=0;
int w=0;
int coul=0;
long p=0;
long p1=0;
long p00=0;
long p01=0;
long p10=0;
long p11=0;
long indicebloc=0;
long compteurresult=0;
double temp=0;
double minib=0,minir=0,minig=0;
double maxib=0,maxig=0,maxir=0;
CString NomFichier = "", ChaineTemp="";
int NombrePixelBloc = LargeurBloc*HauteurBloc;
long TailleImage = LargeurFenetre*HauteurFenetre;
int Nsur2 = LargeurBloc/2;
int m_NParsur3=0;
int deuxm_NParsur3=0;
int HauteurMoinsNsur2 = HauteurFenetre - Nsur2;
int LargeurMoinsNsur2 = LargeurFenetre - Nsur2;
int idebut=0,ifin=0,jdebut=0,jfin=0;

int LargeurBlocSur2 = LargeurBloc/2;
int HauteurBlocSur2 = HauteurBloc/2;

int LargeurBlocSur4 = LargeurBloc/4;
int HauteurBlocSur4 = LargeurBloc/4;

int CCentrageInit=0,LCentrageInit=0,CCentrageFinal=0,LCentrageFinal=0;


int * TableCoordCentrage = (int *)malloc(sizeof(int)*LargeurBlocSur2* HauteurBlocSur2);

rgb * BlocPartieRecopie = (rgb *)malloc(sizeof(rgb)*NombrePixelBloc);
double * ResultatGradientb = (double *)malloc(sizeof(double)*TailleImage);
double * ResultatGradientg = (double *)malloc(sizeof(double)*TailleImage);
double * ResultatGradientr = (double *)malloc(sizeof(double)*TailleImage);

for(i=0;i<TailleImage;i++)
{
	ResultatGradientb[i]=0.0f;
	ResultatGradientg[i]=0.0f;
	ResultatGradientr[i]=0.0f;
	RGBrtn[i].b=0;
	RGBrtn[i].g=0;
	RGBrtn[i].r=0;
}

for(i=0;i<NombrePixelBloc;i++)
{
	BlocPartieRecopie[i].b=0;
	BlocPartieRecopie[i].g=0;
	BlocPartieRecopie[i].r=0;

}

CCentrageInit = LargeurBlocSur2-8;
LCentrageInit = HauteurBlocSur2-8;
CCentrageFinal = LargeurBlocSur2+8;
LCentrageFinal = HauteurBlocSur2+8;


//--------------------------------------
//	On rempli la LUT pour le recentrage
//--------------------------------------
/*i=0;
for(l=LCentrageInit;l<LCentrageFinal;l++) // hauteur
{
	for(c=CCentrageInit;c<CCentrageFinal;c++) //largeur
	{
	TableCoordCentrage[i]=LargeurFenetre*l+c;
	i++;
	}
}*/
//------------------------------------
//	Récupération des caractéristiques
//------------------------------------

LFinal = LargeurBloc;
CFinal = HauteurBloc;

//-1 si pas premier

m_Npar=(LargeurBloc/2)*3;
m_NParsur3 = LargeurBloc/2;
deuxm_NParsur3 = 2*m_NParsur3;
//----------------------
//	HyperImage resultat
//----------------------

double * BlocPartieResultat;

double ** drg = (double **) malloc(sizeof(double*)*TailleImage);
double ** drb = (double **) malloc(sizeof(double*)*TailleImage);
double ** drr = (double **) malloc(sizeof(double*)*TailleImage);

for(i=0;i<TailleImage;i++)
{
		drg[i]= (double*) malloc(sizeof(double)*m_NParsur3);
		drb[i]= (double*) malloc(sizeof(double)*m_NParsur3);
		drr[i]= (double*) malloc(sizeof(double)*m_NParsur3);
}



//------------------------------------
//	Calcul des descripteurs
//------------------------------------


int LFinBalayage = HauteurFenetre-LCentrageInit;
int CFinBalayage = LargeurFenetre-CCentrageInit;

indicebloc=0;
for (l=LCentrageInit;l<LFinBalayage;l++) // hauteur
	{
	for (c=CCentrageInit;c<CFinBalayage;c++) //largeur
		{
		//p1=0;
		w=0;
		idebut=-Nsur2;
		ifin=Nsur2;
		jdebut=-Nsur2;
		jfin=Nsur2;
		for(i=0;i<HauteurBlocSur2;i++)//balayage hauteur sur la petite fenetre
		{
			for(j=0;j<LargeurBlocSur2;j++) //balayage largeur sur la petite fentre
			{
				
					p=(LargeurFenetre*(l+i)+(c+j));
					//p1=TableCoordCentrage[w]; //----> si besoin d'une LUT

					p1=LargeurBloc*(i+LargeurBlocSur4)+(j+LargeurBlocSur4);
					BlocPartieRecopie[p1].r=pImage[p].r;
					BlocPartieRecopie[p1].g=pImage[p].g;
					BlocPartieRecopie[p1].b=pImage[p].b;

					
					//w++; //----------> LUT
				
			
			}
		}
		//calcul de la FFT pour le bloc
		BlocPartieResultat=DescripteursFourier1Complet(BlocPartieRecopie);
		//recopie des descripteurs dans le gros tableau de résultats
		indicebloc=LargeurFenetre*l+c;
		for(i=0;i<m_Npar/3;i++)
		{
			drb[indicebloc][i]=BlocPartieResultat[i];
			drg[indicebloc][i]=BlocPartieResultat[i+m_NParsur3];
			drr[indicebloc][i]=BlocPartieResultat[i+deuxm_NParsur3]; //attention ici
		}
		
		free(BlocPartieResultat);//allocation dans la fonction descripteur
		} // c
	}//l


//--------------------------------
//		Normalisation des blocs
//--------------------------------

// recherche du mini


minib= maxib=drb[0][0];
minig= maxig=drg[0][0];
minir= maxir=drr[0][0];


for(i=0;i<m_NParsur3;i++)
	{
	minib= maxib=drb[0][i];
	minig= maxig=drg[0][i];
	minir= maxir=drr[0][i];

	for(j=0;j<TailleImage;j++)
	{
		if(drb[j][i]<minib)
		{
			minib = drb[j][i];
		}
		if(drb[j][i]>maxib)
		{
			maxib = drb[j][i];
		}
		if(drg[j][i]<minig)
		{
			minig = drg[j][i];
		}
		if(drg[j][i]>maxig)
		{
			maxig = drg[j][i];
		}
		if(drr[j][i]<minir)
		{
			minir = drr[j][i];
		}
		if(drr[j][i]>maxir)
		{
			maxir = drr[j][i];
		}
	}
	for(j=0;j<TailleImage;j++)
	{
		//drb[j][i]-=minib;
		drb[j][i]*=255.0f;
		drb[j][i]/=(maxib);//-minib);

		//drg[j][i]-=minig;
		drg[j][i]*=255.0f;
		drg[j][i]/=(maxig);//-minig);

		//drr[j][i]-=minir;
		drr[j][i]*=255.0f;
		drr[j][i]/=(maxir);//-minir);
	}

}


//--------------------------------
//		Calcul de l'hypergradient
//--------------------------------

for(l=0;l<HauteurFenetre-1;l++)
{
	for(c=0;c<LargeurFenetre-1;c++)
	{
		p00=LargeurFenetre*l+c;
		p01=p00+1;
		p10=p00+LargeurFenetre;
		p11=p10+1;
		compteurresult = p00;
		temp = HyperDistance(drb[p00],drb[p11],m_NParsur3)+HyperDistance(drb[p10],drb[p01],m_NParsur3);
		temp /=2;
		ResultatGradientb[compteurresult]=temp;
		temp = HyperDistance(drg[p00],drg[p11],m_NParsur3)+HyperDistance(drg[p10],drg[p01],m_NParsur3);
		temp /=2;
		ResultatGradientg[compteurresult]=temp;
		temp = HyperDistance(drr[p00],drr[p11],m_NParsur3)+HyperDistance(drr[p10],drr[p01],m_NParsur3);
		temp /=2;
		ResultatGradientr[compteurresult]=temp;
	}
}

//-------------------------------
//	Normalisation et affichage 
//-------------------------------
//on recgerche le mini
minib= maxib=ResultatGradientb[0];
minig= maxig=ResultatGradientg[0];
minir= maxir=ResultatGradientr[0];


	for(i=0;i<m_NParsur3;i++)
	{
		if(ResultatGradientb[i]<minib)
		{
			minib = ResultatGradientb[i];
		}
		if(ResultatGradientb[i]>maxib)
		{
			maxib = ResultatGradientb[i];
		}
		if(ResultatGradientg[i]<minig)
		{
			minig = ResultatGradientg[i];
		}
		if(ResultatGradientg[i]>maxig)
		{
			maxig = ResultatGradientg[i];
		}
		if(ResultatGradientr[i]<minir)
		{
			minir = ResultatGradientr[i];
		}
		if(ResultatGradientr[i]>maxir)
		{
			maxir = ResultatGradientr[i];
		}
	}


//on retranche le mini

	for(i=0;i<m_NParsur3;i++)
	{
	//ResultatGradientb[i]-=minib;
	ResultatGradientb[i]*=255.0f;
	ResultatGradientb[i]/=(maxib);//-minib);

	//ResultatGradientg[i]-=minig;
	ResultatGradientg[i]*=255.0f;
	ResultatGradientg[i]/=(maxig);//-minig);

	//ResultatGradientr[i]-=minir;
	ResultatGradientr[i]*=255.0f;
	ResultatGradientr[i]/=(maxir);//-minir);


	}


//on recopie pour l'affichage

if(m_PImage !=NULL) free(m_PImage);
m_PImage = (rgb*)malloc(sizeof(rgb)*TailleImage);


for(j=0;j<m_NParsur3;j++)
{
	for(i=0;i<TailleImage;i++)
	{
		RGBrtn[i].b=ResultatGradientb[i];
		RGBrtn[i].g=ResultatGradientg[i];
		RGBrtn[i].r=ResultatGradientr[i];

		

		m_PImage[i].b=drr[i][j];
		m_PImage[i].g=drr[i][j];
		m_PImage[i].r=drr[i][j];
		

	}
	NomFichier = "Fichier_num_r";
	ChaineTemp.Format("%d",j);
	NomFichier = NomFichier + ChaineTemp + ".bmp";
	SauveBMP(NomFichier);
	NomFichier = "";


		for(i=0;i<TailleImage;i++)
	{

		
		m_PImage[i].b=drg[i][j];
		m_PImage[i].g=drg[i][j];
		m_PImage[i].r=drg[i][j];
		

	}
	NomFichier = "Fichier_num_g";
	ChaineTemp.Format("%d",j);
	NomFichier = NomFichier + ChaineTemp + ".bmp";
	SauveBMP(NomFichier);
	NomFichier = "";


		for(i=0;i<TailleImage;i++)
	{	

		m_PImage[i].b=drb[i][j];
		m_PImage[i].g=drb[i][j];
		m_PImage[i].r=drb[i][j];
		

		
		}
	NomFichier = "Fichier_num_b";
	ChaineTemp.Format("%d",j);
	NomFichier = NomFichier + ChaineTemp + ".bmp";
	SauveBMP(NomFichier);
	NomFichier = "";
}



//----------------------------
//		Liberation memoire
//----------------------------

free(BlocPartieRecopie);
free(ResultatGradientr);
free(ResultatGradientg);
free(ResultatGradientb);


for(i=0;i<TailleImage;i++)
{
		free(drb[i]);
			free(drg[i]);
				free(drr[i]);
}
free(drb);
free(drg);
free(drr);
}


double CIma::HyperDistanceComplexe(double *V1re, double *V1im, double *V2re, double *V2im, long TailleVecteur)
{


//-----------------------
//	Initialisation
//-----------------------
double Rtn=0;
long i=0;

double RacineTaille = sqrt((double)TailleVecteur);

//double ModulePremElemnt = sqrt(V1re[0]+V1im[0]);

CStdioFile FichierExportDistance;
CString ChaineDistance;

ChaineDistance = "";
//------------------------
// calcul de la distance
//------------------------

//if(FichierExportDistance.Open("EtudeStatDistance.txt",CFile::modeCreate|CFile::modeNoTruncate|CFile::modeWrite))
//{

	for(i=0;i<TailleVecteur;i++)
	{
	Rtn += (sqrt(((V1re[i]-V2re[i])*(V1re[i]-V2re[i]))+((V1im[i]-V2im[i])*(V1im[i]-V2im[i]))));//ModulePremElemnt;
	}

//	Rtn = Rtn/369.3037;
//	ChaineDistance.Format("%20.9f\n",Rtn);
//	FichierExportDistance.WriteString(ChaineDistance);
//	FichierExportDistance.Close();
//}
//--------------------
// retour du resultat
//--------------------
//else
//{
//	MessageBox(NULL,"Ouverture du fichier d'export impossible","Attention",MB_OK|MB_ICONEXCLAMATION);
//}
return Rtn;


}

double CIma::StandardDeviation(double *Tab, long LongueurTab)
{
double Rtn = 0.0f;
int i=0, j=0;
double temp=0.0f;
double * TabCarre= (double *)malloc(sizeof(double)*LongueurTab);
double MoyenneCarres = 0.0f;
double CarreMoyenne = 0.0f;


//copie du tableau
memcpy(TabCarre,Tab,LongueurTab);

//mise au carré
for(i=0;i<LongueurTab;i++)
{
temp += (TabCarre[i]*TabCarre[i]);
CarreMoyenne += (Tab[i]*Tab[i]);
}

MoyenneCarres = temp/LongueurTab;
CarreMoyenne = CarreMoyenne/LongueurTab;
CarreMoyenne = CarreMoyenne*CarreMoyenne;

Rtn = (MoyenneCarres - CarreMoyenne);
Rtn = sqrt(Rtn);

free(TabCarre);
return Rtn;
}

double CIma::HyperDistanceModuleFourier(double * aR, double * aI, double * bR, double * bI, int LargFiltre)
{
	int k;
	double d,dR,dI;
	double Ma,Mb;
	double Distance=0.0;

	for (k=0;k<LargFiltre;k++)
	{
		Ma=sqrt(aR[k]*aR[k]+aI[k]*aI[k]); // est il nécessaire de prendre la racine ?
		Mb=sqrt(bR[k]*bR[k]+bI[k]*bI[k]);
		Distance+=(Ma-Mb)*(Ma-Mb);
	} // for k
	Distance=sqrt(Distance);
	Distance/=LargFiltre;
	return Distance;
}

