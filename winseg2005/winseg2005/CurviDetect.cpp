// CurviDetect.cpp: implementation of the CCurviDetect class.
//
//////////////////////////////////////////////////////////////////////


#include "stdafx.h"
#include "winseg.h"
#include "CurviDetect.h"
#include "svm.h"
#include "progdlg.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


CCurviDetect::CCurviDetect()
{
I_i=I_x=I_y=I_g=NULL;
pi=3.14159f;
pi2=2.0f*pi;
pis2=pi/2.0f;
NPAR=4;
float n=(float)NPAR;
CoefGauss[0]=0.091f;
CoefGauss[1]=0.182f;
CoefGauss[2]=0.454f;
CoefGauss[3]=0.182f;
CoefGauss[4]=0.091f;

CompteurEcretage=0;


EtudeMaxi = 0.0f;
/*
CoefGauss[0]=0.2;
CoefGauss[1]=0.2;
CoefGauss[2]=0.2;
CoefGauss[3]=0.2;
CoefGauss[4]=0.2;*/
/*
CoefFeature[0]=0.2;
CoefFeature[1]=0.1666;
CoefFeature[2]=0.1666;
CoefFeature[3]=0.1666;
CoefFeature[4]=0.1666;
CoefFeature[5]=0.1336;
CoefFeature[6]=0.0; // not used
*/
CoefFeature[0]=1.0f/n;
CoefFeature[1]=1.0f/n;
CoefFeature[2]=1.0f/n;
CoefFeature[3]=1.0f/n;
CoefFeature[4]=1.0f/n;
CoefFeature[5]=1.0f/n;
CoefFeature[6]=1.0f/n; // not used

//------------------------------------------------------------------------------
	LargeurBordFourier  = 16; //Largeur pour le calcul de resemblance texture
	RacineLargeurBordFourier = (int)sqrt((double)LargeurBordFourier);//racine de la ligne d'étude : ---> doit être une puissance de 2
	
	int ValeurModule = LargeurBordFourier;



	
	SFourierGaucheRim = (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierGaucheRre = (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierGaucheGre = (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierGaucheGim = (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierGaucheBre = (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierGaucheBim = (double *)malloc(sizeof(double)*LargeurBordFourier);

	SFourierDroiteRre = (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierDroiteRim = (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierDroiteGre = (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierDroiteGim = (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierDroiteBre = (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierDroiteBim = (double *)malloc(sizeof(double)*LargeurBordFourier);

	SFourierSectionDroiteRre= (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierSectionDroiteRim= (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierSectionDroiteGre= (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierSectionDroiteGim= (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierSectionDroiteBre= (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierSectionDroiteBim= (double *)malloc(sizeof(double)*LargeurBordFourier);

	SFourierSectionGaucheRre= (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierSectionGaucheRim= (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierSectionGaucheGre= (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierSectionGaucheGim= (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierSectionGaucheBre= (double *)malloc(sizeof(double)*LargeurBordFourier);
	SFourierSectionGaucheBim= (double *)malloc(sizeof(double)*LargeurBordFourier);


	SFourierGaucheRModule = (double *)malloc(sizeof(double)*ValeurModule);
	SFourierGaucheGModule = (double *)malloc(sizeof(double)*ValeurModule);
	SFourierGaucheBModule = (double *)malloc(sizeof(double)*ValeurModule);
	SFourierDroiteRModule = (double *)malloc(sizeof(double)*ValeurModule);
	SFourierDroiteGModule = (double *)malloc(sizeof(double)*ValeurModule);
	SFourierDroiteBModule = (double *)malloc(sizeof(double)*ValeurModule);
	SFourierSectionDroiteRModule = (double *)malloc(sizeof(double)*ValeurModule);
	SFourierSectionDroiteGModule = (double *)malloc(sizeof(double)*ValeurModule);
	SFourierSectionDroiteBModule = (double *)malloc(sizeof(double)*ValeurModule);
	SFourierSectionGaucheRModule = (double *)malloc(sizeof(double)*ValeurModule);
	SFourierSectionGaucheGModule = (double *)malloc(sizeof(double)*ValeurModule);
	SFourierSectionGaucheBModule = (double *)malloc(sizeof(double)*ValeurModule);

//------------------------------------------------------------------------------------------
	//ouverture fichier pour l'étude statistique de la distance
if(OpenDistance=FichierExportDistance.Open("EtudeStatDistance.txt",CFile::modeCreate|CFile::modeNoTruncate|CFile::modeWrite))
{
	FichierExportDistance.SeekToEnd();
	
}
	//--------------------
					// retour du resultat
					//--------------------
else
{
	MessageBox(NULL,"Ouverture du fichier d'export impossible","Attention",MB_OK|MB_ICONEXCLAMATION);
}

}

CCurviDetect::~CCurviDetect()
{
	//ChaineEcretage.Format("%6.9f",EtudeMaxi);

	ChaineEcretage.Format("%d",CompteurEcretage);

	MessageBox(NULL,ChaineEcretage,"Info",MB_OK|MB_ICONINFORMATION);

	int i;
	for (i=0;i<Hauteur;i++)
	{
		free(I_i[i]);
		free(I_g[i]);
		free(I_x[i]);
		free(I_y[i]);
		free(I_map[i]);
		free(I_Ext[i]);
		free(I_Qr[i]);
		free(I_Qb[i]);
		free(I_Qg[i]);
		free(Cross[i]);
	}

/*	for (i=0;i<(wmax-wmin);i++) free(Gauss[i]);
	free(Gauss);*/        

	//---------------> modif le 18/10/2006 à remttre juste pour faire une manip rapide

	free(I_i);
	free(I_g);
	free(I_x);
	free(I_y);
	free(I_Qr);
	free(I_Qb);
	free(I_Qg);
	free(I_map);
	free(Cross);
//---------------------------------------------------------------------------------------------------------
//			Desallocation des tableaux gauches et droites de fourier
//---------------------------------------------------------------------------------------------------------
	free(SFourierGaucheRre);
	free(SFourierGaucheRim);
	free(SFourierGaucheGre); 
	free(SFourierGaucheGim);
	free(SFourierGaucheBre); 
	free(SFourierGaucheBim);

	free(SFourierDroiteRre);
	free(SFourierDroiteRim);
	free(SFourierDroiteGre);
	free(SFourierDroiteGim);
	free(SFourierDroiteBre);
	free(SFourierDroiteBim);

	free(SFourierSectionDroiteRre);
	free(SFourierSectionDroiteRim);
	free(SFourierSectionDroiteGre);
	free(SFourierSectionDroiteGim);
	free(SFourierSectionDroiteBre);
	free(SFourierSectionDroiteBim);

	free(SFourierSectionGaucheRre);
	free(SFourierSectionGaucheRim);
	free(SFourierSectionGaucheGre);
	free(SFourierSectionGaucheGim);
	free(SFourierSectionGaucheBre);
	free(SFourierSectionGaucheBim);

	free(SFourierGaucheRModule);
	free(SFourierGaucheGModule);
	free(SFourierGaucheBModule);
	free(SFourierDroiteRModule);
	free(SFourierDroiteGModule);
	free(SFourierDroiteBModule);

	free(SFourierSectionDroiteRModule);
	free(SFourierSectionDroiteGModule);
	free(SFourierSectionDroiteBModule);
	free(SFourierSectionGaucheRModule); 
	free(SFourierSectionGaucheGModule);
	free(SFourierSectionGaucheBModule);


	//-------------------------------------------------------
	//			si SVM
	//-------------------------------------------------------
//	svm_destroy_model(model);

	//-------------------------------------------------------
	// fermeture fichier distance
	//-------------------------------------------------------

	if(OpenDistance)
	{
	FichierExportDistance.Close();
	}
}

void CCurviDetect::Test(rgb *pIma, unsigned char *Result,int Haut, int Larg)
{

pImage=pIma;
Largeur = Larg;
Hauteur=Haut;
int l,c,p;
Init();

  p=0;
  for (l=0;l<Hauteur;l++)
  {
   for (c=0;c<Largeur;c++)
   {
	 // if (I_x[l][c]>0) Result[p]=I_i[l][c]*255;
	  //else Result[p]=0;
		Result[p]=(fabs(I_Ext[l][c]))*255;
	  p++;
   }
  }


}


void CCurviDetect::Go(rgb *pIma, unsigned char *Result,int Haut, int Larg)
{
	pImage=pIma;
	Largeur = Larg;
	Hauteur=Haut;
	int i,p,l,c;


	// initialisation des tableaux
	Init();


	CrossSectionMap();

	//Hough();


	float max,min;

	max=min=I_map[0][0];
	for (l=0;l<Hauteur;l++)
	for (c=0;c<Largeur;c++)
	{
	if (I_map[l][c]>max) max=I_map[l][c];
	if (I_map[l][c]<min) min=I_map[l][c];

	}

	p=0;

	if (max!=min)
	for (l=0;l<Hauteur;l++)
	for (c=0;c<Largeur;c++)
	{
	  Result[p++]=(BYTE)((255.0*(I_map[l][c]-min))/(max-min));
	// 	  Result[p++]=I_g[l][c];
	}
	

/*
   CStdioFile MonFic("sig.xls",CFile::modeWrite|CFile::modeCreate);
	CString Temp;
   for (l=0;l<100;l++)
   {
	   Temp.Format("%f\n",Histo_Sig[l]);
	   MonFic.WriteString(Temp);
   }
   MonFic.Close();*/


}


void CCurviDetect::Go(rgb *pIma, unsigned char *Result,int Haut, int Larg,CPoint point)
{
pImage=pIma;
Largeur = Larg;
Hauteur=Haut;

Init();

CrossSection(point.y,point.x,0);
}

void CCurviDetect::CrossSectionMap()
{

	// init


	int l,c;
	int m_MaxWidth=6;


//	for (l=0;l<100;l++)
//		Histo_Sig[l]=0.0;

 // for each pixel - limits

    for (l=0;l<Hauteur;l++) 
		for (c=0;c<Largeur;c++) 
			I_map[l][c]=0;//I_i[l][c];

	if (wmax>=Largeur) wmax=(Largeur*10)/100;

	CProgressDlg Progress;
	Progress.Create(NULL);
	Progress.SetRange(0,Largeur*Hauteur);//ca c'est la progress bar qui foire


	int pos=0;



	// balayage de l'image

	for (l=dw+1;l<Hauteur-dw;l+=1)
	{
	
		for (c=dw+1;c<Largeur-dw;c+=1)
		{
			Progress.SetPos(pos);
			pos++;
			CrossSection(l,c,1);
		} 
	}

	/*l=218;
	c=251;
	CrossSection(l,c,1);*/


}


void CCurviDetect::RemplitTableauIJCoupe(POINT *TabIJ,int x,int y,float dx,float dy,int w,int dw,int Dir,int *IptInit)
{
	float Pente=0.0F, OrdonneeOrigine=0.0F;
	char equation='0';
	long signe;
	long i,j,l;
	int wx;
	int wy;
	*IptInit=dw;
	
	int NbPointCoupe =(2*dw)+w+1;
	
 
	if (dx==0) equation = '1';   //la perp est verticale
	else if(dy==0) equation = '2';   // la perp est horizontale
	 else equation = '3';   //equation de la forme y=ax+b
    
	 switch (equation)
	{  
	    case '1':
			*IptInit=dw;
			j=y-dw*Dir;  // last change  //-------------------->> A voir
			i=x;
            for (l=0; l<NbPointCoupe; l++)
            {   
				if (i>=Largeur) i=Largeur-1;
				if (j>=Hauteur) j=Hauteur-1;
				if (i<0) i=0;
				if (j<0) j=0;
				TabIJ[l].x = (unsigned short) i;
                TabIJ[l].y = (unsigned short) j;
				j +=Dir;
			}
			break;
		case '2':
			*IptInit=dw;
			i=x-dw*Dir;
			j=y;
			for (l=0; l<NbPointCoupe; l++)
			{   
				if (i>=Largeur) i=Largeur-1;
				if (j>=Hauteur) j=Hauteur-1;
				if (i<0) i=0;
				if (j<0) j=0;
 				TabIJ[l].x = (unsigned short) i;
                TabIJ[l].y = (unsigned short) j;
				i +=Dir;
			}
			break;
		case '3':
			Pente = dy/dx; //initialise la pente a
			OrdonneeOrigine = y - Pente*x;  //initialise l'ordonnee b de y=a*x+b
 			if (Pente>0.0)
			{

				if (Pente>1.0)	
				{

					j = y-dw*Dir;
					i =  (long)  (((float)j-OrdonneeOrigine)/Pente);
					for (l=0; l<NbPointCoupe; l++) 
					{   
						if (i>=Largeur) i=Largeur-1;
						if (j>=Hauteur) j=Hauteur-1;
						if (i<0) i=0;
						if (j<0) j=0;
						TabIJ[l].x = (unsigned short) i;
                		TabIJ[l].y = (unsigned short) j;
						if (TabIJ[l].x==x && TabIJ[l].y==y) *IptInit=l;
						j+=Dir;   
						i = (long)  (((float)j-OrdonneeOrigine)/Pente);   

					}	
				}
				else
				{	
					i =  x-dw*Dir;
					j = Pente*(float)i+OrdonneeOrigine;
					for (l=0; l<NbPointCoupe; l++)         
					{
						if (i>=Largeur) i=Largeur-1;
						if (j>=Hauteur) j=Hauteur-1;
						if (i<0) i=0;
						if (j<0) j=0;
						TabIJ[l].x = (unsigned short) i;
               			TabIJ[l].y = (unsigned short) j;
						if (TabIJ[l].x==x && TabIJ[l].y==y) *IptInit=l;
					 	i +=Dir ;
					 	j = (long)(Pente*(float)i+OrdonneeOrigine);
					}                                   
				}	
			}
			else
			{   
				if (Pente<-1.0)
				{
					j = y+dw*Dir;
					i =  (long)  (((float)j-OrdonneeOrigine)/Pente);
					for (l=0; l<NbPointCoupe; l++)   
					{  
						if (i>=Largeur) i=Largeur-1;
						if (j>=Hauteur) j=Hauteur-1;
						if (i<0) i=0;
						if (j<0) j=0;
						TabIJ[l].x = (unsigned short) i;
                		TabIJ[l].y = (unsigned short) j;
						if (TabIJ[l].x==x && TabIJ[l].y==y) *IptInit=l;
						j-=Dir;   
						i = (long) (((float)j-OrdonneeOrigine)/Pente);   
					}	
				}
				else
				{
					i =  x-dw*Dir;
					j = Pente*(float)i+OrdonneeOrigine;
					for (l=0; l<NbPointCoupe; l++)     
					{
						if (i>=Largeur) i=Largeur-1;
						if (j>=Hauteur) j=Hauteur-1;
						if (i<0) i=0;
						if (j<0) j=0;
						TabIJ[l].x = (unsigned short) i;
						TabIJ[l].y = (unsigned short) j;
						if (TabIJ[l].x==x && TabIJ[l].y==y) *IptInit=l;
					 	i+=Dir ;
					 	j = (long)(Pente*(float)i+OrdonneeOrigine);
					}	
				}
			
			}	
			break;
			
		default:
			break;
	}
}

float CCurviDetect::Diff(float H2, float H1)
{
 float s,t,d;
/* t=H2-H1;
 if (t>=0) s=1; else s=-1;
 if (fabs(t)<=pis2) d=fabs(H2-H1);
 else d=fabs(H2-H1-pi*s);*/
 t=fabs(H2-H1);
 if (t>pi) d=2.0*pi-t;
 else d=t;
 return d;
}

float CCurviDetect::DiffAngle(POINTF p1, POINTF p2)
{

	float dotp=p1.x*p2.x+p1.y*p2.y;
	float normp=sqrt(p1.x*p1.x+p1.y*p1.y)*sqrt(p2.x*p2.x+p2.y*p2.y);
	float d;
	if (normp!=0.0) d=(1.0+dotp/normp)/2.0;
	else d=1.0;

	return d;
}


float CCurviDetect::Hue(rgb pImage)
{
	float num,denom,H;
	float fr,fg,fb;
	fr=(float)pImage.r;
	fg=(float)pImage.g;
	fb=(float)pImage.b;

	num=0.5*(2.0*fr-fg-fb);
	denom=sqrt((fr-fg)*(fr-fg)+(fr-fb)*(fg-fb));
	if (denom==0) H=0.0;
	else H=acos(num/denom);

	return H;
}

float CCurviDetect::AngleGradient(int l, int c)
{
	float a;
	a= atan2((double)I_y[l][c],(double)I_x[l][c]);
//	if (a<0) a=2.0*pi-a;
	return a;

}

void CCurviDetect::Quality(float *FeatureFin,float *pQ, int *pW, int *piw, POINT *Tab,int IptInit,int l, int c)
{

	float Angle, AngleInitial;
	float GradientProduct;
	float Qpipj=0.0;
	float fr,fg,fb;
	float MaxTest=pi;
	float H,m,HLeft, HRight,HMean;
	float AngleMax=pi/2.0;
	float DeltaAngle,DeltaAngleMax=0.5;//0.00001;
	float rac3=1.733;
	float Vrc,Vlc;
	float pti;
	float ptj;
	
	float NonBackground;
	float SimilarBorder;
	double TwoHighGrad;
	double TwoHighGradR;
	double TwoHighGradG;
	double TwoHighGradB;
	float OppositeAngle;

	CString ChaineTemp="";
	CStdioFile FichierExport;

	

	POINTF GradFin,GradInit;
	POINT Axis;
	
	int p,k,w,NbPts;
	int IptFinal;
	int Compteur=0,i=0;


	LargeurBordFourierSur2 = LargeurBordFourier/2;

	

	

//-----------------------------------------------------------------------------------------------------------
	// Pour les distances
	double   DistanceFourierGaucheR=0.0f ;
	double   DistanceFourierGaucheG=0.0f ;
	double   DistanceFourierGaucheB=0.0f ;
	double   DistanceFourierGauche =0.0f;

	double   DistanceFourierDroiteR=0.0f ;
	double   DistanceFourierDroiteG=0.0f ;
	double   DistanceFourierDroiteB=0.0f ;
	double   DistanceFourierDroite=0.0f ;

//------------------------------------------------------------------------------------------------------------

	//NbPts=(2*dw)+wmax+1+2*17;///------------------->>> ICI

//	AngleInitial=AngleGradient(l,c);
//	if (AngleInitial<0) pi2-AngleInitial;

//	for (k=0;k<NbPts;k++) I_map[Tab[k].y][Tab[k].x]=k*5; 
  	struct svm_node *x;
	x = (struct svm_node *) malloc(NPAR*sizeof(struct svm_node)*2);

	
	*pQ=1.0;
	int kw=0;

/*
	CStdioFile Fic("test.xls",CFile::modeWrite|CFile::modeCreate);
	CString Temp;
*/	
	for (w=wmin;w<wmax;w++)// on balaye sur la largeur
	{
		Qpipj=1;
		NbPts=2*dw+w+1;
		m_TailleTabCoupe = NbPts;
		//IptInit=IptInit+16;
		IptFinal=w+IptInit+1;
		CompteurFFT++;
		if (I_Ext[Tab[IptFinal].y][Tab[IptFinal].x]!=0 || Mode!=3)
		{
		/*	Angle=AngleGradient(Tab[IptFinal].y,Tab[IptFinal].x);
			if (Angle<0) pi2-Angle;
			DeltaAngle=Diff(Angle,AngleInitial);
			if (DeltaAngle>AngleMax || Mode!=3) // if opposite direction
			{*/
				// compute standard deviation inside the cross section
			GradInit.x=I_x[Tab[IptInit].y][Tab[IptInit].x];//gradient en x et y ci dessous
			GradInit.y=I_y[Tab[IptInit].y][Tab[IptInit].x];
			GradFin.x=I_x[Tab[IptFinal].y][Tab[IptFinal].x];
			GradFin.y=I_y[Tab[IptFinal].y][Tab[IptFinal].x];
			Axis.x=(Tab[IptInit].x+Tab[IptFinal].x)/2;
			Axis.y=(Tab[IptInit].y+Tab[IptFinal].y)/2;

			DeltaAngle=DiffAngle(GradInit,GradFin);//deviation?
			
			OppositeAngle=1.0-DeltaAngle;

			if (DeltaAngle<DeltaAngleMax || Mode!=3) // if opposite direction
			{

				float mean=0.0;
				int lar=0;
				float sig,sigi,sigQr,sigQb;
				float afi,afr,afb;
				float afiL,afrL,afbL;
				float afiR,afrR,afbR;
				rgb temp;

				afi=afr=afb=0.0;
				sigi=sigQr=sigQb=0.0;

				int nbiter=IptFinal-1-(IptInit+1);
				//ComputeCoefGauss(nbiter);
				int kg=nbiter-wmin+1;

//------------------------------------------------------------------------------------------------
				//calcul de l'écart type sur la section 
				for (k=IptInit+1;k<IptFinal-1;k++)
				{
	/*				afi+=CoGauss[lar]*I_i[Tab[k].y][Tab[k].x];
					afb+=CoGauss[lar]*I_Qb[Tab[k].y][Tab[k].x];
					afr+=CoGauss[lar]*I_Qr[Tab[k].y][Tab[k].x];*/
					
					afi+=Gauss[kg][lar]*I_Qg[Tab[k].y][Tab[k].x];
					afb+=Gauss[kg][lar]*I_Qb[Tab[k].y][Tab[k].x];
					afr+=Gauss[kg][lar]*I_Qr[Tab[k].y][Tab[k].x];
					lar++;
					sigi=sigi+ (I_Qg[Tab[k].y][Tab[k].x]*I_Qg[Tab[k].y][Tab[k].x]);
					sigQr=sigQr+ (I_Qr[Tab[k].y][Tab[k].x]*I_Qr[Tab[k].y][Tab[k].x]);
					sigQb=sigQb+ (I_Qb[Tab[k].y][Tab[k].x]*I_Qb[Tab[k].y][Tab[k].x]);
				}
			//	free(CoGauss);
				sigi=(lar*sigi-afi*afi)/(lar*(lar-1));
				sigQr=(lar*sigQr-afr*afr)/(lar*(lar-1)) ;
				sigQb=(lar*sigQb-afb*afb)/(lar*(lar-1)) ;

				if (sigi>0.0) sigi=sqrt( sigi );  else sigi=0.0;
				if (sigQr>0.0) sigQr=sqrt( sigQr);else sigQr=0.0;
				if (sigQb>0.0) sigQb=sqrt( sigQb );else sigQb=0.0;
				sig=(sigi+sigQr+sigQb)/3.0;
//-----------------------------------------------------------------------------------------------------
//				Calcul de la FFT sur la section
//-----------------------------------------------------------------------------------------------------
				//on vérifie que le gabarit tient dans l'image

				//---------------------------------------------
				//			Pour debuggage
				//---------------------------------------------
				/*if(DistanceEuclPoint(Tab[IptInit],Tab[IptFinal])>16)
					int hhhh=0;*/

			


					/*FichierExport.Open("ExportFFTnormliser.txt",CFile::modeCreate|CFile::modeWrite);
					for(k=0;k<16;k++)
					{
						ChaineTemp.Format("%d",IptInit);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\n");
					}
					FichierExport.Close();*/
					/*(IptInit-17>0)&&*/

				//if ((IptFinal>IptInit+dw)&&(IptFinal+dw<m_TailleTabCoupe))
				//{
					//allocation des tableaux de Fourier pour la section
				
					Compteur = 0;
	//------------------------------------------------------------------------------------------------------
					//Pour la section centrale coté gauche
					for (k=IptInit+1;k<IptInit+17;k++)
					{
						SFourierSectionGaucheRre[Compteur]= pImage[Tab[k].y*Largeur+Tab[k].x].r;
						SFourierSectionGaucheRim[Compteur]= 0.0f;
						SFourierSectionGaucheGre[Compteur]= pImage[Tab[k].y*Largeur+Tab[k].x].b;
						SFourierSectionGaucheGim[Compteur]= 0.0f;
						SFourierSectionGaucheBre[Compteur]= pImage[Tab[k].y*Largeur+Tab[k].x].g;
						SFourierSectionGaucheBim[Compteur]= 0.0f;
						Compteur++;
					}
					Compteur=0;

					//Pour la section centrale coté droit
					for (k=IptFinal-16;k<IptFinal;k++)
					{
						SFourierSectionDroiteRre[Compteur]= pImage[Tab[k].y*Largeur+Tab[k].x].r;
						SFourierSectionDroiteRim[Compteur]= 0.0f;
						SFourierSectionDroiteGre[Compteur]= pImage[Tab[k].y*Largeur+Tab[k].x].b;
						SFourierSectionDroiteGim[Compteur]= 0.0f;
						SFourierSectionDroiteBre[Compteur]= pImage[Tab[k].y*Largeur+Tab[k].x].g;
						SFourierSectionDroiteBim[Compteur]= 0.0f;
						Compteur++;
					}
					Compteur=0;
	//-------------------------------------------------------------------------------------------------------
					//pour le bord gauche
					for (k=IptInit-16;k<IptInit;k++)
					{
						SFourierGaucheRre[Compteur] = pImage[Tab[k].y*Largeur+Tab[k].x].r;
						SFourierGaucheRim[Compteur] = 0.0f;
						SFourierGaucheGre[Compteur] = pImage[Tab[k].y*Largeur+Tab[k].x].b;
						SFourierGaucheGim[Compteur] = 0.0f;
						SFourierGaucheBre[Compteur] = pImage[Tab[k].y*Largeur+Tab[k].x].g;
						SFourierGaucheBim[Compteur] = 0.0f;
						Compteur++;
					}
					Compteur=0;
	//-------------------------------------------------------------------------------------------------------
					//pour le bord droit

					for (k=IptFinal+1;k<IptFinal+17;k++)
					{
						SFourierDroiteRre[Compteur] = pImage[Tab[k].y*Largeur+Tab[k].x].r;
						SFourierDroiteRim[Compteur] = 0.0f;
						SFourierDroiteGre[Compteur] = pImage[Tab[k].y*Largeur+Tab[k].x].b;
						SFourierDroiteGim[Compteur] = 0.0f;
						SFourierDroiteBre[Compteur] = pImage[Tab[k].y*Largeur+Tab[k].x].g;
						SFourierDroiteBim[Compteur] = 0.0f;
						Compteur++;
					}
					Compteur=0;
	//--------------------------------------------------------------------------------------------------------
					//Calcul des FFT 1D

					ObjetCIma.FFT(1,RacineLargeurBordFourier,SFourierGaucheRre,SFourierGaucheRim);
					ObjetCIma.FFT(1,RacineLargeurBordFourier,SFourierDroiteRre,SFourierDroiteRim);

					ObjetCIma.FFT(1,RacineLargeurBordFourier,SFourierGaucheGre,SFourierGaucheGim);
					ObjetCIma.FFT(1,RacineLargeurBordFourier,SFourierDroiteGre,SFourierDroiteGim);

					ObjetCIma.FFT(1,RacineLargeurBordFourier,SFourierGaucheBre,SFourierGaucheBim);
					ObjetCIma.FFT(1,RacineLargeurBordFourier,SFourierDroiteBre,SFourierDroiteBim);

					ObjetCIma.FFT(1,RacineLargeurBordFourier,SFourierSectionGaucheRre,SFourierSectionGaucheRim);
					ObjetCIma.FFT(1,RacineLargeurBordFourier,SFourierSectionGaucheGre,SFourierSectionGaucheGim);
					ObjetCIma.FFT(1,RacineLargeurBordFourier,SFourierSectionGaucheBre,SFourierSectionGaucheBim);

					ObjetCIma.FFT(1,RacineLargeurBordFourier,SFourierSectionDroiteRre,SFourierSectionDroiteRim);
					ObjetCIma.FFT(1,RacineLargeurBordFourier,SFourierSectionDroiteGre,SFourierSectionDroiteGim);
					ObjetCIma.FFT(1,RacineLargeurBordFourier,SFourierSectionDroiteBre,SFourierSectionDroiteBim);

				


	//---------------------------------------------------------------------------------------------------------
					//-------------------------
					// calcul des modules
					//-------------------------

					/*ModuleComplex(SFourierGaucheRModule,SFourierGaucheRre,SFourierGaucheRim,RacineLargeurBordFourier);
					ModuleComplex(SFourierGaucheGModule,SFourierGaucheGre,SFourierGaucheGim,RacineLargeurBordFourier);
					ModuleComplex(SFourierGaucheBModule,SFourierGaucheBre,SFourierGaucheBim,RacineLargeurBordFourier);

					ModuleComplex(SFourierDroiteRModule,SFourierDroiteRre,SFourierDroiteRim,RacineLargeurBordFourier);
					ModuleComplex(SFourierDroiteGModule,SFourierDroiteGre,SFourierDroiteGim,RacineLargeurBordFourier);
					ModuleComplex(SFourierDroiteBModule,SFourierDroiteBre,SFourierDroiteBim,RacineLargeurBordFourier);

					ModuleComplex(SFourierSectionGaucheRModule,SFourierSectionGaucheRre,SFourierSectionGaucheRim,RacineLargeurBordFourier);
					ModuleComplex(SFourierSectionGaucheGModule,SFourierSectionGaucheGre,SFourierSectionGaucheGim,RacineLargeurBordFourier);
					ModuleComplex(SFourierSectionGaucheBModule,SFourierSectionGaucheBre,SFourierSectionGaucheBim,RacineLargeurBordFourier);

					ModuleComplex(SFourierSectionDroiteRModule,SFourierSectionDroiteRre,SFourierSectionDroiteRim,RacineLargeurBordFourier);
					ModuleComplex(SFourierSectionDroiteGModule,SFourierSectionDroiteGre,SFourierSectionDroiteGim,RacineLargeurBordFourier);
					ModuleComplex(SFourierSectionDroiteBModule,SFourierSectionDroiteBre,SFourierSectionDroiteBim,RacineLargeurBordFourier);*/

					

					/*FichierExport.Open("ExportFFT.txt",CFile::modeCreate|CFile::modeWrite);
					for(k=0;k<16;k++)
					{
						ChaineTemp.Format("%f",SFourierSectionDroiteRModule[k]);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\n");
					}
					FichierExport.Close();*/
					

	//---------------------------------------------------------------------------------------------------------
					//------------------------
					// Normalisation
					//------------------------

					/*NormalisationSpectre(SFourierGaucheRModule,LargeurBordFourier);
					NormalisationSpectre(SFourierGaucheGModule,LargeurBordFourier);
					NormalisationSpectre(SFourierGaucheBModule,LargeurBordFourier);

					NormalisationSpectre(SFourierDroiteRModule,LargeurBordFourier);
					NormalisationSpectre(SFourierDroiteGModule,LargeurBordFourier);
					NormalisationSpectre(SFourierDroiteBModule,LargeurBordFourier);

					NormalisationSpectre(SFourierSectionGaucheRModule,LargeurBordFourier);
					NormalisationSpectre(SFourierSectionGaucheGModule,LargeurBordFourier);
					NormalisationSpectre(SFourierSectionGaucheBModule,LargeurBordFourier);

					NormalisationSpectre(SFourierSectionDroiteRModule,LargeurBordFourier);
					NormalisationSpectre(SFourierSectionDroiteGModule,LargeurBordFourier);
					NormalisationSpectre(SFourierSectionDroiteBModule,LargeurBordFourier);*/

					/*NormalisationSpectreComplexe(SFourierGaucheRre,SFourierGaucheRim,LargeurBordFourier);
					NormalisationSpectreComplexe(SFourierGaucheGre,SFourierGaucheGim,LargeurBordFourier);
					NormalisationSpectreComplexe(SFourierGaucheBre,SFourierGaucheBim,LargeurBordFourier);

					NormalisationSpectreComplexe(SFourierDroiteRre,SFourierDroiteRim,LargeurBordFourier);
					NormalisationSpectreComplexe(SFourierDroiteGre,SFourierDroiteGim,LargeurBordFourier);
					NormalisationSpectreComplexe(SFourierDroiteBre,SFourierDroiteBim,LargeurBordFourier);


					NormalisationSpectreComplexe(SFourierSectionGaucheRre,SFourierSectionGaucheRim,LargeurBordFourier);
					NormalisationSpectreComplexe(SFourierSectionGaucheGre,SFourierSectionGaucheGim,LargeurBordFourier);
					NormalisationSpectreComplexe(SFourierSectionGaucheBre,SFourierSectionGaucheBim,LargeurBordFourier);


					NormalisationSpectreComplexe(SFourierSectionDroiteRre,SFourierSectionDroiteRim,LargeurBordFourier);
					NormalisationSpectreComplexe(SFourierSectionDroiteGre,SFourierSectionDroiteGim,LargeurBordFourier);
					NormalisationSpectreComplexe(SFourierSectionDroiteBre,SFourierSectionDroiteBim,LargeurBordFourier);*/


					/*FichierExport.Open("ExportFFTnormliser.txt",CFile::modeCreate|CFile::modeWrite);
					for(k=0;k<16;k++)
					{
						ChaineTemp.Format("%f",SFourierSectionDroiteRModule[k]);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\n");
					}
					FichierExport.Close();*/


					
	//---------------------------------------------------------------------------------------------------------
					//Calcul des hyper-distances 
					
					DistanceFourierDroiteR = ObjetCIma.HyperDistanceComplexe(SFourierDroiteRre,SFourierDroiteRim,SFourierSectionDroiteRre,SFourierSectionDroiteRim,LargeurBordFourierSur2);
					DistanceFourierDroiteG = ObjetCIma.HyperDistanceComplexe(SFourierDroiteGre,SFourierDroiteGim,SFourierSectionDroiteGre,SFourierSectionDroiteGim,LargeurBordFourierSur2);
					DistanceFourierDroiteB = ObjetCIma.HyperDistanceComplexe(SFourierDroiteBre,SFourierDroiteBim,SFourierSectionDroiteBre,SFourierSectionDroiteBim,LargeurBordFourierSur2);

					DistanceFourierGaucheR = ObjetCIma.HyperDistanceComplexe(SFourierGaucheRre,SFourierGaucheRim,SFourierSectionGaucheRre,SFourierSectionGaucheRim,LargeurBordFourierSur2);
					DistanceFourierGaucheG = ObjetCIma.HyperDistanceComplexe(SFourierGaucheGre,SFourierGaucheGim,SFourierSectionGaucheGre,SFourierSectionGaucheGim,LargeurBordFourierSur2);
					DistanceFourierGaucheB = ObjetCIma.HyperDistanceComplexe(SFourierGaucheBre,SFourierGaucheBim,SFourierSectionGaucheBre,SFourierSectionGaucheBim,LargeurBordFourierSur2);

					
					//écrêtage

					if(DistanceFourierDroiteR>1) 
					{
						DistanceFourierDroiteR=1.0f;
						CompteurEcretage++;
					}
					if(DistanceFourierDroiteG>1) 
					{
						DistanceFourierDroiteG=1.0f;
						CompteurEcretage++;
					}
					if(DistanceFourierDroiteB>1) 
					{
						DistanceFourierDroiteB=1.0f;
						CompteurEcretage++;
					}

					if(DistanceFourierGaucheR>1)
					{
						DistanceFourierGaucheR=1.0f;
						CompteurEcretage++;
					}
					if(DistanceFourierGaucheG>1) 
					{
						DistanceFourierGaucheG=1.0f;
						CompteurEcretage++;
					}
					if(DistanceFourierGaucheB>1) 
					{
						DistanceFourierGaucheB=1.0f;
						CompteurEcretage++;
					}

					
					/*if(DistanceFourierDroiteR>EtudeMaxi) 
					{
						EtudeMaxi=DistanceFourierDroiteR;
					}
					if(DistanceFourierDroiteG>EtudeMaxi) 
					{
						EtudeMaxi=DistanceFourierDroiteG;
					}
					if(DistanceFourierDroiteB>EtudeMaxi) 
					{
						EtudeMaxi=DistanceFourierDroiteB;
					}

					if(DistanceFourierGaucheR>EtudeMaxi) 
					{
						EtudeMaxi=DistanceFourierGaucheR;
					}
					if(DistanceFourierGaucheG>EtudeMaxi) 
					{
						EtudeMaxi=DistanceFourierGaucheG;
					}
					if(DistanceFourierGaucheB>EtudeMaxi) 
					{
						EtudeMaxi=DistanceFourierGaucheB;
					}*/


				

					





					/*if(OpenDistance)
					{
					

						
						ChaineDistance.Format("%20.9f\t",DistanceFourierDroiteR);
						FichierExportDistance.WriteString(ChaineDistance);

						ChaineDistance.Format("%20.9f\t",DistanceFourierDroiteG);
						FichierExportDistance.WriteString(ChaineDistance);

						ChaineDistance.Format("%20.9f\t",DistanceFourierDroiteB);
						FichierExportDistance.WriteString(ChaineDistance);

						ChaineDistance.Format("%20.9f\t",DistanceFourierGaucheR);
						FichierExportDistance.WriteString(ChaineDistance);

						ChaineDistance.Format("%20.9f\t",DistanceFourierGaucheG);
						FichierExportDistance.WriteString(ChaineDistance);

						ChaineDistance.Format("%20.9f\t",DistanceFourierGaucheB);
						FichierExportDistance.WriteString(ChaineDistance);

						FichierExportDistance.WriteString("\n");

					}*/
					//--------------------
					// retour du resultat
					//--------------------
					
					/*DistanceFourierGaucheR = ObjetCIma.HyperDistance(SFourierGaucheRModule,SFourierSectionGaucheRModule,LargeurBordFourierSur2);
					DistanceFourierGaucheG = ObjetCIma.HyperDistance(SFourierGaucheGModule,SFourierSectionGaucheGModule,LargeurBordFourierSur2);
					DistanceFourierGaucheB = ObjetCIma.HyperDistance(SFourierGaucheBModule,SFourierSectionGaucheBModule,LargeurBordFourierSur2);

					DistanceFourierDroiteR = ObjetCIma.HyperDistance(SFourierDroiteRModule,SFourierSectionDroiteRModule,LargeurBordFourierSur2);
					DistanceFourierDroiteG = ObjetCIma.HyperDistance(SFourierDroiteGModule,SFourierSectionDroiteGModule,LargeurBordFourierSur2);
					DistanceFourierDroiteB = ObjetCIma.HyperDistance(SFourierDroiteBModule,SFourierSectionDroiteBModule,LargeurBordFourierSur2);*/
					
				/*	if (FichierExport.Open("ExportDistancetextureGD.txt",CFile::modeCreate|CFile::modeWrite|CFile::modeNoTruncate))
					{
						FichierExport.SeekToEnd();
						ChaineTemp.Format("%f",DistanceFourierGaucheR);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\t");
						ChaineTemp.Format("%f",DistanceFourierGaucheG);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\t");
						ChaineTemp.Format("%f",DistanceFourierGaucheB);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\t");
						ChaineTemp.Format("%f",DistanceFourierDroiteR);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\t");
						ChaineTemp.Format("%f",DistanceFourierDroiteG);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\t");
						ChaineTemp.Format("%f",DistanceFourierDroiteB);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\t");
						FichierExport.WriteString("\n");
					
						FichierExport.Close();
					}
					else
						int stop=0;*/

					int ahah=0;
					
				//}
//----------------------------------------------------------------------------------------------------------
				// compute mean on the left
				afiL=afrL=afbL=0.0;
				float sgauss;
				sgauss=0;
				for (k=0;k<dw;k++)
				{
					afiL+=CoefGauss[k]*I_Qg[Tab[k].y][Tab[k].x];
					afbL+=CoefGauss[k]*I_Qb[Tab[k].y][Tab[k].x];
					afrL+=CoefGauss[k]*I_Qr[Tab[k].y][Tab[k].x];
					sgauss+=CoefGauss[k];
				}
//				afiL/=dw;	afbL/=dw;	afrL/=dw;
				afiL/=sgauss;	afbL/=sgauss;	afrL/=sgauss;
//-----------------------------------------------------------------------------------------------------
				// compute mean on the right
				afiR=afbR=afrR=0.0;
				int k2=0;
				sgauss=0;
				for (k=IptFinal;k<NbPts;k++)
				{
					afiR+=CoefGauss[k2]*I_Qg[Tab[k].y][Tab[k].x];
					afbR+=CoefGauss[k2]*I_Qb[Tab[k].y][Tab[k].x];
					afrR+=CoefGauss[k2]*I_Qr[Tab[k].y][Tab[k].x];
					sgauss+=CoefGauss[k2];
					k2++;
				}
//				afiR/=dw;	afbR/=dw;	afrR/=dw;
				afiR/=sgauss;	afbR/=sgauss;	afrR/=sgauss;
//------------------------------------------------------------------------------------------------------
				//proba d'avoir les bords similaires
				m=(afiR-afiL)*(afiR-afiL);
				m+=(afbR-afbL)*(afbR-afbL);
				m+=(afrR-afrL)*(afrR-afrL);
				m=sqrt(m)/3;
				SimilarBorder=1.0-m;

//------------------------------------------------------------------------------------------------------
				//proba d'avoir deux grand gradient attention couleur voir pluôt sur du niveau de gris
				// abs(VL-Vc)
			//	afi/=lar; afb/=lar; afr/=lar;
				//double temp =((afi+afb+afr)/3)-((afiL+afbL+afrL)/3)*((afi+afb+afr)/3)-((afiL+afbL+afrL)/3);
				Vlc=(afi-afiL)*(afi-afiL);
				Vlc+=(afb-afbL)*(afb-afbL);
				Vlc+=(afr-afrL)*(afr-afrL);
				//Vlc=temp;
				Vlc=(sqrt(Vlc)/rac3);

				// abs(VR-Vc)
				//temp =((afi+afb+afr)/3)-((afiR+afbR+afrR)/3)*((afi+afb+afr)/3)-((afiR+afbR+afrR)/3);
				Vrc=(afi-afiR)*(afi-afiR);
				Vrc+=(afb-afbR)*(afb-afbR);
				Vrc+=(afr-afrR)*(afr-afrR);
				//Vrc=temp;
				Vrc=(sqrt(Vrc)/rac3);

				float gpi=I_g[Tab[IptInit].y][Tab[IptInit].x];
				float gpj=I_g[Tab[IptFinal].y][Tab[IptFinal].x];
				TwoHighGrad=gpi*gpj;

				double TempMax[3];
	
				TwoHighGradR = DistanceFourierGaucheR ;//* DistanceFourierDroiteR;
				TwoHighGradG = DistanceFourierGaucheG ;//* DistanceFourierDroiteG;
				TwoHighGradB = DistanceFourierGaucheB ;//* DistanceFourierDroiteB;

				TempMax[0]=TwoHighGradR;
				TempMax[1]=TwoHighGradG;
				TempMax[2]=TwoHighGradB;

				

				/*if (TempMax[0]!=0)
				{
					int arretdebug=0;
				}*/

				
				double MaxiRespComp = 0.0f;
				for(k=0;k<3;k++)
				{
					if(TempMax[i]>MaxiRespComp)
					{
						MaxiRespComp=TempMax[i];
					}
				}

				TwoHighGrad=MaxiRespComp;

				/*if (FichierExport.Open("ExportTwoGrid.txt",CFile::modeCreate|CFile::modeWrite|CFile::modeNoTruncate))
					{
						FichierExport.SeekToEnd();
						ChaineTemp.Format("%f",TwoHighGradR);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\t");
						ChaineTemp.Format("%f",TwoHighGradG);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\t");
						ChaineTemp.Format("%f",TwoHighGradB);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\t");
						ChaineTemp.Format("%f",TwoHighGrad);
						FichierExport.WriteString(ChaineTemp);
						FichierExport.WriteString("\n");					
						FichierExport.Close();
					}
					else
						int stop=0;*/

//				TwoHighGrad=Vlc*Vrc;

			//	Temp.Format("%f;%f\n",pti,ptj);
			//	Fic.WriteString(Temp);

//-------------------------------------------------------------------------------------------------------
				// pas de contours
				//float NoEdge=((1.0-Vlc)+(1.0-Vrc))/2;
				float NoEdge=((1.0-Vlc)*(1.0-Vrc)); 

//-------------------------------------------------------------------------------------------------------
				// le fond?
				float Background=SimilarBorder*NoEdge;
				

				NonBackground=1.0-Background;

				feature[0]=SimilarBorder;
				feature[1]=TwoHighGrad;
				feature[2]=OppositeAngle;
				feature[3]=NonBackground;
				
				/*
				//feature[3]=fabs(pti-ptj);
				//feature[3]=feature[2]*(fabs(gpi-gpj));
				if ((IptFinal-1)-(IptInit+1)>1)
				{
				pti=I_Qg[Tab[IptInit+1].y][Tab[IptInit+1].x];
				ptj=I_Qg[Tab[IptFinal-1].y][Tab[IptFinal-1].x];
				feature[3]=(pti-ptj)*(pti-ptj);
				pti=I_Qr[Tab[IptInit+1].y][Tab[IptInit+1].x];
				ptj=I_Qr[Tab[IptFinal-1].y][Tab[IptFinal-1].x];
				feature[3]+=(pti-ptj)*(pti-ptj);
				pti=I_Qb[Tab[IptInit+1].y][Tab[IptInit+1].x];
				ptj=I_Qb[Tab[IptFinal-1].y][Tab[IptFinal-1].x];
				feature[3]+=(pti-ptj)*(pti-ptj);
				feature[3]=sqrt(feature[3])/rac3;


				pti=I_i[Tab[IptInit+1].y][Tab[IptInit+1].x];
				ptj=I_i[Tab[IptFinal-1].y][Tab[IptFinal-1].x];


				} 
				else feature[3]=0.0;*/


			/*	if (I_g[Tab[IptInit].y][Tab[IptInit].x]==0.0 || I_g[Tab[IptFinal].y][Tab[IptFinal].x]==0.0)	feature[5]=1.0;
				else*/ 
				//feature[5]=sig;
				//feature[5]=(1-fabs(Vrc-Vlc))*(Vrc+Vlc)/2.0;

				//GradientProduct=I_g[l][c]*I_g[Tab[IptFinal].y][Tab[IptFinal].x];

				//Qpipj=GradientProduct;//*exp(-(m+sig));
				//Qpipj=(1-m);//*(1-sig);
				//BYTE bi=sig*100;
				/*if (I_i[l][c]==1.0) */
				//Histo_Sig[bi]++;
		
				if (Mode==3)
				{
					/*for(int par=0;par<NPAR;par++)
 					{ 
						x[par].index=par+1;
						x[par].value=(double)feature[par];
					}
					x[NPAR].index = -1;
					int v = svm_predict(model,x);
					Qpipj=1;
					if (v==1) 
					{
						Qpipj=0 ;
						break;
					}*/
					Qpipj=0.0;

					Qpipj=TwoHighGrad;//*OppositeAngle*SimilarBorder;//*NoEdge;//ajoute de Noedge mais toujours autant de bruit
					Qpipj=1.0-Qpipj;
					if (Qpipj>1.0) Qpipj=1.0;

				}
				if (Mode==4)
				{
					Qpipj=TwoHighGrad;//*OppositeAngle*SimilarBorder*NoEdge;
					Qpipj=1.0-Qpipj;
					//Qpipj/=NPAR;
					if (Qpipj>1.0) Qpipj=1.0;
				}

				
				if (Qpipj<*pQ) // min search
				{
					*pW=w;*pQ=Qpipj;*piw=IptInit;AxisFin=Axis;
					for(int par=0;par<NPAR;par++) FeatureFin[par]=feature[par];
				}
			} // if opposite grad
		} // if grad pj!=0
		kw++;
	} // for w

	//Fic.Close();
	//	if (Qpipj<*pQ && Mode==3) {*pW=w;*pQ=Qpipj;*piw=IptInit;} // for SVM TEST

free(x);


}

void CCurviDetect::WriteSample(CString Nom, float *feature, int t)
{
	CStdioFile Fic;
	CString ChaineVal,Texte,ChainePar,ChaineDonnee,ValClasse,ChaineTemp;
    int nbd;
	int par;

	ChainePar.Format("%d\n",t);
	ChaineDonnee.Format("%010d\n",1);
	Fic.Open(Nom,CFile::modeCreate|CFile::modeReadWrite|CFile::modeNoTruncate);
	if (Fic.GetLength()!=0) // le fichier existait
	{
		Fic.ReadString(ChaineTemp);
		Fic.ReadString(ChaineDonnee);
		nbd=atoi(ChaineDonnee);
		nbd+=1;//2;
		ChaineDonnee.Format("%010d\n",nbd);
		Fic.SeekToBegin();
		Fic.WriteString(ChainePar);
		Fic.WriteString(ChaineDonnee);
		Fic.SeekToEnd();
	}
	else
	{
		nbd=1;//2;
		ChaineDonnee.Format("%010d\n",nbd);
		Fic.WriteString(ChainePar);
		Fic.WriteString(ChaineDonnee);
	}

	for(par=0;par<t;par++)
 	{ 
		ChaineVal.Format("%f ",feature[par]);
		Fic.WriteString(ChaineVal);
	}//for par
	Fic.WriteString("\n");
	Fic.Close();
}


void CCurviDetect::CrossSection(int l, int c, int Boucle)
{
	POINT *Tab1;
	POINT *Tab2;

	int w;
	int NbPts;
	int iw,Width,IptInit;
	int iw2,Width2;

	float Qpipj=1;
	float Qpipj2=1;
	int Widthmin=0;
	int Widthmin2=0;
	int Oldw;
	int par;

	float SeuilNL;
	
	if (I_Ext[l][c]==0.0 && Mode!=4) 
	{
		I_map[l][c]=0.0;
		return;
	}


	CString Aff, Temp;

	QpipjMin=1;
	QpipjMin2=1;
	SeuilNL=0.99999;

	NbPts=2*dw+wmax+1;
	m_TailleTabCoupe = NbPts;
	for( par=0;par<NPAR;par++) featureQmax[par]=0;;
	Tab1= (POINT*) malloc(NbPts*sizeof(POINT));
	Tab2= (POINT*) malloc(NbPts*sizeof(POINT));
	float gx,gy;
	float Alpha,Angle,AngleInit,AngleFinal;
	float AlphaMin,AlphaMin2,AlphaMinSec;
	float Pas;
	int Dir;
	float OldQ=Qpipj;
	float sq,sw;
	int WidthMinSec;

	Angle=AngleGradient(l,c);//Argument du gradient
	AngleInit=Angle-(pi/16.0);//on parcours + ou - pi/16
	AngleFinal=Angle+(pi/16.0);
	AlphaMin=AlphaMinSec=Angle;//recopie de la valeur du gradient 
	Pas=fabs(AngleInit-AngleFinal)/5;

	// init step
	WidthMinSec=wmax;
	gx=cos(AngleInit);//cos et sinus de l'angle initial
	gy=sin(AngleInit);
	if (gx>=0) Dir=1 ; else Dir=-1;//direction en fonction du cosinus

	RemplitTableauIJCoupe(Tab1,c,l,gx,gy,wmax,dw,Dir,&IptInit); //remplie les tableaux avec les coordonnées des directions

	Quality(featureQmax,&Qpipj,&Width,&iw,Tab1,IptInit,l,c); //fonction de cout en quelque sorte

	if (Qpipj<QpipjMin) {QpipjMin=Qpipj;AlphaMin=AngleInit;}
	if (Width<WidthMinSec) {WidthMinSec=Width;AlphaMinSec=AngleInit;}
	OldQ=Qpipj;
	Oldw=Width;
	sq=0;
	sw=0;
	AxisFin.x=0;
	AxisFin.y=0;

	for (Alpha=AngleInit+Pas;Alpha<=AngleFinal;Alpha+=Pas) //parcours pour tous les angles
	{
		gx=cos(Alpha);
		gy=sin(Alpha);
		if (gx>=0) Dir=1 ; else Dir=-1;
		RemplitTableauIJCoupe(Tab1,c,l,gx,gy,wmax,dw,Dir,&IptInit); //rempli tableau des droites
// for testing
		/*
		for (w=0;w<NbPts;w++)
			{
				pImage[Tab1[w].y*Largeur+Tab1[w].x].r=255;
				pImage[Tab1[w].y*Largeur+Tab1[w].x].g=255;
				pImage[Tab1[w].y*Largeur+Tab1[w].x].b=255;
			}*/
		Quality(featureQmax,&Qpipj,&Width,&iw,Tab1,IptInit,l,c); //fonction de cout
		if (Qpipj<QpipjMin) {QpipjMin=Qpipj;AlphaMin=Alpha;}
		if (Width<WidthMinSec) {WidthMinSec=Width;AlphaMinSec=Alpha;}
		sq=sq+fabs(OldQ-Qpipj);
		sw=sw+fabs((double)(Oldw-Width));
		OldQ=Qpipj;
		Oldw=Width;
		if (Qpipj==0 && Mode==3) break;
	}

	if (sq==0.0) AlphaMin=AlphaMinSec; // if QPiPj constant
	if (sq==0.0 && sw==0) AlphaMin=Angle; // if Width and Qpipj constant
	
	gx=cos(AlphaMin);
	gy=sin(AlphaMin);
	if (gx>=0) Dir=1 ; else Dir=-1;
	RemplitTableauIJCoupe(Tab1,c,l,gx,gy,wmax,dw,Dir,&IptInit); 
	Quality(featureQmax,&QpipjMin,&Width,&iw,Tab1,IptInit,l,c);

	// non linearisation
	/*
	for(par=0;par<NPAR;par++)
 	{
		if (featureQmax[par]>SeuilNL) {QpipjMin=1.0; break;}
	}*/
	AxisFin1=AxisFin;


	for(par=0;par<NPAR;par++) featureQmax2[par]=0;
	Qpipj=QpipjMin2=1.0;
	AlphaMin2=AlphaMinSec=Angle;
	WidthMinSec=wmax;

	gx=cos(AngleInit);
	gy=sin(AngleInit);
	if (gx>=0) Dir=-1 ; else Dir=1;
	RemplitTableauIJCoupe(Tab2,c,l,gx,gy,wmax,dw,Dir,&IptInit); 
	Quality(featureQmax2,&Qpipj,&Width2,&iw,Tab2,IptInit,l,c);
	if (Qpipj<QpipjMin2) {QpipjMin2=Qpipj;AlphaMin2=AngleInit;}
	if (Width2<WidthMinSec) {WidthMinSec=Width2;AlphaMinSec=AngleInit;}
	OldQ=Qpipj;
	Oldw=Width2;
	sq=sw=0;
	for (Alpha=AngleInit+Pas;Alpha<AngleFinal;Alpha+=Pas)
	{
		gx=cos(Alpha);
		gy=sin(Alpha);
		if (gx>=0) Dir=-1 ; else Dir=1;
		RemplitTableauIJCoupe(Tab2,c,l,gx,gy,wmax,dw,Dir,&IptInit); 
		Quality(featureQmax2,&Qpipj,&Width2,&iw,Tab2,IptInit,l,c);
		if (Qpipj<QpipjMin2) {QpipjMin2=Qpipj;AlphaMin2=Alpha;}
		if (Width2<WidthMinSec) {WidthMinSec=Width2;AlphaMinSec=Alpha;}
		sq=sq+fabs(OldQ-Qpipj);
		sw=sw+fabs((double)(Oldw-Width2));
		OldQ=Qpipj;
		Oldw=Width2;
		if (Qpipj==0 && Mode==3) break;
	}

	if (sq==0.0) AlphaMin2=AlphaMinSec; // if QPiPj constant
	if (sq==0.0 && sw==0) AlphaMin2=Angle; // if Width and Qpipj constant


	QpipjMin2=1.0;
	gx=cos(AlphaMin2);
	gy=sin(AlphaMin2);
	if (gx>=0) Dir=-1 ; else Dir=1;
	RemplitTableauIJCoupe(Tab2,c,l,gx,gy,wmax,dw,Dir,&IptInit); // ajouter faire varier l'angle
	Quality(featureQmax2,&QpipjMin2,&Width2,&iw2,Tab2,IptInit,l,c);


	// non linearisation
	/*
	for(par=0;par<NPAR;par++)
 	{
		if (featureQmax2[par]>SeuilNL) {QpipjMin2=1.0; break;}
	}*/

	if (Boucle)
	{
		if (QpipjMin<QpipjMin2) 
		{
/*			for (w=iw;w<iw+Width;w++)
			{
				I_map[Tab1[w].y][Tab1[w].x]+=(1-QpipjMin);
				//else I_map[Tab1[w].y][Tab1[w].x]=0;
			}*/
//			if (I_map[AxisFin1.y][AxisFin1.x]<(1.0-QpipjMin)) I_map[AxisFin1.y][AxisFin1.x]=1.0-QpipjMin;
			if (I_map[AxisFin1.y][AxisFin1.x]<(1.0-QpipjMin)) 
			{		
				Cross[AxisFin1.y][AxisFin1.x].Width=Width;
				Cross[AxisFin1.y][AxisFin1.x].Angle=AlphaMin;//+pis2;
				I_map[AxisFin1.y][AxisFin1.x]=1.0-QpipjMin;
			}
		}
		else
		{
/*			for (w=iw2;w<iw2+Width2;w++)
			{
				I_map[Tab2[w].y][Tab2[w].x]+=(1-QpipjMin2);
				//else I_map[Tab2[w].y][Tab2[w].x]=0;
			}*/
	//		if (I_map[l][c]<(1.0-QpipjMin2)) I_map[l][c]=1.0-QpipjMin2;
			if (I_map[AxisFin.y][AxisFin.x]<(1.0-QpipjMin2)) 
			{
				Cross[AxisFin1.y][AxisFin1.x].Width=Width2;
				Cross[AxisFin1.y][AxisFin1.x].Angle=AlphaMin2;//+pis2;
				I_map[AxisFin.y][AxisFin.x]=1.0-QpipjMin2;
			}

		}
	}
	else
	{
		if (QpipjMin<QpipjMin2) 
		{
			for (w=iw;w<iw+Width;w++)
			{
				pImage[Tab1[w].y*Largeur+Tab1[w].x].r=255;
				pImage[Tab1[w].y*Largeur+Tab1[w].x].g=255;
				pImage[Tab1[w].y*Largeur+Tab1[w].x].b=255;
			}
			//	I_map[l][c]=QpipjMax;

		}
		else if (QpipjMin<1.0)
		{
			for (w=iw2;w<iw2+Width2;w++)
			{
				pImage[Tab2[w].y*Largeur+Tab2[w].x].r=127;
				pImage[Tab2[w].y*Largeur+Tab2[w].x].g=127;
				pImage[Tab2[w].y*Largeur+Tab2[w].x].b=127;
			}
		}
	}

 	free(Tab1);
	free(Tab2);
}

void CCurviDetect::Init()
{
 int i,p,l,c;

  wmin=16;
  m_MaxWidth=12;//largeur maximale de la forme curvilinéaire
  wmax=(Largeur*m_MaxWidth)/100;//pourcentage de l'image
  if (wmax>=Largeur) wmax=(Largeur*10)/100; //on prend 10%
  if (wmax<wmin) wmax=wmin+2;//sinon largeur mini + 2
  dw=17;//marge à droite et à gauche  

  int NbW=wmax*2;
  int Taille=wmin;
  Gauss=(float **) malloc (sizeof(float*)*NbW);// tableau + init
  for (i=0;i<NbW;i++)
  {
	  Gauss[i]=(float *)malloc(sizeof(float)*Taille);
	  Taille++;
  }

  int w;
  double x,arg;
  double deuxpi=6.28318;
  double s=0;

 // CStdioFile Fic("test.xls",CFile::modeCreate | CFile ::modeWrite);
 // CString Temp;
  int wfin;
  wfin=wmin-1;//largeur de fin
  for (i=0;i<NbW;i++)//on boucle sur deux fois la largeurmax
  {	  
	  s=0.0;
	  double sigma=0.5*wfin;
	  double unsdeuxsig=1.0/(2.0*sigma);
	  double norm=1.0/(sqrt(deuxpi)*sigma);
	  for (w=0;w<wfin;w++)
	  {
		x=double(w-wfin/2);
		arg=(x*x)*unsdeuxsig;
		Gauss[i][w]=(float)(norm*exp(-arg));
		s=s+Gauss[i][w];
	  }
	  for (w=0;w<wfin;w++)
	  {
		Gauss[i][w]/=s;
//		Temp.Format("%f;",Gauss[i][w]);
//		Fic.WriteString(Temp);
	  }
	  wfin++;
//	  Fic.WriteString("\n");

  }
//  Fic.Close();

  I_i=(float **) malloc(sizeof(float*)*Hauteur);
  I_g=(float **) malloc(sizeof(float)*Hauteur);
  I_x=(float **) malloc(sizeof(float*)*Hauteur);
  I_y=(float **) malloc(sizeof(float*)*Hauteur);
  I_Qr=(float **) malloc(sizeof(float*)*Hauteur);
  I_Qb=(float **) malloc(sizeof(float*)*Hauteur);
  I_Qg=(float **) malloc(sizeof(float*)*Hauteur);
  I_map=(float **) malloc(sizeof(float*)*Hauteur);
  I_Ext=(float **) malloc(sizeof(float*)*Hauteur);
  Cross=(SCross **) malloc(sizeof(SCross*)*Hauteur);

  for (i=0;i<Hauteur;i++)
  {
    I_i[i]=(float *) malloc(Largeur*sizeof(float));
    I_g[i]=(float *) malloc(Largeur*sizeof(float));
    I_x[i]=(float *) malloc(Largeur*sizeof(float));
    I_y[i]=(float *) malloc(Largeur*sizeof(float));
    I_Qr[i]=(float *) malloc(Largeur*sizeof(float));
    I_Qb[i]=(float *) malloc(Largeur*sizeof(float));
    I_Qg[i]=(float *) malloc(Largeur*sizeof(float));
    I_map[i]=(float *) malloc(Largeur*sizeof(float));
    I_Ext[i]=(float *) malloc(Largeur*sizeof(float));
    Cross[i]=(SCross *) malloc(Largeur*sizeof(SCross));
 }


  /*CString ChaineModel("c:\\svm.model");
  model=svm_load_model(ChaineModel);*/

  // compute image intensity
  p=0;
  float H,fr,fg,fb,fi,fi2,num,denom;
  for (l=0;l<Hauteur;l++)
  {
	  for (c=0;c<Largeur;c++)
	  {
		  I_g[l][c]=0.0;
		  I_x[l][c]=0.0;
		  I_y[l][c]=0.0;
		  I_Ext[l][c]=0.0;
		  I_map[l][c]=0.0;
		  fr=(float)pImage[p].r;
		  fg=(float)pImage[p].g;
		  fb=(float)pImage[p].b;
		  fi2=(0.30*fr+0.59*fg+0.11*fb);
		  fi=fr+fg+fb;
//		  I_i[l][c]=fi/(3.0*255.0);
		  I_i[l][c]=fi2/(255.0);
//--------------------------------------------------
//		modif du 10/04/2006 plus de normalisation
//--------------------------------------------------
//      modif du 12/04/2006 retour avec la normalisation
		  if (fi>=1.0)
		  {
				I_Qr[l][c]=fr/fi;
				I_Qb[l][c]=fb/fi;
				I_Qg[l][c]=fg/fi;
		  }
		  else
		  {
				I_Qr[l][c]=0.0;
				I_Qb[l][c]=0.0;
				I_Qg[l][c]=0.0;
		  }
		  // computation of HUE REMOVED
/*		  H=Hue(pImage[p]);
		  if (I_i[l][c]<10.0) I_Ext[l][c]=0.0;
		  else I_Ext[l][c]=H;

		  /*I_Qr[l][c]=(float)pImage[p].r;
		  I_Qb[l][c]=(float)pImage[p].g;
		  I_Qg[l][c]=(float)pImage[p].b;*/
	      p++;
	  }
  }

// dynamique
  float max,min;
  max=min=I_i[0][0];
  for (l=0;l<Hauteur;l++)
   for (c=0;c<Largeur;c++)
   {
	if (I_i[l][c]>max) max=I_i[l][c];
	if (I_i[l][c]<min) min=I_i[l][c];

   }
   if (max!=min)
	for (l=0;l<Hauteur;l++)
	 for (c=0;c<Largeur;c++)
	 {
	  I_i[l][c]=(I_i[l][c]-min)/(max-min);
	 }

  // compute image gradient X, Y and G
/*
  for (c=1;c<Largeur-1;c++)
  {
	  for (l=1;l<Hauteur-1;l++)
	  {
		  I_map[l][c]=(I_i[l-1][c]+2.0*I_i[l][c]+I_i[l+1][c])/3; // lisse vertical
		  I_Ext[l][c]=(I_i[l][c-1]+2.0*I_i[l][c]+I_i[l][c+1])/3; // lisse horizontal
	  }
  }
  for (c=1;c<Largeur-1;c++)
  {
	  for (l=1;l<Hauteur-1;l++)
	  {
		  I_y[l][c]=I_Ext[l-1][c]-I_Ext[l+1][c];
		  I_x[l][c]=I_map[l][c-1]-I_map[l][c+1];
		  I_g[l][c]=fabs(I_x[l][c])+fabs(I_y[l][c]);
	  }
  }
*/
  for (c=1;c<Largeur-1;c++)
  {
	  for (l=1;l<Hauteur-1;l++)
	  {
		  I_x[l][c]=(+I_i[l-1][c-1]+2.0*I_i[l][c-1]+I_i[l+1][c]-
			  I_i[l-1][c+1]-2.0*I_i[l][c+1]-I_i[l+1][c+1])/3.0; 
		  I_y[l][c]=(+I_i[l-1][c-1]+2.0*I_i[l-1][c]+I_i[l][c+1]-
			  I_i[l+1][c-1]-2.0*I_i[l+1][c]-I_i[l+1][c+1])/3.0; 
 		  I_g[l][c]=sqrt(I_x[l][c]*I_x[l][c]+I_y[l][c]*I_y[l][c]);

	  }
  }

  // dynamique for gradient
  max=min=I_g[2][2];
  for (l=2;l<Hauteur-2;l++)
   for (c=2;c<Largeur-2;c++)
   {
	if (I_g[l][c]>max) max=I_g[l][c];
	if (I_g[l][c]<min) min=I_g[l][c];

   }

  if (max!=min)
	for (l=2;l<Hauteur-2;l++)
	 for (c=2;c<Largeur-2;c++)
	 {
	  I_g[l][c]=(I_g[l][c]-min)/(max-min);
	 }
 
	 extrema(0.05);

}




void CCurviDetect::extrema(float seuil)
{
	//calcul les extrema locaux
  int i,j;
  int m=0,n=0;
  int q;
  float tangente;
  float a=0.414;
  float b=2.414;
  int d=2;


  for (i=0;i<Hauteur;i++)
    for (j=0;j<Largeur;j++)
      {
		I_Ext[i][j]=I_g[i][j];
      }

  for (i=d;i<Hauteur-d;i++)
  {
    for (j=d;j<Largeur-d;j++)
      {
	
    	if (I_x[i][j]==0.0) 
		{
			if (I_y[i][j]>=0.0) {m=1;n=0;} else {m=-1;n=0;}
		}
		else 
		{
			tangente=I_y[i][j]/I_x[i][j];
			if (tangente>-a && tangente<a) {m=0; if (I_x[i][j]>=0) n=1; else n=-1;}
			if (tangente>a && tangente<b && I_x[i][j]>=0 && I_y[i][j]>=0) {m=1;n=1;}
			if (tangente>a && tangente<b && I_x[i][j]<0 && I_y[i][j]<0) {m=-1;n=-1;}
			if (fabs(tangente)>b) {n=0; if(I_y[i][j]>=0) m=1;else m=-1;}
			if (tangente>-b && tangente<-a && I_x[i][j]<=0 && I_y[i][j]>=0) {m=1;n=-1;}
			if (tangente>-b && tangente<-a && I_x[i][j]>=0 && I_y[i][j]<=0) {m=-1;n=1;}

		}

	
		for (q=1;q<d+1;q++)
		{  
			if (I_g[i][j]<I_g[i+q*m][j+q*n] 
				|| I_g[i][j]<I_g[i-q*m][j-q*n])
			  I_Ext[i][j]=0; 
		}

		if (I_Ext[i][j]>seuil) 
		{
			I_Ext[i][j]=1.0;
		}
		else I_Ext[i][j]=0;

	} // fin for j
  }// for i


}

void CCurviDetect::Info()
{
	int par;
	CString Temp,Aff;
	Aff="";

	Temp.Format("SimilarBorder=%1.3f \n",featureQmax[0]);
	Aff+=Temp;
	Temp.Format("TwoHighGrad=%1.3f \n",featureQmax[1]);
	Aff+=Temp;
	Temp.Format("OppositeAngle=%1.3f \n",featureQmax[2]);
	Aff+=Temp;
	Temp.Format("NonBackground=%1.3f \n",featureQmax[3]);
	Aff+=Temp;
	Temp.Format("Qmin1=%1.3f \n\n",QpipjMin);
	Aff+=Temp;
	Aff+="\n\n";
/*
	for(par=0;par<NPAR;par++)
 	{
		Temp.Format("%1.3f \n",featureQmax[par]);
		Aff+=Temp;
	}*/
	Temp.Format("SimilarBorder=%1.3f \n",featureQmax2[0]);
	Aff+=Temp;
	Temp.Format("TwoHighGrad=%1.3f \n",featureQmax2[1]);
	Aff+=Temp;
	Temp.Format("OppositeAngle=%1.3f \n",featureQmax2[2]);
	Aff+=Temp;
	Temp.Format("NonBackground=%1.3f \n",featureQmax2[3]);
	Aff+=Temp;
	Temp.Format("Qmin2=%1.3f \n",QpipjMin2);
	Aff+=Temp;
/*
	for(par=0;par<NPAR;par++)
 	{
		Temp.Format("%1.3f \n",featureQmax2[par]);
		Aff+=Temp;
	}*/
	int cla=MessageBox(NULL,Aff,"",MB_OKCANCEL);
	if (cla==IDOK) // 	// learning step : this is a curvi region
	{
				if (QpipjMin<QpipjMin2) WriteSample("C:\\johel\\recherche\\line detection\\Curvi.dat",featureQmax,NPAR);
				else WriteSample("C:\\johel\\recherche\\line detection\\Curvi.dat",featureQmax2,NPAR);
	}
	else 
	{

				if (QpipjMin<QpipjMin2) WriteSample("C:\\johel\\recherche\\line detection\\Not Curvi.dat",featureQmax,NPAR);
				else WriteSample("C:\\johel\\recherche\\line detection\\Not Curvi.dat",featureQmax2,NPAR);
	}


}

void CCurviDetect::ComputeCoefGauss(int nbiter)
{
	double x,arg;
	double deuxpi=6.28318;
	int w;
	float s=0.0;
	double sigma=0.5*nbiter;
	double unsdeuxsig=1.0/(2.0*sigma);
	double norm=1.0/(sqrt(deuxpi)*sigma);
	
	CoGauss=(float*) malloc(sizeof(float)*nbiter);
	for (w=0;w<nbiter;w++)
	{
		x=double(w-nbiter/2);
		arg=(x*x)*unsdeuxsig;
		CoGauss[w]=(float)(norm*exp(-arg));
		s=s+CoGauss[w];
	}
	for (w=0;w<nbiter;w++)
	{
		CoGauss[w]/=s;
	}
}

void CCurviDetect::Hough()
{
	int l,c;
	CProgressDlg Progress;
	Progress.Create(NULL);
	Progress.SetRange(0,Largeur*Hauteur);
	int TailleRho,i,j,pos=0;
	float RhoInit,RhoFin,Rho,PasTeta,Teta,AngleInit,AngleFin,PasRho,Teta2;
	float RhoTest,SeuilRho=1;
	
	float **AccuHough;
	float x,y;
	float deuxpi=2.0*pi;

	RhoFin=sqrt((float)(Largeur*Largeur+Hauteur*Hauteur));
	RhoInit=-RhoFin;
	PasRho=1.0;
	TailleRho=(RhoFin-RhoInit)/PasRho;
	
	float TailleAngle=360;

	PasTeta=(deuxpi)/TailleAngle;

	AccuHough=(float **) malloc(sizeof(float*)*TailleRho);
	for (i=0;i<TailleRho;i++) 
	{
		AccuHough[i]=(float *) malloc(sizeof(float)*TailleAngle);
		for (j=0;j<TailleAngle;j++) AccuHough[i][j]=0.0;
	}



	for (l=dw+1;l<Hauteur-dw;l+=1)
	{
		Progress.SetPos(pos);
		for (c=dw+1;c<Largeur-dw;c+=1)
		{
			pos++;
			if (I_map[l][c]!=0.0)
			{
				AngleInit=Cross[l][c].Angle-(pi/16.0);
				AngleFin=Cross[l][c].Angle+(pi/16.0);
				if (AngleFin>deuxpi) AngleFin-=deuxpi;
				if (AngleInit<0.0)	 
				//AngleInit=0;
				//AngleFin=2.0*pi;
				x=(float)(c);
				y=(float)(l);
				i=0;
				for (Rho=RhoInit;Rho<RhoFin;Rho+=PasRho)
				{
					for(Teta=AngleInit;Teta<AngleFin;Teta+=PasTeta)
					{
						Teta2=Teta;
						if (Teta<0.0) Teta2=deuxpi+Teta;
						if (Teta>deuxpi) Teta2=Teta-deuxpi;
						j=((360*Teta2)/(deuxpi))-1;
						if (j>TailleAngle) 
						{
							MessageBox(NULL,"","rate",MB_OK);
						}
						RhoTest=x*cos(Teta)+y*sin(Teta);
						if (fabs(RhoTest-Rho)<SeuilRho) AccuHough[i][j]++;//=I_map[l][c];
					//	j++;
					}
					i++;
				}
			}
		} // l
	}// c

  float max;
  max=AccuHough[0][0];
  for (l=0;l<TailleRho;l++)
   for (c=0;c<TailleAngle;c++)
   {
	if (AccuHough[l][c]>max) max=AccuHough[l][c];

   }
	FILE *fp=fopen("test.raw","wb");
	for (l=0;l<TailleRho;l++)
	 for (c=0;c<TailleAngle;c++)
	 {
	  AccuHough[l][c]=(AccuHough[l][c]*255)/(max);
	  BYTE a=(BYTE)AccuHough[l][c];
	  fwrite(&a,1,1,fp);
	 }
	fclose(fp);

	for (i=0;i<TailleRho;i++) free(AccuHough[i]);
	free(AccuHough);

}


//-----------------------------------------------
// Calcul le module de 2 tableaux de complexes
// -----> résultat de la FFT Cima
// dont la taille est passé en paramètres 
// -----> Taille DeuxPuissanceNFFT
// Le tableau résultat doit être alloué avant
//-----------------------------------------------

void CCurviDetect::ModuleComplex(double * resultat, double * Reel, double * Imaginaire, int DeuxPuissanceNFFT)
{
	int i=0;
	int NombreValeur = (int)pow(2.0,DeuxPuissanceNFFT);

	for(i=0;i<NombreValeur;i++)
	{
		resultat[i] = sqrt((Reel[i]*Reel[i])+(Imaginaire[i]*Imaginaire[i]));
	}
}



double CCurviDetect::DistanceEuclPoint(POINT Val1, POINT Val2)
{
	double rtn=0;

	rtn = sqrt((double)((Val2.x-Val1.x)*(Val2.x-Val1.x))+((Val2.y-Val1.y)*(Val2.y-Val1.y)));


	return rtn;
}

void CCurviDetect::NormalisationSpectre(double *Spectre, int NbValeur)
{
int i=0;
double Temp = Spectre[0];
for(i=0;i<NbValeur;i++)
{
	if (Spectre[i]==0.0f)
	{
	Spectre[i]=0.0f;
	}
	else
	{
	Spectre[i]/=Temp;
	}
}
}

void CCurviDetect::NormalisationSpectreComplexe(double *SpectreRe, double *SpectreIm, int NbValeur)
{
//-----------------------
//			init
//-----------------------
long i = 0;
double MaxiRe = 0.0f;
double MaxiIm = 0.0f;


//-----------------------------
//			Normalisation
//-----------------------------
for(i=0;i<NbValeur;i++)
{
	if (SpectreRe[i]>=MaxiRe)
	{
	SpectreRe[i]=MaxiRe;
	}
	if (SpectreIm[i]>=MaxiIm)
	{
	SpectreIm[i]=MaxiIm;
	}
}

for(i=0;i<NbValeur;i++)
{
	if (SpectreRe[i]==0.0f)
	{
	SpectreRe[i]=0.0f;
	}
	else
	{
	SpectreRe[i]/=MaxiRe;
	}


	if (SpectreIm[i]==0.0f)
	{
	SpectreIm[i]=0.0f;
	}
	else
	{
	SpectreIm[i]/=MaxiIm;
	}
}

}
