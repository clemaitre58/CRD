// CurviDetect1D.cpp: implementation of the CCurviDetect1D class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "winseg.h"
#include "CurviDetect1D.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define DISTGAUCHE

const double gauss_table_sigma_1[91]=
{1.59838260850225e-005,     2.4942603897572e-005,    3.85354016158044e-005,
 5.89433811294151e-005,    8.92621317398461e-005,     0.000133830937277948,
  0.000198656527551514,     0.000291948477937472,     0.000424782528909873,
  0.000611905183311721,      0.00087268733469374,       0.0012322257196092,
   0.00172257809715151,       0.0023841008765529,      0.00326683642435735,
   0.00443187197401463,      0.00595256406662203,      0.00791549366575433,
    0.0104209902176859,       0.0135830414480152,       0.0175283936833893,
    0.0223946493561441,       0.0283271883432851,       0.0354747814481413,
     0.043983829820749,       0.0539912535580442,       0.0656161636234438,
    0.0789505780421638,       0.0940495773920863,        0.110921424393941,
     0.129518284250705,        0.149728261667037,        0.171369503134903,
     0.194187087380145,        0.217853335251246,        0.241972010964938,
     0.266086664550264,        0.289693092916724,        0.312255593475686,
      0.33322637449202,        0.352067198531885,        0.368272098224365,
     0.381389843122144,        0.391044772967574,        0.396954657889072,
     0.398944401391977,        0.396954657889072,        0.391044772967574,
     0.381389843122144,        0.368272098224365,        0.352067198531885,
      0.33322637449202,        0.312255593475686,        0.289693092916724,
     0.266086664550264,        0.241972010964938,        0.217853335251246,
     0.194187087380145,        0.171369503134903,        0.149728261667037,
     0.129518284250705,        0.110921424393941,       0.0940495773920863,
    0.0789505780421638,       0.0656161636234438,       0.0539912535580442,
     0.043983829820749,       0.0354747814481413,       0.0283271883432851,
    0.0223946493561441,       0.0175283936833893,       0.0135830414480152,
    0.0104209902176859,      0.00791549366575433,      0.00595256406662203,
   0.00443187197401463,      0.00326683642435735,       0.0023841008765529,
   0.00172257809715151,       0.0012322257196092,      0.00087268733469374,
  0.000611905183311721,     0.000424782528909873,     0.000291948477937472,
  0.000198656527551514,     0.000133830937277948,    8.92621317398461e-005,
 5.89433811294151e-005,    3.85354016158044e-005,     2.4942603897572e-005,
 1.59838260850225e-005
};

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCurviDetect1D::CCurviDetect1D()
: Fenetre(0)
, DemiLargeurFourier(0)
, DemiLongueurSignal(0)
, IptInit(0)
, IptFinal(0)
{

}

CCurviDetect1D::~CCurviDetect1D()
{


if(FGaucheRE!=NULL) free(FGaucheRE); 
if(FDroiteRE != NULL) free(FDroiteRE);
if(FCentreGRE != NULL) free(FCentreGRE);
if(FCentreDRE != NULL) free(FCentreDRE);


if(FGaucheIM != NULL) free(FGaucheIM); 
if(FDroiteIM != NULL) free(FDroiteIM); 
if(FCentreGIM != NULL) free(FCentreGIM); 
if(FCentreDIM != NULL) free(FCentreDIM);

if(FSDFixeRE != NULL) free(FSDFixeRE);
if(FSDMobileRE != NULL) free(FSDMobileRE);

if(FSDFixeIM != NULL) free(FSDFixeIM);
if(FSDMobileIM != NULL) free(FSDMobileIM);



if(TabDistanceGauche != NULL) free(TabDistanceGauche);
if(TabDistanceDroite != NULL) free(TabDistanceDroite);
if(TabHomogSection != NULL) free(TabHomogSection);
if(TabBorSim != NULL) free(TabBorSim);
//if(TabDisGlis!= NULL) free(TabDisGlis);
if(MuLTabDis != NULL) free(MuLTabDis);
//if(TabMaxVisu != NULL) free(TabMaxVisu); //---> Pas bon

//if(TabStoreTemp != NULL) free(TabStoreTemp);
if(TabStoreLarg != NULL) free(TabStoreLarg);
//if(TabLargeur != NULL) free(TabLargeur);
if(Fenetre != NULL) free(Fenetre);

}

void CCurviDetect1D::Quality(CString PathSignal, CString NFichierVisu, CString NFichierTotale, int LargeurFourier, int LongueurSignal)
{

int i=0,j=0,k=0,l=0;

//int LargeurPorte = 48;

int DemiLargeurFourier = LargeurFourier/2;



//int DemiLargeurPorte = LargeurPorte/2;
int DemiLongueurSignal = LongueurSignal/2;

double ValDist = 0.0f;

//int CondGauche = DemiLongueurSignal-DemiLargeurPorte;
//int CondDroite = DemiLongueurSignal+DemiLargeurPorte;

int AccuVisu=0;
double MaxVisu = 0.0f;
double MaxIter = 0.0f;

CStdioFile FichierVisu;
CString ChaineVisu;


int IptInit;
int IptFinal;
int LargeurGabaritMin;
int LargeurGabaritMaxi;
int IptInitMax = LongueurSignal-LargeurFourier;//-LargeurPortemodif cédric 01/02/2007
int IptFinalMax = LongueurSignal-LargeurFourier-1; //attention 
int IptInitMin = LargeurFourier;
int IptFinalMin = IptInitMin+LargeurFourier ;//+ LargeurPorte; //modif 01/02/2007

int cpt=0;
int cptTabG=0;
int cptTabD=0;
int cptSection=0;
int cptBord=0;
double max = 0.0f;
int PosMax=0;
int cpti=0;

int * SignalPorte = (int *) malloc(sizeof(int)*LongueurSignal);
int RacineLargeurFourier = 0;//(int)sqrt((double)LargeurFourier);

while((1<<RacineLargeurFourier)<LargeurFourier) RacineLargeurFourier++;

double * FGaucheRE = (double *) malloc(sizeof(double)*LargeurFourier);
double * FDroiteRE = (double *) malloc(sizeof(double)*LargeurFourier);
double * FCentreGRE = (double *) malloc(sizeof(double)*LargeurFourier);
double * FCentreDRE = (double *) malloc(sizeof(double)*LargeurFourier);


double * FGaucheIM = (double *) malloc(sizeof(double)*LargeurFourier);
double * FDroiteIM = (double *) malloc(sizeof(double)*LargeurFourier);
double * FCentreGIM = (double *) malloc(sizeof(double)*LargeurFourier);
double * FCentreDIM = (double *) malloc(sizeof(double)*LargeurFourier);

double * FSDFixeRE = (double *) malloc(sizeof(double)*LargeurFourier);
double * FSDMobileRE = (double *) malloc(sizeof(double)*LargeurFourier);

double * FSDFixeIM = (double *) malloc(sizeof(double)*LargeurFourier);
double * FSDMobileIM = (double *) malloc(sizeof(double)*LargeurFourier);

double MaxGauche = 0.0f,MaxDroite = 0.0f,MaxHomoSec = 0.0f, MaxBordSim = 0.0f;

int TailleTabDist = (IptInitMax-IptInitMin)*(IptFinalMax-IptFinalMin);

double * TabDistanceGauche = (double *) malloc(sizeof(double)*TailleTabDist);
double * TabDistanceDroite = (double *) malloc(sizeof(double)*TailleTabDist);
double * TabHomogSection = (double *) malloc(sizeof(double)*TailleTabDist);
double * TabBorSim = (double *) malloc(sizeof(double)*TailleTabDist);

double * TabDisGlis;
double * MuLTabDis = (double *) malloc(sizeof(double)*TailleTabDist);

double * TabMaxVisu = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

//double * TabStoreTemp = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

double * TabStoreTemp = (double *)malloc(sizeof(double) *(LongueurSignal));
double * TabStoreLarg = (double *)malloc(sizeof(double) *(LongueurSignal));
int * TabLargeur = (int *)malloc(sizeof(int) *(IptInitMax-IptInitMin));

// Allocation du tableau qui contiendra la fenetre de hamming
Fenetre = (double*)malloc(sizeof(double)*LargeurFourier);
// Initialisation de la Fenetre de Haming
InitFenetreHaming(LargeurFourier);

for(i=0;i<LargeurFourier;i++)
{

	FGaucheRE[i]=0.0f;
	FDroiteRE[i]=0.0f;
	FCentreGRE[i]=0.0f;
	FCentreDRE[i]=0.0f;

	FGaucheIM[i]=0.0f;
	FDroiteIM[i]=0.0f;
	FCentreGIM[i]=0.0f;
	FCentreDIM[i]=0.0f;

	FSDFixeRE[i]=0.0f;
	FSDMobileRE[i]=0.0f;
	FSDFixeIM[i]=0.0f;
	FSDMobileIM[i]=0.0f;

}

for (i=0;i<(LongueurSignal);i++)
{
	TabStoreTemp[i] = 0.0f;
	TabStoreLarg[i] = 0.0f;
}

for (i=0;i<(IptInitMax-IptInitMin);i++)
{
	TabMaxVisu[i] = 0.0f;
	TabLargeur[i] = 0;
}

for(i=0;i<LargeurFourier;i++)
{

	TabDistanceGauche[i] = 0.0f;
	TabDistanceDroite[i] = 0.0f;
	

}

for(i=0;i<TailleTabDist;i++)
{
	MuLTabDis[i] = 0.0f;
	TabDistanceGauche[i]=0.0f;
	TabDistanceDroite[i]=0.0f;
	TabHomogSection[i]=0.0f;
	TabBorSim[i]=0.0f;


}

if(Fichier.Open(PathSignal,CFile::modeRead))
{
	
	Fichier.Read(SignalPorte,sizeof(int)*LongueurSignal);
	Fichier.Close();
}
else
{
	MessageBox(NULL,"Ouverture du fichier <<"+PathSignal+">> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
}



//------------------------
//			RAZ 
//------------------------

cptTabG=0;
cptTabD=0;
cptSection=0;
cptBord=0;


//----------------------------------------------------
// recopie des tab pour fourier + calcul distance
//----------------------------------------------------
for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
{
	for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) //idem modif 21/11/2006
	{
	//copie Fgauche 
	
	cpt =0;
	for (i = IptInit-LargeurFourier;i<IptInit;i++)
	{
	FGaucheRE[cpt]= (double)SignalPorte[i];
	cpt++;
	}
	
	cpt =0;
	for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)
	{
	FDroiteRE[cpt]= (double)SignalPorte[i];
	cpt++;
	}

	cpt =0;
	for (i = IptInit+1;i<IptInit+LargeurFourier+1;i++)//modif 09/01/2007
	{
	FCentreGRE[cpt]= (double)SignalPorte[i];
	cpt++;
	}
	
	cpt =0;
	for (i = IptFinal-LargeurFourier;i<IptFinal;i++)//modif 09/01/2007
	{
	FCentreDRE[cpt]= (double)SignalPorte[i];
	cpt++;
	}

	cpt=0;

	/*
	
	//remplissage tableau fixe :
	for(i=IptInit;i<(IptInit+LargeurFourier);i++)
	{
		FSDFixeRE[cpt]=SignalPorte[i];
		cpt++;

	}
	//calcul de la FFT du tableau fixe
	ObjetCIma.FFT(1,RacineLargeurFourier,FSDFixeRE,FSDFixeIM);
	cpt=0;

	TabDisGlis = (double*) malloc(sizeof(double)*(IptFinal-LargeurFourier-IptInit));
	cpti=0;
	//Glissement pour ecart-type :
	for(i=IptInit;i<(IptFinal-LargeurFourier);i++)
	{
		for(j=0;j<LargeurFourier;j++)
		{
			FSDMobileRE[j]=SignalPorte[i+j];
			
		}

		// calcul de la FFT du tableau glissant
		ObjetCIma.FFT(1,RacineLargeurFourier,FSDMobileRE,FSDMobileIM);
		//calcul de la distance entre fixe et mobile
		TabDisGlis[cpti] = ObjetCIma.HyperDistanceComplexe(FSDFixeRE,FSDFixeIM,
																	FSDMobileRE,FSDMobileIM,DemiLargeurFourier);

	cpti++;
	
	}
	TabHomogSection[cptSection] = ObjetCIma.StandardDeviation(TabDisGlis,IptFinal-LargeurFourier-IptInit);
	//
	free(TabDisGlis);*/
	// calcul des FFT

	//Correction de Haming
	CorrectionHaming(FGaucheRE,Fenetre,LargeurFourier);
	CorrectionHaming(FGaucheIM,Fenetre,LargeurFourier);
	CorrectionHaming(FDroiteRE,Fenetre,LargeurFourier);
	CorrectionHaming(FDroiteIM,Fenetre,LargeurFourier);
	CorrectionHaming(FCentreDRE,Fenetre,LargeurFourier);
	CorrectionHaming(FCentreDIM,Fenetre,LargeurFourier);
	CorrectionHaming(FCentreGRE,Fenetre,LargeurFourier);
	CorrectionHaming(FCentreGIM,Fenetre,LargeurFourier);
	//calcul FFT
	ObjetCIma.FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);
	ObjetCIma.FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);
	ObjetCIma.FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);
	ObjetCIma.FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);

					
	TabDistanceGauche[cptTabG] = ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,LargeurFourier);
	TabDistanceDroite[cptTabD] = ObjetCIma.HyperDistanceModuleFourier(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,LargeurFourier);
	TabHomogSection[cptSection] = ObjetCIma.HyperDistanceModuleFourier(FCentreGRE,FCentreGIM,FCentreDRE,FCentreDIM,LargeurFourier);
	
	
	
	TabBorSim[cptBord]=ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,LargeurFourier);
	
	cptTabG++;
	cptTabD++;
	cptSection++;
	cptBord++;
	}
}









//-----------------------------------------------------------------------------------------------------
// recherche des max pour l'inversion des fonctions de probalilité "section homogène" et bord similaire
//-----------------------------------------------------------------------------------------------------
for(i=0;i<TailleTabDist;i++)
{
	
	
	if(TabHomogSection[i]>MaxHomoSec)
	{
		MaxHomoSec = TabHomogSection[i];
	}

	if(TabBorSim[i]>MaxBordSim)
	{
		MaxBordSim = TabBorSim[i];
	}

	if(TabDistanceGauche[i]>MaxGauche)
	{
		MaxGauche = TabDistanceGauche[i];
	}

	if(TabDistanceDroite[i]>MaxDroite)
	{
		MaxDroite = TabDistanceDroite[i];
	}

}




//--------------------------------------------------------
//			Fonction export + calcul des couts
//--------------------------------------------------------

i=0; // compteur
j=0;



if (FichierStdIO.Open(NFichierTotale,CFile::modeCreate|CFile::modeWrite))
{
if (FichierVisu.Open(NFichierVisu,CFile::modeCreate|CFile::modeWrite)) 
{
	j=0;
	for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
	{
		
		//MaxIter = 0.0f;
		PosMax = IptInit;
		for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) // modif 21/11/2006 rajout pour limiter les surbalayage
		{
		
			//MuLTabDis[i] = ((MaxHomoSec-TabHomogSection[i])/MaxHomoSec) * (TabDistanceGauche[i]/MaxGauche) * (TabDistanceDroite[i]/MaxDroite)*(MaxBordSim-TabBorSim[i])/MaxBordSim;
		
			
			// toutes les contraintes


			/*ValDist = ((MaxHomoSec-TabHomogSection[i])/MaxHomoSec) * 
				(TabDistanceGauche[i]/MaxGauche) * 
				(TabDistanceDroite[i]/MaxDroite)*
				(MaxBordSim-TabBorSim[i])/MaxBordSim;
			ChaineDist.Format("%f\t",MuLTabDis[i]);*/


			// Sans la section homogène

			ValDist =(TabDistanceGauche[i]/MaxGauche) * 
				(TabDistanceDroite[i]/MaxDroite)*
				(MaxBordSim-TabBorSim[i])/MaxBordSim;

			ChaineDist.Format("%f\t",MuLTabDis[i]);

			FichierStdIO.WriteString(ChaineDist);
			
			ChaineDist.Format("%d\t",(IptInit));
			FichierStdIO.WriteString(ChaineDist);

			ChaineDist.Format("%d\t",(IptFinal));
			FichierStdIO.WriteString(ChaineDist);

			ChaineDist.Format("%d\t",(IptFinal-IptInit));
			FichierStdIO.WriteString(ChaineDist);
			FichierStdIO.WriteString("\n");

			
			//modif le 09/01/2007

			j=(IptInit+IptFinal)/2;

			if(TabStoreTemp[j]<ValDist)
			{
				/*MaxIter = MuLTabDis[i];
				PosMax = IptFinal;*/
				TabStoreTemp[j] = ValDist;
				TabStoreLarg[j] = IptFinal-IptInit;
			}
		
			i++;
		}//fin IptFinal
		/*j= (PosMax+IptInit)/2;
		if(TabStoreTemp[j]<MaxIter)TabStoreTemp[j]=MaxIter;// a sauver*/
		
		//j++;
     //rechercdu max pour la ligne stockage dans le tableau pour la visu
	}//fin IptInit
	
	for(j=0;j<(LongueurSignal);j++)
	{
		ChaineVisu.Format("%f\t",TabStoreTemp[j]);
		FichierVisu.WriteString(ChaineVisu);
		ChaineVisu.Format("%f\n",TabStoreLarg[j]);
		FichierVisu.WriteString(ChaineVisu);
	}

	FichierVisu.Close();
} // fin si fichier visu
FichierStdIO.Close();
} // fin if ouverture premier fichier






free(FGaucheRE); 
free(FDroiteRE);//---> debordement dans les boucles
free(FCentreGRE); 
free(FCentreDRE); 


free(FGaucheIM);
free(FDroiteIM);
free(FCentreGIM);
free(FCentreDIM);

free(FSDFixeRE);
free(FSDMobileRE);
free(FSDFixeIM);
free(FSDMobileIM);

free(TabDistanceGauche);
free(TabDistanceDroite);

free(TabHomogSection);
free(TabBorSim);
	

free(SignalPorte);
free(MuLTabDis);
free(TabMaxVisu);
free(TabStoreTemp);
free(TabLargeur);

}


// Initialisation du tableau contenant les valeur de la correction de Haming
void CCurviDetect1D::InitFenetreHaming(int LargeurFiltre)
{
	int k;
	double x,deuxpi=6.28318530717959;
	double pi=deuxpi/2.0;
	double pisur2=deuxpi/4.0;
	//int w = 25*LargeurFiltre;// -1
	//for (k=0;k<w;k++)
	//{
	//	x=(deuxpi*k)/(w);
	//	Fenetre[k/25]=(0.5-0.5*cos(x))/*/25.0*/; // bug += 0104
	//} // k	
	int w=LargeurFiltre-1;
	for (k=0;k<LargeurFiltre;k++)
	{
		x=((pi*(double)k)/(double)w)+pisur2;
		Fenetre[k]=(0.5-0.5*cos(x)); 
	} // k	



}

void CCurviDetect1D::CorrectionHaming(double* Signal, double* FenetreHaming, long Largeur)
{
	long i=0;

	for(i=0;i<Largeur;i++)
	{
		Signal[i]=Signal[i]*FenetreHaming[i];
	}
}

void CCurviDetect1D::RespCurviMaxi(BYTE *  TabCoupe)
{
int i=0,j=0;
double ValDistMax=0.0f;
int PosValDis = 0;

double TabAllMaxi[4];
double MaxNorma=0.0f;

//------------------------
//			RAZ 
//------------------------
for(i=0;i<LargeurFourier;i++)
{

	FGaucheRE[i]=0.0f;
	FDroiteRE[i]=0.0f;
	FCentreGRE[i]=0.0f;
	FCentreDRE[i]=0.0f;

	FGaucheIM[i]=0.0f;
	FDroiteIM[i]=0.0f;
	FCentreGIM[i]=0.0f;
	FCentreDIM[i]=0.0f;

	FSDFixeRE[i]=0.0f;
	FSDMobileRE[i]=0.0f;
	FSDFixeIM[i]=0.0f;
	FSDMobileIM[i]=0.0f;

}

for (i=0;i<(LongueurSignal);i++)
{
	TabStoreTemp[i].ValDist = 0.0f;
	TabStoreTemp[i].Centre = 0;
	TabStoreTemp[i].IptFinal = 0;
	TabStoreTemp[i].IptInit = 0;
	TabStoreTemp[i].Largeur = 0;
}

for (i=0;i<(IptInitMax-IptInitMin);i++)
{
	TabMaxVisu[i] = 0.0f;
	TabLargeur[i] = 0;
}

for(i=0;i<LargeurFourier;i++)
{

	TabDistanceGauche[i] = 0.0f;
	TabDistanceDroite[i] = 0.0f;
	

}

for(i=0;i<TailleTabDist;i++)
{
	MuLTabDis[i] = 0.0f;
	TabDistanceGauche[i]=0.0f;
	TabDistanceDroite[i]=0.0f;
	TabHomogSection[i]=0.0f;
	TabBorSim[i]=0.0f;


}

cptTabG=0;
cptTabD=0;
cptSection=0;
cptBord=0;


//----------------------------------------------------
// recopie des tab pour fourier + calcul distance
//----------------------------------------------------
for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
{
	for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) //idem modif 21/11/2006
	{
	//copie Fgauche 
	
	cpt =0;
	for (i = IptInit-LargeurFourier;i<IptInit;i++)
	{
	FGaucheRE[cpt]= (double)TabCoupe[i];
	cpt++;
	}
	
	cpt =0;
	for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)
	{
	FDroiteRE[cpt]= (double)TabCoupe[i];
	cpt++;
	}

	cpt =0;
	for (i = IptInit+1;i<IptInit+LargeurFourier+1;i++)//modif 09/01/2007
	{
	FCentreGRE[cpt]= (double)TabCoupe[i];
	cpt++;
	}
	
	cpt =0;
	for (i = IptFinal-LargeurFourier;i<IptFinal;i++)//modif 09/01/2007
	{
	FCentreDRE[cpt]= (double)TabCoupe[i];
	cpt++;
	}

	cpt=0;

	//Correction de Haming
	CorrectionHaming(FGaucheRE,Fenetre,LargeurFourier);
	CorrectionHaming(FGaucheIM,Fenetre,LargeurFourier);
	CorrectionHaming(FDroiteRE,Fenetre,LargeurFourier);
	CorrectionHaming(FDroiteIM,Fenetre,LargeurFourier);
	CorrectionHaming(FCentreDRE,Fenetre,LargeurFourier);
	CorrectionHaming(FCentreDIM,Fenetre,LargeurFourier);
	CorrectionHaming(FCentreGRE,Fenetre,LargeurFourier);
	CorrectionHaming(FCentreGIM,Fenetre,LargeurFourier);
	//calcul FFT
	ObjetCIma.FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);
	ObjetCIma.FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);
	ObjetCIma.FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);
	ObjetCIma.FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);

					
	TabDistanceGauche[cptTabG] = ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,LargeurFourier);
	TabDistanceDroite[cptTabD] = ObjetCIma.HyperDistanceModuleFourier(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,LargeurFourier);
	TabHomogSection[cptSection] = ObjetCIma.HyperDistanceModuleFourier(FCentreGRE,FCentreGIM,FCentreDRE,FCentreDIM,LargeurFourier);
	
	
	
	TabBorSim[cptBord]=ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,LargeurFourier);
	
	cptTabG++;
	cptTabD++;
	cptSection++;
	cptBord++;
	}
}


//-----------------------------------------------------------------------------------------------------
// recherche des max pour l'inversion des fonctions de probalilité "section homogène" et bord similaire
//-----------------------------------------------------------------------------------------------------
for(i=0;i<TailleTabDist;i++)
{
	
	
	if(TabHomogSection[i]>MaxHomoSec)
	{
		MaxHomoSec = TabHomogSection[i];
		
	}

	if(TabBorSim[i]>MaxBordSim)
	{
		MaxBordSim = TabBorSim[i];
		
	}

	if(TabDistanceGauche[i]>MaxGauche)
	{
		MaxGauche = TabDistanceGauche[i];
		
	}

	if(TabDistanceDroite[i]>MaxDroite)
	{
		MaxDroite = TabDistanceDroite[i];
		
	}

}
TabAllMaxi[0] = MaxHomoSec;
TabAllMaxi[1] = MaxBordSim;
TabAllMaxi[2] = MaxGauche;
TabAllMaxi[3] = MaxDroite;

for(i=0;i<4;i++)
{
	if (TabAllMaxi[i]>MaxNorma)
	{
		MaxNorma = TabAllMaxi[i];
	}
}

if (MaxNorma == 0.0f)
{
	MaxNorma=1;
}
//--------------------------------------------------------
//			Fonction export + calcul des couts
//--------------------------------------------------------

i=0; // compteur
j=0;



	j=0;
	for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
	{
		
		//MaxIter = 0.0f;
//		PosMax = IptInit; ---> modif cedric 18 juin 2007
		for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) // modif 21/11/2006 rajout pour limiter les surbalayage
		{
			// Sans la section homogène

			/*ValDist =(TabDistanceGauche[i]/MaxNorma) * 
				(TabDistanceDroite[i]/MaxNorma);*/

			/*ValDist =(TabDistanceGauche[i]/MaxNorma) * 
				(TabDistanceDroite[i]/MaxNorma)*
				((MaxNorma-TabBorSim[i])/MaxNorma);*/

			ValDist =(TabDistanceGauche[i]/MaxNorma) * 
				(TabDistanceDroite[i]/MaxNorma)*
				((MaxNorma-TabBorSim[i])/MaxNorma)*
				/*(MaxNorma-TabBorSim[i])/MaxNorma **/
				(MaxHomoSec-TabHomogSection[i])/MaxHomoSec;


			
			//modif le 09/01/2007

			j=(IptInit+IptFinal)/2;

			//modif temporaire (18/10/2007) pour tracer les PFD
			if(IptFinal-IptInit <= LargeurFourier+(LargeurFourier/2)+2)ValDist=0.0f;

			if(TabStoreTemp[j].ValDist<ValDist)
			{
				/*MaxIter = MuLTabDis[i];
				PosMax = IptFinal;*/
				TabStoreTemp[j].ValDist = ValDist;
				TabStoreTemp[j].Largeur = IptFinal-IptInit;
				TabStoreTemp[j].IptInit = IptInit;
				TabStoreTemp[j].IptFinal = IptFinal;
				TabStoreTemp[j].Centre = (IptFinal+IptInit)/2;

				//TabStoreLarg[j] = IptFinal-IptInit;
			}
		
			i++;
		}//fin IptFinal
//		if(TabStoreTemp[j]<MaxIter)TabStoreTemp[j]=MaxIter;// a sauver*/
		
		//j++;
	}//fin IptInit
	
	for(j=0;j<(LongueurSignal);j++)
	{
		if(TabStoreTemp[j].ValDist>ValDistMax)
		{
			ValDistMax = TabStoreTemp[j].ValDist;
			PosValDis = j;

		}
		
	}

	MaxiCurvi = TabStoreTemp[PosValDis];


}

CCurviDetect1D::CCurviDetect1D(int ValLongueurSignal, int ValLargeurFourier)
{
int i=0;
LongueurSignal = ValLongueurSignal;
LargeurFourier = ValLargeurFourier;
DemiLargeurFourier = LargeurFourier/2;
DemiLongueurSignal = LongueurSignal/2;



IptInitMax = LongueurSignal-LargeurFourier;//-LargeurPortemodif cédric 01/02/2007
IptFinalMax = LongueurSignal-LargeurFourier-1; //attention 
IptInitMin = LargeurFourier;
IptFinalMin = IptInitMin+LargeurFourier ;//+ LargeurPorte; //modif 01/02/2007

RacineLargeurFourier = 0 ;//(int)sqrt((double)LargeurFourier); bug 0204
while((1<<RacineLargeurFourier)<LargeurFourier) RacineLargeurFourier++;
FGaucheRE = (double *) malloc(sizeof(double)*LargeurFourier);
FDroiteRE = (double *) malloc(sizeof(double)*LargeurFourier);
FCentreGRE = (double *) malloc(sizeof(double)*LargeurFourier);
FCentreDRE = (double *) malloc(sizeof(double)*LargeurFourier);


FGaucheIM = (double *) malloc(sizeof(double)*LargeurFourier);
FDroiteIM = (double *) malloc(sizeof(double)*LargeurFourier);
FCentreGIM = (double *) malloc(sizeof(double)*LargeurFourier);
FCentreDIM = (double *) malloc(sizeof(double)*LargeurFourier);

FSDFixeRE = (double *) malloc(sizeof(double)*LargeurFourier);
FSDMobileRE = (double *) malloc(sizeof(double)*LargeurFourier);

FSDFixeIM = (double *) malloc(sizeof(double)*LargeurFourier);
FSDMobileIM = (double *) malloc(sizeof(double)*LargeurFourier);

MaxGauche = 0.0f,MaxDroite = 0.0f,MaxHomoSec = 0.0f, MaxBordSim = 0.0f;

TailleTabDist = (IptInitMax-IptInitMin)*(IptFinalMax-IptFinalMin);

TabDistanceGauche = (double *) malloc(sizeof(double)*TailleTabDist);
TabDistanceDroite = (double *) malloc(sizeof(double)*TailleTabDist);
TabHomogSection = (double *) malloc(sizeof(double)*TailleTabDist);
TabBorSim = (double *) malloc(sizeof(double)*TailleTabDist);


MuLTabDis = (double *) malloc(sizeof(double)*TailleTabDist);

TabMaxVisu = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

TabStoreTemp = (RespCurvi *)malloc(sizeof(RespCurvi) *(LongueurSignal));
TabStoreLarg = (RespCurvi *)malloc(sizeof(RespCurvi) *(LongueurSignal));
TabLargeur = (int *)malloc(sizeof(int) *(IptInitMax-IptInitMin));

// Allocation du tableau qui contiendra la fenetre de hamming
Fenetre = (double*)malloc(sizeof(double)*LargeurFourier);
// Initialisation de la Fenetre de Haming
InitFenetreHaming(LargeurFourier);

for(i=0;i<LargeurFourier;i++)
{

	FGaucheRE[i]=0.0f;
	FDroiteRE[i]=0.0f;
	FCentreGRE[i]=0.0f;
	FCentreDRE[i]=0.0f;

	FGaucheIM[i]=0.0f;
	FDroiteIM[i]=0.0f;
	FCentreGIM[i]=0.0f;
	FCentreDIM[i]=0.0f;

	FSDFixeRE[i]=0.0f;
	FSDMobileRE[i]=0.0f;
	FSDFixeIM[i]=0.0f;
	FSDMobileIM[i]=0.0f;

}

for (i=0;i<(LongueurSignal);i++)
{
	TabStoreTemp[i].ValDist = 0.0f;
	TabStoreTemp[i].Centre = 0;
	TabStoreTemp[i].IptFinal = 0;
	TabStoreTemp[i].IptInit = 0;
	TabStoreTemp[i].Largeur = 0;
}

for (i=0;i<(IptInitMax-IptInitMin);i++)
{
	TabMaxVisu[i] = 0.0f;
	TabLargeur[i] = 0;
}

for(i=0;i<LargeurFourier;i++)
{

	TabDistanceGauche[i] = 0.0f;
	TabDistanceDroite[i] = 0.0f;
	

}

for(i=0;i<TailleTabDist;i++)
{
	MuLTabDis[i] = 0.0f;
	TabDistanceGauche[i]=0.0f;
	TabDistanceDroite[i]=0.0f;
	TabHomogSection[i]=0.0f;
	TabBorSim[i]=0.0f;


}

}

CCurviDetect1D::CCurviDetect1D(int ValLongueurSignal, int ValLargeurFourier, int ValIptInit)
{
int i=0;

LongueurSignal = ValLongueurSignal;
LargeurFourier = ValLargeurFourier;
DemiLargeurFourier = LargeurFourier/2;
DemiLongueurSignal = LongueurSignal/2;


IptInit = ValIptInit;

//IptInitMax = LongueurSignal-LargeurFourier;//-LargeurPortemodif cédric 01/02/2007
IptFinalMax = LongueurSignal-LargeurFourier-1; //attention 
IptInitMin = LargeurFourier;
//IptFinalMin = IptInitMin+LargeurFourier ;//+ LargeurPorte; //modif 01/02/2007

RacineLargeurFourier = 0;//(int)sqrt((double)LargeurFourier);

while((1<<RacineLargeurFourier)<LargeurFourier) RacineLargeurFourier++;

FGaucheRE = (double *) malloc(sizeof(double)*LargeurFourier);
FDroiteRE = (double *) malloc(sizeof(double)*LargeurFourier);
FCentreGRE = (double *) malloc(sizeof(double)*LargeurFourier);
FCentreDRE = (double *) malloc(sizeof(double)*LargeurFourier);


FGaucheIM = (double *) malloc(sizeof(double)*LargeurFourier);
FDroiteIM = (double *) malloc(sizeof(double)*LargeurFourier);
FCentreGIM = (double *) malloc(sizeof(double)*LargeurFourier);
FCentreDIM = (double *) malloc(sizeof(double)*LargeurFourier);

FSDFixeRE = (double *) malloc(sizeof(double)*LargeurFourier);
FSDMobileRE = (double *) malloc(sizeof(double)*LargeurFourier);

FSDFixeIM = (double *) malloc(sizeof(double)*LargeurFourier);
FSDMobileIM = (double *) malloc(sizeof(double)*LargeurFourier);

MaxGauche = 0.0f,MaxDroite = 0.0f,MaxHomoSec = 0.0f, MaxBordSim = 0.0f;

//TailleTabDist = (IptInitMax-IptInitMin)*(IptFinalMax-IptFinalMin);
TailleTabDist = (IptFinalMax-IptInit); // 13_11_2007 pas sur

TabDistanceGauche = (double *) malloc(sizeof(double)*TailleTabDist);
TabDistanceDroite = (double *) malloc(sizeof(double)*TailleTabDist);
TabHomogSection = (double *) malloc(sizeof(double)*TailleTabDist);
TabBorSim = (double *) malloc(sizeof(double)*TailleTabDist);


MuLTabDis = (double *) malloc(sizeof(double)*TailleTabDist);

//TabMaxVisu = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin)); // 13_11_2007

//TabStoreTemp = (RespCurvi *)malloc(sizeof(RespCurvi) *(LongueurSignal));
TabStoreLarg = (RespCurvi *)malloc(sizeof(RespCurvi) *(LongueurSignal));
// TabLargeur = (int *)malloc(sizeof(int) *(IptInitMax-IptInitMin)); // 13_11_2007

// Allocation du tableau qui contiendra la fenetre de hamming
Fenetre = (double*)malloc(sizeof(double)*LargeurFourier);
// Initialisation de la Fenetre de Haming
InitFenetreHaming(LargeurFourier);

for(i=0;i<LargeurFourier;i++)
{

	FGaucheRE[i]=0.0f;
	FDroiteRE[i]=0.0f;
	FCentreGRE[i]=0.0f;
	FCentreDRE[i]=0.0f;

	FGaucheIM[i]=0.0f;
	FDroiteIM[i]=0.0f;
	FCentreGIM[i]=0.0f;
	FCentreDIM[i]=0.0f;

	FSDFixeRE[i]=0.0f;
	FSDMobileRE[i]=0.0f;
	FSDFixeIM[i]=0.0f;
	FSDMobileIM[i]=0.0f;

}

//for (i=0;i<(LongueurSignal);i++)
//{
//	TabStoreTemp[i].ValDist = 0.0f;
//	TabStoreTemp[i].Centre = 0;
//	TabStoreTemp[i].IptFinal = 0;
//	TabStoreTemp[i].IptInit = 0;
//	TabStoreTemp[i].Largeur = 0;
//}

//for (i=0;i<(IptInitMax-IptInitMin);i++)
//{
//	TabMaxVisu[i] = 0.0f;
//	TabLargeur[i] = 0;
//} -->// 13_11_2007

for(i=0;i<LargeurFourier;i++)
{

	TabDistanceGauche[i] = 0.0f;
	TabDistanceDroite[i] = 0.0f;
	

}

for(i=0;i<TailleTabDist;i++)
{
	MuLTabDis[i] = 0.0f;
	TabDistanceGauche[i]=0.0f;
	TabDistanceDroite[i]=0.0f;
	TabHomogSection[i]=0.0f;
	TabBorSim[i]=0.0f;


}
}

void CCurviDetect1D::RespCurviMaxiFast(BYTE * TabCoupe, CPoint *TabCoupeXY, double WidthPrior, int IptInitFixe,int *ImageEdge,int ImageWidth)
{

IptInit = IptInitFixe;
int i=0,j=0;
double ValDistMax=0.0f;
int PosValDis = 0;

double TabAllMaxi[4];
double MaxNorma=0.0f;

long TailleTabDistance2Passe;

double * TabDistance2Passe;
double * MuLTabDis2Passe;
double * TabDistanceGauche2Passe;
double * TabDistanceDroite2Passe;
double * TabHomogSection2Passe;
double * TabBorSim2Passe;
double TailleGauss=91.0;
TempValDis2Passe=0.0f;

//------------------------
//			RAZ 
//------------------------
for(i=0;i<LargeurFourier;i++)
{

	FGaucheRE[i]=0.0f;
	FDroiteRE[i]=0.0f;
	FCentreGRE[i]=0.0f;
	FCentreDRE[i]=0.0f;

	FGaucheIM[i]=0.0f;
	FDroiteIM[i]=0.0f;
	FCentreGIM[i]=0.0f;
	FCentreDIM[i]=0.0f;

	FSDFixeRE[i]=0.0f;
	FSDMobileRE[i]=0.0f;
	FSDFixeIM[i]=0.0f;
	FSDMobileIM[i]=0.0f;

}

SecondePasse.ValDist = 0.0f;
SecondePasse.IptInit = 0;
SecondePasse.IptFinal = 0;
SecondePasse.Largeur = 0;
SecondePasse.Centre=0;


MaxiCurvi.ValDist = 0.0f;
MaxiCurvi.IptInit=0;
MaxiCurvi.Centre=0;
MaxiCurvi.IptFinal=0;
MaxiCurvi.Largeur=0;

//for (i=0;i<(IptInitMax-IptInitMin);i++)
//{
//	TabMaxVisu[i] = 0.0f;
//	TabLargeur[i] = 0;
//} // 13_11_2007

for(i=0;i<LargeurFourier;i++)
{

	TabDistanceGauche[i] = 0.0f;
	TabDistanceDroite[i] = 0.0f;
	

}

for(i=0;i<TailleTabDist;i++)
{
	MuLTabDis[i] = 0.0f;
	TabDistanceGauche[i]=0.0f;
	TabDistanceDroite[i]=0.0f;
	TabHomogSection[i]=0.0f;
	TabBorSim[i]=0.0f;


}

cptTabG=0;
cptTabD=0;
cptSection=0;
cptBord=0;

	
double DistanceGaucheFixe=1.0;
	
#ifdef DISTGAUCHE

cpt =0;
for (i = IptInit-LargeurFourier;i<IptInit;i++)
{
FGaucheRE[cpt]= (double)TabCoupe[i]; FGaucheIM[cpt]= 0.0f; cpt++; // bug RE 0104
}

cpt =0;
for (i = IptInit+1;i<IptInit+LargeurFourier+1;i++)//modif 09/01/2007
{
FCentreGRE[cpt]= (double)TabCoupe[i]; FCentreGIM[cpt]= 0.0f; cpt++;
}

CorrectionHaming(FGaucheRE,Fenetre,LargeurFourier);
CorrectionHaming(FGaucheIM,Fenetre,LargeurFourier);

CorrectionHaming(FCentreGRE,Fenetre,LargeurFourier);
CorrectionHaming(FCentreGIM,Fenetre,LargeurFourier);

ObjetCIma.FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);
ObjetCIma.FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);	

DistanceGaucheFixe = ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,LargeurFourier);	
#endif
//----------------------------------------------------
// recopie des tab pour fourier + calcul distance
//----------------------------------------------------

for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) //idem modif 21/11/2006
{
	int Edge=ImageEdge[TabCoupeXY[IptFinal].x+TabCoupeXY[IptFinal].y*ImageWidth];
	//Edge=1;
	if (Edge!=0)
	{
		cpt =0;
		for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)
		{
			FDroiteRE[cpt]= (double)TabCoupe[i]; FDroiteIM[cpt]= 0.0f; cpt++;
		}
		cpt=0;
		for (i = IptFinal-LargeurFourier;i<IptFinal;i++)//modif 09/01/2007
		{
			FCentreDRE[cpt]= (double)TabCoupe[i]; FCentreDIM[cpt]= 0.0f; cpt++;
		}
		cpt =0;


		//Correction de Haming
		CorrectionHaming(FDroiteRE,Fenetre,LargeurFourier);
		CorrectionHaming(FDroiteIM,Fenetre,LargeurFourier);
		CorrectionHaming(FCentreDRE,Fenetre,LargeurFourier);
		CorrectionHaming(FCentreDIM,Fenetre,LargeurFourier);
		
		//calcul FFT
		ObjetCIma.FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);
		ObjetCIma.FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);
					

		//TabHomogSection[cptSection] = ObjetCIma.HyperDistanceModuleFourier(FCentreGRE,FCentreGIM,FCentreDRE,FCentreDIM,LargeurFourier);
		TabDistanceDroite[cptTabD]=ObjetCIma.HyperDistanceModuleFourier(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,LargeurFourier);
		TabBorSim[cptBord]=ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,LargeurFourier);
	}
	else
	{
		TabDistanceDroite[cptTabD]=0.0;
		TabBorSim[cptBord]=0.0;
	}

	cptTabG++;
	cptTabD++;
	cptSection++;
	cptBord++;

} // fin for Iptfinal
//CStdioFile FichierVisu;
//CString ChaineVisu;
//FichierVisu.Open("c:\\cerclesbruite8.txt",CFile::modeCreate|CFile::modeWrite);
//
//for(i=0;i<TailleTabDist;i++)
//{
//	ChaineVisu.Format("%f\t",DistanceGaucheFixe);
//	FichierVisu.WriteString(ChaineVisu);
//	ChaineVisu.Format("%f\t",TabBorSim[i]);
//	FichierVisu.WriteString(ChaineVisu);
//	ChaineVisu.Format("%f\n",TabDistanceDroite[i]);
//	FichierVisu.WriteString(ChaineVisu);
//}
//FichierVisu.Close();

//-----------------------------------------------------------------------------------------------------
// recherche des max pour l'inversion des fonctions de probalilité "section homogène" et bord similaire
//-----------------------------------------------------------------------------------------------------
for(i=0;i<TailleTabDist;i++)
{	
	//if(TabHomogSection[i]>MaxHomoSec)MaxHomoSec = TabHomogSection[i];	
	if(TabBorSim[i]>MaxBordSim)MaxBordSim = TabBorSim[i];	
	//if(TabDistanceGauche[i]>MaxGauche)MaxGauche = TabDistanceGauche[i];	
	if(TabDistanceDroite[i]>MaxDroite)MaxDroite = TabDistanceDroite[i];
}
TabAllMaxi[0] = DistanceGaucheFixe;
TabAllMaxi[1] = MaxDroite;
TabAllMaxi[2] = MaxBordSim;
TabAllMaxi[3] = 0.0f;

for(i=0;i<4;i++)
	if (TabAllMaxi[i]>MaxNorma)	MaxNorma = TabAllMaxi[i];

if (MaxNorma == 0.0f)MaxNorma=1;
#ifndef DISTGAUCHE
DistanceGaucheFixe=MaxNorma;
#endif
//--------------------------------------------------------
//			Fonction export + calcul des couts
//--------------------------------------------------------

i=0; // compteur
	for(IptFinal=(IptInit+(LargeurFourier)+1),i=0;IptFinal<IptFinalMax;IptFinal++,i++) // modif 21/11/2006 rajout pour limiter les surbalayage
	{
		double width = sqrt(
		((double)TabCoupeXY[IptInitFixe].x - (double)TabCoupeXY[IptFinal].x)*
		((double)TabCoupeXY[IptInitFixe].x - (double)TabCoupeXY[IptFinal].x)
		+
		((double)TabCoupeXY[IptInitFixe].y - (double)TabCoupeXY[IptFinal].y)*
		((double)TabCoupeXY[IptInitFixe].y - (double)TabCoupeXY[IptFinal].y));
		int g;
		if (WidthPrior>0)
		{
		 g = int( ( (width-WidthPrior)/TailleGauss)*10.0 )+45;  // 45
		 g = __max(0, __min(g, 90));
		} 
		else 
		g = 45;

		ValDist =(TabDistanceDroite[i]/MaxNorma) * 
			(DistanceGaucheFixe/MaxNorma)*
			((MaxNorma-TabBorSim[i])/MaxNorma)
				*gauss_table_sigma_1[g];///gauss_table_sigma_1[45];

		int Edge=ImageEdge[TabCoupeXY[IptFinal].x+TabCoupeXY[IptFinal].y*ImageWidth];
		//Edge=1;
		if(MaxiCurvi.ValDist<ValDist && Edge!=0)
		{
			MaxiCurvi.ValDist = ValDist;
			MaxiCurvi.Largeur = width;
			MaxiCurvi.IptInit = IptInit;
			MaxiCurvi.IptFinal = IptFinal;
			MaxiCurvi.Centre = (IptFinal+IptInit)/2;
		}

	}//fin IptFinal

//-------------------------------------------------------------
//2° passe
//-------------------------------------------------------------

SecondePasse=MaxiCurvi ; // init

if(MaxiCurvi.IptFinal-IptInitFixe > LargeurFourier+1)
{
	MaxGauche2Passe=0.0f;
	MaxDroite2Passe=0.0f;
	MaxHomoSec2Passe=0.0f;
	MaxBordSim2Passe=0.0f;
	ValDist=0.0f;

	TailleTabDistance2Passe = MaxiCurvi.IptFinal-LargeurFourier -IptInitFixe -1;
	TabDistance2Passe = (double *)malloc(sizeof(double)*TailleTabDistance2Passe);
	MuLTabDis2Passe  = (double *)malloc(sizeof(double)*TailleTabDistance2Passe);
	TabDistanceGauche2Passe  = (double *)malloc(sizeof(double)*TailleTabDistance2Passe);
	TabDistanceDroite2Passe = (double *)malloc(sizeof(double)*TailleTabDistance2Passe);
	TabHomogSection2Passe = (double *)malloc(sizeof(double)*TailleTabDistance2Passe);
	TabBorSim2Passe = (double *)malloc(sizeof(double)*TailleTabDistance2Passe);

	//------------------------
	//			RAZ 
	//------------------------
	for(i=0;i<LargeurFourier;i++)
	{
		FGaucheRE[i]=0.0f;
		FDroiteRE[i]=0.0f;
		FCentreGRE[i]=0.0f;
		FCentreDRE[i]=0.0f;
		FGaucheIM[i]=0.0f;
		FDroiteIM[i]=0.0f;
		FCentreGIM[i]=0.0f;
		FCentreDIM[i]=0.0f;
		FSDFixeRE[i]=0.0f;
		FSDMobileRE[i]=0.0f;
		FSDFixeIM[i]=0.0f;
		FSDMobileIM[i]=0.0f;
	}

	for(i=0;i<LargeurFourier;i++)
	{
		TabDistanceGauche[i] = 0.0f; TabDistanceDroite[i] = 0.0f;
	}

	for(i=0;i<TailleTabDistance2Passe;i++)
	{
		TabDistance2Passe[i]=0.0f;
		MuLTabDis2Passe[i]=0.0f;
		TabDistanceGauche2Passe[i]=0.0f;
		TabDistanceDroite2Passe[i]=0.0f;
		TabHomogSection2Passe[i]=0.0f;
		TabBorSim2Passe[i]=0.0f;
	}

	cptTabG=0;
	cptTabD=0;
	cptSection=0;
	cptBord=0;


	cpt=0;
	for (i = MaxiCurvi.IptFinal-LargeurFourier;i<MaxiCurvi.IptFinal;i++)//modif 09/01/2007
	{
		FCentreDRE[cpt]= (double)TabCoupe[i]; FCentreDIM[cpt] = 0.0f; cpt++;
	}
	
	cpt =0;
	for (i = MaxiCurvi.IptFinal+1;i<MaxiCurvi.IptFinal+LargeurFourier+1;i++)
	{
		FDroiteRE[cpt]= (double)TabCoupe[i]; FDroiteIM[cpt]= 0.0f; cpt++;
	}

	CorrectionHaming(FDroiteRE,Fenetre,LargeurFourier);
	CorrectionHaming(FDroiteIM,Fenetre,LargeurFourier);
	CorrectionHaming(FCentreDRE,Fenetre,LargeurFourier);
	CorrectionHaming(FCentreDIM,Fenetre,LargeurFourier);
	
	ObjetCIma.FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);
	ObjetCIma.FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);

	double DistanceDroiteFixe = ObjetCIma.HyperDistanceModuleFourier(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,LargeurFourier);
//----------------------------------------------------
// recopie des tab pour fourier + calcul distance
//----------------------------------------------------
	for(IptFinal=IptInitFixe;IptFinal<MaxiCurvi.IptFinal-LargeurFourier -1;IptFinal++) //idem modif 21/11/2006
	{
		cpt =0;
		for (i = IptFinal-LargeurFourier;i<IptFinal;i++)
		{
			FGaucheRE[cpt]= (double)TabCoupe[i]; FGaucheIM[cpt]= 0.0f; cpt++;
		}

		cpt =0;
		for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)//modif 09/01/2007
		{
			FCentreGRE[cpt]= (double)TabCoupe[i]; FCentreGIM[cpt]= 0.0f; cpt++;
		}
		
		cpt =0;
		
		//Correction de Haming
		CorrectionHaming(FGaucheRE,Fenetre,LargeurFourier);
		CorrectionHaming(FGaucheIM,Fenetre,LargeurFourier);
		CorrectionHaming(FCentreGRE,Fenetre,LargeurFourier);
		CorrectionHaming(FCentreGIM,Fenetre,LargeurFourier);

		//calcul FFT
		ObjetCIma.FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);
		ObjetCIma.FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);

		TabDistanceGauche2Passe[cptTabG] = ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,LargeurFourier);
		//TabHomogSection2Passe[cptSection] = ObjetCIma.HyperDistanceModuleFourier(FCentreGRE,FCentreGIM,FCentreDRE,FCentreDIM,LargeurFourier);
		TabBorSim2Passe[cptBord]=ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,LargeurFourier);
	
		cptTabG++;
		cptTabD++;
		cptSection++;
		cptBord++;
	}



//-----------------------------------------------------------------------------------------------------
// recherche des max pour l'inversion des fonctions de probalilité "section homogène" et bord similaire
//-----------------------------------------------------------------------------------------------------



	for(i=0;i<TailleTabDistance2Passe;i++)
	{
		//if(TabHomogSection[i]>MaxHomoSec2Passe)MaxHomoSec2Passe = TabHomogSection2Passe[i];
		if(TabBorSim2Passe[i]>MaxBordSim2Passe)MaxBordSim2Passe = TabBorSim2Passe[i];
		if(TabDistanceGauche2Passe[i]>MaxGauche2Passe)MaxGauche2Passe = TabDistanceGauche2Passe[i];
		//if(TabDistanceDroite2Passe[i]>MaxDroite2Passe)MaxDroite2Passe = TabDistanceDroite2Passe[i];
	}

	TabAllMaxi[0]=DistanceDroiteFixe;
	TabAllMaxi[1]=MaxBordSim2Passe;
	TabAllMaxi[2]=MaxGauche2Passe;
	TabAllMaxi[3]=0.0f;

	for(i=0;i<4;i++)
	{
		if (TabAllMaxi[i]>MaxNorma)
		{
			MaxNorma = TabAllMaxi[i];
		}
	}

	if (MaxNorma == 0.0f)
	{
		MaxNorma=1.0;
	}
//--------------------------------------------------------
//			Fonction export + calcul des couts
//--------------------------------------------------------
	SecondePasse.ValDist=0.0f; // bug 0204
	int Borne=__min(MaxiCurvi.IptInit+LargeurFourier,MaxiCurvi.IptFinal-LargeurFourier -1);
	for(IptFinal=IptInitFixe,i=0;IptFinal<Borne/*MaxiCurvi.IptFinal-LargeurFourier -1*/;IptFinal++,i++) // modif 21/11/2006 rajout pour limiter les surbalayage
	{
		double width = sqrt(
			((double)TabCoupeXY[IptFinal].x - (double)TabCoupeXY[MaxiCurvi.IptFinal].x)*
			((double)TabCoupeXY[IptFinal].x - (double)TabCoupeXY[MaxiCurvi.IptFinal].x)
			+
			((double)TabCoupeXY[IptFinal].y - (double)TabCoupeXY[MaxiCurvi.IptFinal].y)*
			((double)TabCoupeXY[IptFinal].y - (double)TabCoupeXY[MaxiCurvi.IptFinal].y));
		int g;
		if (WidthPrior>0)
		{
		//			 g = int(((width/WidthPrior-1.0)/0.1))+45; 
		 g = int( ( (width-WidthPrior)/TailleGauss)*10.0 )+45;  
		 g = __max(0, __min(g, 90));
		} else 
		 g = 45;

		// Sans la section homogène
		
		TempValDis2Passe =(TabDistanceGauche2Passe[i]/MaxNorma) * 
			(DistanceDroiteFixe/MaxNorma)*
			((MaxNorma-TabBorSim2Passe[i])/MaxNorma)
			 *gauss_table_sigma_1[g];///gauss_table_sigma_1[45]; 

		if(SecondePasse.ValDist<TempValDis2Passe)
		{
			SecondePasse.ValDist = TempValDis2Passe;
			SecondePasse.IptInit = IptFinal;
			SecondePasse.IptFinal = MaxiCurvi.IptFinal;
			SecondePasse.Largeur = width;
			//SecondePasse.Centre = (SecondePasse.IptFinal+SecondePasse.IptInit)/2;
			//ValDist=TempValDis2Passe;
		}
	}//fin IptFinal


	SecondePasse.Centre = (SecondePasse.IptFinal+SecondePasse.IptInit)/2;

	//MaxiCurvi = SecondePasse; // test 0104

	free(TabDistance2Passe);
	free(MuLTabDis2Passe);
	free(TabDistanceGauche2Passe);
	free(TabDistanceDroite2Passe);
	free(TabHomogSection2Passe);
	free(TabBorSim2Passe);
	} // fin si

}



void CCurviDetect1D::RespCurviMaxiFastCoul(rgb * TabCoupeCoul, CPoint *TabCoupeXY, double WidthPrior, int IptInitFixe,int *ImageEdge,int ImageWidth)
{

IptInit = IptInitFixe;
int i=0,j=0;
double ValDistMax=0.0f;
int PosValDis = 0;

double TabAllMaxi[4];
double MaxNorma=0.0f;

long TailleTabDistance2Passe;

double * TabDistance2Passe;
double * MuLTabDis2Passe;
double * TabDistanceGauche2Passe;
double * TabDistanceDroite2Passe;
double * TabHomogSection2Passe;
double * TabBorSim2Passe;
double TailleGauss=91.0;
TempValDis2Passe=0.0f;

//------------------------
//			RAZ 
//------------------------
for(i=0;i<LargeurFourier;i++)
{

	FGaucheRE[i]=0.0f;
	FDroiteRE[i]=0.0f;
	FCentreGRE[i]=0.0f;
	FCentreDRE[i]=0.0f;

	FGaucheIM[i]=0.0f;
	FDroiteIM[i]=0.0f;
	FCentreGIM[i]=0.0f;
	FCentreDIM[i]=0.0f;

	FSDFixeRE[i]=0.0f;
	FSDMobileRE[i]=0.0f;
	FSDFixeIM[i]=0.0f;
	FSDMobileIM[i]=0.0f;

}

SecondePasse.ValDist = 0.0f;
SecondePasse.IptInit = 0;
SecondePasse.IptFinal = 0;
SecondePasse.Largeur = 0;
SecondePasse.Centre=0;


MaxiCurvi.ValDist = 0.0f;
MaxiCurvi.IptInit=0;
MaxiCurvi.Centre=0;
MaxiCurvi.IptFinal=0;
MaxiCurvi.Largeur=0;

//for (i=0;i<(IptInitMax-IptInitMin);i++)
//{
//	TabMaxVisu[i] = 0.0f;
//	TabLargeur[i] = 0;
//} // 13_11_2007

for(i=0;i<LargeurFourier;i++)
{
	TabDistanceGauche[i] = 0.0f;
	TabDistanceDroite[i] = 0.0f;
}

for(i=0;i<TailleTabDist;i++)
{
	MuLTabDis[i] = 0.0f;
	TabDistanceGauche[i]=0.0f;
	TabDistanceDroite[i]=0.0f;
	TabHomogSection[i]=0.0f;
	TabBorSim[i]=0.0f;
}

cptTabG=0;
cptTabD=0;
cptSection=0;
cptBord=0;

double DistanceGaucheFixe=1.0;

#ifdef DISTGAUCHE
cpt =0;
for (i = IptInit-LargeurFourier;i<IptInit;i++)
{
	FGaucheRE[cpt]= (double)TabCoupeCoul[i].r; FGaucheIM[cpt]= 0.0f; cpt++; // bug RE 0104
}
cpt =0;
for (i = IptInit+1;i<IptInit+LargeurFourier+1;i++)//modif 09/01/2007
{
	FCentreGRE[cpt]= (double)TabCoupeCoul[i].r; FCentreGIM[cpt]= 0.0f; cpt++;
}
CorrectionHaming(FGaucheRE,Fenetre,LargeurFourier);
CorrectionHaming(FCentreGRE,Fenetre,LargeurFourier);
ObjetCIma.FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);
ObjetCIma.FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);	
double DistanceGaucheFixeRed = ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,LargeurFourier);	


cpt =0;
for (i = IptInit-LargeurFourier;i<IptInit;i++)
{
	FGaucheRE[cpt]= (double)TabCoupeCoul[i].g; FGaucheIM[cpt]= 0.0f; cpt++; // bug RE 0104
}
cpt =0;
for (i = IptInit+1;i<IptInit+LargeurFourier+1;i++)//modif 09/01/2007
{
	FCentreGRE[cpt]= (double)TabCoupeCoul[i].g; FCentreGIM[cpt]= 0.0f; cpt++;
}
CorrectionHaming(FGaucheRE,Fenetre,LargeurFourier);
CorrectionHaming(FCentreGRE,Fenetre,LargeurFourier);
ObjetCIma.FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);
ObjetCIma.FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);	
double DistanceGaucheFixeGreen = ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,LargeurFourier);	

cpt =0;
for (i = IptInit-LargeurFourier;i<IptInit;i++)
{
	FGaucheRE[cpt]= (double)TabCoupeCoul[i].b; FGaucheIM[cpt]= 0.0f; cpt++; // bug RE 0104
}
cpt =0;
for (i = IptInit+1;i<IptInit+LargeurFourier+1;i++)//modif 09/01/2007
{
	FCentreGRE[cpt]= (double)TabCoupeCoul[i].b; FCentreGIM[cpt]= 0.0f; cpt++;
}
CorrectionHaming(FGaucheRE,Fenetre,LargeurFourier);
CorrectionHaming(FCentreGRE,Fenetre,LargeurFourier);
ObjetCIma.FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);
ObjetCIma.FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);	
double DistanceGaucheFixeBlue = ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,LargeurFourier);	
// DistanceGaucheFixe=DistanceGaucheFixeRed+DistanceGaucheFixeGreen+DistanceGaucheFixeBlue;
 DistanceGaucheFixe=DistanceGaucheFixeBlue;
 if (DistanceGaucheFixeRed>DistanceGaucheFixe) DistanceGaucheFixe=DistanceGaucheFixeRed;
 if (DistanceGaucheFixeGreen>DistanceGaucheFixe) DistanceGaucheFixe=DistanceGaucheFixeGreen;

#endif
//----------------------------------------------------
// recopie des tab pour fourier + calcul distance
//----------------------------------------------------

for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) //idem modif 21/11/2006
{
	int Edge=ImageEdge[TabCoupeXY[IptFinal].x+TabCoupeXY[IptFinal].y*ImageWidth];
	//Edge=1; // acceleration !!
	if (Edge!=0)
	{
		cpt =0;
		for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)
		{
			FDroiteRE[cpt]= (double)TabCoupeCoul[i].r; FDroiteIM[cpt]= 0.0f; cpt++;
		}
		cpt=0;
		for (i = IptFinal-LargeurFourier;i<IptFinal;i++)//modif 09/01/2007
		{
			FCentreDRE[cpt]= (double)TabCoupeCoul[i].r; FCentreDIM[cpt]= 0.0f; cpt++;
		}
		//Correction de Haming
		CorrectionHaming(FDroiteRE,Fenetre,LargeurFourier);
		CorrectionHaming(FCentreDRE,Fenetre,LargeurFourier);
		ObjetCIma.FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);
		ObjetCIma.FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);
		double DistanceDroiteRed=ObjetCIma.HyperDistanceModuleFourier(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,LargeurFourier);
		double DistanceBorSimRed=ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,LargeurFourier);
		cpt =0;
		for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)
		{
			FDroiteRE[cpt]= (double)TabCoupeCoul[i].g; FDroiteIM[cpt]= 0.0f; cpt++;
		}
		cpt=0;
		for (i = IptFinal-LargeurFourier;i<IptFinal;i++)//modif 09/01/2007
		{
			FCentreDRE[cpt]= (double)TabCoupeCoul[i].g; FCentreDIM[cpt]= 0.0f; cpt++;
		}
		//Correction de Haming
		CorrectionHaming(FDroiteRE,Fenetre,LargeurFourier);
		CorrectionHaming(FCentreDRE,Fenetre,LargeurFourier);
		ObjetCIma.FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);
		ObjetCIma.FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);
		double DistanceDroiteGreen=ObjetCIma.HyperDistanceModuleFourier(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,LargeurFourier);
		double DistanceBorSimGreen=ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,LargeurFourier);
		cpt =0;
		for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)
		{
			FDroiteRE[cpt]= (double)TabCoupeCoul[i].b; FDroiteIM[cpt]= 0.0f; cpt++;
		}
		cpt=0;
		for (i = IptFinal-LargeurFourier;i<IptFinal;i++)//modif 09/01/2007
		{
			FCentreDRE[cpt]= (double)TabCoupeCoul[i].b; FCentreDIM[cpt]= 0.0f; cpt++;
		}
		//Correction de Haming
		CorrectionHaming(FDroiteRE,Fenetre,LargeurFourier);
		CorrectionHaming(FCentreDRE,Fenetre,LargeurFourier);
		ObjetCIma.FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);
		ObjetCIma.FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);
		double DistanceDroiteBlue=ObjetCIma.HyperDistanceModuleFourier(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,LargeurFourier);
		double DistanceBorSimBlue=ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,LargeurFourier);

	//	TabDistanceDroite[cptTabD]=DistanceDroiteRed+DistanceDroiteGreen+DistanceDroiteBlue;
		TabDistanceDroite[cptTabD]=DistanceDroiteBlue;
		if (TabDistanceDroite[cptTabD]<DistanceDroiteRed) TabDistanceDroite[cptTabD]=DistanceDroiteRed;
		if (TabDistanceDroite[cptTabD]<DistanceDroiteGreen) TabDistanceDroite[cptTabD]=DistanceDroiteGreen;
			//ObjetCIma.HyperDistanceModuleFourier(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,LargeurFourier);
	//	TabBorSim[cptBord]=DistanceBorSimRed+DistanceBorSimGreen+DistanceBorSimBlue;
		TabBorSim[cptBord]=DistanceBorSimBlue;
		if (TabBorSim[cptBord]<DistanceBorSimRed)TabBorSim[cptBord]=DistanceBorSimRed;
		if (TabBorSim[cptBord]<DistanceBorSimGreen)TabBorSim[cptBord]=DistanceBorSimGreen;
		//	ObjetCIma.HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,LargeurFourier);
	} 
	else
	{
		TabDistanceDroite[cptTabD]=0.0;
		TabBorSim[cptBord]=0.0;
	}
	cptTabG++;
	cptTabD++;
	cptSection++;
	cptBord++;

} // fin for Iptfinal
//CStdioFile FichierVisu;
//CString ChaineVisu;
//FichierVisu.Open("c:\\cerclesbruite8.txt",CFile::modeCreate|CFile::modeWrite);
//
//for(i=0;i<TailleTabDist;i++)
//{
//	ChaineVisu.Format("%f\t",DistanceGaucheFixe);
//	FichierVisu.WriteString(ChaineVisu);
//	ChaineVisu.Format("%f\t",TabBorSim[i]);
//	FichierVisu.WriteString(ChaineVisu);
//	ChaineVisu.Format("%f\n",TabDistanceDroite[i]);
//	FichierVisu.WriteString(ChaineVisu);
//}
//FichierVisu.Close();

//-----------------------------------------------------------------------------------------------------
// recherche des max pour l'inversion des fonctions de probalilité "section homogène" et bord similaire
//-----------------------------------------------------------------------------------------------------
for(i=0;i<TailleTabDist;i++)
{	
	//if(TabHomogSection[i]>MaxHomoSec)MaxHomoSec = TabHomogSection[i];	
	if(TabBorSim[i]>MaxBordSim)MaxBordSim = TabBorSim[i];	
	//if(TabDistanceGauche[i]>MaxGauche)MaxGauche = TabDistanceGauche[i];	
	if(TabDistanceDroite[i]>MaxDroite)MaxDroite = TabDistanceDroite[i];
}
TabAllMaxi[0] = DistanceGaucheFixe;
TabAllMaxi[1] = MaxDroite;
TabAllMaxi[2] = MaxBordSim;
TabAllMaxi[3] = 0.0f;

for(i=0;i<4;i++)
	if (TabAllMaxi[i]>MaxNorma)	MaxNorma = TabAllMaxi[i];

if (MaxNorma == 0.0f)MaxNorma=1.0;
#ifndef DISTGAUCHE
DistanceGaucheFixe=MaxNorma;
#endif;
//--------------------------------------------------------
//			Fonction export + calcul des couts
//--------------------------------------------------------

i=0; // compteur
	for(IptFinal=(IptInit+(LargeurFourier)+1),i=0;IptFinal<IptFinalMax;IptFinal++,i++) // modif 21/11/2006 rajout pour limiter les surbalayage
	{
		double width = sqrt(
		((double)TabCoupeXY[IptInitFixe].x - (double)TabCoupeXY[IptFinal].x)*
		((double)TabCoupeXY[IptInitFixe].x - (double)TabCoupeXY[IptFinal].x)
		+
		((double)TabCoupeXY[IptInitFixe].y - (double)TabCoupeXY[IptFinal].y)*
		((double)TabCoupeXY[IptInitFixe].y - (double)TabCoupeXY[IptFinal].y));
		int g;
		if (WidthPrior>0)
		{
		 g = int( ( (width-WidthPrior)/TailleGauss)*10.0 )+45;  // 45
		 g = __max(0, __min(g, 90));
		} 
		else 
		g = 45;

		ValDist =(TabDistanceDroite[i]/MaxNorma) * 
			(DistanceGaucheFixe/MaxNorma)*
			((MaxNorma-TabBorSim[i])/MaxNorma)
				*gauss_table_sigma_1[g];///gauss_table_sigma_1[45];
		
		int Edge=ImageEdge[TabCoupeXY[IptFinal].x+TabCoupeXY[IptFinal].y*ImageWidth];
		//Edge=1;
		if(MaxiCurvi.ValDist<ValDist && Edge!=0)
		{
			MaxiCurvi.ValDist = ValDist;
			MaxiCurvi.Largeur = width;
			MaxiCurvi.IptInit = IptInit;
			MaxiCurvi.IptFinal = IptFinal;
			MaxiCurvi.Centre = (IptFinal+IptInit)/2;
		}

	}//fin IptFinal


}
