// CurviDetect1D.h: interface for the CCurviDetect1D class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CURVIDETECT1D_H__F29E1257_B2C6_46E2_A051_2ABDD0C76E57__INCLUDED_)
#define AFX_CURVIDETECT1D_H__F29E1257_B2C6_46E2_A051_2ABDD0C76E57__INCLUDED_

#include "Ima.h"	// Added by ClassView
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CCurviDetect1D  
{
public:
	CString ChaineDist;
	CStdioFile FichierStdIO;
	CIma ObjetCIma;
	CFile Fichier;
	void Quality(CString PathSignal, CString NFichierVisu, CString NFichierTotale, int LargeurFourier, int LongueurSignal);
	CCurviDetect1D();
	virtual ~CCurviDetect1D();

	// Initialisation du tableau contenant les valeur de la correction de Haming
	void InitFenetreHaming(int LargeurFiltre);
	// //Fenetre de correction de Haming
	double * Fenetre;
	void CorrectionHaming(double* Signal, double* FenetreHaming, long Largeur);
	void RespCurviMaxi(BYTE *  TabCoupe);
	CCurviDetect1D(int ValLongueurSignal, int ValLargeurFourier);
	int LongueurSignal;
	int LargeurFourier;
	int DemiLargeurFourier;
	int DemiLongueurSignal;
	int IptInit;
	int IptFinal;
	int LargeurGabaritMin;
	int LargeurGabaritMaxi;
	int IptInitMax;
	int IptFinalMax;
	int IptInitMin;
	int IptFinalMin;
	int cptTabG;
	int cptTabD;
	int cptSection;
	int cptBord;
	int cpt;
	double ValDist;
	

	int RacineLargeurFourier;
	double * FGaucheRE;
	double * FDroiteRE;
	double * FCentreGRE;
	double * FCentreDRE;


	double * FGaucheIM;
	double * FDroiteIM;
	double * FCentreGIM;
	double * FCentreDIM;

	double * FSDFixeRE;
	double * FSDMobileRE;

	double * FSDFixeIM;
	double * FSDMobileIM;

	double MaxGauche;
	double MaxDroite;
	double MaxHomoSec;
	double MaxBordSim;
	double MaxGauche2Passe;
	double MaxDroite2Passe;
	double MaxHomoSec2Passe;
	double MaxBordSim2Passe;
	double TempValDis2Passe;

	int TailleTabDist;

	double * TabDistanceGauche;
	double * TabDistanceDroite;
	double * TabHomogSection;
	double * TabBorSim;

	double * TabDisGlis;
	double * MuLTabDis;

	double * TabMaxVisu;

	
	RespCurvi * TabStoreTemp;
	RespCurvi * TabStoreLarg ;

	RespCurvi MaxiCurvi;
	RespCurvi SecondePasse;

	int * TabLargeur;
	CCurviDetect1D(int ValLongueurSignal, int ValLargeurFourier, int ValIptInit);
	void RespCurviMaxiFast(BYTE * TabCoupe, CPoint *TabCoupeXY, double WidthPrior, int IptInitFixe,int *ImageEdge,int ImageWidth);
	void RespCurviMaxiFastCoul(rgb * TabCoupeCoul, CPoint *TabCoupeXY, double WidthPrior, int IptInitFixe,int *ImageEdge,int ImageWidth);
};

#endif // !defined(AFX_CURVIDETECT1D_H__F29E1257_B2C6_46E2_A051_2ABDD0C76E57__INCLUDED_)
