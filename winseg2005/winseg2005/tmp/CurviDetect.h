// CurviDetect.h: interface for the CCurviDetect class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CURVIDETECT_H__7B6AAE6A_B4F6_47CA_B323_DC75F9EA6CEB__INCLUDED_)
#define AFX_CURVIDETECT_H__7B6AAE6A_B4F6_47CA_B323_DC75F9EA6CEB__INCLUDED_

#include "Ima.h"	// Added by ClassView
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "struct.h"
#include "svm.h"

/*struct SComplexe
{
	float Re;
	float Im;
};*/

struct SCross
{
	float Width;
	float Angle;
};

class CCurviDetect  
{
public :
	long CompteurFFT;
	float CoefGauss[5];
	rgb *pImage;
	int Largeur;
	int Hauteur;
	float pi,pis2,pi2;
	float **I_i;
	float **I_x;
	float **I_y;
	float **I_g;
	float **I_Qr;
	float **I_Qb;
	float **I_Qg;
	float **I_map;
	float **I_Ext;

	SCross **Cross;
	
	float **Gauss;

	float *CoGauss;

	float feature[7];
	float featureQmax[7];
	float featureQmax2[7];
	float CoefFeature[7];

	float QpipjMin;
	float QpipjMin2;

	POINT AxisFin,AxisFin1;


	int Mode;
	int NPAR;
	int wmin,wmax;
	int m_MaxWidth;
	int dw;


	struct svm_model* model;


	double Histo_Sig[100];
public:
	double EtudeMaxi;
	CString ChaineEcretage;
	long CompteurEcretage;
	CString ChaineDistance;
	CStdioFile FichierExportDistance;
	BOOL OpenDistance;
	void NormalisationSpectreComplexe(double * SpectreRe, double * SpectreIm, int NbValeur);
	long m_TailleTabCoupe;
	void NormalisationSpectre(double * Spectre, int NbValeur);
	int LargeurBordFourierSur2;
	double * SFourierGaucheRModule;
	double * SFourierGaucheGModule;
	double * SFourierGaucheBModule;
	double * SFourierDroiteRModule;
	double * SFourierDroiteGModule;
	double * SFourierDroiteBModule;
	double * SFourierSectionDroiteRModule;
	double * SFourierSectionDroiteGModule;
	double * SFourierSectionDroiteBModule;
	double * SFourierSectionGaucheRModule;
	double * SFourierSectionGaucheGModule;
	double * SFourierSectionGaucheBModule;
	double DistanceEuclPoint(POINT Val1, POINT Val2);
	int RacineLargeurBordFourier;
	CIma ObjetCIma;
	int LargeurBordFourier;
	double * SFourierSectionGaucheBim;
	double * SFourierSectionGaucheBre;

	double * SFourierSectionGaucheGim;
	double * SFourierSectionGaucheGre;

	double * SFourierSectionGaucheRim;
	double * SFourierSectionGaucheRre;

	double * SFourierSectionDroiteBim;
	double * SFourierSectionDroiteBre;

	double * SFourierSectionDroiteGim;
	double * SFourierSectionDroiteGre;

	double * SFourierSectionDroiteRim;
	double * SFourierSectionDroiteRre;

	double * SFourierDroiteBim;
	double * SFourierDroiteBre;

	double * SFourierDroiteGim;
	double * SFourierDroiteGre;

	double * SFourierDroiteRim;
	double * SFourierDroiteRre;
	
	double * SFourierGaucheBim;
	double * SFourierGaucheBre;

	double * SFourierGaucheGim;
	double * SFourierGaucheGre;
	
	double * SFourierGaucheRim;
	double * SFourierGaucheRre;
	//double Pi;
	//long puissance2(long p);
	//void FFT1D(long q, SComplexe  *TFFT1D_S);
	void ModuleComplex(double * resultat, double * Reel, double * Imaginaire, int DeuxPuissanceNFFT);
	//double ModuleComplexe(float a,float b);
	void Hough();
	void ComputeCoefGauss(int nbiter);
	void Info();
	void Init();
	void extrema(float seuil);
	void CrossSection(int l, int c, int Boucle);
	void WriteSample(CString Nom,float *feature, int t);
	void Quality(float *FeatureFin,float *pQ, int *pW, int *piw, POINT *Tab,int IPtInit,int l, int c);
	float Hue(rgb pImage);
	float Diff(float H2,float H1);
	float DiffAngle(POINTF p1, POINTF p2);
	void CrossSectionMap();
	void Test(rgb *pIma, unsigned char *Result,int Haut, int Larg);
	void Go(rgb *pIma, unsigned char *Result,int Haut, int Larg);
	void Go(rgb *pIma, unsigned char *Result,int Haut, int Larg,CPoint point);
	void RemplitTableauIJCoupe(POINT *TabIJ,int x,int y,float dx,float dy,int w,int dw,int Dir,int *IptInit);
	float AngleGradient(int l, int c);

	CCurviDetect();
	virtual ~CCurviDetect();

};

#endif // !defined(AFX_CURVIDETECT_H__7B6AAE6A_B4F6_47CA_B323_DC75F9EA6CEB__INCLUDED_)
