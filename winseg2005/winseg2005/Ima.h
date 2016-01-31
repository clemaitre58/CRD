// Ima.h: interface for the CIma class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_IMA_H__F12E6E8E_B6CD_4E1E_BA07_4F5D79A9B5CE__INCLUDED_)
#define AFX_IMA_H__F12E6E8E_B6CD_4E1E_BA07_4F5D79A9B5CE__INCLUDED_
#include <math.h>
#include "FFTReal.h"	// Added by ClassView
#include "BoostMcBis.h"
#include <afxtempl.h>
#include "struct.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000



// quelques structures
/////////////////////////////////////////////////////////////////////////////
struct BARY
{	
	double x;
	double y;
};

/////////////////////////////////////////////////////////////////////////////
typedef struct s_indicemin
{

	double min;
	unsigned long indice;

} indicemin;
/////////////////////////////////////////////////////////////////////////////
/**	\brief	Structure permettant de stocké une image couleur à trois composante : RGB.

*/
//struct rgb	{BYTE b;BYTE g;BYTE r;};

 /*! \class CIma
  *  \brief Fonction pour le traitement d'image
  *  \author Johel Miteran 
  *	 \author Fethi Smach 
  *  \author Cedric Lemaitre
  *  \version 1.0
  *  \date    2005
  *
  *  Fonction de traitement d'images, de caractésiation via des descripteurs (Fourier, Zernike) et encapsulation de la libsvm
  */
// j'ajoute structure complexe


//struct rgb	{BYTE b;BYTE g;BYTE r;};

class CIma  
{

public:

	rgb * m_PImage;
	BITMAPFILEHEADER * m_Tete;
	BITMAPINFO * m_Info;
	rgb *m_PImageReel;


	rgb* PImageCarre;

	_complex **fct;
	int N; // taille de l'image pour la FFT

	
	CList <CPoint,CPoint> listePointContour;
	CString CheminDepart;
	CPoint currentPoint;
	CPoint firstPoint;
	CPoint bufferPoint;
	CPoint Barycentre;
	CPoint Centre;
	CStringArray TabRefNom;
	CString m_Chaine;
	CString m_CheminDesModeles;

	int  nbelement;
	struct svm_model* model;

	double *m_TabDistance;

	int Distance;
	POSITION pos;
	BOOL bContour;

	int *m_TailleModele;
	double *TabC;
	double *m_TabModule;
	double **m_Tab_Modele;	// Tableau de pointeur
	BYTE	*m_Tab;
	double	*m_di_min;
	double	*m_db_final;
	int HistoSurf[256]; // pour les surfaces des zones (256 zones maxi)

	float m_min;
	float m_minMax;
	int m_identite;

	int m_Taille;
	int m_Hauteur;
	int m_Largeur;
	int NbEtiq;

	int m_i_Nb_Modeles;
	int m_TailleAnalyse;
	int Lmin,Lmax;
	int Cmin,Cmax;
	double Rayon;
	//FFTReal *pfft;
	//FFTReal *pfft2D;
	int m_TailleFFT;
	int LCarre,CCarre;
	int LFinal,CFinal;

	float * Tmin;
	float * Tmax;

	_complex cor(double lambda);
	int Nb_desc;
	int TailleMasque;
	int m_TailleFFT2D;
	int m_Npar;

	CBoostMcBis * pBoostMc;


public:
	double StandardDeviation(double * Tab, long LongueurTab);
	double HyperDistanceComplexe(double * V1re, double * V1im, double *V2re, double *V2im, long TailleVecteur);
	void HyperGradientGlissantZeroPadding(rgb *RGBrtn, rgb *pImage, int LargeurBloc, int HauteurBloc, int LargeurFenetre, int HauteurFenetre);
	double * DescripteursFourier1CompletGaussien(rgb *pImage, double sigma);
	void convolution_2D(float  *source,float  *dest,SImage info,SMasque masq);
	void HyperGradientDescripteurGlissant(rgb *RGBrtn,rgb *pImage,int LargeurBloc, int HauteurBloc,int LargeurFenetre,int HauteurFenetre);
	double * DescripteursFourier1Complet(rgb *pImage);
	void HyperGradientDescripteurBloc(rgb *RGBrtn,rgb *pImage,int LargeurBloc, int HauteurBloc,int LargeurFenetre,int HauteurFenetre);
	double HyperDistance(double *Vecteur1, double *Vecteur2, long TailleVecteur);
	rgb * HyperGradientDescripteurBloc(rgb *pImage, int LargeurBloc, int HauteurBloc,int LargeurFenetre,int HauteurFenetre);
	void FourierSurPointsContour(int SeuilG, int Apprend);
	int LireDats(CString Chemin);
	double * DescripteursFourier1(rgb *pImage);
	double * DescripteursFourier2(rgb *pImage);
	
	void Sobel(rgb *PImagel,long Largl,long Hautl);

	void InitPourFFT();
	void CalculCoefCorr();
	BYTE* DitherFloydSteinberg();
	void ReSampleMixteZernikeFourier(BYTE *Image, int Larg, int Haut, CClientDC *pdc,CString Nom,BOOL dec,BOOL Pos,BOOL Sauve,CString NomFich);
	void ReSampleZernike(BYTE *Image, int Larg, int Haut, CClientDC *pdc,CString Nom,BOOL dec,BOOL Pos,BOOL Sauve,CString NomFich);
	void zer_pol_R(int n, int m_in, double x, double y, double *res);
	void zer_pol(int n, int m, double x, double y, _complex *res);
	double * MomentsZernike(rgb *pImage);
	void ReSampleCarreBMPOffset();
	bool FFT(int dir,int m,double *x,double *y);
	void ReSampleCarreBMP();
	void SauveBMP(CString NomFich);
	void VisuBMP(CClientDC *pdc);
	void LitBMP(CString Nom);
	void LireCmc(CString NomCmc);
	void ChargeModelSVM(CString ChaineModel);
	void NormaliseAvecTable(float *tab, int Taille);
	void LitMinMax(CString NomFichier);
	void RecentrageGrilleCouleur(_complex ***grid, long HautIm, long LargIm);
	void ReSampleCarre(BYTE *Image, int Larg, int Haut,CClientDC *pdc,CString Nom,BOOL dec,BOOL Pos,BOOL Sauve,CString NomFich);
	void ReCherche(BYTE *Image, int Larg, int Haut);
	void IncrusteChaine(BYTE *pIm,int x, int y, CString ch,int num);
	void Masquage(CPoint p1, CPoint p2);
	void Fin();
	void AnalyseAvecBoosting();
	BOOL ChargeBoost(CString NomBoost);
	void AjouterDonneeDansRef(CString nomfichier);
	void CalculFFTProfile();
	BOOL m_Echelle;
	void LibereModeles();
	void Normalise(double *tab1, unsigned long t1);
	void Derive(float *tab1, unsigned long t1);
	void IncrusteBoite(BYTE *pIm,int num);
	void Expand(BYTE *pIm,long Taille);
	void Negatif(BYTE *pIm, long Taille);
	void InitAvecEtiquette(BYTE *pImage, int lar, int hau,int num);
	void LibereMemoire();
	void Init(BYTE *pImage, int lar, int hau,int Seuil);
	void VisuNB(BYTE *pImage);
	void BordAZero();
	void Analyse();
	void TabRotation(double *tab, unsigned long t);
	double sum_fabs(double *tab1, double *tab2, unsigned long t1, unsigned long t2);
	BOOL SuiviContour();
	void AjouterRef(CString nomfichier);
	indicemin Identification();
	indicemin db_min(double *tab, unsigned long t);
	void CalculDistance();
	BOOL ChargeModeles();
	void Masquage(BYTE *pImage,int Larg,int Haut, int c1,int l1,int c2, int l2);
	void Erosion(BYTE *pImage,int Larg,int Haut);
	void Seuillage(BYTE *pIm,long Taille,BYTE Seuil);
	long Etiquette(BYTE *pImage, int Larg, int Haut,int Seuil);
	void Contour(BYTE *pImage,int Larg,int Haut);
	CIma(BYTE *pImage,int lar,int hau,int Seuil);
	CIma();
	virtual ~CIma();

	double HyperDistanceModuleFourier(double * aR, double * aI, double * bR, double * bI, int LargFiltre);
};

#endif // !defined(AFX_IMA_H__F12E6E8E_B6CD_4E1E_BA07_4F5D79A9B5CE__INCLUDED_)
