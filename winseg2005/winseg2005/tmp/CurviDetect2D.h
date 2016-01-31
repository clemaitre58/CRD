#include "stdafx.h"
#include "image.h"
#include "ima.h"


class CCurviDetect2D
{
public:
	CCurviDetect2D(void);
	~CCurviDetect2D(void);
	BOOL SuiviContour(BYTE *  pImage,BYTE* pImageDest, long Largeur, long Hauteur);
	void BordAZero(BYTE* pImage, long Largeur, long Hauteur);
	CList <CPoint,CPoint> * ListePointContour;
	CArray <CList<CPoint,CPoint>*,CList<CPoint,CPoint>*> TabListePoint;
	int NbListeInArray ;
	CPoint CurentPoint;
	CPoint FirstPoint;
	double* m_TabDistance;
	long nbelement;
	long NbElementBack;
	POSITION pos;
	long Cmin;
	long Cmax;
	long Lmin;
	long Lmax;
	long Taille;
	BOOL SuiviFreeman(BYTE * pImage, long Largeur, long Hauteur,CMDIFrameWnd *pmd);
	CPoint NouveauPoint;
	void RemplitTableauCoupe(ImagePetrou A, CPoint * TabCoupe, int x, int y, float dx, float dy, int w, int * Centre, long Hauteur, long Largeur);
	void RemplitTableauIJCoupeBiDir(CPoint *TabIJ,int x,int y,float dx,float dy,int w,int dw,int Dir,int *IptInit, long Hauteur, long Largeur);
	BOOL TestCroisement(CPoint x1, CPoint x2, CPoint xp1, CPoint xp2);
	void InvertionPoint(CPoint * xp1, CPoint * xp2);
	void RechercheRectEnglobant(CPoint * TabPoint, long * Xmin, long * Xmax, long * Ymin, long * Ymax);
	void AfficheQuad(CImageWnd * pImage,CRgn * RegionAffichage,long Larg, long * Xmin, long * Xmax, long * Ymin, long * Ymax, int R, int G, int B);
	void MesureRepetabilite(rgb * Img1, rgb * Img2, long Largeur,long Hauteur,long NBRegionCurviImg1,long NBRegionCurviImg2,double * H);
	void MesureRepetJ(CIma * Ima1 , CIma * Ima2 , CString NomFichier, long NbReg1, long NbReg2);
	void RempliQuad(int* pPlan,CRgn * RegionAffichage,long Larg, long * Xmin, long * Xmax, long * Ymin, long * Ymax, int NumReg);
	int MesureRecouv(const int * Ima1 ,const int * Ima2 , int Larg,int Haut,long NbReg1, double SeuilRecouv, int *pi,int *TabListeOk);
	void FermetureContourSimple(BYTE * pImage, long Largeur, long Hauteur);
	double AngleCourbureLocale(Sfpoint P1,Sfpoint P2,Sfpoint P3);
};
