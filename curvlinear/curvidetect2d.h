#ifndef __CURVIDETECT2D_H__
#define __CURVIDETECT2D_H__

#include <list>
#include <vector>
#include "petrou.h"

typedef std::list <IPoint2D>  PointList;
typedef std::vector <PointList*>  TablePointList;

class CCurviDetect2D
{

public:
   TablePointList     TabListePoint;

	CCurviDetect2D(void);
	~CCurviDetect2D(void);
	bool SuiviFreeman(utls::BAry *Image);
	void RemplitTableauIJCoupeBiDirOld(IPoint2D *TabIJ,int x,int y,float dx,float dy,int w,int dw,int Dir,int *IptInit, long Hauteur, long Largeur);
   void RemplitTableauIJCoupeBiDir(IPoint2D *TabIJ,int x,int y,float dx,float dy,int w,int dw,int Dir,int *IptInit, long Hauteur, long Largeur);

private:
	void BordAZero(byte* pImage, int Largeur, int Hauteur);
   bool SuiviContour(byte* pImage, byte* pImageDest, int Largeur, int Hauteur);
   bool TestCroisement(IPoint2D x1, IPoint2D x2, IPoint2D xp1, IPoint2D xp2);
   void InvertionPoint(IPoint2D * xp1, IPoint2D * xp2);
};

#endif // __CURVIDETECT2D_H__
