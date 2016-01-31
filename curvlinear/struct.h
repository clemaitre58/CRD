#ifndef __STRUCT_H__     
#define __STRUCT_H__     

#include <list>
#include <vector>
#include <point.h>
#include <dtypes.h>
#include <ary.h>

enum NBCOLOR {RIEN=0,NB,COLOR};
enum TYPEIMAGE {BMP, MEM, FLOAT32};

enum LETYPE {IMAGE_FIXE,IMAGE_CONTINUE};
       
typedef unsigned char byte;

struct IPoint2D
{
   int x, y;
   IPoint2D() {};
   IPoint2D(int v): x(v), y(v) {}
   IPoint2D(int ax, int ay): x(ax), y(ay) {}
   bool operator ==(IPoint2D &other) { return x==other.x && y==other.y; }
   bool operator !=(IPoint2D &other) { return !(x==other.x && y==other.y); }
};

struct FPoint2D
{
   float x, y;

};

struct DPoint2D
{
   double x, y;

};

struct rgb
{
	unsigned char b;
	unsigned char g;
	unsigned char r;
};

struct drgb	
{
	double b;
	double g;
	double r;
};    

struct SMasque
{  
	float *Filtre;
	int Larg;
	int Haut;
	int Centre_col;
	int Centre_lig;
};

struct SImage
{ 
 	long Larg;
 	long Haut;                   
 	long Lmin,Lmax,Hmin,Hmax;
  	unsigned char  * Image; 
};

struct SFiltre
{
 float *Filtre;
 int Centre;
 int Fin;
}; 

struct SImageFloat
{ 
 	long Larg;
 	long Haut;                   
 	long Lmin,Lmax,Hmin,Hmax;
  	float  * Image; 
};

struct SImageRGB
{ 
 	long Larg;
 	long Haut;                   
 	long Lmin,Lmax,Hmin,Hmax;
  	rgb  * Image; 
};

struct SComplexe
{
	float Re;
	float Im;
};

struct STabIndice
{
	double r;
	int i;
};

struct ConnexityDegreeAFreemanCode
{
	unsigned char ConnexityDegree;
	bool FreemanCode[8];
};

struct ParcoursFreemanHoraire
{
	char Code[8];
};

struct ParcoursFreemanAntiHoraire
{
	char Code[8];
};

struct DirectPrecedent
{
	char Code[8];
};

struct PosDetectCurvi
{
	IPoint2D PositionMaxi;
	IPoint2D PosInit;
	IPoint2D PosFinal;
	double Largeur;
	double ValDist;
   DPoint2D PosAxe;
};

struct CurvilinearRegion;
typedef std::vector<CurvilinearRegion *> TableOfCurvilinearRegions;
typedef std::vector<size_t> IndexVector;

struct CurvilinearRegion : public std::list<PosDetectCurvi>
{   
   DPoint2D bbmin, bbmax; 
   double awidth, bbsize;
   bool merged;
   
   std::vector<size_t> OverlappingRegions;
   
   // computes bounding box out of axis points in the region
   void UpdateAxisBoundingBox()
   {
      iterator it = begin();
      if (it!=end())
      {
         bbmin.x = bbmax.x = it->PosAxe.x;
         bbmin.y = bbmax.y = it->PosAxe.y;
         awidth = it->Largeur;
      } else {
         // nothing to do, invalid BB
         awidth = bbmin.x = bbmin.y = bbmax.x = bbmax.y = bbsize = 0;
         return;
      }

      for (++it; it != end(); ++it)
      {
         if (bbmin.x > it->PosAxe.x) bbmin.x = it->PosAxe.x;
         if (bbmin.y > it->PosAxe.y) bbmin.y = it->PosAxe.y;
         if (bbmax.x < it->PosAxe.x) bbmax.x = it->PosAxe.x;
         if (bbmax.y < it->PosAxe.y) bbmax.y = it->PosAxe.y;
         awidth += it->Largeur;
      }
      
      awidth /= size();      
      bbsize = utls::max( bbmax.x - bbmin.x, bbmax.y - bbmin.y);
   }
};

struct RespCurvi
{
	int IptInit;
	int IptFinal;
	double ValDist;
	double Largeur;
	int Centre;
};

/* small template hack to deduce double from float and RGBValue<double> from RGBValue<float> as accumulator types */
template <typename Accumulated>
struct AccumulatorDeduction
{
   typedef Accumulated item_value;
};

template <> struct AccumulatorDeduction<float>
{
   typedef double value;
};

template <> struct AccumulatorDeduction<utls::RGBValue<float> >
{
   typedef utls::RGBValue<double> value;
};


typedef utls::Ary<float, float, double> Image;
typedef utls::Ary<utls::RGBValue<float>, float, utls::RGBValue<double> > RGBImage;

#endif // __STRUCT_H__