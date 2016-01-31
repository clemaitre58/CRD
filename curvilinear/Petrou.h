// Petrou.h: interface for the CPetrou class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PETROU_H__60C675A3_CBC5_4990_9F78_F12812CEB13D__INCLUDED_)
#define AFX_PETROU_H__60C675A3_CBC5_4990_9F78_F12812CEB13D__INCLUDED_

#include <ary.h>
#include "struct.h"
/*
typedef float Pixelp;
typedef Pixelp* Pixelpp;
struct ImagePetrou 
{
  Pixelp **ima;
  int dim_x;
  int dim_y;
} ;
*/

typedef Image    ImagePetrou;
typedef RGBImage ImagePetrouRGB;

class CPetrouParams
{

};

class CPetrou  
{
public:
   CPetrou(const ImagePetrou* Input, int w, int seuil, int d);
   CPetrou(const ImagePetrouRGB* Input, int w, int seuil, int d);
   ~CPetrou();

   ImagePetrou * I_x;
	ImagePetrou * I_y;
   ImagePetrou * Output;
	
private:
   void compute(const ImagePetrou* A, int w, int seuil, int d);
   void compute(const ImagePetrouRGB* RGB, int w, int seuil, int d);

   void extrema(ImagePetrou* B, ImagePetrou* I_x,ImagePetrou* I_y, int seuil, int d, int w);
	void dynamique(ImagePetrou *ImagePetrou,int H,int L);
};

#endif // !defined(AFX_PETROU_H__60C675A3_CBC5_4990_9F78_F12812CEB13D__INCLUDED_)
