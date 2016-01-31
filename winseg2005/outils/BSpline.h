//------------------------
// outils BSpline de base
//------------------------
#include <windows.h>
#include "struct.h"
#include <math.h>

void OnConvBspline(SImageFloat im,SImageFloat dest);
void OnConvBsplineFloat(SImage im,SImageFloat dest);
void OnGradientBspline(SImageFloat im,SImageFloat dest);
void OnDerDirSec(SImageFloat im,SImageFloat dest);
void OnLigneZero(SImageFloat im,SImage dest);
void OnDecompo(SImageFloat im,SImageFloat dest,int prof);
void OnRecomp(SImageFloat im,SImageFloat dest, int prof);
void OnTraitOnd(SImageFloat im,SImageFloat dest,int prof,float *m_a);