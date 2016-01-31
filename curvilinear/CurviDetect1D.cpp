// CurviDetect1D.cpp: implementation of the CCurviDetect1D class.
//
//////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <cmath>
#include <elemfun.h>
#include <timeutls.h>
#include "curvidetect1d.h"
#include <assert.h>

// a=gauss1d(1,0,10,0.1); fprintf('%25.15g,%25.15g,%25.15g,\n',a);
// spacing 0.1, center at 45
// v = v/0.1 + 45; v = max(0, min(v, 90));
// p = gauss_table_sigma_1[v];
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

double HyperDistanceModuleFourier(double * aR, double * aI, double * bR, double * bI, int LargFiltre)
{
	int k;
	double Ma,Mb;
	double Distance=0.0;

	for (k=0;k<LargFiltre;k++)
	{
		Ma = sqrt(aR[k]*aR[k]+aI[k]*aI[k]);
		Mb = sqrt(bR[k]*bR[k]+bI[k]*bI[k]);
		Distance+=(Ma-Mb)*(Ma-Mb);
	} // for k
	Distance = sqrt(Distance);
	Distance /= LargFiltre;
	return Distance;
}

double HyperDistanceModuleFourier(utls::RGBValue<double> * aR, utls::RGBValue<double> * aI, 
                                  utls::RGBValue<double> * bR, utls::RGBValue<double> * bI, int LargFiltre)
{
	int k;
	utls::RGBValue<double> Ma, Mb;
	utls::RGBValue<double> Distance(0.0);

	for (k=0;k<LargFiltre;k++)
	{
		Ma = sqrt(aR[k]*aR[k]+aI[k]*aI[k]);
		Mb = sqrt(bR[k]*bR[k]+bI[k]*bI[k]);
		Distance+=(Ma-Mb)*(Ma-Mb);
	} // for k
   double MaxDistance = 0;
   for (int i=0;i<3;i++)
   {
      if (Distance[i]>MaxDistance)
         MaxDistance = Distance[i];
   }

	MaxDistance=sqrt(MaxDistance);
	MaxDistance/=LargFiltre;
	return MaxDistance;
}


const int OversampleHamming = 25;

// Initialisation du tableau contenant les valeur de la correction de Haming
void InitFenetreHaming(double *Fenetre, int LargeurFiltre)
{
   int k, w = LargeurFiltre * OversampleHamming - 1;
	for(k=0;k<LargeurFiltre;k++)
	   Fenetre[k] = 0;
	double x, deuxpi = M_PI*2.0;
   // oversample Hamming window and average afterwards
	for (k=0; k <= w;k++)
	{
		x=(deuxpi*k)/double(w);
		Fenetre[k/OversampleHamming] += (0.5-0.5*cos(x))/double(OversampleHamming);
	} // k	
}
