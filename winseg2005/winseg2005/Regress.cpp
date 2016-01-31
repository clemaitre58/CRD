// Regress.cpp: implementation of the CRegress class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "winseg.h"
#include "Regress.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CRegress::CRegress(POINT *Tab,int Taille)
{
	float Sxy;
	float Sx2;
	float Sy2;
	float mx,my,Num,Denom;
	int i;
	mx=0.0;
	my=0.0;
	for (i=0;i<Taille;i++)
	{
		mx+=Tab[i].x;
		my+=Tab[i].y;
	}
	mx/=Taille;
	my/=Taille;

	Sxy=0.0;
	Sx2=0.0;
	Sy2=0.0;

	for (i=0;i<Taille;i++)
	{
		Sxy+=(Tab[i].x-mx)*(Tab[i].y-my);
		Sx2+=(Tab[i].x-mx)*(Tab[i].x-mx);
		Sy2+=(Tab[i].y-my)*(Tab[i].y-my);
	}
	if (Sx2!=0.0) a=Sxy/Sx2; else a=1000.0;
	if(Sxy==0.0) {a=0.0;b=0.0;}
	else b=my-a*mx;

	Num=0;Denom=0;
	for (i=0;i<Taille;i++)
	{
		Num+=(b+a*Tab[i].x-my)*(b+a*Tab[i].x-my);
		Denom+=(Tab[i].y-my)*(Tab[i].y-my);
	}
	if (Denom!=0.0) r=Num/Denom; else r=1.0;
	if (Sy2==0.0) r=1.0;

}

CRegress::~CRegress()
{

}
