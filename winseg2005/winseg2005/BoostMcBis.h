// BoostMcBis.h: interface for the CBoostMcBis class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BOOSTMCBIS_H__49C54486_83C2_46DA_B5FE_F29FF32BF3DA__INCLUDED_)
#define AFX_BOOSTMCBIS_H__49C54486_83C2_46DA_B5FE_F29FF32BF3DA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "DonneeFloat.h"
#include "Boosting.h"	// Added by ClassView

class AFX_EXT_CLASS CBoostMcBis   
{
public:
	int TypeVariation;
	long ntot;
	SBoostOutput FinalClass(float *x);
	void Run(int iter);
	CBoostMcBis(int mod,int Tmaxi,CDonneeFloat *pData);
	virtual ~CBoostMcBis();
	int NbCombi;
	void Lit(CString Nom);
	void Sauve(CString Nom);
	float TestError();
	CBoostMcBis(CDonneeFloat *pData=NULL);
	int Npar,NbClass;
	int Tmax;
	int IterFinal;
	int Model,HardCost;
	CBoosting **pBoost;
	CDonneeFloat *pMesDatas;
	long *erreurT;


};

#endif // !defined(AFX_BOOSTMCBIS_H__49C54486_83C2_46DA_B5FE_F29FF32BF3DA__INCLUDED_)
