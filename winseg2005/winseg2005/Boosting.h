// Boosting.h: interface for the CBoosting class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BOOSTING_H__11CB6184_E7CE_428F_BA67_8925181495ED__INCLUDED_)
#define AFX_BOOSTING_H__11CB6184_E7CE_428F_BA67_8925181495ED__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "donneefloat.h"
#include "polmcfloat.h"

struct SHyperplan
{
	float Borne;
	int Par;
	float Error;
	int Type;
};

struct SHyperInter
{
	float BorneInf;
	float BorneSup;
	float Error;
	int Par;
	int Type;
};

struct SBoostOutput
{
	int cf;
	float Foutput;
};


class CBoosting  
{
public:
	int TypeVariation;
	long ntot;
	float Complex;
	void InitMem(int m,int Tmaxi,int nbpar);
	float ComputeApparenteErrorForOneHyper(SHyperInter *pHr);
	void FindBestHyperrectangleWithMerging(int t);
	void FindBestHyperrectangle(int t);
	float ComputeApparenteError(int t);
	SHyperInter FindInterval(int par);
	void Lit(CString Nom);
	void Sauve(CString Nom);
	float TestError();
	SBoostOutput FinalClass(float *x);
	void Redistrib(CDonneeFloat *pDistrib,int t);
	void FindSetOfHyperrectangles(int t);
	void Smooth(float *p);
	void FindBestHyperplan(int t);
	void FindBestHyperInter(int t);
	void Iteration(int t);
	int Affine(float *p,float m,int indice);
	int Run(int iter,int cla=0);
	int MyPredict(int t,float *x);
	int model;
 	float RealError(int t);

 SHyperplan FindHyperplan(int par);
 float **Poids;
 float **Histo;
 float **HistoCumul;
 float **Sortie;

 BYTE **ht;

 float *DiffP;
 float *DiffM;
 float *MonAlpha;
 long *	erreurT;
 int *ModelRetenu;
 float *TabHardCost;
 float Error;
 float ErrorMini;
 float Epsilon;
 float Zt;

 int SMot;
 int Npar,NbClass;
 int Tmax;
 int IterFinal;
 int IterOptimum;
 int HardCost;

 SHyperInter *TheInter;
 SHyperInter *TheHr;
SHyperplan *TheHp;
 CPolyMcFloat **TabPoly;
 CPolyMcFloat Poly;

 CDonneeFloat *pMesDatas;

	CBoosting(int mod=0,int Tmaxi=1,CDonneeFloat *pData=NULL);
	CBoosting(int mod,int Tmaxi,int nbpar);
	virtual ~CBoosting();

};

#endif // !defined(AFX_BOOSTING_H__11CB6184_E7CE_428F_BA67_8925181495ED__INCLUDED_)
