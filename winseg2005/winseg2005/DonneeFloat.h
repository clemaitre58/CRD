// DonneeFloat.h: interface for the CDonneeFloat class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DONNEE_H__4CB45F95_3EF7_11D5_879D_0050043AFEB2__INCLUDED_)
#define AFX_DONNEE_H__4CB45F95_3EF7_11D5_879D_0050043AFEB2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CDonneeFloat  
{
public:
	long NbMesTot;
	void SauveMinMax(CString NomFichier);
	void LitMinMax(CString NomFichier);
	float * Tmin;
	float * Tmax;
	void Normalise(float binf=0.1F,float bsup=0.9F);
	void NormaliseSimple();
	CDonneeFloat();
	long kmax[256];
	int LireDonnee(CString NomFichier);
	int Erreur;
	long NbMes[256];
	BYTE NPAR;
	BYTE NbClass;
	float ** Tab;
	CDonneeFloat(CString NomFichier);
	virtual ~CDonneeFloat();

};

#endif // !defined(AFX_DONNEE_H__4CB45F95_3EF7_11D5_879D_0050043AFEB2__INCLUDED_)
