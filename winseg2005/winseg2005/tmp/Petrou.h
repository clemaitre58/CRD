// Petrou.h: interface for the CPetrou class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PETROU_H__60C675A3_CBC5_4990_9F78_F12812CEB13D__INCLUDED_)
#define AFX_PETROU_H__60C675A3_CBC5_4990_9F78_F12812CEB13D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "stdafx.h"






class CPetrou  
{
public:
	float m_MaxReg;
	int IsEdge(int l,int c,int dir);
	void Nettoie();
	float AngleGradient(int l,int c);
	float DiffAngle(POINTF p1, POINTF p2);
	bool SuiviInitial(int *pl, int *pc,int lmax);
	float Distance(float *ta,float *tb, int n);
	void ProfilOblique();
	ImagePetrou * I_x;
	ImagePetrou * I_y;
	ImagePetrou * A;
    ImagePetrou * B;

	int		m_TE;
	int		m_MaxLength;
	int		m_MaxWidth;
	int		m_MaxDw;
	int 	m_MinLength;

	BOOL	m_CheckDw;
	BOOL	m_Check_Reg;


	CPetrou();
	virtual ~CPetrou();
	void derive_X(const ImagePetrou* A, ImagePetrou* B);
	void derive_Y(const ImagePetrou* A, ImagePetrou* B);
	void lissage(const ImagePetrou* A, ImagePetrou* B, int w, int d);
	void filtre_petrou(ImagePetrou* A, ImagePetrou* B, int w, int seuil, int d);
	void extrema(ImagePetrou* B, ImagePetrou* I_x,ImagePetrou* I_y, int seuil, int d, int w);
	void dynamique(ImagePetrou *ImagePetrou,int H,int L);

	void Go(BYTE *pImagePetrou, BYTE *dest,int Hauteur,int Largeur,int seuil,int d);
	void RemplitTableauIJCoupe(POINT *TabIJ,int x, int y,float dx,float dy,int w,int *centre);
	void filtre_petrou_couleur(ImagePetrou* R,ImagePetrou* G,ImagePetrou* B, ImagePetrou* I, int w, int seuil, int d);


	CPetrou(ImagePetrou * ImageA, ImagePetrou * ImageB,ImagePetrou * ImageGradX, ImagePetrou * ImageGradY);
};

#endif // !defined(AFX_PETROU_H__60C675A3_CBC5_4990_9F78_F12812CEB13D__INCLUDED_)
