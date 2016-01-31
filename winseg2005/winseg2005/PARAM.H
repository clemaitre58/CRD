// param.h : header file
//
#ifndef __PARAM_H__     
#define __PARAM_H__  

#include "stdafx.h"

typedef unsigned char byte;
typedef unsigned short int word;
typedef unsigned int  dword;



int  comparaison(const void *,const void*);

class CParam
{	    

public:
	CParam();			// protected constructor used by dynamic creation

// Attributes
public:
	
// Operations
public:

 	BOOL LumPix(CString chaine, byte   *mem, byte   *imi, long Larg, long Haut);
 	BOOL LumMoyGlis(CString chaine, byte   *mem, byte   *imi, long Larg, long Haut);
	BOOL Lumirapide(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
 	BOOL LumMoyFix(CString chaine, byte   *mem, byte   *imi, long Larg, long Haut);
 	BOOL GradPix(CString chaine, byte   *mem, byte   *imi, long Larg, long Haut);
 	BOOL GradMoyGlis(CString chaine, byte   *mem, byte   *imi, long Larg, long Haut);
	BOOL GradMoyGlis2(CString chaine, byte   *mem, byte   *imi, long Larg, long Haut);
    BOOL GradMoyFix(CString chaine, byte   *mem, byte   *imi, long Larg, long Haut);
	BOOL Gradientrapide(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL GradientrapideNonNorm(CString chaine, byte  *mem, byte   *imi, long Largl, long Hautl);
	BOOL LumMoyDecLigFix(CString chaine, byte   *mem, byte   *imi, long Larg, long Haut);
	BOOL Entropie(CString chaine, byte   *mem, byte   *imi, long Larg, long Haut);
	BOOL EntropieGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL ContrastLoc(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL Max_Min(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL ContrastGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL ContrastrapideGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL Max_MinGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL MedianGlissant(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL MedianFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL SobelFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL SobelGlissant(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL RehaussementFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL RehaussementGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL Ouverture(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL OuvertureRapide(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL Fermeture(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL FermetureRapide(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL EcartmedFen(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL EcartmedFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL MoyTronqueglis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL MoyTronqueFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL MoyAdap(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL MoyAdapAuto(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL Maxrepglis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL NbMaxrepglis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL ImageEnregistre(byte *mem, long Largl, long Hautl);
	BOOL EcartTypeFenGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL MaxFenGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	BOOL DensMax33(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl);
	void MaximumGlissant(unsigned char *maximumc, unsigned char *imi,int Largl,int Hautl,int,int);
	void MinimumGlissant(unsigned char *Minimumc,unsigned char *imi,int Largl,int Hautl,int,int);
	void GradientRap(unsigned __int32 *lram,int large,int haute);
	void FenetreGlissanteRap(unsigned char* bram,unsigned __int16* cufen,unsigned __int16* lufen,int large,int haute,int LF,int HF);

	void Decode(CString,long *);
 	
// Implementation
public:
	~CParam();



	
};

#endif // __PARAM_H__