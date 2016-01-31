// polymcfloat.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CPolyMcFloat document version FLOAT

#ifndef __POLMCF_H__           
#define __POLMCF_H__


typedef unsigned char byte;
typedef unsigned short int word;
typedef unsigned int  dword;

#include "donneefloat.h"



class CPolyMcFloat 
{

public:
	CPolyMcFloat();			// protected constructor used by dynamic creation

// Attributes
public:
void	sauve_pol_avant(CString NomF);
int LirePol(CString NomIm);
void sauve_pol(CString NomF);
void calcb(int cl,long kmes)  ;
void calcb_all() ;
int interne(long k1,long kpol,int) ;
void etend(long k1,long kpol,int) ;
int externe(long k0,int) ;
int libre(int) ;
void remplace(long kpol,int);
void ajoute(long k1,int);
void calch(long k1,int);
void calch_all();
void polytope();
int Init(CDonneeFloat *);
void FinAlloc();
void Fin();
void AffPol();
void AffMes();

// Operations
public:
  byte NbClass;
  byte NPAR;
  byte texte[256];
  CString  NomPol;//,NomMes;
  short BORD,coef;
  unsigned short  NPOL;//,NMS0,NMS1; 
  long  kmax[256],NbMes[256],Npol[256];//,kmax1;
  
  float **Va;
  float **poli;
  float **pols;
  float **politemp;
  float **polstemp;

  float    *hypi;
  float    *hyps;
  
 
    
// Implementation
public:
	int Verifrecouvrtemp();
	int Verif_recouvr();
//	void VerifSauvMMC();
	BOOL Verif_redondance();
	~CPolyMcFloat();

		// Generated message map functions
protected:
	//{{AFX_MSG(CPolyMcFloat)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
//	DECLARE_MESSAGE_MAP()

private:
	int externfusion(long k0,int cl);
};

#endif