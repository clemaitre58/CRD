#ifndef __ProfilWnd_H__           
#define __ProfilWnd_H__

// ProfilWnd.h : header file
//
#include "graphwnd.h"

/////////////////////////////////////////////////////////////////////////////
// CProfilWnd frame
class CProfilWnd : public CGraphWnd
{
	DECLARE_DYNCREATE(CProfilWnd)
protected:
	CProfilWnd();           // protected constructor used by dynamic creation
public:
	CProfilWnd(int  * InitTab, int Nb, char * InitChaineX,char * InitChaineY, CWnd * InitFeneAppel = NULL, CPoint Ori = NULL, CPoint Fin = NULL, BOOL InitGraphProfil = FALSE);
	CProfilWnd( long  * InitTab, int Nb,char * InitChaineX,char * InitChaineY, CWnd * InitFeneAppel = NULL, CPoint Ori = NULL, CPoint Fin = NULL, BOOL InitGraphProfil = FALSE);
	CProfilWnd( float  * InitTab, int Nb,char * InitChaineX,char * InitChaineY, CWnd * InitFeneAppel = NULL, CPoint Ori = NULL, CPoint Fin = NULL, BOOL InitGraphProfil = FALSE);

private:
	CWnd * FeneAppel;
	CPoint PointOri,PointFin;
// Attributes
public:
	void SetPositionCoupe(long ,long ,long ,long );  
	void GetInfo(long * ,long * ,long * ,long * ,long * );   
private:
	virtual void AfficheDimCoupe();

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CProfilWnd)
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CProfilWnd();

	// Generated message map functions
	//{{AFX_MSG(CProfilWnd)
	afx_msg void OnDestroy();
	afx_msg LONG OnFinProfil(UINT uiParam,LONG lparam );
	afx_msg void OnChildActivate(); 
	afx_msg void OnMinMaxNonAuto();
	afx_msg void OnUpdateMinMaxAuto(CCmdUI* pCmdUI);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#endif //ProfilWnd.h
/////////////////////////////////////////////////////////////////////////////
