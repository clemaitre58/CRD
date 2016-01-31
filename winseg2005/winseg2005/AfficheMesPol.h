// AfficheMesPol.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CAfficheMesPol frame
#pragma once;

#include "PolMc.h"
#include "dialpara.h"
#include "image.h"
#include "mainfrm.h"

//class CGraph2Wnd;
#include "Graph2Wnd.h"
#include "dialparz.h"

extern class CDialParz; 
class CAfficheMesPol : public CMDIChildWnd
{
	DECLARE_DYNCREATE(CAfficheMesPol)
public:
	CAfficheMesPol();           // protected constructor used by dynamic creation
	CWnd * m_pView;
// Attributes
public:
	CGraph2Wnd * pParzen0,* pHisto0,* pHisto1,* pParzen1;
	long histo00[256];
	long parzen00[256];
	long histo01[256];
	long parzen01[256];
	long histo10[256];
	long parzen10[256];
	long histo11[256];
	long parzen11[256];
	CMainFrame* pFrame;

	CMenu menu; 
	CPolyMc PolyMc;
	BOOL Visu,ValidPol;
	CDialParz *DialParz;
	int param0, param1,Hparam;
	CString NomFichierPol;
	CDialPara *DialPara;
	BOOL Parzen,Graph;
	//////affichage
	CRect rcClient;
	int WGrD,WGrG,HGrB,HGrH;
	float Ax,Ay,Bx,By;

	BYTE class1,class2;
	double probc1,probc2;
	BOOL ValidCalcul;
// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAfficheMesPol)
	public:
	virtual BOOL Create( LPCTSTR lpszWindowName);
	virtual BOOL DestroyWindow();
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CAfficheMesPol();
	int Init(CWnd* mere,CString NomImil="defaut.mes",CString NomFichPol="defaut.pol",CDialPara *CPara=NULL);
	// Generated message map functions
	//{{AFX_MSG(CAfficheMesPol)
	afx_msg void OnFichierQuitter();
	afx_msg void OnVisualisationChoixdesparametres();
	afx_msg void OnVisualisationVisualiserlespolytopes();
	afx_msg void OnPaint();
	afx_msg void OnVisualisationInfomesures();
	afx_msg void OnVisualisationInfopolytopes();
	afx_msg void OnParzen();
	afx_msg void OnWindowClosegraphics();
	afx_msg void OnUpdateWindowClosegraphics(CCmdUI* pCmdUI);
	afx_msg void OnParzenCalculdetous();
	afx_msg void OnUpdateVisualisationVisualiserlespolytopes(CCmdUI* pCmdUI);
	afx_msg void OnCorrelationCalcul();
	afx_msg void OnUpdateCorrelationCalcul(CCmdUI* pCmdUI);
	afx_msg void OnCorrelationCalculdans();
	afx_msg void OnUpdateCorrelationCalculdans(CCmdUI* pCmdUI);
	afx_msg void OnVisualisationToolbar();
	afx_msg void OnUpdateVisualisationToolbar(CCmdUI* pCmdUI);
	afx_msg void OnDestroy();
	afx_msg void OnVisualisationChoixdunfichierdepolytopes();
	afx_msg void OnUpdateVisualisationInfopolytopes(CCmdUI* pCmdUI);
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnCalculsChoixdeclasses();
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnUpdateParzen(CCmdUI* pCmdUI);
	afx_msg void OnUpdateParzenCalculdetous(CCmdUI* pCmdUI);
	afx_msg void OnMDIActivate(BOOL bActivate, CWnd* pActivateWnd, CWnd* pDeactivateWnd);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
private:
	double Estimateur (int,int);
	void CalculParzen(int,long *,long *,char * );
	BOOL Definirh();
	double Covariance(int,int);
	double Esperance_xy(int,int);
	double Esperance(int);
	double Variance(int);
	double Esperance2(int);
	double Correlation(int,int);
	//affichage
	void ParametreGraphique();
	int PosGrX(int);
	int PosGrY(int);
protected:
	void AxeGraph(CPaintDC&);


};

/////////////////////////////////////////////////////////////////////////////
