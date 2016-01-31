#pragma once
#include "afxcmn.h"


// Boîte de dialogue CFenetreProgressCurvi

class CFenetreProgressCurvi : public CDialog
{
	DECLARE_DYNAMIC(CFenetreProgressCurvi)

public:
	CFenetreProgressCurvi(CWnd* pParent = NULL);   // constructeur standard
	virtual ~CFenetreProgressCurvi();

// Données de boîte de dialogue
	enum { IDD = IDD_DIALOG_PROGRESSCURVI };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // Prise en charge de DDX/DDV

	DECLARE_MESSAGE_MAP()
public:
	double m_TPetrou;
	double m_TCTracking;
	double m_T1DProcess;
	double m_T2DProcess;
	double m_Affichage;
	CString m_TexteValProgress;
	CProgressCtrl m_ControlBarProgressCurvi;
	afx_msg void OnDestroy();
protected:
	virtual void PostNcDestroy();
public:
	CString m_TxtNbReg;
};
