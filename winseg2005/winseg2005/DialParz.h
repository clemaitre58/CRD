// DialParz.h : header file
//
#pragma once
#include "resource.h"
#include "AfficheMesPol.h"

/////////////////////////////////////////////////////////////////////////////
// CDialParz dialog
//class CAfficheMesPol;


class CDialParz : public CDialog
{
// Construction
CAfficheMesPol * m_pView;

public:
	CDialParz(CWnd* pParent = NULL);   // standard constructor
	CDialParz(CAfficheMesPol * );


// Dialog Data
	//{{AFX_DATA(CDialParz)
	enum { IDD = IDD_DIALPARZ };
	CString	m_recouvr0;
	CString	m_recouvr1;
	CString	m_param0;
	CString	m_param1;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialParz)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialParz)
	virtual void OnCancel();
	virtual BOOL OnInitDialog();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	BOOL Create(CString , char* ,CString , char* );
	void Rafraichissement(CString ,char *,CString , char* );
};
