
// dialpara.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialPara dialog
#ifndef __DIALPARA_H__           
#define __DIALPARA_H__
#include "resource.h"
class CDialPara : public CDialog
{
// Construction
public:
	CDialPara(CWnd* pParent = NULL);	// standard constructor
	
// Dialog Data
	//{{AFX_DATA(CDialPara)
	enum { IDD = IDD_PARAM };
	//}}AFX_DATA
	CString	m_dispo[256];
	CString	m_choisi[256];
	int nbpara;
	LRESULT IChoisi[256];
// Implementation
protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

	// Generated message map functions
	//{{AFX_MSG(CDialPara)
	afx_msg void OnDblclkList1();
	virtual void OnOK();
	virtual void OnCancel();
	virtual BOOL OnInitDialog();
	afx_msg void OnDblclkList2();
	afx_msg void OnAjouter();
	afx_msg void OnEnlever();
	//}}AFX_MSG

	DECLARE_MESSAGE_MAP()
};
/////////////////////////////////////////////////////////////////////////////
// CDialVal dialog

class CDialVal : public CDialog
{
	CString* Titre;
// Construction
public:
	CDialVal(CString* = NULL ,CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	//{{AFX_DATA(CDialVal)
	enum { IDD = IDD_DIALVAL };
	CString	m_val;
	//}}AFX_DATA

// Implementation
protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

	// Generated message map functions
	//{{AFX_MSG(CDialVal)
	virtual void OnOK();
	virtual void OnCancel();
	virtual BOOL OnInitDialog() ;
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
/////////////////////////////////////////////////////////////////////////////
// CDial2Val dialog
class CAfficheMesPol;
class CDial2Val : public CDialog
{
// Construction
public:
	CDial2Val(CDialPara*,int,int,CWnd* pParent = NULL);
	CDialPara * m_pDialPara;
	CComboBox * pCombo1,* pCombo2;
	int m_val0,m_val1;
// Dialog Data
	//{{AFX_DATA(CDial2Val)
	enum { IDD = IDD_DIAL2VAL };
	//}}AFX_DATA

// Implementation
protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

	// Generated message map functions
	//{{AFX_MSG(CDial2Val)
	virtual void OnCancel();
	virtual BOOL OnInitDialog();
	afx_msg void OnSelchangeCombo1();
	afx_msg void OnSelchangeCombo2();
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#endif
