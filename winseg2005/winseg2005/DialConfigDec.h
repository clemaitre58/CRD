#if !defined(AFX_DIALCONFIGDEC_H__52BF7E52_C007_11D1_A4CD_00A0242E2BA0__INCLUDED_)
#define AFX_DIALCONFIGDEC_H__52BF7E52_C007_11D1_A4CD_00A0242E2BA0__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
// DialConfigDec.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialConfigDec dialog

class CDialConfigDec : public CDialog
{
// Construction
public:
	CDialConfigDec(CWnd* pParent = NULL,CString initpmc="",CString initParam="",BOOL initcheck=FALSE,int initclasse=0);   // standard constructor

public:
	CString NomFichPMC,NomFichParam;
	BOOL ValidPMC,ValidPARAM;
	BYTE NBclass,NBpar;

// Dialog Data
	//{{AFX_DATA(CDialConfigDec)
	enum { IDD = IDD_CONFIDEC };
	CButton	m_checkbox;
	CEdit	m_edit1;
	CButton	m_ok;
	CSpinButtonCtrl	m_spin;
	CButton	m_param;
	CButton	m_pmc;
	BOOL	m_check;
	int		m_classeaffect;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialConfigDec)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialConfigDec)
	afx_msg void OnButton1();
	afx_msg void OnButton2();
	afx_msg void OnCheck1();
	virtual BOOL OnInitDialog();
	afx_msg void OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnChangeEdit1();
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DIALCONFIGDEC_H__52BF7E52_C007_11D1_A4CD_00A0242E2BA0__INCLUDED_)
