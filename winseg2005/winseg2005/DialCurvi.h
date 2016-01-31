#if !defined(AFX_DIALCURVI_H__9E5AA350_CEC4_4DA6_A1F4_FE4C5A59E632__INCLUDED_)
#define AFX_DIALCURVI_H__9E5AA350_CEC4_4DA6_A1F4_FE4C5A59E632__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// DialCurvi.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialCurvi dialog

class CDialCurvi : public CDialog
{
// Construction
public:
	int Larg;
	CDialCurvi(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialCurvi)
	enum { IDD = IDD_DIALOG_CURVI };
	int		m_TE;
	int		m_MaxLength;
	int		m_MaxWidth;
	int		m_MaxDw;
	BOOL	m_CheckDw;
	int		m_MinLength;
	CString	m_TxtMaxW;
	float	m_Reg;
	BOOL	m_Check_Reg;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialCurvi)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialCurvi)
	afx_msg void OnChangeEditMw();
	virtual BOOL OnInitDialog();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DIALCURVI_H__9E5AA350_CEC4_4DA6_A1F4_FE4C5A59E632__INCLUDED_)
