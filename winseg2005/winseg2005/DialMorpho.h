#if !defined(AFX_DIALMORPHO_H__F5D63E31_CB08_11D2_A59A_0020AF2B5B4E__INCLUDED_)
#define AFX_DIALMORPHO_H__F5D63E31_CB08_11D2_A59A_0020AF2B5B4E__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
// DialMorpho.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialMorpho dialog

class CDialMorpho : public CDialog
{
// Construction
public:
	CDialMorpho(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialMorpho)
	enum { IDD = IDD_DIALMORPHO };
	int		m_taillestructurant;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialMorpho)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialMorpho)
	virtual void OnOK();
	virtual void OnCancel();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DIALMORPHO_H__F5D63E31_CB08_11D2_A59A_0020AF2B5B4E__INCLUDED_)
