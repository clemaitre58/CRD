#if !defined(AFX_DIALRELIEF_H__1C009B01_CA70_11D2_A59A_0020AF2B5B4E__INCLUDED_)
#define AFX_DIALRELIEF_H__1C009B01_CA70_11D2_A59A_0020AF2B5B4E__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
// Dialrelief.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialrelief dialog

class CDialrelief : public CDialog
{
// Construction
public:
	CDialrelief(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialrelief)
	enum { IDD = IDD_DIALRELIEF };
	float	m_dist;
	float	m_foc;
	float	m_tet;
	float	m_xo;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialrelief)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialrelief)
	virtual void OnOK();
	virtual void OnCancel();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DIALRELIEF_H__1C009B01_CA70_11D2_A59A_0020AF2B5B4E__INCLUDED_)
