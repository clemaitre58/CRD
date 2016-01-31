#if !defined(AFX_DIALNBPOINTAPPR_H__35D71FB1_BA6B_11D1_A4C4_00A0242E2BA0__INCLUDED_)
#define AFX_DIALNBPOINTAPPR_H__35D71FB1_BA6B_11D1_A4C4_00A0242E2BA0__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
// DialNbPointAppr.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialNbPointAppr dialog

class CDialNbPointAppr : public CDialog
{
// Construction
public:
	CDialNbPointAppr(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialNbPointAppr)
	enum { IDD = IDD_NBPOINTAPPR };
	long	m_nbpoint;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialNbPointAppr)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialNbPointAppr)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DIALNBPOINTAPPR_H__35D71FB1_BA6B_11D1_A4C4_00A0242E2BA0__INCLUDED_)
