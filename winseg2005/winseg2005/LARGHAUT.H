#ifndef __larghaut_H__           
#define __larghaut_H__
// larghaut.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CLargHaut dialog

class CLargHaut : public CDialog
{                
	
// Construction
public:
	CLargHaut(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	//{{AFX_DATA(CLargHaut)
	enum { IDD = IDD_LARGHAUT };
	CString	m_TITRE;
	int		m_HAUT;
	int		m_LARG;
	int		COLOR;
	//}}AFX_DATA
	//BOOL COLOR;

// Implementation
protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

	// Generated message map functions
	//{{AFX_MSG(CLargHaut)
	virtual void OnOK();
	afx_msg void OnColor();
	afx_msg void OnNb();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#endif //larghaut