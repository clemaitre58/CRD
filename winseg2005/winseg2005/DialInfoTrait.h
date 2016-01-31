#ifndef __DialInfoTrait_H__           
#define __DialInfoTrait_H__

// DialInfoTrait.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialInfoTrait dialog

class CDialInfoTrait : public CDialog
{
// Construction
public:
	CDialInfoTrait(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialInfoTrait)
	enum { IDD = IDD_INFO_TRAITEMENT };
	CString	m_Comment;
	CString	m_Nom;
	//}}AFX_DATA


public:
	void UpdateInfo(CString, CString);
	virtual BOOL Create(CString InitComment, CString InitNom);
// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialInfoTrait)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialInfoTrait)
	virtual BOOL OnInitDialog();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#endif //DialInfoTrait.h