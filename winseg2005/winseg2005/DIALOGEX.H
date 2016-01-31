
//#include "resource.h"
// dialogex.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialogExtrac dialog

class CDialogExtrac : public CDialog
{
// Construction
public:
	CDialogExtrac(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	//{{AFX_DATA(CDialogExtrac)
	enum { IDD = IDD_DIALOG3 };
	int		m_profondeur;
	//}}AFX_DATA

// Implementation
protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

	// Generated message map functions
	//{{AFX_MSG(CDialogExtrac)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
