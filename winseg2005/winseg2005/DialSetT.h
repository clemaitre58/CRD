// DialSetT.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialSetT dialog

class CDialSetT : public CDialog
{
// Construction
public:
	CDialSetT(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialSetT)
	enum { IDD = IDD_SetT };
	long	m_S;
	long	m_T;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialSetT)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialSetT)
	virtual void OnCancel();
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
