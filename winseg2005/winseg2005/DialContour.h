// DialContour.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialContour dialog

class CDialContour : public CDialog
{
// Construction
public:
	CDialContour(unsigned char,unsigned char, long,long,CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialContour)
	enum { IDD = IDD_Hysteresis };
	long	m_seuilbas;
	long	m_seuilhaut;
	long	m_longmax;
	long	m_seuilcont;
	BOOL	m_check;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialContour)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialContour)
	virtual void OnCancel();
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
