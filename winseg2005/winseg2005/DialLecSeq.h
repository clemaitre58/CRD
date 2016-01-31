// DialLecSeq.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialLecSeq dialog
class CImageWnd;
class CDialLecSeq : public CDialog
{

	CImageWnd* m_pView;
	int Nbimseq;
	CString m_titre;
// Construction
public:
//	CDialLecSeq(CWnd* pParent = NULL);   // standard constructor
	CDialLecSeq(CImageWnd* pParent);
// Dialog Data
	//{{AFX_DATA(CDialLecSeq)
	enum { IDD = IDD_DIALSEQ };
	int		m_edit;
	//}}AFX_DATA
	
// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialLecSeq)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialLecSeq)
	afx_msg void OnButton1();
	afx_msg void OnButton2();
	afx_msg void OnButton3();
	afx_msg void OnTimer(UINT nIDEvent);
	afx_msg void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	virtual void OnOK();
	virtual void OnCancel();
	virtual BOOL OnInitDialog();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()

	void Finlect();
public:
	BOOL Create(int,CString);
	void Redefinition(int);
};
