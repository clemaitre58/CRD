// dialeche.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialEchelle dialog


class CDialEchelle : public CDialog
{
// Construction
public:
	CDialEchelle(int InitEchelle = 1,CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	//{{AFX_DATA(CDialEchelle)
	enum { IDD = IDD_ECHELLE };
	int		m_Echelle;
	//}}AFX_DATA

// Implementation
protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

	// Generated message map functions
	//{{AFX_MSG(CDialEchelle)
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
