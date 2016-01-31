// Vue3d.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CVue3d frame

class CVue3d : public CMDIChildWnd
{
	DECLARE_DYNCREATE(CVue3d)
public:
	CVue3d(unsigned char *PImage=NULL, long Larg=256, long Haut=256,CWnd * pView=NULL);           // protected constructor used by dynamic creation

// Attributes
protected:
CWnd * m_pView;
CMenu menu;
unsigned char * PIm;
long Largeur,Hauteur;
int	PosScroll ;
int PosHFen ;	

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CVue3d)
	public:
	virtual BOOL Create(LPCTSTR szTitle);
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CVue3d();

	// Generated message map functions
	//{{AFX_MSG(CVue3d)
	afx_msg LONG OnFinVue(UINT uiParam,LONG lparam );
	afx_msg void OnPaint();
	afx_msg void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnDestroy();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////
