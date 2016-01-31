#include "acqui.h" 

          
// cette classe hérite de CAcqImageWnd qui elle même hérite de CImageWnd                  
class CTraitAcqWnd : public CAcqImageWnd 
{  
    public:
	// constructeur
	CTraitAcqWnd(int,NBCOLOR);   
	// destructeur
	~CTraitAcqWnd();  

	protected:
	DECLARE_MESSAGE_MAP();
	
	public:  
	afx_msg void OnAutreTraitement();      
	// autre traitement   
	// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CTraitAcqWnd)
	//}}AFX_VIRTUAL

// Implementation
protected:
	// Generated message map functions
	//{{AFX_MSG(CTraitAcqWnd)
	afx_msg void OnAcqurirSquence3aquisitionsquencecalcul();
	afx_msg void OnUpdateAcqurirSquence3aquisitionsquencecalcul(CCmdUI* pCmdUI);
	//}}AFX_MSG

};  
          
