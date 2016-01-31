#ifndef __DialFactZoom_H__           
#define __DialFactZoom_H__

#include "image.h"
/////////////////////////////////////////////////////////////////////////////
// CDialFactZoom dialog

class CDialFactZoom : public CDialog
{
// Construction
public:
	CDialFactZoom(CImageWnd * pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialFactZoom)
	enum { IDD = IDD_FACTEUR_ZOOM };
	int		m_Facteur_Zoom;
	CImageWnd*	m_pSet;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialFactZoom)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialFactZoom)
	virtual void OnCancel();
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#endif //CDialFactZoom