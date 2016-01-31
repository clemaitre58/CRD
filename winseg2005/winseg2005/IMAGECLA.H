// imagecla.h : header file
//
#ifndef __IMAGEC_H__     
#define __IMAGEC_H__  
#include "timage.h"

/////////////////////////////////////////////////////////////////////////////
// CImageClas frame
class CImageClas : public CTraitImageWnd
{
	DECLARE_DYNAMIC(CImageClas)
public:
	//CImageClas();			// protected constructor used by dynamic creation
	CImageClas(HWND pere,LONG InitLarg,LONG InitHaut, TYPEIMAGE InitLeType);     
	CImageClas(HWND pere,CString,long,long,TYPEIMAGE);

// Attributes
public:
   HWND HFenetre;
// Operations
public:

// Implementation
protected:
	virtual ~CImageClas();

	// Generated message map functions
	//{{AFX_MSG(CImageClas)
	afx_msg void OnDestroy();
	afx_msg void OnClose();
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////
#endif
