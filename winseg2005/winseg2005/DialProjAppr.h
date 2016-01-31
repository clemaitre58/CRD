//{{AFX_INCLUDES()
#include "imgedit.h"
//}}AFX_INCLUDES
#if !defined(AFX_DIALPROJAPPR_H__68E186F1_B4FC_11D1_A4BD_00A0242E2BA0__INCLUDED_)
#define AFX_DIALPROJAPPR_H__68E186F1_B4FC_11D1_A4BD_00A0242E2BA0__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
// DialProjAppr.h : header file
//
#include "Timage.h"

/////////////////////////////////////////////////////////////////////////////
// CDialProjAppr dialog

class CDialProjAppr : public CDialog
{
public:
	CString FileProj,ImageAppr,FileAppr,FileParam;
protected:
	BOOL ValidImageAppr,ValidFileAppr,ValidFileParam;
	CTraitImageWnd* m_Parent;
// Construction
public:
	CDialProjAppr(CTraitImageWnd* pParent,CString initFileProj="");   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialProjAppr)
	enum { IDD = IDD_PROJET_APPR };
	CButton	m_sauvegarde;
	CButton	m_FichierProj;
	CButton	m_ImageAppr;
	CButton	m_FichierParam;
	CButton	m_FichierAppr;
//	CImgEdit	m_editimage;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialProjAppr)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialProjAppr)
	virtual void OnCancel();
	virtual void OnOK();
	virtual BOOL OnInitDialog();
	afx_msg void OnNomfichierappr();
	afx_msg void OnNomfichierparam();
	afx_msg void OnNomimageappr();
	afx_msg void OnNomprojet();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DIALPROJAPPR_H__68E186F1_B4FC_11D1_A4BD_00A0242E2BA0__INCLUDED_)
