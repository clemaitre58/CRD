#if !defined(AFX_TREECONFIG_H__0ED37A82_B7F9_11D1_A4BF_00A0242E2BA0__INCLUDED_)
#define AFX_TREECONFIG_H__0ED37A82_B7F9_11D1_A4BF_00A0242E2BA0__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
// TreeConfig.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CTreeConfig window

class CTreeConfig : public CTreeCtrl
{
// Construction
public:
	BYTE Nbpar;
	CString FichParam;
	BOOL DRAG;
	CImageList *imDrag;
	HTREEITEM Moving;
	CStringArray TabString;
	CDWordArray  TabDWord;

public:
	CTreeConfig();

// Attributes
public:

// Operations
	void AddDataItem(HTREEITEM myItem,CString initData);
	void RemoveDataItem(HTREEITEM myItem);
	void TriClasse();
	void AjouterClass();
	void Ajouterfich();
	void Supprimclas();
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CTreeConfig)
	//}}AFX_VIRTUAL

// Implementation
public:
	CString GetNameItem(HTREEITEM);
	int GetNBchild(HTREEITEM);

	virtual ~CTreeConfig();

	// Generated message map functions
protected:
	//{{AFX_MSG(CTreeConfig)
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnBegindrag(NMHDR* pNMHDR, LRESULT* pResult);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnNcMouseMove(UINT nHitTest, CPoint point);
	afx_msg void OnDblclk(NMHDR* pNMHDR, LRESULT* pResult);
	afx_msg void OnRButtonDblClk(UINT nFlags, CPoint point);
	//}}AFX_MSG

	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_TREECONFIG_H__0ED37A82_B7F9_11D1_A4BF_00A0242E2BA0__INCLUDED_)
