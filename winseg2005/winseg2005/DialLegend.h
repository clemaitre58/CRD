#if !defined(AFX_DIALLEGEND_H__352A3654_C550_11D1_A4D2_00A0242E2BA0__INCLUDED_)
#define AFX_DIALLEGEND_H__352A3654_C550_11D1_A4D2_00A0242E2BA0__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
// DialLegend.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDialLegend dialog

class CDialLegend : public CDialog
{
	BOOL NonClass;
	int Nb;
// Construction
public:
	void Init(int , BOOL);
	CDialLegend(CWnd* pParent);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CDialLegend)
	enum { IDD = IDD_LEGENDE };
		// NOTE: the ClassWizard will add data members here
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDialLegend)
	public:
	virtual BOOL Create(int,BOOL);
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CDialLegend)
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnPaint();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
	
	COLORREF  Palette_Couleurs(int);
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DIALLEGEND_H__352A3654_C550_11D1_A4D2_00A0242E2BA0__INCLUDED_)
