
#ifndef __LHFLOAT_H__     
#define __LHFLOAT_H__   

// lhfloat.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// LHfloat dialog
#include "resource.h"
class LHfloat : public CDialog
{  
	DWORD TailleFich;
// Construction
public:
	LHfloat(DWORD InitTailleFich,CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	//{{AFX_DATA(LHfloat)
	enum { IDD = IDD_LHFLOAT };
	long	m_HAUT;
	long	m_LARG;
	//}}AFX_DATA

// Implementation
protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

	// Generated message map functions
	//{{AFX_MSG(LHfloat)
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
#endif // __LHFLOAT_H__

#ifndef __LHBYTE_H__     
#define __LHBYTE_H__   

class LHByte : public CDialog
{  
	DWORD TailleFich;
// Construction
public:
	LHByte(DWORD InitTailleFich,CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	//{{AFX_DATA(LHfloat)
	enum { IDD = IDD_LHFLOAT };
	long	m_HAUT;
	long	m_LARG;
	//}}AFX_DATA

// Implementation
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

	// Generated message map functions
	//{{AFX_MSG(LHfloat)
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#endif // __LHBYTE_H__
