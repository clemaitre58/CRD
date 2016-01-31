// DialNbPointAppr.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "DialNbPointAppr.h"


/////////////////////////////////////////////////////////////////////////////
// CDialNbPointAppr dialog


CDialNbPointAppr::CDialNbPointAppr(CWnd* pParent /*=NULL*/)
	: CDialog(CDialNbPointAppr::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialNbPointAppr)
	m_nbpoint = 1000;
	//}}AFX_DATA_INIT
}


void CDialNbPointAppr::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialNbPointAppr)
	DDX_Text(pDX, IDC_EDIT1, m_nbpoint);
	DDV_MinMaxLong(pDX, m_nbpoint, 1, 32767);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialNbPointAppr, CDialog)
	//{{AFX_MSG_MAP(CDialNbPointAppr)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialNbPointAppr message handlers
