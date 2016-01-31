// DialSetT.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "DialSetT.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDialSetT dialog


CDialSetT::CDialSetT(CWnd* pParent /*=NULL*/)
	: CDialog(CDialSetT::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialSetT)
	m_S = 16;
	m_T = 16;
	//}}AFX_DATA_INIT
}


void CDialSetT::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialSetT)
	DDX_Text(pDX, IDC_EDIT1, m_S);
	DDX_Text(pDX, IDC_EDIT3, m_T);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialSetT, CDialog)
	//{{AFX_MSG_MAP(CDialSetT)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialSetT message handlers

void CDialSetT::OnCancel() 
{
	// TODO: Add extra cleanup here
	
	CDialog::OnCancel();
}

void CDialSetT::OnOK() 
{
	// TODO: Add extra validation here
	UpdateData();
	if((m_S<=0)||(m_T<0)) MessageBox("Valeurs impossibles");
	else CDialog::OnOK();
}

