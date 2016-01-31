// DialParz.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "DialParz.h"
#include "AfficheMesPol.h"


/////////////////////////////////////////////////////////////////////////////
// CDialParz dialog


CDialParz::CDialParz(CWnd* pParent /*=NULL*/)
	: CDialog(CDialParz::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialParz)
	m_recouvr0 = _T("");
	m_recouvr1 = _T("");
	m_param0 = _T("");
	m_param1 = _T("");
	//}}AFX_DATA_INIT
}
CDialParz::CDialParz(CAfficheMesPol * pView)
{
	m_pView=pView;
}

void CDialParz::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialParz)
	DDX_Text(pDX, IDC_RECOUVR0, m_recouvr0);
	DDX_Text(pDX, IDC_RECOUVR1, m_recouvr1);
	DDX_Text(pDX, IDC_PARAM0, m_param0);
	DDX_Text(pDX, IDC_PARAM1, m_param1);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialParz, CDialog)
	//{{AFX_MSG_MAP(CDialParz)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialParz message handlers
BOOL CDialParz::Create(CString para0,char* ch0,CString para1,char* ch1)
{
	m_param0=para0;
	m_param1=para1;
	m_recouvr0=ch0;
	m_recouvr1=ch1;
	return CDialog::Create(CDialParz::IDD);
}

void CDialParz::OnCancel() 
{
	// TODO: Add extra cleanup here
	
	ShowWindow(SW_HIDE);
}

BOOL CDialParz::OnInitDialog() 
{	
	// TODO: Add extra initialization here
	
	return CDialog::OnInitDialog();
}
void CDialParz::Rafraichissement(CString para0,char* ch0,CString para1,char* ch1)
{
	m_param0=para0;
	m_param1=para1;
	m_recouvr0=ch0;
	m_recouvr1=ch1;
	OnInitDialog();
	
}