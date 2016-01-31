// Dialrelief.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "Dialrelief.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDialrelief dialog


CDialrelief::CDialrelief(CWnd* pParent /*=NULL*/)
	: CDialog(CDialrelief::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialrelief)
	m_dist = 2.5f;
	m_foc = 2.0f;
	m_tet = 45.0f;
	m_xo = 0.0f;
	//}}AFX_DATA_INIT
}


void CDialrelief::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialrelief)
	DDX_Text(pDX, IDC_DISTANCE, m_dist);
	DDX_Text(pDX, IDC_FOCALE, m_foc);
	DDX_Text(pDX, IDC_TETA, m_tet);
	DDX_Text(pDX, IDC_X0, m_xo);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialrelief, CDialog)
	//{{AFX_MSG_MAP(CDialrelief)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialrelief message handlers

void CDialrelief::OnOK() 
{
	// TODO: Add extra validation here
	
	CDialog::OnOK();
}

void CDialrelief::OnCancel() 
{
	// TODO: Add extra cleanup here
	
	CDialog::OnCancel();
}
