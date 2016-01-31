// dialeche.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h" 
#include "dialeche.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDialEchelle dialog


CDialEchelle::CDialEchelle(int InitEchelle,CWnd* pParent)
	: CDialog(CDialEchelle::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialEchelle)
	m_Echelle = InitEchelle;
	//}}AFX_DATA_INIT
}

void CDialEchelle::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialEchelle)
	DDX_Text(pDX, IDC_EDIT1, m_Echelle);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CDialEchelle, CDialog)
	//{{AFX_MSG_MAP(CDialEchelle)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CDialEchelle message handlers

void CDialEchelle::OnOK()
{
	
	UpdateData();
	if (m_Echelle>=1) CDialog::OnOK();
}
