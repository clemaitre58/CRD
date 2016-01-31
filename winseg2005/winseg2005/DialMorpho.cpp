// DialMorpho.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "DialMorpho.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDialMorpho dialog


CDialMorpho::CDialMorpho(CWnd* pParent /*=NULL*/)
	: CDialog(CDialMorpho::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialMorpho)
	m_taillestructurant = 1;
	//}}AFX_DATA_INIT
}


void CDialMorpho::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialMorpho)
	DDX_Text(pDX, IDC_TAILLE, m_taillestructurant);
	DDV_MinMaxInt(pDX, m_taillestructurant, 1, 99);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialMorpho, CDialog)
	//{{AFX_MSG_MAP(CDialMorpho)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialMorpho message handlers

void CDialMorpho::OnOK() 
{
	UpdateData(true);
//	if (!(m_taillestructurant%2))
//	{
//		MessageBox("Nombre Pair\nTrès Bien\nEssaie encore !","CDialMorpho");
//		return;
//	}
	CDialog::OnOK();
}

void CDialMorpho::OnCancel() 
{
	// TODO: Add extra cleanup here
	
	CDialog::OnCancel();
}
