// lhfloat.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h" 
#include "lhfloat.h"
//#include "graphi.h"                       
/////////////////////////////////////////////////////////////////////////////
// LHfloat dialog

BEGIN_MESSAGE_MAP(LHfloat, CDialog)
	//{{AFX_MSG_MAP(LHfloat)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

LHfloat::LHfloat(DWORD InitTailleFich,CWnd* pParent /*=NULL*/)
	: CDialog(LHfloat::IDD, pParent)
{    
	TailleFich = InitTailleFich;
	//{{AFX_DATA_INIT(LHfloat)
	m_HAUT = 256;
	m_LARG = 256;
	//}}AFX_DATA_INIT
}

void LHfloat::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(LHfloat)
	DDX_Text(pDX, IDC_HAUT, m_HAUT);
	DDX_Text(pDX, IDC_LARG, m_LARG);
	//}}AFX_DATA_MAP
}




/////////////////////////////////////////////////////////////////////////////
// LHfloat message handlers

void LHfloat::OnOK()
{
	// TODO: Add extra validation here
	UpdateData();
	if ((DWORD)m_LARG*m_HAUT*sizeof(float)!=TailleFich)
	MessageBox("Erreur dans les dimensions de l'image","Attention",MB_OK|MB_ICONEXCLAMATION);
	else CDialog::OnOK();
}


/////////////////////////////////////////////////////////////////////////////
// LHByte dialog


LHByte::LHByte(DWORD InitTailleFich,CWnd* pParent /*=NULL*/)
	: CDialog(LHfloat::IDD, pParent)
{    
	TailleFich = InitTailleFich;
	//{{AFX_DATA_INIT(LHfloat)
	m_HAUT = 256;
	m_LARG = 256;
	//}}AFX_DATA_INIT
}

void LHByte::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(LHfloat)
	DDX_Text(pDX, IDC_HAUT, m_HAUT);
	DDX_Text(pDX, IDC_LARG, m_LARG);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(LHByte, CDialog)
	//{{AFX_MSG_MAP(LHfloat)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// LHfloat message handlers

void LHByte::OnOK()
{
	// TODO: Add extra validation here
	UpdateData();
	if ((DWORD)m_LARG*m_HAUT*sizeof(unsigned char)!=TailleFich)
	MessageBox("Erreur dans les dimensions de l'image","Attention",MB_OK|MB_ICONEXCLAMATION);
	else CDialog::OnOK();

}


