// larghaut.cpp : implementation file
//

#include "stdafx.h"
#include "image.h"
#include "larghaut.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CLargHaut dialog


CLargHaut::CLargHaut(CWnd* pParent /*=NULL*/)
	: CDialog(CLargHaut::IDD, pParent)
{
	//{{AFX_DATA_INIT(CLargHaut)
	m_TITRE = "sans titre";
	m_HAUT = 256;
	m_LARG = 256;
	COLOR = 0;
	//}}AFX_DATA_INIT
	COLOR = FALSE;
}

void CLargHaut::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CLargHaut)
	DDX_Text(pDX, IDC_TITRE, m_TITRE);
	DDV_MaxChars(pDX, m_TITRE, 20);
	DDX_Text(pDX, IDC_HAUT, m_HAUT);
	DDX_Text(pDX, IDC_LARG, m_LARG);
	DDX_Radio(pDX, IDC_NB, COLOR);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CLargHaut, CDialog)
	//{{AFX_MSG_MAP(CLargHaut)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CLargHaut message handlers

void CLargHaut::OnOK()
{
	// TODO: Add extra validation here
	BOOL res = UpdateData();  
	if ((div(m_LARG,4).rem!=0) || (m_LARG<1) || (m_LARG>1024))
	{   
		MessageBox("La largeur doit être :\n* multiple de 4\n* < ou = à 1024\n* > à 0","Attention",MB_OK|MB_ICONEXCLAMATION);

	}
	else 
	if ((m_HAUT<1) || (m_HAUT>1024))
	{   
		MessageBox("La hauteur doit être :\n* < ou = à 1024\n* supérieure à 0","Attention",MB_OK|MB_ICONEXCLAMATION);

	}
	else 
		CDialog::OnOK();
}
