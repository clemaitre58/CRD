// DialInfoTrait.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "DialInfoTrait.h"

/////////////////////////////////////////////////////////////////////////////
// CDialInfoTrait dialog


CDialInfoTrait::CDialInfoTrait(CWnd* pParent /*=NULL*/)
	: CDialog(CDialInfoTrait::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialInfoTrait)
	m_Comment = _T("");
	m_Nom = _T("");
	//}}AFX_DATA_INIT
}
void CDialInfoTrait::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialInfoTrait)
	DDX_Text(pDX, IDC_COMMENTAIRES, m_Comment);
	DDX_Text(pDX, IDC_Nom_Traitement, m_Nom);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialInfoTrait, CDialog)
	//{{AFX_MSG_MAP(CDialInfoTrait)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialInfoTrait message handlers
BOOL CDialInfoTrait::Create(CString InitComment, CString InitNom)
{
	// TODO: Add your specialized code here and/or call the base class
	m_Comment = InitComment;
	m_Nom = InitNom;
	return CDialog::Create(CDialInfoTrait::IDD);
}

BOOL CDialInfoTrait::OnInitDialog() 
{
	return CDialog::OnInitDialog();
}

void CDialInfoTrait::UpdateInfo(CString InitComment, CString InitNom)
{
	m_Comment = InitComment;
	m_Nom = InitNom;
	if(GetSafeHwnd()!=0) UpdateData(FALSE);
	
}