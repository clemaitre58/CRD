// DialogAffine.cpp : fichier d'implémentation
//

#include "stdafx.h"
#include "winseg.h"
#include "DialogAffine.h"


// Boîte de dialogue CDialogAffine

IMPLEMENT_DYNAMIC(CDialogAffine, CDialog)

CDialogAffine::CDialogAffine(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogAffine::IDD, pParent)
	, m_M1(0)
	, m_M2(0)
	, m_M3(0)
	, m_M4(0)
	, m_M5(0)
	, m_M6(0)
	, m_M7(0)
{

}

CDialogAffine::~CDialogAffine()
{
}

void CDialogAffine::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDITM1, m_M1);
	DDX_Text(pDX, IDC_EDITM2, m_M2);
	DDX_Text(pDX, IDC_EDITM3, m_M3);
	DDX_Text(pDX, IDC_EDITM4, m_M4);
	DDX_Text(pDX, IDC_EDITM5, m_M5);
	DDX_Text(pDX, IDC_EDITM6, m_M6);
	DDX_Text(pDX, IDC_EDITM7, m_M7);
}


BEGIN_MESSAGE_MAP(CDialogAffine, CDialog)
END_MESSAGE_MAP()


// Gestionnaires de messages de CDialogAffine
