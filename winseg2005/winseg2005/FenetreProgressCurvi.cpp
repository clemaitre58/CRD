// FenetreProgressCurvi.cpp : fichier d'implémentation
//

#include "stdafx.h"
#include "winseg.h"
#include "FenetreProgressCurvi.h"


// Boîte de dialogue CFenetreProgressCurvi

IMPLEMENT_DYNAMIC(CFenetreProgressCurvi, CDialog)

CFenetreProgressCurvi::CFenetreProgressCurvi(CWnd* pParent /*=NULL*/)
	: CDialog(CFenetreProgressCurvi::IDD, pParent)
	, m_TPetrou(0)
	, m_TCTracking(0)
	, m_T1DProcess(0)
	, m_T2DProcess(0)
	, m_Affichage(0)
	, m_TexteValProgress(_T(""))
	, m_TxtNbReg(_T(""))
{

}

CFenetreProgressCurvi::~CFenetreProgressCurvi()
{
}

void CFenetreProgressCurvi::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_STATIC_Petrou, m_TPetrou);
	DDX_Text(pDX, IDC_STATIC_1DProcess, m_TCTracking);
	DDX_Text(pDX, IDC_STATIC_Display, m_T1DProcess);
	DDX_Text(pDX, IDC_STATIC_2DProcess, m_T2DProcess);
	DDX_Text(pDX, IDC_STATIC_Tracking, m_Affichage);
	DDX_Text(pDX, IDC_STATIC_POURCENTAGE, m_TexteValProgress);
	DDX_Control(pDX, IDC_PROGRESSCURVI, m_ControlBarProgressCurvi);
	DDX_Text(pDX, IDC_NBREG, m_TxtNbReg);
}


BEGIN_MESSAGE_MAP(CFenetreProgressCurvi, CDialog)
	ON_WM_DESTROY()
END_MESSAGE_MAP()


// Gestionnaires de messages de CFenetreProgressCurvi

void CFenetreProgressCurvi::OnDestroy()
{
	CDialog::OnDestroy();

	//delete this;
}

void CFenetreProgressCurvi::PostNcDestroy()
{
	// TODO : ajoutez ici votre code spécialisé et/ou l'appel de la classe de base

	CDialog::PostNcDestroy();
	delete this;
}
