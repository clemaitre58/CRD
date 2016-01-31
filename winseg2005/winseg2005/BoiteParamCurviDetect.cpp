// BoiteParamCurviDetect.cpp : fichier d'implémentation
//

#include "stdafx.h"
#include "winseg.h"
#include "BoiteParamCurviDetect.h"


// Boîte de dialogue BoiteParamCurviDetect

IMPLEMENT_DYNAMIC(BoiteParamCurviDetect, CDialog)

BoiteParamCurviDetect::BoiteParamCurviDetect(CWnd* pParent /*=NULL*/)
	: CDialog(BoiteParamCurviDetect::IDD, pParent)
	, ValLargeur(80)
	, ValLargeurFourier(4)
	, NbListeMini(20)
	, SeuilReponse(0.008)
	, ValLargeurMini(6)
	, m_EditVariationLocaleLargeur(0.15)
	, m_BDebugEdges(FALSE)
	, m_BDebugTracking(FALSE)
	, m_BDebugCurviReg(FALSE)
{

}

BoiteParamCurviDetect::~BoiteParamCurviDetect()
{
}

void BoiteParamCurviDetect::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_SLIDERLARGEUR, SliderLargeur);
	DDX_Text(pDX, IDC_VALLARGEUR, ValLargeur);
	DDX_Control(pDX, IDC_SLIDERLARGEURFOURIER, SliderLargeurFourier);
	DDX_Text(pDX, IDC_VALLARGEURFOURIER, ValLargeurFourier);
	DDX_Text(pDX, IDC_VALLARGEURNBLIST, NbListeMini);
	DDX_Control(pDX, IDC_SLIDERNBLIST, SliderLongueurForme);
	DDX_Text(pDX, IDC_VALSEUILREPONSE, SeuilReponse);
	DDX_Control(pDX, IDC_SLIDERSEUILREPONSE, SliderSeuilReponse);
	DDX_Control(pDX, IDC_SLIDERLARGEURMINI, SliderLargeurMini);
	DDX_Text(pDX, IDC_VALLARGEURMINI, ValLargeurMini);
	//DDV_MinMaxInt(pDX, ValLargeurMini, 4, 256);
	DDX_Control(pDX, IDC_SLIDERVARIATIONLOCALELARGEUR, m_SliderVariationLocaleLargeur);
	DDX_Text(pDX, IDC_VARIATIONLOCALELARGEUR, m_EditVariationLocaleLargeur);
	DDX_Check(pDX, IDC_CHECKDEBUGCONTOURS, m_BDebugEdges);
	DDX_Check(pDX, IDC_CHECK2, m_BDebugTracking);
	DDX_Check(pDX, IDC_CHECK3, m_BDebugCurviReg);
}


BEGIN_MESSAGE_MAP(BoiteParamCurviDetect, CDialog)
	ON_WM_HSCROLL()
END_MESSAGE_MAP()


// Gestionnaires de messages de BoiteParamCurviDetect


void BoiteParamCurviDetect::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar)
{
	// TODO : ajoutez ici le code de votre gestionnaire de messages et/ou les paramètres par défaut des appels
	UpdateData(TRUE);
	SeuilReponse = ((double)SliderSeuilReponse.GetPos())/2000;
	ValLargeur = (double)SliderLargeur.GetPos();
	ValLargeurFourier = (double)SliderLargeurFourier.GetPos();
	NbListeMini = (double)SliderLongueurForme.GetPos();
	ValLargeurMini = SliderLargeurMini.GetPos();
	m_EditVariationLocaleLargeur = (double)m_SliderVariationLocaleLargeur.GetPos()/100;
	UpdateData(FALSE);

	CDialog::OnHScroll(nSBCode, nPos, pScrollBar);
}

BOOL BoiteParamCurviDetect::OnInitDialog()
{
	CDialog::OnInitDialog();

	CFile FichierSaveParam;

	SliderLargeur.SetRange(16,500);
	SliderLargeurFourier.SetRange(4,256);
	SliderLongueurForme.SetRange(0,1000);
	SliderSeuilReponse.SetRange(0,1000);
	SliderLargeurMini.SetRange(4,245);
	m_SliderVariationLocaleLargeur.SetRange(0,100);


	//relecture dernière valeur sur disque sinon valeur par defaut
	
	if(FichierSaveParam.Open("c:\\SavParam.par",CFile::modeRead))
	{
		UpdateData(TRUE);
		FichierSaveParam.Read(&SeuilReponse,sizeof(double));
		FichierSaveParam.Read(&ValLargeur,sizeof(int));
		FichierSaveParam.Read(&ValLargeurFourier,sizeof(int));
		FichierSaveParam.Read(&NbListeMini,sizeof(int));
		FichierSaveParam.Read(&ValLargeurMini,sizeof(int));
		FichierSaveParam.Read(&m_EditVariationLocaleLargeur,sizeof(double));
		UpdateData(FALSE);
		FichierSaveParam.Close();
	}
	else
	{
		SeuilReponse = 0.008;
		ValLargeur = 80;
		ValLargeurFourier = 4;
		NbListeMini = 20;
		ValLargeurMini = 6;
		m_EditVariationLocaleLargeur = 0.15;
	}

	UpdateData(TRUE);
	SliderLargeur.SetPos(ValLargeur);
	SliderLargeurFourier.SetPos(ValLargeurFourier);
	SliderLongueurForme.SetPos(NbListeMini);
	SliderSeuilReponse.SetPos(SeuilReponse*2000);
	SliderLargeurMini.SetPos(ValLargeurMini);
	m_SliderVariationLocaleLargeur.SetPos(m_EditVariationLocaleLargeur*100);
	m_BDebugEdges=FALSE;
	m_BDebugTracking=FALSE;
	m_BDebugCurviReg=TRUE;
	UpdateData(FALSE);

	

	//placement des curseurs pour ne pas avoir de problèmes de raf avec le onhscroll


	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION : les pages de propriétés OCX devraient retourner FALSE
}

void BoiteParamCurviDetect::OnOK()
{
	// TODO : ajoutez ici votre code spécialisé et/ou l'appel de la classe de base
	CFile FichierSaveParam;
	if(FichierSaveParam.Open("c:\\SavParam.par",CFile::modeCreate|CFile::modeWrite))
	{
		FichierSaveParam.Write(&SeuilReponse,sizeof(double));
		FichierSaveParam.Write(&ValLargeur,sizeof(int));
		FichierSaveParam.Write(&ValLargeurFourier,sizeof(int));
		FichierSaveParam.Write(&NbListeMini,sizeof(int));
		FichierSaveParam.Write(&ValLargeurMini,sizeof(int));
		FichierSaveParam.Write(&m_EditVariationLocaleLargeur,sizeof(double));
		FichierSaveParam.Close();
	}


	CDialog::OnOK();
}
