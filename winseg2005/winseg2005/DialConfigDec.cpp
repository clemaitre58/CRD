// DialConfigDec.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "DialConfigDec.h"

/////////////////////////////////////////////////////////////////////////////
// CDialConfigDec dialog


CDialConfigDec::CDialConfigDec(CWnd* pParent /*=NULL*/,CString initpmc,CString initParam,BOOL initcheck,int initclasse)
	: CDialog(CDialConfigDec::IDD, pParent)
{
	NomFichPMC=initpmc;
	NomFichParam=initParam;
	ValidPMC=ValidPARAM=FALSE;
	//{{AFX_DATA_INIT(CDialConfigDec)
	m_check = initcheck;
	m_classeaffect = initclasse;
	//}}AFX_DATA_INIT
}


void CDialConfigDec::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialConfigDec)
	DDX_Control(pDX, IDC_CHECK1, m_checkbox);
	DDX_Control(pDX, IDC_EDIT1, m_edit1);
	DDX_Control(pDX, IDOK, m_ok);
	DDX_Control(pDX, IDC_SPIN1, m_spin);
	DDX_Control(pDX, IDC_BUTTON2, m_param);
	DDX_Control(pDX, IDC_BUTTON1, m_pmc);
	DDX_Check(pDX, IDC_CHECK1, m_check);
	DDX_Text(pDX, IDC_EDIT1, m_classeaffect);
	DDV_MinMaxInt(pDX, m_classeaffect, 0, 255);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialConfigDec, CDialog)
	//{{AFX_MSG_MAP(CDialConfigDec)
	ON_BN_CLICKED(IDC_BUTTON1, OnButton1)
	ON_BN_CLICKED(IDC_BUTTON2, OnButton2)
	ON_BN_CLICKED(IDC_CHECK1, OnCheck1)
	ON_WM_VSCROLL()
	ON_EN_CHANGE(IDC_EDIT1, OnChangeEdit1)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialConfigDec message handlers
BOOL CDialConfigDec::OnInitDialog() 
{
	CDialog::OnInitDialog();
	

	CFile Test;
	BYTE compnbpar;

	if(Test.Open(NomFichPMC,CFile::modeRead))
	{
		Test.Read(&NBclass,1);
		Test.Read(&NBpar,1);
		m_pmc.SetWindowText(Test.GetFileName());
		ValidPMC=TRUE;
		Test.Close();
	
		if(Test.Open(NomFichParam,CFile::modeRead))
		{
			Test.Read(&compnbpar,1);
			if(compnbpar!=NBpar)
			{
				MessageBox("Fichier  "+Test.GetFileName()+"  incompatible!","Nombre de paramètres",MB_OK|MB_ICONSTOP);
				Test.Close();
			}
			else
			{
				m_param.SetWindowText(Test.GetFileName());
				ValidPARAM=TRUE;
				Test.Close();
			}
		}
	}
	m_ok.EnableWindow(ValidPMC&&ValidPARAM);
	m_param.EnableWindow(ValidPMC);
	m_edit1.EnableWindow(m_check);
	m_spin.EnableWindow(m_check);
	m_checkbox.EnableWindow(ValidPMC&&ValidPARAM);
	m_spin.SetRange(0,NBclass-1);
	m_spin.SetPos(m_classeaffect);
	return TRUE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}

void CDialConfigDec::OnButton1() 
{
	CFile Fichier; // variable fichier 	                                        
	BYTE compnbpar;
	CFileDialog * FileData = new CFileDialog (TRUE,"*.pmc",NomFichPMC,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"Polytopes multi-classes (*.pmc)|*.pmc|",this);

	if (FileData->DoModal()==IDOK)
	{
		if(FileData->GetFileExt()!="pmc")
		{
			MessageBox("Le fichier selectionné n'est pas valide!","Erreur",MB_OK|MB_ICONSTOP);
			m_pmc.SetWindowText("*.pmc");
			ValidPMC=FALSE;
		}
		else if (Fichier.Open(FileData->GetPathName(), CFile::modeRead))
		{ 
			Fichier.Read(&NBclass,1);
			Fichier.Read(&compnbpar,1);
			Fichier.Close(); 
			NomFichPMC=FileData->GetPathName();
			m_pmc.SetWindowText(FileData->GetFileName());
			ValidPMC=TRUE;
			if(compnbpar!=NBpar)
			{
				NBpar=compnbpar;
				m_param.SetWindowText("*.param");
				ValidPARAM=FALSE;
			}
			m_classeaffect=0;
			m_spin.SetPos(0);
			m_spin.SetRange(0,NBclass-1);
			if(!m_check) m_edit1.SetWindowText("");
		}
		else 
		{
			MessageBox("Le fichier  "+FileData->GetFileName()+"  n'existe pas!","Erreur",MB_OK|MB_ICONSTOP);
			ValidPMC=FALSE;
		}
	}

	m_ok.EnableWindow(ValidPARAM&&ValidPMC);
	m_checkbox.EnableWindow(ValidPARAM&&ValidPMC);
	m_param.EnableWindow(ValidPMC);
	UpdateData(FALSE);
	delete FileData;	
}

void CDialConfigDec::OnButton2() 
{
	CFile Fichier; // variable fichier 	                                        
	BYTE compnbpar;
	CFileDialog * FileData = new CFileDialog (TRUE,"*.param",NomFichParam,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"Paramètres (*.param)|*.param|",this);

	if (FileData->DoModal()==IDOK)
	{
		if(FileData->GetFileExt()!="param")
		{
			MessageBox("Le fichier selectionné n'est pas valide!","Erreur",MB_OK|MB_ICONSTOP);
			m_pmc.SetWindowText("*.param");
			ValidPARAM=FALSE;
		}
		else if (Fichier.Open(FileData->GetPathName(), CFile::modeRead))
		{ 
			Fichier.Read(&compnbpar,1);

			if(compnbpar!=NBpar)
			{
				MessageBox("Fichier  "+Fichier.GetFileName()+"  incompatible!","Nombre de paramètres",MB_OK|MB_ICONSTOP);
				m_param.SetWindowText("*.param");
				ValidPARAM=FALSE;
				Fichier.Close(); 
			}
			else
			{
				NomFichParam=FileData->GetPathName();
				m_param.SetWindowText(FileData->GetFileName());
				ValidPARAM=TRUE;
				Fichier.Close();
			}
		}
		else 
		{
			MessageBox("Le fichier  "+FileData->GetFileName()+"  n'existe pas!","Erreur",MB_OK|MB_ICONSTOP);
			ValidPARAM=FALSE;
		}
	}

	m_ok.EnableWindow(ValidPARAM&&ValidPMC);
	m_checkbox.EnableWindow(ValidPARAM&&ValidPMC);
	m_param.EnableWindow(ValidPMC);
	UpdateData(FALSE);
	delete FileData;
}

void CDialConfigDec::OnCheck1() 
{
	UpdateData();
	m_edit1.EnableWindow(m_check);
	m_spin.EnableWindow(m_check);
}


void CDialConfigDec::OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar) 
{
	m_classeaffect=nPos;
	UpdateData(FALSE);	
	CDialog::OnVScroll(nSBCode, nPos, pScrollBar);
}


void CDialConfigDec::OnChangeEdit1() 
{
	if(m_check) 
	{
		UpdateData();	
		if(m_classeaffect>(NBclass-1)) 
		{
			CString text;
			text.Format("Valeur comprise entre 0 et %d ",NBclass-1);
			MessageBox(text,NULL,MB_OK|MB_ICONEXCLAMATION);
		}
	}
}

void CDialConfigDec::OnOK() 
{
	// TODO: Add extra validation here
	UpdateData();
	if(m_classeaffect>(NBclass-1)) 
	{
		CString text;
		text.Format("Valeur comprise entre 0 et %d ",NBclass-1);
		MessageBox(text,NULL,MB_OK|MB_ICONEXCLAMATION);
	}else CDialog::OnOK();
}
