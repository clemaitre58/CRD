// DialProjAppr.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "DialProjAppr.h"



/////////////////////////////////////////////////////////////////////////////
// CDialProjAppr dialog


CDialProjAppr::CDialProjAppr(CTraitImageWnd* pParent,CString initFileProj)
	: CDialog(CDialProjAppr::IDD, (CWnd*)pParent)
{
	//{{AFX_DATA_INIT(CDialProjAppr)
	//}}AFX_DATA_INIT
	m_Parent=pParent;
	FileProj=initFileProj;
	ValidImageAppr=ValidFileAppr=ValidFileParam=FALSE;
}


void CDialProjAppr::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialProjAppr)
	DDX_Control(pDX, IDOK, m_sauvegarde);
	DDX_Control(pDX, IDC_NOMPROJET, m_FichierProj);
	DDX_Control(pDX, IDC_NOMIMAGEAPPR, m_ImageAppr);
	DDX_Control(pDX, IDC_NOMFICHIERPARAM, m_FichierParam);
	DDX_Control(pDX, IDC_NOMFICHIERAPPR, m_FichierAppr);
//	DDX_Control(pDX, IDC_EDITCTRL1, m_editimage);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialProjAppr, CDialog)
	//{{AFX_MSG_MAP(CDialProjAppr)
	ON_BN_CLICKED(IDC_NOMFICHIERAPPR, OnNomfichierappr)
	ON_BN_CLICKED(IDC_NOMFICHIERPARAM, OnNomfichierparam)
	ON_BN_CLICKED(IDC_NOMIMAGEAPPR, OnNomimageappr)
	ON_BN_CLICKED(IDC_NOMPROJET, OnNomprojet)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialProjAppr message handlers
BOOL CDialProjAppr::OnInitDialog() 
{
	CDialog::OnInitDialog();
	
	CStdioFile Fichier;// variable fichier proj
	CFile Test; 	                                        
    CString Nom;	
	SetWindowPos(NULL,0,0,306,360,SWP_NOMOVE);
	m_FichierProj.SetWindowText(FileProj.Right(FileProj.GetLength()-FileProj.ReverseFind('\\')-1));

	if (Fichier.Open(FileProj, CFile::modeRead))
  	{ 	
		
		Fichier.ReadString(ImageAppr);

		Fichier.ReadString(FileAppr);

		Fichier.ReadString(FileParam);
  	 	
		Fichier.Close(); 
		
		if(Test.Open(ImageAppr,CFile::modeRead))
		{
			long large,haute;
			WORD BitCount;
			DWORD DWComp;

			Test.Seek(sizeof(BITMAPFILEHEADER)+sizeof(DWORD),CFile::begin);
			Test.Read(&large,sizeof(LONG)); 
			Test.Read(&haute,sizeof(LONG));
			Test.Seek(sizeof(WORD),CFile::current);
			Test.Read(&BitCount,sizeof(WORD));
			Test.Read(&DWComp,sizeof(DWORD));
			
			if ((DWComp != BI_RGB)||(BitCount!=8 && BitCount!=24))
			{
				MessageBox("Impossible de lire les fichiers BMP compressés","Message",MB_OK|MB_ICONEXCLAMATION);
				Test.Close();
			}
  	 		else if((large!=m_Parent->Larg)||(haute!=m_Parent->Haut))
			{
				MessageBox("Image des zones d'apprentissage de taille incompatible","Erreur",MB_OK|MB_ICONSTOP);
				Test.Close();
			}
			else
			{
				ImageAppr=Test.GetFilePath();
				m_ImageAppr.SetWindowText(Test.GetFileName());
				ValidImageAppr=TRUE;
				Test.Close();
				SetWindowPos(NULL,0,0,590,360,SWP_NOMOVE);
			//	m_editimage.SetImage(ImageAppr);
			//	m_editimage.SetZoom(68);
			//	m_editimage.Display();
			}
		}

		if(Test.Open(FileAppr,CFile::modeRead))
		{
			m_FichierAppr.SetWindowText(Test.GetFileName());
			ValidFileAppr=TRUE;
			Test.Close();
		}

		if(Test.Open(FileParam,CFile::modeRead))
		{
			m_FichierParam.SetWindowText(Test.GetFileName());
			ValidFileParam=TRUE;
			Test.Close();
		}
	}
	m_ImageAppr.EnableWindow(TRUE);
	m_FichierAppr.EnableWindow(TRUE);
	m_FichierParam.EnableWindow(TRUE);
	m_sauvegarde.EnableWindow(ValidImageAppr&ValidFileAppr&ValidFileParam);

	UpdateData();
			
	return TRUE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}

void CDialProjAppr::OnCancel() 
{
	// TODO: Add extra cleanup here
	
	CDialog::OnCancel();
}

void CDialProjAppr::OnOK() 
{
	CStdioFile Fichier; // variable fichier BITMAP	                                        
   

	if(Fichier.Open(FileProj, CFile::modeCreate|CFile::modeWrite))
	{
		Fichier.WriteString(ImageAppr);
  	 	Fichier.WriteString("\n");
		Fichier.WriteString(FileAppr);
  	 	Fichier.WriteString("\n");
		Fichier.WriteString(FileParam);
  	 	Fichier.WriteString("\n");
		Fichier.Close(); 
		CDialog::OnOK();
	}
	else MessageBox("Le fichier  "+FileProj.Right(FileProj.GetLength()-FileProj.ReverseFind('\\')-1)+"  n'est pas créé!","Erreur",MB_OK|MB_ICONSTOP);

}



void CDialProjAppr::OnNomfichierappr() 
{
	CFile Fichier; // variable fichier 	                                        

	CFileDialog * FileData = new CFileDialog (FALSE,"*.appr",FileAppr/*.Left(FileAppr.ReverseFind('\\'))+"\\*.appr"*/,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"Apprentissage (*.appr)|*.appr|",this);

	if (FileData->DoModal()==IDOK)
	{
		if (Fichier.Open(FileData->GetPathName(), CFile::modeRead))
		{ 
			Fichier.Close(); 
		}
		else MessageBox("Le fichier  "+FileData->GetFileName()+"  sera créé lors de la génération de la base .appr","Message",MB_OK|MB_ICONEXCLAMATION);
	
		FileAppr=FileData->GetPathName();
		m_FichierAppr.SetWindowText(FileData->GetFileName());
		ValidFileAppr=TRUE;
	}

	m_sauvegarde.EnableWindow(ValidImageAppr&ValidFileAppr&ValidFileParam);
	UpdateData();
	delete FileData;
}

void CDialProjAppr::OnNomfichierparam() 
{
	CFile Fichier; // variable fichier 	                                        

	CFileDialog * FileData = new CFileDialog (TRUE,"*.param",FileParam/*.Left(FileParam.ReverseFind('\\'))+"\\*.param"*/,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"Paramètres (*.param)|*.param|",this);

	if (FileData->DoModal()==IDOK)
	{
		if (Fichier.Open(FileData->GetPathName(), CFile::modeRead))
		{ 
			Fichier.Close(); 
		}
		else MessageBox("Le fichier  "+FileData->GetFileName()+"  sera créé lors de la définition des paramètres","Message",MB_OK|MB_ICONEXCLAMATION);
	
		FileParam=FileData->GetPathName();
		m_FichierParam.SetWindowText(FileData->GetFileName());
		ValidFileParam=TRUE;
	}

	m_sauvegarde.EnableWindow(ValidImageAppr&ValidFileAppr&ValidFileParam);
	UpdateData();
	delete FileData;
}

void CDialProjAppr::OnNomimageappr() 
{
	CFile Fichier; // variable fichier 	                                        
    CString Nom;
	long large,haute;
	DWORD DWComp;
	WORD BitCount;

	CFileDialog * FileData = new CFileDialog (TRUE,"*.bmp",ImageAppr,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"BITMAP (*.bmp)|*.bmp|BRUTE (*.mem)|*.mem||",this);

	if (FileData->DoModal()==IDOK)
	{
		if (Fichier.Open(FileData->GetPathName(), CFile::modeRead))
		{
			Fichier.Seek(sizeof(BITMAPFILEHEADER)+sizeof(DWORD),CFile::begin);
			Fichier.Read(&large,sizeof(LONG)); 
			Fichier.Read(&haute,sizeof(LONG));
			Fichier.Seek(sizeof(WORD),CFile::current);
			Fichier.Read(&BitCount,sizeof(WORD));
			Fichier.Read(&DWComp,sizeof(DWORD));
			if ((DWComp != BI_RGB)||(BitCount!=8 && BitCount!=24))
			{
				MessageBox("Impossible de lire les fichiers BMP compressés","Message",MB_OK|MB_ICONEXCLAMATION);
				Fichier.Close();
				return;
			}
  	 		Fichier.Close(); 
			if((large!=m_Parent->Larg)||(haute!=m_Parent->Haut))
			{
				MessageBox("Image des zones d'apprentissage de taille incompatible","Erreur",MB_OK|MB_ICONSTOP);
			}
			else
			{
				ImageAppr=FileData->GetPathName();
				m_ImageAppr.SetWindowText(FileData->GetFileName());
				ValidImageAppr=TRUE;
				SetWindowPos(NULL,0,0,590,360,SWP_NOMOVE);
				//m_editimage.SetImage(ImageAppr);
				//m_editimage.SetZoom(68);
				//m_editimage.Display();
			}
		}//if fichier.open
		else MessageBox("Ce fichier n'existe pas","Message",MB_OK|MB_ICONEXCLAMATION);
    }//if filedata

	m_sauvegarde.EnableWindow(ValidImageAppr&ValidFileAppr&ValidFileParam);
	UpdateData();
	delete FileData;
}

void CDialProjAppr::OnNomprojet() 
{
	CStdioFile Fichier; // variables fichiers 	 
	CFile Test;
    CString Nom;
	ValidImageAppr=ValidFileAppr=ValidFileParam=FALSE;

	CFileDialog * FileData = new CFileDialog (TRUE,"*.prj",FileProj,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"PROJET (*.prj)|*.prj|",this);
	if (FileData->DoModal()==IDOK)
	{
		if(FileData->GetFileExt()!="prj") {MessageBox("Ce fichier n'est pas un .prj!","Erreur ",MB_OK|MB_ICONSTOP);return;}
		 if (Fichier.Open(FileData->GetPathName(), CFile::modeRead)) 
		 {		
			Fichier.ReadString(ImageAppr);

			Fichier.ReadString(FileAppr);

			Fichier.ReadString(FileParam);
			
			Fichier.Close(); 
			
			if(Test.Open(ImageAppr,CFile::modeRead))
			{
				long large,haute;
				WORD BitCount;
				DWORD DWComp;

				Test.Seek(sizeof(BITMAPFILEHEADER)+sizeof(DWORD),CFile::begin);
				Test.Read(&large,sizeof(LONG)); 
				Test.Read(&haute,sizeof(LONG));
				Test.Seek(sizeof(WORD),CFile::current);
				Test.Read(&BitCount,sizeof(WORD));
				Test.Read(&DWComp,sizeof(DWORD));
				
				if ((DWComp != BI_RGB)||(BitCount!=8 && BitCount!=24))
				{
					MessageBox("Impossible de lire les fichiers BMP compressés","Message",MB_OK|MB_ICONEXCLAMATION);
					Test.Close();
					
				}
  	 			else if((large!=m_Parent->Larg)||(haute!=m_Parent->Haut))
				{
					MessageBox("Image des zones d'apprentissage de taille incompatible","Erreur",MB_OK|MB_ICONSTOP);
					m_ImageAppr.SetWindowText(Test.GetFileName());
					Test.Close();
					SetWindowPos(NULL,0,0,306,360,SWP_NOMOVE);	
				}
				else
				{
					ImageAppr=Test.GetFilePath();
					m_ImageAppr.SetWindowText(Test.GetFileName());
					ValidImageAppr=TRUE;
					Test.Close();
					SetWindowPos(NULL,0,0,590,360,SWP_NOMOVE);
					//m_editimage.SetImage(ImageAppr);
					//m_editimage.SetZoom(68);
					//m_editimage.Display();
				}
			}
			else 
			{
				m_ImageAppr.SetWindowText("*.bmp");
				SetWindowPos(NULL,0,0,306,360,SWP_NOMOVE);
			}

			if(Test.Open(FileAppr,CFile::modeRead))
			{
				m_FichierAppr.SetWindowText(FileAppr.Right(FileAppr.GetLength()-FileAppr.ReverseFind('\\')-1));
				ValidFileAppr=TRUE;
				Test.Close();
			}
			else m_FichierAppr.SetWindowText("*.appr");

			if(Test.Open(FileParam,CFile::modeRead))
			{
				m_FichierParam.SetWindowText(FileParam.Right(FileParam.GetLength()-FileParam.ReverseFind('\\')-1));
				ValidFileParam=TRUE;
				Test.Close();
			}
			else m_FichierParam.SetWindowText("*.param");

		 }//if fichier.open
		 else MessageBox("Le fichier  "+FileData->GetFileName()+"  n'existe pas, il sera créé","Message",MB_OK|MB_ICONEXCLAMATION);


		FileProj=FileData->GetPathName();
		m_FichierProj.SetWindowText(FileData->GetFileName());
		m_ImageAppr.EnableWindow(TRUE);
		m_FichierAppr.EnableWindow(TRUE);
		m_FichierParam.EnableWindow(TRUE);
		m_sauvegarde.EnableWindow(ValidImageAppr&ValidFileAppr&ValidFileParam);
	}//if filedata
	UpdateData();
	delete FileData;		
}
