// DialConfigPoly.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "DialConfigPoly.h"



/////////////////////////////////////////////////////////////////////////////
// CDialConfigPoly dialog


CDialConfigPoly::CDialConfigPoly(CWnd* pParent /*=NULL*/,CString initFich)
	: CDialog(CDialConfigPoly::IDD, pParent)
{

	NomFich=initFich;
	//{{AFX_DATA_INIT(CDialConfigPoly)
	//}}AFX_DATA_INIT
}


void CDialConfigPoly::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialConfigPoly)
	DDX_Control(pDX, IDC_NOMMMC, m_cmc);
	DDX_Control(pDX, IDOK, m_sauve);
	DDX_Control(pDX, IDC_AJOUTERCLAS, m_ajoutclass);
	DDX_Control(pDX, IDC_AJOUTERFICH, m_ajoutfich);
	DDX_Control(pDX, IDC_SUPPRIMCLAS, m_supprimclass);
	DDX_Control(pDX, IDC_TREE1, m_tree);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialConfigPoly, CDialog)
	//{{AFX_MSG_MAP(CDialConfigPoly)
	ON_NOTIFY(TVN_ITEMEXPANDING, IDC_TREE1, OnItemexpandingTree1)
	ON_WM_MOUSEMOVE()
	ON_WM_LBUTTONUP()
	ON_BN_CLICKED(IDC_NOMMMC, OnNommmc)
	ON_BN_CLICKED(IDC_AJOUTERCLAS, OnAjouterclas)
	ON_BN_CLICKED(IDC_AJOUTERFICH, OnAjouterfich)
	ON_BN_CLICKED(IDC_SUPPRIMCLAS, OnSupprimclas)
	ON_NOTIFY(TVN_SELCHANGED, IDC_TREE1, OnSelchangedTree1)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialConfigPoly message handlers

BOOL CDialConfigPoly::OnInitDialog() 
{
	CDialog::OnInitDialog();
	int i;

	CImageList *imageTree=new CImageList();
	imageTree->Create( IDB_BITMAPTREE,16,5, RGB(255,255,255));
	m_tree.SetImageList(imageTree,TVSIL_NORMAL);

	CStdioFile Fichier; // variable fichier
	CFile prem;
	HTREEITEM root,courant;
	CString text;
	int NBclass,*NBFichParclass;
	if (Fichier.Open(NomFich, CFile::modeRead))
	{ 
		m_cmc.SetWindowText(Fichier.GetFileName());
		root=m_tree.InsertItem(Fichier.GetFileName(),0,1,TVI_ROOT,TVI_FIRST );
		m_tree.AddDataItem(root,Fichier.GetFilePath());

		Fichier.ReadString(text);
		NBclass=atoi(LPCSTR(text));
		NBFichParclass=new int[NBclass];
		for(i=0;i<NBclass;i++)
		{
			Fichier.ReadString(text);
			NBFichParclass[i]=atoi(LPCSTR(text));

		}
		
		for(i=0;i<NBclass;i++)
		{
			m_tree.AjouterClass();
			if(i==0) courant=m_tree.GetNextItem(root,TVGN_CHILD);
			else courant=m_tree.GetNextSiblingItem(courant);
			for(int j=0;j<NBFichParclass[i];j++)
			{
				Fichier.ReadString(text);		
				m_tree.AddDataItem(m_tree.InsertItem(text.Right(text.GetLength()-text.ReverseFind('\\')-1),4,5,courant),text);
				if(!j&&!i) 
				{
					if(prem.Open(text,CFile::modeRead))
					{
						prem.Read(m_tree.FichParam.GetBuffer(1024),1024);
						m_tree.FichParam.ReleaseBuffer();
						prem.Read(&m_tree.Nbpar,1);
						prem.Close();
					}
				}
			}					
		}
		m_tree.Select(root,TVGN_CARET);
		m_tree.Expand(root,TVE_EXPAND);
		Fichier.Close(); 
		m_ajoutclass.EnableWindow(TRUE);
		m_supprimclass.EnableWindow(FALSE);
		m_ajoutfich.EnableWindow(FALSE);
		m_sauve.EnableWindow(TRUE);

	}
	else
	{
		m_ajoutclass.EnableWindow(FALSE);
		m_supprimclass.EnableWindow(FALSE);
		m_ajoutfich.EnableWindow(FALSE);
		m_sauve.EnableWindow(FALSE);
	}
	return TRUE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}

void CDialConfigPoly::OnItemexpandingTree1(NMHDR* pNMHDR, LRESULT* pResult) 
{
	NM_TREEVIEW* pNMTreeView = (NM_TREEVIEW*)pNMHDR;
	int im, selim;
	m_tree.GetItemImage(m_tree.GetSelectedItem(),im, selim);
	if(im==2)
	{
		im=0;
		selim=1;
	}
	else
	{
		im=2;
		selim=3;
	}
	m_tree.SetItemImage(m_tree.GetSelectedItem(),im, selim);
	*pResult = 0;
}

void CDialConfigPoly::OnMouseMove(UINT nFlags, CPoint point) 
{
	if (m_tree.DRAG)
	{
		SetCursor(LoadCursor(0,IDC_NO));
	}
	
	CDialog::OnMouseMove(nFlags, point);
}


void CDialConfigPoly::OnLButtonUp(UINT nFlags, CPoint point) 
{
	if (m_tree.DRAG)
	{
		m_tree.imDrag->EndDrag();
		m_tree.DRAG=FALSE;
		m_tree.RedrawWindow();
	}

	CDialog::OnLButtonUp(nFlags, point);
}

void CDialConfigPoly::OnNommmc() 
{
	CStdioFile Fichier; // variable fichier
	CFile prem;
	HTREEITEM root,courant;
	CString text;
	int NBclass,*NBFichParclass,i;

	CFileDialog * FileData = new CFileDialog (TRUE,"*.cmc",NomFich,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"Configuration  multi-classes (*.cmc)|*.cmc|",this);

	if (FileData->DoModal()==IDOK)
	{
		if (Fichier.Open(FileData->GetPathName(), CFile::modeRead))
		{ 
			NomFich=FileData->GetPathName();
			m_tree.DeleteAllItems();
			m_tree.TabString.RemoveAll();
			m_tree.TabDWord.RemoveAll();
			root=m_tree.InsertItem(FileData->GetFileName(),0,1,TVI_ROOT,TVI_FIRST );
			m_tree.AddDataItem(root,FileData->GetPathName());

			Fichier.ReadString(text);
			NBclass=atoi(LPCSTR(text));
			NBFichParclass=new int[NBclass];
			for(i=0;i<NBclass;i++)
			{
				Fichier.ReadString(text);
				NBFichParclass[i]=atoi(LPCSTR(text));

			}
			
			for(i=0;i<NBclass;i++)
			{
				m_tree.AjouterClass();
				if(i==0) courant=m_tree.GetNextItem(root,TVGN_CHILD);
				else courant=m_tree.GetNextSiblingItem(courant);
				for(int j=0;j<NBFichParclass[i];j++)
				{
					Fichier.ReadString(text);		
					m_tree.AddDataItem(m_tree.InsertItem(text.Right(text.GetLength()-text.ReverseFind('\\')-1),4,5,courant),text);
					if(!j&&!i) 
					{	
						if(prem.Open(text,CFile::modeRead))
						{	
							prem.Read(m_tree.FichParam.GetBuffer(1024),1024);
							m_tree.FichParam.ReleaseBuffer();
							prem.Read(&m_tree.Nbpar,1);
							prem.Close();
						}
					}
				}					
			}
			m_tree.Select(root,TVGN_CARET);
			m_tree.Expand(root,TVE_EXPAND);
			Fichier.Close(); 
		}
		else if (m_tree.GetCount())
		{
			MessageBox("Le fichier  "+FileData->GetFileName()+" sera créé !","Message",MB_OK|MB_ICONEXCLAMATION);
			NomFich=FileData->GetPathName();
			root = m_tree.GetRootItem();
			m_tree.RemoveDataItem(root);
			m_tree.SetItemText(root,FileData->GetFileName());
			m_tree.AddDataItem(root,FileData->GetPathName());
		}
		else
		{
			MessageBox("Le fichier  "+FileData->GetFileName()+" sera créé !","Message",MB_OK|MB_ICONEXCLAMATION);
			NomFich=FileData->GetPathName();
			m_tree.DeleteAllItems();
			m_tree.TabString.RemoveAll();
			m_tree.TabDWord.RemoveAll();
			root=m_tree.InsertItem(FileData->GetFileName(),0,1,TVI_ROOT,TVI_FIRST );
			m_tree.AddDataItem(root,FileData->GetPathName());
			m_tree.AjouterClass();
			m_tree.SelectItem(root);
			m_tree.Expand(root,TVE_EXPAND);
		}
		m_ajoutclass.EnableWindow(TRUE);
		m_cmc.SetWindowText(FileData->GetFileName());
		m_sauve.EnableWindow(TRUE);
	}
	delete FileData;
}

void CDialConfigPoly::OnAjouterclas() 
{
	m_tree.AjouterClass();
}

void CDialConfigPoly::OnAjouterfich() 
{
	m_tree.Ajouterfich();
}

void CDialConfigPoly::OnSupprimclas() 
{
	m_tree.Supprimclas();
}


void CDialConfigPoly::OnSelchangedTree1(NMHDR* pNMHDR, LRESULT* pResult) 
{
	NM_TREEVIEW* pNMTreeView = (NM_TREEVIEW*)pNMHDR;
	// TODO: Add your control notification handler code here
	HTREEITEM selected=pNMTreeView->itemNew.hItem;
	if(selected)
	{
		if(m_tree.GetParentItem(selected)==m_tree.GetRootItem()&&!m_tree.ItemHasChildren(selected)) // sur classe vide
		{
			m_supprimclass.EnableWindow(TRUE);
			m_ajoutfich.EnableWindow(TRUE);
		}
		else if (selected == m_tree.GetRootItem()) // sur root
		{
			m_supprimclass.EnableWindow(FALSE);
			m_ajoutfich.EnableWindow(FALSE);
		}
		else if(m_tree.GetParentItem(selected)==m_tree.GetRootItem()&&m_tree.ItemHasChildren(selected)) // sur classe non vide
		{
			m_supprimclass.EnableWindow(FALSE);
			m_ajoutfich.EnableWindow(TRUE);
		}
		else //sur fichier
		{
			m_supprimclass.EnableWindow(TRUE);
			m_ajoutfich.EnableWindow(FALSE);
		}
	}
	*pResult = 0;
}

void CDialConfigPoly::OnOK() 
{
	HTREEITEM root=m_tree.GetRootItem();
	CStdioFile fichierCMC;
	CFile fichierMMC,**fichierAPPR;
	CString text;
	BYTE NBclass;
	int *NBFichierParClass;
	CString** Nomfichierappr;
	DWORD LongFich;
	int i;

	if(fichierCMC.Open((LPCSTR)m_tree.GetNameItem(root),CFile::modeCreate|CFile::modeWrite))
	{
		
		HTREEITEM courant=m_tree.GetNextItem(root,TVGN_CHILD);
		HTREEITEM FichCourant;

		NBclass=m_tree.GetNBchild(root);
		text.Format("%d",NBclass);
		fichierCMC.WriteString(text);
		fichierCMC.WriteString("\n");

		NBFichierParClass=new int[NBclass];
		Nomfichierappr=new CString*[NBclass];
		fichierAPPR=new CFile*[NBclass];

		for(i=0;i<NBclass;i++)
		{	
			NBFichierParClass[i]=m_tree.GetNBchild(courant);
			if(!NBFichierParClass[i]) 
			{
				MessageBox("Votre fichier contient une classe vide!","Sauvegarde",MB_OK|MB_ICONEXCLAMATION);
				fichierCMC.Close();
				return;
			}
			Nomfichierappr[i]=new CString[NBFichierParClass[i]];
			fichierAPPR[i]=new CFile[NBFichierParClass[i]];
			text.Format("%d",NBFichierParClass[i]);
			fichierCMC.WriteString(text);
			fichierCMC.WriteString("\n");
			courant=m_tree.GetNextSiblingItem(courant);
		}

		courant=m_tree.GetNextItem(root,TVGN_CHILD);
		for(i=0;i<NBclass;i++)
		{	
			FichCourant=m_tree.GetNextItem(courant,TVGN_CHILD);
			for(int j=0;j<NBFichierParClass[i];j++)
			{	
				text=m_tree.GetNameItem(FichCourant);
				Nomfichierappr[i][j]=text;
				fichierCMC.WriteString(text);
				fichierCMC.WriteString("\n");

				FichCourant=m_tree.GetNextSiblingItem(FichCourant);

			}
			courant=m_tree.GetNextSiblingItem(courant);
		}

		fichierCMC.Close();

		if(fichierMMC.Open(fichierCMC.GetFilePath().Left(fichierCMC.GetFilePath().ReverseFind('.'))+".mmc",CFile::modeCreate|CFile::modeWrite))
		{
			///attention DWORD avant WORD...
			DWORD NMES=0;
			long NMestemp;
			BYTE * StTemp;
			int j;
			fichierMMC.Write(&NBclass,1);
			fichierMMC.Write(&m_tree.Nbpar,1);
			for(i=0;i<NBclass;i++)
			{	
				NMES=0;
				for(j=0;j<NBFichierParClass[i];j++)
				{
					if(fichierAPPR[i][j].Open(Nomfichierappr[i][j],CFile::modeRead))
					{
						fichierAPPR[i][j].Seek(1025,CFile::begin); // ne prend pas nomfichparam ni Nbpar
						fichierAPPR[i][j].Read(&NMestemp,sizeof(long));
						NMES+=(WORD)NMestemp;
						fichierAPPR[i][j].Close();
					}
				}
				fichierMMC.Write(&NMES,sizeof(DWORD));

				for(j=0;j<NBFichierParClass[i];j++)
				{
					if(fichierAPPR[i][j].Open(Nomfichierappr[i][j],CFile::modeRead))
					{
						fichierAPPR[i][j].Seek(1025+sizeof(long),CFile::begin); // ne prend pas nomfichparam ni Nbpar ni nbmes
						LongFich = fichierAPPR[i][j].GetLength()-1025-sizeof(long);
						StTemp = new BYTE[LongFich];
						fichierAPPR[i][j].Read(StTemp,LongFich);
						fichierMMC.Write(StTemp,LongFich);
						delete StTemp;
						fichierAPPR[i][j].Close();
					}
				}
			}
			fichierMMC.Close();
		}
	}
	
	CDialog::OnOK();
}
