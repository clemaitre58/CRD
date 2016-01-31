// TreeConfig.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "TreeConfig.h"


/////////////////////////////////////////////////////////////////////////////
// CTreeConfig

CTreeConfig::CTreeConfig()
{
	DRAG=FALSE;
	Nbpar=0;
	FichParam="";
}

CTreeConfig::~CTreeConfig()
{
}


BEGIN_MESSAGE_MAP(CTreeConfig, CTreeCtrl)
	//{{AFX_MSG_MAP(CTreeConfig)
	ON_WM_MOUSEMOVE()
	ON_NOTIFY_REFLECT(TVN_BEGINDRAG, OnBegindrag)
	ON_WM_LBUTTONUP()
	ON_WM_NCMOUSEMOVE()
	ON_NOTIFY_REFLECT(NM_DBLCLK, OnDblclk)
	ON_WM_RBUTTONDBLCLK()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CTreeConfig message handlers

void CTreeConfig::OnMouseMove(UINT nFlags, CPoint point) 
{
	if (DRAG)
	{
		RedrawWindow();
		imDrag->Draw( GetDC(),0, point+CPoint(-10,-10), ILD_NORMAL);
	}
	
	CTreeCtrl::OnMouseMove(nFlags, point);
}

void CTreeConfig::OnBegindrag(NMHDR* pNMHDR, LRESULT* pResult) 
{
	NM_TREEVIEW* pNMTreeView = (NM_TREEVIEW*)pNMHDR;

	HTREEITEM hSelected=pNMTreeView->itemNew.hItem;

	if(hSelected)
	{
		if((GetParentItem(hSelected)!=GetRootItem())&&(hSelected!=GetRootItem()))
		{
			DRAG = TRUE;
			CPoint pt;
			GetCursorPos(&pt);
			imDrag = CreateDragImage(hSelected);
			pt.Offset(0,0);
			imDrag->BeginDrag(0,pt);
			Moving=hSelected;
		}
	}
	*pResult = 0;
}

void CTreeConfig::OnLButtonUp(UINT nFlags, CPoint point) 
{
	HTREEITEM select,courant;
	CString text;
	UINT FLAG;
	if (DRAG)
	{
		imDrag->EndDrag();
		DRAG=FALSE;
		RedrawWindow();
		FLAG=TVHT_ONITEM;
		select=HitTest(point,&FLAG);
		if(select)
		{	
			if(GetParentItem(select)==GetRootItem())
			{	
				SelectItem(select);
				text=GetItemText(Moving);
				courant = InsertItem(LPCSTR(text),4,5,select);
				SetItemData(courant,GetItemData(Moving));
				DeleteItem(Moving);
				Expand(GetParentItem(courant),TVE_EXPAND);
				SelectItem(courant);
				SelectItem(GetParentItem(courant));
			}
		}
	}
	CTreeCtrl::OnLButtonUp(nFlags, point);
}



void CTreeConfig::OnNcMouseMove(UINT nHitTest, CPoint point) 
{
	RedrawWindow();
	CTreeCtrl::OnNcMouseMove(nHitTest, point);
}

void CTreeConfig::AddDataItem(HTREEITEM myItem,CString initData)
{
	DWORD temp;
	int test = TabString.Add((LPCSTR)initData);
	if (!test)
	{
		TabDWord.Add(0);
		temp=0;
	}
	else
	{
		temp = TabDWord.GetAt(test-1)+1;
		TabDWord.Add(temp);
	}
	SetItemData(myItem,temp);
}

void CTreeConfig::RemoveDataItem(HTREEITEM myItem)
{
	DWORD test;
	int cpt=0;
	test=GetItemData(myItem);
	while (TabDWord.GetAt(cpt)!=test)
		cpt++;
	TabDWord.RemoveAt(cpt);
	TabString.RemoveAt(cpt);

}

void CTreeConfig::OnDblclk(NMHDR* pNMHDR, LRESULT* pResult) 
{
	NM_TREEVIEW* pNMTreeView = (NM_TREEVIEW*)pNMHDR;

	HTREEITEM hSelected=GetSelectedItem();
	//DWORD test;

	if(hSelected)
	{
		if(GetParentItem(hSelected)!=GetRootItem()&&hSelected!=GetRootItem())
		{
		/*int cpt=0;
			test=GetItemData(hSelected);
			while (TabDWord.GetAt(cpt)!=test)
				cpt++;*/
			MessageBox(GetNameItem(hSelected));
		}
	}
	*pResult = 0;
}

void CTreeConfig::OnRButtonDblClk(UINT nFlags, CPoint point) 
{
	HTREEITEM select,prev;
	UINT FLAG;
	FLAG=TVHT_ONITEM;
	select=HitTest(point,&FLAG);
	if(select)
	{
		if((GetParentItem(select)!=GetRootItem())&&(select!=GetRootItem()))
		{
			SelectItem(select);
			prev=GetPrevVisibleItem(select);
			RemoveDataItem(select);
			DeleteItem(select);
			CTreeCtrl::OnRButtonDblClk(nFlags, point);
			SelectItem(prev);
		}
	}
	

}

void CTreeConfig::TriClasse()
{
	HTREEITEM Root=GetRootItem(), courant;
	if(ItemHasChildren(Root))
	{
		int i=0;
		CString text;
		courant=GetNextItem(Root,TVGN_CHILD);
		while(courant!=NULL)
		{
			text.Format("Classe n°%d",i++);
			SetItemText(courant,LPCSTR(text));
			courant=GetNextSiblingItem(courant);
		}
	}
}

void CTreeConfig::AjouterClass()
{
	HTREEITEM Root=GetRootItem(), courant;
	if(ItemHasChildren(Root))
	{
		int i=0;
		CString text;
		courant=GetNextItem(Root,TVGN_CHILD);
		while(courant!=NULL)
		{
			i++;
			courant=GetNextSiblingItem(courant);
		}
		text.Format("Classe n°%d",i);
		InsertItem(text,0,1,Root);
	}
	else InsertItem("Classe n°0",0,1,Root);
}


void CTreeConfig::Ajouterfich()
{
	CFile Fichier; // variable fichier
	HTREEITEM courant ;
	CString testfichparam;

	CFileDialog * FileData = new CFileDialog (TRUE,"*.appr",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"Apprentissage (*.appr)|*.appr|",this);

	if (FileData->DoModal()==IDOK)
	{
		if(FileData->GetFileExt()!="appr")
		{
			MessageBox("Ce Fichier n'est pas conforme!","Erreur",MB_OK|MB_ICONSTOP);
		}
		else
		if (Fichier.Open(FileData->GetPathName(), CFile::modeRead))
		{ 
			Fichier.Read(testfichparam.GetBuffer(1024),1024);
			
			if(FichParam=="") {FichParam=testfichparam;Fichier.Read(&Nbpar,1);}
			else if(testfichparam!=FichParam) 
			{
				testfichparam.ReleaseBuffer();
				MessageBox("Ce fichier n'est pas compatible!\n\nNom du fichier de paramètres différent!","Erreur",MB_OK|MB_ICONSTOP);
				Fichier.Close();
				delete FileData;
				return;
			}
			Fichier.Close(); 
			courant=InsertItem(FileData->GetFileName(),4,5,GetSelectedItem());
			AddDataItem(courant,FileData->GetPathName());
			Expand(GetParentItem(courant),TVE_EXPAND);
			SelectItem(courant);
			SelectItem(GetParentItem(courant));
		}
		else MessageBox("Le fichier  "+FileData->GetFileName()+" n'existe pas!","Message",MB_OK|MB_ICONEXCLAMATION);
	}
	delete FileData;
}


void CTreeConfig::Supprimclas()
{
	HTREEITEM courant = GetSelectedItem(),prev=GetPrevVisibleItem(courant);
	if(GetParentItem(courant)==GetRootItem())
	{
		DeleteItem(courant);
		TriClasse();
		SelectItem(prev);
	}
	else
	{
		RemoveDataItem(courant);
		DeleteItem(courant);
		SelectItem(prev);
	}

}



int CTreeConfig::GetNBchild(HTREEITEM inititem)
{
	HTREEITEM courant=GetNextItem(inititem,TVGN_CHILD);
	int nb=0;
	while(courant!=NULL)
	{
		nb++;
		courant=GetNextSiblingItem(courant);
	}
	return nb;
}

CString CTreeConfig::GetNameItem(HTREEITEM initem)
{
	DWORD test;
	int cpt=0;
	test=GetItemData(initem);
	while (TabDWord.GetAt(cpt)!=test)
		cpt++;
	return TabString.GetAt(cpt);
}
