
// imagecla.cpp : implementation file
//
  
#include "stdafx.h"
#include "winseg.h"
#include "imagecla.h"
#include "message.h"
#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif



const UINT uiMsg=::RegisterWindowMessage("WM_FINCLASSE");


/////////////////////////////////////////////////////////////////////////////
// CImageClas

IMPLEMENT_DYNAMIC(CImageClas, CTraitImageWnd)
                   
  
CImageClas::CImageClas(HWND pere,CString S, long InitLarg, long InitHaut, TYPEIMAGE InitTypeImage)
					:CTraitImageWnd(S,InitLarg,InitHaut,InitTypeImage)
{
	HFenetre=pere;
}    

CImageClas::CImageClas(HWND pere,LONG InitLarg,LONG InitHaut,  TYPEIMAGE InitLeType)
					:CTraitImageWnd( InitLarg, InitHaut,  InitLeType)
{
	HFenetre=pere;
}
  
  

CImageClas::~CImageClas()
{
	::SendMessage(HFenetre,uiMsg,0,(LPARAM)this);
}


BEGIN_MESSAGE_MAP(CImageClas, CTraitImageWnd)
	//{{AFX_MSG_MAP(CImageClas)
	ON_WM_DESTROY()
	ON_WM_CLOSE()
	ON_WM_CREATE()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CImageClas message handlers

void CImageClas::OnDestroy()
{
	CTraitImageWnd::OnDestroy();
	
	::SendMessage(HFenetre,uiMsg,0,(LPARAM)this);

	// TODO: Add your message handler code here

}


void CImageClas::OnClose() 
{
	// TODO: Add your message handler code here and/or call default
	::SendMessage(HFenetre,uiMsg,0,(LPARAM)this);

	CTraitImageWnd::OnClose();
}

int CImageClas::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CTraitImageWnd::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	// TODO: Add your specialized creation code here
	menu.EnableMenuItem(6,MF_GRAYED|MF_BYPOSITION);
	menu.EnableMenuItem(7,MF_GRAYED|MF_BYPOSITION);

	return 0;
}
