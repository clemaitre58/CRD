// Vue3d.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "Vue3d.h"
#include "resource.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CVue3d
const UINT uiMsg4=::RegisterWindowMessage("WM_FINVUE3D");


IMPLEMENT_DYNCREATE(CVue3d, CMDIChildWnd)

CVue3d::CVue3d(unsigned char *PImage, long Larg, long Haut,CWnd * pView)
{
	m_pView=pView;
	PIm=PImage;
	Largeur=Larg;
	Hauteur=Haut;
	PosScroll = 0;
	PosHFen 	= 0;	

}

CVue3d::~CVue3d()
{
}


BEGIN_MESSAGE_MAP(CVue3d, CMDIChildWnd)
	//{{AFX_MSG_MAP(CVue3d)
	ON_REGISTERED_MESSAGE(uiMsg4, OnFinVue)
	ON_WM_PAINT()
	ON_WM_HSCROLL()
	ON_WM_VSCROLL()
	ON_WM_DESTROY()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CVue3d message handlers

BOOL CVue3d::Create(LPCTSTR szTitle) 
{
	// TODO: Add your specialized code here and/or call the base class
		CRect rect(0,0,660,520); 


	menu.LoadMenu(IDR_VUE3D);
	m_hMenuShared = menu.m_hMenu;

	BOOL LeBool=CMDIChildWnd::Create(
	AfxRegisterWndClass(CS_HREDRAW | CS_VREDRAW,NULL,/*(HBRUSH)(COLOR_WINDOW+1)*/(HBRUSH)(COLOR_BACKGROUND),NULL), 
	szTitle, WS_CHILD|WS_VISIBLE|WS_OVERLAPPED|WS_CAPTION|WS_THICKFRAME|WS_MINIMIZEBOX|WS_MAXIMIZEBOX 
/*|WS_HSCROLL|WS_VSCROLL*/, rect, NULL);

	return LeBool;
}


void CVue3d::OnPaint() 
{
	CPaintDC DevC(this); // device context for painting
	// TODO: Add your message handler code here
	CRect rect;
	GetClientRect(&rect);
	DevC.FillSolidRect(0,0,rect.right,rect.bottom,RGB(0,0,0));
	CPen Pen(PS_SOLID,1,RGB(0,0,255));
	DevC.SelectObject(&Pen);

unsigned char yy0,yy;
int xx,xx0,dl,dc;
int ddl=3;
int ddc=3;
int paff=3;
unsigned char *buff;

long p=0;
int l,c;
long decal=0;
buff=(unsigned char *) malloc(Largeur*Hauteur);
for (c=0;c<Largeur;c++)  buff[c]=0;
dl=0;dc=0;
for (l=0;l<Hauteur;l=l+ddl)
  {
    xx0=0;
    yy0=PIm[Largeur*l];
    if ((yy0+dl)<rect.bottom)
    {
    for (c=0;c<Largeur;c++)
     {
       p=Largeur*l+c;
       xx=c;
       yy=PIm[p];
       if (yy>=buff[c+decal] && yy+dl<rect.bottom)
		{
		if (c>0 && yy0>=buff[(c-1)+decal]) 
		{
			//Pen=new CPen(PS_SOLID,1,RGB((BYTE)PIm[p],(BYTE)PIm[p],(BYTE)PIm[p]));
			//DevC.SelectObject(Pen);
			DevC.MoveTo(xx0+dc,rect.bottom-(int(yy0)+dl));
			DevC.LineTo(xx+dc,rect.bottom-(int(yy)+dl));
			//delete Pen;
			
		}
		else  
		{	
			//Pen=new CPen(PS_SOLID,1,RGB((BYTE)PIm[p],(BYTE)PIm[p],(BYTE)PIm[p]));
			//DevC.SelectObject(Pen);
			DevC.MoveTo(xx0+dc,rect.bottom-(int(buff[c+decal])+dl));
			DevC.LineTo(xx+dc,rect.bottom-(int(yy)+dl));
			//delete Pen;
			
		}
		buff[c+decal]=yy;
		}
       xx0=xx;
       yy0=yy;
     }  // for c
    }
    dl=dl+paff;
    dc=dc+ddc;
    //buff=buff+ddc;
	decal=decal+ddc;
    for (c=Largeur-ddc;c<Largeur;c++) buff[c+decal]=0;
    for (c=0;c<Largeur-ddc;c++)
    {
     if (buff[c+decal]>paff) buff[c+decal]=buff[c+decal]-paff;
     else buff[c+decal]=0;
    }
  } // for l
free(buff);	
	// Do not call CMDIChildWnd::OnPaint() for painting messages
}

void CVue3d::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar) 
{
	// TODO: Add your message handler code here and/or call default
	
	CMDIChildWnd::OnHScroll(nSBCode, nPos, pScrollBar);
}

void CVue3d::OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar) 
{
CRect TheRectClient;
	GetClientRect (TheRectClient);     
 if ((TheRectClient.Height()!=Hauteur))    
 { 
 	int  tempPos = PosScroll;
 	switch( nSBCode )
 	{
  		case SB_PAGEDOWN 		: PosScroll = PosScroll +20; 	break;    
  		case SB_LINEDOWN 		: PosScroll ++;  				break;
  		case SB_PAGEUP 			: PosScroll = PosScroll-20;  	break;   
  		case SB_LINEUP 			: PosScroll --;  				break;   
  		case SB_THUMBPOSITION 	: PosScroll = nPos;  			break;
 	}   
 	if (PosScroll>100) PosScroll = 100;    
 	if (PosScroll<0) PosScroll = 0;   
 	PosHFen = PosScroll*((int)Hauteur-TheRectClient.Height())/100;                                         
 	if (tempPos != PosScroll) Invalidate(FALSE);
 }
 else 
 {  
  PosScroll = 0;
 }
 SetScrollPos(SB_VERT,PosScroll,TRUE); 	

	
	CMDIChildWnd::OnVScroll(nSBCode, nPos, pScrollBar);
}

LONG CVue3d::OnFinVue(UINT uiParam,LONG lparam )
{
	if(	(CWnd*)lparam==m_pView) DestroyWindow();
	return 0;
}

void CVue3d::OnDestroy() 
{
	CMDIChildWnd::OnDestroy();
	
		if(GetMDIFrame()->MDIGetActive()==NULL)
	{
		CMenu menutemp;
		menutemp.LoadMenu(IDR_MAINFRAME);
		GetMDIFrame()->MDISetMenu(&menutemp,NULL);
		AfxGetMainWnd()->DrawMenuBar();
	}

}
