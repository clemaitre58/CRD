// ProfilWnd.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "ProfilWnd.h"
#include "message.h"

/////////////////////////////////////////////////////////////////////////////
// CProfilWnd
const UINT uiMsg3=::RegisterWindowMessage("WM_FINPROFIL");

IMPLEMENT_DYNCREATE(CProfilWnd, CGraphWnd)

CProfilWnd::CProfilWnd()
{
}

CProfilWnd::CProfilWnd(int  * InitTab, int Nb, char * InitChaineX,char * InitChaineY, CWnd * InitFeneAppel,CPoint Ori,CPoint Fin, BOOL InitGraphProfil)
					:CGraphWnd( InitTab,  Nb,  InitChaineX, InitChaineY,NULL,InitGraphProfil)
{
	PointOri = Ori;
	PointFin = Fin;      
	FeneAppel =  InitFeneAppel;
	//FeneAppel->SendMessage(WM_COUPE,NULL,(LONG) (HWND *) this );
    MIN_MAX_AUTO = TRUE;

}

CProfilWnd::CProfilWnd(long  * InitTab,int Nb, char * InitChaineX,char * InitChaineY, CWnd * InitFeneAppel,CPoint Ori,CPoint Fin, BOOL InitGraphProfil)
					:CGraphWnd((long *) InitTab,  Nb,  InitChaineX, InitChaineY,NULL,InitGraphProfil)
{
	PointOri = Ori;
	PointFin = Fin;      
	FeneAppel =  InitFeneAppel;
	//FeneAppel->SendMessage(WM_COUPE,NULL,(LONG) (HWND *) this );
    MIN_MAX_AUTO = TRUE;

}


CProfilWnd::CProfilWnd(float  * InitTab,int Nb, char * InitChaineX,char * InitChaineY, CWnd * InitFeneAppel,CPoint Ori,CPoint Fin, BOOL InitGraphProfil)
					:CGraphWnd((float *) InitTab,  Nb,  InitChaineX, InitChaineY,NULL,InitGraphProfil)
{
	PointOri = Ori;
	PointFin = Fin;      
	FeneAppel =  InitFeneAppel;
	//FeneAppel->SendMessage(WM_COUPE,NULL,(LONG) (HWND *) this );
    MIN_MAX_AUTO = TRUE;

}

CProfilWnd::~CProfilWnd()
{
}


BEGIN_MESSAGE_MAP(CProfilWnd, CGraphWnd)
	//{{AFX_MSG_MAP(CProfilWnd)
	ON_WM_DESTROY()
	ON_REGISTERED_MESSAGE(uiMsg3, OnFinProfil)
	ON_WM_CHILDACTIVATE()
	ON_COMMAND(ID_MIN_MAX_AUTO, OnMinMaxNonAuto)
	ON_UPDATE_COMMAND_UI(ID_MIN_MAX_AUTO, OnUpdateMinMaxAuto)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CProfilWnd message handlers

void CProfilWnd::OnDestroy() 
{
	CGraphWnd::OnDestroy();
	
	FeneAppel->SendMessage(WM_COUPE,NULL,(LONG) (HWND *) NULL );
	if (GetMDIFrame()->MDIGetActive() == NULL)
	{
		CMenu toto;
		toto.LoadMenu(IDR_MAINFRAME);
		GetMDIFrame()->MDISetMenu(&toto,NULL);
		AfxGetMainWnd()->DrawMenuBar();
	}
	
}

afx_msg void CProfilWnd::OnChildActivate()
{                                               
	FeneAppel->SendMessage(WM_COUPE,NULL,(LONG) (HWND *) this );	
	CMDIChildWnd::OnChildActivate();
}

void CProfilWnd :: AfficheDimCoupe()
{
	char dim[100];
	sprintf(dim,"Coupe: (%d,%d) - (%d,%d)", PointFin.x, PointFin.y, PointOri.x, PointOri.y);
	AfficheInfo2(dim);
}


void CProfilWnd ::SetPositionCoupe(long i1,long j1,long i2,long j2)
{
	PointOri.x = (int)i1;
	PointOri.y = (int)j1;
	PointFin.x = (int)i2;
	PointFin.y = (int)j2; 
}

void CProfilWnd ::GetInfo(long * i1,long * j1,long * i2,long * j2,long * N)
{
	*i1 = PointOri.x;
	*j1 = PointOri.y; 
	*i2 = PointFin.x;
	*j2 = PointFin.y;
	*N	= NbPoint;		
}

void CProfilWnd::OnMinMaxNonAuto()
{
	if (MIN_MAX_AUTO) 
	{
		MIN_MAX_AUTO = FALSE; 
	}	
	else 
	{
		MIN_MAX_AUTO = TRUE;
	}
	Invalidate();
		
}

LONG CProfilWnd::OnFinProfil(UINT uiParam,LONG lparam )
{
	if(	(CWnd*)lparam==FeneAppel) DestroyWindow();
	return 0;
}

void CProfilWnd::OnUpdateMinMaxAuto(CCmdUI* pCmdUI) 
{
	pCmdUI->SetCheck(MIN_MAX_AUTO);	
}
