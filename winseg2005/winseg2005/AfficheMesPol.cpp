// AfficheMesPol.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "AfficheMesPol.h"
#include "resource.h"
#include "wingdi.h"
#include "DialParz.h"
#include "Graph2Wnd.h"
#include "DialChoixClass.h"
#define _CRT_SECURE_NO_DEPRECATE
#pragma warning( disable : 4996 )

const UINT uiMsg2=::RegisterWindowMessage("WM_FINVISUMES");
/////////////////////////////////////////////////////////////////////////////
// CAfficheMesPol

//CMenu NEAR CAfficheMesPol::menu;     

IMPLEMENT_DYNCREATE(CAfficheMesPol, CMDIChildWnd)

CAfficheMesPol::CAfficheMesPol()
{
	param0=0;
	param1=1;
	Hparam=0;
	Visu=FALSE;
	DialParz=new CDialParz(this);
	pParzen0=NULL;
	pHisto0=NULL;
	pHisto1=NULL;
	pParzen1=NULL;
	Parzen=FALSE;
	Graph=FALSE;
	class1=0;
	class2=1;
	probc1=1.0;
	probc2=1.0;
	ValidCalcul=FALSE;
}

CAfficheMesPol::~CAfficheMesPol()
{
	delete DialParz;
	::SendMessage(m_pView->m_hWnd,uiMsg2,0,0);
}
	


BEGIN_MESSAGE_MAP(CAfficheMesPol, CMDIChildWnd)
	//{{AFX_MSG_MAP(CAfficheMesPol)
	ON_COMMAND(ID_FICHIER_QUITTER, OnFichierQuitter)
	ON_COMMAND(ID_VISUALISATION_CHOIXDESPARAMETRES, OnVisualisationChoixdesparametres)
	ON_COMMAND(ID_VISUALISATION_VISUALISERLESPOLYTOPES, OnVisualisationVisualiserlespolytopes)
	ON_WM_PAINT()
	ON_COMMAND(ID_VISUALISATION_INFOMESURES, OnVisualisationInfomesures)
	ON_COMMAND(ID_VISUALISATION_INFOPOLYTOPES, OnVisualisationInfopolytopes)
	ON_COMMAND(ID_PARZEN, OnParzen)
	ON_COMMAND(ID_WINDOW_CLOSEGRAPHICS, OnWindowClosegraphics)
	ON_UPDATE_COMMAND_UI(ID_WINDOW_CLOSEGRAPHICS, OnUpdateWindowClosegraphics)
	ON_COMMAND(ID_PARZEN_CALCULDETOUS, OnParzenCalculdetous)
	ON_UPDATE_COMMAND_UI(ID_VISUALISATION_VISUALISERLESPOLYTOPES, OnUpdateVisualisationVisualiserlespolytopes)
	ON_COMMAND(ID_CORRELATION_CALCUL, OnCorrelationCalcul)
	ON_UPDATE_COMMAND_UI(ID_CORRELATION_CALCUL, OnUpdateCorrelationCalcul)
	ON_COMMAND(ID_CORRELATION_CALCULDANS, OnCorrelationCalculdans)
	ON_UPDATE_COMMAND_UI(ID_CORRELATION_CALCULDANS, OnUpdateCorrelationCalculdans)
	ON_COMMAND(ID_VISUALISATION_TOOLBAR, OnVisualisationToolbar)
	ON_UPDATE_COMMAND_UI(ID_VISUALISATION_TOOLBAR, OnUpdateVisualisationToolbar)
	ON_WM_DESTROY()
	ON_COMMAND(ID_VISUALISATION_CHOIXDUNFICHIERDEPOLYTOPES, OnVisualisationChoixdunfichierdepolytopes)
	ON_UPDATE_COMMAND_UI(ID_VISUALISATION_INFOPOLYTOPES, OnUpdateVisualisationInfopolytopes)
	ON_WM_SIZE()
	ON_COMMAND(ID_CALCULS_CHOIXDECLASSES, OnCalculsChoixdeclasses)
	ON_WM_CREATE()
	ON_UPDATE_COMMAND_UI(ID_PARZEN, OnUpdateParzen)
	ON_UPDATE_COMMAND_UI(ID_PARZEN_CALCULDETOUS, OnUpdateParzenCalculdetous)
	ON_WM_MDIACTIVATE()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CAfficheMesPol message handlers
BOOL CAfficheMesPol::Create(LPCSTR szTitle)
{
  //dimension de la fenetre
  
	CRect rect(0,0,300,330); 


	menu.LoadMenu(IDR_MESPOL);
	m_hMenuShared = menu.m_hMenu;

	BOOL LeBool=CMDIChildWnd::Create(
	AfxRegisterWndClass( CS_DBLCLKS |CS_HREDRAW | CS_VREDRAW | CS_BYTEALIGNWINDOW | CS_SAVEBITS,NULL,/*(HBRUSH)(COLOR_WINDOW+1)*/(HBRUSH)(COLOR_BACKGROUND),NULL), 
	szTitle, WS_CHILD|WS_VISIBLE|WS_OVERLAPPEDWINDOW, rect, NULL);
	
	//ToolBar
	CToolBarCtrl& pToolBarCtrl=(&((CMainFrame*) GetMDIFrame())->bar)->GetToolBarCtrl();
	pToolBarCtrl.HideButton(ID_VISUALISATION_CHOIXDESPARAMETRES,FALSE);
	pToolBarCtrl.HideButton(ID_VISUALISATION_VISUALISERLESPOLYTOPES,FALSE);
	pToolBarCtrl.HideButton(ID_PARZEN,FALSE);
//	pToolBarCtrl.UpdateWindow( );

//	((CMainFrame*)GetMDIFrame())->RecalcLayout();
	return LeBool;
}
 
int CAfficheMesPol::Init(CWnd* mere,CString NomImil,CString NomFich,CDialPara *Dialp)
{
	m_pView=mere;
	DialPara=Dialp;
	NomFichierPol=NomFich;
	ValidPol=FALSE;
//	BYTE NbPar;

	if (!PolyMc.Init(NomImil))  
	{
	  MessageBox("erreur init pol");
	 // Poly.FinAlloc(); 
	  return 0;
	}
	else 
	 {
	/*	NbPar=PolyMc.NPAR;
		ValidPol=PolyMc.LirePol(NomFichierPol);	
		if(PolyMc.NPAR!=NbPar) 
		{
			ValidPol=FALSE;
			PolyMc.NPAR=NbPar;
			MessageBox("Fichier de polytopes incompatible!","Erreur!",MB_OK|MB_ICONSTOP);
			return 0;
		}
		else */
		return 1;
	 }
} // init

/*void CAfficheMesPol::OnFichierOuvrir() 
{
	CFileDialog * FileData1 = new CFileDialog (TRUE,"*.mmc",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"Mesures (*.mmc)|*.mmc||",this);
	if (FileData1->DoModal()==IDOK)
	{
	CFileDialog * FileData2 = new CFileDialog (TRUE,"*.pmc",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"Mesures (*.pmc)|*.pmc||",this);
	if (FileData2->DoModal()==IDOK)
	{
	if (PolyMc.Va!=NULL) PolyMc.FinAlloc();
	PolyMc.Init(FileData1->GetFileName(),FileData2->GetFileName());
	PolyMc.LirePol(	FileData2->GetFileName());
	Invalidate(FALSE);
   	}
	delete FileData2;

	}

	delete FileData1;

}
*/
void CAfficheMesPol::OnFichierQuitter() 
{
	PolyMc.FinAlloc();
	DestroyWindow();
}

void CAfficheMesPol::OnVisualisationChoixdesparametres() 
{
	if (PolyMc.NPAR==1) {param1=0;param0=0;}
	CDial2Val Dial2Val(DialPara,param0,param1);

	if (Dial2Val.DoModal()==IDOK)
	{
		if (Dial2Val.m_val0<0 || Dial2Val.m_val0>(PolyMc.NPAR-1)) {MessageBox("valeur impossible");return;}	
		if (Dial2Val.m_val1<0 || Dial2Val.m_val1>(PolyMc.NPAR-1)) {MessageBox("valeur impossible");return;}	
		param0=Dial2Val.m_val0;
		param1=Dial2Val.m_val1;
	}
	if(Parzen) OnParzen();
	Invalidate(TRUE);	
}

void CAfficheMesPol::OnVisualisationVisualiserlespolytopes() 
{
	if(Visu==TRUE) 
	{
		Visu=FALSE;
	}
	else 
	{
		Visu=TRUE;   
	}
	Invalidate(TRUE);
}


void CAfficheMesPol::ParametreGraphique()
{
	int decalage = 20; 
	if(IsWindow(GetSafeHwnd())) GetClientRect(rcClient);
	WGrD = rcClient.Width()-decalage;   
	WGrG = decalage;    
	HGrB = rcClient.Height()-decalage;   
	HGrH = decalage+20; //20 pour légende                    
	Ay = (float)(HGrH-HGrB)/255.0F; 
	By = (float)HGrB;//- Ay*(float)MinY;
		  
	Bx = (float)WGrG;
	Ax = (float)(WGrD-Bx)/255.0F;
}


int CAfficheMesPol::PosGrX(int X)
{
	return (int) (Ax*X+Bx);
}

int CAfficheMesPol::PosGrY(int Y)
{
	return (int) (Ay*Y+By);
}

void CAfficheMesPol::AxeGraph(CPaintDC &dc)
{
	dc.FillSolidRect(PosGrX(0),PosGrY(255),PosGrX(255)-20,PosGrY(0)-40,RGB(255,255,255));

	CPen PenAxe(PS_SOLID,1,RGB(0,0,0));
	dc.SelectObject(&PenAxe);
	dc.MoveTo(PosGrX(0),PosGrY(0));
	dc.LineTo(PosGrX(0),PosGrY(255));
	dc.MoveTo(PosGrX(0),PosGrY(0));
	dc.LineTo(PosGrX(255),PosGrY(0));

	for(int i=5;i<255;i+=5)
	{
		dc.MoveTo(PosGrX(0),PosGrY(i));
		dc.LineTo(PosGrX(3),PosGrY(i));
		dc.MoveTo(PosGrX(i),PosGrY(0));
		dc.LineTo(PosGrX(i),PosGrY(3));
	}
	DeleteObject(&PenAxe);
}

void CAfficheMesPol::OnSize(UINT nType, int cx, int cy) 
{
	if(IsWindow(GetSafeHwnd())) GetClientRect (rcClient); 
	if ((rcClient.Height()<230) || (rcClient.Width()<200))	
	{ 	
		SetWindowPos(NULL ,0,0,200,230,SWP_NOMOVE );       // SWP_NOMOVE        
		cx = 200;cy = 230;
	}
	else CMDIChildWnd::OnSize(nType, cx, cy);
}

void CAfficheMesPol::OnPaint() 
{
	CPaintDC dc(this); // device context for painting
	int x,y,cl;
	long i;
	int Rouge,Vert,Bleu;
	dc.SetBkMode(TRANSPARENT);
	dc.TextOut(20,0,"N° de classe:");
	dc.SetBkMode(OPAQUE);
	int x1,x2,y1,y2,anpar;
	ParametreGraphique();
	AxeGraph(dc);
	///////////////////////////////////
	CPen Pen(PS_SOLID,1,RGB(0,255,0));
	dc.SelectObject(&Pen);
	if (PolyMc.NPAR==1) {param1=0;param0=0;}

	if (Visu==TRUE)
	{
		anpar=PolyMc.NPAR;

		if (anpar!=PolyMc.NPAR) 
		{
			MessageBox("erreur sur le nombre de paramètres");
			return;
		}
		
		for (cl=0;cl<PolyMc.NbClass;cl++)
		{
			if((cl%3)==0) Rouge=255-(16*cl);else Rouge=0;
			if(((cl+1)%3)==0) Vert=255-(16*cl);  else Vert=0;
			if(((cl+2)%3)==0) Bleu=255-(16*cl);else Bleu=0;


			CPen Pen(PS_SOLID,1,RGB(Rouge,Vert,Bleu));

			dc.SelectObject(&Pen);

			for (i=0;i<PolyMc.Npol[cl];i++)
			{
				x1=PolyMc.poli[cl][i*PolyMc.NPAR+param0];
				x2=PolyMc.pols[cl][i*PolyMc.NPAR+param0];
				y1=PolyMc.poli[cl][i*PolyMc.NPAR+param1];
				y2=PolyMc.pols[cl][i*PolyMc.NPAR+param1];
				dc.MoveTo(PosGrX(x1),PosGrY(y1));
				dc.LineTo(PosGrX(x2),PosGrY(y1));
				dc.LineTo(PosGrX(x2),PosGrY(y2));
				dc.LineTo(PosGrX(x1),PosGrY(y2));
				dc.LineTo(PosGrX(x1),PosGrY(y1));
			}
		}

	} // si visu
	CString TextClass;
	dc.SetTextColor(RGB(255,255,255));
	for (cl=0;cl<PolyMc.NbClass;cl++)
	{
		if((cl%3)==0) Rouge=255-(16*cl);else Rouge=0;
		if(((cl+1)%3)==0) Vert=255-(16*cl);  else Vert=0;
		if(((cl+2)%3)==0) Bleu=255-(16*cl);else Bleu=0;

		dc.SetBkColor(RGB(Rouge,Vert,Bleu));
		TextClass.Format("%d",cl);
		dc.TextOut(30+cl*20,20,TextClass);
		
		for (i=0;i<PolyMc.NbMes[cl];i++)
		{
			x=(int)PolyMc.Va[cl][i*PolyMc.NPAR+param0];
			y=(int)PolyMc.Va[cl][i*PolyMc.NPAR+param1];
			dc.SetPixel(PosGrX(x),PosGrY(y),RGB(Rouge,Vert,Bleu));
			//dc.FillSolidRect(x,300-y,2,2,RGB(Rouge,Vert,Bleu));
		}

	} // cl
	// Do not call CMDIChildWnd::OnPaint() for painting messages
}

void CAfficheMesPol::OnVisualisationInfomesures() 
{
	PolyMc.AffMes();	
}

void CAfficheMesPol::OnVisualisationInfopolytopes() 
{
	PolyMc.AffPol();
}

BOOL CAfficheMesPol::Definirh() 
{
	CString* titre;
	titre=new CString("Entrez une valeur >0");
	CDialVal DialVal(titre);

	if(DialVal.DoModal()==IDOK)
	{
		Hparam=atoi(DialVal.m_val);
		if (Hparam<=0) {MessageBox("Valeur impossible","Attention");return Definirh();}
		Parzen=TRUE;
		return TRUE;
	}
	else return FALSE;
	

}

void CAfficheMesPol::OnParzen() 
{
	// TODO: Add your command handler code here
	long i;
	char ch[50],ch0[15],ch1[15],info[50];
	
	CString titre;
	titre.Format("Recouvrement histog. classes n°%d, n°%d",class1,class2);
	
	Graph=TRUE;
	if(!Parzen) 
	{
		if(!Definirh()) return;
	}
		
	if (PolyMc.NPAR==1) {param1=0;param0=0;}

	if(pParzen0!=NULL) { pParzen0->DestroyWindow();pParzen0=NULL;}
	if(pHisto0!=NULL)  {pHisto0->DestroyWindow();pHisto0=NULL;}
	if(pParzen1!=NULL) {pParzen1->DestroyWindow();pParzen1=NULL;}
	if(pHisto1!=NULL)  {pHisto1->DestroyWindow();pHisto1=NULL;}

	sprintf(info,"c0 en rouge ; c1 en bleu ; h=%d",Hparam);

	CalculParzen( param0,parzen00,parzen10,ch0);

	for (i=0;i<256;i++)
	{
		histo00[i]=0;
		histo10[i]=0;
	}

	
	for ( i=0;i<PolyMc.NbMes[class1];i++)
	{
		histo00[PolyMc.Va[class1][i*PolyMc.NPAR+param0]]++;
	}
	

	for ( i=0;i<PolyMc.NbMes[class2];i++)
	{
		histo10[PolyMc.Va[class2][i*PolyMc.NPAR+param0]]++;
	}
	
	if (param0!=param1) 
	{
		
		CalculParzen( param1,parzen01,parzen11,ch1);
		long i;
		for (i=0;i<256;i++)
		{
			histo01[i]=0;
			histo11[i]=0;
		}


		for ( i=0;i<PolyMc.NbMes[class1];i++)
		{
			histo01[PolyMc.Va[class1][i*PolyMc.NPAR+param1]]++;
		}
	
		for ( i=0;i<PolyMc.NbMes[class2];i++)
		{
			histo11[PolyMc.Va[class2][i*PolyMc.NPAR+param1]]++;
		}

		sprintf(ch,"Parzen :"+DialPara->m_choisi[param1]+" classe n°%d et n°%d",class1,class2);
		
		pParzen1 = (CGraph2Wnd*) new CGraph2Wnd(parzen01,parzen11,256,"Niveau","Nombre",info);  
		if (!pParzen1->Create(ch,GetMDIFrame())) return;

		sprintf(ch,"Histo. :"+DialPara->m_choisi[param1]+" classe n°%d et n°%d",class1,class2);
		
		pHisto1 = (CGraph2Wnd*) new CGraph2Wnd(histo01,histo11,256,"Niveau","Nombre","c0 en rouge ; c1 en bleu");  
		if (!pHisto1->Create(ch,GetMDIFrame())) return;


		if(!IsWindow(DialParz->GetSafeHwnd())) 
		{
			DialParz->Create(DialPara->m_choisi[param0]+":",ch0,DialPara->m_choisi[param1]+":",ch1);
			DialParz->SetWindowText(titre);
		}
		else
		{
			DialParz->Rafraichissement(DialPara->m_choisi[param0]+":",ch0,DialPara->m_choisi[param1]+":",ch1);
			DialParz->SetWindowText(titre);
			DialParz->ShowWindow(SW_SHOW);
		}
	}//if param0 != param1
	else 
	{
		if(!IsWindow(DialParz->GetSafeHwnd())) 
		{
			DialParz->Create(DialPara->m_choisi[param0]+":",ch0,"pas de deuxième paramètre","");
			DialParz->SetWindowText(titre);
		}
		else
		{
			DialParz->Rafraichissement(DialPara->m_choisi[param0]+":",ch0,"pas de deuxième paramètre","");
			DialParz->SetWindowText(titre);
			DialParz->ShowWindow(SW_SHOW);
		}	
	}
	
	sprintf(ch,"Parzen :"+DialPara->m_choisi[param0]+" classe n°%d et n°%d",class1,class2);
	
	pParzen0 = (CGraph2Wnd*) new CGraph2Wnd(parzen00,parzen10,256,"Niveau","Nombre",info);  
	if (!pParzen0->Create(ch,GetMDIFrame())) return;
	

	sprintf(ch,"Histo.:"+DialPara->m_choisi[param0]+" classe n°%d et n°%d",class1,class2);
	
	pHisto0 = (CGraph2Wnd*) new CGraph2Wnd(histo00,histo10,256,"Niveau","Nombre","c0 en rouge ; c1 en bleu");  
	if (!pHisto0->Create(ch,GetMDIFrame())) return;

	

	GetMDIFrame()->MDITile(MDITILE_VERTICAL);
	this->SetFocus();
}//OnParzen

void CAfficheMesPol::CalculParzen(int par,long * parzen0,long *parzen1 ,char * ch)
{
	
	long recouvr;
	double temp;
	int i;
	for (i=0;i<256;i++)
	{
		parzen0[i]=0;
		parzen1[i]=0;
	}
	
	for ( i=0;i<256;i++)
	{
		temp=0;
		for (long j=0;j<PolyMc.NbMes[class1];j++)
		{
			temp+=Estimateur(i,PolyMc.Va[class1][j*PolyMc.NPAR+par]);
		}
		temp/=PolyMc.NbMes[class1]*Hparam;
		parzen0[i]=(long)(probc1*temp);
		
	}
		
	for ( i=0;i<256;i++)
	{
		temp=0;
		for (long j=0;j<PolyMc.NbMes[class2];j++)
		{
			temp+=Estimateur(i,PolyMc.Va[class2][j*PolyMc.NPAR+par]);
		}
		temp/=PolyMc.NbMes[class2]*Hparam;
		parzen1[i]=(long)(probc2*temp);
	}

	recouvr=0;
	for(i=0;i<256;i++)
	{
		recouvr+=__min(parzen0[i],parzen1[i]);
	}

	double rapport= (double) recouvr/1000000.0;
	sprintf(ch,"%2.4f",rapport);
	

}//CalculParzen

double CAfficheMesPol::Estimateur (int x0,int xn)

{
//modif	
/*	if (abs(x0-xn)>(double) (Hparam*Hparam)/2.0) return 0;
	else return  25500*(2.0/(double)Hparam - (double)(4*abs(x0-xn))/(double)(Hparam*Hparam*Hparam));*/
	if (abs(x0-xn)> Hparam) return 0;
	else return  100000000*(1 - (double)(abs(x0-xn))/(double)(Hparam));
}

void CAfficheMesPol::OnParzenCalculdetous() 
{
	long parzen0[256];
	long parzen1[256];
	char ch[15];
	CString ch0;
	FILE *Fichier;
	
	if(!Parzen) 
	{
		if(!Definirh()) return;
	}

	CFileDialog FileData(FALSE,"*.xls",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"EXCEL (*.xls)|*.xls||",this);
	if (FileData.DoModal()==IDOK)
	{
		if (Fichier=fopen(FileData.GetPathName(),"w"))
		{
			SetCursor(LoadCursor(0,IDC_WAIT));
			fprintf(Fichier,"recouvrement des classes n°%d et n°%d en %% \n",class1,class2);
			fprintf(Fichier,"paramètre de parzen: %d\n",Hparam);
			for (int i=0;i<PolyMc.NPAR;i++)
			{
				CalculParzen(i,parzen0,parzen1,ch);
				
				ch0=ch;
				
				ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
				fprintf(Fichier,DialPara->m_choisi[i]+";"+ch0+"\n");
			}
			fclose(Fichier);
			SetCursor(LoadCursor(0,IDC_ARROW));
		}
	}   
	
}
/////////////////////////////////////////////////
////coeff corrélation
/////////////////////////////////////////////

void CAfficheMesPol::OnCorrelationCalcul() 
{
	
	long nbp,p=0;
	double  r;
	char ch[300];
	CString titre;
	nbp=PolyMc.NbMes[class1]+PolyMc.NbMes[class2];
	
	r=Correlation(param0,param1);
	sprintf(ch,DialPara->m_choisi[param0]+"  &  "+DialPara->m_choisi[param1]+"\n\n Coeffient de correlation = %1.4f",r);
	titre.Format("Correlation sur les classes n°%d et n°%d",class1,class2);
	MessageBox(ch,titre);
}
double CAfficheMesPol::Correlation(int par1,int par2)
{
	return (Covariance(par1,par2)/(pow(Variance(par1)*Variance(par2),0.5)));
}

double CAfficheMesPol::Covariance(int par1,int par2)
{
	double e,exy;
	exy=Esperance_xy(par1,par2);
	e=Esperance(par1)*Esperance(par2);
	return (exy-e);
}

double CAfficheMesPol::Esperance_xy(int par1, int par2)
{
	double e=0;
	long i;
	for (i=0;i<PolyMc.NbMes[class1];i++)
	 {
		e+=double(PolyMc.Va[class1][i*PolyMc.NPAR+par1])*double(PolyMc.Va[class1][i*PolyMc.NPAR+par2]);
	 }
	for (i=0;i<PolyMc.NbMes[class2];i++)
	 {
		e+=double(PolyMc.Va[class2][i*PolyMc.NPAR+par1])*double(PolyMc.Va[class2][i*PolyMc.NPAR+par2]);
	 }
	return (e/double(PolyMc.NbMes[class1]+PolyMc.NbMes[class2]));
}

double CAfficheMesPol::Esperance(int par)
{
	double e=0;
	long i;

	for (i=0;i<PolyMc.NbMes[class1];i++)
	{
		e+=double(PolyMc.Va[class1][i*PolyMc.NPAR+par]);
	}

	for (i=0;i<PolyMc.NbMes[class2];i++)
	{
		e+=double(PolyMc.Va[class2][i*PolyMc.NPAR+par]);
	}
	return (e/double(PolyMc.NbMes[class1]+PolyMc.NbMes[class2]));
}

double CAfficheMesPol::Variance(int par)
{
	double e,ex2;
	e=Esperance2(par);
	ex2=Esperance(par);
	ex2=ex2*ex2;
	return (e-ex2);
}

double CAfficheMesPol::Esperance2(int par)
{
	double e=0;
	long i;

	for (i=0;i<PolyMc.NbMes[class1];i++)
	 {
		e+=double(PolyMc.Va[class1][i*PolyMc.NPAR+par])*double(PolyMc.Va[class1][i*PolyMc.NPAR+par]);
	 }
	for (i=0;i<PolyMc.NbMes[class2];i++)
	 {
		e+=double(PolyMc.Va[class2][i*PolyMc.NPAR+par])*double(PolyMc.Va[class2][i*PolyMc.NPAR+par]);
	 }
	return (e/double(PolyMc.NbMes[class1]+PolyMc.NbMes[class2]));
}

void CAfficheMesPol::OnCorrelationCalculdans() 
{
	char ch[30];
	CString ch0;
	FILE *Fichier;
	
	CFileDialog FileData(FALSE,"*.xls",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"EXCEL (*.xls)|*.xls||",this);
	if (FileData.DoModal()==IDOK)
	{
		if (Fichier=fopen(FileData.GetPathName(),"w"))
		{
			SetCursor(LoadCursor(0,IDC_WAIT));
			fprintf(Fichier,"Coefficient de correlation des paramètres sur les classe n°%d et n°%d\n",class1,class2);
			for (int i=0;i<PolyMc.NPAR-1;i++)
			{
				for(int j=i+1;j<PolyMc.NPAR;j++)
				{
				sprintf(ch,"%1.4f",Correlation(i,j));
				
				ch0=ch;
				
				ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
				fprintf(Fichier,DialPara->m_choisi[i]+" & "+DialPara->m_choisi[j]+";"+ch0+"\n");
				}
			}
			fclose(Fichier);
			SetCursor(LoadCursor(0,IDC_ARROW));
		}
	}   
	
}
/////////////////////////////////////////////////
////***********************************************
/////////////////////////////////////////////
void CAfficheMesPol::OnWindowClosegraphics() 
{
	if(IsWindow(pParzen0->GetSafeHwnd()))  pParzen0->DestroyWindow();
	if(IsWindow(pHisto0->GetSafeHwnd()))  pHisto0->DestroyWindow();
	if(IsWindow(pParzen1->GetSafeHwnd()))  pParzen1->DestroyWindow();
	if(IsWindow(pHisto1->GetSafeHwnd()))  pHisto1->DestroyWindow();
	pParzen0=NULL;
	pHisto0=NULL;
	pParzen1=NULL;
	pHisto1=NULL;
	Parzen=FALSE;	
}


///////////////////////gestion menu/////////////////////////////
void CAfficheMesPol::OnUpdateWindowClosegraphics(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Graph);
}

void CAfficheMesPol::OnUpdateVisualisationVisualiserlespolytopes(CCmdUI* pCmdUI) 
{
	pCmdUI->SetCheck(Visu);
	pCmdUI->Enable(ValidPol);
}

void CAfficheMesPol::OnUpdateCorrelationCalcul(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable((PolyMc.NPAR!=1)&&ValidCalcul);	
}

BOOL CAfficheMesPol::DestroyWindow() 
{
	OnWindowClosegraphics();
	CToolBarCtrl& pToolBarCtrl=(&((CMainFrame*) GetMDIFrame())->bar)->GetToolBarCtrl();
	pToolBarCtrl.HideButton(ID_VISUALISATION_CHOIXDESPARAMETRES,TRUE);
	pToolBarCtrl.HideButton(ID_VISUALISATION_VISUALISERLESPOLYTOPES,TRUE);
	pToolBarCtrl.HideButton(ID_PARZEN,TRUE);
	((CMainFrame*)GetMDIFrame())->RecalcLayout();
	

	return CMDIChildWnd::DestroyWindow();
}

void CAfficheMesPol::OnUpdateCorrelationCalculdans(CCmdUI* pCmdUI) 
{
		pCmdUI->Enable((PolyMc.NPAR!=1)&&ValidCalcul);
}

void CAfficheMesPol::OnVisualisationToolbar() 
{
	CToolBar* pToolBar=& ((CMainFrame*)GetMDIFrame())->bar;
	ShowControlBar(pToolBar,(pToolBar->GetStyle() & WS_VISIBLE)==0, FALSE);
	
}

void CAfficheMesPol::OnUpdateVisualisationToolbar(CCmdUI* pCmdUI) 
{
	// TODO: Add your command update UI handler code here
	pCmdUI->SetRadio((((CMainFrame*)GetMDIFrame())->bar.GetStyle() & WS_VISIBLE)!=0);
}


void CAfficheMesPol::OnDestroy() 
{
	CMDIChildWnd::OnDestroy();
	
	// TODO: Add your message handler code here
	if(GetMDIFrame()->MDIGetActive()==NULL)
	{
		CMenu menutemp;
		menutemp.LoadMenu(IDR_MAINFRAME);
		GetMDIFrame()->MDISetMenu(&menutemp,NULL);
		AfxGetMainWnd()->DrawMenuBar();
	}
	(((CMainFrame*) GetMDIFrame())->bar.GetToolBarCtrl()).SetOwner( (CWnd*)GetMDIFrame() );
	((CMainFrame*) GetMDIFrame())->bar.OnUpdateCmdUI( (CFrameWnd*) GetMDIFrame(), TRUE);
}

void CAfficheMesPol::OnVisualisationChoixdunfichierdepolytopes() 
{
	BYTE NbPar=PolyMc.NPAR;
	CPolyMc test;
	CFileDialog * FileData = new CFileDialog (TRUE,"*.pmc",NomFichierPol,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"Polytopes multi-classes (*.pmc)|*.pmc||",this);
	if(FileData->DoModal()==IDOK)
	{
		test.LirePol(FileData->GetPathName());
		test.Fin();
		if(PolyMc.NPAR!=test.NPAR||PolyMc.NbClass!=test.NbClass)
		{
			MessageBox("Fichier polytope incompatible!","Erreur!",MB_OK|MB_ICONSTOP);
			ValidPol=FALSE;
			Visu=FALSE;
		}
		else 
		{
			ValidPol=TRUE;
			NomFichierPol=FileData->GetPathName();
			PolyMc.Fin();
			ValidPol=PolyMc.LirePol(NomFichierPol);
		}
	}
	Invalidate(TRUE);
	delete FileData;
}

void CAfficheMesPol::OnUpdateVisualisationInfopolytopes(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValidPol);	
}



void CAfficheMesPol::OnCalculsChoixdeclasses() 
{
	CDialChoixClass Dlg(this,class1,class2,PolyMc.NbClass,probc1,probc2);
	if(Dlg.DoModal()==IDOK)
	{
		class1=Dlg.m_class1;
		class2=Dlg.m_class2;
		probc1=Dlg.m_probc1;
		probc2=Dlg.m_probc2;
		menu.GetSubMenu(1)->EnableMenuItem(2,MF_ENABLED|MF_BYPOSITION);	
		menu.GetSubMenu(1)->EnableMenuItem(3,MF_ENABLED|MF_BYPOSITION);
		ValidCalcul=TRUE;
		if(Parzen) OnParzen();
	}

/*	CString ch;
	ch.Format("probc1=%1.4f\nprobc2=%1.4f",probc1,probc2);
	MessageBox(ch);*/
}

int CAfficheMesPol::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CMDIChildWnd::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	menu.GetSubMenu(1)->EnableMenuItem(2,MF_GRAYED|MF_BYPOSITION);	
	menu.GetSubMenu(1)->EnableMenuItem(3,MF_GRAYED|MF_BYPOSITION);
	return 0;
}

void CAfficheMesPol::OnUpdateParzen(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValidCalcul);		
}

void CAfficheMesPol::OnUpdateParzenCalculdetous(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValidCalcul);	
}

void CAfficheMesPol::OnMDIActivate(BOOL bActivate, CWnd* pActivateWnd, CWnd* pDeactivateWnd) 
{
	CMDIChildWnd::OnMDIActivate(bActivate, pActivateWnd, pDeactivateWnd);
	
	(((CMainFrame*) GetMDIFrame())->bar.GetToolBarCtrl()).SetOwner((CWnd*)this );
	((CMainFrame*) GetMDIFrame())->bar.OnUpdateCmdUI( (CFrameWnd*) GetMDIFrame(), FALSE);
	
}
