// DialLegend.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "DialLegend.h"


/////////////////////////////////////////////////////////////////////////////
// CDialLegend dialog


CDialLegend::CDialLegend(CWnd* pParent /*=NULL*/)
	: CDialog(CDialLegend::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialLegend)
		// NOTE: the ClassWizard will add member initialization here
	//}}AFX_DATA_INIT
}


void CDialLegend::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialLegend)
		// NOTE: the ClassWizard will add DDX and DDV calls here
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialLegend, CDialog)
	//{{AFX_MSG_MAP(CDialLegend)
	ON_WM_CREATE()
	ON_WM_PAINT()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialLegend message handlers

int CDialLegend::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CDialog::OnCreate(lpCreateStruct) == -1)
		return -1;
		SetWindowPos(NULL ,0,0,150,70+Nb*20,SWP_NOMOVE ); 
	// TODO: Add your specialized creation code here
	
	return 0;
}



void CDialLegend::OnPaint() 
{
	CPaintDC dc(this); // device context for painting
	
	CString text;
	dc.SetBkMode(TRANSPARENT);
	int i;
	for(i=0;i<Nb;i++)
	{
		text.Format("Classe n°%d",i);
		dc.TextOut(40,20+20*i,text);
	}
	if(NonClass) 
	{
		text.Format("Non classés");
		dc.TextOut(40,20+20*i,text);
	}
	dc.SetBkMode(OPAQUE);
	for( i=0;i<Nb;i++)
	{
		dc.FillSolidRect(20,22+20*i,10,10,Palette_Couleurs(i));
	}
	if(NonClass) 
	{
		dc.FillSolidRect(20,22+20*i,10,10,Palette_Couleurs(i));
	}

	
	// TODO: Add your message handler code here
	
	// Do not call CDialog::OnPaint() for painting messages
}

COLORREF  CDialLegend::Palette_Couleurs(int i)
{
	 int vert,bleu,rouge;
	 unsigned char l;
	 int nombre=Nb;
	 if(NonClass)nombre++;
	 l=(unsigned char) (i*((double)255/nombre));

	if ((l>=0)&&(l<=33))
	{
					vert  = (int)0;
					bleu  = (int)(255*(0-l)/(0-33));
					rouge = (int)(255*(0-l)/(0-33)); 
	}else
	if ((l>=34)&&(l<=67))
	{
				 vert = 0;
				 bleu = 255;
				 rouge = (int)(255*(77-l)/(77-34));
	}else   
	if ((l>=68)&&(l<=100) )
	{
				 vert  = (int)(255*(68-l)/(68-100));
				 bleu  = (int)255;
				 rouge = (int)0;
	}else 
	if( (l>=101)&&(l<=134) )
	{
				 vert = 255;
				 bleu = (int)(255.0*(float)(134-l)/(float)(134-101));
				 rouge = 0;
	}else
	if ((l>=135)&&(l<=177))
	{ 
				 vert  = 255;
				bleu  = 0;
				 rouge = (int)(255.0*(float)(135-l)/(float)(135-177));
	}else                 
	if ((l>=178)&&(l<=201))
	{
		  vert = (int)(255.0*(float)(201-l)/(float)(201-178));
		  bleu = 0;
		  rouge =255;
	}else	  
	if ((l>=202)&&(l<=255))
	{
		  vert = (int)(255.0*(float)(202-l)/(float)(202-255));
		  bleu = (int)(255.0*(float)(202-l)/(float)(202-255));
		  rouge = 255;
	}


	return RGB(rouge,vert,bleu);
 }

BOOL CDialLegend::Create(int initnb,BOOL initaffectnclas) 
{
	// TODO: Add your specialized code here and/or call the base class
	Nb=initnb;
	NonClass=initaffectnclas;

	return CDialog::Create(IDD);
}

void CDialLegend::Init(int initnb,BOOL initnnclas)
{
	Nb=initnb;
	NonClass=initnnclas;
	SetWindowPos(NULL ,0,0,150,70+Nb*20,SWP_NOMOVE ); 
	Invalidate(TRUE);
}
