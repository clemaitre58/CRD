// DialContour.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "DialContour.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDialContour dialog


CDialContour::CDialContour(unsigned char seuilhaut,unsigned char seuilbas,long longmax,long seuilcont,CWnd* pParent /*=NULL*/)
	: CDialog(CDialContour::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialContour)
	m_seuilbas = seuilbas;
	m_seuilhaut = seuilhaut;
	m_longmax = longmax;
	m_seuilcont = seuilcont;
	m_check = FALSE;
	//}}AFX_DATA_INIT
}


void CDialContour::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialContour)
	DDX_Text(pDX, IDC_EDIT3, m_seuilbas);
	DDX_Text(pDX, IDC_EDIT1, m_seuilhaut);
	DDX_Text(pDX, IDC_EDIT4, m_longmax);
	DDX_Text(pDX, IDC_EDIT5, m_seuilcont);
	DDX_Check(pDX, IDC_CHECK1, m_check);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialContour, CDialog)
	//{{AFX_MSG_MAP(CDialContour)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialContour message handlers

void CDialContour::OnCancel() 
{
	// TODO: Add extra cleanup here
	
	CDialog::OnCancel();
}

void CDialContour::OnOK() 
{
	UpdateData();
	if((m_seuilbas<1)||(m_seuilbas>255)||(m_seuilhaut<1)||(m_seuilhaut>255)||
		(m_longmax<1)||(m_longmax>4000)||(m_seuilcont<1)||(m_seuilcont>4000))
		MessageBox("Valeurs impossibles !");
	else CDialog::OnOK();
}
