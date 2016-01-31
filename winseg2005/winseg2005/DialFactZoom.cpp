#include "stdafx.h"
#include "winseg.h"
#include "DialFactZoom.h"


/////////////////////////////////////////////////////////////////////////////
// CDialFactZoom dialog


CDialFactZoom::CDialFactZoom(CImageWnd* pParent /*=NULL*/)
	: CDialog(CDialFactZoom::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialFactZoom)
	m_Facteur_Zoom = 1;
	m_pSet = pParent;
	//}}AFX_DATA_INIT
}


void CDialFactZoom::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialFactZoom)
	DDX_Text(pDX, IDC_EDIT_FACTEUR_ZOOM, m_Facteur_Zoom);
	DDV_MinMaxInt(pDX, m_Facteur_Zoom, 1, 10);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialFactZoom, CDialog)
	//{{AFX_MSG_MAP(CDialFactZoom)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialFactZoom message handlers

void CDialFactZoom::OnCancel() 
{
	// TODO: Add extra cleanup here
	
	CDialog::OnCancel();
}

void CDialFactZoom::OnOK() 
{
	// TODO: Add extra validation here
	UpdateData();
	m_pSet->m_FacteurZoom = m_Facteur_Zoom;
	CDialog::OnOK();
}
