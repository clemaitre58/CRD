// DialCurvi.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "DialCurvi.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDialCurvi dialog


CDialCurvi::CDialCurvi(CWnd* pParent /*=NULL*/)
	: CDialog(CDialCurvi::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialCurvi)
	m_TE = 3;
	m_MaxLength = 3000;
	m_MaxWidth = 6;
	m_MaxDw = 6;
	m_CheckDw = FALSE;
	m_MinLength = 10;
	m_TxtMaxW = _T("");
	m_Reg = 0.8f;
	m_Check_Reg = FALSE;
	//}}AFX_DATA_INIT
}


void CDialCurvi::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialCurvi)
	DDX_Text(pDX, IDC_EDIT_TE, m_TE);
	DDX_Text(pDX, IDC_EDIT_ML, m_MaxLength);
	DDX_Text(pDX, IDC_EDIT_MW, m_MaxWidth);
	DDX_Text(pDX, IDC_EDIT_DW, m_MaxDw);
	DDX_Check(pDX, IDC_CHECK1, m_CheckDw);
	DDX_Text(pDX, IDC_EDIT_ML2, m_MinLength);
	DDX_Text(pDX, IDC_TXTMAXW, m_TxtMaxW);
	DDX_Text(pDX, IDC_EDIT_REG, m_Reg);
	DDV_MinMaxFloat(pDX, m_Reg, 0.f, 1.f);
	DDX_Check(pDX, IDC_CHECK_REG, m_Check_Reg);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialCurvi, CDialog)
	//{{AFX_MSG_MAP(CDialCurvi)
	ON_EN_CHANGE(IDC_EDIT_MW, OnChangeEditMw)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialCurvi message handlers

void CDialCurvi::OnChangeEditMw() 
{
	// TODO: If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.
	
	// TODO: Add your control notification handler code here
	UpdateData(TRUE);
	if (Larg==0) return;
	int a=(Larg*m_MaxWidth)/100;
	if (a>(Larg/2)) m_TxtMaxW="Forbidden";
	else m_TxtMaxW.Format("%d pixels",a);
	UpdateData(FALSE);

}

BOOL CDialCurvi::OnInitDialog() 
{
	CDialog::OnInitDialog();
	
	// TODO: Add extra initialization here
	if (Larg!=0) 
	{
	 int a=(Larg*m_MaxWidth)/100;
	 m_TxtMaxW.Format("%d pixels",a);
	 UpdateData(FALSE);
	}
	return TRUE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}
