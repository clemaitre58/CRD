#include "stdafx.h"
#include "winseg.h"
#include "dialogex.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDialogExtrac dialog


CDialogExtrac::CDialogExtrac(CWnd* pParent /*=NULL*/)
	: CDialog(CDialogExtrac::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialogExtrac)
	m_profondeur = 2;
	//}}AFX_DATA_INIT
}

void CDialogExtrac::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialogExtrac)
	DDX_Text(pDX, IDC_REPONSE, m_profondeur);
	DDV_MinMaxInt(pDX, m_profondeur, 1, 5);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CDialogExtrac, CDialog)
	//{{AFX_MSG_MAP(CDialogExtrac)
		// NOTE: the ClassWizard will add message map macros here
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CDialogExtrac message handlers
