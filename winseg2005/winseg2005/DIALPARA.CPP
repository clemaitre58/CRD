// dialpara.cpp : implementation file
//

#include "stdafx.h"
//#include "timage.h"
#include "dialpara.h"



/////////////////////////////////////////////////////////////////////////////
// CDialPara dialog


CDialPara::CDialPara(CWnd* pParent /*=NULL*/)
	: CDialog(CDialPara::IDD, pParent)
{
 int i;		
 nbpara=0;                              
 for (i=0;i<256;i++) IChoisi[i]=0xff;
   	m_dispo[0] = "Luminance pixel:n"; 
	m_dispo[1] = "Lum moy fen glis:n,n";
	m_dispo[2] = "Lum moy fen fix:n,n";
	m_dispo[3] = "Gradient pixel:n"; 
	m_dispo[4] = "Grad moy fen glis:n,n";
	m_dispo[5] = "Grad moy fen fix:n,n";
	m_dispo[6] = "Lum moy décal lign fen fix:n,n";
	m_dispo[7] = "Entropie fix :n,n,n";
	m_dispo[8] = "Entropie glis :n,n,n";
	m_dispo[9] = "Contrast fix :n,n";
	m_dispo[10]= "Max_min fix :n,n";
	m_dispo[11]= "Contrast glis :n,n";
	m_dispo[12]= "Max_min glis :n,n";
	m_dispo[13]= "Median glis :n,n,n";
	m_dispo[14]= "Median fix :n,n,n";
	m_dispo[15]= "Image Enregistrée:";
	m_dispo[16]= "Sobel moy fen fix:n,n";
	m_dispo[17]= "Sobel moy fen glis:n,n";
	m_dispo[18]= "Rehaussement fix:n,n";
	m_dispo[19]= "Rehaussement glis:n,n";
	m_dispo[20]= "Ouverture :n,n";
	m_dispo[21]= "Fermeture :n,n";
	m_dispo[22]= "Ecart/Median fen glis:n,n,n";
	m_dispo[23]= "Ecart/Median fen fix:n,n,n";
	m_dispo[24]= "Moy. adap.& seuil:n,n,n";
	m_dispo[25]= "Moy. adap.auto.:n,n";
	m_dispo[26]= "Moy. tronquée glis:n,n,n";
	m_dispo[27]= "Moy. tronquée fix:n,n,n";
	m_dispo[28]= "Max histo glis:n,n";
	m_dispo[29]= "NbMax histo glis:n,n";
	m_dispo[30]= "Grad moy fen glis2:n,n";
	m_dispo[31]= "Fermeture rapide:n,n";
	m_dispo[32]= "Ouverture rapide:n,n";
	m_dispo[33]= "Contrast rapide glis:n,n";
	m_dispo[34]= "Grad rapide glis:n,n";
	m_dispo[35]= "Lum rapide glis:n,n";
	m_dispo[36]= "Grad rapide NonNorm:n,n";
	m_dispo[37]= "Ecart typ fen glis:n,n";
	m_dispo[38]= "Max fen glis:n,n";
	m_dispo[39]= "Densite max loc 3,3:n,n";
	//{{AFX_DATA_INIT(CDialPara)
	//}}AFX_DATA_INIT
}

void CDialPara::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialPara)
//	DDX_LBString(pDX, IDC_LIST1, m_dispo);
//	DDX_LBString(pDX, IDC_LIST2, m_choisi);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CDialPara, CDialog)
	//{{AFX_MSG_MAP(CDialPara)
	ON_LBN_DBLCLK(IDC_LIST1, OnDblclkList1)
	ON_LBN_DBLCLK(IDC_LIST2, OnDblclkList2)
	ON_BN_CLICKED(IDC_AJOUTER, OnAjouter)
	ON_BN_CLICKED(IDC_ENLEVER, OnEnlever)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CDialPara message handlers



void CDialPara::OnOK()
{
	// TODO: Add extra validation here
	
	CDialog::OnOK();
}

void CDialPara::OnCancel()
{
	// TODO: Add extra cleanup here
	
	CDialog::OnCancel();
}

BOOL CDialPara::OnInitDialog()
{
	CDialog::OnInitDialog();
	
	// TODO: Add extra initialization here
	DWORD dwIndex[256];
    int i;
    
//    nbpara=0;
    	

    
//    for (i=0;i<256;i++) IChoisi[i]=0xff;
   
    for (i=0;i<40;i++)
    {
	dwIndex[i] = SendDlgItemMessage(IDC_LIST1,LB_ADDSTRING, 0, (LPARAM) ((LPSTR) m_dispo[i].GetBuffer(40)));
    }
	for (i=0;i<nbpara;i++)
    {
	dwIndex[i] = SendDlgItemMessage(IDC_LIST2,LB_ADDSTRING, 0, (LPARAM) ((LPSTR) m_choisi[i].GetBuffer(40)));
    }

	return TRUE;  // return TRUE  unless you set the focus to a control

}


void CDialPara::OnDblclkList1()
{
	// TODO: Add your control notification handler code here
	CDialVal *DialVal;
	//CDial2Val *Dial2Val;
	CString Cval;
	int Pos;
	LRESULT lr;
	
	lr=SendDlgItemMessage(IDC_LIST1,LB_GETCURSEL, 0, 0);
   	if(lr!=LB_ERR)
	{
	Pos=m_dispo[lr].Find(":");
	if (Pos==-1) {MessageBox("erreur dans le parametre","Message",MB_OK);return;}
	Cval=m_dispo[lr].Right(m_dispo[lr].GetLength()-Pos);
	//MessageBox(Cval,"",MB_OK);
	m_choisi[nbpara]=m_dispo[lr].Left(Pos+1);

	do
	{
		DialVal=new CDialVal();
		if (DialVal->DoModal()==IDOK) 
   		{
   			if (DialVal->m_val.Left(1)!="+" && DialVal->m_val.Left(1)!="-") 	m_choisi[nbpara]+="+";
			m_choisi[nbpara]+=DialVal->m_val;
			
	 	}  // si OK
		else return;

		delete DialVal;
		
		Pos=Cval.Find(",");
		if (Pos!=-1) m_choisi[nbpara]+=",";
		Cval=Cval.Right(Cval.GetLength()-Pos-1);
		//MessageBox(Cval,"",MB_OK);
	} while(Pos!=-1);
   	SendDlgItemMessage(IDC_LIST2,LB_ADDSTRING, 0, (LPARAM) ((LPSTR) m_choisi[nbpara].GetBuffer(40)));
	IChoisi[nbpara]=lr;
	nbpara++;
	}
	
} //   si doubleclk



void CDialPara::OnDblclkList2()
{
	// TODO: Add your control notification handler code here
	LRESULT lr,i;
	lr=SendDlgItemMessage(IDC_LIST2,LB_GETCURSEL, 0, 0);
	if(lr!=LB_ERR)
	{
  // 	IChoisi[lr]=0xff;
    if (nbpara>0)
    {
    	for (i=lr;i<nbpara;i++)
    	{
    		m_choisi[i]=m_choisi[i+1];
    		IChoisi[i]=IChoisi[i+1];
    	}
    	nbpara--;
    	SendDlgItemMessage(IDC_LIST2,LB_DELETESTRING, (WPARAM)lr, 0);
	} // si reste para
	}
}

void CDialPara::OnAjouter()
{
	// TODO: Add your control notification handler code here
	OnDblclkList1();
}

void CDialPara::OnEnlever()
{
	OnDblclkList2();	
}
/////////////////////////////////////////////////////////////////////////////
// CDialVal dialog


CDialVal::CDialVal(CString* initTitre,CWnd* pParent /*=NULL*/)
	: CDialog(CDialVal::IDD, pParent)
{
	Titre=initTitre;
	//{{AFX_DATA_INIT(CDialVal)
	m_val = "";
	//}}AFX_DATA_INIT	
}

void CDialVal::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialVal)
	DDX_Text(pDX, IDC_EDIT1, m_val);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CDialVal, CDialog)
	//{{AFX_MSG_MAP(CDialVal)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

BOOL CDialVal::OnInitDialog() 
{
	if(Titre!=NULL) SetWindowText(*Titre);
	return CDialog::OnInitDialog();
}
void CDialVal::OnOK()
{
	// TODO: Add extra validation here
	
	CDialog::OnOK();
}

void CDialVal::OnCancel()
{
	// TODO: Add extra cleanup here
	
	CDialog::OnCancel();
}

/////////////////////////////////////////////////////////////////////////////
// CDial2Val dialog
/////////////////////////////////////////////////////////////////////////////

CDial2Val::CDial2Val(CDialPara* initdialpara,int param0,int param1,CWnd* pParent /*=NULL*/)
	: CDialog(CDial2Val::IDD, pParent)
{
	m_pDialPara=initdialpara;
	m_val0=param0;
	m_val1=param1;
}

void CDial2Val::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDial2Val)
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CDial2Val, CDialog)
	//{{AFX_MSG_MAP(CDial2Val)
	ON_CBN_SELCHANGE(IDC_COMBO1, OnSelchangeCombo1)
	ON_CBN_SELCHANGE(IDC_COMBO2, OnSelchangeCombo2)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


/////////////////////////////////////////////////////////////////////////////
// CDial2Val message handlers



void CDial2Val::OnCancel()
{
	// TODO: Add extra cleanup here
	CDialog::OnCancel();
}


BOOL CDial2Val::OnInitDialog() 
{
	
	pCombo1=(CComboBox *) GetDlgItem(IDC_COMBO1);
	pCombo2=(CComboBox *) GetDlgItem(IDC_COMBO2);
    for (int i=0;i<m_pDialPara->nbpara;i++)
    {
	pCombo1->InsertString (i, m_pDialPara->m_choisi[i]);
	pCombo2->InsertString (i, m_pDialPara->m_choisi[i]);
    }
	pCombo1->SelectString(m_val0, m_pDialPara->m_choisi[m_val0]);
	pCombo2->SelectString(m_val1, m_pDialPara->m_choisi[m_val1]);
	return CDialog::OnInitDialog();
}

void CDial2Val::OnSelchangeCombo1() 
{
	m_val0=pCombo1->GetCurSel();
}

void CDial2Val::OnSelchangeCombo2() 
{
	m_val1=pCombo2->GetCurSel();
}

void CDial2Val::OnOK() 
{
	// TODO: Add extra validation here
	
	CDialog::OnOK();
}
