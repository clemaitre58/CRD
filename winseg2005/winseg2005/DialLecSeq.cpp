// DialLecSeq.cpp : implementation file
//

#include "stdafx.h"
#include "winseg.h"
#include "DialLecSeq.h"
#include "image.h"



/////////////////////////////////////////////////////////////////////////////
// CDialLecSeq dialog

/*
CDialLecSeq::CDialLecSeq(CWnd* pParent /*=NULL)
	: CDialog(CDialLecSeq::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDialLecSeq)
	m_edit = 0;
	//}}AFX_DATA_INIT
}*/
CDialLecSeq::CDialLecSeq(CImageWnd* pParent)
{
	m_pView=pParent;
	m_edit = 1;
}

void CDialLecSeq::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDialLecSeq)
	DDX_Text(pDX, IDC_EDIT1, m_edit);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDialLecSeq, CDialog)
	//{{AFX_MSG_MAP(CDialLecSeq)
	ON_BN_CLICKED(IDC_BUTTON1, OnButton1)
	ON_BN_CLICKED(IDC_BUTTON2, OnButton2)
	ON_BN_CLICKED(IDC_BUTTON3, OnButton3)
	ON_WM_TIMER()
	ON_WM_HSCROLL()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDialLecSeq message handlers

void CDialLecSeq::OnButton1() 
{
	m_pView->AfficheSequence(0,Nbimseq);
	m_edit=Nbimseq;
	CSliderCtrl *pSL=(CSliderCtrl *) GetDlgItem(IDC_SLIDER1);
	pSL->SetPos(m_edit);
	UpdateData(FALSE);
}

void CDialLecSeq::OnButton2() 
{
	SetTimer(0,100,NULL);	
	CButton* pBut=(CButton *) GetDlgItem(IDC_BUTTON2);
	pBut->SetState( 1 );
}

void CDialLecSeq::OnButton3() 
{
	m_edit =1;
	UpdateData(FALSE);
	CSliderCtrl *pSL=(CSliderCtrl *) GetDlgItem(IDC_SLIDER1);
	pSL->SetPos(m_edit);
	m_pView->AfficheSequence(m_edit-1,m_edit);
}

void CDialLecSeq::OnTimer(UINT nIDEvent) 
{
	if(m_edit==Nbimseq) Finlect();
	else 
	{
		m_edit ++;
		UpdateData(FALSE);
		CSliderCtrl *pSL=(CSliderCtrl *) GetDlgItem(IDC_SLIDER1);
		pSL->SetPos(m_edit);
		m_pView->AfficheSequence(m_edit-1,m_edit);
	}
	CDialog::OnTimer(nIDEvent);
}

void CDialLecSeq::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar) 
{
	// TODO: Add your message handler code here and/or call default
	switch(nSBCode)
	{
	case SB_THUMBTRACK:
		m_edit = nPos;
		UpdateData(FALSE);
		break;
	case SB_PAGEDOWN:
		m_edit = ((CSliderCtrl*)pScrollBar)->GetPos();
		UpdateData(FALSE);
		break;
	case TB_THUMBPOSITION :
		m_edit = nPos;
		UpdateData(FALSE);
		break;
	case SB_PAGEUP:
		m_edit = ((CSliderCtrl*)pScrollBar)->GetPos();
		UpdateData(FALSE);
		break;
	case SB_LINEUP:
		m_edit = ((CSliderCtrl*)pScrollBar)->GetPos();
		UpdateData(FALSE);
		break;
	case SB_LINEDOWN:
		m_edit = ((CSliderCtrl*)pScrollBar)->GetPos();
		UpdateData(FALSE);
		break;

	}
	m_pView->AfficheSequence(m_edit-1,m_edit);
	CDialog::OnHScroll(nSBCode, nPos, pScrollBar);
}

void CDialLecSeq::OnOK() 
{
	Finlect();
//	CDialog::OnOK();

}

void CDialLecSeq::OnCancel() 
{
	ShowWindow(SW_HIDE);
//	CDialog::OnCancel();
}

BOOL CDialLecSeq::OnInitDialog() 
{
	CDialog::OnInitDialog();
	SetWindowText(LPCTSTR (m_titre));
	CSliderCtrl *pSL=(CSliderCtrl *) GetDlgItem(IDC_SLIDER1);
	pSL->SetRange(1,Nbimseq);
	if(Nbimseq>40)pSL->SetTicFreq(2);
	else pSL->SetTicFreq(1);
	pSL->SetLineSize(1);
	pSL->SetPageSize(5);
	pSL->SetSelection( 1, Nbimseq);
	CString ch;
	ch.Format("%d images",Nbimseq);
	SetDlgItemText(IDC_TEXT,ch);
	
	return TRUE;  // return TRUE unless you set the focus to a control
	              // EXCEPTION: OCX Property Pages should return FALSE
}

BOOL CDialLecSeq::Create(int nbseq,CString Titre) 
{
	Nbimseq=nbseq;
	m_edit=1;
	m_titre="Lecture de "+Titre;
	return CDialog::Create(CDialLecSeq::IDD);
}

void CDialLecSeq::Finlect()
 {
	 KillTimer(0);
	 CButton* pBut=(CButton *) GetDlgItem(IDC_BUTTON2);
	 pBut->SetState( 0);
 }

void CDialLecSeq::Redefinition(int NewNbimseq)
{
	if(NewNbimseq!=Nbimseq)
	{
		Nbimseq=NewNbimseq;
		CSliderCtrl *pSL=(CSliderCtrl *) GetDlgItem(IDC_SLIDER1);
		pSL->SetRange(1,Nbimseq);
		if(Nbimseq>40)pSL->SetTicFreq(2);
		else pSL->SetTicFreq(1);
		pSL->SetLineSize(1);
		pSL->SetPageSize(5);
		pSL->SetSelection( 1, Nbimseq);
		CString ch;
		ch.Format("%d images",Nbimseq);
		SetDlgItemText(IDC_TEXT,ch);
	}
}