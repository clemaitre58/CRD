   
#include "stdafx.h"
#include "resource.h"
#include "dialogue.h"        
     
               
// *************************************************************************************
// *********          choix d'une zone rectangle                  **********************
// *************************************************************************************
BEGIN_MESSAGE_MAP(CDialZone, CDialog)     
   	//{{AFX_MSG_MAP(CDialZone)
   	ON_COMMAND(IDOK, OnOK) 
   	ON_COMMAND(IDCANCEL, OnCancel)
   	//}}AFX_MSG_MAP
END_MESSAGE_MAP()                            


CDialZone::CDialZone(long *initi1,long *initj1,long *initi2,long *initj2,CWnd *pParent): CDialog(IDD_ZONE,NULL)
{ 
  LeParent=pParent;
  i1=initi1; j1=initj1;i2=initi2;j2=initj2; 
}

BOOL CDialZone::OnInitDialog()
{
  char ch1[10],ch2[10],ch3[10],ch4[10];  
  sprintf(ch1,"%ld",*i1);                           
  SetDlgItemText(IDC_i1,(LPSTR) ch1); 
  sprintf(ch2,"%ld",*j1);      
  SetDlgItemText(IDC_j1,(LPSTR) ch2);  
  sprintf(ch3,"%ld",*i2);  
  SetDlgItemText(IDC_i2,(LPSTR) ch3);
  sprintf(ch4,"%ld",*j2);
  SetDlgItemText(IDC_j2,(LPSTR) ch4); 
 return CDialog::OnInitDialog(); 
}
  
  
void CDialZone::OnOK( )
{ 
 
  long x1,x2,y1,y2;
  char ch1[10],ch2[10],ch3[10],ch4[10];                             
  GetDlgItemText(IDC_i1,(LPSTR) ch1,sizeof(ch1)); 
  x1=atol(ch1);
  GetDlgItemText(IDC_j1,(LPSTR) ch2,sizeof(ch2)); 
  y1=atol(ch2);
  GetDlgItemText(IDC_i2,(LPSTR) ch3,sizeof(ch3)); 
  x2=atol(ch3);
  GetDlgItemText(IDC_j2,(LPSTR) ch4,sizeof(ch4)); 
  y2=atol(ch4);  
  div_t res;  
  res = div((int)(x2-x1+1),4);
  if ((x1>x2)||(y1>y2)||(res.rem!=0)) MessageBox("Sélection de la zone\nou la largeur n'est pas multiple de 4","Erreur");
  else
  {                     
   *i1 = x1; *j1= y1;
   *i2 = x2; *j2= y2;
   LeParent->Invalidate(FALSE);        
   EndDialog(IDOK);
  }                                      
} 

void CDialZone::OnCancel( )
{
   LeParent->Invalidate(FALSE);
   EndDialog(IDOK);
}



// *********************************************************************************
// ****          choix d'un seuil de binarisation **********************************
// *********************************************************************************
BEGIN_MESSAGE_MAP(CDialSeuil, CDialog)     
    //{{AFX_MSG_MAP(CDialSeuil)
   	ON_COMMAND(IDOK, OnOK) 
   	ON_COMMAND(IDCANCEL, OnCancel)
   	ON_COMMAND(IDC_CHECK,OnCheck) 
   	//}}AFX_MSG_MAP
END_MESSAGE_MAP()                            


CDialSeuil::CDialSeuil(unsigned char *initval,BOOL *initcheck): CDialog(IDD_SEUIL,NULL)
{     
	val=initval;    
	check = initcheck;
	*check = TRUE;
}                         
                     
                     
void CDialSeuil::OnCheck()
{
  if (*check) *check = FALSE;else *check = TRUE;
}
                     
BOOL CDialSeuil::OnInitDialog()
{
 char ch1[10];
 sprintf(ch1,"%d",*val);  
 SetDlgItemText(IDC_SEUIL,(LPSTR) ch1);
 CButton * p = ( CButton *) GetDlgItem(IDC_CHECK); 
 p->SetCheck(1);    
 return CDialog::OnInitDialog(); 
}

void CDialSeuil::OnOK( )
{ 
  int x1;
  char ch1[10];                             
  GetDlgItemText(IDC_SEUIL,(LPSTR) ch1,sizeof(ch1)); 
  x1 = (int) atoi(ch1);
  if ((x1<0)||(x1>255)) MessageBox("Erreur","Choix du seuil");
  else
  {                     
   *val = (unsigned char) x1; 
   EndDialog(IDOK);
  } 
} 


// *********************************************************************************
// ****          choix des limites de la sequence **********************************
// *********************************************************************************
BEGIN_MESSAGE_MAP(CDialSeqlimite, CDialog)     
   	//{{AFX_MSG_MAP(CDialSeqlimite)
   	ON_COMMAND(IDOK, OnOK) 
   	ON_COMMAND(IDCANCEL, OnCancel)
   	//}}AFX_MSG_MAP
END_MESSAGE_MAP()                            


CDialSeqlimite::CDialSeqlimite(int *Initi1,int *Initi2): CDialog(IDD_LIMITE_SEQUENCE,NULL)
{     
	i1 = Initi1;
	i2 = Initi2; 
	depas = 256;
}

BOOL CDialSeqlimite::OnInitDialog()
{
 char ch1[10],ch2[10];
 sprintf(ch1,"%d",*i1);  
 SetDlgItemText(IDC_NUM1,(LPSTR) ch1); 
 sprintf(ch2,"%d",*i2);  
 SetDlgItemText(IDC_NUM2,(LPSTR) ch2);   
 return CDialog::OnInitDialog(); 
}

void CDialSeqlimite::OnOK( )
{ 
  int x1,x2;
  char ch1[10],ch2[10];                             
  GetDlgItemText(IDC_NUM1,(LPSTR) ch1,sizeof(ch1)); 
  x1 = (int) atoi(ch1); 
  GetDlgItemText(IDC_NUM2,(LPSTR) ch2,sizeof(ch2)); 
  x2 = (int) atoi(ch2);
  if ((x1<0)||(x1>x2)||(x2>depas)) MessageBox("Erreur","Choix des limites");
  else
  {                     
   *i1 = (int) x1; 
   *i2 = (int) x2; 
   EndDialog(IDOK);       
  } 
}                              


// *********************************************************************************
// ****          choix des caractéristique de la sequence **********************************
// *********************************************************************************
BEGIN_MESSAGE_MAP(CDialSeqCarac, CDialog)     
	//{{AFX_MSG_MAP(CDialSeqCarac)
   	ON_COMMAND(IDOK, OnOK) 
   	ON_COMMAND(IDCANCEL, OnCancel)
   	//}}AFX_MSG_MAP
END_MESSAGE_MAP()                            


CDialSeqCarac::CDialSeqCarac(int *Initi1,unsigned char *Initi2): CDialog(IDD_CARAC_ZONE,NULL)
{     
	i1 = Initi1;
	i2 = Initi2; 
}

BOOL CDialSeqCarac::OnInitDialog()
{
 char ch1[10],ch2[10];
 sprintf(ch1,"%d",*i1);  
 SetDlgItemText(IDC_NBIMAGE,(LPSTR) ch1); 
 sprintf(ch2,"%d",*i2);  
 SetDlgItemText(IDC_TEMPS,(LPSTR) ch2);   
 return CDialog::OnInitDialog(); 
}

void CDialSeqCarac::OnOK( )
{ 
  int x1,x2;
  char ch1[10],ch2[10];                             
  GetDlgItemText(IDC_NBIMAGE,(LPSTR) ch1,sizeof(ch1)); 
  x1 = (int) atoi(ch1); 
  GetDlgItemText(IDC_TEMPS,(LPSTR) ch2,sizeof(ch2)); 
  x2 = (int) atoi(ch2);
  if ((x1<0)) MessageBox("Erreur","Choix des caractéristiques");
  else
  {                     
   *i1 = (int) x1; 
   *i2 = (unsigned char) x2; 
   EndDialog(IDOK);
  } 
} 

// *********************************************************************************
// ****          choix de deux seuils  de binarisation **********************************
// *********************************************************************************
BEGIN_MESSAGE_MAP(CDialDoubleSeuil, CDialog)     
   	//{{AFX_MSG_MAP(CDialDoubleSeuil)
   	ON_COMMAND(IDOK, OnOK) 
   	ON_COMMAND(IDCANCEL, OnCancel)
   	//}}AFX_MSG_MAP
END_MESSAGE_MAP()                            


CDialDoubleSeuil::CDialDoubleSeuil(unsigned int IDDIAL,unsigned char *initval1, unsigned char *initval2): CDialog(IDDIAL,NULL)
{     
	val1=initval1; 
	val2=initval2;
}

BOOL CDialDoubleSeuil::OnInitDialog()
{
 char ch1[10];
 char ch2[10];
 sprintf(ch1,"%d",*val1); 
 sprintf(ch2,"%d",*val2); 
 SetDlgItemText(IDC_SEUIL1,(LPSTR) ch1);
 SetDlgItemText(IDC_SEUIL2,(LPSTR) ch2);   
 return CDialog::OnInitDialog(); 
}

void CDialDoubleSeuil::OnOK( )
{ 
  int x1, x2;
  char ch1[10], ch2[10]; 
                              
  GetDlgItemText(IDC_SEUIL1,(LPSTR) ch1,sizeof(ch1)); 
  x1 = (int) atoi(ch1); 
  GetDlgItemText(IDC_SEUIL2,(LPSTR) ch2,sizeof(ch2)); 
  x2 = (int) atoi(ch2);
  if((x1<0)||(x2>255)) MessageBox("Erreur","Choix des seuils");
  else
  {                     
   *val1 = (unsigned char) x1; 
   *val2 = (unsigned char) x2;
   EndDialog(IDOK);
  } 
} 




//********************************************