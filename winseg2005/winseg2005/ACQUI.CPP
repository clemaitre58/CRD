 
#include "stdafx.h"
#include "mainfrm.h"
#include "acqui.h"  
//#include "carte.h"
#include "dialogue.h"
#include "ProfilWnd.h"

// ouverture d'une fenêtre d'acquisition  
BOOL AcquisitionImage = FALSE;

BEGIN_MESSAGE_MAP(CAcqImageWnd,CTraitImageWnd)      
   	//{{AFX_MSG_MAP(CAcqImageWnd)
	ON_WM_PAINT()   
	ON_COMMAND(ID_SNAP, OnSnap)   
	ON_COMMAND(ID_GRAB, OnGrab)    
	ON_COMMAND(ID_FREEZE, OnFreeze) 
	ON_COMMAND(ID_CREER_SEQUENCE, OnCreerSequence )   
	ON_WM_TIMER()
	ON_UPDATE_COMMAND_UI(ID_GRAB, OnUpdateGrab)
	ON_UPDATE_COMMAND_UI(ID_SNAP, OnUpdateSnap)
	ON_UPDATE_COMMAND_UI(ID_FREEZE, OnUpdateFreeze)
	ON_UPDATE_COMMAND_UI(ID_CREER_SEQUENCE, OnUpdateCreerSequence)
	ON_COMMAND(ID_ACQUI_SEQUENCE, OnAcquerirSequence)
	ON_UPDATE_COMMAND_UI(ID_ACQUI_SEQUENCE, OnUpdateSequence)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()                            


BOOL CAcqImageWnd::Create(LPCSTR szTitle,LONG style,CMDIFrameWnd* pParent)
{ 
	BOOL b = CImageWnd::Create(szTitle,style,pParent);
  	return b;
}


CAcqImageWnd::CAcqImageWnd(int InitType, NBCOLOR initColor = NB)
	:CTraitImageWnd(ReturnLarg(),ReturnHaut(),initColor)
{
 	// initialisation de constante de la sequence 
	TypeAcq = InitType;  
	Init();
} 


CAcqImageWnd::~CAcqImageWnd()
{
	// destruction de la sequence si pas detruite
	AcquisitionImage = FALSE;
}  


void CAcqImageWnd::Init()
{
 	Temps  = 0; 
 	LeTimerStop = TRUE;	
 	NbImage = 10;
}
 
void CAcqImageWnd::ChargeMenu()
{   
	if (menu.m_hMenu==NULL) 
	{
		if (TypeColor == COLOR)  menu.LoadMenu(IDR_ACQIMAGECOLOR);
		else menu.LoadMenu(IDR_ACQIMAGENB);
	}
	m_hMenuShared = menu.m_hMenu;  
	AcquisitionImage = TRUE;       
}	

afx_msg void CAcqImageWnd::OnPaint(void)
{

	if (!LeTimerStop)
	{
		Acquerir(2,(unsigned char  *) Image,(int)0,(int)0,(int)Larg-1,(int)Haut-1);
	}	

 CImageWnd::OnPaint();
} 


// ****************************************
// Routines liées à l'affichage en continu 
// ****************************************
void CAcqImageWnd::MettreLeContinu()
{                                                                                 
 	if (LeTimerStop)
	{ 
		TypeAcq=IMAGE_CONTINUE;
		LeTimerStop = FALSE;
		SetTimer(ID_TIMER,0,NULL);   // branchement à OnTimer toutes les 0 ms en théorie !!
 	}	
	((CMainFrame*) AfxGetApp()->m_pMainWnd)->bar.OnUpdateCmdUI( (CFrameWnd*) GetMDIFrame(), TRUE);
}
     
afx_msg void CAcqImageWnd::OnTimer(UINT  nIDEvent)
{
	 CImageWnd::OnTimer(nIDEvent);  
	 if (IsWindow(pCoupeVariable->GetSafeHwnd())) MettreAJourLaCoupeActive();
	 Invalidate(FALSE);   // envoi d'un message WM_PAINT
}                  

      
void CAcqImageWnd::ArreteLeContinu() 
{ 
	if (!LeTimerStop)
	{   
		TypeAcq=IMAGE_FIXE;
	/*	menu.CheckMenuItem(ID_SNAP,MF_CHECKED);
		menu.CheckMenuItem(ID_FREEZE,MF_UNCHECKED);
		menu.CheckMenuItem(ID_GRAB,MF_CHECKED);
 	*/	LeTimerStop = TRUE; 
 		KillTimer(ID_TIMER);                              
	} 
}
 


// *****************************************************
// Routine de Snap
// ****************************************************
afx_msg void CAcqImageWnd::OnSnap()               
{  
	ArreteLeContinu();      
/*	menu.CheckMenuItem(ID_SNAP,MF_CHECKED);
	menu.CheckMenuItem(ID_FREEZE,MF_UNCHECKED);
	menu.CheckMenuItem(ID_GRAB,MF_UNCHECKED);    
*/	
	Acquerir(0,(unsigned char  *) Image,(int)0,(int)0,(int)Larg-1,(int)Haut-1);
	Invalidate(FALSE);
	
	Memorise();
}

// *****************************************************
// Routine de Freeze
// ****************************************************
afx_msg void CAcqImageWnd::OnFreeze()               
{  
	ArreteLeContinu();  
/*	menu.CheckMenuItem(ID_SNAP,MF_UNCHECKED);
	menu.CheckMenuItem(ID_FREEZE,MF_CHECKED);
	menu.CheckMenuItem(ID_GRAB,MF_UNCHECKED); 
*/	 
	Acquerir(1,(unsigned char   *) Image,(int)0,(int)0,(int)Larg-1,(int)Haut-1);
	Invalidate(FALSE);
	Memorise();
}


// *****************************************************
// Routine Grab
// ****************************************************
afx_msg void CAcqImageWnd::OnGrab()               
{
	if (LeTimerStop)
	{ 
		MettreLeContinu();
	/*	menu.CheckMenuItem(ID_SNAP,MF_UNCHECKED);
		menu.CheckMenuItem(ID_FREEZE,MF_UNCHECKED);
		menu.CheckMenuItem(ID_GRAB,MF_CHECKED);        
	*/	Invalidate(TRUE);
	} 
}
        

// **********************************************************************************
// *************** GESTION DE LA SEQUENCE 
// **********************************************************************************
afx_msg  void CAcqImageWnd::OnAcquerirSequence()
{
	if (SEQUENCE_CREE)
	{     
		ArreteLeContinu();
	/*	menu.CheckMenuItem(ID_SNAP,MF_CHECKED);
		menu.CheckMenuItem(ID_FREEZE,MF_UNCHECKED);
		menu.CheckMenuItem(ID_GRAB,MF_UNCHECKED);
	*/	           
		SetCursor(LoadCursor(0,IDC_WAIT));	   
		DWORD dt = SequenceCarte(NbMaxImageSeq,Bpile,Temps,(int)debut_i,(int)debut_j,(int)fin_i,(int)fin_j);  
		SetCursor(LoadCursor(0,IDC_ARROW));
		char chaine[50];
		sprintf(chaine,"vitesse de transfert = %lu images/s",(DWORD)NbMaxImageSeq*1000L/dt);
		MessageBox(chaine,"* Résultat *");
		// souris normal
		SEQ_EXISTE=TRUE;
	}
	else MessageBox("Créer la séquence avant de l'acquérir","",MB_OK|MB_ICONEXCLAMATION); 
}

afx_msg  void CAcqImageWnd::OnCreerSequence()
{ 
	ArreteLeContinu();
	if (!SEQUENCE_CREE)
	{ 
		CDialSeqCarac dial(&NbImage,&Temps); 
		if (dial.DoModal()==IDOK)
			AlloueSequence(&NbImage);
	}
	else MessageBox("Détruire la séquence précédente\navant d'en créer une nouvelle","",MB_OK|MB_ICONEXCLAMATION);
}                                

	  
void CAcqImageWnd::OnUpdateGrab(CCmdUI* pCmdUI) 
{
	// TODO: Add your command update UI handler code here
	pCmdUI->SetCheck(TypeAcq==IMAGE_CONTINUE);
}

void CAcqImageWnd::OnUpdateSnap(CCmdUI* pCmdUI) 
{
	// TODO: Add your command update UI handler code here
//	pCmdUI->SetCheck(TypeAcq!=IMAGE_CONTINUE);
}

void CAcqImageWnd::OnUpdateFreeze(CCmdUI* pCmdUI) 
{
	// TODO: Add your command update UI handler code here
	//pCmdUI->SetCheck(TypeAcq==IMAGE_FIXE);
}

void CAcqImageWnd::OnUpdateSequence(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(SEQUENCE_CREE);	
}

void CAcqImageWnd::OnUpdateCreerSequence(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(!SEQUENCE_CREE);	
}

