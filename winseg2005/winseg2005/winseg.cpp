// 05/01/98   
// version 8.0

#include "stdafx.h"
#include "resource.h"
// include du projet
#include "winseg.h" 
#include "mainfrm.h"
#include "Splash.h"
#include <process.h>
                       

/////////////////////////////////////////////////////////////////////////////
// CMdiApp
BEGIN_MESSAGE_MAP(CMdiApp, CWinApp)  
//{{AFX_MSG_MAP(CMdiApp)
//}}AFX_MSG_MAP
END_MESSAGE_MAP()                  

CMdiApp  theApp;                          


CMdiApp::~CMdiApp()
{  
}
   
                           
BOOL CMdiApp::InitInstance()
{   

	///donne la classe temps réel au process

/*	HANDLE myProcess = OpenProcess(PROCESS_ALL_ACCESS,TRUE,_getpid());

	if(!SetPriorityClass(myProcess,REALTIME_PRIORITY_CLASS))
		MessageBox(NULL,"Impossible de modifier la priorité du process !","pb2",MB_OK); 

	if(!SetThreadPriority( THREAD_PRIORITY_TIME_CRITICAL)) 
		MessageBox(NULL,"Impossible de modifier la priorité du thread principal","pb",MB_OK); 

	if((GetThreadPriority()!=THREAD_PRIORITY_TIME_CRITICAL)||(GetPriorityClass(myProcess)!=REALTIME_PRIORITY_CLASS)) MessageBox(NULL,"PB de priorité","Erreur",MB_OK);
	
	CloseHandle(myProcess);
*/	// CG: This line was added by the ActiveX Control Containment component
	
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	// À définir pour inclure toutes les classes de contrôles communs à utiliser
	// dans votre application.
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinApp::InitInstance();

	// Initialiser les bibliothèques OLE
	if (!AfxOleInit())
	{
//		AfxMessageBox(IDP_OLE_INIT_FAILED);
		return FALSE;
	}
	AfxEnableControlContainer();
	
	
	
	AfxEnableControlContainer();

	// CG: The following block was added by the Splash Screen component.
	{
		CCommandLineInfo cmdInfo;
		ParseCommandLine(cmdInfo);

		CSplashWnd::EnableSplashScreen(cmdInfo.m_bShowSplash);
	}

	LoadStdProfileSettings();  // Load standard INI file options (including MRU)

	// create main MDI Frame window
	CMainFrame * pMainFrame = new CMainFrame;  
	if (!pMainFrame->Create(AfxRegisterWndClass(NULL,0,0,LoadIcon(IDI_LE2I)), "Le2i-Winseg 12", 
		WS_OVERLAPPEDWINDOW | FWS_ADDTOTITLE|WS_SYSMENU,CFrameWnd::rectDefault,NULL,
		MAKEINTRESOURCE(IDR_MAINFRAME),  0 , NULL )) return FALSE;
//	if (!pMainFrame->LoadFrame(IDR_MAINFRAME)) return FALSE;
	pMainFrame->ShowWindow(SW_HIDE);
	pMainFrame->LoadAccelTable(MAKEINTRESOURCE(IDR_MAINFRAME));
	m_pMainWnd = pMainFrame;
	return TRUE;
}

BOOL CMdiApp::PreTranslateMessage(MSG* pMsg)
{
	// CG: The following lines were added by the Splash Screen component.
	if (CSplashWnd::PreTranslateAppMessage(pMsg))
		return TRUE;

	return CWinApp::PreTranslateMessage(pMsg);
}
