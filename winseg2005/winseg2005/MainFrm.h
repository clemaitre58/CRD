#ifndef __MainFrame_H__           
#define __MainFrame_H__

#include <afxadv.h>

class CMainFrame : public CMDIFrameWnd
{
	
	DECLARE_DYNAMIC(CMainFrame)

	public:
	BOOL premier;

	int m_nDatePaneNo;
	int m_nTimePaneNo;
	CToolBar bar;


	/////////Parzen 3D
	byte EspClas[256*256*256];
	BOOL EspValid;
/*	unsigned __int32 EspProb1[256*256*256];
	unsigned __int32 EspProb2[256*256*256];
	unsigned __int32 EspProb3[256*256*256];
	unsigned __int32 EspProb4[256*256*256];
*/	public:
	CMainFrame(); 
	~CMainFrame(); 
	  
	public:
	CStatusBar  m_wndStatusBar;
	afx_msg void OnUpdateDate(CCmdUI* pCmdUI);
	afx_msg void OnUpdateTime(CCmdUI* pCmdUI);
	//{{AFX_MSG(CMainFrame)
	afx_msg void OnCloseAll();   
	afx_msg void OnNewImage();   
	afx_msg void OnOuvrirSequence();
	afx_msg void OnImagevide();
	afx_msg void OnAcquerirImage();
	afx_msg void OnAcquerirContinue();
	afx_msg void OnFileMruFile1();
	afx_msg void OnFileMruFile2();
	afx_msg void OnFileMruFile3();
	afx_msg void OnFileMruFile4();
	afx_msg void OnUpdateFileMruFile2(CCmdUI* pCmdUI);
	afx_msg void OnUpdateFileMruFile3(CCmdUI* pCmdUI);
	afx_msg void OnUpdateFileMruFile4(CCmdUI* pCmdUI);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnEditToolbar();
	afx_msg void OnUpdateEditToolbar(CCmdUI* pCmdUI);
	afx_msg void OnFileExit();
	afx_msg void OnClose();
	afx_msg void OnAbout();
	afx_msg void OnUpdateGrab(CCmdUI* pCmdUI);
	afx_msg void OnUpdateSnap(CCmdUI* pCmdUI);
	afx_msg void OnParzen3dCalculdesespacesdeproba();
	afx_msg void OnParzen3dCrerunespacededcision();
	afx_msg void OnParzen3dOuvrirunespacededcision();
	afx_msg void OnParzen3dSauvegarderunespacededcision();
	afx_msg void OnParzen3dVisualiserlespacededcision();
	afx_msg void OnUpdateParzen3dSauvegarderunespacededcision(CCmdUI* pCmdUI);
	afx_msg void OnUpdateParzen3dVisualiserlespacededcision(CCmdUI* pCmdUI);
	afx_msg void OnParzen3dMesrurerecouvrementsaveccl0();
	afx_msg void OnParzen3dMesuredevolumedesclasses();
	afx_msg void OnUpdateParzen3dMesuredevolumedesclasses(CCmdUI* pCmdUI);
	afx_msg void OnTestOuvrir();
	afx_msg void OnZernikeBatchmoments();
	afx_msg void OnEtudecurvi1dGenerationfichierprofilPorte();
	afx_msg void OnEtudecurvi1dGenerationfichierprofilPorte32();
	afx_msg void OnEtudecurvi1dGenerationfichierprofilPorte16();
	afx_msg void OnEtudecurvi1dGenerationfichierprofilPorte12();
	afx_msg void OnEtudecurvi1dEtudeavecobjetEtudeporte32();
	afx_msg void OnEtudecurvi1dEtudeavecobjetEtudeporte();
	afx_msg void OnEtudecurvi1dEtudeavecobjetEtudeporte16();
	afx_msg void OnEtudecurvi1dEtudeavecobjetEtudeporte12();
	afx_msg void OnEtudecurvi1dEtudeavecobjetEtudeporte12f4();
	afx_msg void OnEtudecurvi1dGenerationfichierprofilBitexture48();
	afx_msg void OnEtudecurvi1dEtudeavecobjetEtudebitexture48();
	afx_msg void OnEtudecurvi1dEtudeavecobjetEtudebitexture48f4();
	afx_msg void OnEtudecurvi1dGenerationfichierprofilPortebordnonsim48();
	afx_msg void OnEtudecurvi1dEtudeavecobjetEtudeportenonsim48f16();
	afx_msg void OnEtudecurvi1dGenerationfichierprofilMultisignaux1d();
	afx_msg void OnEtudecurvi1dEtudeavecobjetEtudesignauxmultiple48f16();
	afx_msg void OnEtudecurvi1dEtudeavecobjetEtudesignauxrels();
	afx_msg void OnEtudecurvi1dEtudeavecobjetTestnew();
	afx_msg void OnEtudecurvi1dEtudeavecobjetEtudemultifichier();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
private:
	BOOL InitStatusBar(UINT *pIndicators, int nSize, int nSeconds);
public:
	afx_msg void OnStegerConvolution1d();
	afx_msg void OnConversionSignal();
	afx_msg void OnOutilsAjoutbruit();
};
    
      
      


#endif // CMainFrame
