#ifndef __ACQUI_H__     
#define __ACQUI_H__     

#include "timage.h" 

          
class CAcqImageWnd : public CTraitImageWnd 
{   
	protected:
	DECLARE_MESSAGE_MAP();
	 
	 // variable du timer
	public:
	BOOL LeTimerStop;
	int TypeAcq;
	    
	unsigned char Temps;   
	
	// constructeur
	CAcqImageWnd(int, NBCOLOR);	
	void Init();
	
	// destructeur
	~CAcqImageWnd();                               
		
	BOOL Create(LPCSTR szTitle,LONG style = 0,CMDIFrameWnd* pParent = NULL);	 


	//{{AFX_MSG(CAcqImageWnd) 	
	afx_msg void OnPaint(void);   
	afx_msg void OnSnap();               
 	afx_msg void OnGrab() ; 
	afx_msg void OnFreeze();
	afx_msg void OnCreerSequence();   
	afx_msg void OnTimer(UINT  nIDEvent); 
	afx_msg void OnUpdateGrab(CCmdUI* pCmdUI);
	afx_msg void OnUpdateSnap(CCmdUI* pCmdUI);
	afx_msg void OnUpdateFreeze(CCmdUI* pCmdUI);
	afx_msg void OnUpdateCreerSequence(CCmdUI* pCmdUI);
	afx_msg void OnAcquerirSequence();
	afx_msg void OnUpdateSequence(CCmdUI* pCmdUI);
	//}}AFX_MSG
		
	virtual void ChargeMenu();
	


	// fonctions du timer
	void MettreLeContinu();
	void ArreteLeContinu();    
	    
	protected:
	LONG OnSetCoupeVariable(UINT uiParam,LONG lparam );
		
	// acquisition de l'image de la carte dans Image
};

#endif // __ACQUI_H__
