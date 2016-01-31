#ifndef __INFO_H__     
#define __INFO_H__     


class CInfoWnd : public CWnd
{   
 	public:
 	char chaine1[200];
 	char chaine2[200];
 	char chaine3[200]; 
 	char chaine4[200];
 	 
 	BOOL TEXTE;
    protected:
	DECLARE_MESSAGE_MAP();
	
    public:    
    LONG OnAfficheTexte1(UINT uiParam,LONG lParam);    
    LONG OnAfficheTexte2(UINT uiParam,LONG lParam);
    LONG OnAfficheTexte3(UINT uiParam,LONG lParam);
    LONG OnAfficheTexte4(UINT uiParam,LONG lParam);
	// constructeur
	CInfoWnd();     
	
	BOOL Create(CRect rect,CWnd * pParent);   
	//{{AFX_MSG(CInfoWnd)
	afx_msg void OnPaint();        
	//}}AFX_MSG
	
};  

#endif // __INFO_H__