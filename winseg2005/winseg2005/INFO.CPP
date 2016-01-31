 /////////////////////////////////////////////////////////////////////////////   
#include "stdafx.h" 
#include "info.h"
#include "resource.h"
#include "message.h"
       
#pragma warning( disable : 4996 )


BEGIN_MESSAGE_MAP(CInfoWnd, CWnd)
	//{{AFX_MSG_MAP(CInfoWnd)
	// fonctions déclarées dans CImage	
	ON_MESSAGE(WM_TEXTE1,OnAfficheTexte1)
	ON_MESSAGE(WM_TEXTE2,OnAfficheTexte2)     
	ON_MESSAGE(WM_TEXTE3,OnAfficheTexte3)  
	ON_MESSAGE(WM_TEXTE4,OnAfficheTexte4)  
    ON_WM_PAINT()           
    //}}AFX_MSG_MAP
END_MESSAGE_MAP()                            

BOOL CInfoWnd::Create(CRect rect, CWnd * pParent)
{   
   return CWnd::Create(AfxRegisterWndClass( CS_HREDRAW | CS_VREDRAW ,LoadCursor(NULL, IDC_ARROW ),(HBRUSH)(16),NULL)
   		,"",WS_VISIBLE,rect,pParent,NULL);
}   
 
LONG CInfoWnd::OnAfficheTexte1(UINT uiParam,LONG lParam)
{
  TEXTE = TRUE;
  sprintf(chaine1,"%s ",(LPSTR)lParam);
  Invalidate(FALSE);
  return 0;
}    

 
LONG CInfoWnd::OnAfficheTexte2(UINT uiParam,LONG lParam)
{
  TEXTE = TRUE;   
  sprintf(chaine2,"%s   ",(LPSTR)lParam);
  Invalidate(FALSE);
  return 0;
}             

LONG CInfoWnd::OnAfficheTexte3(UINT uiParam,LONG lParam)
{
  TEXTE = TRUE;   
  sprintf(chaine3,"%s   ",(LPSTR)lParam);
  Invalidate(FALSE);
  return 0;
}    

LONG CInfoWnd::OnAfficheTexte4(UINT uiParam,LONG lParam)
{
  TEXTE = TRUE;   
  sprintf(chaine4,"%s   ",(LPSTR)lParam);
  Invalidate(FALSE);
  return 0;
}    


CInfoWnd::CInfoWnd()
{
   sprintf(chaine1,"");  
   sprintf(chaine2,"");
   sprintf(chaine3,"");
   sprintf(chaine4,"");
   TEXTE = FALSE;
}


afx_msg void CInfoWnd::OnPaint()
{  
//    char ch[200];
	int taille1, taille2, taille3;
	CPaintDC dc(this);    
	CClientDC DC(this); 
	CFont ft1,ft2,ft3,ft4;
	CRect rcClient;
	GetClientRect(&rcClient);
	CBrush bFond;
	
	bFond.CreateSolidBrush(RGB(180,180,180));
	dc.SelectObject(&bFond);
	dc.FillRect(&rcClient,&bFond);
//	DeleteObject(&bFond);
	if (TEXTE)
	{	
/*	    strcpy(ch,chaine1);
	    strcat(ch,chaine2); 
	    strcat(ch,chaine3);
	    strcat(ch,chaine4);*/
	    taille1 = strlen(chaine1);
	    taille2 = strlen(chaine2);
	    taille3 = strlen(chaine3);
		DC.SetBkMode(TRANSPARENT);
 		ft1.CreateFont(11,5,0,0,FW_NORMAL,0,0,0,0,OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH,NULL);
  		DC.SelectObject(&ft1);
  		DC.TextOut(5,5,chaine1);
  		ft1.DeleteObject();  

		ft2.CreateFont(12,6,0,0,FW_BOLD,0,0,0,0,OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH,NULL);
  		DC.SelectObject(&ft2);
  		DC.TextOut(5+taille1*5,4,chaine2);
  		ft2.DeleteObject();  

		ft3.CreateFont(14,7,0,0,FW_SEMIBOLD,0,0,0,0,OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH,NULL);
  		DC.SelectObject(&ft3);
  		DC.TextOut(5+taille1*5+taille2*6,4,chaine3);
  		ft3.DeleteObject();  

 		ft4.CreateFont(16,8,0,0,FW_DEMIBOLD,0,0,0,0,OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,DEFAULT_QUALITY,DEFAULT_PITCH,NULL);
  		DC.SelectObject(&ft4);
  		DC.TextOut(5+taille1*5+taille2*6+taille3*7,4,chaine4);
  		ft4.DeleteObject();  
	}	              
	ReleaseDC(&DC);
	ReleaseDC(&dc);   
} 