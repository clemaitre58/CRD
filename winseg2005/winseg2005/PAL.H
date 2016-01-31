      
#ifndef _PAL_H__     
#define _PAL_H__     
                      
class CPaletteWnd : public CWnd
{   
	public:
	unsigned char  * Image;
	LOGPALETTE		* 	MyLogPalette; 
	    
	private: 
	//pour l'affichage du Bitmap 
	BITMAPINFO	 			*	MyInitInfo;
	BITMAPINFOHEADER	*	MyInfoHeader;   
	HGLOBAL id;
	
	CPalette 		*	MyPalette; 
	BOOL			Appuie,Appuie2;    
	CEdit 			*	EditSeuil;
	unsigned char	Seuil,Seuil2;
	CRect				rcClient;
	CButton 		*   MyButSeuiller,*MyButFermer; 
	
	CWnd * 	pImage; 	
	public:
	// constructeur
	CPaletteWnd(CPalette * InitMyPalette,CWnd *  pointeur);      
	// destructeur
	~CPaletteWnd();   
		
	private:
	void Mapping_Palette();
	void Modifie_Color(int i);
	void Modifie_palette();
		
	public:
	LONG OnMemorise(UINT ,LONG );
	    
	protected:
	DECLARE_MESSAGE_MAP(); 

		          
	protected: 
	// gestion des événements souris
    //{{AFX_MSG(CPaletteWnd)
	afx_msg void OnLButtonDown(UINT nFlags,CPoint Point);  
	afx_msg void OnLButtonUp(UINT nFlags,CPoint Point);
	afx_msg void OnMouseMove(UINT nFlags,CPoint Point); 
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	//}}AFX_MSG
	  	
	  
		
	public:    
	BOOL Create(LPCSTR szTitle, LONG style ,CWnd * parent);   
//	BOOL Create(LPCSTR szTitle, LONG style = 0,CMDIFrameWnd* pParent = NULL);
		
	// reponses au message du menu
	afx_msg void OnPaint(void);  
};  

#endif // __PAL_H__