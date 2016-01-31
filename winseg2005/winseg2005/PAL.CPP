#include "resource.h"   
#include "stdafx.h"
#include "pal.h"        
#include "message.h"  
 #pragma warning( disable : 4996 )       
 
BEGIN_MESSAGE_MAP(CPaletteWnd, CWnd)      
    //{{AFX_MSG_MAP(CPaletteWnd)
    ON_WM_LBUTTONDOWN()     
    ON_WM_LBUTTONUP()           
    ON_WM_MOUSEMOVE()          
	ON_WM_KEYDOWN()
    ON_WM_PAINT()    
    ON_MESSAGE(WM_TEXTE1,OnMemorise)       
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONUP()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()


BOOL CPaletteWnd::Create(LPCSTR szTitle, LONG style ,CWnd * parent)
{ 
	BOOL LeBool =CWnd::CreateEx(/*style */WS_EX_PALETTEWINDOW , AfxRegisterWndClass( CS_HREDRAW | CS_VREDRAW ,NULL,(HBRUSH)(COLOR_BACKGROUND),NULL),
	szTitle, style,100,100,300, 125, parent->m_hWnd, NULL);
	return LeBool;  					
}
 

CPaletteWnd::CPaletteWnd(CPalette * InitMyPalette,CWnd *  pointeur)
{   
	pImage = pointeur;
	DWORD Taille;   
	Taille  = (DWORD)256*(DWORD)20;
	id = GlobalAlloc(GMEM_MOVEABLE,Taille);            
	Image =  (unsigned char  *) GlobalLock(id);    
	Appuie = FALSE;
	Appuie2 = FALSE;	
	Seuil = 0;
	Seuil2 = 255;
	int i,j;
	for (j=0;j<256;j++) for (i=0;i<20;i++) Image[i*256+j] = (unsigned  char) j;
    
    // création du DIB                           
	MyInitInfo = (BITMAPINFO*) new unsigned char[ sizeof(BITMAPINFOHEADER) + sizeof(RGBQUAD)*256];
	MyInfoHeader = (BITMAPINFOHEADER*) new unsigned char[ sizeof(BITMAPINFOHEADER)];
	MyInfoHeader-> biSize = sizeof(BITMAPINFOHEADER);
	MyInfoHeader-> biPlanes = 1;
	MyInfoHeader-> biBitCount = 8;                    
	MyInfoHeader-> biWidth = 256;
	MyInfoHeader-> biHeight = 20;
	MyInfoHeader-> biCompression = BI_RGB;
	MyInfoHeader-> biSizeImage = 0;
	MyInfoHeader-> biXPelsPerMeter = 0;
	MyInfoHeader-> biYPelsPerMeter = 0;
	MyInfoHeader-> biClrUsed = 0;
	MyInfoHeader-> biClrImportant = 0;
	MyInitInfo->  bmiHeader = *MyInfoHeader; 
	
	// Création de la palette 
	MyLogPalette =(LOGPALETTE*) new unsigned char[sizeof(LOGPALETTE) + sizeof(PALETTEENTRY)*256 ];
	MyLogPalette->palVersion = 0x300;
	MyLogPalette->palNumEntries = 256;
	MyPalette = InitMyPalette; 
		
}                         

LONG CPaletteWnd::OnMemorise(UINT uiParam,LONG lParam)
{
	MyPalette->GetPaletteEntries(0, 256, MyLogPalette->palPalEntry );    
	Invalidate(TRUE);   
	return 0;
}

CPaletteWnd::~CPaletteWnd()
{   
	delete MyInitInfo;                                  
	delete MyInfoHeader;    
	delete MyLogPalette;   
	// libération de la mémoire de l'image palette
	GlobalUnlock(id);
	free(id); 
                        
}
 
afx_msg void CPaletteWnd::OnPaint(void)   
{       
	CPaintDC dc(this);                                        
	Mapping_Palette();
	dc.Rectangle(19,19,21+256,21+20);
	dc.Rectangle(53+10,44+10,205+20,65+20);
	SetDIBitsToDevice(dc.GetSafeHdc(),20,20, (int) 256,(int) 20,0,0,0,(int) 20,Image,MyInitInfo,DIB_RGB_COLORS);
	dc.SetBkMode(TRANSPARENT);
	dc.TextOut(17,46,"0");   dc.TextOut(8+255,46,"255"); 
	dc.SetBkMode(OPAQUE);
	dc.TextOut(55+10,46+20,"                      "); 
	if (Seuil!=0 ||Seuil2!=0) 
	{
		char tabSeuil[25];
		sprintf(tabSeuil,"Seuil1=%3d Seuil2=%3d",Seuil,Seuil2);     
		dc.TextOut(55+15,46+20,tabSeuil);      
	}
	ReleaseDC(&dc); // le dc est relaché          WS_
} 

void CPaletteWnd::Mapping_Palette()
{   
	MyPalette->GetPaletteEntries(0, 256,   MyLogPalette->palPalEntry ); 
	for (int i = 0; i<256; i++)
	{
		MyInitInfo->bmiColors[i].rgbBlue     =   MyLogPalette->palPalEntry[i].peBlue;
		MyInitInfo->bmiColors[i].rgbGreen    =   MyLogPalette->palPalEntry[i].peGreen;
		MyInitInfo->bmiColors[i].rgbRed      =   MyLogPalette->palPalEntry[i].peRed;
		MyInitInfo->bmiColors[i].rgbReserved =   MyLogPalette->palPalEntry[i].peFlags;
	}	
}	



// ****************************************************************************
// 							gestion de la souris                         
// ****************************************************************************
afx_msg void CPaletteWnd::OnLButtonDown(UINT nFlags,CPoint Point) 
{ 
		if ((Point.x>=20)&&(Point.x <Seuil2+20)&&(Point.y>20 )&&(Point.y<40)) 
		{ 
			Seuil = Point.x-20; 
			Appuie = TRUE;
			Modifie_palette();    
			pImage->Invalidate(FALSE); 
			Invalidate(FALSE);
		} 
}      
void CPaletteWnd::OnRButtonDown(UINT nFlags, CPoint point) 
{
	
	if ((point.x>=20+Seuil)&&(point.x <20+256)&&(point.y>20 )&&(point.y<40)) 
	{ 
		Seuil2 = point.x-20; 
		Appuie2 = TRUE;
		Modifie_palette();    
		pImage->Invalidate(FALSE); 
		Invalidate(FALSE);
	} 

	CWnd::OnRButtonDown(nFlags, point);
}

void CPaletteWnd::OnRButtonUp(UINT nFlags, CPoint point) 
{
	Appuie2 = FALSE;

	CWnd::OnRButtonUp(nFlags, point);
}

   
afx_msg void CPaletteWnd::OnLButtonUp(UINT nFlags,CPoint Point)
{ 
	Appuie = FALSE;
}                  

afx_msg void CPaletteWnd::OnMouseMove(UINT nFlags,CPoint Point)
{
	if((Point.x>=20)&&(Point.x <20+256)&&(Point.y>20 )&&(Point.y<40))
			SetCursor(LoadCursor(NULL,IDC_SIZEWE));
	else	SetCursor(LoadCursor(NULL,IDC_ARROW));

	if (Appuie)
	{    
		if ((Point.x>=20)&&(Point.x <20+256)&&(Point.y>20 )&&(Point.y<40)) 
		{
			Seuil = Point.x-20;  
			Modifie_palette();    
			pImage->Invalidate(FALSE); 
			Invalidate(FALSE);
		} 
	}	  
	if (Appuie2)
	{    
		if ((Point.x>=Seuil+20)&&(Point.x <20+256)&&(Point.y>20 )&&(Point.y<40)) 
		{
			Seuil2 = Point.x-20;  
			Modifie_palette();    
			pImage->Invalidate(FALSE); 
			Invalidate(FALSE);
		} 
	}	  
}

    
    
void CPaletteWnd::Modifie_Color(int i)
{                           
	if (i>Seuil2)
	{
		MyLogPalette->palPalEntry[i].peBlue = 0  ; 
		MyLogPalette->palPalEntry[i].peGreen = 0;  
		MyLogPalette->palPalEntry[i].peRed = 255;
	}  
	else if(i>Seuil)
	{
		MyLogPalette->palPalEntry[i].peBlue = 255  ; 
		MyLogPalette->palPalEntry[i].peGreen = 255;  
		MyLogPalette->palPalEntry[i].peRed = 255;
	}
	else
	{
		MyLogPalette->palPalEntry[i].peBlue = 0  ; 
		MyLogPalette->palPalEntry[i].peGreen = 0;  
		MyLogPalette->palPalEntry[i].peRed = 0;
	}
}
    
void CPaletteWnd::Modifie_palette()
{
	for (int i = 0; i<256; i++)
	{
		Modifie_Color(i);    
		MyInitInfo->bmiColors[i].rgbBlue = MyLogPalette->palPalEntry[i].peBlue ;
		MyInitInfo->bmiColors[i].rgbGreen = MyLogPalette->palPalEntry[i].peGreen;    
		MyInitInfo->bmiColors[i].rgbRed = MyLogPalette->palPalEntry[i].peRed ;
	}	
	MyPalette->SetPaletteEntries(0, 256,MyLogPalette->palPalEntry ); 
}    


// ************** FIN DE GESTION DE LA SOURIS *******************************************
// **************************************************************************************


void CPaletteWnd::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	
	if ((nChar == 37) || (nChar == 39) ) 
	{
		if(nChar==37)  Seuil--; //Fleche à gauche <-                    
		if(nChar==39)  Seuil++; //Fleche à droite ->
				
		if (Seuil <0) Seuil = 0;
		if (Seuil > 255) Seuil = 255;
		Modifie_palette();    
		pImage->Invalidate(FALSE); 
		Invalidate(FALSE);                        
	}		
	CWnd::OnKeyDown(nChar, nRepCnt, nFlags);   
}



