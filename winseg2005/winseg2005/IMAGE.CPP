#include "stdafx.h"
#include "ProfilWnd.h"
#include "GraphWnd.h"
#include "bitmap.h"
#include "dialogue.h"  
#include "image.h"            
#include "timage.h"
#include "dialeche.h"
#include "DialFactZoom.h"
#include "dialogex.h"
#include "DialLecSeq.h"
#include "dialmorpho.h"
#include "dialrelief.h"
#include "mainfrm.h"
#include "petrou.h"
#include "dialcurvi.h"
#include "curvidetect.h"
#include "regress.h"

#include "FenetreProgressCurvi.h"
#include "curvlinear.h"
#include <conio.h>

//int NbCoupe = 0;
ldiv_t res;
#pragma warning( disable : 4996 )

IMPLEMENT_DYNAMIC(CImageWnd, CMDIChildWnd)

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                           
// NE PAS TOUCHER A LA CLASSE CImageWnd SANS AUTORISATION !!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
const UINT uiMsg3=::RegisterWindowMessage("WM_FINPROFIL");

BEGIN_MESSAGE_MAP(CImageWnd, CMDIChildWnd)
   //{{AFX_MSG_MAP(CImageWnd)
   ON_COMMAND(ID_FILE_MRU_FILE1, OnFileMruFile1)
   ON_COMMAND(ID_FILE_MRU_FILE2, OnFileMruFile2)
   ON_COMMAND(ID_FILE_MRU_FILE3, OnFileMruFile3)
   ON_COMMAND(ID_FILE_MRU_FILE4, OnFileMruFile4)
   ON_WM_PAINT()  
   ON_COMMAND(ID_SAUVER, OnSauve) 
   ON_COMMAND(ID_ZONE,   OnSelectionneZone)
   ON_COMMAND(ID_ZONE_MANUELLE, OnSelectionneZoneClavier)
   ON_COMMAND(ID_RESTAURE_IMAGE,OnRestaureImage)
   ON_COMMAND(ID_MIROIR_VERT, OnMiroirVert)
   ON_COMMAND(ID_MIROIR_HORI, OnMiroirH)
   ON_COMMAND(ID_ROBERTS,OnRoberts)  
   ON_COMMAND(ID_PREWITT,OnPrewitt)   
   ON_COMMAND(ID_SOBEL,OnSobel) 
   ON_COMMAND(ID_MULTI_ECHELLE,OnMultiEchelle)
   ON_COMMAND(ID_FFT,OnFFT)
   ON_COMMAND(ID_SEUILLAGE1, OnSeuillerImage1) 	
   ON_COMMAND(ID_SEUILLAGE2, OnSeuillerImage2) 	 
   ON_COMMAND(ID_INVERSION_VIDEO, OnInversionVideo) 
   ON_COMMAND(ID_COULEUR,OnCouleur)     
   ON_COMMAND(ID_AFFICHE_PALETTE,OnPaletteAffiche) 
   ON_COMMAND(ID_NB,OnNoirEtBlanc)
   ON_COMMAND(ID_INFO_GENERALE, OnInformationGenerale)   
   ON_COMMAND(ID_VERTICALE,OnSelectionneCoupeV)
   ON_COMMAND(ID_OBLIQUE,OnSelectionneCoupeOblique)        
   ON_COMMAND(ID_HORIZONTALE,OnSelectionneCoupeH) 
   ON_COMMAND(ID_HISTOGRAMME,OnHistogramme)   
   ON_COMMAND(ID_FONCTION_REPARTITION,OnFonctionRepartition)
   ON_COMMAND(ID_LIRE_SEQUENCE, OnLireSequence )    
   ON_COMMAND(ID_SAUVE_SEQUENCE, OnSauveSequence ) 
   ON_COMMAND(ID_DETRUIRE_SEQUENCE, OnDetruireSequence )  
   ON_WM_LBUTTONDOWN()               
   ON_WM_LBUTTONUP()
   ON_WM_MOUSEMOVE()
   ON_WM_KEYDOWN()
   ON_WM_SIZE() 
   ON_WM_VSCROLL()   
   ON_COMMAND(ID_ETIQUETTAGE, OnEtiquettage)
   ON_COMMAND(ID_ROTATIONS_DROITE, OnRotationsDroite)
   ON_COMMAND(ID_ROTATIONS_GAUCHE, OnRotationsGauche)
   ON_COMMAND(ID_COLOR, OnSeparationRGB)
   ON_COMMAND(ID_NORMALISATION, OnNormalisation)
   ON_COMMAND(ID_HAAR_FLOAT, OnHaarFloat)
   ON_COMMAND(ID_SPLINE_FLOAT, OnSplineFloat)
   ON_COMMAND(ID_GRAD_ROBERTS, OnGradRoberts)
   ON_COMMAND(ID_OP_PREWITT, OnOpPrewitt)
   ON_COMMAND(ID_LAPLACIEN1, OnLaplacien1)
   ON_COMMAND(ID_LAPLACIEN_2, OnLaplacien2)
   ON_COMMAND(ID_ECHANT_LIGNE_COL, OnEchantLigneCol)
   ON_COMMAND(ID_BASSE, OnBasse)
   ON_COMMAND(ID_COPIE_EN_FLOAT, OnCopieEnFloat)
   ON_COMMAND(ID_COPIE_EN_BYTE, OnCopieEnByte)
   ON_COMMAND(ID_MISE_A_ZERO, OnMiseAZero)    
   ON_COMMAND(ID_CHROMINANCE, OnLumchro)
   ON_COMMAND(ID_ZOOM_100, OnZoom100)
   ON_COMMAND(ID_ZOOM_25, OnZoom25)
   ON_COMMAND(ID_ZOOM_50, OnZoom50)
   ON_COMMAND(ID_HAAR, OnHaar)
   ON_COMMAND(ID_BSPLINE, OnBSpline)
   ON_COMMAND(ID_RECONSTRUCT, OnReconstruct)
   ON_COMMAND(ID_CONT, OnCont)
   ON_COMMAND(ID_SOUS_ECHAN_COL, OnSousEchanCol)
   ON_COMMAND(ID_SOUS_ECHAN_LIGNE, OnSousEchanLigne)
   ON_COMMAND(ID_AFFECTER_ECHELLE, OnAffecterEchelle)
   ON_WM_CREATE()
   ON_UPDATE_COMMAND_UI(ID_AFFICHE_PALETTE, OnUpdateAffichePalette)
   ON_WM_MDIACTIVATE()
   ON_UPDATE_COMMAND_UI(ID_NB, OnUpdateNb)
   ON_UPDATE_COMMAND_UI(ID_COULEUR, OnUpdateCouleur)
   ON_WM_DESTROY()
   ON_COMMAND(ID_INFO_IMAGE, OnInfoImage)
   ON_UPDATE_COMMAND_UI(ID_INFO_IMAGE, OnUpdateInfoImage)
   ON_COMMAND(ID_SAUVE_EXCEL, OnSauveExcel)
   ON_UPDATE_COMMAND_UI(ID_ZONE, OnUpdateZone)
   ON_WM_HSCROLL()
   ON_COMMAND(ID_EDIT_COPY, OnEditCopy)
   ON_COMMAND(ID_FILE_PRINT, OnFilePrint)
   ON_COMMAND(ID_EDIT_CUT, OnEditCut)
   ON_COMMAND(ID_EDIT_PASTE, OnEditPaste)
   ON_COMMAND(ID_MOYENNEURXY121, OnMoyenneurxy121)
   ON_COMMAND(ID__MOYENNEURXY111, OnMoyenneurxy111)
   ON_COMMAND(ID_MOYENNEURX111, OnMoyenneurx111)
   ON_COMMAND(ID_MOYENNEURY111, OnMoyenneury111)
   ON_COMMAND(ID_ZOOM_X, OnZoomX)
   ON_COMMAND(ID_EXTRA_PAUL, OnExtraPaul)
   ON_UPDATE_COMMAND_UI(ID_HORIZONTALE, OnUpdateHorizontale)
   ON_UPDATE_COMMAND_UI(ID_OBLIQUE, OnUpdateOblique)
   ON_UPDATE_COMMAND_UI(ID_VERTICALE, OnUpdateVerticale)
   ON_COMMAND(ID_SAUVER_SOUS, OnSauverSous)
   ON_WM_CLOSE()
   ON_UPDATE_COMMAND_UI(ID_LIRE_SEQUENCE, OnUpdateLireSequence)
   ON_UPDATE_COMMAND_UI(ID_DETRUIRE_SEQUENCE, OnUpdateDetruireSequence)
   ON_UPDATE_COMMAND_UI(ID_SAUVE_SEQUENCE, OnUpdateSauveSequence)
   ON_COMMAND(ID_GRAB, OnGrab)
   ON_UPDATE_COMMAND_UI(ID_GRAB, OnUpdateGrab)
   ON_COMMAND(ID_SNAP, OnSnap)
   ON_UPDATE_COMMAND_UI(ID_SNAP, OnUpdateSnap)
   ON_WM_SHOWWINDOW()
   ON_COMMAND(ID_TRAITEMENT_MORPHOLOGIE_EROSION, OnTraitementMorphologieErosion)
   ON_COMMAND(ID_TRAITEMENT_MORPHOLOGIE_DILATATION, OnTraitementMorphologieDilatation)
   ON_COMMAND(ID_TRAITEMENT_MORPHOLOGIE_OUVERTURE, OnTraitementMorphologieOuverture)
   ON_COMMAND(ID_TRAITEMENT_MORPHOLOGIE_FERMETURE, OnTraitementMorphologieFermeture)
   ON_COMMAND(ID_OUTILS_PSEUDORELIEF, OnOutilsPseudorelief)
   ON_COMMAND(ID_TRANSFORMATION_CONVERSIONENBSPLINE, OnTransformationConversionenbspline)
   ON_COMMAND(ID_TRAITEMENT_NORMALISATIONDESCOMPOSANTES, OnTraitementNormalisationdescomposantes)
   ON_COMMAND(ID_TRAITEMENT_PETROU, OnTraitementPetrou)
   ON_COMMAND(ID_TRAITEMENT_CONTOURPETROU, OnTraitementContourpetrou)
   ON_COMMAND(ID_TRAITEMENT_CONVOLUTIONGAUSSIEN, OnTraitementConvolutiongaussien)
   ON_COMMAND(ID_TRAITEMENT_CONVOLUTIONGAUSSIEN_GAUSSIEN1DX, OnTraitementConvolutiongaussienGaussien1dx)
   ON_COMMAND(ID_TRAITEMENT_CONVOLUTIONGAUSSIEN_GAUSSIEN1DY, OnTraitementConvolutiongaussienGaussien1dy)
   ON_COMMAND(ID_TRAITEMENT_CONVOLUTIONGAUSSIEN_CURVILINEARDETECTION, OnTraitementConvolutiongaussienCurvilineardetection)
   ON_COMMAND(ID_TRAITEMENT_CURVILINEAR, OnTraitementCurvilinear)
   ON_COMMAND(ID_COLOR_CURVI, OnColorCurvi)
   ON_WM_RBUTTONDBLCLK()
   ON_COMMAND(ID_GRAD_Y, OnGradY)
   ON_COMMAND(ID_TRAITEMENT_ESTIMATIONTEMPSFFT, OnTraitementEstimationtempsfft)
   ON_COMMAND(ID_TRAITEMENT_HYPERGRADIENT_PARBLOC, OnTraitementHypergradientParbloc)
   ON_COMMAND(ID_TRAITEMENT_HYPERGRADIENT_FENETREGLISSANTE, OnTraitementHypergradientFenetreglissante)
   ON_COMMAND(ID_TRAITEMENT_HYPERGRADIENT_FENETREGLISSANTEZEROPADDING, OnTraitementHypergradientFenetreglissantezeropadding)
   ON_WM_SYSKEYDOWN()
   ON_MESSAGE(WM_COUPE, OnSetCoupeVariable)  
   ON_COMMAND(ID_TRAITEMENT_OUTILS_ETUDEDISTANCEMAXI, OnTraitementOutilsEtudedistancemaxi)
   //}}AFX_MSG_MAP
   ON_COMMAND(ID_OUTILS_SUIVIDECONTOURS, &CImageWnd::OnOutilsSuividecontours)
   ON_COMMAND(ID_OUTILS_PARAMCURVIDETECT, &CImageWnd::OnOutilsParamcurvidetect)
   ON_COMMAND(ID_MATCHING_SELECTION3POINTS, &CImageWnd::OnMatchingSelection3points)
   ON_COMMAND(ID_TRANSFORMATIONIMAGE_PARAMTRANSFORMATION, &CImageWnd::OnTransformationimageParamtransformation)
   ON_COMMAND(ID_TRANSFORMATIONIMAGE_RUN, &CImageWnd::OnTransformationimageRun)
   ON_COMMAND(ID_MATCHING_CALCULMATRICEPASSAGE, &CImageWnd::OnMatchingCalculmatricepassage)
   ON_COMMAND(ID_MATCHING_CALCULMATRICEPASSAGE504, &CImageWnd::OnMatchingCalculmatricepassage504)
   ON_COMMAND(ID_MATCHING_SELECTION4POINTS, &CImageWnd::OnMatchingSelection4points)
   ON_COMMAND(ID_TRANSFORMATIONIMAGE_RUN506, &CImageWnd::OnTransformationimageRun506)
   ON_COMMAND(ID_MATCHING_REPETABILIT507, &CImageWnd::OnMatchingRepetabilit507)
   ON_COMMAND(ID_TRAITEMENT_PETROUCOULEUR, &CImageWnd::OnTraitementPetroucouleur)
   ON_COMMAND(ID_TRAITEMENT_CURVICOULEUR, &CImageWnd::OnTraitementCurvicouleur)
   ON_COMMAND(ID_TRAITEMENT_GAUSSIENCOULEUR, &CImageWnd::OnTraitementGaussiencouleur)
END_MESSAGE_MAP()                            



/////////////////////////////////////////////////////////////////////////////
// CImageWnd creation
BOOL CImageWnd::Create(LPCSTR szTitle, LONG style ,CMDIFrameWnd* parent)
{ 
   BOOL LeBool = FALSE;
   if (PImage!=NULL)
   {
      CRect rectMain;
      (AfxGetApp()->m_pMainWnd)->GetClientRect(rectMain);
      // 27 hauteur du menu, 22 Hauteur de la status bar
      LFene = (int) (Larg+12); //8
      HFene = (int) (Haut+31+HautWndInfo);

      if (HFene>rectMain.Height()-100)
      {
         style = style | WS_VSCROLL;	  
         HFene = rectMain.Height()-100;
      }
      if (LFene>rectMain.Width()-50)
      {
         style = style | WS_HSCROLL;	  
         LFene = rectMain.Width()-50;
      }

      // dimension de la fenetre
      CRect rect(0,0,LFene,HFene); 
      ChargeMenu();

      if (SEQ_EXISTE) menu.InsertMenu(0,MF_BYPOSITION,ID_LIRE_SEQUENCE,"Lire la séquence");    

      // creation de la fenetre                                           
      if (TypeImage == BMP)
      {
         LeBool = CMDIChildWnd::Create	( 
            AfxRegisterWndClass( CS_DBLCLKS |CS_HREDRAW | CS_VREDRAW | CS_BYTEALIGNWINDOW | CS_SAVEBITS,/*LoadCursor(NULL, IDC_ARROW )*/NULL,/*(HBRUSH)(COLOR_WINDOW+1)*/(HBRUSH)(COLOR_BACKGROUND),LoadIcon(AfxGetInstanceHandle(),MAKEINTRESOURCE(IDI_BMP)))
            , szTitle, style ,rect, parent);
      }
      else if (TypeImage == MEM)
      {
         LeBool = CMDIChildWnd::Create	( 
            AfxRegisterWndClass( CS_DBLCLKS |CS_HREDRAW | CS_VREDRAW | CS_BYTEALIGNWINDOW | CS_SAVEBITS,/*LoadCursor(NULL, IDC_ARROW )*/NULL,/*(HBRUSH)(COLOR_WINDOW+1)*/(HBRUSH)(COLOR_BACKGROUND),LoadIcon(AfxGetInstanceHandle(),MAKEINTRESOURCE(IDI_MEM)))
            , szTitle, style ,rect, parent);
      }
      else if (TypeImage == FLOAT32)
      {
         FloatToImage();
         LeBool = CMDIChildWnd::Create	( 
            AfxRegisterWndClass( CS_DBLCLKS |CS_HREDRAW | CS_VREDRAW | CS_BYTEALIGNWINDOW | CS_SAVEBITS,/*LoadCursor(NULL, IDC_ARROW )*/NULL,/*(HBRUSH)(COLOR_WINDOW+1)*/(HBRUSH)(COLOR_BACKGROUND),LoadIcon(AfxGetInstanceHandle(),MAKEINTRESOURCE(IDI_F32)))
            , szTitle, style ,rect, parent);
      }

      // fenetre d'information	
      GetClientRect(TheRectClient);  
      AfficheDimensions(); // fenetre info n°1
      GetWindowText(Titre);   
      Memorise();
      IMAGE_MODIFIEE = FALSE;
      menu.CheckMenuItem(ID_ZOOM_100,MF_CHECKED);  

   }
   return LeBool;
}

void	CImageWnd::ChargeMenu()
{   
   if (menu.m_hMenu==NULL) 
   {
      if (TypeImage == MEM || TypeImage == BMP)
      {
         if (TypeColor == COLOR)  menu.LoadMenu(IDR_COLOR);
         else menu.LoadMenu(IDR_IMAGE);
      }
      if (TypeImage == FLOAT32)
      {
         menu.LoadMenu(IDR_FLOAT);
      }
      m_hMenuShared = menu.m_hMenu;
   }
}

// *********************************************************
//  Constructeur de CImage
// **********************************************************
CImageWnd::CImageWnd()
: m_LargeurEtude(0)
, m_LargeurMini(0)
, VLocaleLargeur(0)
, MATCHINGCURVI(FALSE)
, CptManipMatching(0)
{
}

CImageWnd::CImageWnd(LONG InitLarg,LONG InitHaut,NBCOLOR InitColor, BOOL SEQ)
{   
   TypeColor = InitColor;  
   TypeImage = BMP; 
   InitGlobal(InitLarg,InitHaut);  
   SEQ_EXISTE =  SEQ;
   PImage = AlloueMemoire(Larg,Haut);
   PImageTemp = AlloueMemoire(Larg,Haut);      
   ImageRGB = (rgb  *) PImage;  
   Image = (unsigned char  *) PImage;

   TheTete = (BITMAPFILEHEADER*) new unsigned char[sizeof(BITMAPFILEHEADER)];    
   SetTeteBitmap(TheTete,Larg,Haut,TypeColor);    
   if (TypeColor == NB) MyInitInfo = (BITMAPINFO*) new unsigned char[ sizeof(BITMAPINFOHEADER) + sizeof(RGBQUAD)*256];
   else MyInitInfo = (BITMAPINFO*) new unsigned char[ sizeof(BITMAPINFOHEADER)];
   SetInfoBitmap(MyInitInfo,Larg,Haut,TypeColor);   
   if  (TypeColor == NB)
   {
      MyPalette->CreatePalette(MyLogPaletteNB);
      Mapping_Palette(MyLogPaletteNB);
      TypePalette = PALNB;
   }	
   // structure de l'image 
   SetSImage();
}  

CImageWnd::CImageWnd(long InitLarg, long InitHaut, TYPEIMAGE InitLeType)
{
   TypeColor = NB;  
   TypeImage = InitLeType;
   InitGlobal(InitLarg,InitHaut);  
   PImage = AlloueMemoire(Larg,Haut);
   PImageTemp = AlloueMemoire(Larg,Haut);      
   ImageRGB = (rgb  *) PImage;  
   Image = (unsigned char  *) PImage;

   TheTete = (BITMAPFILEHEADER*) new unsigned char[sizeof(BITMAPFILEHEADER)];    
   SetTeteBitmap(TheTete,Larg,Haut,TypeColor);    
   if (TypeColor == NB) MyInitInfo = (BITMAPINFO*) new unsigned char[ sizeof(BITMAPINFOHEADER) + sizeof(RGBQUAD)*256];
   else MyInitInfo = (BITMAPINFO*) new unsigned char[ sizeof(BITMAPINFOHEADER)];
   SetInfoBitmap(MyInitInfo,Larg,Haut,TypeColor);   

   if (TypeImage == FLOAT32)
   {
      ImageFloat = AlloueMemoireFloat(Larg,Haut);
      SetSImageFloat();
   }
   // structure de l'image 
   SetSImage();
   MyPalette->CreatePalette(MyLogPaletteNB);
   Mapping_Palette(MyLogPaletteNB);
}

CImageWnd::CImageWnd(CString InitFileName, long InitLarg, long InitHaut, TYPEIMAGE InitTypeImage)
{  
   TypeImage = InitTypeImage;
   FileName = InitFileName;
   InitGlobal(InitLarg,InitHaut);
   if (TypeImage == BMP)
      InitialisationPourBMP();
   else
   {
      TypeColor = NB;
      PImage = AlloueMemoire(Larg,Haut);
      PImageTemp = AlloueMemoire(Larg,Haut);      
      ImageRGB = (rgb  *) PImage;  
      Image = (unsigned char  *) PImage;
      if (TypeImage == FLOAT32)
      {
         ImageFloat = AlloueMemoireFloat(Larg,Haut);
         InitialisationPourFLOAT();
         SetSImageFloat();
      }
      else if (TypeImage == MEM)
         InitialisationPourMEM();

      TheTete = (BITMAPFILEHEADER*) new unsigned char[sizeof(BITMAPFILEHEADER)];    
      SetTeteBitmap(TheTete,Larg,Haut,TypeColor);    
      if (TypeColor == NB) MyInitInfo = (BITMAPINFO*) new unsigned char[ sizeof(BITMAPINFOHEADER) + sizeof(RGBQUAD)*256];
      else MyInitInfo = (BITMAPINFO*) new unsigned char[ sizeof(BITMAPINFOHEADER)];
      SetInfoBitmap(MyInitInfo,Larg,Haut,TypeColor);   
   }
   // structure de l'image 
   SetSImage();

   MyPalette->CreatePalette(MyLogPaletteInit);
   Mapping_Palette(MyLogPaletteInit);
}  

void CImageWnd::InitGlobal(LONG InitLarg,LONG InitHaut)
{
   LaZone = Appuie = LaCoupe = zoneSelectionnee = CoupeSelectionnee = PALETTE_AFFICHEE = FALSE;
   INFO_IMAGE = FALSE;
   TRACER_COUPE = FALSE;
   IMAGE_MODIFIEE = FALSE;

   NbProfil=0;
   pCoupeVariable=NULL;
   PImage = NULL;
   ImageFloat = NULL;
   mini = 0.F; maxi = 255.F;

   dialInfo = new CDialInfoTrait(this);
   NomTraitement = " ";
   CommentTraitement = " ";

   SEQUENCE_CREE = FALSE;  
   NbMaxImageSeq = 100; 
   SEQ_EXISTE =  FALSE;  
   DialLecSeq=new CDialLecSeq((CImageWnd*)this);
   // initialisation des seuils
   Seuil = 128;
   Seuil2 = 250;               
   // Initialisation de l'échelle par défaut
   Echelle =1;

   // pointeur sur la fenetre palette
   pPaletteWnd	=	NULL;

   // position de la souris dans l'image
   I=J=0;
   // position verticale du scrollbar
   PosVScroll = 0; PosHScroll = 0;
   PosHFen = 0; PosVFen = 0;	
   // Hauteur des fenêtres informations
   HautWndInfo =20;
   // visualisation de l'image 
   ZOOM = 1;

   // initialisation des dimensions de l'image et réservation mémoire
   Larg = InitLarg;    
   Haut = InitHaut; 
   if (div((int)Larg,4).rem!=0) Larg = Larg-div((int)Larg,4).rem;
   // zone rectangle
   debut_i = 0;fin_j = Haut-1;   
   debut_j = 0;fin_i = Larg-1;   
   // zone pour la sequence
   debut_i_seq = 0;fin_j_seq = Haut-1;   
   debut_j_seq = 0;fin_i_seq = Larg-1;	
   // Création de la palette 
   MyLogPaletteNB = (LOGPALETTE*) new BYTE [sizeof(LOGPALETTE)+sizeof(PALETTEENTRY)*256];
   MyLogPaletteNB->palVersion = 0x300;
   MyLogPaletteNB->palNumEntries = 256;
   Palette_NB();

   MyLogPaletteCOUL = (LOGPALETTE*) new BYTE [sizeof(LOGPALETTE)+sizeof(PALETTEENTRY)*256];
   MyLogPaletteCOUL->palVersion = 0x300;
   MyLogPaletteCOUL->palNumEntries = 256;
   Palette_Couleurs();

   MyLogPaletteInit =(LOGPALETTE*) new unsigned char[sizeof(LOGPALETTE) + sizeof(PALETTEENTRY)*256];
   MyLogPaletteInit->palVersion = 0x300;
   MyLogPaletteInit->palNumEntries = 256;
   for(int i=0;i<255;i++)
   {
      MyLogPaletteInit->palPalEntry[i].peBlue  = i;
      MyLogPaletteInit->palPalEntry[i].peGreen = i;
      MyLogPaletteInit->palPalEntry[i].peRed   = i;
      MyLogPaletteInit->palPalEntry[i].peFlags = PC_NOCOLLAPSE;
   }
   MyLogPaletteInit->palPalEntry[255].peBlue  = 255;
   MyLogPaletteInit->palPalEntry[255].peGreen = 255;
   MyLogPaletteInit->palPalEntry[255].peRed   = 255;
   MyLogPaletteInit->palPalEntry[255].peFlags = PC_NOCOLLAPSE;

   MyPalette = new CPalette();
   TypePalette=PALINIT;
   MATCHINGCURVI = FALSE;
   MATCHINGCURVI4 = FALSE;
}

void CImageWnd::InitialisationPourBMP()
{
   if (LeFichier.Open(FileName, CFile::modeRead))
   {   
      TheTete = (BITMAPFILEHEADER*) new BYTE[sizeof(BITMAPFILEHEADER)];    
      MyInitInfo = (BITMAPINFO*) new BYTE[ sizeof(BITMAPINFOHEADER) + sizeof(RGBQUAD)*256];

      if (LeFichier.Read(TheTete,sizeof(BITMAPFILEHEADER))) 
         if (TheTete->bfType==19778)
         {
            LeFichier.Read(MyInitInfo,sizeof(BITMAPINFOHEADER) + sizeof(RGBQUAD)*256);
            switch (MyInitInfo->bmiHeader.biBitCount)
            {
            case 8 		: TypeColor = NB;		break;    // BMP N&B
            case 24		: TypeColor = COLOR;	break;    // BMP COLOR
            }	
            if((TypeColor!=NB)&&(TypeColor!=COLOR))
            {
               MessageBox("Impossible de lire ce format de BMP","winseg 8.2",MB_OK|MB_ICONSTOP);
               return;
            }

            //allocation memoire du plan image
            PImage = AlloueMemoire(Larg,Haut);
            PImageTemp = AlloueMemoire(Larg,Haut);      
            ImageRGB = (rgb  *) PImage;  
            Image = (unsigned char  *) PImage;

            for(int i=0;i<256;i++)
            {
               MyLogPaletteInit->palPalEntry[i].peBlue  = MyInitInfo->bmiColors[i].rgbBlue;
               MyLogPaletteInit->palPalEntry[i].peGreen = MyInitInfo->bmiColors[i].rgbGreen;
               MyLogPaletteInit->palPalEntry[i].peRed   = MyInitInfo->bmiColors[i].rgbRed;
               MyLogPaletteInit->palPalEntry[i].peFlags = PC_NOCOLLAPSE;
            }
            switch (TypeColor)
            {
            case NB : // lecture de la palette du fichier
               LeFichier.Seek(TheTete->bfOffBits,CFile::begin);
               SetImage(Image);
               break;
            case COLOR : 
               LeFichier.Seek(TheTete->bfOffBits,CFile::begin);
               SetImage(ImageRGB); 
               break;
            }
         }	 
         LeFichier.Close();
   }	
}
void CImageWnd::InitialisationPourFLOAT()
{                                

   LeFichier.Open(FileName, CFile::modeRead);     
   LeFichier.Read(ImageFloat,Larg*Haut*sizeof(float));		    				
   LeFichier.Close();
   FloatToImage();
}

void CImageWnd::InitialisationPourMEM()
{
   LeFichier.Open(FileName, CFile::modeRead);     
   LeFichier.Read(Image,Larg*Haut);		    				
   LeFichier.Close();
}

// ******************** fin constructeur **********************************

void * CImageWnd::AlloueMemoire(long L, long H)
{ 
   void * PBuf;
   int t=L*H;;
   switch (TypeColor)
   {   
   case NB 		: 	{ 
      Taille = t;
      PBuf 	=  (unsigned char   *) calloc(Taille,1);
      break;
                  }	
   case COLOR 		: 	{ 

      Taille = sizeof(rgb)*t;
      PBuf =  (rgb   *) calloc(t,sizeof(rgb));
      break;
                     }   
   } 
   return PBuf;
}

float  * CImageWnd::AlloueMemoireFloat(long L,long H)
{                                         
   float  * p;
   p = (float  *) malloc((DWORD)L*(DWORD)H*(DWORD)sizeof(float));
   return p;
}


void CImageWnd::SetImage(unsigned char  *  PBuf)
{
   for (long j=0;j<Haut;j++)
   {   
      LeFichier.Read(&PBuf[j*Larg],Larg);
      if (res.rem != 0) LeFichier.Seek(res.rem,CFile::current); 
   }
}

void CImageWnd::SetImage(rgb  *  PBuf)
{
   for (long j=0;j<Haut;j++)
   {   
      LeFichier.Read(&PBuf[j*Larg],Larg*sizeof(rgb));
      if (res.rem != 0) LeFichier.Seek((res.rem)*sizeof(rgb), CFile::current);
   }
}		      	

void CImageWnd::ImageToFloat()
{
   if (ImageFloat == NULL)	ImageFloat = AlloueMemoireFloat(Larg,Haut);
   for (long i=0; i<Larg*Haut; i++)	ImageFloat[i] = (float) Image[i];
}

void  CImageWnd::FloatToImage()
{
   if (Echelle==1) Normalisation(0,0,Larg-1,Haut-1);
   else 
   {  
      Normalisation(0,0,Larg/Echelle-1,Haut/Echelle-1); //1         //34
      Normalisation(Larg/Echelle,0,Larg-1,Haut/Echelle-1);//2       //12
      Normalisation(0,Haut/Echelle,Larg/Echelle,Haut-1);//3
      Normalisation(Larg/Echelle,Haut/Echelle,Larg-1,Haut-1);//4
   } 
   /*
   // ------------ Variables locales (Dimensions de  l'image float ... ----
   float temp;	 	
   float norm;

   //-------------------- Création de l'image destination ------------
   // -------------------- Module de conversion float -> 256 niveaux de gris ------
   // ----- et affectation du buffer ImageFloat ----------------------

   //-- recherche du maximun et du minimun dans l'image float ---

   // NORMALISATION ENTRE MIN ET MAX
   mini=maxi=ImageFloat[0];		
   for (long l = 0 ; l<Larg*Haut;l++) 
   {
   temp=floor(ImageFloat[l]+0.5F);
   if (temp>maxi) maxi=temp;
   if (temp<mini) mini=temp;
   }

   // -- Création de la nouvelle image	normalisée                 
   norm = (maxi-mini)/255;	
   for ( l = 0 ; l<Larg*Haut;l++) 
   {
   Image[l]= (unsigned char) ( (floor(ImageFloat[l]+0.5F)-mini)/norm );
   }			

   /*
   // NORMALISATION ENTRE 0 ET 4000
   mini = 0;
   maxi = 4000;
   norm = (maxi-mini)/255;	
   for (long l = 0 ; l<Larg*Haut;l++) 
   {
   if (ImageFloat[l]<=0)	Image[l] = 0;
   else if (ImageFloat[l]>4000)	Image[l] = 255;
   else Image[l]= (unsigned char) ( (ImageFloat[l]-mini)/norm );
   }
   */

   /*
   // NORMALISATION EXPONENTIELLE
   mini=maxi=ImageFloat[0];		
   for (long l = 0 ; l<Larg*Haut;l++) 
   {
   temp=ImageFloat[l];
   if (temp>maxi) maxi=temp;
   if (temp<mini) mini=temp;
   }

   // -- Création de la nouvelle image	normalisée                 
   norm = (maxi-mini);	
   for ( l = 0 ; l<Larg*Haut;l++) 
   {
   Image[l]= (unsigned char) ( 256*(1-exp((-5.5*(ImageFloat[l]-mini))/(norm))) );
   }

   // NORMALISATION ARC TANGENTE (Plage et Centre)
   float Plage = 1500.0F, Centre = 2000.0F;
   mini=maxi=ImageFloat[0];		
   for (long l = 0 ; l<Larg*Haut;l++) 
   {
   temp=ImageFloat[l];
   if (temp>maxi) maxi=temp;
   if (temp<mini) mini=temp;
   }

   // -- Création de la nouvelle image	normalisée                 
   norm = (maxi-mini);	
   for ( l = 0 ; l<Larg*Haut;l++) 
   {
   Image[l]= (BYTE)(255*(atan((2*(ImageFloat[l]-Centre))/Plage)-atan((2*(mini-Centre))/Plage))/(atan((2*(maxi-Centre))/Plage)-atan((2*(mini-Centre))/Plage)));
   }*/
}

float CImageWnd::GetPixelNB(long index)
{   
   if (TypeColor == NB) return Image[index];   
   else return 0.0F;
}


void CImageWnd::SetSImage()
{
   // structure de l'image 
   Im.Larg 	= Larg;
   Im.Haut 	= Haut;
   Im.Image 	= Image;
   Im.Lmin		= 0;
   Im.Lmax		= Larg-1;
   Im.Hmin		= 0;
   Im.Hmax		= Haut-1;
}

void CImageWnd::SetSImageFloat()
{
   // structure de l'image float
   ImFloat.Larg 	= Larg;
   ImFloat.Haut 	= Haut;
   ImFloat.Image 	= ImageFloat;
   ImFloat.Lmin	= 0;
   ImFloat.Lmax	= Larg-1;
   ImFloat.Hmin	= 0;
   ImFloat.Hmax	= Haut-1;
}

void CImageWnd::SetSImageRGB()
{
   // structure de l'image 
   ImRGB.Larg 	= Larg;
   ImRGB.Haut 	= Haut;
   ImRGB.Image = ImageRGB;
   ImRGB.Lmin	= 0;
   ImRGB.Lmax	= Larg-1;
   ImRGB.Hmin	= 0;
   ImRGB.Hmax	= Haut-1;
}

SImage  CImageWnd::GetSImage()
{
   return Im;
}

SImageRGB  CImageWnd::GetSImageRGB()
{
   return ImRGB;
}


SImageFloat  CImageWnd::GetSImageFloat()
{
   return ImFloat;
}

// *************************************************                                                     
// routine appelée lorsque la fenêtre est retaillée 
// *************************************************
afx_msg void CImageWnd::OnSize( UINT nType, int cx, int cy )
{   
   if (cx<Larg) //Gestion Scroll Horizontal
      ModifyStyle(NULL, WS_HSCROLL, SWP_SHOWWINDOW);
   else
      ModifyStyle(WS_HSCROLL, NULL, SWP_SHOWWINDOW);
   if (cy<Haut+HautWndInfo)	//Gestion Scroll Vertical
      ModifyStyle(NULL, WS_VSCROLL, SWP_SHOWWINDOW);
   else
      ModifyStyle(WS_VSCROLL, NULL, SWP_SHOWWINDOW);

   CMDIChildWnd::OnSize(  nType,  LFene,  HFene ); 
   GetClientRect(TheRectClient); 

   UINT nID, nStyle;
   int nWidth;
   m_wndStatusBarImage.GetPaneInfo(0, nID, nStyle, nWidth);
   m_wndStatusBarImage.SetPaneInfo(0, nID, nStyle, LFene/3);
}



afx_msg void CImageWnd::OnVScroll(UINT nSBCode, UINT nPos, CScrollBar * pScrollBar)
{ 

   if (ZOOM==1)
   {				                        
      GetClientRect (TheRectClient);     
      if ((TheRectClient.Height()-HautWndInfo<Haut))    
      { 
         int  tempPos = PosVScroll;
         switch( nSBCode )
         {
         case SB_PAGEDOWN 		: PosVScroll = PosVScroll +20; 	break;    
         case SB_LINEDOWN 		: PosVScroll ++;  				break;
         case SB_PAGEUP 			: PosVScroll = PosVScroll-20;  	break;   
         case SB_LINEUP 			: PosVScroll --;  				break;   
         case SB_THUMBPOSITION 	: PosVScroll = nPos;  			break;
         }   
         if (PosVScroll>100) PosVScroll = 100;    
         if (PosVScroll<0) PosVScroll = 0;   
         PosHFen = PosVScroll*((int)Haut-TheRectClient.Height()+(int)HautWndInfo)/100;                                         
         if (tempPos != PosVScroll) Invalidate(FALSE);
      }
      else 
      {  
         PosVScroll = 0;
      }
      SetScrollPos(SB_VERT,PosVScroll,TRUE); 	
   }
}

void CImageWnd::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar) 
{
   if (ZOOM==1)
   {				                        
      GetClientRect (TheRectClient);     
      if ((TheRectClient.Width()<Larg))    
      { 
         int  tempPos = PosHScroll;
         switch( nSBCode )
         {
         case SB_PAGERIGHT 		: PosHScroll = PosHScroll +20; 	break;    
         case SB_LINERIGHT 		: PosHScroll ++;  				break;
         case SB_PAGELEFT		: PosHScroll = PosHScroll-20;  	break;   
         case SB_LINELEFT 		: PosHScroll --;  				break;   
         case SB_THUMBPOSITION 	: PosHScroll = nPos;  			break;
         }   
         if (PosHScroll>100) PosHScroll = 100;    
         if (PosHScroll<0) PosHScroll = 0;   
         PosVFen = PosHScroll*((int)Larg-TheRectClient.Width())/100;                                         
         if (tempPos != PosHScroll) Invalidate(FALSE);
      }
      else 
      {  
         PosHScroll = 0;
      }
      SetScrollPos(SB_HORZ,PosHScroll,TRUE); 	
   }
}

// *********************************************
//   DESTRUCTEUR DE CImage
// *********************************************
CImageWnd::~CImageWnd()
{   
   int i=0;
   MyPalette->DeleteObject();                             
   if (PALETTE_AFFICHEE) pPaletteWnd->DestroyWindow();
   delete TheTete;
   delete MyInitInfo;                                  
   delete MyLogPaletteInit;    
   delete MyLogPaletteCOUL;    
   delete MyLogPaletteNB;
   if (INFO_IMAGE) delete dialInfo;

   free(PImage);
   free(PImageTemp);      
   if (TypeImage == FLOAT32) free(ImageFloat);      
   OnDetruireSequence();
   delete DialLecSeq;


   //liberation pour des listes pour la detection curvilinéaire

   //----> A voir ça plante 

   /*if(TableauListesPointsFiltrees.GetCount()!=0)
   {
   for(i=0;i<TableauListesPointsFiltrees.GetCount();i++)
   {
   delete TableauListesPointsFiltrees[i];
   }
   }*/

}
// ****************************************************************

// ****************************************************************************************
// 								ONPAINT réponse à WM_PAINT
// ****************************************************************************************

afx_msg void CImageWnd::OnPaint(void)
{
   CPaintDC dc(this);                                       
   GetClientRect (TheRectClient);       
   if (PALETTE_AFFICHEE) 	Mapping_Palette(pPaletteWnd->MyLogPalette);

   if (ZOOM==1) 
      SetDIBitsToDevice(dc.GetSafeHdc(),0,HautWndInfo, Larg,Haut,PosVFen,-PosHFen,0,Haut,PImage,MyInitInfo,DIB_RGB_COLORS);
   else StretchDIBits(dc.GetSafeHdc(),0,HautWndInfo,(int)(Larg/ZOOM),(int)(Haut/ZOOM),
      PosVFen, -PosHFen,(int) Larg, (int)Haut ,PImage, MyInitInfo,DIB_RGB_COLORS,SRCCOPY);
   ReleaseDC(&dc); // le dc est relaché    

   // dessin du rectangle si zone
   if ((debut_i!=0)||(debut_j!=0)||(fin_i!=Larg-1)||(fin_j!=Haut-1)) DessineRectangle();
   if(SEQUENCE_CREE&&((debut_i_seq!=0)||(debut_j_seq!=0)||(fin_i_seq!=Larg-1)||(fin_j_seq!=Haut-1))) DessineRectangle(debut_i_seq-1,debut_j_seq-2,fin_i_seq+2,fin_j_seq+1,RGB(255,255,0));
   if ((CoupeSelectionnee)&&(pCoupeVariable!=NULL)) DessineLigne(Lg_i1,Lg_j1,Lg_i2,Lg_j2);
   else CoupeSelectionnee = FALSE;

} // fin OnPaint                                  

// *****************************************************************************
// routine d'affichage par exemple pour le continu (for) independamment du timer
// *****************************************************************************
afx_msg void CImageWnd::Affichage()
{
   CClientDC dc(this);                                        
   GetClientRect (TheRectClient);       

   if (ZOOM==1) 
      SetDIBitsToDevice(dc.GetSafeHdc(),0,HautWndInfo, Larg,Haut,PosVFen,-PosHFen,0,Haut,PImage,MyInitInfo,DIB_RGB_COLORS);
   else StretchDIBits(dc.GetSafeHdc(),0,HautWndInfo,(int)(Larg/ZOOM),(int)(Haut/ZOOM),
      PosVFen, -PosHFen,(int) Larg, (int)Haut ,PImage, MyInitInfo,DIB_RGB_COLORS,SRCCOPY);
   ReleaseDC(&dc); // le dc est relaché  
} // fin OnPaint                                  

void CImageWnd::Memorise()
{   
   long l;
   IMAGE_MODIFIEE = TRUE;
   switch (TypeColor)
   {   
   case NB:
      for (l=0;l<Larg*Haut;l++) ((unsigned char  *)PImageTemp)[l] = ((unsigned char  *) PImage)[l]; 
      break;
   case COLOR:
      for (l=0;l<Larg*Haut;l++) ((rgb  * )PImageTemp)[l] = ((rgb  *) PImage)[l];     
      break;
   }
}

// *****************************************************************************
// Routine de sauegarde de l'image en format BITMAP 256 N&B
// *****************************************************************************

void CImageWnd::OnSauverSous() 
{
   if(SEQUENCE_CREE) {OnSauveSequence(); return;}
   ArreteLeContinu();
   CString extension, liste;
   CFile Fichier; // variable fichier BITMAP
   int i1,i2,j1,j2;
   i1=debut_i;
   i2=fin_i;
   j1=debut_j;
   j2=fin_j;

   if ((debut_i!=0) || (fin_i!=Larg-1) || (debut_j!=0) || (fin_j!=Haut-1))
      switch(MessageBox("Voulez-vous sauvegardez que la sélection?","Sauvegarde",MB_YESNOCANCEL))
   {
      case IDNO:	// sauvegarde de la zone sélectionnée de l'image
         debut_i=0;debut_j=0;fin_i=Larg-1;fin_j=Haut-1;
         break;
      case IDCANCEL:
         return;
      default:
         break;
   }

   if (TypeImage == BMP)
   {
      liste = "BITMAP (*.bmp)|*.bmp|BRUTE (*.mem)|*.mem|Float (*.f32)|*.f32||";
      extension = "*.bmp";
   }
   else if(TypeImage == FLOAT32)
   {
      liste = "Float (*.f32)|*.f32|BITMAP (*.bmp)|*.bmp|BRUTE (*.mem)|*.mem||";
      extension = "*.f32";
   }
   else if(TypeImage == MEM)
   {
      liste = "BRUTE (*.mem)|*.mem|Float (*.f32)|*.f32|BITMAP (*.bmp)|*.bmp||";
      extension = "*.mem";
   }

   CFileDialog FileData(FALSE,extension,NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,liste,this);
   if (FileData.DoModal()==IDOK)
      if (Fichier.Open(FileData.GetPathName(),CFile::modeCreate | CFile::modeWrite))
      {   
         if (FileData.GetFileExt() == "bmp" || FileData.GetFileExt() == "f32" || FileData.GetFileExt() == "mem")
         {
            SetCursor(LoadCursor(0,IDC_WAIT));  

            // sauvegarde de l'entête si BMP
            if (FileData.GetFileExt() == "bmp")
            {
               BITMAPFILEHEADER * Tete =  (BITMAPFILEHEADER *) new BITMAPFILEHEADER[sizeof(BITMAPFILEHEADER)] ;
               SetTeteBitmap(Tete,fin_i-debut_i+1,fin_j-debut_j+1,TypeColor);

               BITMAPINFO * Info;
               if (TypeColor == NB) 
               {
                  Info = (BITMAPINFO *) new BITMAPINFO[sizeof(BITMAPINFOHEADER)+256*sizeof(RGBQUAD)];
                  SetInfoBitmap(Info,fin_i-debut_i+1,fin_j-debut_j+1,MyInitInfo); 
               }	
               else
               {
                  Info = (BITMAPINFO *) new BITMAPINFO[sizeof(BITMAPINFOHEADER)]; 
                  SetInfoBitmap(Info,fin_i-debut_i+1,fin_j-debut_j+1,COLOR);         
               }	

               Fichier.Write(Tete,(UINT) sizeof(BITMAPFILEHEADER)); 
               if (TypeColor==NB) Fichier.Write(Info,sizeof(BITMAPINFOHEADER)+256*sizeof(RGBQUAD));      
               else Fichier.Write(Info,sizeof(BITMAPINFOHEADER));   
               delete Tete;
               delete Info;                   

            }	
            if (FileData.GetFileExt() == "f32" && (TypeImage == MEM || TypeImage == BMP))
               ImageToFloat();
            long l;                          
            if ((debut_i==0) && (fin_i==Larg-1) && (debut_j==0) && (fin_j==Haut-1))
            {
               // sauvegarde de l'image
               if (FileData.GetFileExt() == "bmp" || FileData.GetFileExt() == "mem") 
               {
                  if (TypeColor==NB) Fichier.Write((unsigned char  *)PImage,Taille);
                  else Fichier.Write(ImageRGB,Taille);	
               }
               if (FileData.GetFileExt() == "f32") Fichier.Write((float  *)ImageFloat,Larg*Haut*sizeof(float)); 
            } 
            else    
            {
               // sauvegarde de la zone sélectionnée de l'image
               for (long j=debut_j;j<=fin_j;j++)
               { 
                  l = j*Larg+debut_i;
                  if (FileData.GetFileExt() == "bmp" || FileData.GetFileExt() == "mem")
                  {
                     if (TypeColor==NB) Fichier.Write(&Image[l],sizeof(unsigned char)*(fin_i-debut_i+1));	
                     else Fichier.Write(&ImageRGB[l],sizeof(rgb)*(fin_i-debut_i+1));	
                  }
                  if (FileData.GetFileExt() == "f32") Fichier.Write(&ImageFloat[l],sizeof(float)*(fin_i-debut_i+1));	
               } 
            }	
            Fichier.Close();    
            if (FileData.GetFileExt() == "f32" && (TypeImage == MEM || TypeImage == BMP))
            { free(ImageFloat); ImageFloat = NULL;}
            SetCursor(LoadCursor(0,IDC_ARROW));
            IMAGE_MODIFIEE=FALSE;
         }
         else MessageBox("Impossible d'ouvrir \n ce type de fichier image","Fichier autre que bmp, mem ou f32",MB_OK|MB_ICONEXCLAMATION);
      }
      else MessageBox("Sauvegarde impossible","Erreur",MB_OK|MB_ICONEXCLAMATION);
      debut_i=i1;
      fin_i=i2;
      debut_j=j1;
      fin_j=j2;

}
void CImageWnd::OnSauve()
{             
   int i1,i2,j1,j2;
   ArreteLeContinu();
   //if (MessageBox("Enregistrer les modifications"+FileName,"Enregistrement",MB_DEFBUTTON2|MB_YESNO|MB_ICONEXCLAMATION) == IDYES)
   //{
   CFile Fichier; // variable fichier BITMAP
   if (Fichier.Open(FileName,CFile::modeReadWrite))
   {   
      SetCursor(LoadCursor(0,IDC_WAIT));  
      // sauvegarde de l'entête si BMP
      i1=debut_i;
      i2=fin_i;
      j1=debut_j;
      j2=fin_j;

      if ((debut_i!=0) || (fin_i!=Larg-1) || (debut_j!=0) || (fin_j!=Haut-1))
         switch(MessageBox("Voulez-vous sauvegardez que la sélection?","Sauvegarde",MB_YESNOCANCEL))
      {
         case IDNO:	// sauvegarde de la zone sélectionnée de l'image
            debut_i=0;debut_j=0;fin_i=0;fin_j=0;
            break;
         case IDCANCEL:
            Fichier.Close();
            return;
         default:
            break;
      }

      if (TypeImage == BMP)
      {
         BITMAPFILEHEADER * Tete =  (BITMAPFILEHEADER *) new BITMAPFILEHEADER[sizeof(BITMAPFILEHEADER)] ;
         SetTeteBitmap(Tete,fin_i-debut_i+1,fin_j-debut_j+1,TypeColor);

         BITMAPINFO * Info;
         if (TypeColor == NB) 
         {
            Info = (BITMAPINFO *) new BITMAPINFO[sizeof(BITMAPINFOHEADER)+256*sizeof(RGBQUAD)];
            SetInfoBitmap(Info,fin_i-debut_i+1,fin_j-debut_j+1,MyInitInfo); 
         }	
         else
         {
            Info = (BITMAPINFO *) new BITMAPINFO[sizeof(BITMAPINFOHEADER)]; 
            SetInfoBitmap(Info,fin_i-debut_i+1,fin_j-debut_j+1,COLOR);         
         }	

         Fichier.Write(Tete,(UINT) sizeof(BITMAPFILEHEADER)); 
         if (TypeColor==NB) Fichier.Write(Info,sizeof(BITMAPINFOHEADER)+256*sizeof(RGBQUAD));      
         else Fichier.Write(Info,sizeof(BITMAPINFOHEADER));   
         delete Tete;
         delete Info;                   

      }	
      long l;                          
      if ((debut_i==0) && (fin_i==Larg-1) && (debut_j==0) && (fin_j==Haut-1))
      {
         // sauvegarde de l'image
         if (TypeImage == MEM || TypeImage == BMP)
         {
            if (TypeColor==NB) Fichier.Write((unsigned char  *)PImage,Taille);
            else Fichier.Write(ImageRGB,Taille);	
         }
         if (TypeImage == FLOAT32) Fichier.Write((float  *)ImageFloat,Larg*Haut*sizeof(float)); 
      } 
      else    
      {
         for (long j=debut_j;j<=fin_j;j++)
         { 
            l = j*Larg+debut_i;
            if (TypeImage == MEM || TypeImage == BMP)
            {
               if (TypeColor==NB) Fichier.Write(&Image[l],sizeof(unsigned char)*(fin_i-debut_i+1));	
               else Fichier.Write(&ImageRGB[l],sizeof(rgb)*(fin_i-debut_i+1));	
            }
            if (TypeImage == FLOAT32) Fichier.Write(&ImageFloat[l],sizeof(float)*(fin_i-debut_i+1));	
         } 
      }	
      Fichier.Close();    
      debut_i=i1;
      fin_i=i2;
      debut_j=j1;
      fin_j=j2;

      SetCursor(LoadCursor(0,IDC_ARROW));
      IMAGE_MODIFIEE=FALSE;
   }
   else 
      OnSauverSous();
   //}
} 

// *******************************************************************************
// ******************************* PALETTE ***************************************
// *******************************************************************************
// gestion de la palette
afx_msg void CImageWnd::OnNoirEtBlanc()
{
   if(TypePalette!=PALNB)
   {
      TypePalette=PALNB;
      MyPalette->SetPaletteEntries(0, 256, MyLogPaletteNB->palPalEntry );
      Mapping_Palette(MyLogPaletteNB);
      if (PALETTE_AFFICHEE) pPaletteWnd->SendMessage(WM_TEXTE1,0,0);
      Invalidate(FALSE);	
   }
   else
   {
      TypePalette=PALINIT;
      MyPalette->SetPaletteEntries(0, 256, MyLogPaletteInit->palPalEntry );
      Mapping_Palette(MyLogPaletteInit);
      if (PALETTE_AFFICHEE) pPaletteWnd->SendMessage(WM_TEXTE1,0,0);
      Invalidate(FALSE);
   }
   ((CMainFrame*) GetMDIFrame())->bar.OnUpdateCmdUI( (CFrameWnd*) GetMDIFrame(), TRUE);

}

afx_msg void CImageWnd::OnCouleur()
{
   if(TypePalette!=PALCOUL)
   {
      TypePalette=PALCOUL;
      MyPalette->SetPaletteEntries(0, 256, MyLogPaletteCOUL->palPalEntry );
      Mapping_Palette(MyLogPaletteCOUL);
      if (PALETTE_AFFICHEE) pPaletteWnd->SendMessage(WM_TEXTE1,0,0);
      Invalidate(FALSE); 
   }
   else
   {
      TypePalette=PALINIT;;
      MyPalette->SetPaletteEntries(0, 256, MyLogPaletteInit->palPalEntry );
      Mapping_Palette(MyLogPaletteInit);
      if (PALETTE_AFFICHEE) pPaletteWnd->SendMessage(WM_TEXTE1,0,0);
      Invalidate(FALSE);
   }
   ((CMainFrame*) GetMDIFrame())->bar.OnUpdateCmdUI( (CFrameWnd*) GetMDIFrame(), TRUE);

}

afx_msg void  CImageWnd::OnPaletteAffiche() 
{ 
   if (!PALETTE_AFFICHEE)                                                       
   {
      PALETTE_AFFICHEE = TRUE;
      pPaletteWnd = (CPaletteWnd * ) new CPaletteWnd(MyPalette,this);   
      pPaletteWnd->Create("",WS_POPUP|WS_SYSMENU| WS_VISIBLE | WS_CAPTION,this);   
      pPaletteWnd->SetWindowText("Palette de "+Titre);

   }
   else
   {
      Mapping_Palette(MyLogPaletteInit);
      MyPalette->SetPaletteEntries(0, 256, MyLogPaletteInit->palPalEntry );
      TypePalette=PALINIT;
      Invalidate(FALSE);	

      pPaletteWnd->DestroyWindow();
      PALETTE_AFFICHEE = FALSE;
   }
}

void CImageWnd::OnUpdateAffichePalette(CCmdUI* pCmdUI) 
{
   if (TypeColor == NB) pCmdUI->SetCheck(PALETTE_AFFICHEE);
   else pCmdUI->Enable(FALSE);
}

void CImageWnd::OnUpdateNb(CCmdUI* pCmdUI) 
{
   if (TypeColor == NB) pCmdUI->SetCheck(TypePalette==PALNB);
   else pCmdUI->Enable(FALSE);
}

void CImageWnd::OnUpdateCouleur(CCmdUI* pCmdUI) 
{
   if (TypeColor == NB) pCmdUI->SetCheck(TypePalette==PALCOUL);
   else pCmdUI->Enable(FALSE);
}

// *******************************************************************************
// *************************** FIN PALETTE ***************************************

// **********************************************************************
// Transforme les coordonnées écran en coordonnées image et son inverse
// **********************************************************************
CPoint CImageWnd::Get_ij(CPoint point)
{  
   point.x  = (int) (point.x)+PosVFen; 
   point.y  = (int) (Haut-1-(int)(point.y)+HautWndInfo) - PosHFen;
   return point;
} 

CPoint CImageWnd::Get_ij(int x,int y)
{ 
   CPoint point; 
   point.x  = (int) (x)+PosVFen; 
   point.y  = (int) (Haut-1-y+HautWndInfo) - PosHFen;
   return point;
} 

CPoint CImageWnd::Inv_Get_ij(int i,int j)
{
   CPoint point;
   point.x  = i-PosVFen;  
   point.y  = (int)Haut-1-j+HautWndInfo - PosHFen; 
   return point;
}


// *********************************************************************/
////////////////////////// Gestion du clipboard  ////////////////////////
// *********************************************************************/
void CImageWnd::OnEditCopy() 
{
   HBITMAP bitm;
   HDC hMemDC;

   CClientDC dc(this);
   VERIFY(::OpenClipboard(GetSafeHwnd()));
   bitm = CreateCompatibleBitmap(dc.GetSafeHdc(),fin_i-debut_i+1,fin_j-debut_j+1);
   hMemDC 	= CreateCompatibleDC(dc.GetSafeHdc());
   SelectObject(hMemDC, bitm);                 
   BitBlt(hMemDC,0,0,fin_i-debut_i+1,fin_j-debut_j+1,dc.GetSafeHdc(),Inv_Get_ij(debut_i,fin_j).x,Inv_Get_ij(debut_i,fin_j).y,SRCCOPY);

   VERIFY(::EmptyClipboard());
   VERIFY(::SetClipboardData(CF_BITMAP,bitm));
   VERIFY(::CloseClipboard());

   DeleteDC(hMemDC);
   DeleteObject(bitm);  
   ReleaseDC(&dc); // le dc est relaché
}

void CImageWnd::OnEditCut() 
{
   OnEditCopy();
   Memorise();
   for (long i=debut_i;i<=fin_i;i++) 	
      for(long j=debut_j;j<=fin_j;j++) 
         Image[i+j*Larg]=255;

   Invalidate(FALSE);
}

void CImageWnd::OnEditPaste() 
{
   HBITMAP bitm;
   HDC hMemDC;
   BITMAP bmp;

   SetCursor(LoadCursor(0,IDC_WAIT));

   CClientDC dc(this);
   VERIFY(::OpenClipboard(GetSafeHwnd()));
   VERIFY(bitm = (HBITMAP)::GetClipboardData(CF_BITMAP));

   hMemDC 	= CreateCompatibleDC(dc.GetSafeHdc());
   SelectObject(hMemDC, bitm);
   GetObject(bitm,sizeof(BITMAP),(LPSTR)&bmp);

   CTraitImageWnd * pImagePaste =  (CTraitImageWnd * ) new CTraitImageWnd(bmp.bmWidth,bmp.bmHeight,COLOR);
   GetDIBits( hMemDC, bitm, 0, bmp.bmHeight, pImagePaste->Image, pImagePaste->MyInitInfo, DIB_RGB_COLORS );
   pImagePaste->Create(" Image Presse-papier",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());


   VERIFY(::CloseClipboard());

   DeleteDC(hMemDC);
   DeleteObject(bitm);  
   ReleaseDC(&dc); // le dc est relaché  
   SetCursor(LoadCursor(0,IDC_ARROW));  
}

// *********************************************************************/
/////////////////////  Gestion de l'impression  /////////////////////////
// *********************************************************************/
void CImageWnd::OnFilePrint() 
{
   HDC    hdcPrn ;

   // Instantiate a CPrintDialog.
   CPrintDialog *printDlg =
      new CPrintDialog(FALSE, PD_ALLPAGES | PD_RETURNDC | PD_ENABLESETUPHOOK, NULL);

   // Initialize some of the fields in PRINTDLG structure.
   printDlg->m_pd.nMinPage = printDlg->m_pd.nMaxPage = 1;
   printDlg->m_pd.nFromPage = printDlg->m_pd.nToPage = 1;

   // Display Windows print dialog box.
   printDlg->DoModal();

   // Obtain a handle to the device context.
   hdcPrn = printDlg->GetPrinterDC();
   if (hdcPrn != NULL)
   {
      CDC *pDC = new CDC;
      pDC->Attach (hdcPrn);      // attach a printer DC

      pDC->StartDoc("test");     // begin a new print job

      SetPrintAlign(pDC, hdcPrn);// Set the printing alignment

      pDC->StartPage();          // begin a new page
      StretchDIBits(pDC->GetSafeHdc(),0,HautWndInfo, 2*Larg,2*Haut,PosVFen,-PosHFen,Larg,Haut,PImage,MyInitInfo,DIB_RGB_COLORS,SRCCOPY);

      pDC->EndPage();            // end a page

      pDC->EndDoc();             // end a print job

      pDC->Detach();             // detach the printer DC
      delete pDC;
   }

   delete printDlg;

}

void CImageWnd::SetPrintAlign(CDC *pDC, HDC hdcPrn)
{
   int nHsize = (int) (pDC->GetDeviceCaps(HORZSIZE) * 1700 / 254L);
   int nVsize = (int) (pDC->GetDeviceCaps(VERTSIZE) * 1700 / 254L);
   pDC->SetMapMode(MM_ANISOTROPIC);
   pDC->SetWindowExt(nHsize, nVsize);
   pDC->SetViewportExt(pDC->GetDeviceCaps(HORZRES),
      pDC->GetDeviceCaps(VERTRES));
   /*    short cxPage, cyPage, cxSize, cySize;

   cxPage = ::GetDeviceCaps (hdcPrn, HORZRES) ;
   cyPage = ::GetDeviceCaps (hdcPrn, VERTRES) ;
   cxSize = ::GetDeviceCaps (hdcPrn, PHYSICALWIDTH);
   cySize = ::GetDeviceCaps (hdcPrn, PHYSICALHEIGHT);
   pDC->SetMapMode (MM_ISOTROPIC) ;
   pDC->SetWindowExt ( cxSize, cySize) ;

   pDC->SetViewportExt (cxSize, cySize) ;
   pDC->SetViewportOrg (0,  cyPage / 2) ;

   //pDC->SetViewportExt (cxPage / 2, -cyPage / 2) ;
   //pDC->SetViewportOrg (cxPage / 2,  cyPage / 2) ;
   pDC->SetTextAlign (TA_BASELINE | TA_CENTER) ;*/
}

// ***********************************************************
// Restaure la zone rectangle sélectionnée
// ***********************************************************
void CImageWnd::OnRestaureZone()
{
   debut_i = 0;fin_j = Haut-1;   
   debut_j = 0;fin_i = Larg-1;  
   // remplissage de la structure zone
   Im.Lmin		= debut_i;
   Im.Lmax		= fin_i;
   Im.Hmin		= debut_j;
   Im.Hmax		= fin_j;   
   zoneSelectionnee = FALSE;  
   Invalidate(FALSE);
}  

afx_msg void CImageWnd::OnRestaureImage()
{   
   long l;
   switch (TypeColor)
   {   
   case NB:
      for (l=0;l<Larg*Haut;l++) ( (unsigned char  *)PImage)[l] = ( (unsigned char  *) PImageTemp)[l];
      break;
   case COLOR:
      for (l=0;l<Larg*Haut;l++) ( (rgb  *)PImage)[l] = ( (rgb  *) PImageTemp)[l];	               		
      break;
   } 
   NomTraitement = "Aucun";
   CommentTraitement = "Image restaurée";
   dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
   Invalidate(FALSE);
}


afx_msg void CImageWnd::OnSelectionneZoneClavier()
{    
   CDialZone dial(&debut_i,&debut_j,&fin_i,&fin_j,this);   
   if (dial.DoModal()==IDOK) 
   {      
      TestMinMax(debut_i,debut_j,fin_i,fin_j); 
      Troncature();
      Im.Lmin		= debut_i;
      Im.Lmax		= fin_i;
      Im.Hmin		= debut_j;
      Im.Hmax		= fin_j;  
      if((debut_i!=0)||(debut_j!=0)||(fin_i!=Larg-1)||(fin_j!=Haut-1))zoneSelectionnee=TRUE;  
   }
} 


// **********************************************************
// Mets des Booléens pour gestion du rectangle par la souris
// **********************************************************
afx_msg void CImageWnd::OnSelectionneZone()
{    
   ArreteLeContinu();
   if (zoneSelectionnee)
   {
      LaZone = FALSE;
      Appuie  = FALSE;
      first = FALSE;
      zoneSelectionnee = FALSE;
      OnRestaureZone();
   }
   else
   {
      LaZone = TRUE;
      Appuie  = FALSE;
      first = TRUE;
      zoneSelectionnee = TRUE;
   }
}

void CImageWnd::OnUpdateZone(CCmdUI* pCmdUI) 
{
   pCmdUI->SetCheck(zoneSelectionnee);
}

// **************************************************
// Coupe horizontale
// **************************************************
afx_msg void CImageWnd::OnSelectionneCoupeH()
{
   ArreteLeContinu();
   LaCoupe = TRUE;
   Appuie  = FALSE;
   first = TRUE;
   CoupeSelectionnee = FALSE;
   COUPE_VER = FALSE;COUPE_HOR = TRUE;COUPE_OBLI = FALSE;     
   AfficheCoupe = TRUE;
}

// **********************************************
// Coupe verticale
// **********************************************
afx_msg void CImageWnd::OnSelectionneCoupeV()
{    
   ArreteLeContinu();
   LaCoupe = TRUE;
   Appuie  = FALSE;
   first = TRUE;
   CoupeSelectionnee = FALSE;
   COUPE_VER = TRUE;COUPE_HOR = FALSE;COUPE_OBLI = FALSE;          
   AfficheCoupe = TRUE;
}

// ****************************************************
// Coupe oblique
// ***************************************************
afx_msg void CImageWnd::OnSelectionneCoupeOblique()
{    
   ArreteLeContinu();
   LaCoupe = TRUE;
   Appuie  = FALSE;
   first = TRUE;
   CoupeSelectionnee = FALSE;
   COUPE_VER = FALSE;COUPE_HOR = FALSE;COUPE_OBLI = TRUE;      
   AfficheCoupe = TRUE;
}


afx_msg void CImageWnd::OnHistogramme()
{
   long TabHisto[256];
   Histo(TabHisto);   
   CGraphWnd * pHisto = (CGraphWnd*) new CGraphWnd(TabHisto,256,"Nombre","Niveau");  
   if (!pHisto->Create("Histogramme : "+Titre,GetMDIFrame())) return;
}

afx_msg void CImageWnd::OnFonctionRepartition()
{    
   long TabFRepar[256];    
   long TabHisto[256];   
   Histo(TabHisto);  
   FonctionRepartition(TabFRepar,TabHisto);  
   CGraphWnd * pFRepart = (CGraphWnd*) new CGraphWnd(TabFRepar,256,"Nombre","Niveau");
   if (!pFRepart->Create("Fonction de répartition : "+Titre,GetMDIFrame())) return; 
}


// ************************************************
// Routine d'information moyenne, écart type, mini, maxi sur zone                       
// **********************************************
afx_msg void CImageWnd::OnInformationGenerale()
{
   char chaine[100];
   SetCursor(LoadCursor(0,IDC_WAIT)); 

   long l,cpt=0;
   float n,sig=0.0F;
   float Moy,moy=0.0F,Sig;
   float Min,Max;
   if (TypeImage == MEM || TypeImage == BMP) Min = Max = (float) Image[0];
   if (TypeImage == FLOAT32) Min = Max = (float) ImageFloat[0];

   for (long j=debut_j;j<=fin_j;j++)
   { 
      l = j*Larg+debut_i;
      for (long i=debut_i;i<=fin_i;i++)
      {  
         if (TypeImage == MEM || TypeImage == BMP) n = (float) Image[l++];
         if (TypeImage == FLOAT32) n = (float) ImageFloat[l++];
         moy += n; 
         sig += n*n;  
         if (n>Max) Max = n; 
         if (n<Min) Min = n;
         cpt ++;
      } 
   } 
   Moy = (float) moy / (float)cpt;
   Sig = (float) sig / (float)cpt;
   Sig = (float) sqrt((double)(Sig-Moy*Moy));

   SetCursor(LoadCursor(0,IDC_ARROW));
   if (TypeImage == MEM || TypeImage == BMP) 	sprintf(chaine,"niveau moyen : %3.1f\nécart type : %3.1f\nmin : %d\nmax : %d",Moy,Sig,(unsigned char)Min,(unsigned char)Max);
   if (TypeImage == FLOAT32) 	sprintf(chaine,"moyenne : %3.1f\nécart type : %3.1f\nmin : %3.1f\nmax : %3.1f",Moy,Sig,Min,Max);
   MessageBox(chaine,"Information générale");

}



// ****************************************************************************
// 							gestion de la souris                         
// ****************************************************************************
afx_msg void CImageWnd::OnLButtonDown(UINT nFlags,CPoint Point) 
{
   CString ChaineTemp ="";
   if (ZOOM==1)
      if ((LaZone) || (LaCoupe))
      { 

         C1.x =  Point.x;
         C1.y =  Point.y;    
         C2.x =  Point.x;
         C2.y =  Point.y;
         Appuie = TRUE;
      }
      if(zoneSelectionnee && !LaZone)
      {
         int l=fin_i-debut_i;
         int h=fin_j-debut_j;
         /*debut_i = Get_ij(Point).x;
         fin_i   = Get_ij(Point).x+l;  
         debut_j = Get_ij(Point).y;
         fin_j   = Get_ij(Point).y+h;
         */	C1.x=Point.x;
         C1.y=Point.y;
         C2.x=Point.x+l;
         C2.y=Point.y-h;
         Appuie=TRUE;
         first=TRUE;
      }

      if(MATCHINGCURVI)
      {

         if(CptManipMatching<3)
         {
            ChaineTemp.Format("%d",Point.x);
            FichierMatching.WriteString(ChaineTemp);
            FichierMatching.WriteString("\n");
            ChaineTemp.Format("%d",Haut-Point.y+19);//car inversé dans l'espace de l'image et largeur bordure
            FichierMatching.WriteString(ChaineTemp);
            FichierMatching.WriteString("\n");
            if(CptManipMatching==2)
            {
               FichierMatching.Close();
               MATCHINGCURVI=FALSE;
            }
         }
         CptManipMatching++;
      }

      if(MATCHINGCURVI4)
      {

         if(CptManipMatching<4)
         {
            ChaineTemp.Format("%d",Point.x);
            FichierMatching4.WriteString(ChaineTemp);
            FichierMatching4.WriteString("\n");
            ChaineTemp.Format("%d",Haut-Point.y+19);//car inversé dans l'espace de l'image et largeur bordure
            FichierMatching4.WriteString(ChaineTemp);
            FichierMatching4.WriteString("\n");
            if(CptManipMatching==3)
            {
               FichierMatching4.Close();
               MATCHINGCURVI4=FALSE;
            }
         }
         CptManipMatching++;
      }
}      

afx_msg void CImageWnd::OnLButtonUp(UINT nFlags,CPoint Point)
{    
   if (ZOOM==1)
      if (Appuie) 
      {  

         Appuie = FALSE; 
         if (LaZone)  SourisRelachePourZone();
         if (LaCoupe) SourisRelachePourCoupe();
         if(zoneSelectionnee && !LaZone) SourisRelachePourBougeZone();
      }
} 

afx_msg void CImageWnd::OnMouseMove(UINT nFlags,CPoint Point)
{ 
   if (ZOOM==1)
   {
      AfficheNiveau(Point); 
      if((LaZone)||(LaCoupe)) SetCursor(LoadCursor(0,IDC_CROSS));
      else SetCursor(LoadCursor(0,IDC_ARROW));
      if (Appuie)
      {  
         if (LaZone) SourisBougePourZone(Point);   
         if (LaCoupe) SourisBougePourCoupe(Point); 
         if(zoneSelectionnee && !LaZone) SourisBougePourBougeZone(Point);
         first = FALSE;
      } 
   }	
   if(MATCHINGCURVI)
   {
      SetCursor(LoadCursor(0,IDC_CROSS));
   }

   if(MATCHINGCURVI4)
   {
      SetCursor(LoadCursor(0,IDC_CROSS));
   }
}


void CImageWnd::MettreAJourLaCoupeActive()
{   
   Pixel TabIJ[1000];
   RemplitTableauIJCoupe(TabIJ);          

   void * GraphPoint;                         
   if (TypeImage == MEM || TypeImage == BMP)
   {
      GraphPoint = (int *) new int[NbPointCoupe*sizeof(int)];  
      for (long l=0;l<NbPointCoupe;l++) ((int *)GraphPoint)[l] = (int) Image[(long)TabIJ[l].j*Larg+(long)TabIJ[l].i];
      pCoupeVariable->SetValeurTab((int *)GraphPoint); 
   }
   if (TypeImage==FLOAT32)
   {

      GraphPoint = (float *) new float[NbPointCoupe*sizeof(float)];  
      for (long l=0;l<NbPointCoupe;l++) ((float *)GraphPoint)[l] =  ImageFloat[(long)TabIJ[l].j*Larg+(long)TabIJ[l].i];
      pCoupeVariable->SetValeurTab((float *)GraphPoint); 
   }    

   pCoupeVariable->SetPositionCoupe(Lg_i1,Lg_j1,Lg_i2,Lg_j2); 
   pCoupeVariable->InvalidateRect(pCoupeVariable->rcInflate,TRUE);
   delete GraphPoint;
}

afx_msg void CImageWnd::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{
   if ((CoupeSelectionnee) && (pCoupeVariable!=NULL))
   {
      if ((Lg_i1==Lg_i2) || (Lg_j1==Lg_j2))
      {  
         CPoint p1((int)Lg_i1,(int)Lg_j1),p2((int)Lg_i2,(int)Lg_j2);

         if 	((nChar==39) 	&& (Lg_i2<Larg-1))	{ Lg_i1++; Lg_i2++;	}
         if	((nChar==37)  && (Lg_i1>1)) 			{ Lg_i1 --; Lg_i2--;}
         if	((nChar==38) 	&& (Lg_j1<Haut-1)) 	{	Lg_j1++; Lg_j2++; }
         if 	((nChar==40)  && (Lg_j2>1)) 			{	Lg_j1--; Lg_j2--; }    
         CPoint Point1,Point2;
         Point1 = Get_ij(Inv_Get_ij((int)Lg_i1,(int)Lg_j1)); 
         Point2 = Get_ij(Inv_Get_ij((int)Lg_i2,(int)Lg_j2));  
         if ((Point1.x!=(int)Lg_i1) || (Point1.y!=(int)Lg_j1)  || (Point2.x!=(int)Lg_i2) || (Point2.y!=(int)Lg_j2) )   
         {
            Lg_i1 = (long) p1.x; 
            Lg_j1 = (long) p1.y; 
            Lg_i2 = (long) p2.x; 
            Lg_j2 = (long) p2.y; 
         }
         else
         {                   
            MettreAJourLaCoupeActive();
            Invalidate(FALSE);
         }	
      }	
   }	
}
// ************** FIN DE GESTION DE LA SOURIS *******************************************
// **************************************************************************************




//  TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
void CImageWnd::OnSeuillerImage1()
{               
   ArreteLeContinu();
   BOOL flag;
   CDialSeuil dial (&Seuil,&flag); 
   if (dial.DoModal()==IDOK)
   {
      SetCursor(LoadCursor(0,IDC_WAIT));
      Seuillage1(Im,Seuil,flag);      
      SetCursor(LoadCursor(0,IDC_ARROW));
      Invalidate(FALSE);
   }
}    

void CImageWnd::OnSeuillerImage2()
{               
   ArreteLeContinu();
   CDialDoubleSeuil  dial(IDD_SEUIL2,&Seuil,&Seuil2); 
   if (dial.DoModal()==IDOK)
   {
      SetCursor(LoadCursor(0,IDC_WAIT));
      Seuillage2(Im,Seuil,Seuil2);      
      SetCursor(LoadCursor(0,IDC_ARROW));
      Invalidate(FALSE);
   }
}     


// fin routine de traitement



// *********************************************************************************
// µµµµµµµµµµµµµµµµµµµµµµµµµµµµ   ROUTINES PRIVEES  µµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµ
// *********************************************************************************
void CImageWnd::DessineLigne(long i1, long j1, long i2,long j2)
{     
   CClientDC DC(this);     
   CPen stylo(PS_SOLID,1,RGB(255,0,0));
   DC.SelectObject(&stylo); 
   CPoint h1,h2;
   h1 = Inv_Get_ij((int)i1,(int)j1); 
   h2 = Inv_Get_ij((int)i2,(int)j2);
   DC.MoveTo(h1.x,h1.y);
   DC.LineTo(h2.x,h2.y);
   ReleaseDC(&DC);   
   TRACER_COUPE = TRUE;   
} 

void CImageWnd::DessineLigne(long i1, long j1, long i2,long j2, COLORREF Color)
{     
   //   if (!TRACER_COUPE)
   { 
      CClientDC DC(this);     
      CPen stylo(PS_SOLID,1,Color);
      DC.SelectObject(&stylo); 
      CPoint h1,h2;
      h1 = Inv_Get_ij((int)i1,(int)j1); 
      h2 = Inv_Get_ij((int)i2,(int)j2);
      DC.MoveTo(h1.x,h1.y);
      DC.LineTo(h2.x,h2.y);
      ReleaseDC(&DC);   
      TRACER_COUPE = TRUE;   
   }  
} 

void CImageWnd::DessineRectangle(int i1, int j1, int i2,int j2, COLORREF Color)
{     
   CClientDC DC(this);     
   CPen stylo(PS_SOLID,1,Color);
   DC.SelectObject(&stylo);   
   CPoint h1,h2;        
   h1 = Inv_Get_ij(i1,j1); 
   h2 = Inv_Get_ij(i2,j2); 
   DC.SelectStockObject(NULL_BRUSH);
   DC.Rectangle(CRect(h1, h2));
   ReleaseDC(&DC);   
} 


void CImageWnd::DessineRectangle()
{     
   CClientDC DC(this);     
   CPen stylo(PS_SOLID,1,RGB(0,0,255));
   DC.SelectObject(&stylo);
   CPoint h1,h2;
   h1 = Inv_Get_ij((int)debut_i-1,(int)debut_j-2); 
   h2 = Inv_Get_ij((int)fin_i+2,(int)fin_j+1);
   DC.SelectStockObject(NULL_BRUSH);
   DC.Rectangle(CRect(h1, h2));
   ReleaseDC(&DC);
} 


// ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// test du min et max des limites du rectangle
// ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
void CImageWnd::TestMinMax(long& i1,long& j1,long& i2,long& j2)
{  
   long temp;
   if (i1>i2) 
   { 
      temp =  i2;
      i2 = i1;
      i1 = temp;
   }  
   if (j1>j2) 
   { 
      temp =  j2;
      j2 = j1;
      j1 = temp;
   }             
}

void CImageWnd::Troncature()
{
   if (debut_i<0) debut_i = 0;
   if (debut_j<0) debut_j = 0; 
   if (fin_j>=Haut) fin_j = Haut-1;
   if (fin_i>=Larg) fin_i = Larg-1;   
   if (fin_i<0) fin_i = 0;
   if (fin_j<0) fin_j = 0; 
   if (debut_j>=Haut) debut_j = Haut-1;
   if (debut_i>=Larg) debut_i = Larg-1; 
}

// ***************************** LA ZONE *********************************************
//  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// Routine appelée lorsque la zone et BoutonUp de la souris
// ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
void  CImageWnd::SourisRelachePourZone()
{
   if (C1!=C2)
   { 
      zoneSelectionnee = TRUE;   
      debut_i = Get_ij(C1).x;
      fin_i   = Get_ij(C2).x;  
      debut_j = Get_ij(C1).y;
      fin_j   = Get_ij(C2).y;    
      TestMinMax(debut_i,debut_j,fin_i,fin_j); 
      Troncature();  
      div_t res = div((int)(fin_i-debut_i+1),4);   
      if (res.rem !=0) fin_i -= res.rem;   
      res = div((int)(fin_j-debut_j+1),4);   
      if (res.rem !=0) fin_j -= res.rem;
      // remplissage de la structure zone
      Im.Lmin		= debut_i;
      Im.Lmax		= fin_i;
      Im.Hmin		= debut_j;
      Im.Hmax		= fin_j;  
      Invalidate(TRUE);  
   } 
   LaZone = FALSE;
}

void  CImageWnd::SourisRelachePourBougeZone()
{ 
   debut_i = Get_ij(C1).x;
   fin_i   = Get_ij(C2).x;  
   debut_j = Get_ij(C1).y;
   fin_j   = Get_ij(C2).y;    
   TestMinMax(debut_i,debut_j,fin_i,fin_j); 
   Troncature();  
   // remplissage de la structure zone
   Im.Lmin		= debut_i;
   Im.Lmax		= fin_i;
   Im.Hmin		= debut_j;
   Im.Hmax		= fin_j;  
   Invalidate(TRUE);  
}  
//  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// Routine appelée lorsque la zone et Move de la souris
// ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
void CImageWnd::SourisBougePourZone(CPoint Point)     
{
   if (C1!=Point)
   { 
      CClientDC dc(this);         
      dc.SetROP2(R2_XORPEN);
      CPen stylo(PS_SOLID,1,RGB(0,0,255));
      dc.SelectObject(&stylo); 
      dc.SelectStockObject(NULL_BRUSH);
      if (!first)
      {    
         dc.Rectangle(CRect(C1, C2));
      }
      dc.Rectangle(CRect(C1, Point));
      ReleaseDC(&dc); 
      C2 = Point; 
      debut_i = Get_ij(C1).x;
      fin_i   = Get_ij(C2).x;  
      debut_j = Get_ij(C1).y;
      fin_j   = Get_ij(C2).y; 
      TestMinMax(debut_i,debut_j,fin_i,fin_j);         
      AfficheDimensions();
   }  
}

void CImageWnd::SourisBougePourBougeZone(CPoint Point)     
{

   CClientDC dc(this);         
   dc.SetROP2(R2_XORPEN);
   CPen stylo(PS_SOLID,1,RGB(0,0,255));
   dc.SelectObject(&stylo); 
   dc.SelectStockObject(NULL_BRUSH);
   if(first)dc.Rectangle(CRect(Inv_Get_ij(debut_i-1,debut_j-2), Inv_Get_ij(fin_i+2,fin_j+1)));
   else dc.Rectangle(CRect(C1,C2));

   int l=fin_i-debut_i;
   int h=fin_j-debut_j;
   C1.x=Point.x;
   C1.y=Point.y;
   C2.x=Point.x+l;
   C2.y=Point.y-h;
   dc.Rectangle(CRect(C1,C2));
   ReleaseDC(&dc);
   debut_i = Get_ij(C1).x;
   fin_i   = Get_ij(C2).x;  
   debut_j = Get_ij(C1).y;
   fin_j   = Get_ij(C2).y; 
   TestMinMax(debut_i,debut_j,fin_i,fin_j);         
   AfficheDimensions(); 
}

// ******************************* FIN ZONE *****************************************************



// ***************************** LA COUPE *********************************************
//  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// Routine appelée lorsque la Coupe et BoutonUp de la souris
// ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
void CImageWnd::SourisBougePourCoupe(CPoint Point)     
{
   if (C1!=Point)
   { 
      CClientDC dc(this);         
      dc.SetROP2(R2_XORPEN);
      CPen stylo(PS_SOLID,1,RGB(255,0,0));
      dc.SelectObject(&stylo); 
      if (!first)
      {    
         dc.MoveTo(C1.x,C1.y);                                  
         if (COUPE_VER) dc.LineTo(C1.x,C2.y);
         else if (COUPE_HOR) dc.LineTo(C2.x,C1.y);
         else dc.LineTo(C2.x,C2.y);
      }
      if (COUPE_VER) {C1.x = Point.x;}
      if (COUPE_HOR) {C1.y = Point.y;}   
      dc.MoveTo(C1.x,C1.y);
      if (COUPE_VER) dc.LineTo(C1.x,Point.y);
      else if (COUPE_HOR) dc.LineTo(Point.x,C1.y);
      else dc.LineTo(Point.x,Point.y); 
      ReleaseDC(&dc); 
      C2 = Point;
      AffichePositionCoupe();
   }

}


//  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// Routine appelée lorsque la zone et Move de la souris
// ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
void  CImageWnd::SourisRelachePourCoupe()
{   
   LaCoupe = FALSE; 
   if (C1!=C2)
   {  
      CoupeSelectionnee = TRUE;   
      if(!COUPE_OBLI)TestMinMax(C1.x,C1.y,C2.x,C2.y);
      Lg_i1 = Get_ij(C1).x;
      Lg_i2 = Get_ij(C2).x;  
      Lg_j1 = Get_ij(C1).y;
      Lg_j2 = Get_ij(C2).y;   
      if(Lg_i1>Larg || Lg_i2>Larg || Lg_j1>Haut||Lg_j2>Haut||Lg_i1<0|| Lg_i2<0 || Lg_j1<0||Lg_j2<0)
      {
         MessageBox("Coupe en dehors de l'image!","Erreur!",MB_OK|MB_ICONSTOP);
         CoupeSelectionnee=FALSE;
         Invalidate(TRUE);
         return;
      }
      DessineLigne(Lg_i1,Lg_j1,Lg_i2,Lg_j2);
      Pixel TabIJ[1000];
      RemplitTableauIJCoupe(TabIJ);
      void * GraphPoint;                         
      CProfilWnd  * pCoupe;
      if (TypeImage == MEM || TypeImage == BMP)
      {
         GraphPoint = (int *) new int[NbPointCoupe*sizeof(int)];  
         for (long l=0;l<NbPointCoupe;l++) ((int *)GraphPoint)[l] = (int) Image[(long)TabIJ[l].j*Larg+(long)TabIJ[l].i];
         if (AfficheCoupe) pCoupe = (CProfilWnd*)(new CProfilWnd((int *)GraphPoint,(int)NbPointCoupe, "Niveau","Position",this,Get_ij(C1),Get_ij(C2),TRUE));
      }
      if (TypeImage==FLOAT32)
      {
         GraphPoint = (float *) new float[NbPointCoupe*sizeof(float)];  
         for (long l=0;l<NbPointCoupe;l++) ((float *)GraphPoint)[l] =  ImageFloat[(long)TabIJ[l].j*Larg+(long)TabIJ[l].i];
         if (AfficheCoupe) pCoupe = (CProfilWnd*)(new CProfilWnd((float *)GraphPoint,(int)NbPointCoupe, "Niveau","Position",this,Get_ij(C1),Get_ij(C2),TRUE));
      }    
      delete GraphPoint;
      if (AfficheCoupe)
      {	  
         if(COUPE_VER)  if (!pCoupe->Create("Coupe Verticale : "+Titre,GetMDIFrame())) return;
         if(COUPE_HOR)  if (!pCoupe->Create("Coupe Horizontale : "+Titre, GetMDIFrame())) return;
         if(COUPE_OBLI) if (!pCoupe->Create("Coupe Oblique : "+Titre, GetMDIFrame())) return;
         AfficheCoupe = FALSE;
         NbProfil++;
      }         
   }
}



void CImageWnd::RemplitTableauIJCoupe(Pixel TabIJ[1000])
{
   float Pente=0.0F, OrdonneeOrigine=0.0F;
   char equation='0';
   long signe;
   long i,j,l;

   i = Lg_i1;
   j = Lg_j1;          

   TabIJ[0].i = (unsigned short) i;
   TabIJ[0].j = (unsigned short) j;

   if ((Lg_i2-Lg_i1)==0) equation = '1';   //equation de la forme x=cste
   else if((Lg_j2-Lg_j1)==0) equation = '2';   //equation de la forme y=cste
   else equation = '3';   //equation de la forme y=ax+b
   switch (equation)
   {  
   case '1':       
      if (Lg_j1>Lg_j2) signe=-1;else signe = 1;  
      NbPointCoupe = labs(Lg_j2-Lg_j1)+1;         
      for (l=1; l<NbPointCoupe; l++)
      {   j +=signe;
      TabIJ[l].i = (unsigned short) i;
      TabIJ[l].j = (unsigned short) j;
      }
      break;
   case '2':
      if (Lg_i1>Lg_i2) signe=-1;else signe = 1;    
      NbPointCoupe = labs(Lg_i2-Lg_i1)+1;
      for (l=1; l<NbPointCoupe; l++)
      {   
         i += signe;
         TabIJ[l].i = (unsigned short) i;
         TabIJ[l].j = (unsigned short) j;
      }
      break;
   case '3':
      if (!((float)Lg_i2-(float)Lg_i1))
         break;
      Pente = (float)(Lg_j2-Lg_j1)/((float)Lg_i2-(float)Lg_i1); //initialise la pente a
      OrdonneeOrigine = (float)Lg_j2 - Pente*(float)Lg_i2;  //initialise l'ordonnee b de y=a*x+b
      if (Pente>0.0)
      {

         if (Pente>1.0)	
         {
            if (Lg_j1>Lg_j2) signe=-1;else signe = 1; 
            NbPointCoupe = labs(Lg_j2-Lg_j1)+1;
            for (l=1; l<NbPointCoupe; l++) 
            {   
               j+=signe;   
               i = (long)  (((float)j*1.0-OrdonneeOrigine)/Pente);   
               TabIJ[l].i = (unsigned short) i;
               TabIJ[l].j = (unsigned short) j;

            }	
         }
         else
         {	
            if (Lg_i1>Lg_i2) signe=-1;else signe = 1; 
            NbPointCoupe = labs(Lg_i2-Lg_i1)+1;
            for (l=1; l<NbPointCoupe; l++)         
            {
               i +=signe ;
               j = (long)(Pente*(float)i+OrdonneeOrigine);
               TabIJ[l].i = (unsigned short) i;
               TabIJ[l].j = (unsigned short) j;
            }                                   
         }	
      }
      else
      {   
         if (Pente<-1.0)
         {
            if (Lg_j1>Lg_j2) signe=-1;else signe = 1;  
            NbPointCoupe = labs(Lg_j2-Lg_j1)+1;
            for (l=1; l<NbPointCoupe; l++)   
            {  
               j+=signe;   
               i = (long) (((float)j*1.0-OrdonneeOrigine)/Pente);   
               TabIJ[l].i = (unsigned short) i;
               TabIJ[l].j = (unsigned short) j;
            }	
         }
         else
         {
            if (Lg_i1>Lg_i2) signe=-1;else signe = 1; 
            NbPointCoupe = labs(Lg_i2-Lg_i1)+1;
            for (l=1; l<NbPointCoupe; l++)     
            {
               i +=signe ;
               j = (long)(Pente*(float)i+OrdonneeOrigine);
               TabIJ[l].i = (unsigned short) i;
               TabIJ[l].j = (unsigned short) j;
            }	
         }

      }	
      break;

   default:
      MessageBox ("Pas d'initialisation","",MB_OK);
      break;
   }
}

void  CImageWnd::AffichePositionCoupe()
{
   char titre[100];
   sprintf(titre," coupe entre (%3ld,%3ld) et (%3ld, %3ld)",Lg_i1,Lg_j1,Lg_i2,Lg_j2);   
   CStatusBar *pStatus = (CStatusBar *)
      (GetMDIFrame()->GetDescendantWindow(AFX_IDW_STATUS_BAR));
   if (pStatus)
   {
      wsprintf(titre,"Coupe entre (%3ld,%3ld) et (%3ld, %3ld)",C1.x, C1.y, C2.x, C2.y);
      pStatus->SetPaneText(1, titre);	// premier panneau de la barre d'etat
   }
}


// ********************************* FIN COUPE **********************************************



void  CImageWnd::AfficheDimensions()
{
   char titre[40];
   sprintf(titre,"LxH (%3ldx%3ld)",Larg/ZOOM,Haut/ZOOM);   
   m_wndStatusBarImage.SetPaneText(0, titre);
   char text[100];
   CStatusBar *pStatus = (CStatusBar *)
      (GetMDIFrame()->GetDescendantWindow(AFX_IDW_STATUS_BAR));
   if (pStatus)
   {
      wsprintf(text,"Dim ZONE = (%d, %d)",fin_i-debut_i+1,fin_j-debut_j+1);
      pStatus->SetPaneText(1, text);	// premier panneau de la barre d'etat
   }
}

void  CImageWnd::Histo( long TabHisto[256])
{
   long i;
   long l;
   for (i=0;i<256;i++) TabHisto[i] = 0;  

   for (long j=debut_j;j<=fin_j;j++)
   { 
      l = j*Larg+debut_i;
      for (i=debut_i;i<=fin_i;i++) TabHisto[Image[l++]]++;
   } 
}

void  CImageWnd::FonctionRepartition(long TabFRepar[256],long TabHisto[256])
{ 
   long i;
   for (i=0;i<256;i++) TabFRepar[i] = 0;
   TabFRepar[0] = TabHisto[0];   
   for (i=1;i<256;i++) TabFRepar[i] = TabFRepar[i-1] + TabHisto[i]; 
}                  


// *******************************************
// routine d'affichage du niveau de gris et i,j
void CImageWnd::AfficheNiveau(CPoint Point)
{   
   char chaine[50];          
   I = Get_ij(Point).x;
   J = Get_ij(Point).y;  
   if (I<0) I = 0;
   if (J<0) J = 0; 
   if (J>=(int)Haut) J = (int)Haut-1;
   if (I>=(int)Larg) I = (int)Larg-1;
   if (TypeImage == MEM || TypeImage == BMP)
   {       
      if (TypeColor == NB)
      {  
         unsigned char  niv;
         niv = Image[J*Larg+I]; 
         sprintf(chaine,"NB de (%3d,%3d) = %3d",I,J,niv);
      }	
      else
      {   
         rgb RGB;	
         RGB = ImageRGB[J*Larg+I]; 
         sprintf(chaine,"RGB de (%3d,%3d) = %3d,%3d,%3d",I,J,RGB.r,RGB.g,RGB.b);
      }
   }
   if (TypeImage == FLOAT32)
   { 
      float v;
      v = ImageFloat[J*Larg+I]; 
      sprintf(chaine,"float de (%3d,%3d) = %f",I,J,v);
   }
   m_wndStatusBarImage.SetPaneText(1,chaine);
}   

// ******************************LA PALETTE******************************
void CImageWnd::Mapping_Palette(LOGPALETTE * log)
{    
   if (TypeColor == NB)
      for (int i = 0; i<256; i++)
      {
         MyInitInfo->bmiColors[i].rgbBlue     =   log->palPalEntry[i].peBlue;
         MyInitInfo->bmiColors[i].rgbGreen    =   log->palPalEntry[i].peGreen;
         MyInitInfo->bmiColors[i].rgbRed      =   log->palPalEntry[i].peRed;
         MyInitInfo->bmiColors[i].rgbReserved =   log->palPalEntry[i].peFlags;
      }

}	

void CImageWnd::Palette_Couleurs()
{
   int vert,bleu,rouge;
   for (int l=0;l<256;l++)
   {
      if ((l>=0)&&(l<=33))
      {
         vert  = (int)0;
         bleu  = (int)(255*(0-l)/(0-33));
         rouge = (int)(255*(0-l)/(0-33)); 
      }else
         if ((l>=34)&&(l<=67))
         {
            vert = 0;
            bleu = 255;
            rouge = (int)(255*(77-l)/(77-34));

         }
         else   
            if ((l>=68)&&(l<=100) )
            {
               vert  = (int)(255*(68-l)/(68-100));
               bleu  = (int)255;
               rouge = (int)0;
            }else 
               if( (l>=101)&&(l<=134) )
               {
                  vert = 255;
                  bleu = (int)(255.0*(float)(134-l)/(float)(134-101));
                  rouge = 0;
               }else
                  if ((l>=135)&&(l<=177))
                  { 
                     vert  = 255;
                     bleu  = 0;
                     rouge = (int)(255.0*(float)(135-l)/(float)(135-177));
                  }else                 
                     if ((l>=178)&&(l<=201))
                     {
                        vert = (int)(255.0*(float)(201-l)/(float)(201-178));
                        bleu = 0;
                        rouge =255;
                     }else	  
                        if ((l>=202)&&(l<=255))
                        {
                           vert = (int)(255.0*(float)(202-l)/(float)(202-255));
                           bleu = (int)(255.0*(float)(202-l)/(float)(202-255));
                           rouge = 255;
                        }
                        MyLogPaletteCOUL->palPalEntry[l].peRed   =  rouge;
                        MyLogPaletteCOUL->palPalEntry[l].peGreen =  vert;
                        MyLogPaletteCOUL->palPalEntry[l].peBlue  =  bleu;
                        MyLogPaletteCOUL->palPalEntry[l].peFlags = PC_NOCOLLAPSE;
   }
}



void CImageWnd::Palette_NB()   
{
   for (int l=0;l<255;l++)
   {
      MyLogPaletteNB->palPalEntry[l].peRed=l;
      MyLogPaletteNB->palPalEntry[l].peGreen=l;
      MyLogPaletteNB->palPalEntry[l].peBlue=l;
      MyLogPaletteNB->palPalEntry[l].peFlags=PC_NOCOLLAPSE;
   }  
   MyLogPaletteNB->palPalEntry[255].peRed=255;
   MyLogPaletteNB->palPalEntry[255].peGreen=255;
   MyLogPaletteNB->palPalEntry[255].peBlue=255;
   MyLogPaletteNB->palPalEntry[255].peFlags=PC_NOCOLLAPSE;
}
// **************************** FIN LA PALETTE*********************************



//  TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
//               Routines de traitement 
//  TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
//  TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Inversion vidéo
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void CImageWnd::OnInversionVideo()
{ 
   ArreteLeContinu();
   SetCursor(LoadCursor(0,IDC_WAIT)); 
   Memorise(); 
   Inversion(GetSImage());
   SetCursor(LoadCursor(0,IDC_ARROW));
   Invalidate(FALSE);                                           
}       

// MIROIR VERTICAL 
afx_msg void CImageWnd:: OnMiroirVert()
{      
   ArreteLeContinu();
   SetCursor(LoadCursor(0,IDC_WAIT));   
   Memorise();  
   MiroirV(GetSImage());
   SetCursor(LoadCursor(0,IDC_ARROW));  
   Invalidate(FALSE);
   NomTraitement = "Miroir Vertical";
   CommentTraitement = "Symetrie par rapport a un plan vertical";
   dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
}


// MIROIR HORIZONTAL 
afx_msg	void CImageWnd::OnMiroirH()
{
   ArreteLeContinu();
   SetCursor(LoadCursor(0,IDC_WAIT)); 
   Memorise();  
   MiroirH(GetSImage());
   SetCursor(LoadCursor(0,IDC_ARROW));
   Invalidate(FALSE); 
   NomTraitement = "Miroir Horizontal";
   CommentTraitement = "Symetrie par rapport a un plan horizontal";
   dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
} 


afx_msg void CImageWnd::OnRoberts()
{   
   ArreteLeContinu();
   CDialDoubleSeuil  dial(IDD_SEUIL2,&Seuil,&Seuil2);  
   if (dial.DoModal()==IDOK)
   {
      SetCursor(LoadCursor(0,IDC_WAIT)); 
      Memorise();  
      Roberts(GetSImage(),Seuil,Seuil2); 
      SetCursor(LoadCursor(0,IDC_ARROW)); 
      Invalidate(FALSE);
      NomTraitement = "Gradient de Roberts";
      CommentTraitement = "Détection non optimal des contours";
      dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
   }
}

afx_msg void CImageWnd:: OnPrewitt()
{   
   ArreteLeContinu();
   CDialDoubleSeuil dial(IDD_SEUIL2,&Seuil,&Seuil2);  
   if (dial.DoModal()==IDOK)
   {  
      SetCursor(LoadCursor(0,IDC_WAIT)); 
      Memorise(); 
      Prewitt(GetSImage(),Seuil,Seuil2); 
      SetCursor(LoadCursor(0,IDC_ARROW)); 
      Invalidate(FALSE);
      NomTraitement = "Gradient de Prewitt";
      CommentTraitement = "Détection non optimal des contours";
      dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
   }	
}

afx_msg void CImageWnd:: OnSobel()
{   
   ArreteLeContinu();
   CDialDoubleSeuil dial(IDD_SEUIL2,&Seuil,&Seuil2);  
   if (dial.DoModal()==IDOK)
   {
      Memorise();
      SetCursor(LoadCursor(0,IDC_WAIT));  
      Sobel(GetSImage(),Seuil,Seuil2); 
      SetCursor(LoadCursor(0,IDC_ARROW)); 
      Invalidate(FALSE);
      NomTraitement = "Gradient de Sobel";
      CommentTraitement = "Détection non optimal des contours";
      dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
   }
}


afx_msg void CImageWnd:: OnFFT()
{   
   ArreteLeContinu();
   SetCursor(LoadCursor(0,IDC_WAIT));  
   long l = TestPuissance2(fin_i-debut_i+1);
   long h = TestPuissance2(fin_j-debut_j+1); 
   fin_i = debut_i + puissance2(l)-1; 
   fin_j = debut_j + puissance2(h)-1;
   CTraitImageWnd * pImageTransModule = (CTraitImageWnd * ) new CTraitImageWnd((fin_i-debut_i+1),(fin_j-debut_j+1),NB);   
   CTraitImageWnd * pImageTransPhase =  (CTraitImageWnd * ) new CTraitImageWnd((fin_i-debut_i+1),(fin_j-debut_j+1),NB);   
   // ** 
   FFT2D(GetSImage(),pImageTransModule->GetSImage(),pImageTransPhase->GetSImage()); 
   // **
   SetCursor(LoadCursor(0,IDC_ARROW));    
   pImageTransModule	->	Create(Titre+" : FFT - MODULE",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   pImageTransPhase	->	Create(Titre+" : FFT - PHASE",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   Invalidate(FALSE);
   NomTraitement = "FFT";
   CommentTraitement = "Détection des variations de frequences ";
   dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
}

//******************************************************************************************//

void CImageWnd::OnRotationsDroite()
{
   ArreteLeContinu();  
   // pour bitmap largeur multiple de 4

   long Fin_j=fin_j;
   div_t res = div((int)(fin_j-debut_j+1),4);   
   if (res.rem !=0) Fin_j -= res.rem;

   CTraitImageWnd * pImageRot =  (CTraitImageWnd * ) new CTraitImageWnd(Fin_j-debut_j+1,fin_i-debut_i+1,NB);   
   SetCursor(LoadCursor(0,IDC_WAIT));
   RotationD(GetSImage() , pImageRot->GetSImage() ); 
   SetCursor(LoadCursor(0,IDC_ARROW));  
   pImageRot->Create(Titre+" : Rotation Droite",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   pImageRot->NomTraitement = NomTraitement = "Rotation a droite de 90°";
   pImageRot->CommentTraitement = CommentTraitement = " ";
   pImageRot->dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
}

void CImageWnd::OnRotationsGauche()
{
   ArreteLeContinu();   
   // pour bitmap largeur multiple de 4
   long Fin_j=fin_j;
   div_t res = div((int)(fin_j-debut_j+1),4);   
   if (res.rem !=0) Fin_j -= res.rem;

   CTraitImageWnd * pImageRot =  (CTraitImageWnd * ) new CTraitImageWnd(Fin_j-debut_j+1,fin_i-debut_i+1,NB);   
   SetCursor(LoadCursor(0,IDC_WAIT));
   RotationG(GetSImage() , pImageRot->GetSImage()); 
   SetCursor(LoadCursor(0,IDC_ARROW));  

   pImageRot->Create(Titre+" :  Rotation Gauche",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   NomTraitement = "Rotation a gauche de 90°";
   CommentTraitement = " ";
   dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
}

void CImageWnd::OnBSpline()
{
   // --- Filtres ---
   SetCursor(LoadCursor(0,IDC_WAIT)); 
   // --- Filtres ---
   SFiltre h,g;
   SImage SSource,SDest,SCont,SReconst;
   BOOL ok_reconst=FALSE,ok_cont=FALSE,ok_basse=FALSE,ok_extra=FALSE;   

   CTraitImageWnd * pImage = NULL;	
   CTraitImageWnd * rImage = NULL;
   CTraitImageWnd * cImage = NULL;

   // ---------------- Initialisation des filtres ------------------
   //-- filtre associé à la fonction d'échelle --
   h.Filtre=new float[64*sizeof(float)];     
   h.Fin=63;h.Centre=32;

   h.Filtre[0] = (float)0.000004;h.Filtre[1] = (float)-0.000001;h.Filtre[2] = (float)-0.000005;h.Filtre[3] = (float)0.000005;h.Filtre[4] = (float)0.000008;h.Filtre[5] = (float)-0.000011;h.Filtre[6] = (float)-0.000014;
   h.Filtre[7] = (float)0.000021;h.Filtre[8] = (float)0.000027;h.Filtre[9] = (float)-0.000042;h.Filtre[10] = (float)-0.000053;h.Filtre[11] = (float)0.000083;h.Filtre[12] = (float)0.000103;h.Filtre[13] = (float)-0.000164;
   h.Filtre[14] = (float)-0.000202;h.Filtre[15] = (float)0.000327;h.Filtre[16] = (float)0.000396;h.Filtre[17] = (float)-0.000656;h.Filtre[18] = (float)-0.000780;h.Filtre[19] = (float)0.001331;h.Filtre[20] = (float)0.001546;
   h.Filtre[21] = (float)-0.002745;h.Filtre[22] = (float)-0.003079;h.Filtre[23] = (float)0.005800;h.Filtre[24] = (float)0.006142;h.Filtre[25] = (float)-0.012715;h.Filtre[26] = (float)-0.012146;h.Filtre[27] = (float)0.029747;
   h.Filtre[28] = (float)0.022685;h.Filtre[29] = (float)-0.077808;	h.Filtre[30] = (float)-0.035498;h.Filtre[31] = (float)0.306830;h.Filtre[32] = (float)0.541736;

   h.Filtre[63] = (float)0.000001;h.Filtre[62] = (float)-0.000005;h.Filtre[61] = (float)0.000005;h.Filtre[60] = (float)0.000008;h.Filtre[59] = (float)-0.000011; h.Filtre[58] = (float)-0.000014;
   h.Filtre[57] = (float)0.000021;h.Filtre[56] = (float)0.000027;
   h.Filtre[55] = (float)-0.000042;h.Filtre[54] = (float)-0.000053;h.Filtre[53] = (float)0.000083;
   h.Filtre[52] = (float)0.000103;h.Filtre[51] = (float)-0.000164;h.Filtre[50] = (float)-0.000202;
   h.Filtre[49] = (float)0.000327; h.Filtre[48] = (float)0.000396;h.Filtre[47] = (float)-0.000656;
   h.Filtre[46] = (float)-0.000780;h.Filtre[45] = (float)0.001331;h.Filtre[44] = (float)0.001546;
   h.Filtre[43] = (float)-0.002745;h.Filtre[42] = (float)-0.003079;
   h.Filtre[41] = (float)0.005799;h.Filtre[40] = (float)0.006142;h.Filtre[39] = (float)-0.012715;
   h.Filtre[38] = (float)-0.012146;h.Filtre[37] = (float)0.029747;h.Filtre[36] = (float)0.022685;h.Filtre[35] = (float)-0.077808;
   h.Filtre[34] = (float)-0.035498;h.Filtre[33] = (float)0.306830;


   //-- filtre associé à l'ondelette --	
   g.Filtre=new float[64*sizeof(float)];	
   g.Fin=63;g.Centre=32;
   for (int i=0;i<=h.Fin;i++) g.Filtre[i] = (float)pow((double)-1,(double)i)*h.Filtre[h.Fin-i];

   // --------- Controle du menu -----------
   // affectation par defaut
   SSource=GetSImage();
   SCont=GetSImage();
   SReconst=GetSImage();	
   pImage =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);  
   pImage->Echelle = Echelle*2;
   SDest=pImage->GetSImage();  
   if (menu.GetMenuState(ID_MULTI_ECHELLE,MF_BYCOMMAND) ==  MF_CHECKED)	
   { 
      Extraction(h,g,SSource,SDest,"B-Splines cubiques");	
      SetCursor(LoadCursor(0,IDC_WAIT));
   }

   if (menu.GetMenuState( ID_RECONSTRUCT,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      rImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut); 
      rImage->Echelle = Echelle; 
      ok_reconst=TRUE;    
      SReconst=rImage->GetSImage();
   }

   if (menu.GetMenuState( ID_CONT,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);  
      ok_cont=TRUE;      
      SCont=cImage->GetSImage();
   }

   if (menu.GetMenuState( ID_BASSE,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);  
      ok_basse=TRUE;      
      SCont=cImage->GetSImage();
   }           

   if (menu.GetMenuState( ID_EXTRA_PAUL,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);  
      ok_extra=TRUE;      
      SCont=cImage->GetSImage();
   }                        

   TOD(h,g,SSource,SDest,SReconst,SCont,ok_reconst,ok_cont,ok_basse,ok_extra,Echelle);	

   pImage->Create(Titre+" Analyse Multi-résolutions : B-Splines cubiques ",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 

   if (menu.GetMenuState( ID_RECONSTRUCT,MF_BYCOMMAND) ==  MF_CHECKED)	
      rImage->Create(Titre+" Reconstruction : B-Splines cubiques",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 

   if (menu.GetMenuState( ID_CONT,MF_BYCOMMAND) ==  MF_CHECKED)	
      cImage->Create(Titre+" Extraction des détails : B-Splines cubiques",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 

   if (menu.GetMenuState( ID_BASSE,MF_BYCOMMAND) ==  MF_CHECKED)	
      cImage->Create(Titre+" Extraction des Vij : B-Splines cubiques",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 

   if (menu.GetMenuState( ID_EXTRA_PAUL,MF_BYCOMMAND) ==  MF_CHECKED)	
      cImage->Create(Titre+" Extrapolation : B-Splines cubiques",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 

   SetCursor(LoadCursor(0,IDC_ARROW));

   delete g.Filtre;
   delete h.Filtre;	
   Invalidate(FALSE);

}

void CImageWnd::OnReconstruct()
{
   if (menu.GetMenuState( ID_RECONSTRUCT,MF_BYCOMMAND) ==  MF_UNCHECKED)
      menu.CheckMenuItem(ID_RECONSTRUCT,MF_CHECKED);  
   else  menu.CheckMenuItem(ID_RECONSTRUCT,MF_UNCHECKED);  
}


void CImageWnd::OnCont()
{ 
   if (menu.GetMenuState( ID_CONT,MF_BYCOMMAND) ==  MF_UNCHECKED)
      menu.CheckMenuItem(ID_CONT,MF_CHECKED);  
   else  menu.CheckMenuItem(ID_CONT,MF_UNCHECKED);  

}

void CImageWnd::OnHaar()
{
   SetCursor(LoadCursor(0,IDC_WAIT)); 

   // TODO: Add your command handler code here
   // --- Filtres ---
   SFiltre h,g;
   SImage SCont,SReconst,SSource,SDest;
   BOOL ok_reconst=FALSE,ok_cont=FALSE,ok_basse=FALSE,ok_extra=FALSE;
   CTraitImageWnd * pImage = NULL;	
   CTraitImageWnd * rImage = NULL;
   CTraitImageWnd * cImage = NULL;

   // ---------------- Initialisation des filtres ------------------

   //-- filtre associé à la fonction d'échelle --

   h.Filtre=new float[2];     
   h.Fin=1;h.Centre=1;
   h.Filtre[0] = (float)(1.0/2.0);h.Filtre[1] = (float)(1.0/2.0);    	   
   //-- filtre associé à l'ondelette --	

   g.Filtre=new float[2];	
   g.Fin=1;g.Centre=1;
   g.Filtre[0] = (float)(-1.0/2.0);g.Filtre[1] = (float)(1.0/2.0);	                
   // --------- Controle du menu -----------
   //
   // affectation par defaut
   SSource=GetSImage();	// valeur exacte
   SCont=GetSImage();
   SReconst=GetSImage();    
   pImage =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);
   pImage->Echelle = Echelle*2;

   SDest=pImage->GetSImage();  

   if (menu.GetMenuState( ID_MULTI_ECHELLE,MF_BYCOMMAND) ==  MF_CHECKED)	
   {                                           
      Extraction(h,g,SSource,SDest,"HAAR");	
      SetCursor(LoadCursor(0,IDC_WAIT)); 

   }


   if (menu.GetMenuState( ID_RECONSTRUCT,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      rImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);  
      rImage->Echelle = Echelle;
      ok_reconst=TRUE;    
      SReconst=rImage->GetSImage();
   }

   if (menu.GetMenuState( ID_CONT,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);  
      cImage->Echelle = 1;
      ok_cont=TRUE;      
      SCont=cImage->GetSImage();
   }                        
   if (menu.GetMenuState( ID_BASSE,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);  
      cImage->Echelle = 1;
      ok_basse=TRUE;      
      SCont=cImage->GetSImage();
   }                        
   if (menu.GetMenuState( ID_EXTRA_PAUL,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);  
      cImage->Echelle = 1;
      ok_extra=TRUE;      
      SCont=cImage->GetSImage();
   }
   TOD(h,g,SSource,SDest,SReconst,SCont,ok_reconst,ok_cont,ok_basse,ok_extra,Echelle);	

   pImage->Create(Titre+" Analyse Multi-résolutions : HAAR",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 

   if (menu.GetMenuState( ID_RECONSTRUCT,MF_BYCOMMAND) ==  MF_CHECKED)	
      rImage->Create(Titre+" Reconstruction : HAAR",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 

   if (menu.GetMenuState( ID_CONT,MF_BYCOMMAND) ==  MF_CHECKED)	
      cImage->Create(Titre+" Extraction des détails : HAAR",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 

   if (menu.GetMenuState( ID_BASSE,MF_BYCOMMAND) ==  MF_CHECKED)	
      cImage->Create(Titre+" Extraction des Vij : HAAR",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 

   if (menu.GetMenuState( ID_EXTRA_PAUL,MF_BYCOMMAND) ==  MF_CHECKED)	
      cImage->Create(Titre+" Extrapolation : HAAR",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 

   SetCursor(LoadCursor(0,IDC_ARROW));


   Invalidate(FALSE);
   delete h.Filtre;
   delete g.Filtre;
}


//************************************** FIN DES TRAITEMENTS ************************
//************************************** FIN DES TRAITEMENTS ************************
//************************************** FIN DES TRAITEMENTS ************************

void CImageWnd::ArreteLeContinu()
{
} 

void CImageWnd::MettreLeContinu()
{
} 

// **********************************************************************************
// *************** GESTION DE LA SEQUENCE 
// **********************************************************************************

afx_msg  void CImageWnd::OnSauveSequence()  
{
   ArreteLeContinu();   

   int  num1=1,num2; 
   num2=NbMaxImageSeq;
   CDialSeqlimite dial(&num1,&num2); 
   if (dial.DoModal()==IDOK)
   {
      CFile Fichier; // variable fichier BITMAP	 
      int l; 
      CString NomFichier,s;
      char NumExt[20];

      CFileDialog  FileData(FALSE,"*.bmp",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"BITMAP (*.bmp)|*.bmp||",this);
      if (FileData.DoModal()==IDOK) 
      {
         SetCursor(LoadCursor(0,IDC_WAIT));  
         MyInitInfo->bmiHeader.biWidth =  (fin_i_seq-debut_i_seq+1);
         MyInitInfo->bmiHeader.biHeight = (fin_j_seq-debut_j_seq+1);
         l = FileData.GetPathName().GetLength();
         s = FileData.GetPathName().Left(l-4); 
         unsigned char  * P;
         for (int i=num1;i<=num2;i++)
         {    
            sprintf(NumExt,"%.3d.bmp",i); 
            NomFichier = s + NumExt;  
            Fichier.Open(NomFichier,CFile::modeCreate | CFile::modeWrite) ;     
            Fichier.Write(TheTete,(UINT) sizeof(BITMAPFILEHEADER)); 
            Fichier.Write(MyInitInfo,sizeof(BITMAPINFOHEADER)+256*sizeof(RGBQUAD));  
            P = (unsigned char  *) Bpile[i-1];
            Fichier.Write(P,TailleSeq);   
            Fichier.Close();
         } 
         MyInitInfo->bmiHeader.biWidth =  Larg;
         MyInitInfo->bmiHeader.biHeight = Haut;  
         SetCursor(LoadCursor(0,IDC_ARROW));
      }
   }  
}

void CImageWnd::OnUpdateSauveSequence(CCmdUI* pCmdUI) 
{
   pCmdUI->Enable(SEQUENCE_CREE);
}


void  CImageWnd::AlloueSequence(int * Nb)
{     
   Bpile = new unsigned char  * [*Nb];  
   SEQUENCE_CREE = TRUE;
   debut_i_seq=debut_i;
   fin_i_seq=fin_i;
   debut_j_seq=debut_j;
   fin_j_seq=fin_j;
   TailleSeq = (fin_i_seq-debut_i_seq+1)*(fin_j_seq-debut_j_seq+1); 
   NbImage = *Nb;  
   for (int i=0;i<NbImage;i++)
      Bpile[i] = (unsigned char   *) malloc( TailleSeq); 
   //*Nb = NbImage;  
   NbMaxImageSeq = NbImage; 
}               

void CImageWnd::AfficheSequence(int debut,int fin)
{   
   long c,p;
   char chiffre[30];
   for (int l=debut;l<fin;l++)
   {
      c= 0;
      Buf = Bpile[l];  
      for (int j=(int) debut_j_seq;j<=(int)fin_j_seq;j++)
      {
         p = (long)j*Larg+debut_i_seq; 
         for (int i=(int)debut_i_seq;i<=(int)fin_i_seq;i++)
            ((unsigned char  *)Image)[p++] =  Buf[c++];
      }
      Affichage();
      if ((CoupeSelectionnee) && (pCoupeVariable!=NULL)) MettreAJourLaCoupeActive();
      SetWindowText(Titre+": "+_itoa(l+1,chiffre,10));

   }

   Invalidate(FALSE);
}   


afx_msg  void CImageWnd::OnDetruireSequence()      
{
   if (SEQUENCE_CREE)
   {
      for (int i=0;i<NbMaxImageSeq;i++)
      {                 
         free(Bpile[i]);
      }
      SEQUENCE_CREE = FALSE;   
      SEQ_EXISTE=FALSE;
      delete Bpile;  
      if(IsWindow(DialLecSeq->GetSafeHwnd())) DialLecSeq->ShowWindow(SW_HIDE);
   }
}  	

void CImageWnd::OnUpdateDetruireSequence(CCmdUI* pCmdUI) 
{
   pCmdUI->Enable(SEQUENCE_CREE);
}


afx_msg  void CImageWnd::OnLireSequence()
{
   ArreteLeContinu();  
   if(!IsWindow(DialLecSeq->GetSafeHwnd())) 
      DialLecSeq->Create(NbMaxImageSeq,Titre);
   else 
   {	
      DialLecSeq->Redefinition(NbMaxImageSeq);
      DialLecSeq->ShowWindow(SW_SHOW);
   }

}

void CImageWnd::OnUpdateLireSequence(CCmdUI* pCmdUI) 
{
   pCmdUI->Enable(SEQUENCE_CREE);
}

// *************************** FIN DE GESTION DE LA SEQUENCE **********************


LONG CImageWnd::OnSetCoupeVariable(UINT uiParam,LONG lparam )
{ 
   CProfilWnd * pCoupeVariableAvant = pCoupeVariable;                       
   pCoupeVariable = (CProfilWnd *) lparam; 
   if(lparam==NULL) 
      NbProfil--;
   if (pCoupeVariableAvant!=pCoupeVariable)
   { 
      if (pCoupeVariable == NULL && !NbProfil) CoupeSelectionnee = FALSE;
      else if (pCoupeVariable != NULL)
      {   pCoupeVariable->GetInfo(&Lg_i1,&Lg_j1,&Lg_i2,&Lg_j2,&NbPointCoupe);
      CoupeSelectionnee = TRUE; 
      }
      else if (pCoupeVariableAvant != NULL && NbProfil) {
         pCoupeVariable = pCoupeVariableAvant;
         pCoupeVariable->GetInfo(&Lg_i1,&Lg_j1,&Lg_i2,&Lg_j2,&NbPointCoupe);
         CoupeSelectionnee = TRUE;
      }
   }	
   Invalidate(TRUE);  
   return 0;
}    


void CImageWnd::OnSeparationRGB()
{
   rgb RGB;
   long l = 0;
   CTraitImageWnd * r = new CTraitImageWnd((fin_i-debut_i+1),(fin_j-debut_j+1),NB);
   CTraitImageWnd * v = new CTraitImageWnd((fin_i-debut_i+1),(fin_j-debut_j+1),NB);
   CTraitImageWnd * b = new CTraitImageWnd((fin_i-debut_i+1),(fin_j-debut_j+1),NB);          
   unsigned char   * ImR = (unsigned char   *) r->Image;
   unsigned char   * ImV = (unsigned char   *) v->Image;
   unsigned char   * ImB = (unsigned char   *) b->Image;
   for (long j= debut_j;j<=fin_j;j++) 
   {                              
      for (long i=debut_i;i<=fin_i;i++) 
      {       
         RGB = ImageRGB[j*Larg+i];
         ((unsigned char   * )ImR)[l] = RGB.r; 
         ((unsigned char   * )ImV)[l] = RGB.g; 
         ((unsigned char   * )ImB)[l] = RGB.b; 
         l++;
      } 
   }	
   r->Create(Titre + " : rouge",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   v->Create(Titre + " : vert",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   b->Create(Titre + " : bleu",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
}

void CImageWnd::OnLumchro()
{
   rgb RGB;                                                        
   long l = 0;
   CTraitImageWnd * Y = new CTraitImageWnd((fin_i-debut_i+1),(fin_j-debut_j+1),NB);
   CTraitImageWnd * C0 = new CTraitImageWnd((fin_i-debut_i+1),(fin_j-debut_j+1),NB);
   CTraitImageWnd * C1 = new CTraitImageWnd((fin_i-debut_i+1),(fin_j-debut_j+1),NB);
   unsigned char   * ImY = (unsigned char   *) Y->Image;   
   unsigned char   * ImC0 = (unsigned char   *) C0->Image;   
   unsigned char   * ImC1 = (unsigned char   *) C1->Image;   
   int  * Chromi0  = (int  * ) malloc( (fin_i-debut_i+1)*(fin_j-debut_j+1)*sizeof(int));
   int  * Chromi1  = (int  *) malloc( (fin_i-debut_i+1)*(fin_j-debut_j+1)*sizeof(int));
   int lum;
   long i,j;  
   int MaxRouge=-255,MaxBleu=-255;
   int MinRouge=255,MinBleu=255;
   for (j= debut_j;j<=fin_j;j++) 
   {                              
      for (i=debut_i;i<=fin_i;i++) 
      {                              
         RGB = ImageRGB[j*Larg+i]; 
         lum =  (int)(0.3*RGB.r + 0.59*RGB.g + 0.11*RGB.b);
         ((unsigned char   * )ImY)[l] = (unsigned char) lum; 
         Chromi0[l] = (int)RGB.r-lum;
         Chromi1[l] = (int)RGB.b-lum;
         if (Chromi0[l] 	> MaxRouge) MaxRouge =Chromi0[l];
         if (Chromi1[l]  > MaxBleu) MaxBleu = Chromi1[l] ;
         if (Chromi0[l]  < MinRouge) MinRouge = Chromi0[l] ;
         if (Chromi1[l]  < MinBleu) MinBleu = Chromi1[l] ;
         l++;
      }
   } 
   l = 0;   
   // normalisation des chrominances 
   float normR,normB,norm; 
   if (MaxRouge==MinRouge) normR= 1; else normR = (MaxRouge-MinRouge)/255;
   if (MaxBleu==MinBleu) normB= 1; else normB = (MaxBleu-MinBleu)/255;      
   if (normR>normB) norm = normR; else norm = normB;   
   for (j= debut_j;j<=fin_j;j++) 
   {                              
      for (i=debut_i;i<=fin_i;i++) 
      {                              
         RGB = ImageRGB[j*Larg+i]; 
         ((unsigned char   * )ImC0)[l] = (unsigned char) (Chromi0[l]-MinRouge)/norm;
         ((unsigned char   * )ImC1)[l] = (unsigned char) (Chromi1[l] -MinBleu)/norm;
         l++;
      }
   }

   free(Chromi0);	
   free(Chromi1);	
   Y->Create(Titre + " : Luminance Y R(30%) G(59%) B(11%)",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   C0->Create(Titre + " : Chrominance R-Y",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   C1->Create(Titre + " : Chrominance B-Y",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
}


void CImageWnd::OnZoom100()
{
   ZOOM = 1; 
   menu.CheckMenuItem(ID_ZOOM_100,MF_CHECKED);
   menu.CheckMenuItem(ID_ZOOM_50,MF_UNCHECKED);
   menu.CheckMenuItem(ID_ZOOM_25,MF_UNCHECKED);
   Invalidate(TRUE);
   NomTraitement = "Image initiale";
   CommentTraitement = " Pas de zoom";
   dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
}

void CImageWnd::OnZoom25()
{                                                 
   ZOOM = 4;
   menu.CheckMenuItem(ID_ZOOM_100,MF_UNCHECKED);
   menu.CheckMenuItem(ID_ZOOM_50,MF_UNCHECKED);
   menu.CheckMenuItem(ID_ZOOM_25,MF_CHECKED);
   Invalidate(TRUE);
   NomTraitement = "Zoom 1/4";
   CommentTraitement = "Reduction de la taille de l'image";
   dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
}

void CImageWnd::OnZoom50()
{
   ZOOM = 2;
   menu.CheckMenuItem(ID_ZOOM_100,MF_UNCHECKED);
   menu.CheckMenuItem(ID_ZOOM_50,MF_CHECKED);
   menu.CheckMenuItem(ID_ZOOM_25,MF_UNCHECKED);
   Invalidate(TRUE);
   NomTraitement = "Zoom 1/2";
   CommentTraitement = "Reduction de la taille de l'image";
   dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
}


void CImageWnd::OnSousEchanCol()
{    
   CTraitImageWnd * PImageEchant = new CTraitImageWnd(Larg,Haut/2,NB);
   SousEchantillonnageColonne(GetSImage(),PImageEchant->GetSImage());
   PImageEchant->Create(Titre + " Sous-échantillonnage Colonne",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
}

void CImageWnd::OnSousEchanLigne()
{
   CTraitImageWnd * PImageEchant = new CTraitImageWnd(Larg/2,Haut,NB);
   SousEchantillonnageLigne(GetSImage(),PImageEchant->GetSImage());
   PImageEchant->Create(Titre + " Sous-échantillonnage Ligne",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
}

void CImageWnd::OnAffecterEchelle()
{
   CDialEchelle dial(Echelle);
   if (dial.DoModal()==IDOK) Echelle = dial.m_Echelle;
}

int CImageWnd::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
   if (CMDIChildWnd::OnCreate(lpCreateStruct) == -1)
      return -1;

   // Create an array for status bar indicators
   UINT pIndicators[2] = {ID_SEPARATOR, ID_SEPARATOR};
   if (!m_wndStatusBarImage.Create(this, 
      WS_CHILD | WS_VISIBLE | CBRS_TOP , ID_IMAGE_STATUS_BAR) ||
      !m_wndStatusBarImage.SetIndicators(pIndicators,
      2))
      TRACE0("Failed to create Status Bar\n");
   Invalidate(TRUE);

   m_wndStatusBarImage.GetStatusBarCtrl().ModifyStyle(SBARS_SIZEGRIP,CCS_NORESIZE,SWP_SHOWWINDOW);

   UINT nID, nStyle;
   int nWidth;
   m_wndStatusBarImage.GetPaneInfo(0, nID, nStyle, nWidth);
   m_wndStatusBarImage.SetPaneInfo(0, nID, nStyle, LFene/3);
   return 0;
}

void CImageWnd::OnMDIActivate(BOOL bActivate, CWnd* pActivateWnd, CWnd* pDeactivateWnd) 
{
   CMDIChildWnd::OnMDIActivate(bActivate, pActivateWnd, pDeactivateWnd);

   AfficheDimensions();
   (((CMainFrame*) GetMDIFrame())->bar.GetToolBarCtrl()).SetOwner((CWnd*)this );
   ((CMainFrame*) GetMDIFrame())->bar.OnUpdateCmdUI( (CFrameWnd*) GetMDIFrame(), FALSE);

}

void CImageWnd::OnClose() 
{
   if (IMAGE_MODIFIEE)
   {
      int test = MessageBox("Voulez-vous enregistrez les modifications avant de fermer le fichier :\n" + FileName,
         "Sauvegarde du fichier",MB_DEFBUTTON2|MB_YESNOCANCEL|MB_ICONQUESTION);
      switch (test)
      {
      case IDYES:
         OnSauve();
         CMDIChildWnd::OnClose();
         break;
      case IDNO:
         CMDIChildWnd::OnClose();
         break;
      case IDCANCEL:
         break;
      }
   }
   else
      CMDIChildWnd::OnClose();
}
int CImageWnd::CloseFen()
{
   int test=IDCANCEL;
   if (IMAGE_MODIFIEE)
   {
      test = MessageBox("Voulez-vous enregistrez les modifications avant de fermer le fichier :\n" + FileName,
         "Sauvegarde du fichier",MB_DEFBUTTON2|MB_YESNOCANCEL|MB_ICONQUESTION);
      switch (test)
      {
      case IDYES:
         OnSauve();
         IMAGE_MODIFIEE=FALSE;
         OnClose();
         break;
      case IDNO:
         IMAGE_MODIFIEE=FALSE;
         OnClose();
         break;
      case IDCANCEL:
         break;
      }
   }
   else OnClose();
   return test;
}

void CImageWnd::OnDestroy() 
{
   CMDIChildWnd::OnDestroy();

   char titre[100];
   CStatusBar *pStatus = (CStatusBar *)
      AfxGetApp()->m_pMainWnd->GetDescendantWindow(AFX_IDW_STATUS_BAR);
   if (pStatus)
   {
      wsprintf(titre," ");
      pStatus->SetPaneText(1, titre);	// Barre d'etat de MainFrame
   }

   if (GetMDIFrame()->MDIGetActive() == NULL)
   {
      CMenu toto;
      toto.LoadMenu(IDR_MAINFRAME);
      GetMDIFrame()->MDISetMenu(&toto,NULL);
      AfxGetMainWnd()->DrawMenuBar();
   }

   // Suppression des coupes associees à l'image en partant du haut de la pile de coupes
   while(NbProfil>0)
   {
      GetMDIFrame()->SendMessageToDescendants( uiMsg3,0,(LPARAM)this,TRUE,FALSE);
   }

   (((CMainFrame*) GetMDIFrame())->bar.GetToolBarCtrl()).SetOwner( (CWnd*)GetMDIFrame() );
   ((CMainFrame*) GetMDIFrame())->bar.OnUpdateCmdUI( (CFrameWnd*) GetMDIFrame(), TRUE);
}

/////////////////////////////////////////////////////////////////////
///			GESTION AFFICHAGE INFO			/////////////////////////
/////////////////////////////////////////////////////////////////////
void CImageWnd::OnInfoImage() 
{
   if (!INFO_IMAGE)
   {
      INFO_IMAGE = TRUE;
      if (dialInfo->GetSafeHwnd()==0)
         dialInfo->Create(CommentTraitement,NomTraitement);
      else
      {
         dialInfo->ShowWindow(SW_SHOW);
         dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
      }
   }
   else 
   {
      INFO_IMAGE = FALSE;
      dialInfo->ShowWindow(SW_HIDE);
   }
}

CString CImageWnd::ConcateneInfo(CString C1,CString C2,CString C3,CString C4,CString C5,CString C6)
{
   CString ChaineRetour = _T("");
   ChaineRetour += C1;
   ChaineRetour += char(13); ChaineRetour += char(10);
   ChaineRetour += C2;
   ChaineRetour += char(13); ChaineRetour += char(10);
   ChaineRetour += C3;
   ChaineRetour += char(13); ChaineRetour += char(10);
   ChaineRetour += C4;
   ChaineRetour += char(13); ChaineRetour += char(10);
   ChaineRetour += C5;
   ChaineRetour += char(13); ChaineRetour += char(10);
   ChaineRetour += C6;
   return ChaineRetour;
}
void CImageWnd::OnUpdateInfoImage(CCmdUI* pCmdUI) 
{
   pCmdUI->SetRadio(INFO_IMAGE);
}

/////////////////////////////////////////////////////////////////////
///			FIN GESTION AFFICHAGE INFO			/////////////////////
/////////////////////////////////////////////////////////////////////

void CImageWnd::OnSauveExcel() 
{
   CFile Fichier; // variable fichier BITMAP
   char temp[10];
   int nb;
   CFileDialog FileData(FALSE,"*.xls",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"Fichier Excel Text (*.xls)|*.xls||",this);
   if (FileData.DoModal()==IDOK)
      if (Fichier.Open(FileData.GetPathName(),CFile::modeCreate | CFile::modeWrite))
      {   
         SetCursor(LoadCursor(0,IDC_WAIT));  

         // sauvegarde de l'entête si BMP
         if (FileData.GetFileExt() == "xls")
         {
            for (long j=debut_j;j<fin_j;j++)
            {
               for (long i=debut_i; i<fin_i;i++)
               {
                  nb = sprintf(temp,"%d ",Image[i+Larg*(j)]);
                  Fichier.Write(temp,nb); 
               }
               Fichier.Write("\n",1); 
            }
         }	
         Fichier.Close();    
         SetCursor(LoadCursor(0,IDC_ARROW));
      }  
}

//*****************************************************
////////////  REDUCTEURS de BRUITS  type MOYENNEUR  ///
//*****************************************************
void CImageWnd::OnMoyenneurxy121() 
{
   float tab[3][3];
   SetCursor(LoadCursor(0,IDC_WAIT));
   tab[0][1]=tab[1][0]=tab[1][2]=tab[2][1]=2.0F;
   tab[1][1]=4.0F;
   tab[0][0]=tab[0][2]=tab[2][0]=tab[2][2]=1.0F;
   CalculLissage(tab);
   SetCursor(LoadCursor(0,IDC_ARROW));
   Invalidate(FALSE);
   NomTraitement = _T("Moyenneur XY 121");
   s1 = _T("Atténuation du bruit");
   s2 = _T("Utilisation d'un filtre passe bas type moyenneur.");
   s3 = _T("1 2 1");
   s4 = _T("2 4 2");
   s5 = _T("1 2 1");
   CommentTraitement = ConcateneInfo(s1,s2,s3,s4,s5);
   dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
}

void CImageWnd::OnMoyenneurxy111() 
{
   float tab[3][3];
   SetCursor(LoadCursor(0,IDC_WAIT));
   tab[0][1]=tab[1][0]=tab[1][2]=tab[2][1]=1.0F;
   tab[1][1]=1.0F;
   tab[0][0]=tab[0][2]=tab[2][0]=tab[2][2]=1.0F;
   CalculLissage(tab);
   SetCursor(LoadCursor(0,IDC_ARROW));
   Invalidate(FALSE);
   NomTraitement = _T("Moyenneur XY 111");
   s1 = _T("Atténuation du bruit");
   s2 = _T("Utilisation d'un filtre passe bas type moyenneur.");
   s3 = _T("1 1 1");
   s4 = _T("1 1 1");
   s5 = _T("1 1 1");
   CommentTraitement = ConcateneInfo(s1,s2,s3,s4,s5);
   dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
}

void CImageWnd::OnMoyenneurx111() 
{
   float tab[3][3];
   SetCursor(LoadCursor(0,IDC_WAIT));
   tab[1][0]=tab[1][1]=tab[1][2]=1.0F;
   tab[0][0]=tab[0][1]=tab[0][2]=tab[2][0]=tab[2][2]=tab[2][1]=0.0F;
   CalculLissage(tab);
   SetCursor(LoadCursor(0,IDC_ARROW));
   NomTraitement = _T("Moyenneur X 111");
   s1 = _T("Atténuation du bruit");
   s2 = _T("Utilisation d'un filtre passe bas type moyenneur.");
   s3 = _T("0 0 0");
   s4 = _T("1 1 1");
   s5 = _T("0 0 0");
   CommentTraitement = ConcateneInfo(s1,s2,s3,s4,s5);
   dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
   Invalidate(FALSE);
}

void CImageWnd::OnMoyenneury111() 
{
   float tab[3][3];
   SetCursor(LoadCursor(0,IDC_WAIT));
   tab[0][1]=tab[2][1]=1.0F;
   tab[1][1]= 1.0F;
   tab[0][0]=tab[1][0]=tab[0][2]=tab[2][0]=tab[2][2]=tab[1][2]=0.0F;
   CalculLissage(tab);
   SetCursor(LoadCursor(0,IDC_ARROW));
   NomTraitement = _T("Moyenneur Y 111");
   s1 = _T("Atténuation du bruit");
   s2 = _T("Utilisation d'un filtre passe bas type moyenneur.");
   s3 = _T("0 1 0");
   s4 = _T("0 1 0");
   s5 = _T("0 1 0");
   CommentTraitement = ConcateneInfo(s1,s2,s3,s4,s5);
   dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
   Invalidate(FALSE);
}

///////////////////////////////////////////////////////////////////////////
void CImageWnd :: CalculLissage(float CoefL[3][3])
{
   HGLOBAL ic;
   float  * ImageBuf;
   DWORD TailleBuf = Larg * Haut * sizeof(float);
   ic = GlobalAlloc(GHND,TailleBuf);            
   ImageBuf = (float  *) GlobalLock(ic);

   float TotalCoef = CoefL[0][0]+CoefL[0][1]+CoefL[0][2]+CoefL[1][0]+CoefL[1][1]+CoefL[1][2]+
      CoefL[2][0]+CoefL[2][1]+CoefL[2][2];                      
   long cpte=0,l, i, j;
   for (j = debut_j+1; j < fin_j; j++)
   {
      l = debut_i + 1 + Larg * j;
      for (i = debut_i+1; i < fin_i; i++)
      {
         ImageBuf[l] = CoefL[1][0]*Image[l-1] + CoefL[1][1]*Image[l] + CoefL[1][2]*Image[l+1] + CoefL[0][0]*Image[l+Larg-1];
         ImageBuf[l] += CoefL[0][1]*Image[l+Larg]+CoefL[0][2]*Image[l-Larg+1]+CoefL[2][0]*Image[l-Larg-1]+CoefL[2][1]*Image[l-Larg]+CoefL[2][2]*Image[l-Larg+1];
         ImageBuf[l] =(1/TotalCoef)*ImageBuf[l];
         l++;
      }
   }
   for (j = debut_j+1; j < fin_j; j++)
   {
      l = debut_i + 1 + Larg * j;
      for (i = debut_i+1; i < fin_i; i++)
         Image[l++] = (unsigned char)ImageBuf[l];
   } 

   GlobalUnlock(ic);
   free(ic);
}

void CImageWnd::OnZoomX() 
{
   CDialFactZoom dial(this);
   if (dial.DoModal() == IDOK)
   {
      CTraitImageWnd * pImageWnd =  new CTraitImageWnd((fin_i-debut_i+1)*m_FacteurZoom,(fin_j-debut_j+1)*m_FacteurZoom,BMP);
      for(long j=0;j<fin_j-debut_j+1;j++)
         for (long i=0;i<fin_i-debut_i+1;i++)
         {
            for(long k=0;k<m_FacteurZoom;k++)
            {
               for(long h=0;h<m_FacteurZoom;h++)
                  pImageWnd->Image[(m_FacteurZoom*i+k)+pImageWnd->Larg*(m_FacteurZoom*j+h)] = Image[(i+debut_i)+Larg*(j+debut_j)];
            }
         }

         pImageWnd->Create(Titre,WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());
   }
}

void CImageWnd::OnEtiquettage()
{
   ArreteLeContinu();

   BOOL flag;                        
   CDialSeuil dial (&Seuil,&flag); 
   if (dial.DoModal()==IDOK)
   {              

      Memorise(); 
      CTraitImageWnd * pImageEtiquette =  (CTraitImageWnd * ) new CTraitImageWnd((fin_i-debut_i+1),(fin_j-debut_j+1),NB);  
      CTraitImageWnd * pImageCentreGravite =  (CTraitImageWnd * ) new CTraitImageWnd((fin_i-debut_i+1),(fin_j-debut_j+1),NB);  
      SetCursor(LoadCursor(0,IDC_WAIT)); 
      long nbetiqu = Etiquett(GetSImage(), pImageEtiquette->GetSImage(),Seuil); 
      CentresGravites(pImageEtiquette->GetSImage(),pImageCentreGravite->GetSImage(), nbetiqu);
      SetCursor(LoadCursor(0,IDC_ARROW));  

      Invalidate(FALSE);
      UpdateWindow();  

      pImageEtiquette->Create(Titre+" : étiquettée", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
      pImageCentreGravite->Create(Titre+" : centre de gravité",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
      char chaine[30];
      sprintf(chaine,"Nombre d'étiquettes : %d",nbetiqu);
      MessageBox(chaine,"Résultat");
   }	
}

void CImageWnd::OnNormalisation()
{
   Normalisation(debut_i,debut_j,fin_i,fin_j);
   Invalidate(FALSE);
}
void CImageWnd::Normalisation(long i1,long j1,long i2,long j2)
{   

   float maxi,mini,temp;	 	
   float norm; 
   long i,j,l;

   //-- recherche du maximun et du minimun dans l'image float ---

   if (TypeImage==FLOAT32) { mini=maxi=ImageFloat[j1*Larg+i1];}
   if (TypeImage==MEM || TypeImage==BMP){ mini=maxi=Image[j1*Larg+i1];}

   for ( j = j1 ; j<=j2;j++) 	
   {                               
      l = j*Larg+i1;
      for ( i = i1 ; i<=i2;i++) 
      {
         if (TypeImage==FLOAT32) temp=ImageFloat[l++];
         if (TypeImage==MEM || TypeImage==BMP) temp=(float) Image[l++];
         if (temp>maxi) maxi=temp;
         if (temp<mini) mini=temp;
      }
   }
   // -- Création de la nouvelle image	normalisée                 
   if (maxi==mini) norm=1.F;
   else norm = (maxi-mini)/255;	
   for ( j = j1 ; j<=j2;j++) 	
   {                               
      l = j*Larg+i1;
      for ( i = i1 ; i<=i2;i++) 
      {
         if (TypeImage==FLOAT32) Image[l]= (unsigned char) ( (ImageFloat[l++]-mini)/norm ); 
         if (TypeImage==MEM || TypeImage==BMP) Image[l]= (unsigned char) ( (Image[l++]-mini)/norm );   
      }
   }				

}


void CImageWnd::OnHaarFloat()
{
   SetCursor(LoadCursor(0,IDC_WAIT)); 
   // --- Filtres ---
   SFiltre h,g;
   SImageFloat SCont,SReconst,SSource,SDest; // Type inconnu deux structures possibles
   long i;
   BOOL ok_reconst=FALSE,ok_cont=FALSE,ok_basse=FALSE;		 
   CTraitImageWnd * pImage = NULL;	
   CTraitImageWnd * rImage = NULL;
   CTraitImageWnd * cImage = NULL;

   // ---------------- Initialisation des filtres ------------------

   //-- filtre associé à la fonction d'échelle --

   h.Filtre=new float[2];     
   h.Fin=1;h.Centre=1;
   h.Filtre[0]=(1.F/2.F);h.Filtre[1]=(1.F/2.F);    	   
   //-- filtre associé à l'ondelette --	

   g.Filtre=new float[2];	
   g.Fin=1;g.Centre=1;
   g.Filtre[0]=(-1.F/2.F);g.Filtre[1]=(1.F/2.F);	                
   // --------- Controle du menu -----------
   //
   // affectation par defaut
   SSource=GetSImageFloat();	// valeur exacte
   SCont=GetSImageFloat();
   SReconst=GetSImageFloat();    
   pImage =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);
   pImage->Echelle = Echelle*2;

   SDest=pImage->GetSImageFloat(); 
   // --- recopie du buffer dans pImage afin de conserver le passé ... c'est beau, non ?
   // ImageFloat 
   for (i=0;i<(long)Taille;i++)
      SDest.Image[i]=SSource.Image[i];     	
   // Image         
   for (i=0;i<(long)Taille;i++)
      pImage->Image[i]=Image[i];     	

   if (menu.GetMenuState( ID_MULTI_ECHELLE,MF_BYCOMMAND) ==  MF_CHECKED)	
   {                                           
      Extraction(h,g,SSource,SDest,"HAAR (Float)");	
      SetCursor(LoadCursor(0,IDC_WAIT)); 

   }


   if (menu.GetMenuState( ID_RECONSTRUCT,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      rImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);  
      rImage->Echelle = Echelle;
      ok_reconst=TRUE;    
      SReconst=rImage->GetSImageFloat();
   }

   if (menu.GetMenuState( ID_CONT,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);  
      cImage->Echelle = 1;
      ok_cont=TRUE;      
      SCont=cImage->GetSImageFloat();
   }                        
   if (menu.GetMenuState( ID_BASSE,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);  
      cImage->Echelle = 1;
      ok_basse=TRUE;      
      SCont=cImage->GetSImageFloat();
   }                        


   TOD(h,g,SSource,SDest,SReconst,SCont,ok_reconst,ok_cont,ok_basse,Echelle);	
   pImage->FloatToImage();
   pImage->Create(Titre+" Analyse Multi-résolutions : HAAR (Float)",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 

   if (menu.GetMenuState(ID_RECONSTRUCT,MF_BYCOMMAND) ==  MF_CHECKED)		
   {
      rImage->FloatToImage();	
      rImage->Create(Titre+" Reconstruction : HAAR (Float)",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 
   }


   if (menu.GetMenuState( ID_BASSE,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage->FloatToImage();	
      cImage->Create(Titre+" Extraction des Vij : HAAR (Float)",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 
   }
   if (menu.GetMenuState( ID_CONT,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage->FloatToImage();	
      cImage->Create(Titre+" Extraction des détails : HAAR (Float)",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 
   }

   SetCursor(LoadCursor(0,IDC_ARROW));


   Invalidate(FALSE);
   delete h.Filtre;
   delete g.Filtre;		
}

void CImageWnd::OnSplineFloat()
{
   // TODO: Add your command handler code here
   SetCursor(LoadCursor(0,IDC_WAIT)); 

   // --- Filtres ---
   SFiltre h,g;
   SImageFloat SCont,SReconst,SSource,SDest; // Type inconnu deux structures possibles
   long i;
   BOOL ok_reconst=FALSE,ok_cont=FALSE,ok_basse=FALSE;		 
   CTraitImageWnd * pImage = NULL;	
   CTraitImageWnd * rImage = NULL;
   CTraitImageWnd * cImage = NULL;

   // ---------------- Initialisation des filtres ------------------
   //-- filtre associé à la fonction d'échelle --
   h.Filtre=new float[64*sizeof(float)];     
   h.Fin=63;h.Centre=32;

   h.Filtre[0]=0.000004F;h.Filtre[1]=-0.000001F;h.Filtre[2]=-0.000005F;h.Filtre[3]=0.000005F;h.Filtre[4]=0.000008F;h.Filtre[5]=-0.000011F;h.Filtre[6]=-0.000014F;
   h.Filtre[7]=0.000021F;h.Filtre[8]=0.000027F;h.Filtre[9]=-0.000042F;h.Filtre[10]=-0.000053F;h.Filtre[11]=0.000083F;h.Filtre[12]=0.000103F;h.Filtre[13]=-0.000164F;
   h.Filtre[14]=-0.000202F;h.Filtre[15]=0.000327F;h.Filtre[16]=0.000396F;h.Filtre[17]=-0.000656F;h.Filtre[18]=-0.000780F;h.Filtre[19]=0.001331F;h.Filtre[20]=0.001546F;
   h.Filtre[21]=-0.002745F;h.Filtre[22]=-0.003079F;h.Filtre[23]=0.005800F;h.Filtre[24]=0.006142F;h.Filtre[25]=-0.012715F;h.Filtre[26]=-0.012146F;h.Filtre[27]=0.029747F;
   h.Filtre[28]=0.022685F;h.Filtre[29]=-0.077808F;	h.Filtre[30]=-0.035498F;h.Filtre[31]=0.306830F;h.Filtre[32]=0.541736F;

   h.Filtre[63]=0.000001F;h.Filtre[62]=-0.000005F;h.Filtre[61]=0.000005F;h.Filtre[60]=0.000008F;h.Filtre[59]=-0.000011F; h.Filtre[58]=-0.000014F;h.Filtre[57]=0.000021F;h.Filtre[56]=0.000027F;
   h.Filtre[55]=-0.000042F;h.Filtre[54]=-0.000053F;h.Filtre[53]=0.000083F;h.Filtre[52]=0.000103F;h.Filtre[51]=-0.000164F;h.Filtre[50]=-0.000202F;h.Filtre[49]=0.000327F;
   h.Filtre[48]=0.000396F;h.Filtre[47]=-0.000656F;h.Filtre[46]=-0.000780F;h.Filtre[45]=0.001331F;h.Filtre[44]=0.001546F;h.Filtre[43]=-0.002745F;h.Filtre[42]=-0.003079F;
   h.Filtre[41]=0.005799F;h.Filtre[40]=0.006142F;h.Filtre[39]=-0.012715F;h.Filtre[38]=-0.012146F;h.Filtre[37]=0.029747F;h.Filtre[36]=0.022685F;h.Filtre[35]=-0.077808F;
   h.Filtre[34]=-0.035498F;h.Filtre[33]=0.306830F;


   //-- filtre associé à l'ondelette --	
   g.Filtre=new float[64*sizeof(float)];	
   g.Fin=63;g.Centre=32;
   for (i=0;i<=h.Fin;i++) g.Filtre[i] = (float)pow((double)-1,(double)i)*h.Filtre[h.Fin-i];


   // --------- Controle du menu -----------
   //
   // affectation par defaut
   SSource=GetSImageFloat();	// valeur exacte
   SCont=GetSImageFloat();
   SReconst=GetSImageFloat();    
   pImage =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);
   pImage->Echelle = Echelle*2;

   SDest=pImage->GetSImageFloat(); 
   // --- recopie du buffer dans pImage afin de conserver le passé ... c'est beau, non ?
   // ImageFloat 
   for (i=0;i<(long)Taille;i++)
      SDest.Image[i]=SSource.Image[i];     	
   // Image         
   for (i=0;i<(long)Taille;i++)
      pImage->Image[i]=Image[i];     	

   if (menu.GetMenuState( ID_MULTI_ECHELLE,MF_BYCOMMAND) ==  MF_CHECKED)	
   {                                           
      Extraction(h,g,SSource,SDest,"B-Splines cubiques (Float)");	
      SetCursor(LoadCursor(0,IDC_WAIT)); 

   }


   if (menu.GetMenuState( ID_RECONSTRUCT,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      rImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);  
      rImage->Echelle = Echelle;
      ok_reconst=TRUE;    
      SReconst=rImage->GetSImageFloat();
   }

   if (menu.GetMenuState( ID_CONT,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);  
      cImage->Echelle = 1;
      ok_cont=TRUE;      
      SCont=cImage->GetSImageFloat();
   }                        
   if (menu.GetMenuState( ID_BASSE,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);  
      cImage->Echelle = 1;
      ok_basse=TRUE;      
      SCont=cImage->GetSImageFloat();
   }                        


   TOD(h,g,SSource,SDest,SReconst,SCont,ok_reconst,ok_cont,ok_basse,Echelle);	
   pImage->FloatToImage();
   pImage->Create(Titre+" Analyse Multi-résolutions : B-Splines cubiques  (Float)",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 

   if (menu.GetMenuState( ID_RECONSTRUCT,MF_BYCOMMAND) ==  MF_CHECKED)		
   {
      rImage->FloatToImage();	
      rImage->Create(Titre+" Reconstruction : B-Splines cubiques  (Float)",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 
   }


   if (menu.GetMenuState( ID_CONT,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage->FloatToImage();	
      cImage->Create(Titre+" Extraction des détails : B-Splines cubiques  (Float)",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 
   }
   if (menu.GetMenuState( ID_BASSE,MF_BYCOMMAND) ==  MF_CHECKED)	
   {
      cImage->FloatToImage();	
      cImage->Create(Titre+" Extraction des Vij : B-Splines cubiques  (Float)",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 
   }

   SetCursor(LoadCursor(0,IDC_ARROW));


   Invalidate(FALSE);
   delete h.Filtre;
   delete g.Filtre;			
}


void CImageWnd::OnGradRoberts()
{
   SMasque masqGX,masqGY;
   // --- parametres du masque ----
   // GX
   masqGX.Filtre=new float[9];
   masqGX.Centre_col=1;
   masqGX.Centre_lig=1;
   masqGX.Haut=3;
   masqGX.Larg=3;
   masqGX.Filtre[0]=0.F;		masqGX.Filtre[1]=0.F;	masqGX.Filtre[2]=0.F;
   masqGX.Filtre[3]=0.F;		masqGX.Filtre[4]=1.F;	masqGX.Filtre[5]=0.F;
   masqGX.Filtre[6]=0.F;		masqGX.Filtre[7]=0.F;	masqGX.Filtre[8]=-1.F;
   // GY
   masqGY.Filtre=new float[9];
   masqGY.Centre_col=1;
   masqGY.Centre_lig=1;
   masqGY.Haut=3;
   masqGY.Larg=3;
   masqGY.Filtre[0]=0.F;		masqGY.Filtre[1]=0.F;	masqGY.Filtre[2]=0.F;
   masqGY.Filtre[3]=0.F;		masqGY.Filtre[4]=0.F;	masqGY.Filtre[5]=1.F;
   masqGY.Filtre[6]=0.F;		masqGY.Filtre[7]=-1.F;masqGY.Filtre[8]=0.F;


   Filtre_2D_GX_ET_GY(masqGX,masqGY,"Gradient de ROBERTS");   

   delete masqGX.Filtre;
   delete masqGY.Filtre;	

}

void CImageWnd::OnOpPrewitt()
{
   SMasque masqGX,masqGY;

   // --- parametres du masque ----
   // GX
   masqGX.Filtre=new float[9];
   masqGX.Centre_col=1;
   masqGX.Centre_lig=1;
   masqGX.Haut=3;
   masqGX.Larg=3;
   masqGX.Filtre[0]=-1.F;		masqGX.Filtre[1]=0.F;	masqGX.Filtre[2]=1.F;
   masqGX.Filtre[3]=-1.F;		masqGX.Filtre[4]=0.F;	masqGX.Filtre[5]=1.F;
   masqGX.Filtre[6]=-1.F;		masqGX.Filtre[7]=0.F;	masqGX.Filtre[8]=1.F;
   // GY
   masqGY.Filtre=new float[9];
   masqGY.Centre_col=1;
   masqGY.Centre_lig=1;
   masqGY.Haut=3;
   masqGY.Larg=3;
   masqGY.Filtre[0]=-1.F;	masqGY.Filtre[1]=-1.F;masqGY.Filtre[2]=-1.F;
   masqGY.Filtre[3]=0.F;		masqGY.Filtre[4]=0.F;	masqGY.Filtre[5]=0.F;
   masqGY.Filtre[6]=1.F;		masqGY.Filtre[7]=1.F;	masqGY.Filtre[8]=1.F;


   Filtre_2D_GX_ET_GY(masqGX,masqGY,"Opérateur de PREWITT");   

   delete masqGX.Filtre;
   delete masqGY.Filtre;		
}

void CImageWnd::OnLaplacien1()
{
   SMasque masq;

   // --- parametres du masque ----
   masq.Filtre=new float[9];
   masq.Centre_col=1;
   masq.Centre_lig=1;
   masq.Haut=3;
   masq.Larg=3;
   masq.Filtre[0]=0.F;		masq.Filtre[1]=1.F;	masq.Filtre[2]=0.F;
   masq.Filtre[3]=1.F;		masq.Filtre[4]=-4.F;	masq.Filtre[5]=1.F;
   masq.Filtre[6]=0.F;		masq.Filtre[7]=1.F;	masq.Filtre[8]=0.F;


   Filtre_2D(masq,"Laplacien 1");   

   delete masq.Filtre;

}

void CImageWnd::OnLaplacien2()
{
   SMasque masq;
   // --- parametres du masque ----
   masq.Filtre=new float[9];
   masq.Centre_col=1;
   masq.Centre_lig=1;
   masq.Haut=3;
   masq.Larg=3;
   masq.Filtre[0]=2.F;		masq.Filtre[1]=-1.F;	masq.Filtre[2]=2.F;
   masq.Filtre[3]=-1.F;		masq.Filtre[4]=-4.F;	masq.Filtre[5]=-1.F;
   masq.Filtre[6]=2.F;		masq.Filtre[7]=-1.F;	masq.Filtre[8]=2.F;


   Filtre_2D(masq,"Laplacien 2");   

   delete masq.Filtre;
}



void CImageWnd::OnBasse()
{
   if (menu.GetMenuState( ID_BASSE,MF_BYCOMMAND) ==  MF_UNCHECKED)
      menu.CheckMenuItem(ID_BASSE,MF_CHECKED);  
   else  menu.CheckMenuItem(ID_BASSE,MF_UNCHECKED);  	
}


void CImageWnd::OnCopieEnFloat()
{   
   long l=0;
   CTraitImageWnd * pImage = new CTraitImageWnd(fin_i-debut_i+1,fin_j-debut_j+1,FLOAT32);  
   for (long j=debut_j;j<=fin_j;j++)
      for (long i=debut_i;i<=fin_i;i++)		
      {   
         if (TypeImage==MEM || TypeImage==BMP) pImage->ImageFloat[l++] = (float) Image[j*Larg+i];
         else pImage->ImageFloat[l++] = ImageFloat[j*Larg+i];
      }	
      pImage->Create(Titre+" copie",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
}

void CImageWnd::OnCopieEnByte()
{  
   long l=0;
   CTraitImageWnd * pImage = new CTraitImageWnd(fin_i-debut_i+1,fin_j-debut_j+1,BMP);  
   for (long j=debut_j;j<=fin_j;j++)
      for (long i=debut_i;i<=fin_i;i++)		
      {
         pImage->Image[l++] = Image[j*Larg+i];        
      }	
      pImage->Create(Titre+" copie",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
}

void CImageWnd::OnMiseAZero()
{    
   SetCursor(LoadCursor(0,IDC_WAIT));
   Memorise();
   for (long j=debut_j;j<=fin_j;j++)
      for (long i=debut_i;i<=fin_i;i++)		
      {
         if (TypeImage==MEM || TypeImage==BMP) Image[j*Larg+i]=0;
         else  ImageFloat[j*Larg+i]=(float)0;
      }	
      if (TypeImage==FLOAT32) FloatToImage();
      SetCursor(LoadCursor(0,IDC_ARROW));
      Invalidate(FALSE);
}

void CImageWnd::OnEchantLigneCol()
{
   CTraitImageWnd * PImageEchant = new CTraitImageWnd(Larg/2,Haut/2,NB);
   SousEchantillonnageLigneColonne(GetSImage(),PImageEchant->GetSImage());
   PImageEchant->Create(Titre + " Sous-échantillonnage",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
}


void CImageWnd::OnMultiEchelle()
{   

   if (menu.GetMenuState(ID_MULTI_ECHELLE,MF_BYCOMMAND) ==  MF_UNCHECKED)
      menu.CheckMenuItem(ID_MULTI_ECHELLE,MF_CHECKED);  
   else  menu.CheckMenuItem(ID_MULTI_ECHELLE,MF_UNCHECKED);	

}

void CImageWnd::Filtre_2D(SMasque masq,CString Chaine)
{                                    
   // --- Déclarations des variables -----	



   float  * source;	
   float mini,maxi;	
   SImage info;

   // ---------- informations relatives a l'image -----
   info=GetSImage();

   CTraitImageWnd * pImage=NULL;	
   source = (float   *) malloc( Larg*Haut*sizeof(float));   

   float  * dest;
   dest = (float   *) malloc( Larg*Haut*sizeof(float));   
   if (TypeImage==FLOAT32) {  
      pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);
      recop_image(ImageFloat,source,Taille);
   }
   if (TypeImage==MEM || TypeImage==BMP){ 
      pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);

      recop_image(Image,source,Taille);	
   }
   SetCursor(LoadCursor(0,IDC_WAIT));
   // ---- Attention le filtre n'est tjrs pas retourné  hé hé ---

   masque_convolution_2D(source,dest,info,masq);  

   if (TypeImage==MEM || TypeImage==BMP){ 
      min_max(dest,&mini,&maxi,pImage->Taille);
      r_affichage(dest,pImage->Image,maxi,mini,Haut,Larg,1);	                    
   }                     
   if (TypeImage==FLOAT32) {
      recop_image(dest,pImage->ImageFloat,Taille);
      pImage->FloatToImage(); }
   pImage->Create(Titre+" : "+Chaine,WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   Invalidate(FALSE);
   free( source);
   free( dest);		
   SetCursor(LoadCursor(0,IDC_ARROW));	

}
void CImageWnd::Filtre_2D_GX_ET_GY(SMasque masqGX,SMasque masqGY,CString Chaine)
{                                    
   // --- Déclarations des variables -----	


   long i;
   float  * source;	
   float mini,maxi;	
   SImage info;

   // ---------- informations relatives a l'image -----
   info=GetSImage();

   CTraitImageWnd * pImage=NULL;	
   source = (float   *) malloc( Larg*Haut*sizeof(float));   

   float  * destGX;
   destGX = (float   *) malloc( Larg*Haut*sizeof(float));   
   float  * destGY;
   destGY = (float   *) malloc( Larg*Haut*sizeof(float));   
   float  * dest;
   dest = (float   *) malloc( Larg*Haut*sizeof(float));   

   if (TypeImage==FLOAT32) {  
      pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);
      recop_image(ImageFloat,source,Taille);
   }
   if (TypeImage==MEM || TypeImage==BMP){ 
      pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);

      recop_image(Image,source,Taille);	
   }
   SetCursor(LoadCursor(0,IDC_WAIT));  
   // ---- Attention le filtre n'est tjrs pas retourné  hé hé ---
   masque_convolution_2D(source,destGX,info,masqGX); 
   masque_convolution_2D(source,destGY,info,masqGY); 
   for(i=0;i<(long)Taille;i++)
      dest[i]=(float)fabs((double)destGX[i])+(float)fabs((double)destGY[i]);

   if (TypeImage==MEM || TypeImage==BMP){ 
      min_max(dest,&mini,&maxi,pImage->Taille);
      r_affichage(dest,pImage->Image,maxi,mini,Haut,Larg,1);	                    
   }                     
   if (TypeImage==FLOAT32) {
      recop_image(dest,pImage->ImageFloat,Taille);
      pImage->FloatToImage(); }
   pImage->Create(Titre+" : "+Chaine,WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   Invalidate(FALSE);	     

   free( source);
   free( dest);	
   SetCursor(LoadCursor(0,IDC_ARROW));	
}
void CImageWnd::Filtre_2D_GX_PUIS_GY(SMasque masqGX,SMasque masqGY,CString Chaine)
{                                    
   // --- Déclarations des variables -----	


   float  * source;	
   float mini,maxi;	
   SImage info;

   // ---------- informations relatives a l'image -----
   info=GetSImage();

   CTraitImageWnd * pImage=NULL;	
   source = (float   *) malloc( Larg*Haut*sizeof(float));   
   float  * dest;
   dest = (float   *) malloc( Larg*Haut*sizeof(float));   

   if (TypeImage==FLOAT32) {  
      pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);
      recop_image(ImageFloat,source,Taille);
   }
   if (TypeImage==MEM || TypeImage==BMP){ 
      pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);

      recop_image(Image,source,Taille);	
   }
   SetCursor(LoadCursor(0,IDC_WAIT));	
   // ---- Attention le filtre n'est tjrs pas retourné  hé hé ---
   masque_convolution_2D(source,dest,info,masqGX); 
   recop_image(dest,source,Taille);
   masque_convolution_2D(source,dest,info,masqGY);
   if (TypeImage==MEM || TypeImage==BMP){ 
      min_max(dest,&mini,&maxi,pImage->Taille);
      r_affichage(dest,pImage->Image,maxi,mini,Haut,Larg,1);	                    
   }                     
   if (TypeImage==FLOAT32) {
      recop_image(dest,pImage->ImageFloat,Taille);
      pImage->FloatToImage(); }
   pImage->Create(Titre+" : "+Chaine,WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   Invalidate(FALSE);	     

   free( source);
   free( dest);		
   SetCursor(LoadCursor(0,IDC_ARROW));	

}

void CImageWnd::Extraction(SFiltre h,SFiltre g,SImage SSource,SImage SDest,CString chaine)
{
   int cpt=0,prof,echelle;
   long i;	
   CDialogExtrac * Boite = new CDialogExtrac;
   SetCursor(LoadCursor(0,IDC_ARROW));
   Boite->DoModal(); 
   prof=Boite->m_profondeur;
   delete Boite;

   SetCursor(LoadCursor(0,IDC_WAIT));
   SImage SExtrac,SDecomp;
   // --- Allocation de la mémoire ---- 

   CTraitImageWnd * eImage = NULL;   
   CTraitImageWnd * sImage = NULL;
   eImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);
   eImage->Echelle = Echelle;	
   sImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);
   sImage->Echelle = Echelle;
   eImage->Create(Titre+" Extraction Multi-Echelles :"+chaine,WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 
   sImage->Create(Titre+" Décomposition (Extraction) : "+chaine,WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());   
   SExtrac=eImage->GetSImage();
   SDecomp=sImage->GetSImage();
   // ---- recopie de la source dans la SDecmpo ---
   for (i=0;i<(long)Taille;i++)
      SDecomp.Image[i]=SSource.Image[i];

   // --- Boite de dialogue ----

   float  * * tab_v_j = new float  * [prof*sizeof(float  *)];
   float  * * tab_d_h = new float  * [prof*sizeof(float  *)];
   float  * * tab_d_v = new float  * [prof*sizeof(float  *)];
   float  * * tab_d_o = new float  * [prof*sizeof(float  *)];
   for (i=0;i<prof;i++)
   {
      tab_d_h[i] = (float  *) malloc( Larg*Haut*sizeof(float));
      tab_d_v[i] = (float  *) malloc( Larg*Haut*sizeof(float));   
      tab_d_o[i] = (float  *) malloc( Larg*Haut*sizeof(float));   
      tab_v_j[i] = (float  *) malloc( Larg*Haut*sizeof(float));
   }             

   // ----- boucle d'Extraction des renseignements ---
   echelle=Echelle;
   for (cpt=0;cpt<prof;cpt++)   
   {
      TOD_EX (h,g,SDecomp,SDest,SExtrac,tab_v_j,tab_d_h,tab_d_v,tab_d_o,echelle,cpt);

      //SSource.Image=SDest.Image;
      for (i=0;i<(long)Taille;i++)
      {
         SDecomp.Image[i]=SDest.Image[i];
         SDest.Image[i]=0;
      }
      echelle=echelle*2;	
      sImage->Echelle=echelle; // --- affecter	
   }

   // ---- fin de la boucle d'Extraction ----


   // ---- Libération de la mémoire ---
   for (i=0;i<prof;i++)
   {
      free(tab_v_j[i]);
      free(tab_d_o[i]);
      free(tab_d_v[i]);
      free(tab_d_h[i]);
   }
   delete tab_v_j;
   delete tab_d_h;
   delete tab_d_o;	
   delete tab_d_v;
   eImage->Echelle=1;
   SetCursor(LoadCursor(0,IDC_ARROW));

}
void CImageWnd::Extraction(SFiltre h,SFiltre g,SImageFloat SSource,SImageFloat SDest,CString chaine)
{
   int cpt=0,prof,echelle;
   long i;	
   CDialogExtrac * Boite = new CDialogExtrac;
   SetCursor(LoadCursor(0,IDC_ARROW));
   Boite->DoModal(); 
   prof=Boite->m_profondeur;
   delete Boite;

   SetCursor(LoadCursor(0,IDC_WAIT));
   SImageFloat SExtrac,SDecomp;
   // --- Allocation de la mémoire ---- 

   CTraitImageWnd * eImage = NULL;   
   CTraitImageWnd * sImage = NULL;
   eImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);
   eImage->Echelle = Echelle;	
   sImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);
   sImage->Echelle = Echelle;
   eImage->Create(Titre+" Extraction Multi-Echelles :"+chaine,WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame()); 
   sImage->Create(Titre+" Décomposition (Extraction) : "+chaine,WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());   
   SExtrac=eImage->GetSImageFloat();
   SDecomp=sImage->GetSImageFloat();
   // ---- recopie de la source dans la SDecmpo ---
   for (i=0;i<(long)Taille;i++)
      SDecomp.Image[i]=SSource.Image[i];

   // --- Boite de dialogue ----

   float  * * tab_v_j = new float  * [prof*sizeof(float  *)];
   float  * * tab_d_h = new float  * [prof*sizeof(float  *)];
   float  * * tab_d_v = new float  * [prof*sizeof(float  *)];


   float  * * tab_d_o = new float  * [prof*sizeof(float  *)];
   for (i=0;i<prof;i++)
   {
      tab_d_h[i] = (float  *) malloc( Larg*Haut*sizeof(float));
      tab_d_v[i] = (float  *) malloc( Larg*Haut*sizeof(float));   
      tab_d_o[i] = (float  *) malloc( Larg*Haut*sizeof(float));   
      tab_v_j[i] = (float  *) malloc( Larg*Haut*sizeof(float));
   }             

   // ----- boucle d'Extraction des renseignements ---
   echelle=Echelle;
   for (cpt=0;cpt<prof;cpt++)   
   {
      TOD_EX (h,g,SDecomp,SDest,SExtrac,tab_v_j,tab_d_h,tab_d_v,tab_d_o,echelle,cpt);
      //SSource.Image=SDest.Image;	
      for (i=0;i<(long)Taille;i++)
      {
         SDecomp.Image[i]=SDest.Image[i];
         SDest.Image[i]=0.F;
      }
      echelle=echelle*2;
      sImage->Echelle=echelle;	
      sImage->FloatToImage();

   }
   //sImage->Echelle=echelle/2;	
   // ---- fin de la boucle d'Extraction ----


   // ---- Libération de la mémoire ---
   for (i=0;i<prof;i++)
   {
      free(tab_v_j[i]);
      free(tab_d_o[i]);
      free(tab_d_v[i]);
      free(tab_d_h[i]);
   }
   delete tab_v_j;
   delete tab_d_h;
   delete tab_d_o;
   delete tab_d_v;     
   eImage->Echelle=1;   		
   eImage->FloatToImage(); // -- On  normalise l'Extraction --
   SetCursor(LoadCursor(0,IDC_ARROW));

}

void CImageWnd::OnTraitementMorphologieErosion() 
{
   if ((Larg<3)||(Haut<3))
   {
      MessageBox("Image trop petite !","Erosion",MB_OK|MB_ICONSTOP);
      return;
   }
   CDialMorpho boite;
   if(boite.DoModal()==IDOK)
   {
      SetCursor(LoadCursor(0,IDC_WAIT));
      Erode(GetSImage(),boite.m_taillestructurant);
      SetCursor(LoadCursor(0,IDC_ARROW));  
      Invalidate();
      NomTraitement = "Erosion";
      CommentTraitement = "Erosion";
      dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
   }
}

void CImageWnd::OnTraitementMorphologieDilatation() 
{
   if ((Larg<3)||(Haut<3))
   {
      MessageBox("Image trop petite !","Dilatation",MB_OK|MB_ICONSTOP);
      return;
   }
   CDialMorpho boite;
   if(boite.DoModal()==IDOK)
   {
      SetCursor(LoadCursor(0,IDC_WAIT));
      Dilate(GetSImage(),boite.m_taillestructurant);
      SetCursor(LoadCursor(0,IDC_ARROW));
      Invalidate();
      NomTraitement = "Dilatation";
      CommentTraitement = "Dilatation";
      dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
   }
}

void CImageWnd::OnTraitementMorphologieOuverture() 
{
   if ((Larg<3)||(Haut<3))
   {
      MessageBox("Image trop petite !","Ouverture",MB_OK|MB_ICONSTOP);
      return;
   }
   CDialMorpho boite;
   if(boite.DoModal()==IDOK)
   {
      SetCursor(LoadCursor(0,IDC_WAIT));
      Erode(GetSImage(),boite.m_taillestructurant);
      Dilate(Im,boite.m_taillestructurant);
      SetCursor(LoadCursor(0,IDC_ARROW));
      Invalidate();
      NomTraitement = "Ouverture";
      CommentTraitement = "Ouverture";
      dialInfo->UpdateInfo(CommentTraitement,NomTraitement);
   }
}

void CImageWnd::OnTraitementMorphologieFermeture() 
{
   if ((Larg<3)||(Haut<3))
   {
      MessageBox("Image trop petite !","Fermeture",MB_OK|MB_ICONSTOP);
      return;
   }
   CDialMorpho boite;
   if(boite.DoModal()==IDOK)
   {
      SetCursor(LoadCursor(0,IDC_WAIT));

      Dilate(GetSImage(),boite.m_taillestructurant);
      Erode(Im,boite.m_taillestructurant);
      SetCursor(LoadCursor(0,IDC_ARROW));
      Invalidate();
      NomTraitement = "Fermeture";
      CommentTraitement = "Fermeture";
      dialInfo->UpdateInfo(CommentTraitement,NomTraitement);	
   }
}

void CImageWnd::OnOutilsPseudorelief() 
{
   float x,y,z,xp,yp;
   float f,d,al,be,ga,tet,co,si;
   long larg,haut;
   long i,j,k,l,u,v,w,kw;
   float X,Y,Z,XP,YP,ZP;
   float uo,vo,xo;

   CDialrelief Boite;   


   if(Boite.DoModal()==IDOK)
   {
      CTraitImageWnd *PImageFene = new CTraitImageWnd(fin_i-debut_i+1, fin_j-debut_j+1, NB);

      SetCursor(LoadCursor(0,IDC_WAIT));


      larg=fin_i-debut_i+1;    
      haut=fin_j-debut_j+1;

      al=0.3f;	// coefficient de hauteur
      be=2.3f;	// coefficient de distance
      ga=2.0f;	// coefficient de focale

      f=ga*(float)haut;
      d=be*(float)haut;

      xo = float((Boite.m_xo+0.5)*(float)larg);
      tet =float(Boite.m_tet*3.141592/180.0);
      be = Boite.m_dist;
      ga = Boite.m_foc;
      co=float(cos(tet));si=float(sin(tet));

      vo=float((float)haut/2.0);
      uo=xo;
      l=0;
      for(j=debut_j;j<=fin_j;j++)
      {
         for (i=debut_i;i<=fin_i;i++)
         {
            PImageFene->Image[l++]=0;
         }
      }

      for(j=fin_j;j>=debut_j;j--)
      {
         for(i=debut_i;i<=fin_i;i++)
         {
            w=Image[j*Larg+i];
            kw=(int)(al*w);
            for(k=0;k<=kw;k++)
            {
               X=(float)i;Y=(float)j;Z=(float)k;
               XP=X-xo;
               YP=Y+(float)haut*(float)(be*co-0.5);
               ZP=Z-(float)haut*be*si;
               x=XP;
               y=(ZP*si+YP*co);
               z=(ZP*co-YP*si);
               xp=-x*f/z;
               yp=-y*f/z;
               u=(int)(xp+uo);
               v=(int)(yp+vo);
               if((u>=0)&&(u<larg)&&(v>=0)&&(v<haut))
               {
                  PImageFene->Image[v*larg+u]=(int)(k/al);	
               }
            }
         }
      }
      PImageFene->Create(Titre+" Pseudo Relief", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   }	
}
void CImageWnd::OnTransformationConversionenbspline() 
{
   CTraitImageWnd *PImageFene = new CTraitImageWnd(Larg, Haut, FLOAT32);
   OnConvBsplineFloat(GetSImage(),PImageFene->GetSImageFloat());
   PImageFene->Create(Titre+"  B-Spline", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   SetCursor(LoadCursor(0,IDC_ARROW));	
}

void CImageWnd::OnExtraPaul() 
{
   if (menu.GetMenuState( ID_EXTRA_PAUL,MF_BYCOMMAND) ==  MF_UNCHECKED)
      menu.CheckMenuItem(ID_EXTRA_PAUL,MF_CHECKED);  
   else  menu.CheckMenuItem(ID_EXTRA_PAUL,MF_UNCHECKED); 
}

void CImageWnd::OnFileMruFile1() 
{
   CString NomFich;
   long Larg, Haut;
   CFile Fichier;
   int pos;
   NomFich = AfxGetApp()->GetProfileString( "Recent File List",  "File1");
   if (Fichier.Open(NomFich, CFile::modeRead))
   { 
      CTraitImageWnd	* pImageWnd;
      Fichier.Seek(sizeof(BITMAPFILEHEADER)+sizeof(DWORD), CFile::begin);
      Fichier.Read(&Larg, sizeof(long));
      Fichier.Read(&Haut, sizeof(long));
      Fichier.Close();
      pImageWnd = new CTraitImageWnd(NomFich, Larg, Haut, BMP);
      pos=NomFich.ReverseFind(0x5c);
      NomFich=NomFich.Right(NomFich.GetLength()-pos-1);
      pImageWnd->Create(NomFich.Left(NomFich.GetLength()-4), WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());

   }
   else MessageBox("Ce fichier n'existe pas","Message",MB_OK|MB_ICONEXCLAMATION);
}

void CImageWnd::OnFileMruFile2() 
{
   // TODO: Add your command handler code here
   CString NomFich;
   long Larg, Haut;
   CFile Fichier;
   int pos;
   NomFich = AfxGetApp()->GetProfileString( "Recent File List",  "File2");
   if (Fichier.Open(NomFich, CFile::modeRead))
   { 
      CString Temp = AfxGetApp()->GetProfileString( "Recent File List",  "File1");
      AfxGetApp()->AddToRecentFileList(NomFich);
      AfxGetApp()->WriteProfileString( "Recent File List",  "File1", NomFich);
      AfxGetApp()->WriteProfileString( "Recent File List",  "File2", Temp);
      CTraitImageWnd	* pImageWnd;
      Fichier.Seek(sizeof(BITMAPFILEHEADER)+sizeof(DWORD), CFile::begin);
      Fichier.Read(&Larg, sizeof(long));
      Fichier.Read(&Haut, sizeof(long));
      Fichier.Close();
      pImageWnd = new CTraitImageWnd(NomFich, Larg, Haut, BMP);
      pos=NomFich.ReverseFind(0x5c);
      NomFich=NomFich.Right(NomFich.GetLength()-pos-1);
      pImageWnd->Create(NomFich.Left(NomFich.GetLength()-4), WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());

   }
   else MessageBox("Ce fichier n'existe pas","Message",MB_OK|MB_ICONEXCLAMATION);
}

void CImageWnd::OnFileMruFile3() 
{
   CString NomFich;
   long Larg, Haut;
   CFile Fichier;
   int pos;
   NomFich = AfxGetApp()->GetProfileString( "Recent File List",  "File3");
   if (Fichier.Open(NomFich, CFile::modeRead))
   { 
      CString Temp = AfxGetApp()->GetProfileString( "Recent File List",  "File1");
      AfxGetApp()->AddToRecentFileList(NomFich);
      AfxGetApp()->WriteProfileString( "Recent File List",  "File1", NomFich);
      AfxGetApp()->WriteProfileString( "Recent File List",  "File3", AfxGetApp()->GetProfileString( "Recent File List",  "File2"));
      AfxGetApp()->WriteProfileString( "Recent File List",  "File2", Temp);
      CTraitImageWnd	* pImageWnd;
      Fichier.Seek(sizeof(BITMAPFILEHEADER)+sizeof(DWORD), CFile::begin);
      Fichier.Read(&Larg, sizeof(long));
      Fichier.Read(&Haut, sizeof(long));
      Fichier.Close();
      pImageWnd = new CTraitImageWnd(NomFich, Larg, Haut, BMP);
      pos=NomFich.ReverseFind(0x5c);
      NomFich=NomFich.Right(NomFich.GetLength()-pos-1);
      pImageWnd->Create(NomFich.Left(NomFich.GetLength()-4), WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());

   }
   else MessageBox("Ce fichier n'existe pas","Message",MB_OK|MB_ICONEXCLAMATION);

}

void CImageWnd::OnFileMruFile4() 
{
   CString NomFich;
   long Larg, Haut;
   CFile Fichier;
   int pos;
   NomFich = AfxGetApp()->GetProfileString( "Recent File List",  "File4");
   if (Fichier.Open(NomFich, CFile::modeRead))
   { 
      AfxGetApp()->AddToRecentFileList(NomFich);
      CString Temp = AfxGetApp()->GetProfileString( "Recent File List",  "File1");
      AfxGetApp()->WriteProfileString( "Recent File List",  "File1", NomFich);
      AfxGetApp()->WriteProfileString( "Recent File List",  "File4", AfxGetApp()->GetProfileString( "Recent File List",  "File3"));
      AfxGetApp()->WriteProfileString( "Recent File List",  "File3", AfxGetApp()->GetProfileString( "Recent File List",  "File2"));
      AfxGetApp()->WriteProfileString( "Recent File List",  "File2", Temp);
      CTraitImageWnd	* pImageWnd;
      Fichier.Seek(sizeof(BITMAPFILEHEADER)+sizeof(DWORD), CFile::begin);
      Fichier.Read(&Larg, sizeof(long));
      Fichier.Read(&Haut, sizeof(long));
      Fichier.Close();
      pImageWnd = new CTraitImageWnd(NomFich, Larg, Haut, BMP);
      pos=NomFich.ReverseFind(0x5c);
      NomFich=NomFich.Right(NomFich.GetLength()-pos-1);
      pImageWnd->Create(NomFich.Left(NomFich.GetLength()-4), WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());

   }
   else MessageBox("Ce fichier n'existe pas","Message",MB_OK|MB_ICONEXCLAMATION);

}

void CImageWnd::OnUpdateHorizontale(CCmdUI* pCmdUI) 
{
   if (TypeColor != NB) pCmdUI->Enable(FALSE);
}

void CImageWnd::OnUpdateOblique(CCmdUI* pCmdUI) 
{
   if (TypeColor != NB) pCmdUI->Enable(FALSE);
}

void CImageWnd::OnUpdateVerticale(CCmdUI* pCmdUI) 
{
   if (TypeColor != NB) pCmdUI->Enable(FALSE);
}

void CImageWnd::OnGrab() 
{
   ((CMainFrame*) GetMDIFrame())->OnAcquerirContinue();	
}

void CImageWnd::OnUpdateGrab(CCmdUI* pCmdUI) 
{
   ((CMainFrame*) GetMDIFrame())->OnUpdateGrab( pCmdUI); 
}

void CImageWnd::OnSnap() 
{
   ((CMainFrame*) GetMDIFrame())->OnAcquerirImage();
}

void CImageWnd::OnUpdateSnap(CCmdUI* pCmdUI) 
{
   ((CMainFrame*) GetMDIFrame())->OnUpdateSnap( pCmdUI); 
}

void CImageWnd::OnShowWindow(BOOL bShow, UINT nStatus) 
{
   CMDIChildWnd::OnShowWindow(bShow, nStatus);

   if(PALETTE_AFFICHEE) OnPaletteAffiche();	
}

void CImageWnd::OnTraitementNormalisationdescomposantes() 
{
   float som;
   for(long i=0;i<Larg*Haut;i++)
   {
      if(ImageRGB[i].r!=0 ||ImageRGB[i].g!=0||ImageRGB[i].b!=0)
      {
         som=(float)ImageRGB[i].r+(float)ImageRGB[i].g+(float)ImageRGB[i].b;
         ImageRGB[i].r=(unsigned char)(255.0*(float)ImageRGB[i].r/som);
         ImageRGB[i].g=(unsigned char)(255.0*(float)ImageRGB[i].g/som);
         ImageRGB[i].b=(unsigned char)(255.0*(float)ImageRGB[i].b/som);
      }
   }
   Invalidate();
}

void CImageWnd::OnTraitementPetrou() 
{

}

void CImageWnd::OnTraitementContourpetrou() 
{
#ifdef MP
   ImagePetrou IR,B,I_x, I_y;
   rgb RGB;
   int SeuilPetrou=4;
   int SigmaGaussien=3;
   BYTE * ImageGauss = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);

   CTraitImageWnd *pImagePetrou=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);


   //--------------------------------------------------------
   //			Init pour Détecteur contour Petrou
   //--------------------------------------------------------
   int i,j;
   IR.dim_x =Larg;
   IR.dim_y =Haut;

   B.dim_x = Larg;
   B.dim_y = Haut;
   I_x.dim_x = Larg;
   I_x.dim_y = Haut;
   I_y.dim_x = Larg;
   I_y.dim_y = Haut;

   IR.ima = (Pixelp**) malloc(sizeof(Pixelp*)*IR.dim_y);

   B.ima = (Pixelp**) malloc(sizeof(Pixelp*)*B.dim_y);
   I_x.ima = (Pixelp**) malloc(sizeof(Pixelp*)*I_x.dim_y);
   I_y.ima = (Pixelp**) malloc(sizeof(Pixelp*)*I_y.dim_y);

   for(i=0;i<IR.dim_y;i++)
   {
      IR.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*IR.dim_x);
      B.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*B.dim_x);
      I_x.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*I_x.dim_x);
      I_y.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*I_y.dim_x);
   }

   Convolution2DCedricCurvi(ImageGauss,Image,SigmaGaussien,Larg,Haut);
   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         IR.ima[i][j] = (Pixelp)ImageGauss[i*Larg+j]; // niveau de gris
      }
   }

   CPetrou ContourPetrou(&IR, &B, &I_x,&I_y);
   ContourPetrou.filtre_petrou(&IR,&B,6,SeuilPetrou,1);
   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         pImagePetrou->ImageRGB[i*Larg+j].r = (BYTE) B.ima[i][j];
         pImagePetrou->ImageRGB[i*Larg+j].g = (BYTE) B.ima[i][j];
         pImagePetrou->ImageRGB[i*Larg+j].b = (BYTE) B.ima[i][j];
      }
   }

   pImagePetrou->Create(Titre + " : Petrou",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());

   for(i=0;i<IR.dim_y;i++)
   {
      free(IR.ima[i]);
      free(B.ima[i]);
      free(I_x.ima[i]);
      free(I_y.ima[i]);
   }

   free(IR.ima);
   free(B.ima);
   free(I_x.ima);
   free(I_y.ima);

   free(ImageGauss);

#endif
}




void CImageWnd::OnTraitementConvolutiongaussien() 
{
   // TODO: Add your command handler code here
   CDialVal DialVal;
   double sigma;
   if (DialVal.DoModal() != IDOK) return;
   if (DialVal.m_val=="") sigma=0.1;
   else sigma=atof(DialVal.m_val);
   SMasque masq;
   int l,c,p;
   double deuxpi=6.28318;
   double sigma2=sigma*sigma;
   double unsdeuxsigc=1.0/(2.0*sigma2);
   double norm=1.0/(deuxpi*sigma2);
   double norm2=1.0/(deuxpi*sigma2*sigma2);
   double norm3=1.0/(deuxpi*sigma2*sigma2*sigma2);
   double x,y;
   double arg;

   // --- parametres du masque ----
   masq.Haut=19;
   masq.Larg=19;
   masq.Centre_col=masq.Larg/2;
   masq.Centre_lig=masq.Haut/2;
   masq.Filtre=new float[masq.Haut*masq.Larg];

   for (l=0;l<masq.Haut;l++)
   {
      for (c=0;c<masq.Larg;c++)
      {
         p=masq.Larg*l+c;
         x=double(c-masq.Centre_col);
         y=double(l-masq.Centre_lig);
         arg=(x*x+y*y)*unsdeuxsigc;
         // gaussien
         masq.Filtre[p]=(float)(norm*exp(-arg));
         // dérivée de gaussien suivant c
         //	masq.Filtre[p]=(float)(-x*norm2*exp(-arg));
         //	derivee seconde
         //	masq.Filtre[p]=(float)(-x*norm2*exp(-arg));

      }
   }

   Filtre_2D(masq,"gaussien X");


   delete masq.Filtre;
}

void CImageWnd::OnTraitementConvolutiongaussienGaussien1dx() 
{

}

void CImageWnd::OnTraitementConvolutiongaussienGaussien1dy() 
{
   // TODO: Add your command handler code here

}

void CImageWnd::OnTraitementConvolutiongaussienCurvilineardetection() 
{
   // TODO: Add your command handler code here
   SMasque masq;
   int l,c,p;
   float mini,maxi;	
   SImage info;
   double deuxpi=6.28318;
   double sigma=10;
   double sigma2=sigma*sigma;
   double unsdeuxsigc=1.0/(2.0*sigma2);
   double norm=1.0/(deuxpi*sigma2);
   double norm2=1.0/(deuxpi*sigma2*sigma2);
   double norm3=1.0/(deuxpi*sigma2*sigma2*sigma2);
   double x,y;
   double arg;

   double Lambda1,Lambda2,LambdaMax;
   double Nx,Ny;
   double t;
   double Px,Py;


   // --- parametres du masque ----
   masq.Haut=19;
   masq.Larg=19;
   masq.Centre_col=masq.Larg/2;
   masq.Centre_lig=masq.Haut/2;
   masq.Filtre=new float[masq.Haut*masq.Larg];

   float  * source;


   float  * Rx;	
   float  * Ry;	
   float  * Rxx;	
   float  * Rxy;	
   float  * Ryy;	
   float  * Lbd;	

   info=GetSImage();

   CTraitImageWnd * pImage=NULL;	
   Rx = (float   *) malloc(Larg*Haut*sizeof(float));   
   Ry = (float   *) malloc(Larg*Haut*sizeof(float));   
   Rxx = (float   *) malloc(Larg*Haut*sizeof(float));   
   Rxy = (float   *) malloc(Larg*Haut*sizeof(float));   
   Ryy = (float   *) malloc(Larg*Haut*sizeof(float));   
   Lbd = (float   *) malloc(Larg*Haut*sizeof(float));   
   source = (float   *) malloc(Larg*Haut*sizeof(float));   

   if (TypeImage==MEM || TypeImage==BMP)
   { 
      pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);
      recop_image(Image,source,Taille);	
   }

   for (l=0;l<masq.Haut;l++)
   {
      for (c=0;c<masq.Larg;c++)
      {
         p=masq.Larg*l+c;
         x=double(c-masq.Centre_col);
         y=double(l-masq.Centre_lig);
         arg=(x*x+y*y)*unsdeuxsigc;
         masq.Filtre[p]=(float)(-x*norm2*exp(-arg));

      }
   }

   masque_convolution_2D(source,Rx,info,masq);  
   //	Filtre_2D(masq,"gaussien X");

   for (l=0;l<masq.Haut;l++)
   {
      for (c=0;c<masq.Larg;c++)
      {
         p=masq.Larg*l+c;
         x=double(c-masq.Centre_col);
         y=double(l-masq.Centre_lig);
         arg=(x*x+y*y)*unsdeuxsigc;
         masq.Filtre[p]=(float)(-y*norm2*exp(-arg));

      }
   }
   masque_convolution_2D(source,Ry,info,masq);  

   for (l=0;l<masq.Haut;l++)
   {
      for (c=0;c<masq.Larg;c++)
      {
         p=masq.Larg*l+c;
         x=double(c-masq.Centre_col);
         y=double(l-masq.Centre_lig);
         arg=(x*x+y*y)*unsdeuxsigc;
         masq.Filtre[p]=(float)((x*x-sigma2)*norm3*exp(-arg));

      }
   }
   masque_convolution_2D(source,Rxx,info,masq);  

   for (l=0;l<masq.Haut;l++)
   {
      for (c=0;c<masq.Larg;c++)
      {
         p=masq.Larg*l+c;
         x=double(c-masq.Centre_col);
         y=double(l-masq.Centre_lig);
         arg=(x*x+y*y)*unsdeuxsigc;
         masq.Filtre[p]=(float)(y*x*norm3*exp(-arg));
      }
   }
   masque_convolution_2D(source,Rxy,info,masq);  

   for (l=0;l<masq.Haut;l++)
   {
      for (c=0;c<masq.Larg;c++)
      {
         p=masq.Larg*l+c;
         x=double(c-masq.Centre_col);
         y=double(l-masq.Centre_lig);
         arg=(x*x+y*y)*unsdeuxsigc;
         masq.Filtre[p]=(float)((y*y-sigma2)*norm3*exp(-arg));
      }
   }
   masque_convolution_2D(source,Ryy,info,masq);  

   double temp;
   double denom;

   for (l=0;l<Haut;l++)
   {
      for (c=0;c<Larg;c++)
      {
         p=Larg*l+c;
         temp=(Rxx[p]-Ryy[p])*(Rxx[p]-Ryy[p]);
         Lambda1=(Rxx[p]+Ryy[p]+sqrt(temp+4.0*Rxy[p]*Rxy[p]))*0.5;
         Lambda2=(Rxx[p]+Ryy[p]-sqrt(temp+4.0*Rxy[p]*Rxy[p]))*0.5;
         if(fabs(Lambda1)>fabs(Lambda2)) LambdaMax=Lambda1; else LambdaMax=Lambda2;
         Lbd[p]=fabs(LambdaMax);
         temp=(LambdaMax-Rxx[p])*(LambdaMax-Rxx[p]);
         denom=sqrt(Rxy[p]*Rxy[p]+temp);
         Nx=(Rxy[p])/denom;
         Ny=(LambdaMax-Rxx[p])/denom;
         t=-(Rx[p]*Nx+Ry[p]*Ny);
         t=t/( Rxx[p]*Nx*Nx  +  2.0*Rxy[p]*Nx*Ny  +  Ryy[p]*Ny*Ny);
         Px=t*Ny;
         Py=t*Ny;
         if ((Px>=-0.5 && Px<=0.5) && (Py>=-0.5 && Py<=0.5) && fabs(LambdaMax)>0.8)
         {
            Image[p]=250; 
         }
         //else Image[p]=10;

      }
   }



   if (TypeImage==MEM || TypeImage==BMP)
   { 
      min_max(Lbd,&mini,&maxi,pImage->Taille);
      r_affichage(Lbd,pImage->Image,maxi,mini,Haut,Larg,1);	                    
   }                     
   pImage->Create(Titre+" : "+"Lbd",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   Invalidate(FALSE);


   delete masq.Filtre;	
   free(source);
   free(Rx);
   free(Ry);
   free(Rxx);
   free(Rxy);
   free(Ryy);

}

void CImageWnd::OnTraitementCurvilinear() 
{
#ifdef MP
   // TODO: Add your command handler code here
   SImage info;
   info=GetSImage();
   CDialCurvi DialCurvi;
   DialCurvi.Larg=info.Larg;
   if (DialCurvi.DoModal()==IDCANCEL) return;

   CPetrou Petrou;

   Petrou.m_TE=DialCurvi.m_TE;
   Petrou.m_MaxLength=DialCurvi.m_MaxLength;
   Petrou.m_MaxWidth=DialCurvi.m_MaxWidth;
   Petrou.m_MaxDw=DialCurvi.m_MaxDw;
   Petrou.m_CheckDw=DialCurvi.m_CheckDw;
   Petrou.m_MinLength = DialCurvi.m_MinLength;	
   Petrou.m_MaxReg=DialCurvi.m_Reg;
   Petrou.m_Check_Reg=DialCurvi.m_Check_Reg;
   // ---------- informations relatives a l'image -----

   CTraitImageWnd * pImage=NULL;	
   pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);
   SetCursor(LoadCursor(0,IDC_WAIT));

   Petrou.Go(Image,pImage->Image,info.Haut,info.Larg,DialCurvi.m_TE,1);  

   CTraitImageWnd * pImageMix=NULL;	
   pImageMix=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);
   for (int p=0; p<Larg*Haut;p++)
   {
      if (Image[p]>50) pImageMix->Image[p]=Image[p]-50;
      else pImageMix->Image[p]=0;
      if (pImage->Image[p]!=0) pImageMix->Image[p]=255;//pImage->Image[p];
      //		 else pImageMix->Image[p]=Image[p];

   }


   pImage->Create(Titre+" Curvilinear ",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   pImageMix->Create(Titre+" Curvilinear Mix",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   Invalidate(FALSE);
   SetCursor(LoadCursor(0,IDC_ARROW));			
#endif
}

void CImageWnd::OnColorCurvi() 
{
   // TODO: Add your command handler code here

   CCurviDetect Curvi;
   Curvi.CompteurFFT=0;

   CTraitImageWnd * pImage=NULL;	
   pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);
   Curvi.Mode=4;
   Curvi.Go(ImageRGB,pImage->Image,Haut,Larg);
   /*CString temp;
   temp.Format("%d",Curvi.CompteurFFT);
   MessageBox(temp,"Nombre de FFT",MB_OK);*/

   pImage->Create(Titre+" Curvilinear color",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   Invalidate(FALSE);
}

void CImageWnd::OnRButtonDblClk(UINT nFlags, CPoint point) 
{
   // TODO: Add your message handler code here and/or call default
   CCurviDetect Curvi;
   CPoint pointt;
   pointt=Get_ij(point.x,point.y);

   CTraitImageWnd * pImage=NULL;	
   pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);
   Curvi.Mode=4;
   Curvi.Go(ImageRGB,pImage->Image,Haut,Larg,pointt); 
   Invalidate(FALSE);
   Curvi.Info();


   CMDIChildWnd::OnRButtonDblClk(nFlags, point);
}

void CImageWnd::OnGradY() 
{
   // TODO: Add your command handler code here
   CCurviDetect Curvi;

   CTraitImageWnd * pImage=NULL;	
   pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);
   Curvi.Mode=4;
   Curvi.Test(ImageRGB,pImage->Image,Haut,Larg);  

   pImage->Create(Titre+" Curvilinear color",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   Invalidate(FALSE);
   /*
   POINT point[50];
   for (int i=0;i<50;i++) {point[i].y=-0.34*i+280.6;point[i].x=i;}
   CRegress Regress(point,50);
   CString Temp;
   Temp.Format("a=%f    b=%f   r=%f",Regress.a,Regress.b,Regress.r);
   MessageBox(Temp);
   */

}

void CImageWnd::OnTraitementEstimationtempsfft() 
{
   CString NomImage;
   LARGE_INTEGER Debut,Fin,Freq;

   //---------------------------
   //CIma pour la FFT
   //---------------------------

   CIma MonImage;

   __int64 tps;
   CString Temp;
   int LargeurBordFourier = 16;
   int ifft=9999;
   double * TabEstimeTempsx = (double *)malloc(sizeof(double)*LargeurBordFourier);
   double * TabEstimeTempsy = (double *)malloc(sizeof(double)*LargeurBordFourier);
   for (int i=0;i<LargeurBordFourier;i++)
   {
      TabEstimeTempsx[i] = 255.0f;
      TabEstimeTempsy[i] = 255.0f;
   }
   QueryPerformanceFrequency(&Freq);
   QueryPerformanceCounter(&Debut);

   while (ifft!=0)
   {
      MonImage.FFT(1,4,TabEstimeTempsx,TabEstimeTempsy);
      ifft--;
   }



   QueryPerformanceCounter(&Fin);

   free(TabEstimeTempsx);
   free(TabEstimeTempsy);
   tps=Fin.QuadPart-Debut.QuadPart;
   tps=1000*tps;
   tps=tps/Freq.QuadPart;
   Temp.Format("%d",tps);
   MessageBox(Temp,"temps pour une FFT 16 valeurs",MB_OK);

}

void CImageWnd::OnTraitementHypergradientParbloc() 
{
   BYTE * imagealac = (BYTE *)malloc(sizeof(BYTE)*(Haut/16)*(Larg/16)*3);

   CIma ClasseImage(imagealac,Haut/16,Larg/16,100);
   CImageWnd * pImage=NULL;	


   pImage=(CTraitImageWnd * ) new CImageWnd(Larg/16,Haut/16, COLOR);
   ClasseImage.HyperGradientDescripteurBloc(pImage->ImageRGB,ImageRGB,16,16,Larg,Haut);

   pImage->Create(Titre+" Curvilinear color",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   Invalidate(FALSE);	
}

void CImageWnd::OnTraitementHypergradientFenetreglissante() 
{
   BYTE * imagealac = (BYTE *)malloc(sizeof(BYTE)*Haut*Larg*3);

   CIma ClasseImage(imagealac,Haut,Larg,100);

   CImageWnd * pImage=NULL;	


   pImage=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);

   ClasseImage.HyperGradientDescripteurGlissant(pImage->ImageRGB,ImageRGB,16,16,Larg,Haut);

   pImage->Create(Titre+" Curvilinear color",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   Invalidate(FALSE);
}

void CImageWnd::OnTraitementHypergradientFenetreglissantezeropadding() 
{
   BYTE * imagealac = (BYTE *)malloc(sizeof(BYTE)*Haut*Larg*3);

   CIma ClasseImage(imagealac,Haut,Larg,100);

   CImageWnd * pImage=NULL;	


   pImage=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);

   ClasseImage.HyperGradientGlissantZeroPadding(pImage->ImageRGB,ImageRGB,32,32,Larg,Haut);

   pImage->Create(Titre+" Curvilinear color",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   Invalidate(FALSE);	
}

void CImageWnd::OnTraitementOutilsEtudedistancemaxi() 
{
   // Pour faire une étude statistique sur les distances maxi pour la normalisation
   CCurviDetect Curvi;
   Curvi.CompteurFFT=0;

   CTraitImageWnd * pImage=NULL;	
   pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);
   Curvi.Mode=4;
   Curvi.Go(ImageRGB,pImage->Image,Haut,Larg);


   pImage->Create(Titre+" Curvilinear color",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   Invalidate(FALSE);	
}

void CImageWnd::OnOutilsSuividecontours()
{
   if (BoiteReglageParam.DoModal()==IDCANCEL) return ;
   
   CCurvLinearParams p;

   p.SeuilReponse = BoiteReglageParam.SeuilReponse;
   p.LargeurEtude = BoiteReglageParam.ValLargeur;
   p.LargeurFourier = BoiteReglageParam.ValLargeurFourier;
   p.LongueurMiniForme = BoiteReglageParam.NbListeMini;
   p.LargeurMini = BoiteReglageParam.ValLargeurMini;
   p.VLocaleLargeur = BoiteReglageParam.m_EditVariationLocaleLargeur;

   BOOL DebugEgde = BoiteReglageParam.m_BDebugEdges;
   BOOL DebugTracking = BoiteReglageParam.m_BDebugTracking;

   SImage info = GetSImage();
   
   utls::BAry *tmp = new utls::BAry(Haut, Larg, this->Image);
   typedef utls::Ary<float, float, double> CurvImage;
   typedef utls::Ary<utls::RGBValue<float>, float, utls::RGBValue<double> > CurvRGBImage;

   CurvImage *im = utls::convert<utls::BAry, CurvImage>(tmp);
   utls::BAry *edge;
   TableOfCurvilinearRegions TableauListesPointsFiltrees, TableauListesPointsFiltreesC2D;
   DetectCurvLinearPoints(im, p, 0, &edge, TableauListesPointsFiltrees, TableauListesPointsFiltreesC2D);

   tmp->detach(); delete tmp;
#ifdef OLD
   bool bAff=TRUE; // for progress dialog  false for debug and ANIM

   CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> TableauListesPointsFiltrees;
   CList <PosDetectCurvi,PosDetectCurvi> * ListePointsFiltres;
   CList <PosDetectCurvi,PosDetectCurvi> * ListePointsFiltresComp;
   CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> TableauListesPointsFiltreesC2D;
   CList <PosDetectCurvi,PosDetectCurvi> * ListePointsFiltresC2D;
   CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> TableauListesFinal;

   CArray <CList<int,int>*,CList<int,int>*> TabNumListesConnex;
   CArray <SBoundingBox,SBoundingBox> TableauListesBoundingBox;
   SImage info;

   CCurviDetect2D CurviDetect2D;
   CRgn * ConstrucRegion;
   CPoint * TabTraceRegion;
   int MaxiLongueurListe;
   int IndiceListeMaxi;
   int IndiceReg;

   float SigneDx = 0.0f;
   int Sens=-1;
   float OldSigneDx=0.0f;

   BOOL Init = FALSE;
   PosDetectCurvi PosNM1,PosN,PosNP1;
   double MoyPrev, MoyNext;
   POSITION PositionN, PositionNM1,PositionCListeCopy;
   int CptTableauListe = 0;
   Taille=Larg*Haut;

   CFenetreProgressCurvi * AffichageProgression;
   AffichageProgression = new CFenetreProgressCurvi(this);
   AffichageProgression->Create(IDD_DIALOG_PROGRESSCURVI);
   LARGE_INTEGER  FreqProc, TpsStart, TPetrou, TTracking, T1DProcess, T2DProcess, TAffichage;
   LARGE_INTEGER  tps1,tps2,tps3,tps4,tps5;
   QueryPerformanceFrequency(&FreqProc);

   int m_LongueurMiniRegionFiltrage = m_LongueurMiniForme/2;
   m_LongueurMiniForme=20; // debug
   long CompteurRegionFiltrage = 0;

   int SeuilPetrou=2;
   int SigmaGaussien=3;

   int SeuilDepass=5;
   int SeuilExtrapol=5;

   BOOL BoiteAxe=TRUE;


   int TabHisto[1000];

   CPoint p1,p2,p3;
   Taille=Larg*Haut;


   ImagePetrou A,APasGaussien/*Contient l'image original*/,B,I_x, I_y;
   int NbListe;
   int i=0,j=0,k=0;

   int w = m_LargeurEtude;
   int dw = w;
   int cpti = 0;
   int cptpoint = 0;
   int LargeurFourier = m_LargeurFourier;

   int NBPts = 2*(w+dw)+ 1;
   int Centre;
   POSITION pos;
   CPoint CI1, CI2, CI3;
   BOOL bCI1 = FALSE,bCI2 = FALSE,bCI3 = FALSE;
   float GradMoy_x = 0.0f, GradMoy_y = 0.0f;
   CPoint CurentPoint;

   /*BYTE * TabGauche = (BYTE *)malloc(sizeof(BYTE)*((NBPts/2)+2*LargeurFourier));
   BYTE * TabDroite = (BYTE *)malloc(sizeof(BYTE)*(((NBPts/2)+1+2*LargeurFourier)));*/

   BYTE * TabGauche = (BYTE *)malloc(sizeof(BYTE)*(w+1+LargeurFourier));
   BYTE * TabDroite = (BYTE *)malloc(sizeof(BYTE)*(w+1+LargeurFourier));

   CPoint * TabCoupe = (CPoint*) malloc(sizeof(CPoint)* NBPts);
   CPoint * TabCoupeGauche   = new CPoint[w+1+LargeurFourier];
   for(i=0;i<NBPts;i++)
   {
      TabCoupe[i].x=0;
      TabCoupe[i].y=0;
   }

   BYTE * ImageFreeman = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);
   BYTE * ImageGauss = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);
   int *ImageLabelEdge=(int*) malloc(sizeof(int)*Larg*Haut);

   long p=0;
   int cptCoul=0;
   int MaxListe = 0, PosMaxListe = 0;
   int CptPtsGrad = 0;
   BYTE r,g,b;

   int ValIptInitFixe;


   PosDetectCurvi PosGauche;
   PosDetectCurvi PosDroite;
   PosDetectCurvi PosARetenir;

   CString ChaineVisu;
   //CFile FichierVisu;
   CStdioFile FichierVisu;
   CImageWnd * pImage=NULL;	
   CImageWnd * pImagePetrou=NULL;
   CImageWnd * pImageTracking=NULL;

   CPoint TabPourQuadri[4];




   //--------------------------------------------------------------------
   //			construction object ccurvidetect1d
   //--------------------------------------------------------------------

   CCurviDetect1D CurviDetect1D(w+1+LargeurFourier,LargeurFourier,LargeurFourier);
   // ---------- informations relatives a l'image -----
   info=GetSImage();


   pImage=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);
   pImagePetrou=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);
   pImageTracking=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);
   int l,c,u;
   for (l=0;l<Haut;l++)
   {
      for (c=0;c<Larg;c++)
      {
         u=Larg*l+c;
         pImage->ImageRGB[u].r=0;
         pImage->ImageRGB[u].g=0;
         pImage->ImageRGB[u].b=0;
         ImageLabelEdge[u]=0;
      }
   }

   SetCursor(LoadCursor(0,IDC_WAIT));

   /* pImageAffichePoints=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);
   SetCursor(LoadCursor(0,IDC_WAIT));*/

   CMDIFrameWnd *pmd=GetMDIFrame();// Pour l'affichage du debugage

   pImage->Create(Titre+" Curvilinear Regions",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   pImage->SetActiveWindow();

   if(DebugEgde)
   {
      pImagePetrou->Create(Titre+" Edge detections",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   }

   if(DebugTracking)
   {
      pImageTracking->Create(Titre+" Edge detections",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   }

   //--------------------------------------------------------
   //			Lissage gaussien pour avoir plus de connexité
   //--------------------------------------------------------
   
   Convolution2DCedricCurvi(ImageGauss,Image,SigmaGaussien,Larg,Haut);

   //--------------------------------------------------------
   //			Init pour Détecteur contour Petrou
   //--------------------------------------------------------

   A.dim_x = Larg;
   A.dim_y = Haut;

   APasGaussien.dim_x = Larg;
   APasGaussien.dim_y = Haut;

   B.dim_x = Larg;
   B.dim_y = Haut;

   I_x.dim_x = Larg;
   I_x.dim_y = Haut;

   I_y.dim_x = Larg;
   I_y.dim_y = Haut;

   A.ima = (Pixelp**) malloc(sizeof(Pixelp*)*A.dim_y);
   APasGaussien.ima = (Pixelp**) malloc(sizeof(Pixelp*)*APasGaussien.dim_y);
   B.ima = (Pixelp**) malloc(sizeof(Pixelp*)*B.dim_y);
   I_x.ima = (Pixelp**) malloc(sizeof(Pixelp*)*I_x.dim_y);
   I_y.ima = (Pixelp**) malloc(sizeof(Pixelp*)*I_y.dim_y);

   for(i=0;i<A.dim_y;i++)
   {
      A.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*A.dim_x);
      APasGaussien.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*APasGaussien.dim_x);
      B.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*B.dim_x);
      I_x.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*I_x.dim_x);
      I_y.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*I_y.dim_x);
   }

   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         A.ima[i][j] = (Pixelp)ImageGauss[i*Larg+j];
         APasGaussien.ima[i][j] = (Pixelp)Image[i*Larg+j];
         /*pImage->ImageRGB[i*Larg+j].r = ImageGauss[i*Larg+j];
         pImage->ImageRGB[i*Larg+j].g = ImageGauss[i*Larg+j];
         pImage->ImageRGB[i*Larg+j].b = ImageGauss[i*Larg+j];*/

      }

   }




   //--------------------------------------------------------
   //			Contour Petrou
   //--------------------------------------------------------
   //AffichageProgression->UpdateData(TRUE);
   AffichageProgression->m_TPetrou = 0.0f;
   AffichageProgression->m_TCTracking = 0.0f;
   AffichageProgression->m_T1DProcess = 0.0f;
   AffichageProgression->m_T2DProcess = 0.0f;
   AffichageProgression->m_Affichage = 0.0f;
   AffichageProgression->m_TexteValProgress = "0 %";
   AffichageProgression->m_ControlBarProgressCurvi.SetPos(0);
   AffichageProgression->m_TxtNbReg="";
   AffichageProgression->UpdateData(FALSE);

   AffichageProgression->ShowWindow(bAff); // temp

   CPetrou ContourPetrou(&A, &B, &I_x,&I_y);

   QueryPerformanceCounter(&TpsStart); 
   ContourPetrou.filtre_petrou(&A,&B,6,SeuilPetrou,1);
   QueryPerformanceCounter(&TPetrou);
   tps1.QuadPart = TPetrou.QuadPart-TpsStart.QuadPart;



   AffichageProgression->m_TPetrou = (double)tps1.QuadPart/(double)FreqProc.QuadPart;
   AffichageProgression->m_TCTracking = 0.0f;
   AffichageProgression->m_T1DProcess = 0.0f;
   AffichageProgression->m_T2DProcess = 0.0f;
   AffichageProgression->m_Affichage = 0.0f;
   AffichageProgression->m_TexteValProgress = "20 %";
   AffichageProgression->m_ControlBarProgressCurvi.SetPos(20);
   AffichageProgression->UpdateData(FALSE);

   AffichageProgression->ShowWindow(bAff); //  FALSE temporaire
   //--------------------------------------------------------
   //			Recopie Image Contour Petrou pour suivi
   //--------------------------------------------------------

   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         if(DebugEgde)
         {
            pImagePetrou->ImageRGB[i*Larg+j].r = (BYTE) B.ima[i][j];
            pImagePetrou->ImageRGB[i*Larg+j].g = (BYTE) B.ima[i][j];
            pImagePetrou->ImageRGB[i*Larg+j].b = (BYTE) B.ima[i][j];
         }
         //pImageAffichePoints->ImageRGB[i*Larg+j].r = Image[i*Larg+j];//0;
         //pImageAffichePoints->ImageRGB[i*Larg+j].g = Image[i*Larg+j];//0;
         //pImageAffichePoints->ImageRGB[i*Larg+j].b = Image[i*Larg+j];//0;
         ImageFreeman[i*Larg+j] = (BYTE) B.ima[i][j];
      }

   }

   //--------------------------------------------------------
   //			suivi de contours
   //--------------------------------------------------------
   //MessageBox("avant appel suivi freeman");
   CurviDetect2D.FermetureContourSimple(ImageFreeman,info.Larg,info.Haut);

   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         if(DebugEgde)
         {
            pImagePetrou->ImageRGB[i*Larg+j].r = (BYTE) ImageFreeman[i*Larg+j];
            pImagePetrou->ImageRGB[i*Larg+j].g = (BYTE) ImageFreeman[i*Larg+j];
            pImagePetrou->ImageRGB[i*Larg+j].b = (BYTE) ImageFreeman[i*Larg+j];
         }

      }

   }

   CurviDetect2D.SuiviFreeman(ImageFreeman,info.Larg,info.Haut,pmd);

   QueryPerformanceCounter(&TTracking);

   TableauListesPointsFiltrees.RemoveAll();
   TableauListesPointsFiltreesC2D.RemoveAll();

   tps2.QuadPart = TTracking.QuadPart-TPetrou.QuadPart;

   AffichageProgression->m_TCTracking = (double)tps2.QuadPart/(double)FreqProc.QuadPart;
   AffichageProgression->m_T1DProcess = 0.0f;
   AffichageProgression->m_T2DProcess = 0.0f;
   AffichageProgression->m_Affichage = 0.0f;
   AffichageProgression->m_TexteValProgress = "40 %";
   AffichageProgression->m_ControlBarProgressCurvi.SetPos(40);
   AffichageProgression->UpdateData(FALSE);
   AffichageProgression->Invalidate(TRUE);
   AffichageProgression->ShowWindow(bAff); // temp

   //::MessageBox(NULL,"STOP","Après suivi freeman",MB_OK);
   ////---------------------------------------------------------
   ////				Affichage couleur fenetre résultats
   ////---------------------------------------------------------

   int cptlistetemp=0;

   //pImageAffichePoints->Create("Points Détecteur curvi",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   int OldIptFinalGauche=0;
   int OldIptFinalDroite=0;
   int CptExtrapolGauche,CptExtrapolDroite;



   NbListe = CurviDetect2D.TabListePoint.GetCount();

   if (NbListe !=0)
   {
      for(i=0;i<NbListe;i++)
      {
         pos = CurviDetect2D.TabListePoint[i]->GetHeadPosition();
         while(pos)
         {
            CurentPoint = CurviDetect2D.TabListePoint[i]->GetAt(pos);
            CurviDetect2D.TabListePoint[i]->GetNext(pos);
            u=Larg*CurentPoint.y+CurentPoint.x;
            ImageLabelEdge[u]=i+1;

            //if(DebugTracking)
            //{
            //	pImageTracking->ImageRGB[u].r = 255-ImageLabelEdge[u];
            //	pImageTracking->ImageRGB[u].g = 255-ImageLabelEdge[u];
            //	pImageTracking->ImageRGB[u].b = 255-ImageLabelEdge[u];
            //}
         }
      }
   }

   if (NbListe !=0)
   {
      for(i=0;i<NbListe;i++)
      {

         if(CurviDetect2D.TabListePoint[i]->GetCount()>m_LongueurMiniForme)
         {
            ListePointsFiltres = new CList<PosDetectCurvi,PosDetectCurvi>;// creation pour les listes dépassants la contraintes de largeurs mini.
            ListePointsFiltresComp = new CList<PosDetectCurvi,PosDetectCurvi>;

            double LastWidthGauche = -1.0f,LastWidthDroite = -1.0f;
            cptCoul++;
            pos = CurviDetect2D.TabListePoint[i]->GetHeadPosition();
            CptExtrapolGauche=0;
            CptExtrapolDroite=0;

            //tirage au sort de la couleur pour la forme entière

            r= (rand()*235)+20;
            g= (rand()*235)+20;
            b= (rand()*235)+20;

            CurviDetect1D.MaxiCurvi.IptFinal=-1;
            OldIptFinalGauche=-1;
            OldIptFinalDroite=-1;

            int longueur;
            CPoint PointSuiv;
            while(pos)
            {
               CurentPoint = CurviDetect2D.TabListePoint[i]->GetAt(pos);
               if (CurentPoint.x>LargeurFourier && CurentPoint.x<Larg-LargeurFourier
                  && CurentPoint.y>LargeurFourier && CurentPoint.y<Haut-LargeurFourier)
               {

                  GradMoy_x = I_x.ima[CurentPoint.y][CurentPoint.x];
                  GradMoy_y = I_y.ima[CurentPoint.y][CurentPoint.x];

                  if(GradMoy_x>0)
                     SigneDx = -Sens; 
                  else
                     SigneDx=+Sens;

                  if(GradMoy_x==0 && GradMoy_y>=0.0) SigneDx=-Sens;
                  if(GradMoy_x==0 && GradMoy_y<0.0) SigneDx=+Sens;


                  CurviDetect2D.RemplitTableauIJCoupeBiDir(TabCoupe,CurentPoint.x,CurentPoint.y,GradMoy_x,GradMoy_y,w,dw,SigneDx,&Centre,Haut,Larg);

                  //----------------------------------------------------------------------
                  //			Creation tableau gauche et droite
                  //----------------------------------------------------------------------
                  for(j=Centre-w;j<Centre+1+LargeurFourier;j++)//Gauche
                  {
                     TabGauche[Centre-j+LargeurFourier] = APasGaussien.ima[TabCoupe[j].y][TabCoupe[j].x];
                     TabCoupeGauche[Centre+LargeurFourier-j] = TabCoupe[j];
                     p=TabCoupe[j].x+Larg*TabCoupe[j].y;
                  }
                  cpti=0;
                  for(j=Centre-LargeurFourier;j<Centre+w+1;j++)//droite
                  {
                     TabDroite[cpti++] = APasGaussien.ima[TabCoupe[j].y][TabCoupe[j].x];
                  }//fin for j

                  ValIptInitFixe = LargeurFourier;

                  CurviDetect1D.RespCurviMaxiFast(TabGauche, TabCoupeGauche, LastWidthGauche,LargeurFourier,ImageLabelEdge,Larg);// modif

                  PosGauche.PosInit = TabCoupe[Centre+LargeurFourier-CurviDetect1D.MaxiCurvi.IptInit]; 
                  PosGauche.PosFinal = TabCoupe[Centre+LargeurFourier-CurviDetect1D.MaxiCurvi.IptFinal];
                  PosGauche.PositionMaxi = TabCoupe[Centre+LargeurFourier-CurviDetect1D.MaxiCurvi.Centre];
                  PosGauche.ValDist = CurviDetect1D.MaxiCurvi.ValDist;
                  PosGauche.Largeur = LastWidthGauche = CurviDetect1D.MaxiCurvi.Largeur;
                  PosGauche.PosAxe.x=((double)PosGauche.PosInit.x+(double)PosGauche.PosFinal.x)/2;
                  PosGauche.PosAxe.y=((double)PosGauche.PosInit.y+(double)PosGauche.PosFinal.y)/2;

                  if ((PosGauche.ValDist>=m_SeuilReponse) &&(PosGauche.Largeur>=m_LargeurMini))
                  {
                     ListePointsFiltres->AddTail(PosGauche);// ajout à la liste des candidats gardés
                     CptExtrapolGauche=0;
                  }
                  else
                  {
                     CptExtrapolGauche++;
                     if (CptExtrapolGauche>SeuilExtrapol)
                     {
                        longueur=ListePointsFiltres->GetCount();
                        if (longueur>m_LongueurMiniForme) 
                        { 
                           TableauListesPointsFiltrees.Add(ListePointsFiltres);
                           cptlistetemp++;
                        }
                        else delete ListePointsFiltres;
                        ListePointsFiltres = new CList<PosDetectCurvi,PosDetectCurvi>;// creation pour les listes dépassants la contraintes de largeurs mini.
                        CptExtrapolGauche=0;
                     }
                  }
                  OldIptFinalGauche=CurviDetect1D.MaxiCurvi.IptFinal;

                  CurviDetect1D.RespCurviMaxiFast(TabDroite, TabCoupe+Centre+1-LargeurFourier, LastWidthDroite,LargeurFourier,ImageLabelEdge,Larg);// modif

                  PosDroite.PosInit = TabCoupe[Centre+1-LargeurFourier+CurviDetect1D.MaxiCurvi.IptInit];
                  PosDroite.PosFinal = TabCoupe[Centre+1-LargeurFourier+CurviDetect1D.MaxiCurvi.IptFinal];
                  PosDroite.PositionMaxi = TabCoupe[Centre+1-LargeurFourier+CurviDetect1D.MaxiCurvi.Centre];
                  PosDroite.ValDist = CurviDetect1D.MaxiCurvi.ValDist;
                  PosDroite.Largeur = LastWidthDroite = CurviDetect1D.MaxiCurvi.Largeur;
                  PosDroite.PosAxe.x=((double)PosDroite.PosInit.x+(double)PosDroite.PosFinal.x)/2;
                  PosDroite.PosAxe.y=((double)PosDroite.PosInit.y+(double)PosDroite.PosFinal.y)/2;

                  if ((PosDroite.ValDist>=m_SeuilReponse) && (PosDroite.Largeur>=m_LargeurMini) )
                  {
                     ListePointsFiltresComp->AddTail(PosDroite);// ajout à la liste des candidats gardés
                     CptExtrapolDroite=0;	
                  }
                  else
                  {
                     CptExtrapolDroite++;
                     if (CptExtrapolDroite>SeuilExtrapol)
                     {
                        longueur=ListePointsFiltresComp->GetCount();
                        if (longueur>m_LongueurMiniForme)
                        {
                           TableauListesPointsFiltrees.Add(ListePointsFiltresComp);
                           cptlistetemp++;
                        }
                        else delete ListePointsFiltresComp;
                        ListePointsFiltresComp = new CList<PosDetectCurvi,PosDetectCurvi>;
                        CptExtrapolDroite=0;
                     }

                  }
               } // si on est dans l'image moins LargeurFourier
               OldIptFinalDroite=CurviDetect1D.MaxiCurvi.IptFinal;
               OldSigneDx=SigneDx;
               if(pos)CurviDetect2D.TabListePoint[PosMaxListe]->GetNext(pos);
               CptPtsGrad++;
               cptpoint++;
               /*}*/


            }// while(pos)

            // correction Jo
            longueur=ListePointsFiltres->GetCount();
            if (longueur>m_LongueurMiniForme)
            {
               TableauListesPointsFiltrees.Add(ListePointsFiltres);
               cptlistetemp++;

            } 
            else delete ListePointsFiltres;
            longueur=ListePointsFiltresComp->GetCount();
            if (longueur>m_LongueurMiniForme)
            {
               TableauListesPointsFiltrees.Add(ListePointsFiltresComp);
               cptlistetemp++;
            }
            else delete ListePointsFiltresComp;

         }//fin if nb points liste > 10
         //ajout de la liste en cours au tableau de liste de structure poscurvidetect

      }//for i
   }//fin if liste !=0 

   QueryPerformanceCounter(&T1DProcess);
   tps3.QuadPart = T1DProcess.QuadPart-TTracking.QuadPart;
   AffichageProgression->m_T1DProcess = (double)tps3.QuadPart/(double)FreqProc.QuadPart;
   AffichageProgression->m_T2DProcess = 0.0f;
   AffichageProgression->m_Affichage = 0.0f;
   AffichageProgression->m_TexteValProgress = "60 %";
   AffichageProgression->m_ControlBarProgressCurvi.SetPos(60);
   AffichageProgression->UpdateData(FALSE);
   AffichageProgression->Invalidate(TRUE);
   AffichageProgression->ShowWindow(bAff); // temp

   SectionOptimization(pImageTracking,TableauListesPointsFiltrees);
   // no filter on the length at the 2D stage
   TwoDimConstraints(TableauListesPointsFiltreesC2D,TableauListesPointsFiltrees, 10 /*m_LongueurMiniRegionFiltrage*/, SeuilDepass);

   //MessageBox("2dd termine");


   //free(TabPointTemp);
   // filtrage des carrés et cercles
   //Filter of circles and squared shapes

   PosDetectCurvi PointCourant;
   Sfpoint AxeMin, AxeMax, AxeDif;
   SBoundingBox Box;
   double LargMoy,DifMax;
   double LargMin,LargMax,EtendueLarg;
   int longueur;
   int LongueurTableauListesPointsFiltreesC2D=TableauListesPointsFiltreesC2D.GetCount();
   if (BoiteAxe)
   {
      CList<int,int> *ListeNumConnex;
      int NumConnex;

      for(IndiceReg=0;IndiceReg<LongueurTableauListesPointsFiltreesC2D;IndiceReg++)//pour toutes les régions
      {
         LargMoy=0;
         pos = TableauListesPointsFiltreesC2D[IndiceReg]->GetHeadPosition();
         longueur=TableauListesPointsFiltreesC2D[IndiceReg]->GetCount();
         if(longueur>0 /*m_LongueurMiniRegionFiltrage*/)// si la liste n'est pas vide
         {
            // dermination liste englobante
            PointCourant = TableauListesPointsFiltreesC2D[IndiceReg]->GetAt(pos); // on récupère ce qui sera n-1
            TableauListesPointsFiltreesC2D[IndiceReg]->GetNext(pos);// on pointe sur le suivant
            LargMoy+=PointCourant.Largeur;
            AxeMin.x=AxeMax.x=PointCourant.PosAxe.x;
            AxeMin.y=AxeMax.y=PointCourant.PosAxe.y;
            LargMin=LargMax=PointCourant.Largeur;
            for (k=1;k<longueur;k++)
            {
               PointCourant = TableauListesPointsFiltreesC2D[IndiceReg]->GetAt(pos); // on récupère ce qui sera n-1
               if (PointCourant.PosAxe.x<AxeMin.x) AxeMin.x=PointCourant.PosAxe.x;
               if (PointCourant.PosAxe.y<AxeMin.y) AxeMin.y=PointCourant.PosAxe.y;
               if (PointCourant.PosAxe.x>AxeMax.x) AxeMax.x=PointCourant.PosAxe.x;
               if (PointCourant.PosAxe.y>AxeMax.y) AxeMax.y=PointCourant.PosAxe.y;
               LargMoy+=PointCourant.Largeur;
               if (PointCourant.Largeur>LargMax) LargMax=PointCourant.Largeur;
               if (PointCourant.Largeur<LargMin) LargMin=PointCourant.Largeur;
               TableauListesPointsFiltreesC2D[IndiceReg]->GetNext(pos);// on pointe sur le suivant
            }
            LargMoy/=longueur;
            AxeDif.x=fabs(AxeMin.x-AxeMax.x);
            AxeDif.y=fabs(AxeMin.y-AxeMax.y);
            DifMax=__max(AxeDif.x,AxeDif.y);
            EtendueLarg=fabs(LargMax-LargMin);
            if (DifMax<LargMoy*1.5/* || EtendueLarg>LargMoy*/) 
            {
               TableauListesPointsFiltreesC2D.RemoveAt(IndiceReg);
               IndiceReg--;
               LongueurTableauListesPointsFiltreesC2D--;
            }
            else
            {
               Box.AxeMin=AxeMin;
               Box.AxeMax=AxeMax;
               TableauListesBoundingBox.Add(Box);
               ListeNumConnex =new CList<int,int>;
               TabNumListesConnex.Add(ListeNumConnex);
            }
         } // si longueur ok

      } // for indice
   } // si boite Axe

   LongueurTableauListesPointsFiltreesC2D=TableauListesPointsFiltreesC2D.GetCount();
   int NumConnex;
   SBoundingBox BoxA,BoxB;
   int longlis=TabNumListesConnex.GetCount();

   for(IndiceReg=0;IndiceReg<LongueurTableauListesPointsFiltreesC2D;IndiceReg++)//pour toutes les régions
   {
      BoxA=TableauListesBoundingBox[IndiceReg];
      for( int IndiceRegB=0;IndiceRegB<LongueurTableauListesPointsFiltreesC2D;IndiceRegB++)//pour toutes les régions
      {
         if (IndiceReg!=IndiceRegB) 
         {
            BoxB=TableauListesBoundingBox[IndiceRegB];
            int Ext=0;
            if (BoxA.AxeMin.x>BoxB.AxeMax.x || BoxA.AxeMax.x<BoxB.AxeMin.x) Ext=1;
            if (BoxA.AxeMin.y>BoxB.AxeMax.y || BoxA.AxeMax.y<BoxB.AxeMin.y) Ext=1;
            if (Ext==0) // there is an intersection between boxes
            {
               // add the num of the list to the liste of candidates
               TabNumListesConnex[IndiceReg]->AddHead(IndiceRegB);
            }
         }
      }
   }

   LongueurTableauListesPointsFiltreesC2D=TableauListesPointsFiltreesC2D.GetCount();
   int *TabIndice;
   int *TabListeOk;
   TabIndice=(int *) malloc (LongueurTableauListesPointsFiltreesC2D*sizeof(int));
   TabListeOk=(int *) malloc (LongueurTableauListesPointsFiltreesC2D*sizeof(int));


   for(IndiceReg=0;IndiceReg<LongueurTableauListesPointsFiltreesC2D;IndiceReg++)//pour toutes les régions
   {
      TabListeOk[IndiceReg]=1; // mark all the regions as OK
      TabIndice[IndiceReg]=IndiceReg;
   }

   // fast merge
   TryToMerge(TableauListesPointsFiltreesC2D,TabNumListesConnex,TabListeOk,pImage);

   //MergeConnectedRegions(TableauListesPointsFiltreesC2D,pImage);

   QueryPerformanceCounter(&T2DProcess);
   tps4.QuadPart = T2DProcess.QuadPart-T1DProcess.QuadPart;
   AffichageProgression->m_T2DProcess = (double)tps4.QuadPart/(double)FreqProc.QuadPart;
   AffichageProgression->m_Affichage = 0.0f;
   AffichageProgression->m_TexteValProgress = "80 %";
   AffichageProgression->m_ControlBarProgressCurvi.SetPos(80);
   AffichageProgression->UpdateData(FALSE);
   AffichageProgression->Invalidate(TRUE);
   AffichageProgression->ShowWindow(bAff); // temp

   //BOOL OuvertureFichierOk = FichierExportMatas.Open("ExportCurvipoints.txt",CFile::modeCreate|CFile::modeWrite);


   // remove too short regions

   for(IndiceReg=0;IndiceReg<LongueurTableauListesPointsFiltreesC2D;IndiceReg++)//pour toutes les régions
   {
      longueur=TableauListesPointsFiltreesC2D[IndiceReg]->GetCount();
      if(longueur<m_LongueurMiniRegionFiltrage && TabListeOk[IndiceReg]==1)// si la liste n'est pas vide
      {	
         TabListeOk[IndiceReg]=-4;
      }
   }

   //-----------------------------------------------------------------------------------
   // remove multiple regions
   //-----------------------------------------------------------------------------------
#define REMOVEMULTIPLENB
#ifdef REMOVEMULTIPLENB

   //int LongueurTableauListesPointsFiltreesC2D;
   int cptaff = 1;

   int Intersect=0,Surf1=0,Surf2=0,SurfMin;
   double Recouv;
   double SeuilRecou=0.7;
   int NbReg1=1;
   //CIma Test;
   //Test.LitBMP("f:\\johel\\prog\\bmp\\temp.bmp");
   IncrementeCouleur(&r,&g,&b,256);
   for (u=0;u<Taille;u++) {Plan1[u]=0;}	

   LongueurTableauListesPointsFiltreesC2D=TableauListesPointsFiltreesC2D.GetCount();
   // find the first non merged region
   int IndiceRegInit=0;
   while(TabListeOk[IndiceRegInit]!=1 && IndiceRegInit<LongueurTableauListesPointsFiltreesC2D)
   {
      IndiceRegInit++;
   } 
   if (IndiceRegInit!=LongueurTableauListesPointsFiltreesC2D)
   {
      BuildImageFromList(Plan1,TableauListesPointsFiltreesC2D[IndiceRegInit]);

      //for (u=0;u<Taille;u++)
      //{
      //if (Plan1[u]!=0) pImage->ImageRGB[u].r=255; else pImage->ImageRGB[u].r=0;
      //if (Plan1[u]!=0) pImage->ImageRGB[u].g=255;else pImage->ImageRGB[u].g=0;
      //if (Plan1[u]!=0) pImage->ImageRGB[u].b=255;else pImage->ImageRGB[u].b=0;
      //}
      //pImage->SendMessage(WM_PAINT);
      //pImage->Invalidate(TRUE);
      //Sleep(5000);
      for(IndiceReg=IndiceRegInit+1;IndiceReg<LongueurTableauListesPointsFiltreesC2D;IndiceReg++)//pour toutes les régions
      {
         longueur=TableauListesPointsFiltreesC2D[IndiceReg]->GetCount();
         if(TabListeOk[IndiceReg]==1)// si la liste n'est pas vide
         {	
            BuildImageFromList(Plan2,TableauListesPointsFiltreesC2D[IndiceReg]);
            /*for (u=0;u<Taille;u++)
            {
            if (Plan2[u]!=0) pImage->ImageRGB[u].r=255; else pImage->ImageRGB[u].r=0;
            if (Plan2[u]!=0) pImage->ImageRGB[u].g=255;else pImage->ImageRGB[u].g=0;
            if (Plan2[u]!=0) pImage->ImageRGB[u].b=255;else pImage->ImageRGB[u].b=0;
            }
            pImage->SendMessage(WM_PAINT);
            pImage->Invalidate(TRUE);
            Sleep(5000);*/
            // tester recouvrement p1 et p2, et recopier dans p1
            int IndiceMax=0;
            int Sens;
            //		Sens=CurviDetect2D.MesureRecouv(Plan1,Plan2,Larg,Haut,NbReg1,SeuilRecou,&IndiceMax); // bug 17 7
            Sens=CurviDetect2D.MesureRecouv(Plan1,Plan2,Larg,Haut,IndiceReg+1/*NbReg1*/,SeuilRecou,&IndiceMax,TabListeOk);
            //Sens=0;
            if (Sens==12) // on garde la 1 et on vire la 2
            {
               TabListeOk[IndiceReg]=-1; // mark the region as removed
               //MessageBox("region 2 removed");
            }
            if (Sens==21) //recopie région du plan2 dans plan1// if (Surf2>Surf1)
            {

               //TabListeOk[TabIndice[IndiceMax-1]]=-3; 
               TabListeOk[IndiceMax-1]=-3; // mark the region as removed
               for (u=0;u<Taille;u++)
               {
                  if (Plan1[u]==IndiceMax) Plan1[u]=0;
                  if (Plan2[u]!=0 ) Plan1[u]=IndiceMax;
               }
               //MessageBox("region 1 removed");
               //TabIndice[IndiceMax-1]=IndiceReg;
            } // si recouv

            if (Sens==0) 
            {
               NbReg1++;
               for (u=0;u<Taille;u++)
               {
                  //if (Plan2[u]!=0 && Plan1[u]==0) Plan1[u]=NbReg1;//Plan2[u]; bug 17 7
                  if (Plan2[u]!=0 && Plan1[u]==0) Plan1[u]=IndiceReg+1;//NbReg1;//Plan2[u];

               }
            } // si peu de recou

         }// si les listes ne sont pas vide (fin if)


      }//fin for IndiceReg

   } // if more than one region
#endif 

   // visualization of final results
   int NbRegFinal=0;
   int LngTableauListesPointsFiltrees2D=TableauListesPointsFiltreesC2D.GetCount();
   for(IndiceReg=0;IndiceReg<LngTableauListesPointsFiltrees2D;IndiceReg++)//pour toutes les régions
   {
      if (TabListeOk[IndiceReg]==1)
      {
         NbRegFinal++;
         BuildImageFromList(Plan1,TableauListesPointsFiltreesC2D[IndiceReg]);
         TireCouleurAlea(&r,&g,&b);
         for (u=0;u<Taille;u++)
         {
            if (Plan1[u]==1) pImage->ImageRGB[u].r=r;
            if (Plan1[u]==1) pImage->ImageRGB[u].g=g;
            if (Plan1[u]==1) pImage->ImageRGB[u].b=b;
         }
      }
   }


   // // affichage des régions non filtrées  dans image tracking
   if(DebugTracking)
   {
      int LngTableauListesPointsFiltrees=TableauListesPointsFiltrees.GetCount();
      for(IndiceReg=0;IndiceReg<LngTableauListesPointsFiltrees;IndiceReg++)//pour toutes les régions
      {
         BuildImageFromList(Plan1,TableauListesPointsFiltrees[IndiceReg]);
         TireCouleurAlea(&r,&g,&b);
         for (u=0;u<Taille;u++)
         {
            if (Plan1[u]==1) pImageTracking->ImageRGB[u].r=r;
            if (Plan1[u]==1) pImageTracking->ImageRGB[u].g=g;
            if (Plan1[u]==1) pImageTracking->ImageRGB[u].b=b;
         }
      }
   }


   free(TabIndice);
   free(TabListeOk);	

   QueryPerformanceCounter(&TAffichage);


   tps5.QuadPart = TAffichage.QuadPart-T2DProcess.QuadPart;

   AffichageProgression->m_Affichage = (double)tps5.QuadPart/(double)FreqProc.QuadPart;
   AffichageProgression->m_TexteValProgress = "100 %";
   AffichageProgression->m_ControlBarProgressCurvi.SetPos(100);
   AffichageProgression->m_TxtNbReg.Format("Nb Reg=%d",NbRegFinal);
   AffichageProgression->UpdateData(FALSE);

   //AffichageProgression->Invalidate(TRUE);


   //MessageBox("Export & affichage ok ","Export & affichage ok ",MB_OK);
   //----------------------------------------------------------------------------------------------------
   Invalidate(FALSE);

   SetCursor(LoadCursor(0,IDC_ARROW));

   for(i=0;i<A.dim_y;i++)
   {
      free(A.ima[i]);
      free(APasGaussien.ima[i]);
      free(B.ima[i]);
      free(I_x.ima[i]);
      free(I_y.ima[i]);
   }

   free(A.ima);
   free(APasGaussien.ima);
   free(B.ima);
   free(I_x.ima);
   free(I_y.ima);
   free(TabCoupe);


   free(TabDroite);
   free(TabGauche);
   free(ImageFreeman);
   free(ImageGauss);

   free(ImageLabelEdge);

   delete TabCoupeGauche;

   for(IndiceReg=0;IndiceReg<TableauListesPointsFiltreesC2D.GetCount();IndiceReg++)//pour toutes les régions
   {
      TabNumListesConnex[IndiceReg]->RemoveAll();
      delete TabNumListesConnex[IndiceReg];
      TableauListesPointsFiltreesC2D[IndiceReg]->RemoveAll();
      delete TableauListesPointsFiltreesC2D[IndiceReg];
   }
   for(IndiceReg=0;IndiceReg<TableauListesPointsFiltrees.GetCount();IndiceReg++)//pour toutes les régions
   {
      TableauListesPointsFiltrees[IndiceReg]->RemoveAll();
      delete TableauListesPointsFiltrees[IndiceReg];
   }



#endif

   VisualiseCurvi(DebugEgde, DebugTracking, Larg, Haut, edge, TableauListesPointsFiltrees, TableauListesPointsFiltreesC2D);
   //MessageBox("Export & affichage ok ","Export & affichage ok ",MB_OK);
   //----------------------------------------------------------------------------------------------------
   Invalidate(FALSE);
   SetCursor(LoadCursor(0,IDC_ARROW));

   delete im; delete edge;
}

void CImageWnd::VisualiseCurvi(bool DebugEdge, bool DebugTracking, int Larg, int Haut, utls::BAry *edge, TableOfCurvilinearRegions &TableauListesPointsFiltrees, TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D)
{
   CImageWnd * pImage=NULL;	
   CImageWnd * pImagePetrou=NULL;
   CImageWnd * pImageTracking=NULL;

   int Taille = Larg*Haut;
   pImage=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);
   pImagePetrou=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);
   pImageTracking=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);

   pImage->Create(Titre+" Curvilinear Regions",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   pImage->SetActiveWindow();

   if(DebugEdge)
   {
      unsigned char *ImageFreeman = edge->first();
      pImagePetrou->Create(Titre+" Edge detections",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
      for (int u=0;u<Taille;u++)
      {
         pImagePetrou->ImageRGB[u].r = ImageFreeman[u];
         pImagePetrou->ImageRGB[u].g = ImageFreeman[u];
         pImagePetrou->ImageRGB[u].b = ImageFreeman[u];
      }
   }

   if(DebugTracking)
   {
      pImageTracking->Create(Titre+" Edge detections",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   }

   BYTE r,g,b;

   utls::IAry *p1 = new utls::IAry(Haut, Larg);
   int *Plan1 = p1->first();  p1->clear();

   // visualization of final results
   int NbRegFinal=0;
   int LngTableauListesPointsFiltrees2D=TableauListesPointsFiltreesC2D.size();
   for(int IndiceReg=0;IndiceReg<LngTableauListesPointsFiltrees2D;IndiceReg++)//pour toutes les régions
   {
      IPoint2D bbmin, bbmax;
      BuildImageFromList(p1, TableauListesPointsFiltreesC2D[IndiceReg], bbmin, bbmax);
      TireCouleurAlea(&r,&g,&b);
      for (int u=0;u<Taille;u++)
      {
         if (Plan1[u]==1) pImage->ImageRGB[u].r=r;
         if (Plan1[u]==1) pImage->ImageRGB[u].g=g;
         if (Plan1[u]==1) pImage->ImageRGB[u].b=b;
      }
      p1->clear(bbmin.y, bbmax.y, bbmin.x, bbmax.x);
   }


   // // affichage des régions non filtrées  dans image tracking
   if(DebugTracking)
   {
      int LngTableauListesPointsFiltrees=TableauListesPointsFiltrees.size();
      for(int IndiceReg=0;IndiceReg<LngTableauListesPointsFiltrees;IndiceReg++)//pour toutes les régions
      {
         IPoint2D bbmin, bbmax;
         BuildImageFromList(p1, TableauListesPointsFiltrees[IndiceReg], bbmin, bbmax);
         TireCouleurAlea(&r,&g,&b);
         for (int u=0;u<Taille;u++)
         {
            if (Plan1[u]==1) pImageTracking->ImageRGB[u].r=r;
            if (Plan1[u]==1) pImageTracking->ImageRGB[u].g=g;
            if (Plan1[u]==1) pImageTracking->ImageRGB[u].b=b;
         }
         p1->clear(bbmin.y, bbmax.y, bbmin.x, bbmax.x);
      }
   }
   free(Plan1);
}

void CImageWnd::OnOutilsParamcurvidetect()
{
   // TODO : ajoutez ici le code de votre gestionnaire de commande

   BoiteReglageParam.DoModal();
}

void CImageWnd::OnMatchingSelection3points()
{
   MATCHINGCURVI = TRUE;
   CptManipMatching = 0;
   FichierMatching.Open("FichierMaching.txt",CFile::modeCreate|CFile::modeWrite|CFile::modeNoTruncate);
   FichierMatching.SeekToEnd();

}

void CImageWnd::OnTransformationimageParamtransformation()
{

   CFile FichierSav;

   if (BoiteParamAffine.DoModal()==IDOK)
   {
      ParamAffine[0] = BoiteParamAffine.m_M1;
      ParamAffine[1] = BoiteParamAffine.m_M2;
      ParamAffine[2] = BoiteParamAffine.m_M3;
      ParamAffine[3] = BoiteParamAffine.m_M4;
      ParamAffine[4] = BoiteParamAffine.m_M5;
      ParamAffine[5] = BoiteParamAffine.m_M6;
      ParamAffine[6] = BoiteParamAffine.m_M7;
      if(FichierSav.Open("AffineParam.aff",CFile::modeCreate|CFile::modeWrite))
      {
         FichierSav.Write(&ParamAffine,sizeof(double)*7);
         FichierSav.Close();
      }
   }


}

void CImageWnd::OnTransformationimageRun()
{
   CFile FichierSav;
   CStdioFile FichierPoints;
   CString ChaineRelecture;
   if(FichierSav.Open("AffineParam.aff",CFile::modeRead))
   {
      FichierSav.Read(&ParamAffine,sizeof(double)*7);
      FichierSav.Close();
   }

   //ouverture du fichier comptenant les 3 pts




   double x1,x2,x3,y1,y2,y3;
   double xp1,xp2,xp3,yp1,yp2,yp3;
   int i=0,imodif=0,j=0,jmodif=0;
   int MaxLargeur=0,MinLargeur=0,MaxHauteur=0,MinHauteur=0;
   int LargeurNewImage, HauteurNewImage;

   CPoint * LUTPourTransAff = (CPoint *)malloc(sizeof(CPoint)*Haut*Larg);



   if(FichierPoints.Open("FichierMaching.txt",CFile::modeRead))
   {
      FichierPoints.ReadString(ChaineRelecture);
      x1 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      y1 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      x2 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      y2 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      x3 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      y3 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.Close();
   }





   xp1 = ParamAffine[0]*(double)x1 + ParamAffine[1]*(double)y1 +ParamAffine[2];
   yp1 = ParamAffine[3]*(double)x1 + ParamAffine[4]*(double)y1 +ParamAffine[5];

   xp2 = ParamAffine[0]*(double)x2 + ParamAffine[1]*(double)y2 +ParamAffine[2];
   yp2 = ParamAffine[3]*(double)x2 + ParamAffine[4]*(double)y2 +ParamAffine[5];

   xp3 = ParamAffine[0]*(double)x3 + ParamAffine[1]*(double)y3 +ParamAffine[2];
   yp3 = ParamAffine[3]*(double)x3 + ParamAffine[4]*(double)y3 +ParamAffine[5];


   //init LUT à 0
   for(i=0;i<Larg;i++)
      for(j=0;j<Haut;j++)
      {
         LUTPourTransAff[j*Larg+i].x = 0;
         LUTPourTransAff[j*Larg+i].y = 0;
      }
      //calcul LUT
      for(i=0;i<Larg;i++)
         for(j=0;j<Haut;j++)
         {
            imodif = ParamAffine[0]*(double)i + ParamAffine[1]*(double)j +ParamAffine[2];
            jmodif = ParamAffine[3]*(double)i + ParamAffine[4]*(double)j +ParamAffine[5];

            LUTPourTransAff[j*Larg+i].x = imodif;
            LUTPourTransAff[j*Larg+i].y = jmodif;
         }
         ////recherche taille max
         //for(i=0;i<Larg;i++)
         //	for(j=0;j<Haut;j++)
         //	{
         //		if(LUTPourTransAff[j*Larg+i].x > MaxLargeur )
         //			MaxLargeur = LUTPourTransAff[j*Larg+i].x;

         //		if(LUTPourTransAff[j*Larg+i].x < MinLargeur )
         //			MinLargeur = LUTPourTransAff[j*Larg+i].x;

         //		if(LUTPourTransAff[j*Larg+i].y > MaxHauteur)
         //			MaxHauteur= LUTPourTransAff[j*Larg+i].y;

         //		if(LUTPourTransAff[j*Larg+i].y < MinHauteur)
         //			MinHauteur= LUTPourTransAff[j*Larg+i].y;
         //	}

         //	LargeurNewImage = MaxLargeur - MinLargeur;
         //	HauteurNewImage = MaxHauteur - MinHauteur;

         ////Normalisation entre 0 et MaxLargeur/MaxHauteur
         //	for(i=0;i<Larg;i++)
         //		for(j=0;j<Haut;j++)
         //		{
         //			LUTPourTransAff[j*Larg+i].x = LUTPourTransAff[j*Larg+i].x - MinLargeur;
         //			LUTPourTransAff[j*Larg+i].y = LUTPourTransAff[j*Larg+i].y - MinHauteur;
         //		}



         CImageWnd * pImage=NULL;	
         pImage=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);

         SetCursor(LoadCursor(0,IDC_WAIT));

         CMDIFrameWnd *pmd=GetMDIFrame();
         pImage->Create(Titre+" Curvilinear color",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
         for(i=0;i<Larg;i++)
            for(j=0;j<Haut;j++)
            {

               if((LUTPourTransAff[j*Larg+i].y<Haut)&&(LUTPourTransAff[j*Larg+i].y>0)&&(LUTPourTransAff[j*Larg+i].x<Larg)&&(LUTPourTransAff[j*Larg+i].x>0))
               {
                  pImage->ImageRGB[LUTPourTransAff[j*Larg+i].y *Larg+LUTPourTransAff[j*Larg+i].x].r = Image[j*Larg+i];
                  pImage->ImageRGB[LUTPourTransAff[j*Larg+i].y *Larg+LUTPourTransAff[j*Larg+i].x].g = 0;//Image[j*Larg+i];
                  pImage->ImageRGB[LUTPourTransAff[j*Larg+i].y *Larg+LUTPourTransAff[j*Larg+i].x].b = 0;//Image[j*Larg+i];
               }

            }

            /*if ((xp1<Larg)&& (yp1<Haut))
            {
            pImage->ImageRGB[yp1*Larg+xp1].r = 255;
            pImage->ImageRGB[yp1*Larg+xp1].g = 0;
            pImage->ImageRGB[yp1*Larg+xp1].b = 0;
            }
            if ((xp2<Larg)&& (yp2<Haut))
            {
            pImage->ImageRGB[yp2*Larg+xp2].r = 0;
            pImage->ImageRGB[yp2*Larg+xp2].g = 255;
            pImage->ImageRGB[yp2*Larg+xp2].b = 0;
            }

            if ((xp3<Larg)&& (yp3<Haut))
            {
            pImage->ImageRGB[yp2*Larg+xp2].r = 0;
            pImage->ImageRGB[yp2*Larg+xp2].g = 0;
            pImage->ImageRGB[yp2*Larg+xp2].b = 255;
            }*/


            Invalidate(FALSE);
            SetCursor(LoadCursor(0,IDC_ARROW));


}

void CImageWnd::OnMatchingCalculmatricepassage()
{

   CStdioFile FichierPoints;
   CString ChaineRelecture;


   //ouverture du fichier comptenant les 3 pts




   double x1,x2,x3,x4,y1,y2,y3,y4;
   double xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4;

   int i=0,imodif=0,j=0,jmodif=0;
   int MaxLargeur=0,MinLargeur=0,MaxHauteur=0,MinHauteur=0;
   int LargeurNewImage, HauteurNewImage;

   double m1,m2,m3,m4,m5,m6,Denominateur;

   //relecture 6 pts

   if(FichierPoints.Open("FichierMaching.txt",CFile::modeRead))
   {
      FichierPoints.ReadString(ChaineRelecture);
      x1 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      y1 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      x2 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      y2 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      x3 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      y3 = atof((LPCSTR)ChaineRelecture);
      /*FichierPoints.ReadString(ChaineRelecture);
      x4 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      y4 = atof((LPCSTR)ChaineRelecture);*/
      FichierPoints.ReadString(ChaineRelecture);			
      xp1 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      yp1 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      xp2 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      yp2 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      xp3 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      yp3 = atof((LPCSTR)ChaineRelecture);
      /*FichierPoints.ReadString(ChaineRelecture);
      xp4 = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      yp4 = atof((LPCSTR)ChaineRelecture);*/
      FichierPoints.Close();
   }

   //Denominateur = (x1*y2)-(x3*y2)+(x3*y1)+(x2*y3)-(x2*y1)-(x1*y3);


   m1 = -(y3*xp1-y2*xp1+y2*xp3+xp2*y1-xp2*y3-xp3*y1)/(x3*y1-x3*y2-x2*y1-x1*y3+x1*y2+x2*y3);
   m2 = (-x2*xp1-x3*xp2+x3*xp1-x1*xp3+x1*xp2+x2*xp3)/(x3*y1-x3*y2-x2*y1-x1*y3+x1*y2+x2*y3);
   m3 = (y3*x2*xp1-xp3*x2*y1-y2*x3*xp1+y2*x1*xp3+xp2*x3*y1-y3*x1*xp2)/(x3*y1-x3*y2-x2*y1-x1*y3+x1*y2+x2*y3);


   m4 = -(y3*yp1-y2*yp1+y2*yp3+yp2*y1-yp2*y3-yp3*y1)/(x3*y1-x3*y2-x2*y1-x1*y3+x1*y2+x2*y3);
   m5 = (-x2*yp1-x1*yp3+x3*yp1+x1*yp2+x2*yp3-x3*yp2)/(x3*y1-x3*y2-x2*y1-x1*y3+x1*y2+x2*y3);
   m6 = (y3*x2*yp1-yp3*x2*y1+y2*x1*yp3-y2*x3*yp1+yp2*x3*y1-y3*x1*yp2)/(x3*y1-x3*y2-x2*y1-x1*y3+x1*y2+x2*y3);


   //m2 = (-x2*xp1+x1*xp2+x4*xp1-x4*xp2+x2*xp4-x1*xp4)/(-x2*y1+x2*y4+x4*y1+x1*y2-x1*y4-x4*y2);
   //m6 = (-x1*y4*yp3+x1*y3*yp4+x4*y1*yp3-x3*y1*yp4+y4*x3*yp1-y3*x4*yp1)/(-x1*y4+x3*y4-x3*y1-x4*y3+x4*y1+x1*y3);
   //m4 = -(-y4*yp3+yp3*y1-y1*yp4+yp1*y4-y3*yp1+y3*yp4)/(-x1*y4+x3*y4-x3*y1-x4*y3+x4*y1+x1*y3);
   ////m7 = 0;
   ////m8 = 0;
   //m5 = (x1*yp3+x4*yp1-x4*yp3-x3*yp1-x1*yp4+x3*yp4)/(-x1*y4+x3*y4-x3*y1-x4*y3+x4*y1+x1*y3);
   ////m9 = 1;
   //m3 = (-x1*y4*xp2+x1*y2*xp4+x4*y1*xp2-x2*y1*xp4+y4*x2*xp1-y2*x4*xp1)/(-x2*y1+x2*y4+x4*y1+x1*y2-x1*y4-x4*y2);
   //m1 = -(-y4*xp2+xp2*y1-y1*xp4-y2*xp1+xp1*y4+y2*xp4)/(-x2*y1+x2*y4+x4*y1+x1*y2-x1*y4-x4*y2);

   ParamAffine[0] = m1;
   ParamAffine[1] = m2;
   ParamAffine[2] = m3;
   ParamAffine[3] = m4;
   ParamAffine[4] = m5;
   ParamAffine[5] = m6;
   ParamAffine[6] = 1;


   CFile FichierSav;

   if(FichierSav.Open("AffineParam.aff",CFile::modeCreate|CFile::modeWrite))
   {
      FichierSav.Write(&ParamAffine,sizeof(double)*7);
      FichierSav.Close();
   }
}

void CImageWnd::OnMatchingCalculmatricepassage504()
{
   CStdioFile FichierPoints;
   CString ChaineRelecture;


   //ouverture du fichier comptenant les 3 pts




   double xa,xb,xc,xd,ya,yb,yc,yd;
   double xpa,xpb,xpc,xpd,ypa,ypb,ypc,ypd;
   double k1,k2,k3,kp1,kp2,kp3;
   double R1,R2,R3;

   int i=0,imodif=0,j=0,jmodif=0;
   int MaxLargeur=0,MinLargeur=0,MaxHauteur=0,MinHauteur=0;
   int LargeurNewImage, HauteurNewImage;

   double m1,m2,m3,m4,m5,m6,m7,m8,m9;

   //relecture 8 pts

   if(FichierPoints.Open("FichierMaching4.txt",CFile::modeRead))
   {
      FichierPoints.ReadString(ChaineRelecture);
      xa = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      ya = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      xb = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      yb = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      xc = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      yc = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      xd = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      yd = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);			
      xpa = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      ypa = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      xpb = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      ypb = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      xpc = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      ypc = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      xpd = atof((LPCSTR)ChaineRelecture);
      FichierPoints.ReadString(ChaineRelecture);
      ypd = atof((LPCSTR)ChaineRelecture);
      FichierPoints.Close();
   }



   k3 = (xb*ya-yb*xa-xb*yd+xa*yd-ya*xd+yb*xd)/(xb*ya-xb*yc+xa*yc-ya*xc+yb*xc-yb*xa);
   k2 = -(xa*yd-xa*yc+ya*xc-xc*yd+xd*yc-ya*xd)/(xb*ya-xb*yc+xa*yc-ya*xc+yb*xc-yb*xa); 
   k1 = -(xc*yd-xd*yc-xb*yd+yb*xd+xb*yc-yb*xc)/(xb*ya-xb*yc+xa*yc-ya*xc+yb*xc-yb*xa);



   kp2 = -(-xpa*ypd+xpa*ypc+xpc*ypd-ypa*xpc-xpd*ypc+ypa*xpd)/(-xpb*ypa+xpb*ypc-xpa*ypc+ypa*xpc+ypb*xpa-ypb*xpc);
   kp1 = (xpb*ypc-ypb*xpc+xpc*ypd-xpd*ypc-xpb*ypd+ypb*xpd)/(-xpb*ypa+xpb*ypc-xpa*ypc+ypa*xpc+ypb*xpa-ypb*xpc);
   kp3 = (-xpb*ypa-xpa*ypd+ypa*xpd+xpb*ypd+ypb*xpa-ypb*xpd)/(-xpb*ypa+xpb*ypc-xpa*ypc+ypa*xpc+ypb*xpa-ypb*xpc);

   R1 = kp1/k1;
   R2 = kp2/k2;
   R3 = kp3/k3;


   m2 = -(xb*R3*xpc+xa*R2*xpb-xa*R3*xpc+xc*R1*xpa-xc*R2*xpb-xb*R1*xpa)/(xb*ya-xb*yc+xa*yc-ya*xc+yb*xc-yb*xa);
   m1 = (ya*R2*xpb-ya*R3*xpc+yb*R3*xpc-R2*xpb*yc+R1*xpa*yc-R1*xpa*yb)/(xb*ya-xb*yc+xa*yc-ya*xc+yb*xc-yb*xa);
   m3 = (ya*xb*R3*xpc-R1*xpa*xb*yc-yb*xa*R3*xpc+yb*xc*R1*xpa+R2*xpb*xa*yc-ya*xc*R2*xpb)/(xb*ya-xb*yc+xa*yc-ya*xc+yb*xc-yb*xa);
   m5 = -(-xa*R3*ypc+xc*R1*ypa+xb*R3*ypc-xb*R1*ypa-xc*R2*ypb+xa*R2*ypb)/(xb*ya-xb*yc+xa*yc-ya*xc+yb*xc-yb*xa);
   m4 = (yc*R1*ypa+yb*R3*ypc-yb*R1*ypa-R3*ypc*ya+R2*ypb*ya-R2*ypb*yc)/(xb*ya-xb*yc+xa*yc-ya*xc+yb*xc-yb*xa);
   m6 = (yb*xc*R1*ypa-R2*ypb*ya*xc-yc*xb*R1*ypa+yc*xa*R2*ypb+R3*ypc*xb*ya-yb*xa*R3*ypc)/(xb*ya-xb*yc+xa*yc-ya*xc+yb*xc-yb*xa);
   m7 = (R1*yc-R1*yb+yb*R3+ya*R2-R2*yc-ya*R3)/(xb*ya-xb*yc+xa*yc-ya*xc+yb*xc-yb*xa);
   m8 = -(-xb*R1-xc*R2+xc*R1-xa*R3+xa*R2+xb*R3)/(xb*ya-xb*yc+xa*yc-ya*xc+yb*xc-yb*xa);
   m9 = (-R1*xb*yc+ya*xb*R3+yb*xc*R1-yb*xa*R3+R2*xa*yc-ya*xc*R2)/(xb*ya-xb*yc+xa*yc-ya*xc+yb*xc-yb*xa);

   ParamHomogra[0] = m1;
   ParamHomogra[1] = m2;
   ParamHomogra[2] = m3;
   ParamHomogra[3] = m4;
   ParamHomogra[4] = m5;
   ParamHomogra[5] = m6;
   ParamHomogra[6] = m7;
   ParamHomogra[7] = m8;
   ParamHomogra[8] = m9;

   CFile FichierSav;

   if(FichierSav.Open("HomographParam.aff",CFile::modeCreate|CFile::modeWrite))
   {
      FichierSav.Write(&ParamHomogra,sizeof(double)*9);
      FichierSav.Close();
   }
}

void CImageWnd::OnMatchingSelection4points()
{
   MATCHINGCURVI4 = TRUE;
   CptManipMatching = 0;
   FichierMatching4.Open("FichierMaching4.txt",CFile::modeCreate|CFile::modeWrite|CFile::modeNoTruncate);
   FichierMatching4.SeekToEnd();
}

void CImageWnd::OnTransformationimageRun506()
{
   CFile FichierSav;
   CStdioFile FichierPoints;
   CString ChaineRelecture;
   if(FichierSav.Open("HomographParam.aff",CFile::modeRead))
   {
      FichierSav.Read(&ParamHomogra,sizeof(double)*9);
      FichierSav.Close();
   }

   //ouverture du fichier comptenant les 3 pts




   double x1,x2,x3,y1,y2,y3;
   double xp1,xp2,xp3,yp1,yp2,yp3;
   int i=0,j=0,z=0;
   double imodif=0,jmodif=0,zmodif=0;
   int MaxLargeur=0,MinLargeur=0,MaxHauteur=0,MinHauteur=0;
   int LargeurNewImage, HauteurNewImage;

   CPoint * LUTPourTransAff = (CPoint *)malloc(sizeof(CPoint)*Haut*Larg);



   ///*if(FichierPoints.Open("FichierMaching.txt",CFile::modeRead))
   //{
   //	FichierPoints.ReadString(ChaineRelecture);
   //	x1 = atof((LPCSTR)ChaineRelecture);
   //	FichierPoints.ReadString(ChaineRelecture);
   //	y1 = atof((LPCSTR)ChaineRelecture);
   //	FichierPoints.ReadString(ChaineRelecture);
   //	x2 = atof((LPCSTR)ChaineRelecture);
   //	FichierPoints.ReadString(ChaineRelecture);
   //	y2 = atof((LPCSTR)ChaineRelecture);
   //	FichierPoints.ReadString(ChaineRelecture);
   //	x3 = atof((LPCSTR)ChaineRelecture);
   //	FichierPoints.ReadString(ChaineRelecture);
   //	y3 = atof((LPCSTR)ChaineRelecture);
   //	FichierPoints.Close();
   //}

   //



   //xp1 = ParamAffine[0]*(double)x1 + ParamAffine[1]*(double)y1 +ParamAffine[2];
   //yp1 = ParamAffine[3]*(double)x1 + ParamAffine[4]*(double)y1 +ParamAffine[5];

   //xp2 = ParamAffine[0]*(double)x2 + ParamAffine[1]*(double)y2 +ParamAffine[2];
   //yp2 = ParamAffine[3]*(double)x2 + ParamAffine[4]*(double)y2 +ParamAffine[5];

   //xp3 = ParamAffine[0]*(double)x3 + ParamAffine[1]*(double)y3 +ParamAffine[2];
   //yp3 = ParamAffine[3]*(double)x3 + ParamAffine[4]*(double)y3 +ParamAffine[5];*/


   //init LUT à 0
   for(i=0;i<Larg;i++)
      for(j=0;j<Haut;j++)
      {
         LUTPourTransAff[j*Larg+i].x = 0;
         LUTPourTransAff[j*Larg+i].y = 0;
      }
      //calcul LUT
      for(i=0;i<Larg;i++)
         for(j=0;j<Haut;j++)
         {
            imodif = ParamHomogra[0]*(double)i + ParamHomogra[1]*(double)j +ParamHomogra[2];
            jmodif = ParamHomogra[3]*(double)i + ParamHomogra[4]*(double)j +ParamHomogra[5];
            zmodif = ParamHomogra[6]*(double)i + ParamHomogra[7]*(double)j +ParamHomogra[8];

            /*if (zmodif!=0.0f)
            {*/
            LUTPourTransAff[j*Larg+i].x = imodif/zmodif;
            LUTPourTransAff[j*Larg+i].y = jmodif/zmodif;
            /*}*/
         }
         i=277;j=66;
         imodif = ParamHomogra[0]*(double)i + ParamHomogra[1]*(double)j +ParamHomogra[2];
         jmodif = ParamHomogra[3]*(double)i + ParamHomogra[4]*(double)j +ParamHomogra[5];
         zmodif = ParamHomogra[6]*(double)i + ParamHomogra[7]*(double)j +ParamHomogra[8];

         imodif = imodif/zmodif;
         jmodif = jmodif/zmodif;

         ////recherche taille max
         //for(i=0;i<Larg;i++)
         //	for(j=0;j<Haut;j++)
         //	{
         //		if(LUTPourTransAff[j*Larg+i].x > MaxLargeur )
         //			MaxLargeur = LUTPourTransAff[j*Larg+i].x;

         //		if(LUTPourTransAff[j*Larg+i].x < MinLargeur )
         //			MinLargeur = LUTPourTransAff[j*Larg+i].x;

         //		if(LUTPourTransAff[j*Larg+i].y > MaxHauteur)
         //			MaxHauteur= LUTPourTransAff[j*Larg+i].y;

         //		if(LUTPourTransAff[j*Larg+i].y < MinHauteur)
         //			MinHauteur= LUTPourTransAff[j*Larg+i].y;
         //	}

         //	LargeurNewImage = MaxLargeur - MinLargeur;
         //	HauteurNewImage = MaxHauteur - MinHauteur;

         ////Normalisation entre 0 et MaxLargeur/MaxHauteur
         //	for(i=0;i<Larg;i++)
         //		for(j=0;j<Haut;j++)
         //		{
         //			LUTPourTransAff[j*Larg+i].x = LUTPourTransAff[j*Larg+i].x - MinLargeur;
         //			LUTPourTransAff[j*Larg+i].y = LUTPourTransAff[j*Larg+i].y - MinHauteur;
         //		}



         CImageWnd * pImage=NULL;	
         pImage=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);

         SetCursor(LoadCursor(0,IDC_WAIT));

         CMDIFrameWnd *pmd=GetMDIFrame();
         pImage->Create(Titre+" Curvilinear color",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
         for(i=0;i<Larg;i++)
            for(j=0;j<Haut;j++)
            {

               if((LUTPourTransAff[j*Larg+i].y<Haut)&&(LUTPourTransAff[j*Larg+i].y>0)&&(LUTPourTransAff[j*Larg+i].x<Larg)&&(LUTPourTransAff[j*Larg+i].x>0))
               {
                  pImage->ImageRGB[LUTPourTransAff[j*Larg+i].y *Larg+LUTPourTransAff[j*Larg+i].x].r = Image[j*Larg+i];
                  pImage->ImageRGB[LUTPourTransAff[j*Larg+i].y *Larg+LUTPourTransAff[j*Larg+i].x].g = 0;//Image[j*Larg+i];
                  pImage->ImageRGB[LUTPourTransAff[j*Larg+i].y *Larg+LUTPourTransAff[j*Larg+i].x].b = 0;//Image[j*Larg+i];
               }

            }

            /*if ((xp1<Larg)&& (yp1<Haut))
            {
            pImage->ImageRGB[yp1*Larg+xp1].r = 255;
            pImage->ImageRGB[yp1*Larg+xp1].g = 0;
            pImage->ImageRGB[yp1*Larg+xp1].b = 0;
            }
            if ((xp2<Larg)&& (yp2<Haut))
            {
            pImage->ImageRGB[yp2*Larg+xp2].r = 0;
            pImage->ImageRGB[yp2*Larg+xp2].g = 255;
            pImage->ImageRGB[yp2*Larg+xp2].b = 0;
            }

            if ((xp3<Larg)&& (yp3<Haut))
            {
            pImage->ImageRGB[yp2*Larg+xp2].r = 0;
            pImage->ImageRGB[yp2*Larg+xp2].g = 0;
            pImage->ImageRGB[yp2*Larg+xp2].b = 255;
            }*/


            Invalidate(FALSE);
            SetCursor(LoadCursor(0,IDC_ARROW));
}

void CImageWnd::Convolution2DCedricCurvi(BYTE * ImageDst, BYTE * ImageSrc, double ValSigma, int Largeur, int Hauteur)
{

   double sigma;

   sigma = ValSigma;
   SMasque masq;
   int l,c,p;
   double deuxpi=6.28318;
   double sigma2=sigma*sigma;
   double unsdeuxsigc=1.0/(2.0*sigma2);
   double norm=1.0/(deuxpi*sigma2);
   double norm2=1.0/(deuxpi*sigma2*sigma2);
   double norm3=1.0/(deuxpi*sigma2*sigma2*sigma2);
   double x,y;
   double arg;

   // --- parametres du masque ----
   masq.Haut=19;
   masq.Larg=19;
   masq.Centre_col=masq.Larg/2;
   masq.Centre_lig=masq.Haut/2;
   masq.Filtre=new float[masq.Haut*masq.Larg];

   for (l=0;l<masq.Haut;l++)
   {
      for (c=0;c<masq.Larg;c++)
      {
         p=masq.Larg*l+c;
         x=double(c-masq.Centre_col);
         y=double(l-masq.Centre_lig);
         arg=(x*x+y*y)*unsdeuxsigc;
         // gaussien
         masq.Filtre[p]=(float)(norm*exp(-arg));
         // dérivée de gaussien suivant c
         //	masq.Filtre[p]=(float)(-x*norm2*exp(-arg));
         //	derivee seconde
         //	masq.Filtre[p]=(float)(-x*norm2*exp(-arg));

      }
   }

   Filtre_2D(ImageDst,ImageSrc,masq,Largeur,Hauteur);


   delete masq.Filtre;
}
void CImageWnd::Filtre_2D(BYTE * ImageDst, BYTE * ImageSrc, SMasque masq,int Largeur, int Hauteur)
{                                    
   // --- Déclarations des variables -----	



   float  * source;	
   float mini,maxi;

   SImage info;

   // ---------- informations relatives a l'image -----

   info=GetSImage();

   //CTraitImageWnd * pImage=NULL;	
   source = (float   *) malloc( Largeur*Hauteur*sizeof(float));   

   float  * dest;
   dest = (float   *) malloc( Largeur*Hauteur*sizeof(float)); 

   /*if (TypeImage==FLOAT32) {  
   pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,FLOAT32);
   recop_image(ImageFloat,source,Taille);
   }*/
   /*if (TypeImage==MEM || TypeImage==BMP){ 
   pImage=(CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut);*/

   recop_image(ImageSrc,source,Largeur*Hauteur);	
   /*}
   SetCursor(LoadCursor(0,IDC_WAIT));*/
   // ---- Attention le filtre n'est tjrs pas retourné  hé hé ---

   masque_convolution_2D(source,dest,info,masq);  

   if (TypeImage==MEM || TypeImage==BMP){ 
      min_max(dest,&mini,&maxi,Largeur*Hauteur);

      r_affichage(dest,ImageDst,maxi,mini,Hauteur,Largeur,1);	                    
   }  

   /*if (TypeImage==FLOAT32) {
   recop_image(dest,pImage->ImageFloat,Taille);
   pImage->FloatToImage(); }*/

   /*pImage->Create(Titre+" : "+Chaine,WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
   Invalidate(FALSE);*/
   free(source);
   free(dest);		
   /*SetCursor(LoadCursor(0,IDC_ARROW));*/	

}

void CImageWnd::TireCouleurAlea(BYTE * r, BYTE * g, BYTE * b)
{
   *r= (rand()*235)+20;
   *g= (rand()*235)+20;
   *b= (rand()*235)+20;
}

void CImageWnd::IncrementeCouleur(BYTE * r, BYTE * g, BYTE * b,long cpt)
{
   *r= (BYTE) cpt;
   *g= (BYTE)(cpt>>8);
   *b= (BYTE) (cpt>>16);
}

void CImageWnd::OnMatchingRepetabilit507()
{
#ifdef MP
   CIma Image1, Image2;
   CCurviDetect2D CurviDetect2D;
   CString CheminDossier = "C:\\users\\Cedric Lemaitre\\code curvilineaire\\ImageMichal\\";
   CString Fichier1 = "HRI1.bmp";
   CString Fichier2 = "RI2.bmp";


   /*double H[9];*/

   long NbRegion1, NbRegion2;

   //H[0]=-0.9036;
   //H[1]=0.0263;
   //H[2]=480.3369;
   //H[3]=0.0324;
   //H[4]=-1.0045;
   //H[5]=496.1943;
   //H[6]=0.0003;
   //H[7]=7.3634;
   //H[8]=0.9417;

   NbRegion1=2;
   NbRegion2=2;

   Image1.LitBMP(CheminDossier+Fichier1);
   Image2.LitBMP(CheminDossier+Fichier2);

   //CurviDetect2D.MesureRepetabilite(Image1.m_PImage,Image2.m_PImage,Image1.m_Largeur,Image1.m_Hauteur,NbRegion1,NbRegion2,H);
   CurviDetect2D.MesureRepetJ(&Image1,&Image2,"ExportResultRepet.txt",NbRegion1,NbRegion2);
#endif
}



int CImageWnd::IsAnEdge(int *Image,int Larg,CPoint A)
{
   int l,c,p;
   p=Larg*A.y+A.x;
   if (Image[p]!=0) return Image[p]; // if edge at this point

   // in the 3x3 neigboughood, what is the edge ?
   //for (c=A.x-1;c<=A.x+1;c++)
   //{
   //	for (l=A.y-1;l<=A.y+1;l++)
   //	{
   //		
   //	}
   //}
   return 0;
}

#ifdef MP
void CImageWnd::BuildImageFromList(int *Plan1,CList<PosDetectCurvi,PosDetectCurvi>* ListeCurvi)
{
   CCurviDetect2D CurviDetect2D;
   int u;
   bool Init;
   POSITION pos;
   CPoint TabPourQuadri[4];
   PosDetectCurvi PosNM1,PosN;
   for (u=0;u<Taille;u++) {Plan1[u]=0;}	
   Init = FALSE;
   pos = ListeCurvi->GetHeadPosition();
   while(pos)
   {
      if(!Init)//initialisation pour avoir n-1
      {
         PosNM1 = ListeCurvi->GetAt(pos); // on récupère ce qui sera n-1
         ListeCurvi->GetNext(pos);// on pointe sur le suivant
         Init = TRUE;

      }
      else // phase de traitement sur la longueur
      {
         PosN =  ListeCurvi->GetAt(pos);// On récupère la section n 
         if(pos)// si le pointeur de position ok
         {
            TabPourQuadri[0]=PosNM1.PosInit;
            TabPourQuadri[1]=PosNM1.PosFinal;
            TabPourQuadri[3]=PosN.PosInit;
            TabPourQuadri[2]=PosN.PosFinal;
            TraceRegionCurvi = new CRgn;
            if(TraceRegionCurvi->CreatePolygonRgn(TabPourQuadri,4,ALTERNATE))
            {
               CurviDetect2D.RechercheRectEnglobant(TabPourQuadri,&Xmin,&Xmax,&Ymin,&Ymax);
               CurviDetect2D.RempliQuad(Plan1,TraceRegionCurvi,Larg,&Xmin,&Xmax,&Ymin,&Ymax,1);
            }
            delete TraceRegionCurvi;
            PosNM1 = PosN; // n-1=n;
            ListeCurvi->GetNext(pos);// on passe à la section suivante
            if(pos)PosN = ListeCurvi->GetAt(pos);
         }// si la position est valable
      }// else si init ok
   }// fin while tant qu'on est dans liste

}



void CImageWnd::MergeConnectedRegions(CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> & TableauListesPointsFiltreesC2D,CImageWnd * pImage)
{
   //-----------------------------------------------------------------------------------
   // merge connected regions
   //-----------------------------------------------------------------------------------
   int Taille=Larg*Haut;
   int *Plan1;
   int *Plan2;
   Plan1=(int*) malloc(Taille*sizeof(int));
   Plan2=(int*) malloc(Taille*sizeof(int));
   int uTail,uHead;
   int AddSens,ParcoursSens;

   //CIma Test;
   //Test.LitBMP("f:\\johel\\prog\\bmp\\temp.bmp");

   int LongueurTableauListesPointsFiltreesC2D=TableauListesPointsFiltreesC2D.GetCount();
   int NbOverlap;
   int NbOverlapB;
   double DistGrav;
   double SeuilDistGrav=8.0;
   double SeuilLargMoy=2.0;
   int SeuilOverlapp=3;
   int IndiceReg,IndiceRegB;

   int u;
   PosDetectCurvi PosN,PosNHead,PosNTail;
   POSITION pos;
   BYTE r,g,b;
   double LargMoyA,LargMoyB,LargeurA,LargeurB;
   Sfpoint CentreA,CentreB;
   Sfpoint CentreMaxA,CentreMaxB;
   Sfpoint CentreMinA,CentreMinB;

   //#define ANIM

   for(IndiceReg=0;IndiceReg<LongueurTableauListesPointsFiltreesC2D;IndiceReg++)//pour toutes les régions
   {
      BuildImageFromList(Plan1,TableauListesPointsFiltreesC2D[IndiceReg]);
#ifdef ANIM		
      TireCouleurAlea(&r,&g,&b);
      for (u=0;u<Taille;u++)
      {
         if (Plan1[u]==1) pImage->ImageRGB[u].r=r; else pImage->ImageRGB[u].r=0; 
         if (Plan1[u]==1) pImage->ImageRGB[u].g=g; else pImage->ImageRGB[u].g=0;
         if (Plan1[u]==1) pImage->ImageRGB[u].b=b; else pImage->ImageRGB[u].b=0;
      }
      pImage->SendMessage(WM_PAINT);
      pImage->Invalidate(TRUE);
      Sleep(500);
#endif


      for(IndiceRegB=IndiceReg+1;IndiceRegB<LongueurTableauListesPointsFiltreesC2D;IndiceRegB++)//pour toutes les régions
      {
         BuildImageFromList(Plan2,TableauListesPointsFiltreesC2D[IndiceRegB]);
#ifdef ANIM
         TireCouleurAlea(&r,&g,&b);
         for (u=0;u<Taille;u++)
         {
            if (Plan2[u]==1) pImage->ImageRGB[u].r=r;// else pImage->ImageRGB[u].r=0; 
            if (Plan2[u]==1) pImage->ImageRGB[u].g=0; //else pImage->ImageRGB[u].g=0;
            if (Plan2[u]==1) pImage->ImageRGB[u].b=b; //else pImage->ImageRGB[u].b=0;
         }
         pImage->SendMessage(WM_PAINT);
         pImage->Invalidate(TRUE);
         Sleep(500);
#endif

         // compute the gravity center of each axis points in the overlapping region
         //Sfpoint Grav1,Grav2;
         //Grav1.x=Grav1.y=0.0;
         //Grav2.x=Grav2.y=0.0;
         NbOverlap=0;
         NbOverlapB=0;
         AddSens=-1;
         ParcoursSens=-1;
         LargMoyA=0;
         LargMoyB=0;
         LargeurB=LargeurA=0.0;
         CentreMaxA.x=0;
         CentreMaxA.y=0;
         CentreMinA.x=Larg-1;
         CentreMinA.y=Haut-1;
         CentreMaxB.x=0;
         CentreMaxB.y=0;
         CentreMinB.x=Larg-1;
         CentreMinB.y=Haut-1;

         // find the position of the overlapp : head or tail or neither ?
         pos = TableauListesPointsFiltreesC2D[IndiceReg]->GetTailPosition();
         PosNTail =  TableauListesPointsFiltreesC2D[IndiceReg]->GetAt(pos);
         uTail=Larg*(int)PosNTail.PosAxe.y+(int)PosNTail.PosAxe.x;
         pos = TableauListesPointsFiltreesC2D[IndiceReg]->GetHeadPosition();
         PosNHead =  TableauListesPointsFiltreesC2D[IndiceReg]->GetAt(pos);
         uHead=Larg*(int)PosNHead.PosAxe.y+(int)PosNHead.PosAxe.x;
         if (Plan2[uTail]==1 && Plan2[uHead]==1) // A is in B, no merging with this part
         {
            AddSens=2;
         }
         if (Plan2[uTail]==1 && Plan2[uHead]==0) // l'axe tail de A est dans B, we need add tail during merging
         {
            AddSens=1;
            LargeurA=PosNTail.Largeur;
         }
         if (Plan2[uTail]==0 && Plan2[uHead]==1) // l'axe head de A est dans B, we need add head during merging
         {
            AddSens=0;
            LargeurA=PosNHead.Largeur;

         }
         //if (AddSens==-1 || AddSens==2)  // no merging


         if (AddSens==0 || AddSens==1)
         {
            while (pos)
            {
               PosN =  TableauListesPointsFiltreesC2D[IndiceReg]->GetAt(pos);
               u=Larg*(int)PosN.PosAxe.y+(int)PosN.PosAxe.x;
               LargMoyA+=PosN.Largeur;

               if (Plan1[u]==1 && Plan2[u]==1) 
               {
                  //Grav1.x+=PosN.PosAxe.x;
                  //Grav1.y+=PosN.PosAxe.y;
                  NbOverlap++;
                  if (PosN.PosAxe.x>CentreMaxA.x) CentreMaxA.x=PosN.PosAxe.x;
                  if (PosN.PosAxe.y>CentreMaxA.y) CentreMaxA.y=PosN.PosAxe.y;
                  if (PosN.PosAxe.x<CentreMinA.x) CentreMinA.x=PosN.PosAxe.x;
                  if (PosN.PosAxe.y<CentreMinA.y) CentreMinA.y=PosN.PosAxe.y;

               }
               TableauListesPointsFiltreesC2D[IndiceReg]->GetNext(pos);
            }
            LargMoyA/=TableauListesPointsFiltreesC2D[IndiceReg]->GetCount();

            if (NbOverlap>SeuilOverlapp) //  if  overlap, try to merge
            {
               //Grav1.x=Grav1.x/NbOverlap;
               //Grav1.y=Grav1.y/NbOverlap;
               CentreA.x=(CentreMinA.x+CentreMaxA.x)/2.0;
               CentreA.y=(CentreMinA.y+CentreMaxA.y)/2.0;

               pos = TableauListesPointsFiltreesC2D[IndiceRegB]->GetTailPosition();
               PosNTail =  TableauListesPointsFiltreesC2D[IndiceRegB]->GetAt(pos);
               uTail=Larg*(int)PosNTail.PosAxe.y+(int)PosNTail.PosAxe.x;
               pos = TableauListesPointsFiltreesC2D[IndiceRegB]->GetHeadPosition();
               PosNHead =  TableauListesPointsFiltreesC2D[IndiceRegB]->GetAt(pos);
               uHead=Larg*(int)PosNHead.PosAxe.y+(int)PosNHead.PosAxe.x;	
               if (Plan1[uTail]==1 && Plan1[uHead]==1) // B is in A, no merging with this part
               {
                  ParcoursSens=2;
               }
               if (Plan1[uTail]==1 && Plan1[uHead]==0) // l'axe tail de B est dans A, we need PREV during merging
               {
                  ParcoursSens=1;
                  LargeurB=PosNTail.Largeur;

               }
               if (Plan1[uTail]==0 && Plan1[uHead]==1) // l'axe head de B est dans 1, we need NEXT during merging
               {
                  ParcoursSens=0;
                  LargeurB=PosNHead.Largeur;
               }
               if (ParcoursSens==1 || ParcoursSens==0)
               {
                  pos = TableauListesPointsFiltreesC2D[IndiceRegB]->GetHeadPosition();
                  int numpt=0;
                  while (pos)
                  {
                     PosN =  TableauListesPointsFiltreesC2D[IndiceRegB]->GetAt(pos);
                     LargMoyB+=PosN.Largeur;
                     u=Larg*(int)PosN.PosAxe.y+(int)PosN.PosAxe.x;
                     if (Plan1[u]==1 && Plan2[u]==1) 
                     {
                        //	Grav2.x+=PosN.PosAxe.x;
                        //	Grav2.y+=PosN.PosAxe.y;
                        NbOverlapB++;
                        if (PosN.PosAxe.x>CentreMaxB.x) CentreMaxB.x=PosN.PosAxe.x;
                        if (PosN.PosAxe.y>CentreMaxB.y) CentreMaxB.y=PosN.PosAxe.y;
                        if (PosN.PosAxe.x<CentreMinB.x) CentreMinB.x=PosN.PosAxe.x;
                        if (PosN.PosAxe.y<CentreMinB.y) CentreMinB.y=PosN.PosAxe.y;
                        if (ParcoursSens==0) // on prend la derniere largeur comme LargeurB
                        {
                           LargeurB=PosN.Largeur;
                        }

                        if (ParcoursSens==1) // on prend la premiere largeur comme LargeurB
                        {
                           if (numpt==0) { LargeurB=PosN.Largeur; numpt++;}
                        }

                     }
                     TableauListesPointsFiltreesC2D[IndiceRegB]->GetNext(pos);// on passe à la section suivante
                  }
                  LargMoyB/=TableauListesPointsFiltreesC2D[IndiceRegB]->GetCount();

                  if (NbOverlapB!=0)
                  {
                     //Grav2.x=Grav2.x/NbOverlapB;
                     //Grav2.y=Grav2.y/NbOverlapB;
                     CentreB.x=(CentreMinB.x+CentreMaxB.x)/2.0;
                     CentreB.y=(CentreMinB.y+CentreMaxB.y)/2.0;

                  } else 
                  {int bugggg=0; }
                  // if the distance between the gravity centers is low, merge the two lists

                  //DistGrav=(Grav1.x-Grav2.x)*(Grav1.x-Grav2.x)+(Grav1.y-Grav2.y)*(Grav1.y-Grav2.y);
                  DistGrav=(CentreA.x-CentreB.x)*(CentreA.x-CentreB.x)+(CentreA.y-CentreB.y)*(CentreA.y-CentreB.y);
                  DistGrav=sqrt(DistGrav);
                  //		double DistLargMoy=fabs(LargMoyB-LargMoyA)/log(__max(LargMoyA,LargMoyB));
                  double DistLargMoy=fabs(LargeurB-LargeurA)/log(__max(LargeurA,LargeurB));
#ifdef ANIM
                  u=Larg*(int)CentreA.y+(int)CentreA.x;
                  pImage->ImageRGB[u].r=255; pImage->ImageRGB[u].g=255; pImage->ImageRGB[u].b=255;
                  u=Larg*(int)CentreB.y+(int)CentreB.x;
                  pImage->ImageRGB[u].r=0; pImage->ImageRGB[u].g=255; pImage->ImageRGB[u].b=0;
                  CString Tempo;
                  Tempo.Format("DistLargMoy=%f\n distgrav=%f",DistLargMoy,DistGrav);
                  MessageBox(Tempo,"");
#endif
                  if (DistGrav<__max(LargMoyA,LargMoyB)/4 && DistLargMoy<SeuilLargMoy) //  merge
                  {
#ifdef ANIM 
                     MessageBox("Fusion");
#endif
                     if (ParcoursSens==0)
                     {
                        int NbpAj=0;
                        pos = TableauListesPointsFiltreesC2D[IndiceRegB]->GetHeadPosition();
                        while (pos)
                        {
                           PosN =  TableauListesPointsFiltreesC2D[IndiceRegB]->GetAt(pos);
                           u=Larg*(int)PosN.PosAxe.y+(int)PosN.PosAxe.x;
                           if (Plan1[u]==0)  // if the axis is not already in A region
                           {
                              if (AddSens==1) TableauListesPointsFiltreesC2D[IndiceReg]->AddTail(PosN);
                              if (AddSens==0)  TableauListesPointsFiltreesC2D[IndiceReg]->AddHead(PosN);
                           }
                           TableauListesPointsFiltreesC2D[IndiceRegB]->GetNext(pos);
                           NbpAj++;
                        } // while pos
                        /*CString Tempo;
                        Tempo.Format("Nb Sens 0=%d",NbpAj);
                        MessageBox(Tempo,"");*/
                        //CStdioFile FichierTemp;
                        //FichierTemp.Open("c:\\liste.txt",CFile::modeCreate|CFile::modeWrite);
                        //CString Tempo;
                        //pos = TableauListesPointsFiltreesC2D[IndiceReg]->GetHeadPosition();
                        //while (pos)
                        //{
                        //	PosN =  TableauListesPointsFiltreesC2D[IndiceReg]->GetAt(pos);
                        //	Tempo.Format("%d;%d;%d;%d\n",PosN.PosInit.x,PosN.PosInit.y,PosN.PosFinal.x,PosN.PosFinal.y);
                        //	FichierTemp.WriteString(Tempo);
                        //	TableauListesPointsFiltreesC2D[IndiceReg]->GetNext(pos);
                        //}
                        //			FichierTemp.Close();
                     }

                     if (ParcoursSens==1)
                     {
                        int NbpAj=0;
                        pos = TableauListesPointsFiltreesC2D[IndiceRegB]->GetTailPosition();
                        while (pos)
                        {
                           PosN =  TableauListesPointsFiltreesC2D[IndiceRegB]->GetAt(pos);
                           u=Larg*(int)PosN.PosAxe.y+(int)PosN.PosAxe.x;
                           if (Plan1[u]==0 )  // if the axis is not already in the A region
                           {
                              if (AddSens==1) TableauListesPointsFiltreesC2D[IndiceReg]->AddTail(PosN);
                              if (AddSens==0) TableauListesPointsFiltreesC2D[IndiceReg]->AddHead(PosN);
                           }
                           TableauListesPointsFiltreesC2D[IndiceRegB]->GetPrev(pos);
                           NbpAj++;
                        }
                        /*CString Tempo;
                        Tempo.Format("Nb sens 1=%d",NbpAj);
                        MessageBox(Tempo,"");*/
                     }

                     // update the plan1 from merged regions
                     BuildImageFromList(Plan1,TableauListesPointsFiltreesC2D[IndiceReg]);
#ifdef ANIM
                     for (u=0;u<Taille;u++)
                     {
                        if (Plan1[u]==1) pImage->ImageRGB[u].r=255; else pImage->ImageRGB[u].r=0; 
                        if (Plan1[u]==1) pImage->ImageRGB[u].g=0; else pImage->ImageRGB[u].g=0;
                        if (Plan1[u]==1) pImage->ImageRGB[u].b=0; else pImage->ImageRGB[u].b=0;
                     }
                     pImage->SendMessage(WM_PAINT);
                     pImage->Invalidate(TRUE);
                     Sleep(2000);
#endif
                     // remove the region 2 merged in the region 1
                     TableauListesPointsFiltreesC2D[IndiceRegB]->RemoveAll();
                     delete TableauListesPointsFiltreesC2D[IndiceRegB];
                     TableauListesPointsFiltreesC2D.RemoveAt(IndiceRegB);
                     LongueurTableauListesPointsFiltreesC2D--;
                     IndiceRegB=IndiceReg; // start again because the A region changed ! will be  incremented by the loop
                  } // end if merging is possible
               } // end if B is in A

            } // end if no overlap
         } // end if no merging

      }// fin for IndiceRegB

   }//fin for IndiceReg
   free(Plan1);
   free(Plan2);
#ifdef ANIM
   for (u=0;u<Taille;u++)
   {
      pImage->ImageRGB[u].r=0;
      pImage->ImageRGB[u].g=0; 
      pImage->ImageRGB[u].b=0; 
   }
#endif
}



void CImageWnd::SectionOptimization(CImageWnd * pImageTracking,CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> & TableauListesPointsFiltrees)
{
   POSITION pos,PositionCListeCopy;
   PosDetectCurvi PosNM1,PosN,PosNP1;
   int i,k,IndiceReg,longueur;
   double LargTemp;
   int LengthLocalLine=7;
   PosDetectCurvi *pLocalLine;
   Sfpoint *pPosInit;
   Sfpoint *pPosTemp;
   Sfpoint *pPosFinal;
   Sfpoint GravInit;
   Sfpoint PointTwo;
   Sfpoint PointTwoPerp;

   pLocalLine=(PosDetectCurvi *) malloc(LengthLocalLine*sizeof(PosDetectCurvi));
   pPosInit=(Sfpoint *) malloc(LengthLocalLine*sizeof(Sfpoint));
   pPosTemp=(Sfpoint *) malloc(LengthLocalLine*sizeof(Sfpoint));
   pPosFinal=(Sfpoint *) malloc(LengthLocalLine*sizeof(Sfpoint));

   STabIndice *TabR=(STabIndice*)malloc(LengthLocalLine*sizeof(STabIndice));
   STabIndice *TabRInv=(STabIndice*)malloc(LengthLocalLine*sizeof(STabIndice));
   double VarInit,CovInit;
   double Dist, DistMin;
   int IndiceMin;
   CCurviDetect2D CurviDetect2D;
   CPoint * TabCoupe = (CPoint*) malloc(sizeof(CPoint)* (Haut+Larg));


   int longueurTableauListesPointsFiltrees=TableauListesPointsFiltrees.GetCount();

   for(IndiceReg=0;IndiceReg<longueurTableauListesPointsFiltrees;IndiceReg++)//pour toutes les régions
   {


      // visualisation BEFORE optim
      //pos = TableauListesPointsFiltrees[IndiceReg]->GetHeadPosition();
      //while(pos)
      //{
      //	PosN = TableauListesPointsFiltrees[IndiceReg]->GetAt(pos); 
      //	TableauListesPointsFiltrees[IndiceReg]->GetNext(pos);// on pointe sur le suivant
      //	int SigneDx;
      //	int Sens=-1;
      //	int Centre=0;
      //	int w=PosN.Largeur;
      //	float GradMoy_x=PosN.PosInit.x-PosN.PosFinal.x;
      //	float GradMoy_y=PosN.PosInit.y-PosN.PosFinal.y;
      //	if(GradMoy_x>0)
      //		SigneDx = -Sens; 
      //	else
      //		SigneDx=+Sens;
      //	if(GradMoy_x==0 && GradMoy_y>=0.0) SigneDx=-Sens;
      //	if(GradMoy_x==0 && GradMoy_y<0.0) SigneDx=+Sens;
      //	CurviDetect2D.RemplitTableauIJCoupeBiDir(TabCoupe,PosN.PosAxe.x,PosN.PosAxe.y,GradMoy_x,GradMoy_y,w,0,SigneDx,&Centre,Haut,Larg);
      //	for(int j=Centre-w/2;j<Centre+w/2;j++)//Gauche
      //	{
      //		long p=TabCoupe[j].x+Larg*TabCoupe[j].y;
      //		pImageTracking->ImageRGB[p].r =  (j+20)*2;; 
      //		pImageTracking->ImageRGB[p].g =  0;
      //		pImageTracking->ImageRGB[p].b =  0;
      //	}
      //		long p=(long)PosN.PosInit.x+Larg*(long)PosN.PosInit.y;
      //		pImageTracking->ImageRGB[p].r =  255; 
      //		pImageTracking->ImageRGB[p].g =  255;
      //		pImageTracking->ImageRGB[p].b =  255;
      //		p=(long)PosN.PosAxe.x+Larg*(long)PosN.PosAxe.y;
      //		pImageTracking->ImageRGB[p].r =  0; 
      //		pImageTracking->ImageRGB[p].g =  255;
      //		pImageTracking->ImageRGB[p].b =  0;

      //	pImageTracking->SendMessage(WM_PAINT);
      //	pImageTracking->Invalidate(TRUE);
      //	Sleep(100);
      //}



      pos = TableauListesPointsFiltrees[IndiceReg]->GetHeadPosition();
      //	ListePointsFiltresC2D = new CList <PosDetectCurvi,PosDetectCurvi>;//on alloue une nouvelle liste
      bool fin = false;
      while(pos)
      {
         GravInit.x=GravInit.y=0.0;
         for (i=0;i<LengthLocalLine;i++)
         {
            TabR[i].i=i;
            if (pos) pLocalLine[i] = TableauListesPointsFiltrees[IndiceReg]->GetAt(pos); 
            else {fin =true; break;}
            pPosInit[i].x=(double)pLocalLine[i].PosInit.x;
            pPosInit[i].y=(double)pLocalLine[i].PosInit.y;
            GravInit.x+=pPosInit[i].x;
            GravInit.y+=pPosInit[i].y;
            pPosTemp[i].x=(double)pLocalLine[i].PosFinal.x;
            pPosTemp[i].y=(double)pLocalLine[i].PosFinal.y;
            TableauListesPointsFiltrees[IndiceReg]->GetNext(pos);// on pointe sur le suivant

         }
         if (!fin)
         {
            GravInit.x/=LengthLocalLine;
            GravInit.y/=LengthLocalLine;
            if (pos) for (i=0;i<LengthLocalLine;i++) // rewind the list if we are not at the end
            {
               TableauListesPointsFiltrees[IndiceReg]->GetPrev(pos); 
            }
            VarInit=Variance(pPosInit,LengthLocalLine,GravInit.x);
            CovInit=Covariance(pPosInit,LengthLocalLine,GravInit);
            PointTwo.x=0.0;
            PointTwoPerp.x=0;
            if (VarInit==0.0) // if vertical
            {
               PointTwo.x=GravInit.x;
               PointTwo.y=GravInit.y+100.0;
               PointTwoPerp.y=GravInit.y; // perp is horizontal
               PointTwoPerp.x=GravInit.x+100 ;
            } 
            else
            {
               PointTwo.y=GravInit.y-(CovInit*GravInit.x)/VarInit;
               if (CovInit==0.0)  {PointTwoPerp.y=GravInit.x+100;PointTwoPerp.x=GravInit.x;}
               else
                  PointTwoPerp.y=GravInit.y+(VarInit/CovInit)*GravInit.x;
            }

            //int SigneDx;
            //int Sens=-1;
            //int Centre=0;
            //float GradMoy_x=GravInit.x-PointTwo.x;
            //float GradMoy_y=GravInit.y-PointTwo.y;
            //if(GradMoy_x>0)
            //	SigneDx = -Sens; 
            //else
            //	SigneDx=+Sens;
            //if(GradMoy_x==0 && GradMoy_y>=0.0) SigneDx=-Sens;
            //if(GradMoy_x==0 && GradMoy_y<0.0) SigneDx=+Sens;
            //CurviDetect2D.RemplitTableauIJCoupeBiDir(TabCoupe,(int)GravInit.x,(int)GravInit.y,GradMoy_x,GradMoy_y,40,10,SigneDx,&Centre,Haut,Larg);
            //for(int j=Centre-30;j<Centre+30;j++)//Gauche
            //{
            //	long p=TabCoupe[j].x+Larg*TabCoupe[j].y;
            //	pImageTracking->ImageRGB[p].r =  0; 
            //	pImageTracking->ImageRGB[p].g =  (j+20)*2;
            //	pImageTracking->ImageRGB[p].b =  0;
            //}
            //long p=(long)GravInit.x+Larg*(long)GravInit.y;
            //pImageTracking->ImageRGB[p].r =  255; 
            //pImageTracking->ImageRGB[p].g =  255;
            //pImageTracking->ImageRGB[p].b =  255;
            //pImageTracking->SendMessage(WM_PAINT);
            //pImageTracking->Invalidate(TRUE);
            //Sleep(1);

            // point A = centre de gravite
            // point B = pointtwo
            // point C = point final à projeter

            for (i=0;i<LengthLocalLine;i++)
            {
               TabR[i].r=Scalaire(GravInit,PointTwo,pPosTemp[i]);
            }
            qsort(TabR,LengthLocalLine,sizeof(STabIndice),CImageWnd::CompareDouble);
            InvertTable(TabR,TabRInv,LengthLocalLine);

            for (i=1;i<LengthLocalLine-1;i++)
            {
               //DistMin=__max((double)Larg,(double)Haut);
               //pPosFinal[i].x=pPosTemp[TabRInv[i].i].x; // default value to the ordered point
               //pPosFinal[i].y=pPosTemp[TabRInv[i].i].y; // default value to the ordered point
               pPosFinal[i].x=pPosTemp[i].x; // default value to the original point
               pPosFinal[i].y=pPosTemp[i].y; // default value to the original point
               DistMin=Distance(pPosFinal[i],pPosInit[i]);
               double Scal;
               Sfpoint Final;
               // find the nearest point from ordered final points or projection
               for (k=0;k<LengthLocalLine-1;k++)
               {
                  Final=Projection(pPosTemp[TabRInv[k].i],pPosTemp[TabRInv[k+1].i],pPosInit[i],&Scal);
                  Dist=Distance(Final,pPosInit[i]);
                  if (Dist<DistMin && Scal>=0.0 && Scal <=1.0)  
                  { 
                     DistMin=Dist; 
                     IndiceMin=k;
                     pPosFinal[i].x=Final.x; // optimised value
                     pPosFinal[i].y=Final.y; // optimised value
                  }
               } // for all projections
            } // for initial points

            // replace the modifed points and recompute the new axis position
            if (pos) TableauListesPointsFiltrees[IndiceReg]->GetNext(pos);// on pointe sur le suivant
            for (i=1;i<LengthLocalLine-1;i++)
            {
               pLocalLine[i].PosFinal.x=(int)pPosFinal[i].x;
               pLocalLine[i].PosFinal.y=(int)pPosFinal[i].y;
               pLocalLine[i].PosAxe.x=(pPosInit[i].x+pPosFinal[i].x)/2;
               pLocalLine[i].PosAxe.y=(pPosInit[i].y+pPosFinal[i].y)/2;
               pLocalLine[i].Largeur=Distance(pPosInit[i],pPosFinal[i]);
               if (pos) TableauListesPointsFiltrees[IndiceReg]->SetAt(pos,pLocalLine[i]);
               if (pos) TableauListesPointsFiltrees[IndiceReg]->GetNext(pos);// on pointe sur le suivant
            }
            if (pos) for (i=1;i<LengthLocalLine/4-1;i++) // rewind the list but one if we are not at the end
            {
               TableauListesPointsFiltrees[IndiceReg]->GetPrev(pos); 
            }

         } // if not fin
      } // while pos

      // visualisation
      //pos = TableauListesPointsFiltrees[IndiceReg]->GetHeadPosition();
      //while(pos)
      //{
      //	PosN = TableauListesPointsFiltrees[IndiceReg]->GetAt(pos); 
      //	TableauListesPointsFiltrees[IndiceReg]->GetNext(pos);// on pointe sur le suivant
      //	int SigneDx;
      //	int Sens=-1;
      //	int Centre=0;
      //	int w=PosN.Largeur;
      //	float GradMoy_x=PosN.PosInit.x-PosN.PosFinal.x;
      //	float GradMoy_y=PosN.PosInit.y-PosN.PosFinal.y;
      //	if(GradMoy_x>0)
      //		SigneDx = -Sens; 
      //	else
      //		SigneDx=+Sens;
      //	if(GradMoy_x==0 && GradMoy_y>=0.0) SigneDx=-Sens;
      //	if(GradMoy_x==0 && GradMoy_y<0.0) SigneDx=+Sens;
      //	CurviDetect2D.RemplitTableauIJCoupeBiDir(TabCoupe,PosN.PosAxe.x,PosN.PosAxe.y,GradMoy_x,GradMoy_y,w,0,SigneDx,&Centre,Haut,Larg);
      //	for(int j=Centre-w/2;j<Centre+w/2;j++)//Gauche
      //	{
      //		long p=TabCoupe[j].x+Larg*TabCoupe[j].y;
      //		pImageTracking->ImageRGB[p].r =  (j+20)*2;; 
      //		pImageTracking->ImageRGB[p].g =  0;
      //		pImageTracking->ImageRGB[p].b =  0;
      //	}
      //		long p=(long)PosN.PosInit.x+Larg*(long)PosN.PosInit.y;
      //		pImageTracking->ImageRGB[p].r =  255; 
      //		pImageTracking->ImageRGB[p].g =  255;
      //		pImageTracking->ImageRGB[p].b =  255;
      //		p=(long)PosN.PosAxe.x+Larg*(long)PosN.PosAxe.y;
      //		pImageTracking->ImageRGB[p].r =  0; 
      //		pImageTracking->ImageRGB[p].g =  255;
      //		pImageTracking->ImageRGB[p].b =  0;

      //	pImageTracking->SendMessage(WM_PAINT);
      //	pImageTracking->Invalidate(TRUE);
      //	Sleep(100);
      //}

   } // for all regions

   free(pLocalLine);
   free(TabR);
   free(TabRInv);
   free(pPosTemp);
   free(pPosInit);
   free(pPosFinal);
   free (TabCoupe);

}


void CImageWnd::TwoDimConstraints(CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> & TableauListesPointsFiltreesC2D,CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> & TableauListesPointsFiltrees,int m_LongueurMiniRegionFiltrage,int SeuilDepass)
{
   //------------------------------------------
   // Contraintes 2D
   //------------------------------------------
   int IndiceReg,longueur,cptDepassw=0;
   double LargTemp;
   double ddw,DistAxe;
   double AngleCourb;
   double SeuilCourb=1.0;
   CList <PosDetectCurvi,PosDetectCurvi> * ListePointsFiltresC2D;
   BOOL Init = FALSE;
   POSITION pos,PositionCListeCopy;
   PosDetectCurvi PosNM1,PosN,PosNP1;


   //PosDetectCurvi *TabPointTemp=(PosDetectCurvi *)malloc(sizeof (PosDetectCurvi)*SeuilDepass);
   int longueurTableauListesPointsFiltrees=TableauListesPointsFiltrees.GetCount();
   for(IndiceReg=0;IndiceReg<longueurTableauListesPointsFiltrees;IndiceReg++)//pour toutes les régions
   {
      pos = TableauListesPointsFiltrees[IndiceReg]->GetHeadPosition();
      Init = FALSE;
      if(TableauListesPointsFiltrees[IndiceReg]->GetCount()>m_LongueurMiniRegionFiltrage)// si la liste n'est pas vide
      {	
         ListePointsFiltresC2D = new CList <PosDetectCurvi,PosDetectCurvi>;//on alloue une nouvelle liste

         while(pos)
         {
            if(!Init)//initialisation pour avoir n-1
            {
               PosNM1 = TableauListesPointsFiltrees[IndiceReg]->GetAt(pos); // on récupère ce qui sera n-1
               TableauListesPointsFiltrees[IndiceReg]->GetNext(pos);// on pointure le suivant
               PositionCListeCopy = ListePointsFiltresC2D->AddHead(PosNM1);
               Init = TRUE;
               cptDepassw=0;
               LargTemp=(double)PosNM1.Largeur;

            }
            else // phase de traitement sur la largeur
            {
               PosN =  TableauListesPointsFiltrees[IndiceReg]->GetAt(pos);// On récupère la section n 
               AngleCourb=-10;
               if(pos)// si le pointeur de position ok
               {
                  //PosN=TableauListesPointsFiltrees[IndiceReg]->GetAt(pos);

                  ddw=fabs((double)PosN.Largeur-LargTemp);
                  ddw=ddw/log(LargTemp);


                  TableauListesPointsFiltrees[IndiceReg]->GetNext(pos);
                  if (pos) 
                  {
                     PosNP1=TableauListesPointsFiltrees[IndiceReg]->GetAt(pos);
                     TableauListesPointsFiltrees[IndiceReg]->GetPrev(pos);
                     if (ddw>(VLocaleLargeur*4.0) )
                        //if (0)
                     {
                        cptDepassw++;
                        if (cptDepassw>=SeuilDepass)
                        {
                           if (ListePointsFiltresC2D->GetCount()<m_LongueurMiniRegionFiltrage) delete ListePointsFiltresC2D;
                           else TableauListesPointsFiltreesC2D.Add(ListePointsFiltresC2D);
                           ListePointsFiltresC2D = new CList <PosDetectCurvi,PosDetectCurvi>;//on alloue une nouvelle liste parce que la section en cours ne respecte pas la contrainte
                           for (int tt=0;tt<SeuilDepass;tt++) TableauListesPointsFiltrees[IndiceReg]->GetPrev(pos);
                           Init=FALSE; // retour au debut
                           cptDepassw=0;
                           LargTemp=(double)PosN.Largeur;
                        }
                     }
                     else
                     {
                        // no more need of local curvature neither axe distance

                        //AngleCourb=CurviDetect2D.AngleCourbureLocale(PosNM1.PosAxe,PosN.PosAxe,PosNP1.PosAxe);
                        //DistAxe=(PosN.PosAxe.x-PosNM1.PosAxe.x)*(PosN.PosAxe.x-PosNM1.PosAxe.x);
                        //DistAxe+=(PosN.PosAxe.y-PosNM1.PosAxe.y)*(PosN.PosAxe.y-PosNM1.PosAxe.y);
                        //DistAxe=sqrt(DistAxe)*2.0;
                        //if ((DistAxe>(double)PosNM1.Largeur && AngleCourb<SeuilCourb))
                        if (0)
                        {
                           if (ListePointsFiltresC2D->GetCount()<m_LongueurMiniRegionFiltrage) delete ListePointsFiltresC2D;
                           else TableauListesPointsFiltreesC2D.Add(ListePointsFiltresC2D);
                           ListePointsFiltresC2D = new CList <PosDetectCurvi,PosDetectCurvi>;//on alloue une nouvelle liste parce que la section en cours ne respecte pas la contrainte
                           ListePointsFiltresC2D->AddHead(PosN);//on ajoute la section dans la liste en cours

                           Init=FALSE; // retour au debut
                           cptDepassw=0;
                           LargTemp=(double)PosN.Largeur;
                        }
                        else
                        {
                           ListePointsFiltresC2D->AddHead(PosN);//on ajoute la section dans la liste en cours
                           cptDepassw=0;
                           LargTemp=(double)PosN.Largeur;
                        }

                     }
                     PosNM1 = PosN; // n-1=n;
                     TableauListesPointsFiltrees[IndiceReg]->GetNext(pos);// on passe à la section suivante
                     //if(pos)PosN = TableauListesPointsFiltrees[IndiceReg]->GetAt(pos);

                  }// si la position est valable
               } // si pos
            }// else si init ok
         }// fin while tant qu'on est dans liste


         longueur=ListePointsFiltresC2D->GetCount();
         //if (longueur>m_LongueurMiniForme)	
         TableauListesPointsFiltreesC2D.Add(ListePointsFiltresC2D);//on ajoute la liste en cours dans le tableau de liste
      }// si les listes ne sont pas vide
   }//fin for IndiceReg
}

void CImageWnd::TryToMerge(CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> & TableauListesPointsFiltreesC2D,
                           CArray <CList<int,int>*,CList<int,int>*> &TabNumListesConnex,
                           int *TabListeOk,
                           CImageWnd * pImage)
{
   //-----------------------------------------------------------------------------------
   // merge connected regions
   //-----------------------------------------------------------------------------------
   //#define ANIM


   int *Plan1;
   int *Plan2;
   int Taille=Larg*Haut;
   Plan1=(int*) malloc(Taille*sizeof(int));
   Plan2=(int*) malloc(Taille*sizeof(int));
   int uTail,uHead;
   int AddSens,ParcoursSens;

   //CIma Test;
   //Test.LitBMP("f:\\johel\\prog\\bmp\\temp.bmp");
   int LongueurTableauListesPointsFiltreesC2D=TableauListesPointsFiltreesC2D.GetCount();
   int NbOverlap;
   int NbOverlapB;
   double DistGrav;
   double SeuilDistGrav=8.0;
   double SeuilLargMoy=2.0;
   double SumDis;
   double MinDisPt;
   double DistPtAxe;
   double VarianceA,VarianceB,CovarianceA,CovarianceB;
   int SeuilOverlapp=3;
   int IndiceReg,IndiceRegB;

   int u;
   PosDetectCurvi PosN,PosNHead,PosNTail;
   POSITION pos;
   BYTE r,g,b;
   double LargMoyA,LargMoyB,LargeurA,LargeurB;
   Sfpoint CentreA,CentreB;
   Sfpoint CentreMaxA,CentreMaxB;
   Sfpoint CentreMinA,CentreMinB;
   Sfpoint Grav1,Grav2;

   CList<PosDetectCurvi,PosDetectCurvi>*ListeA;
   CList<PosDetectCurvi,PosDetectCurvi>*ListeB;
   bool MergeOk;
   int NbRegConnex;
   POSITION posA,posB,posC;
   int IndiceRegConnected,IndiceRegConnectedToB;
   int n,k;
   CArray<Sfpoint,Sfpoint> ListeAxeA;
   CArray<Sfpoint,Sfpoint> ListeAxeB;

   for(IndiceReg=0;IndiceReg<LongueurTableauListesPointsFiltreesC2D;IndiceReg++)//pour toutes les régions
   {
      if (TabListeOk[IndiceReg]==1)
      {
         ListeA=TableauListesPointsFiltreesC2D[IndiceReg];
         BuildImageFromList(Plan1,ListeA);
#ifdef ANIM		
         TireCouleurAlea(&r,&g,&b);
         for (u=0;u<Taille;u++)
         {
            if (Plan1[u]==1) pImage->ImageRGB[u].r=0; else pImage->ImageRGB[u].r=0; 
            if (Plan1[u]==1) pImage->ImageRGB[u].g=255; else pImage->ImageRGB[u].g=0;
            if (Plan1[u]==1) pImage->ImageRGB[u].b=0; else pImage->ImageRGB[u].b=0;
         }
         pImage->SendMessage(WM_PAINT);
         pImage->Invalidate(TRUE);
         Sleep(500);
#endif

         NbRegConnex=TabNumListesConnex[IndiceReg]->GetCount();
         posA=TabNumListesConnex[IndiceReg]->GetHeadPosition();
         while(posA) // for all box-connected lists
         {
            IndiceRegB=TabNumListesConnex[IndiceReg]->GetAt(posA);
            ListeB=TableauListesPointsFiltreesC2D[IndiceRegB];
            BuildImageFromList(Plan2,ListeB);
#ifdef ANIM
            TireCouleurAlea(&r,&g,&b);
            for (u=0;u<Taille;u++)
            {
               if (Plan2[u]==1) pImage->ImageRGB[u].r=0;// else pImage->ImageRGB[u].r=0; 
               if (Plan2[u]==1) pImage->ImageRGB[u].g=0; //else pImage->ImageRGB[u].g=0;
               if (Plan2[u]==1) pImage->ImageRGB[u].b=255; //else pImage->ImageRGB[u].b=0;
            }
            pImage->SendMessage(WM_PAINT);
            pImage->Invalidate(TRUE);
            Sleep(500);
#endif

            // try to merge A and B
            // if merge ok, mark the region B as merged
            MergeOk=false;
            if (TabListeOk[IndiceRegB]==1)// if list not already merged
            {
               NbOverlap=0;
               NbOverlapB=0;
               AddSens=-1;
               ParcoursSens=-1;
               LargMoyA=0;
               LargMoyB=0;
               LargeurB=LargeurA=0.0;
               CentreMaxA.x=0;
               CentreMaxA.y=0;
               CentreMinA.x=Larg-1;
               CentreMinA.y=Haut-1;
               CentreMaxB.x=0;
               CentreMaxB.y=0;
               CentreMinB.x=Larg-1;
               CentreMinB.y=Haut-1;
               Grav1.x=Grav1.y=Grav2.x=Grav2.y=0.0;
               SumDis=0;


               // find the position of the overlapp : head or tail or neither ?
               pos = ListeA->GetTailPosition();
               PosNTail =  ListeA->GetAt(pos);
               uTail=Larg*(int)PosNTail.PosAxe.y+(int)PosNTail.PosAxe.x;
               pos = ListeA->GetHeadPosition();
               PosNHead =  ListeA->GetAt(pos);
               uHead=Larg*(int)PosNHead.PosAxe.y+(int)PosNHead.PosAxe.x;
               if (Plan2[uTail]==1 && Plan2[uHead]==1) // A is in B, no merging with this part
               {
                  AddSens=2;
                  //compute the average width
                  // if average width are similar, mark A as removed
                  //LargeurA=PosNTail.Largeur; // bug 17 /7
               }
               if (Plan2[uTail]==1 && Plan2[uHead]==0) // l'axe tail de A est dans B, we need add tail during merging
               {
                  AddSens=1;
                  LargeurA=PosNTail.Largeur;
               }
               if (Plan2[uTail]==0 && Plan2[uHead]==1) // l'axe head de A est dans B, we need add head during merging
               {
                  AddSens=0;
                  LargeurA=PosNHead.Largeur;

               }
               //if (AddSens==-1 || AddSens==2)  // no merging


               if (AddSens==0 || AddSens==1)
                  //			if (AddSens!=-1)
               {
                  ListeAxeA.RemoveAll();
                  while (pos)
                  {
                     PosN =  ListeA->GetAt(pos);
                     u=Larg*(int)PosN.PosAxe.y+(int)PosN.PosAxe.x;
                     LargMoyA+=PosN.Largeur;

                     if (Plan1[u]==1 && Plan2[u]==1) 
                     {
                        Grav1.x+=PosN.PosAxe.x;
                        Grav1.y+=PosN.PosAxe.y;
                        NbOverlap++;
                        if (PosN.PosAxe.x>CentreMaxA.x) CentreMaxA.x=PosN.PosAxe.x;
                        if (PosN.PosAxe.y>CentreMaxA.y) CentreMaxA.y=PosN.PosAxe.y;
                        if (PosN.PosAxe.x<CentreMinA.x) CentreMinA.x=PosN.PosAxe.x;
                        if (PosN.PosAxe.y<CentreMinA.y) CentreMinA.y=PosN.PosAxe.y;
                        ListeAxeA.Add(PosN.PosAxe);
                     }
                     ListeA->GetNext(pos);
                  }
                  LargMoyA/=ListeA->GetCount();

                  if (NbOverlap>SeuilOverlapp) //  if  overlap, try to merge
                  {
                     Grav1.x=Grav1.x/NbOverlap;
                     Grav1.y=Grav1.y/NbOverlap;
                     CentreA.x=(CentreMinA.x+CentreMaxA.x)/2.0;
                     CentreA.y=(CentreMinA.y+CentreMaxA.y)/2.0;
                     ListeAxeB.RemoveAll();
                     pos = ListeB->GetTailPosition();
                     PosNTail =  ListeB->GetAt(pos);
                     uTail=Larg*(int)PosNTail.PosAxe.y+(int)PosNTail.PosAxe.x;
                     pos = ListeB->GetHeadPosition();
                     PosNHead =  ListeB->GetAt(pos);
                     uHead=Larg*(int)PosNHead.PosAxe.y+(int)PosNHead.PosAxe.x;	
                     if (Plan1[uTail]==1 && Plan1[uHead]==1) // B is in A, no merging with this part
                     {
                        ParcoursSens=2;
                        //LargeurB=PosNTail.Largeur;
                     }
                     if (Plan1[uTail]==1 && Plan1[uHead]==0) // l'axe tail de B est dans A, we need PREV during merging
                     {
                        ParcoursSens=1;
                        LargeurB=PosNTail.Largeur;

                     }
                     if (Plan1[uTail]==0 && Plan1[uHead]==1) // l'axe head de B est dans 1, we need NEXT during merging
                     {
                        ParcoursSens=0;
                        LargeurB=PosNHead.Largeur;
                     }

                     if (ParcoursSens==1 || ParcoursSens==0)
                        //					if (ParcoursSens!=-1)
                     {
                        pos = ListeB->GetHeadPosition();
                        int numpt=0;
                        while (pos)
                        {
                           PosN =  ListeB->GetAt(pos);
                           LargMoyB+=PosN.Largeur;
                           u=Larg*(int)PosN.PosAxe.y+(int)PosN.PosAxe.x;
                           if (Plan1[u]==1 && Plan2[u]==1) 
                           {
                              Grav2.x+=PosN.PosAxe.x;
                              Grav2.y+=PosN.PosAxe.y;
                              NbOverlapB++;
                              if (PosN.PosAxe.x>CentreMaxB.x) CentreMaxB.x=PosN.PosAxe.x;
                              if (PosN.PosAxe.y>CentreMaxB.y) CentreMaxB.y=PosN.PosAxe.y;
                              if (PosN.PosAxe.x<CentreMinB.x) CentreMinB.x=PosN.PosAxe.x;
                              if (PosN.PosAxe.y<CentreMinB.y) CentreMinB.y=PosN.PosAxe.y;
                              if (ParcoursSens==0) // on prend la derniere largeur comme LargeurB
                              {
                                 LargeurB=PosN.Largeur;
                              }

                              if (ParcoursSens==1) // on prend la premiere largeur comme LargeurB
                              {
                                 if (numpt==0) { LargeurB=PosN.Largeur; numpt++;}
                              }
                              ListeAxeB.Add(PosN.PosAxe);
                           }
                           ListeB->GetNext(pos);// on passe à la section suivante
                        }
                        LargMoyB/=ListeB->GetCount();

                        if (NbOverlapB!=0)
                        {
                           Grav2.x=Grav2.x/NbOverlapB;
                           Grav2.y=Grav2.y/NbOverlapB;
                           CentreB.x=(CentreMinB.x+CentreMaxB.x)/2.0;
                           CentreB.y=(CentreMinB.y+CentreMaxB.y)/2.0;

                           //SumDis=0.0;
                           //for (n=0;n<NbOverlap;n++)
                           //{
                           //	// find the dist min
                           //	
                           //	MinDisPt=(double)__max(Larg,Haut);
                           //	for (k=0;k<__min(NbOverlapB,NbOverlap);k++)
                           //	{
                           //		DistPtAxe=(ListeAxeA[n].x-ListeAxeB[k].x)*(ListeAxeA[n].x-ListeAxeB[k].x)+
                           //			(ListeAxeA[n].y-ListeAxeB[k].y)*(ListeAxeA[n].y-ListeAxeB[k].y);
                           //		if (DistPtAxe<MinDisPt) MinDisPt=DistPtAxe;
                           //	}
                           //	SumDis+=MinDisPt;
                           //}
                           //SumDis/=(double)NbOverlap;

                           VarianceA=Variance(ListeAxeA,Grav1.x);
                           VarianceB=Variance(ListeAxeB,Grav2.x);
                           CovarianceA=Covariance(ListeAxeA,Grav1);
                           CovarianceB=Covariance(ListeAxeB,Grav2);
                           //						SumDis=fabs(VarianceA*CovarianceB-VarianceB*CovarianceA);
                           double AlphaA,AlphaB;
                           AlphaA=atan2(CovarianceA,VarianceA);
                           AlphaB=atan2(CovarianceB,VarianceB);
                           SumDis=fabs(AlphaA-AlphaB)*180.0/3.14159;

                           // if the distance between the gravity centers is low, merge the two lists

                           //DistGrav=(Grav1.x-Grav2.x)*(Grav1.x-Grav2.x)+(Grav1.y-Grav2.y)*(Grav1.y-Grav2.y);
                           DistGrav=(CentreA.x-CentreB.x)*(CentreA.x-CentreB.x)+(CentreA.y-CentreB.y)*(CentreA.y-CentreB.y);
                           DistGrav=sqrt(DistGrav);
                           //		double DistLargMoy=fabs(LargMoyB-LargMoyA)/log(__max(LargMoyA,LargMoyB));
                           double DistLargMoy=fabs(LargeurB-LargeurA)/log(__max(LargeurA,LargeurB));

#ifdef ANIM
                           u=Larg*(int)CentreA.y+(int)CentreA.x;
                           pImage->ImageRGB[u].r=255; pImage->ImageRGB[u].g=255; pImage->ImageRGB[u].b=255;
                           u=Larg*(int)CentreB.y+(int)CentreB.x;
                           pImage->ImageRGB[u].r=0; pImage->ImageRGB[u].g=255; pImage->ImageRGB[u].b=0;
                           CString Tempo;
                           Tempo.Format("DistLargMoy=%f\n distgrav=%f \n SumDis=%f",DistLargMoy,DistGrav,SumDis);
                           MessageBox(Tempo,"");
#endif
                           if (DistGrav<__max(LargMoyA,LargMoyB)/4 && DistLargMoy<SeuilLargMoy && SumDis<45.0 ) //  merge
                           {
#ifdef ANIM 
                              MessageBox("Fusion");
#endif
                              MergeOk=true;
                              if (ParcoursSens==0)
                              {
                                 int NbpAj=0;
                                 pos = ListeB->GetHeadPosition();
                                 while (pos)
                                 {
                                    PosN =  ListeB->GetAt(pos);
                                    u=Larg*(int)PosN.PosAxe.y+(int)PosN.PosAxe.x;
                                    if (Plan1[u]==0)  // if the axis is not already in A region
                                    {
                                       if (AddSens==1) ListeA->AddTail(PosN);
                                       if (AddSens==0)  ListeA->AddHead(PosN);
                                    }
                                    ListeB->GetNext(pos);
                                    NbpAj++;
                                 } // while pos
                                 /*CString Tempo;
                                 Tempo.Format("Nb Sens 0=%d",NbpAj);
                                 MessageBox(Tempo,"");*/
                                 //CStdioFile FichierTemp;
                                 //FichierTemp.Open("c:\\liste.txt",CFile::modeCreate|CFile::modeWrite);
                                 //CString Tempo;
                                 //pos = TableauListesPointsFiltreesC2D[IndiceReg]->GetHeadPosition();
                                 //while (pos)
                                 //{
                                 //	PosN =  TableauListesPointsFiltreesC2D[IndiceReg]->GetAt(pos);
                                 //	Tempo.Format("%d;%d;%d;%d\n",PosN.PosInit.x,PosN.PosInit.y,PosN.PosFinal.x,PosN.PosFinal.y);
                                 //	FichierTemp.WriteString(Tempo);
                                 //	TableauListesPointsFiltreesC2D[IndiceReg]->GetNext(pos);
                                 //}
                                 //			FichierTemp.Close();
                              }

                              if (ParcoursSens==1)
                              {
                                 int NbpAj=0;
                                 pos = ListeB->GetTailPosition();
                                 while (pos)
                                 {
                                    PosN =  TableauListesPointsFiltreesC2D[IndiceRegB]->GetAt(pos);
                                    u=Larg*(int)PosN.PosAxe.y+(int)PosN.PosAxe.x;
                                    if (Plan1[u]==0 )  // if the axis is not already in the A region
                                    {
                                       if (AddSens==1) ListeA->AddTail(PosN);
                                       if (AddSens==0) ListeA->AddHead(PosN);
                                    }
                                    ListeB->GetPrev(pos);
                                    NbpAj++;
                                 }
                                 /*CString Tempo;
                                 Tempo.Format("Nb sens 1=%d",NbpAj);
                                 MessageBox(Tempo,"");*/
                              }

                              // update the plan1 from merged regions
                              BuildImageFromList(Plan1,ListeA);
#ifdef ANIM
                              for (u=0;u<Taille;u++)
                              {
                                 if (Plan1[u]==1) pImage->ImageRGB[u].r=255; else pImage->ImageRGB[u].r=0; 
                                 if (Plan1[u]==1) pImage->ImageRGB[u].g=0; else pImage->ImageRGB[u].g=0;
                                 if (Plan1[u]==1) pImage->ImageRGB[u].b=0; else pImage->ImageRGB[u].b=0;
                              }
                              pImage->SendMessage(WM_PAINT);
                              pImage->Invalidate(TRUE);
                              CString Ttempo;
                              Ttempo.Format(" %d dans %d",IndiceRegB,IndiceReg);
                              MessageBox("liste rouge fusionnee"+Ttempo);
                              Sleep(2000);
#endif

                           } // end if merging is possible
                        } // if nboverlapp2
                        else 
                        {int bugggg=0; }
                     } // end if B is in A
                  } // end if no overlap
               } // end if no merging


               if (MergeOk)
               {
                  TabListeOk[IndiceRegB]=-2; // mark all the regions B as merged
                  // update the indices of all box connected regions
                  posB=TabNumListesConnex[IndiceRegB]->GetHeadPosition();
                  while(posB) // for all box-connected lists
                  {
                     IndiceRegConnected=TabNumListesConnex[IndiceRegB]->GetAt(posB);
                     posC=TabNumListesConnex[IndiceRegConnected]->GetHeadPosition();
                     while (posC)
                     {
                        IndiceRegConnectedToB=TabNumListesConnex[IndiceRegConnected]->GetAt(posC);
                        if (IndiceRegConnectedToB==IndiceRegB) TabNumListesConnex[IndiceRegConnected]->SetAt(posC,IndiceReg);
                        TabNumListesConnex[IndiceRegConnected]->GetNext(posC);
                     }
                     TabNumListesConnex[IndiceRegB]->GetNext(posB);
                  }
               }
            } // if list not already merged
            TabNumListesConnex[IndiceReg]->GetNext(posA);
         } // while pos
      } // if the region is not yet merged
   }


   free(Plan1);
   free(Plan2);
#ifdef ANIM
   for (u=0;u<Taille;u++)
   {
      pImage->ImageRGB[u].r=0;
      pImage->ImageRGB[u].g=0; 
      pImage->ImageRGB[u].b=0; 
   }
#endif
}



double CImageWnd::Variance(CArray <Sfpoint,Sfpoint> & Tab, double Moy)
{
   int NbPt=Tab.GetCount();
   double Var=0.0;
   for (int n=0;n<NbPt;n++)
   {
      Var+=(Tab[n].x-Moy)*(Tab[n].x-Moy);
   }
   Var/=(double)NbPt;
   return Var;
}

double CImageWnd::Variance(Sfpoint *pTab, long NbPt,double Moy)
{
   double Var=0.0;
   for (int n=0;n<NbPt;n++)
   {
      Var+=(pTab[n].x-Moy)*(pTab[n].x-Moy);
   }
   Var/=(double)NbPt;
   return Var;
}

Sfpoint CImageWnd::Projection(Sfpoint A, Sfpoint B,Sfpoint C,double *pr)
{
   // P is the projection of C onto AB
   Sfpoint P;
   double r;
   double sqnorme;
   sqnorme=(B.x-A.x)*(B.x-A.x)+(B.y-A.y)*(B.y-A.y);
   r=(A.y-C.y)*(A.y-B.y) - (A.x-C.x)*(B.x-A.x);
   r=r/sqnorme;
   P.x=A.x+r*(B.x-A.x);
   P.y=A.y+r*(B.y-A.y);

   *pr=r;
   /*
   r=0 if P=A
   r=1 if p=B
   r<0 if P is on backward extension of AB
   r>0 if P is on forward ext. of AB
   0<r<1 if P in interior to AB
   */
   return P;
}

double CImageWnd::Scalaire(Sfpoint A, Sfpoint B,Sfpoint C)
{
   // compute the position from A to the projection of C on AB 
   Sfpoint P;
   double r;
   double sqnorme;
   sqnorme=(B.x-A.x)*(B.x-A.x)+(B.y-A.y)*(B.y-A.y);
   r=(B.x-A.x)*(C.x-A.x)+(B.y-A.y)*(C.y-A.y);
   r=r/sqnorme;

   /*
   r=0 if P=A
   r=1 if p=B
   r<0 if P is on backward extension of AB
   r>0 if P is on forward ext. of AB
   0<r<1 if P in interior to AB
   */
   return r;
}


double CImageWnd::Covariance(Sfpoint *pTab, long NbPt, Sfpoint Moy)
{
   double Cov=0.0;
   for (int n=0;n<NbPt;n++)
   {
      Cov+=(pTab[n].x-Moy.x)*(pTab[n].y-Moy.y);
   }
   Cov/=(double)NbPt;
   return Cov;

}
double CImageWnd::Covariance(CArray <Sfpoint,Sfpoint> & Tab, Sfpoint Moy)
{
   int NbPt=Tab.GetCount();
   double Cov=0.0;
   for (int n=0;n<NbPt;n++)
   {
      Cov+=(Tab[n].x-Moy.x)*(Tab[n].y-Moy.y);
   }
   Cov/=(double)NbPt;
   return Cov;

}

void CImageWnd::PetrouCouleur(BYTE *pDest,rgb *pSource,int SeuilPetrou, int SigmaGaussien,ImagePetrou *I_x,ImagePetrou *I_y)
{
   ImagePetrou IR,IB,IG,B;
   rgb RGB;
   //int SeuilPetrou=3;
   //int SigmaGaussien=1;
   long l = 0;
   int i,j;

   BYTE   * ImR = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);
   BYTE   * ImV = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);
   BYTE   * ImB = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);
   BYTE   * ImageGauss = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);

   for (long j= 0;j<Haut;j++) 
   {                              
      for (long i=0;i<Larg;i++) 
      {       
         RGB = pSource[j*Larg+i];
         ((unsigned char   * )ImR)[l] = RGB.r; 
         ((unsigned char   * )ImV)[l] = RGB.g; 
         ((unsigned char   * )ImB)[l] = RGB.b; 
         l++;
      } 
   }	


   IR.dim_x = IG.dim_x=IB.dim_x=Larg;
   IR.dim_y = IG.dim_y =IB.dim_y =Haut;

   B.dim_x = Larg;
   B.dim_y = Haut;
   I_x->dim_x = Larg;
   I_x->dim_y = Haut;
   I_y->dim_x = Larg;
   I_y->dim_y = Haut;

   IR.ima = (Pixelp**) malloc(sizeof(Pixelp*)*IR.dim_y);
   IG.ima = (Pixelp**) malloc(sizeof(Pixelp*)*IG.dim_y);
   IB.ima = (Pixelp**) malloc(sizeof(Pixelp*)*IB.dim_y);

   B.ima = (Pixelp**) malloc(sizeof(Pixelp*)*B.dim_y);

   for(i=0;i<IR.dim_y;i++)
   {
      IR.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*IR.dim_x);
      IG.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*IG.dim_x);
      IB.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*IB.dim_x);
      B.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*B.dim_x);
   }

   Convolution2DCedricCurvi(ImageGauss,ImR,SigmaGaussien,Larg,Haut);
   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         IR.ima[i][j] = (Pixelp)ImageGauss[i*Larg+j]; // rouge
      }
   }
   Convolution2DCedricCurvi(ImageGauss,ImB,SigmaGaussien,Larg,Haut);
   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         IB.ima[i][j] = (Pixelp)ImageGauss[i*Larg+j]; // bleu
      }
   }
   Convolution2DCedricCurvi(ImageGauss,ImV,SigmaGaussien,Larg,Haut);
   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         IG.ima[i][j] = (Pixelp)ImageGauss[i*Larg+j]; // vert
      }
   }

   CPetrou ContourPetrou(&IR, &B, I_x,I_y);
   ContourPetrou.filtre_petrou_couleur(&IR,&IB,&IG,&B,6,SeuilPetrou,1);

   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         pDest[i*Larg+j] = (BYTE) B.ima[i][j];
      }
   }

   for(i=0;i<IR.dim_y;i++)
   {
      free(IR.ima[i]);
      free(IG.ima[i]);
      free(IB.ima[i]);
      free(B.ima[i]);

   }

   free(IR.ima);
   free(IG.ima);
   free(IB.ima);
   free(B.ima);


   free(ImageGauss);
   free(ImR);
   free(ImV);
   free(ImB);
}
#endif

void CImageWnd::OnTraitementPetroucouleur()
{
#ifdef MP
   // TODO : ajoutez ici le code de votre gestionnaire de commande
   CTraitImageWnd *pImagePetrou=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, NB);
   int SeuilPetrou=5;
   ImagePetrou I_x,I_y;
   I_x.dim_x=I_y.dim_x=Larg;
   I_x.dim_y=I_y.dim_y=Haut;
   int i;

   I_x.ima = (Pixelp**) malloc(sizeof(Pixelp*)*I_x.dim_y);
   I_y.ima = (Pixelp**) malloc(sizeof(Pixelp*)*I_y.dim_y);

   for(i=0;i<Haut;i++)
   {
      I_x.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*I_x.dim_x);
      I_y.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*I_y.dim_x);

   }

   PetrouCouleur(pImagePetrou->Image,ImageRGB,SeuilPetrou,1,&I_x,&I_y);
   pImagePetrou->Create(Titre + " : Petrou",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());

   for(i=0;i<Haut;i++)
   {
      free(I_x.ima[i]); free (I_y.ima[i]);
   }
   free(I_x.ima);
   free(I_y.ima);
#endif
}

void CImageWnd::OnTraitementCurvicouleur()
{
   if (BoiteReglageParam.DoModal()==IDCANCEL) return ;
   
   CCurvLinearParams p;

   p.SeuilReponse = BoiteReglageParam.SeuilReponse;
   p.LargeurEtude = BoiteReglageParam.ValLargeur;
   p.LargeurFourier = BoiteReglageParam.ValLargeurFourier;
   p.LongueurMiniForme = BoiteReglageParam.NbListeMini;
   p.LargeurMini = BoiteReglageParam.ValLargeurMini;
   p.VLocaleLargeur = BoiteReglageParam.m_EditVariationLocaleLargeur;

   BOOL DebugEgde = BoiteReglageParam.m_BDebugEdges;
   BOOL DebugTracking = BoiteReglageParam.m_BDebugTracking;

   SImage info = GetSImage();
   
   utls::BRGBAry *tmp = new utls::BRGBAry(info.Haut, info.Larg, (utls::RGBValue<unsigned char> *)this->ImageRGB);
   typedef utls::Ary<float, float, double> CurvImage;
   typedef utls::Ary<utls::RGBValue<float>, float, utls::RGBValue<double> > CurvRGBImage;

   CurvRGBImage *im = utls::convert<utls::BRGBAry, CurvRGBImage>(tmp);
   utls::BAry *edge;
   TableOfCurvilinearRegions TableauListesPointsFiltrees, TableauListesPointsFiltreesC2D;
   DetectCurvLinearPoints(im, p, 0, &edge, TableauListesPointsFiltrees, TableauListesPointsFiltreesC2D);

   tmp->detach(); delete tmp;


   VisualiseCurvi(DebugEgde, DebugTracking, info.Larg, info.Haut, edge, TableauListesPointsFiltrees, TableauListesPointsFiltreesC2D);
   //MessageBox("Export & affichage ok ","Export & affichage ok ",MB_OK);
   //----------------------------------------------------------------------------------------------------
   Invalidate(FALSE);
   SetCursor(LoadCursor(0,IDC_ARROW));

   delete im; delete edge;

#ifdef MP
   if (BoiteReglageParam.DoModal()==IDCANCEL) return ;
   bool bAff=TRUE; // for progress


   CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> TableauListesPointsFiltrees;
   CList <PosDetectCurvi,PosDetectCurvi> * ListePointsFiltres;
   CList <PosDetectCurvi,PosDetectCurvi> * ListePointsFiltresComp;
   CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> TableauListesPointsFiltreesC2D;
   CList <PosDetectCurvi,PosDetectCurvi> * ListePointsFiltresC2D;
   CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> TableauListesFinal;

   CArray <CList<int,int>*,CList<int,int>*> TabNumListesConnex;
   CArray <SBoundingBox,SBoundingBox> TableauListesBoundingBox;
   SImage info;

   CCurviDetect2D CurviDetect2D;
   CRgn * ConstrucRegion;
   CPoint * TabTraceRegion;
   int MaxiLongueurListe;
   int IndiceListeMaxi;
   int IndiceReg;

   float SigneDx = 0.0f;
   int Sens=-1;
   float OldSigneDx=0.0f;

   BOOL Init = FALSE;
   PosDetectCurvi PosNM1,PosN,PosNP1;
   double MoyPrev, MoyNext;
   POSITION PositionN, PositionNM1,PositionCListeCopy;
   int CptTableauListe = 0;
   Taille=Larg*Haut;

   CFenetreProgressCurvi * AffichageProgression;
   AffichageProgression = new CFenetreProgressCurvi(this);
   AffichageProgression->Create(IDD_DIALOG_PROGRESSCURVI);
   LARGE_INTEGER  FreqProc, TpsStart, TPetrou, TTracking, T1DProcess, T2DProcess, TAffichage;
   LARGE_INTEGER  tps1,tps2,tps3,tps4,tps5;
   QueryPerformanceFrequency(&FreqProc);

   m_SeuilReponse = BoiteReglageParam.SeuilReponse;
   m_LargeurEtude = BoiteReglageParam.ValLargeur;
   m_LargeurFourier = BoiteReglageParam.ValLargeurFourier;
   m_LongueurMiniForme = BoiteReglageParam.NbListeMini;
   m_LargeurMini = BoiteReglageParam.ValLargeurMini;
   VLocaleLargeur = BoiteReglageParam.m_EditVariationLocaleLargeur;
   BOOL DebugEgde = BoiteReglageParam.m_BDebugEdges;
   BOOL DebugTracking = BoiteReglageParam.m_BDebugTracking;
   int m_LongueurMiniRegionFiltrage = m_LongueurMiniForme/2;
   m_LongueurMiniForme=20; // debug
   long CompteurRegionFiltrage = 0;
   int SeuilPetrou=4;
   int SigmaGaussien=3;

   int SeuilDepass=5;
   int SeuilExtrapol=5;

   BOOL BoiteAxe=TRUE;


   int TabHisto[1000];

   CPoint p1,p2,p3;
   Taille=Larg*Haut;


   ImagePetrou APasGaussien/*Contient l'image original*/,B,I_x, I_y;
   int NbListe;
   int i=0,j=0,k=0;

   int w = m_LargeurEtude;
   int dw = w;
   int cpti = 0;
   int cptpoint = 0;
   int LargeurFourier = m_LargeurFourier;

   int NBPts = 2*(w+dw)+ 1;
   int Centre;
   POSITION pos;
   CPoint CI1, CI2, CI3;
   BOOL bCI1 = FALSE,bCI2 = FALSE,bCI3 = FALSE;
   float GradMoy_x = 0.0f, GradMoy_y = 0.0f;
   CPoint CurentPoint;

   /*BYTE * TabGauche = (BYTE *)malloc(sizeof(BYTE)*((NBPts/2)+2*LargeurFourier));
   BYTE * TabDroite = (BYTE *)malloc(sizeof(BYTE)*(((NBPts/2)+1+2*LargeurFourier)));*/

   rgb * TabGaucheCoul = (rgb *)malloc(sizeof(rgb)*(w+1+LargeurFourier));
   rgb * TabDroiteCoul = (rgb *)malloc(sizeof(rgb)*(w+1+LargeurFourier));

   CPoint * TabCoupe = (CPoint*) malloc(sizeof(CPoint)* NBPts);
   CPoint * TabCoupeGauche   = new CPoint[w+1+LargeurFourier];

   CPoint * TabCoupeEdge = (CPoint*) malloc(sizeof(CPoint)* NBPts);
   CPoint * TabCoupeGaucheEdge   = new CPoint[w+1+LargeurFourier];

   for(i=0;i<NBPts;i++)
   {
      TabCoupe[i].x=0;
      TabCoupe[i].y=0;
   }

   BYTE * ImageFreeman = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);
   BYTE * ImageGauss = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);
   int *ImageLabelEdge=(int*) malloc(sizeof(int)*Larg*Haut);

   long p=0;
   int cptCoul=0;
   int MaxListe = 0, PosMaxListe = 0;
   int CptPtsGrad = 0;
   BYTE r,g,b;

   int ValIptInitFixe;
   PosDetectCurvi PosGauche;
   PosDetectCurvi PosDroite;
   PosDetectCurvi PosARetenir;

   CString ChaineVisu;
   //CFile FichierVisu;
   CStdioFile FichierVisu;
   CImageWnd * pImage=NULL;	
   CImageWnd * pImagePetrou=NULL;
   CImageWnd * pImageTracking=NULL;

   CPoint TabPourQuadri[4];

   //--------------------------------------------------------------------
   //			construction object ccurvidetect1d
   //--------------------------------------------------------------------

   CCurviDetect1D CurviDetect1D(w+1+LargeurFourier,LargeurFourier,LargeurFourier);
   // ---------- informations relatives a l'image -----
   info=GetSImage();


   pImage=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);
   pImagePetrou=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);
   pImageTracking=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);
   int l,c,u;
   for (l=0;l<Haut;l++)
   {
      for (c=0;c<Larg;c++)
      {
         u=Larg*l+c;
         pImage->ImageRGB[u].r=0;
         pImage->ImageRGB[u].g=0;
         pImage->ImageRGB[u].b=0;
         ImageLabelEdge[u]=0;
      }
   }

   SetCursor(LoadCursor(0,IDC_WAIT));

   /* pImageAffichePoints=(CTraitImageWnd * ) new CImageWnd(Larg,Haut, COLOR);
   SetCursor(LoadCursor(0,IDC_WAIT));*/

   CMDIFrameWnd *pmd=GetMDIFrame();// Pour l'affichage du debugage

   pImage->Create(Titre+" Curvilinear Regions",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   pImage->SetActiveWindow();

   if(DebugEgde)
   {
      pImagePetrou->Create(Titre+" Edge detections",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   }

   if(DebugTracking)
   {
      pImageTracking->Create(Titre+" Edge detections",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   }


   APasGaussien.dim_x =I_x.dim_x=I_y.dim_x= Larg;
   APasGaussien.dim_y =I_x.dim_y=I_y.dim_y= Haut;
   APasGaussien.ima = (Pixelp**) malloc(sizeof(Pixelp*)*APasGaussien.dim_y);

   I_x.ima = (Pixelp**) malloc(sizeof(Pixelp*)*I_x.dim_y);
   I_y.ima = (Pixelp**) malloc(sizeof(Pixelp*)*I_y.dim_y);

   for(i=0;i<Haut;i++)
   {
      APasGaussien.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*APasGaussien.dim_x);
      I_x.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*I_x.dim_x);
      I_y.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*I_y.dim_x);

   }
   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         APasGaussien.ima[i][j] = (Pixelp)(ImageRGB[i*Larg+j].r)/3; // attention couleur
         APasGaussien.ima[i][j]+= (Pixelp)(ImageRGB[i*Larg+j].g)/3; // attention couleur
         APasGaussien.ima[i][j]+= (Pixelp)(ImageRGB[i*Larg+j].b)/3; // attention couleur
      }

   }




   //--------------------------------------------------------
   //			Contour Petrou
   //--------------------------------------------------------
   //AffichageProgression->UpdateData(TRUE);
   AffichageProgression->m_TPetrou = 0.0f;
   AffichageProgression->m_TCTracking = 0.0f;
   AffichageProgression->m_T1DProcess = 0.0f;
   AffichageProgression->m_T2DProcess = 0.0f;
   AffichageProgression->m_Affichage = 0.0f;
   AffichageProgression->m_TexteValProgress = "0 %";
   AffichageProgression->m_ControlBarProgressCurvi.SetPos(0);
   AffichageProgression->m_TxtNbReg="";
   AffichageProgression->UpdateData(FALSE);

   AffichageProgression->ShowWindow(bAff); // temp

   //CPetrou ContourPetrou(NULL, &B, &I_x,&I_y);

   QueryPerformanceCounter(&TpsStart); 

   PetrouCouleur(ImageFreeman,ImageRGB,3,1,&I_x,&I_y);

   QueryPerformanceCounter(&TPetrou);
   tps1.QuadPart = TPetrou.QuadPart-TpsStart.QuadPart;



   AffichageProgression->m_TPetrou = (double)tps1.QuadPart/(double)FreqProc.QuadPart;
   AffichageProgression->m_TCTracking = 0.0f;
   AffichageProgression->m_T1DProcess = 0.0f;
   AffichageProgression->m_T2DProcess = 0.0f;
   AffichageProgression->m_Affichage = 0.0f;
   AffichageProgression->m_TexteValProgress = "20 %";
   AffichageProgression->m_ControlBarProgressCurvi.SetPos(20);
   AffichageProgression->UpdateData(FALSE);
   AffichageProgression->ShowWindow(bAff); //  FALSE temporaire
   //--------------------------------------------------------
   //			Recopie Image Contour Petrou pour suivi
   //--------------------------------------------------------

   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         if(DebugEgde)
         {
            pImagePetrou->ImageRGB[i*Larg+j].r = ImageFreeman[i*Larg+j];
            pImagePetrou->ImageRGB[i*Larg+j].g = ImageFreeman[i*Larg+j];
            pImagePetrou->ImageRGB[i*Larg+j].b = ImageFreeman[i*Larg+j];
         }
         //pImageAffichePoints->ImageRGB[i*Larg+j].r = Image[i*Larg+j];//0;
         //pImageAffichePoints->ImageRGB[i*Larg+j].g = Image[i*Larg+j];//0;
         //pImageAffichePoints->ImageRGB[i*Larg+j].b = Image[i*Larg+j];//0;
      }

   }

   //--------------------------------------------------------
   //			suivi de contours
   //--------------------------------------------------------
   //MessageBox("avant appel suivi freeman");
   CurviDetect2D.FermetureContourSimple(ImageFreeman,info.Larg,info.Haut);
   if(DebugEgde)
   {
      for(i=0;i<Haut;i++)
      {
         for(j=0;j<Larg;j++)
         {
            pImagePetrou->ImageRGB[i*Larg+j].r = (BYTE) ImageFreeman[i*Larg+j];
            pImagePetrou->ImageRGB[i*Larg+j].g = (BYTE) ImageFreeman[i*Larg+j];
            pImagePetrou->ImageRGB[i*Larg+j].b = (BYTE) ImageFreeman[i*Larg+j];
         }
      }
   }

   CurviDetect2D.SuiviFreeman(ImageFreeman,info.Larg,info.Haut,pmd);

   QueryPerformanceCounter(&TTracking);

   TableauListesPointsFiltrees.RemoveAll();
   TableauListesPointsFiltreesC2D.RemoveAll();

   tps2.QuadPart = TTracking.QuadPart-TPetrou.QuadPart;

   AffichageProgression->m_TCTracking = (double)tps2.QuadPart/(double)FreqProc.QuadPart;
   AffichageProgression->m_T1DProcess = 0.0f;
   AffichageProgression->m_T2DProcess = 0.0f;
   AffichageProgression->m_Affichage = 0.0f;
   AffichageProgression->m_TexteValProgress = "40 %";
   AffichageProgression->m_ControlBarProgressCurvi.SetPos(40);
   AffichageProgression->UpdateData(FALSE);
   AffichageProgression->Invalidate(TRUE);
   AffichageProgression->ShowWindow(bAff); // temp

   //::MessageBox(NULL,"STOP","Après suivi freeman",MB_OK);
   ////---------------------------------------------------------
   ////				Affichage couleur fenetre résultats
   ////---------------------------------------------------------

   int cptlistetemp=0;

   //pImageAffichePoints->Create("Points Détecteur curvi",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
   int OldIptFinalGauche=0;
   int OldIptFinalDroite=0;
   int CptExtrapolGauche,CptExtrapolDroite;



   NbListe = CurviDetect2D.TabListePoint.GetCount();

   if (NbListe !=0)
   {
      for(i=0;i<NbListe;i++)
      {
         pos = CurviDetect2D.TabListePoint[i]->GetHeadPosition();
         while(pos)
         {
            CurentPoint = CurviDetect2D.TabListePoint[i]->GetAt(pos);
            CurviDetect2D.TabListePoint[i]->GetNext(pos);
            u=Larg*CurentPoint.y+CurentPoint.x;
            ImageLabelEdge[u]=i+1;

            //if(DebugTracking)
            //{
            //	pImageTracking->ImageRGB[u].r = 255-ImageLabelEdge[u];
            //	pImageTracking->ImageRGB[u].g = 255-ImageLabelEdge[u];
            //	pImageTracking->ImageRGB[u].b = 255-ImageLabelEdge[u];
            //}
         }
      }
   }

   if (NbListe !=0)
   {
      for(i=0;i<NbListe;i++)
      {

         if(CurviDetect2D.TabListePoint[i]->GetCount()>m_LongueurMiniForme)
         {
            ListePointsFiltres = new CList<PosDetectCurvi,PosDetectCurvi>;// creation pour les listes dépassants la contraintes de largeurs mini.
            ListePointsFiltresComp = new CList<PosDetectCurvi,PosDetectCurvi>;

            double LastWidthGauche = -1.0f,LastWidthDroite = -1.0f;
            cptCoul++;
            pos = CurviDetect2D.TabListePoint[i]->GetHeadPosition();
            CptExtrapolGauche=0;
            CptExtrapolDroite=0;

            //tirage au sort de la couleur pour la forme entière

            r= (rand()*235)+20;
            g= (rand()*235)+20;
            b= (rand()*235)+20;

            CurviDetect1D.MaxiCurvi.IptFinal=-1;
            OldIptFinalGauche=-1;
            OldIptFinalDroite=-1;

            int longueur;
            CPoint PointSuiv;
            while(pos)
            {
               CurentPoint = CurviDetect2D.TabListePoint[i]->GetAt(pos);
               /*			if(pos)CurviDetect2D.TabListePoint[i]->GetNext(pos);
               if(pos)CurviDetect2D.TabListePoint[i]->GetNext(pos);
               if(pos)CurviDetect2D.TabListePoint[i]->GetNext(pos);
               if(pos)CurviDetect2D.TabListePoint[i]->GetNext(pos);
               if(pos)CurviDetect2D.TabListePoint[i]->GetNext(pos);
               if(pos)CurviDetect2D.TabListePoint[i]->GetNext(pos);
               if (pos) PointSuiv=CurviDetect2D.TabListePoint[i]->GetAt(pos);
               else break;
               CurviDetect2D.TabListePoint[i]->GetPrev(pos);
               if(pos)CurviDetect2D.TabListePoint[i]->GetPrev(pos);
               if(pos)CurviDetect2D.TabListePoint[i]->GetPrev(pos);
               if(pos)CurviDetect2D.TabListePoint[i]->GetPrev(pos);
               if(pos)CurviDetect2D.TabListePoint[i]->GetPrev(pos);
               if(pos)CurviDetect2D.TabListePoint[i]->GetPrev(pos);*/

               if (CurentPoint.x>LargeurFourier && CurentPoint.x<Larg-LargeurFourier
                  && CurentPoint.y>LargeurFourier && CurentPoint.y<Haut-LargeurFourier)
               {

                  GradMoy_x = I_x.ima[CurentPoint.y][CurentPoint.x];
                  GradMoy_y = I_y.ima[CurentPoint.y][CurentPoint.x];

                  /*GradMoy_x=-(CurentPoint.y-PointSuiv.y);
                  GradMoy_y=(CurentPoint.x-PointSuiv.x);
                  if (GradMoy_x==0 && GradMoy_y==0) break;*/

                  //if(DebugTracking)
                  //{
                  //	pImageTracking->ImageRGB[CurentPoint.y*Larg+CurentPoint.x].r = r;
                  //	pImageTracking->ImageRGB[CurentPoint.y*Larg+CurentPoint.x].g = g;
                  //	pImageTracking->ImageRGB[CurentPoint.y*Larg+CurentPoint.x].b = b;
                  //}

                  if(GradMoy_x>0)
                     SigneDx = -Sens; 
                  else
                     SigneDx=+Sens;

                  if(GradMoy_x==0 && GradMoy_y>=0.0) SigneDx=-Sens;
                  if(GradMoy_x==0 && GradMoy_y<0.0) SigneDx=+Sens;

                  CurviDetect2D.RemplitTableauIJCoupeBiDir(TabCoupe,CurentPoint.x,CurentPoint.y,GradMoy_x,GradMoy_y,w,dw,SigneDx,&Centre,Haut,Larg);

                  /* affichage des IptInit*/
                  //p=TabCoupe[Centre].x+Larg*TabCoupe[Centre].y;
                  //pImageTracking->ImageRGB[p].r =  0; 
                  //pImageTracking->ImageRGB[p].g =  255;//(j+20)*2;
                  //pImageTracking->ImageRGB[p].b =  0;

                  //----------------------------------------------------------------------
                  //			Creation tableau gauche et droite
                  //----------------------------------------------------------------------
                  for(j=Centre-w;j<Centre+1+LargeurFourier;j++)//Gauche
                  {
                     p=TabCoupe[j].x+Larg*TabCoupe[j].y;
                     TabGaucheCoul[Centre-j+LargeurFourier].r = ImageRGB[p].r;
                     TabGaucheCoul[Centre-j+LargeurFourier].g = ImageRGB[p].g;
                     TabGaucheCoul[Centre-j+LargeurFourier].b = ImageRGB[p].b;
                     //APasGaussien.ima[TabCoupe[j].y][TabCoupe[j].x];
                     TabCoupeGauche[Centre+LargeurFourier-j] = TabCoupe[j];
                     //// papillons
                     //if (CurentPoint.x>200 && CurentPoint.y >105 && CurentPoint.x<283 && CurentPoint.y<150)
                     //{
                     //pImageTracking->ImageRGB[p].r =  0; 
                     //pImageTracking->ImageRGB[p].g =  (j+20)*2;
                     //pImageTracking->ImageRGB[p].b =  0;
                     //}

                  }

                  //CStdioFile FichierTemp;
                  //FichierTemp.Open("c:\\coupegauche.txt",CFile::modeCreate|CFile::modeWrite);
                  //CString Tempo;
                  //for (j=0;j<=w+LargeurFourier;j++ )
                  //{
                  //	Tempo.Format("%d\n",TabGauche[j]);
                  //	FichierTemp.WriteString(Tempo);
                  //}FichierTemp.Close();


                  cpti=0;
                  for(j=Centre-LargeurFourier;j<Centre+w+1;j++)//droite
                  {
                     p=TabCoupe[j].x+Larg*TabCoupe[j].y;
                     TabDroiteCoul[cpti].r = ImageRGB[p].r;//APasGaussien.ima[TabCoupe[j].y][TabCoupe[j].x];
                     TabDroiteCoul[cpti].g = ImageRGB[p].g;
                     TabDroiteCoul[cpti].b = ImageRGB[p].b;
                     cpti++;
                     /*p=TabCoupe[j].x+Larg*TabCoupe[j].y;
                     // papillons

                     pImagePetrou->ImageRGB[p].r =  (j+20)*2 ; 
                     pImagePetrou->ImageRGB[p].g =  0;
                     pImagePetrou->ImageRGB[p].b =  0;*/

                  }//fin for j

                  //FichierTemp.Open("c:\\coupedroite.txt",CFile::modeCreate|CFile::modeWrite);
                  //for (j=0;j<=w+LargeurFourier;j++ )
                  //{
                  //	Tempo.Format("%d\n",TabDroite[j]);
                  //	FichierTemp.WriteString(Tempo);
                  //}	
                  //FichierTemp.Close();


                  ValIptInitFixe = LargeurFourier;

                  CurviDetect1D.RespCurviMaxiFastCoul(TabGaucheCoul, TabCoupeGauche, LastWidthGauche,LargeurFourier,ImageLabelEdge,Larg);// modif

                  PosGauche.PosInit = TabCoupe[Centre+LargeurFourier-CurviDetect1D.MaxiCurvi.IptInit]; 
                  PosGauche.PosFinal = TabCoupe[Centre+LargeurFourier-CurviDetect1D.MaxiCurvi.IptFinal];
                  PosGauche.PositionMaxi = TabCoupe[Centre+LargeurFourier-CurviDetect1D.MaxiCurvi.Centre];
                  PosGauche.ValDist = CurviDetect1D.MaxiCurvi.ValDist;
                  PosGauche.Largeur = LastWidthGauche = CurviDetect1D.MaxiCurvi.Largeur;
                  PosGauche.PosAxe.x=((double)PosGauche.PosInit.x+(double)PosGauche.PosFinal.x)/2;
                  PosGauche.PosAxe.y=((double)PosGauche.PosInit.y+(double)PosGauche.PosFinal.y)/2;

                  if ((PosGauche.ValDist>=m_SeuilReponse) &&(PosGauche.Largeur>=m_LargeurMini))
                  {
                     ListePointsFiltres->AddTail(PosGauche);// ajout à la liste des candidats gardés
                     CptExtrapolGauche=0;
                  }
                  else
                  {
                     CptExtrapolGauche++;
                     if (CptExtrapolGauche>SeuilExtrapol)
                     {
                        longueur=ListePointsFiltres->GetCount();
                        if (longueur>m_LongueurMiniForme) 
                        { 
                           TableauListesPointsFiltrees.Add(ListePointsFiltres);
                           cptlistetemp++;
                        }
                        else delete ListePointsFiltres;
                        ListePointsFiltres = new CList<PosDetectCurvi,PosDetectCurvi>;// creation pour les listes dépassants la contraintes de largeurs mini.
                        CptExtrapolGauche=0;
                     }
                  }
                  //else
                  //{
                  //	if (OldIptFinalGauche!=-1)
                  //	{
                  //		CptExtrapolGauche++;
                  //		if (CptExtrapolGauche<SeuilExtrapol)
                  //		{
                  //			PosGauche.PosInit=CurentPoint;
                  //			if (OldSigneDx*SigneDx<0) PosGauche.PosFinal=TabCoupe[w+LargeurFourier-OldIptFinalGauche];
                  //			if (GradMoy_x==0 && GradMoy_y>0.0) PosGauche.PosFinal=TabCoupe[w+1-LargeurFourier+OldIptFinalGauche];
                  //			else
                  //			{
                  //			if (GradMoy_x==0 && GradMoy_y>=0.0) PosGauche.PosFinal=TabCoupe[w+1-LargeurFourier+OldIptFinalGauche];
                  //			else PosGauche.PosFinal=TabCoupe[w+LargeurFourier-OldIptFinalGauche];
                  //			}
                  //				PosGauche.Largeur=sqrt(
                  //			( (double)PosGauche.PosInit.x-(double)PosGauche.PosFinal.x )*
                  //			( (double)PosGauche.PosInit.x-(double)PosGauche.PosFinal.x )+
                  //			( (double)PosGauche.PosInit.y-(double)PosGauche.PosFinal.y )*
                  //			( (double)PosGauche.PosInit.y-(double)PosGauche.PosFinal.y ) );
                  //			PosGauche.PosAxe.x=((double)PosGauche.PosInit.x+(double)PosGauche.PosFinal.x)/2;
                  //			PosGauche.PosAxe.y=((double)PosGauche.PosInit.y+(double)PosGauche.PosFinal.y)/2;
                  //			ListePointsFiltres->AddTail(PosGauche);// ajout à la liste des candidats gardés
                  //		}

                  //	}
                  //}
                  OldIptFinalGauche=CurviDetect1D.MaxiCurvi.IptFinal;

                  CurviDetect1D.RespCurviMaxiFastCoul(TabDroiteCoul, TabCoupe+Centre+1-LargeurFourier, LastWidthDroite,LargeurFourier,ImageLabelEdge,Larg);

                  PosDroite.PosInit = TabCoupe[Centre+1-LargeurFourier+CurviDetect1D.MaxiCurvi.IptInit];
                  PosDroite.PosFinal = TabCoupe[Centre+1-LargeurFourier+CurviDetect1D.MaxiCurvi.IptFinal];
                  PosDroite.PositionMaxi = TabCoupe[Centre+1-LargeurFourier+CurviDetect1D.MaxiCurvi.Centre];
                  PosDroite.ValDist = CurviDetect1D.MaxiCurvi.ValDist;
                  PosDroite.Largeur = LastWidthDroite = CurviDetect1D.MaxiCurvi.Largeur;
                  PosDroite.PosAxe.x=((double)PosDroite.PosInit.x+(double)PosDroite.PosFinal.x)/2;
                  PosDroite.PosAxe.y=((double)PosDroite.PosInit.y+(double)PosDroite.PosFinal.y)/2;

                  if ((PosDroite.ValDist>=m_SeuilReponse) && (PosDroite.Largeur>=m_LargeurMini) )
                  {
                     ListePointsFiltresComp->AddTail(PosDroite);// ajout à la liste des candidats gardés
                     CptExtrapolDroite=0;	
                  }
                  else
                  {
                     CptExtrapolDroite++;
                     if (CptExtrapolDroite>SeuilExtrapol)
                     {
                        longueur=ListePointsFiltresComp->GetCount();
                        if (longueur>m_LongueurMiniForme)
                        {
                           TableauListesPointsFiltrees.Add(ListePointsFiltresComp);
                           cptlistetemp++;
                        }
                        else delete ListePointsFiltresComp;
                        ListePointsFiltresComp = new CList<PosDetectCurvi,PosDetectCurvi>;
                        CptExtrapolDroite=0;
                     }

                  }

                  //else
                  //{
                  //	if (OldIptFinalDroite!=-1) 
                  //	{
                  //		CptExtrapolDroite++;
                  //		if (CptExtrapolDroite<SeuilExtrapol)
                  //		{

                  //			PosDroite.PosInit=CurentPoint;
                  //			if (OldSigneDx*SigneDx<0) PosDroite.PosFinal=TabCoupe[w+1-LargeurFourier+OldIptFinalDroite];
                  //			if (GradMoy_x==0 && GradMoy_y>0.0) PosDroite.PosFinal=TabCoupe[w+LargeurFourier-OldIptFinalDroite];
                  //			else
                  //			{
                  //			if (GradMoy_x==0 && GradMoy_y>=0.0) PosDroite.PosFinal=TabCoupe[w+1-LargeurFourier+OldIptFinalDroite];
                  //			else PosDroite.PosFinal=TabCoupe[w+LargeurFourier-OldIptFinalDroite];
                  //			}
                  //			PosDroite.Largeur=sqrt(
                  //			( (double)PosDroite.PosInit.x-(double)PosDroite.PosFinal.x )*
                  //			( (double)PosDroite.PosInit.x-(double)PosDroite.PosFinal.x )+
                  //			( (double)PosDroite.PosInit.y-(double)PosDroite.PosFinal.y )*
                  //			( (double)PosDroite.PosInit.y-(double)PosDroite.PosFinal.y ) );
                  //			ListePointsFiltres->AddTail(PosDroite);// ajout à la liste des candidats gardés
                  //			PosDroite.PosAxe.x=((double)PosDroite.PosInit.x+(double)PosDroite.PosFinal.x)/2;
                  //			PosDroite.PosAxe.y=((double)PosDroite.PosInit.y+(double)PosDroite.PosFinal.y)/2;
                  //		}
                  //	}
                  //}
               } // si on est dans l'image moins LargeurFourier
               OldIptFinalDroite=CurviDetect1D.MaxiCurvi.IptFinal;
               OldSigneDx=SigneDx;
               if(pos)CurviDetect2D.TabListePoint[PosMaxListe]->GetNext(pos);
               CptPtsGrad++;
               cptpoint++;
               /*}*/


            }// while(pos)

            // correction Jo
            longueur=ListePointsFiltres->GetCount();
            if (longueur>m_LongueurMiniForme)
            {
               TableauListesPointsFiltrees.Add(ListePointsFiltres);
               cptlistetemp++;

            } 
            else delete ListePointsFiltres;
            longueur=ListePointsFiltresComp->GetCount();
            if (longueur>m_LongueurMiniForme)
            {
               TableauListesPointsFiltrees.Add(ListePointsFiltresComp);
               cptlistetemp++;
            }
            else delete ListePointsFiltresComp;

         }//fin if nb points liste > 10
         //ajout de la liste en cours au tableau de liste de structure poscurvidetect

      }//for i
   }//fin if liste !=0 

   QueryPerformanceCounter(&T1DProcess);
   tps3.QuadPart = T1DProcess.QuadPart-TTracking.QuadPart;
   AffichageProgression->m_T1DProcess = (double)tps3.QuadPart/(double)FreqProc.QuadPart;
   AffichageProgression->m_T2DProcess = 0.0f;
   AffichageProgression->m_Affichage = 0.0f;
   AffichageProgression->m_TexteValProgress = "60 %";
   AffichageProgression->m_ControlBarProgressCurvi.SetPos(60);
   AffichageProgression->UpdateData(FALSE);
   AffichageProgression->Invalidate(TRUE);
   AffichageProgression->ShowWindow(bAff); // temp

   //MessageBox("1d termine");



   SectionOptimization(pImageTracking,TableauListesPointsFiltrees);

   int *Plan1;
   int *Plan2;
   Plan1=(int*) malloc(Taille*sizeof(int));
   Plan2=(int*) malloc(Taille*sizeof(int));
   if(DebugTracking)
   {
      int LngTableauListesPointsFiltrees=TableauListesPointsFiltrees.GetCount();
      for(IndiceReg=0;IndiceReg<LngTableauListesPointsFiltrees;IndiceReg++)//pour toutes les régions
      {
         BuildImageFromList(Plan1,TableauListesPointsFiltrees[IndiceReg]);
         TireCouleurAlea(&r,&g,&b);
         for (u=0;u<Taille;u++)
         {
            if (Plan1[u]==1) pImageTracking->ImageRGB[u].r=r;
            if (Plan1[u]==1) pImageTracking->ImageRGB[u].g=g;
            if (Plan1[u]==1) pImageTracking->ImageRGB[u].b=b;
         }
      }
   }
#define FAIRELASUITE
#ifdef FAIRELASUITE
   // no filter on the length at the 2D stage
   TwoDimConstraints(TableauListesPointsFiltreesC2D,TableauListesPointsFiltrees, 10 /*m_LongueurMiniRegionFiltrage*/, SeuilDepass);

   //MessageBox("2dd termine");


   //free(TabPointTemp);
   // filtrage des carrés et cercles
   //Filter of circles and squared shapes

   PosDetectCurvi PointCourant;
   Sfpoint AxeMin, AxeMax, AxeDif;
   SBoundingBox Box;
   double LargMoy,DifMax;
   double LargMin,LargMax,EtendueLarg;
   int longueur;
   int LongueurTableauListesPointsFiltreesC2D=TableauListesPointsFiltreesC2D.GetCount();
   if (BoiteAxe)
   {
      CList<int,int> *ListeNumConnex;
      int NumConnex;

      for(IndiceReg=0;IndiceReg<LongueurTableauListesPointsFiltreesC2D;IndiceReg++)//pour toutes les régions
      {
         LargMoy=0;
         pos = TableauListesPointsFiltreesC2D[IndiceReg]->GetHeadPosition();
         longueur=TableauListesPointsFiltreesC2D[IndiceReg]->GetCount();
         if(longueur>0 /*m_LongueurMiniRegionFiltrage*/)// si la liste n'est pas vide
         {
            // dermination liste englobante
            PointCourant = TableauListesPointsFiltreesC2D[IndiceReg]->GetAt(pos); // on récupère ce qui sera n-1
            TableauListesPointsFiltreesC2D[IndiceReg]->GetNext(pos);// on pointe sur le suivant
            LargMoy+=PointCourant.Largeur;
            AxeMin.x=AxeMax.x=PointCourant.PosAxe.x;
            AxeMin.y=AxeMax.y=PointCourant.PosAxe.y;
            LargMin=LargMax=PointCourant.Largeur;
            for (k=1;k<longueur;k++)
            {
               PointCourant = TableauListesPointsFiltreesC2D[IndiceReg]->GetAt(pos); // on récupère ce qui sera n-1
               if (PointCourant.PosAxe.x<AxeMin.x) AxeMin.x=PointCourant.PosAxe.x;
               if (PointCourant.PosAxe.y<AxeMin.y) AxeMin.y=PointCourant.PosAxe.y;
               if (PointCourant.PosAxe.x>AxeMax.x) AxeMax.x=PointCourant.PosAxe.x;
               if (PointCourant.PosAxe.y>AxeMax.y) AxeMax.y=PointCourant.PosAxe.y;
               LargMoy+=PointCourant.Largeur;
               if (PointCourant.Largeur>LargMax) LargMax=PointCourant.Largeur;
               if (PointCourant.Largeur<LargMin) LargMin=PointCourant.Largeur;
               TableauListesPointsFiltreesC2D[IndiceReg]->GetNext(pos);// on pointe sur le suivant
            }
            LargMoy/=longueur;
            AxeDif.x=fabs(AxeMin.x-AxeMax.x);
            AxeDif.y=fabs(AxeMin.y-AxeMax.y);
            DifMax=__max(AxeDif.x,AxeDif.y);
            EtendueLarg=fabs(LargMax-LargMin);
            if (DifMax<LargMoy*1.5/* || EtendueLarg>LargMoy*/) 
            {
               TableauListesPointsFiltreesC2D.RemoveAt(IndiceReg);
               IndiceReg--;
               LongueurTableauListesPointsFiltreesC2D--;
            }
            else
            {
               Box.AxeMin=AxeMin;
               Box.AxeMax=AxeMax;
               TableauListesBoundingBox.Add(Box);
               ListeNumConnex =new CList<int,int>;
               TabNumListesConnex.Add(ListeNumConnex);
            }
         } // si longueur ok

      } // for indice
   } // si boite Axe

   LongueurTableauListesPointsFiltreesC2D=TableauListesPointsFiltreesC2D.GetCount();
   int NumConnex;
   SBoundingBox BoxA,BoxB;
   int longlis=TabNumListesConnex.GetCount();

   for(IndiceReg=0;IndiceReg<LongueurTableauListesPointsFiltreesC2D;IndiceReg++)//pour toutes les régions
   {
      BoxA=TableauListesBoundingBox[IndiceReg];
      for( int IndiceRegB=0;IndiceRegB<LongueurTableauListesPointsFiltreesC2D;IndiceRegB++)//pour toutes les régions
      {
         if (IndiceReg!=IndiceRegB) 
         {
            BoxB=TableauListesBoundingBox[IndiceRegB];
            int Ext=0;
            if (BoxA.AxeMin.x>BoxB.AxeMax.x || BoxA.AxeMax.x<BoxB.AxeMin.x) Ext=1;
            if (BoxA.AxeMin.y>BoxB.AxeMax.y || BoxA.AxeMax.y<BoxB.AxeMin.y) Ext=1;
            if (Ext==0) // there is an intersection between boxes
            {
               // add the num of the list to the liste of candidates
               TabNumListesConnex[IndiceReg]->AddHead(IndiceRegB);
            }
         }
      }
   }

   LongueurTableauListesPointsFiltreesC2D=TableauListesPointsFiltreesC2D.GetCount();
   int *TabIndice;
   int *TabListeOk;
   TabIndice=(int *) malloc (LongueurTableauListesPointsFiltreesC2D*sizeof(int));
   TabListeOk=(int *) malloc (LongueurTableauListesPointsFiltreesC2D*sizeof(int));


   for(IndiceReg=0;IndiceReg<LongueurTableauListesPointsFiltreesC2D;IndiceReg++)//pour toutes les régions
   {
      TabListeOk[IndiceReg]=1; // mark all the regions as OK
      TabIndice[IndiceReg]=IndiceReg;
   }

   // fast merge

   TryToMerge(TableauListesPointsFiltreesC2D,TabNumListesConnex,TabListeOk,pImage);

   //MergeConnectedRegions(TableauListesPointsFiltreesC2D,pImage);

   QueryPerformanceCounter(&T2DProcess);
   tps4.QuadPart = T2DProcess.QuadPart-T1DProcess.QuadPart;
   AffichageProgression->m_T2DProcess = (double)tps4.QuadPart/(double)FreqProc.QuadPart;
   AffichageProgression->m_Affichage = 0.0f;
   AffichageProgression->m_TexteValProgress = "80 %";
   AffichageProgression->m_ControlBarProgressCurvi.SetPos(80);
   AffichageProgression->UpdateData(FALSE);
   AffichageProgression->Invalidate(TRUE);
   AffichageProgression->ShowWindow(bAff); // temp


   // remove too short regions

   for(IndiceReg=0;IndiceReg<LongueurTableauListesPointsFiltreesC2D;IndiceReg++)//pour toutes les régions
   {
      longueur=TableauListesPointsFiltreesC2D[IndiceReg]->GetCount();
      if(longueur<m_LongueurMiniRegionFiltrage && TabListeOk[IndiceReg]==1)// si la liste n'est pas vide
      {	
         TabListeOk[IndiceReg]=-4;
      }
   }



   //-----------------------------------------------------------------------------------
   // remove multiple regions
   //-----------------------------------------------------------------------------------
#define REMOVEMULTIPLECOUL
#ifdef REMOVEMULTIPLECOUL

   //int LongueurTableauListesPointsFiltreesC2D;
   int cptaff = 1;

   int Intersect=0,Surf1=0,Surf2=0,SurfMin;
   double Recouv;
   double SeuilRecou=0.7;
   int NbReg1=1;
   //CIma Test;
   //Test.LitBMP("f:\\johel\\prog\\bmp\\temp.bmp");
   IncrementeCouleur(&r,&g,&b,256);
   for (u=0;u<Taille;u++) {Plan1[u]=0;}	

   LongueurTableauListesPointsFiltreesC2D=TableauListesPointsFiltreesC2D.GetCount();
   // find the first non merged region
   int IndiceRegInit=0;
   while(TabListeOk[IndiceRegInit]!=1 && IndiceRegInit<LongueurTableauListesPointsFiltreesC2D)
   {
      IndiceRegInit++;
   } 
   if (IndiceRegInit!=LongueurTableauListesPointsFiltreesC2D)
   {
      BuildImageFromList(Plan1,TableauListesPointsFiltreesC2D[IndiceRegInit]);

      //for (u=0;u<Taille;u++)
      //{
      //if (Plan1[u]!=0) pImage->ImageRGB[u].r=255; else pImage->ImageRGB[u].r=0;
      //if (Plan1[u]!=0) pImage->ImageRGB[u].g=255;else pImage->ImageRGB[u].g=0;
      //if (Plan1[u]!=0) pImage->ImageRGB[u].b=255;else pImage->ImageRGB[u].b=0;
      //}
      //pImage->SendMessage(WM_PAINT);
      //pImage->Invalidate(TRUE);
      //Sleep(5000);
      for(IndiceReg=IndiceRegInit+1;IndiceReg<LongueurTableauListesPointsFiltreesC2D;IndiceReg++)//pour toutes les régions
      {
         longueur=TableauListesPointsFiltreesC2D[IndiceReg]->GetCount();
         if(TabListeOk[IndiceReg]==1)// si la liste n'est pas vide
         {	
            BuildImageFromList(Plan2,TableauListesPointsFiltreesC2D[IndiceReg]);
            /*for (u=0;u<Taille;u++)
            {
            if (Plan2[u]!=0) pImage->ImageRGB[u].r=255; else pImage->ImageRGB[u].r=0;
            if (Plan2[u]!=0) pImage->ImageRGB[u].g=255;else pImage->ImageRGB[u].g=0;
            if (Plan2[u]!=0) pImage->ImageRGB[u].b=255;else pImage->ImageRGB[u].b=0;
            }
            pImage->SendMessage(WM_PAINT);
            pImage->Invalidate(TRUE);
            Sleep(5000);*/
            // tester recouvrement p1 et p2, et recopier dans p1
            int IndiceMax=0;
            int Sens;
            //		Sens=CurviDetect2D.MesureRecouv(Plan1,Plan2,Larg,Haut,NbReg1,SeuilRecou,&IndiceMax); // bug 17 7
            Sens=CurviDetect2D.MesureRecouv(Plan1,Plan2,Larg,Haut,IndiceReg+1/*NbReg1*/,SeuilRecou,&IndiceMax,TabListeOk);
            //Sens=0;
            if (Sens==12) // on garde la 1 et on vire la 2
            {
               TabListeOk[IndiceReg]=-1; // mark the region as removed
               //MessageBox("region 2 removed");
            }
            if (Sens==21) //recopie région du plan2 dans plan1// if (Surf2>Surf1)
            {

               //TabListeOk[TabIndice[IndiceMax-1]]=-3; 
               TabListeOk[IndiceMax-1]=-3; // mark the region as removed
               for (u=0;u<Taille;u++)
               {
                  if (Plan1[u]==IndiceMax) Plan1[u]=0;
                  if (Plan2[u]!=0 ) Plan1[u]=IndiceMax;
               }
               //MessageBox("region 1 removed");
               //TabIndice[IndiceMax-1]=IndiceReg;
            } // si recouv

            if (Sens==0) 
            {
               NbReg1++;
               for (u=0;u<Taille;u++)
               {
                  //if (Plan2[u]!=0 && Plan1[u]==0) Plan1[u]=NbReg1;//Plan2[u]; bug 17 7
                  if (Plan2[u]!=0 && Plan1[u]==0) Plan1[u]=IndiceReg+1;//NbReg1;//Plan2[u];

               }
            } // si peu de recou

         }// si les listes ne sont pas vide (fin if)
      }//fin for IndiceReg
   } // if more than one region
#endif 



   //for (l=0;l<Haut;l++)
   //{
   //	for (c=0;c<Larg;c++)
   //	{
   //		u=Larg*l+c;
   //		if (Plan1[u]) pImage->ImageRGB[u].r=255-Plan1[u];
   //		if (Plan1[u]) pImage->ImageRGB[u].g=255-Plan1[u]>>8;
   //		if (Plan1[u]) pImage->ImageRGB[u].b=255-Plan1[u]>>16;
   //	}
   //}


   // visualization of final results
   int NbRegFinal=0;
   int LngTableauListesPointsFiltrees2D=TableauListesPointsFiltreesC2D.GetCount();
   for(IndiceReg=0;IndiceReg<LngTableauListesPointsFiltrees2D;IndiceReg++)//pour toutes les régions
   {
      if (TabListeOk[IndiceReg]==1)
      {
         NbRegFinal++;
         BuildImageFromList(Plan1,TableauListesPointsFiltreesC2D[IndiceReg]);
         TireCouleurAlea(&r,&g,&b);
         for (u=0;u<Taille;u++)
         {
            if (Plan1[u]==1) pImage->ImageRGB[u].r=r;
            if (Plan1[u]==1) pImage->ImageRGB[u].g=g;
            if (Plan1[u]==1) pImage->ImageRGB[u].b=b;
         }
      }
   }




   free(TabIndice);
   free(TabListeOk);	


   QueryPerformanceCounter(&TAffichage);
   tps5.QuadPart = TAffichage.QuadPart-T2DProcess.QuadPart;

   AffichageProgression->m_Affichage = (double)tps5.QuadPart/(double)FreqProc.QuadPart;
   AffichageProgression->m_TexteValProgress = "100 %";
   AffichageProgression->m_ControlBarProgressCurvi.SetPos(100);
   AffichageProgression->m_TxtNbReg.Format("Nb Reg=%d",NbRegFinal);
   AffichageProgression->UpdateData(FALSE);
#endif
   //AffichageProgression->Invalidate(TRUE);


   //MessageBox("Export & affichage ok ","Export & affichage ok ",MB_OK);
   //----------------------------------------------------------------------------------------------------
   Invalidate(FALSE);

   SetCursor(LoadCursor(0,IDC_ARROW));

   for(i=0;i<Haut;i++)
   {
      free(APasGaussien.ima[i]);
      free(I_x.ima[i]);
      free(I_y.ima[i]);
   }

   free(I_x.ima);
   free(I_y.ima);
   free(APasGaussien.ima);
   free(TabCoupe);


   free(TabDroiteCoul);
   free(TabGaucheCoul);
   free(ImageFreeman);
   free(ImageGauss);

   free(ImageLabelEdge);

   delete TabCoupeGauche;

   for(IndiceReg=0;IndiceReg<TableauListesPointsFiltreesC2D.GetCount();IndiceReg++)//pour toutes les régions
   {
      TabNumListesConnex[IndiceReg]->RemoveAll();
      delete TabNumListesConnex[IndiceReg];
      TableauListesPointsFiltreesC2D[IndiceReg]->RemoveAll();
      delete TableauListesPointsFiltreesC2D[IndiceReg];
   }
   for(IndiceReg=0;IndiceReg<TableauListesPointsFiltrees.GetCount();IndiceReg++)//pour toutes les régions
   {
      TableauListesPointsFiltrees[IndiceReg]->RemoveAll();
      delete TableauListesPointsFiltrees[IndiceReg];
   }


   free(Plan1);
   free(Plan2);
#endif
}

void CImageWnd::OnTraitementGaussiencouleur()
{
#ifdef MP
   ImagePetrou IR,IB,IG,B;
   rgb RGB;
   //int SeuilPetrou=3;
   int SigmaGaussien=3;
   long l = 0;
   int i,j;

   BYTE   * ImR = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);
   BYTE   * ImV = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);
   BYTE   * ImB = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);
   BYTE   * ImageGauss = (BYTE*) malloc(sizeof(BYTE)*Larg*Haut);

   for (long j= 0;j<Haut;j++) 
   {                              
      for (long i=0;i<Larg;i++) 
      {       
         RGB = ImageRGB[j*Larg+i];
         ((unsigned char   * )ImR)[l] = RGB.r; 
         ((unsigned char   * )ImV)[l] = RGB.g; 
         ((unsigned char   * )ImB)[l] = RGB.b; 
         l++;
      } 
   }	


   IR.dim_x = IG.dim_x=IB.dim_x=Larg;
   IR.dim_y = IG.dim_y =IB.dim_y =Haut;

   B.dim_x = Larg;
   B.dim_y = Haut;


   IR.ima = (Pixelp**) malloc(sizeof(Pixelp*)*IR.dim_y);
   IG.ima = (Pixelp**) malloc(sizeof(Pixelp*)*IG.dim_y);
   IB.ima = (Pixelp**) malloc(sizeof(Pixelp*)*IB.dim_y);

   B.ima = (Pixelp**) malloc(sizeof(Pixelp*)*B.dim_y);

   for(i=0;i<IR.dim_y;i++)
   {
      IR.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*IR.dim_x);
      IG.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*IG.dim_x);
      IB.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*IB.dim_x);
      B.ima[i] = (Pixelp*) malloc(sizeof(Pixelp)*B.dim_x);
   }

   Convolution2DCedricCurvi(ImageGauss,ImR,SigmaGaussien,Larg,Haut);
   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         ImageRGB[i*Larg+j].r = (Pixelp)ImageGauss[i*Larg+j]; // rouge
      }
   }
   Convolution2DCedricCurvi(ImageGauss,ImB,SigmaGaussien,Larg,Haut);
   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         ImageRGB[i*Larg+j].b = (Pixelp)ImageGauss[i*Larg+j]; // bleu
      }
   }
   Convolution2DCedricCurvi(ImageGauss,ImV,SigmaGaussien,Larg,Haut);
   for(i=0;i<Haut;i++)
   {
      for(j=0;j<Larg;j++)
      {
         ImageRGB[i*Larg+j].g = (Pixelp)ImageGauss[i*Larg+j]; // vert
      }
   }



   for(i=0;i<IR.dim_y;i++)
   {
      free(IR.ima[i]);
      free(IG.ima[i]);
      free(IB.ima[i]);
      free(B.ima[i]);

   }

   free(IR.ima);
   free(IG.ima);
   free(IB.ima);
   free(B.ima);


   free(ImageGauss);
   free(ImR);
   free(ImV);
   free(ImB);

   Invalidate(TRUE);
#endif
}

#ifdef MP
int CImageWnd::CompareDouble(const void *a,const void *b)
{
   STabIndice *val1=(STabIndice *)a;
   STabIndice *val2=(STabIndice *)b;
   if (val1->r > val2->r) return 1; // 1
   else if (val2->r > val1->r) return -1; //-1
   else return 0;
}

void CImageWnd::InvertTable(STabIndice *Tabi, STabIndice *Tabf,int Taille)
{
   int i,k;
   for (i=0;i<Taille;i++)
   {
      for (k=0;k<Taille;k++)
      {
         if (Tabi[k].i==i)
         {
            Tabf[i].i=k;
            break;
         }
      }
   }
}

double CImageWnd::Distance(Sfpoint p1,Sfpoint p2)
{
   double d;
   d=(p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y);
   d=sqrt(d);
   return d; 

}
#endif