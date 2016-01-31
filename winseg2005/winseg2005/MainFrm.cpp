// ******************************************************
// Attention : pour les acquistions d'images, il faut 
// enlevées 1  ligne en commentaire 
// ******************************************************             

// include du projet
#include "stdafx.h"
#include "curvlinear.h"
#include "dialogue.h"  
#include "mainfrm.h" 
#include "constdef.h" 
#include "resource.h"
#include "largHaut.h"
#include "lhfloat.h"
#include "timage.h"
#include "Splash.h"
#include "ProgDlg.h"
#include "Steger1D.h"
#include "time.h"
#include "stdlib.h"
#include "math.h"
#include "ima.h"
#include "CurviDetect1D.h"

using namespace utls;
//#include "iostream.h"
// ********************************************                       
// Enlever ou mettre en commentaire pour ne pas faire d'acquisitions
//#define ACQUISITION TRUE

//#ifdef ACQUISITION  
//#include "carte.h"
//#include "tacq.h"
//CTraitAcqWnd * pImageAcqWndContinue;
//#endif

#ifdef _DEBUG
#define new DEBUG_NEW
#endif
#pragma warning( disable : 4996 )

BOOL INSTALL = FALSE;                        
NBCOLOR TypeColor = RIEN;

// Variables d'acquisition  
IMPLEMENT_DYNAMIC(CMainFrame, CMDIFrameWnd)

BEGIN_MESSAGE_MAP(CMainFrame, CMDIFrameWnd)  
   ON_UPDATE_COMMAND_UI(ID_INDICATOR_DATE, OnUpdateDate)
   ON_UPDATE_COMMAND_UI(ID_INDICATOR_TIME, OnUpdateTime)
   //{{AFX_MSG_MAP(CMainFrame)
   ON_COMMAND(ID_CLOSE, OnCloseAll) 
   ON_COMMAND(ID_IMAGE, OnNewImage)
   ON_COMMAND(ID_SEQUENCE,OnOuvrirSequence)
   ON_COMMAND(ID_IMAGEVIDE, OnImagevide)
   ON_COMMAND(ID_SNAP, OnAcquerirImage)
   ON_COMMAND(ID_GRAB, OnAcquerirContinue)
   ON_COMMAND(ID_FILE_MRU_FILE1, OnFileMruFile1)
   ON_COMMAND(ID_FILE_MRU_FILE2, OnFileMruFile2)
   ON_COMMAND(ID_FILE_MRU_FILE3, OnFileMruFile3)
   ON_COMMAND(ID_FILE_MRU_FILE4, OnFileMruFile4)
   ON_UPDATE_COMMAND_UI(ID_FILE_MRU_FILE2, OnUpdateFileMruFile2)
   ON_UPDATE_COMMAND_UI(ID_FILE_MRU_FILE3, OnUpdateFileMruFile3)
   ON_UPDATE_COMMAND_UI(ID_FILE_MRU_FILE4, OnUpdateFileMruFile4)
   ON_WM_CREATE()
   ON_COMMAND(ID_EDIT_TOOLBAR, OnEditToolbar)
   ON_UPDATE_COMMAND_UI(ID_EDIT_TOOLBAR, OnUpdateEditToolbar)
   ON_COMMAND(ID_FILE_EXIT, OnFileExit)
   ON_WM_CLOSE()
   ON_COMMAND(ID_ABOUT, OnAbout)
   ON_UPDATE_COMMAND_UI(ID_GRAB, OnUpdateGrab)
   ON_UPDATE_COMMAND_UI(ID_SNAP, OnUpdateSnap)
   ON_COMMAND(ID_PARZEN3D_CALCULDESESPACESDEPROBA, OnParzen3dCalculdesespacesdeproba)
   ON_COMMAND(ID_PARZEN3D_CRERUNESPACEDEDCISION, OnParzen3dCrerunespacededcision)
   ON_COMMAND(ID_PARZEN3D_OUVRIRUNESPACEDEDCISION, OnParzen3dOuvrirunespacededcision)
   ON_COMMAND(ID_PARZEN3D_SAUVEGARDERUNESPACEDEDCISION, OnParzen3dSauvegarderunespacededcision)
   ON_COMMAND(ID_PARZEN3D_VISUALISERLESPACEDEDCISION, OnParzen3dVisualiserlespacededcision)
   ON_UPDATE_COMMAND_UI(ID_PARZEN3D_SAUVEGARDERUNESPACEDEDCISION, OnUpdateParzen3dSauvegarderunespacededcision)
   ON_UPDATE_COMMAND_UI(ID_PARZEN3D_VISUALISERLESPACEDEDCISION, OnUpdateParzen3dVisualiserlespacededcision)
   ON_COMMAND(ID_PARZEN3D_MESRURERECOUVREMENTSAVECCL0, OnParzen3dMesrurerecouvrementsaveccl0)
   ON_COMMAND(ID_PARZEN3D_MESUREDEVOLUMEDESCLASSES, OnParzen3dMesuredevolumedesclasses)
   ON_UPDATE_COMMAND_UI(ID_PARZEN3D_MESUREDEVOLUMEDESCLASSES, OnUpdateParzen3dMesuredevolumedesclasses)
   ON_COMMAND(ID_TEST_OUVRIR, OnTestOuvrir)
   ON_COMMAND(ID_ZERNIKE_BATCHMOMENTS, OnZernikeBatchmoments)
   ON_COMMAND(ID_ETUDECURVI1D_GENERATIONFICHIERPROFIL_PORTE, OnEtudecurvi1dGenerationfichierprofilPorte)
   ON_COMMAND(ID_ETUDECURVI1D_GENERATIONFICHIERPROFIL_PORTE32, OnEtudecurvi1dGenerationfichierprofilPorte32)
   ON_COMMAND(ID_ETUDECURVI1D_GENERATIONFICHIERPROFIL_PORTE16, OnEtudecurvi1dGenerationfichierprofilPorte16)
   ON_COMMAND(ID_ETUDECURVI1D_GENERATIONFICHIERPROFIL_PORTE12, OnEtudecurvi1dGenerationfichierprofilPorte12)
   ON_COMMAND(ID_ETUDECURVI1D_ETUDEAVECOBJET_ETUDEPORTE32, OnEtudecurvi1dEtudeavecobjetEtudeporte32)
   ON_COMMAND(ID_ETUDECURVI1D_ETUDEAVECOBJET_ETUDEPORTE, OnEtudecurvi1dEtudeavecobjetEtudeporte)
   ON_COMMAND(ID_ETUDECURVI1D_ETUDEAVECOBJET_ETUDEPORTE16, OnEtudecurvi1dEtudeavecobjetEtudeporte16)
   ON_COMMAND(ID_ETUDECURVI1D_ETUDEAVECOBJET_ETUDEPORTE12, OnEtudecurvi1dEtudeavecobjetEtudeporte12)
   ON_COMMAND(ID_ETUDECURVI1D_ETUDEAVECOBJET_ETUDEPORTE12F4, OnEtudecurvi1dEtudeavecobjetEtudeporte12f4)
   ON_COMMAND(ID_ETUDECURVI1D_GENERATIONFICHIERPROFIL_BITEXTURE48, OnEtudecurvi1dGenerationfichierprofilBitexture48)
   ON_COMMAND(ID_ETUDECURVI1D_ETUDEAVECOBJET_ETUDEBITEXTURE48, OnEtudecurvi1dEtudeavecobjetEtudebitexture48)
   ON_COMMAND(ID_ETUDECURVI1D_ETUDEAVECOBJET_ETUDEBITEXTURE48F4, OnEtudecurvi1dEtudeavecobjetEtudebitexture48f4)
   ON_COMMAND(ID_ETUDECURVI1D_GENERATIONFICHIERPROFIL_PORTEBORDNONSIM48, OnEtudecurvi1dGenerationfichierprofilPortebordnonsim48)
   ON_COMMAND(ID_ETUDECURVI1D_ETUDEAVECOBJET_ETUDEPORTENONSIM48F16, OnEtudecurvi1dEtudeavecobjetEtudeportenonsim48f16)
   ON_COMMAND(ID_ETUDECURVI1D_GENERATIONFICHIERPROFIL_MULTISIGNAUX1D, OnEtudecurvi1dGenerationfichierprofilMultisignaux1d)
   ON_COMMAND(ID_ETUDECURVI1D_ETUDEAVECOBJET_ETUDESIGNAUXMULTIPLE48F16, OnEtudecurvi1dEtudeavecobjetEtudesignauxmultiple48f16)
   ON_COMMAND(ID_ETUDECURVI1D_ETUDEAVECOBJET_ETUDESIGNAUXRELS, OnEtudecurvi1dEtudeavecobjetEtudesignauxrels)
   ON_COMMAND(ID_ETUDECURVI1D_ETUDEAVECOBJET_TESTNEW, OnEtudecurvi1dEtudeavecobjetTestnew)
   ON_COMMAND(ID_ETUDECURVI1D_ETUDEAVECOBJET_ETUDEMULTIFICHIER, OnEtudecurvi1dEtudeavecobjetEtudemultifichier)
   //}}AFX_MSG_MAP
   ON_COMMAND(ID_STEGER_CONVOLUTION1D, &CMainFrame::OnStegerConvolution1d)
   ON_COMMAND(ID_CONVERSION_SIGNAL, &CMainFrame::OnConversionSignal)
   ON_COMMAND(ID_OUTILS_AJOUTBRUIT, &CMainFrame::OnOutilsAjoutbruit)
END_MESSAGE_MAP()



CMainFrame::CMainFrame()
{
   premier=TRUE;
   EspValid=FALSE;
}                                                   


CMainFrame::~CMainFrame()
{
#ifdef ACQUISITION  
   if (INSTALL) FermeCarte();
#endif
}

void CMainFrame::OnFileExit() 
{
   OnClose();
}

void CMainFrame::OnClose() 
{
   int ferm=IDYES;
   while (MDIGetActive()!=NULL && ferm!=IDCANCEL)
   {
      if(MDIGetActive()->IsKindOf( RUNTIME_CLASS( CTraitImageWnd)))
      {
         CTraitImageWnd *ptr = (CTraitImageWnd *)MDIGetActive();
         ferm=ptr->CloseFen();
         //		if (IsWindow(ptr->GetSafeHwnd()))
         //			return;
      }
      else MDIGetActive()->MDIDestroy();
   }
   //DestroyWindow();
   if(ferm!=IDCANCEL) CMDIFrameWnd::OnClose();
}

afx_msg void  CMainFrame::OnCloseAll()
{   
   //	while (MDIGetActive()!=NULL) MDIGetActive()-> MDIDestroy();           
   int ferm=IDYES;
   while (MDIGetActive()!=NULL && ferm!=IDCANCEL)
   {
      if(MDIGetActive()->IsKindOf( RUNTIME_CLASS( CTraitImageWnd)))
      {
         CTraitImageWnd *ptr = (CTraitImageWnd *)MDIGetActive();
         ferm=ptr->CloseFen();
         //		if (IsWindow(ptr->GetSafeHwnd()))
         //			return;
      }
      else MDIGetActive()->MDIDestroy();
   }
}

afx_msg void CMainFrame::OnAcquerirImage()
{   
   // **********************************************
   // A ENLEVER POUR LES ACQUISITIONS ET INVERSEMENT...
#ifdef ACQUISITION
   if (!INSTALL) 
   {   
      TypeColor = InitialisationCarte();   
      INSTALL = TRUE;
   }
   if (TypeColor!=RIEN) 
   {
      SetCursor(LoadCursor(0,IDC_WAIT));
      CTraitAcqWnd * pImageAcqWnd =  (CTraitAcqWnd * ) new CTraitAcqWnd(IMAGE_FIXE,TypeColor);  
      pImageAcqWnd->Create("Acquisition de l'image",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,this);
      pImageAcqWnd->OnSnap();
      SetCursor(LoadCursor(0,IDC_ARROW));
   }
   else MessageBox("erreur d'initialisation de la carte","Message");
#endif    
} 

afx_msg void CMainFrame::OnAcquerirContinue()
{  
   // **********************************************
   // A ENLEVER POUR LES ACQUISITIONS ET INVERSEMENT ....
#ifdef ACQUISITION
   if (!INSTALL) 
   {
      TypeColor = InitialisationCarte();   
      INSTALL = TRUE;
   }	
   if (TypeColor!=RIEN) 
   { 
      SetCursor(LoadCursor(0,IDC_WAIT));

      pImageAcqWndContinue = (CTraitAcqWnd * ) new CTraitAcqWnd(IMAGE_CONTINUE,TypeColor);   
      pImageAcqWndContinue->Create("Acquisition de l'image", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,this);
      pImageAcqWndContinue->OnGrab();
      //	bar.OnUpdateCmdUI( (CFrameWnd*) this, FALSE);
      SetCursor(LoadCursor(0,IDC_ARROW));
   }   
   else MessageBox("erreur d'initialisation de la carte","Message");
#endif    
} 

afx_msg void CMainFrame::OnNewImage()           
{
   CFile Fichier; // variable fichier BITMAP	  
   long Larg, Larg2, Haut;
   DWORD DWComp; // test la compression du fichier bitmap
   WORD BitCount;
   CTraitImageWnd	* pImageWnd;
   CFileDialog * FileData ;
   // ouverture de la boite de dialogue d'ouverture d'un fichier
   if (premier) {FileData= new CFileDialog (TRUE,"*.bmp;*.mem;*.f32","c:\\bmp\\*.bmp",OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"Common files (*.bmp;*.mem;*.f32)|*.bmp;*.mem;*.f32|BITMAP (*.bmp)|*.bmp|BRUTE (*.mem)|*.mem|Image Float (*.f32)|*.f32||",this);premier=FALSE;}
   else FileData= new CFileDialog (TRUE,"*.bmp;*.mem;*.f32",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"Common files (*.bmp;*.mem;*.f32)|*.bmp;*.mem;*.f32|BITMAP (*.bmp)|*.bmp|BRUTE (*.mem)|*.mem|Image Float (*.f32)|*.f32||",this);
   if (FileData->DoModal()==IDOK)
   {
      if (FileData->GetFileExt() == "bmp" || FileData->GetFileExt() == "f32" || FileData->GetFileExt() == "mem")
      {
         if (Fichier.Open(FileData->GetPathName(), CFile::modeRead))
         { 
            CString NomTemp;
            if (FileData->GetFileExt() == "bmp")
            {
               Fichier.Seek(sizeof(BITMAPFILEHEADER)+sizeof(DWORD),CFile::begin);
               Fichier.Read(&Larg,sizeof(LONG)); 
               if (Larg%4!=0)
               {
                  MessageBox("la largeur de l'image doit etre multiple de 4","Message",MB_OK|MB_ICONEXCLAMATION);
                  Fichier.Close();
                  return;
               }

               Fichier.Read(&Haut,sizeof(LONG));
               Fichier.Seek(sizeof(WORD),CFile::current);
               Fichier.Read(&BitCount,sizeof(WORD));
               Fichier.Read(&DWComp,sizeof(DWORD));
               if ((DWComp != BI_RGB)||(BitCount!=8 && BitCount!=24))
               {
                  MessageBox("Impossible de lire les fichiers BMP compressés","Message",MB_OK|MB_ICONEXCLAMATION);
                  Fichier.Close();
                  return;
               }
               Fichier.Close();
               pImageWnd = new CTraitImageWnd(FileData->GetPathName(), Larg, Haut, BMP);
               if ((FileData->GetPathName()!=AfxGetApp()->GetProfileString("Recent File List","File1"))
                  &&(FileData->GetPathName()!=AfxGetApp()->GetProfileString("Recent File List","File2"))
                  &&(FileData->GetPathName()!=AfxGetApp()->GetProfileString("Recent File List","File3"))
                  &&(FileData->GetPathName()!=AfxGetApp()->GetProfileString("Recent File List","File4")))
               {
                  NomTemp=AfxGetApp()->GetProfileString("Recent File List",  "File3");
                  AfxGetApp()->WriteProfileString("Recent File List",  "File4",NomTemp);
                  NomTemp=AfxGetApp()->GetProfileString("Recent File List",  "File2");
                  AfxGetApp()->WriteProfileString("Recent File List",  "File3",NomTemp);
                  NomTemp=AfxGetApp()->GetProfileString("Recent File List",  "File1");
                  AfxGetApp()->WriteProfileString("Recent File List",  "File2",NomTemp);
                  AfxGetApp()->WriteProfileString("Recent File List",  "File1",FileData->GetPathName());

                  AfxGetApp()->AddToRecentFileList( FileData->GetPathName());
               }
               //AfxGetApp()->AddToRecentFileList( FileData->GetPathName());
            }
            else if (FileData->GetFileExt() == "mem")
            {
               LHByte dial(Fichier.GetLength());
               if (dial.DoModal()==IDOK)
               {
                  Larg = dial.m_LARG; 
                  Haut = dial.m_HAUT;
                  Larg2 = (long)(div(Larg,(long)4).quot)*4;
                  Fichier.Close();
                  pImageWnd = new CTraitImageWnd(FileData->GetPathName(), Larg2, Haut, MEM);
               }
               else return;
            }
            else if (FileData->GetFileExt() == "f32")
            {
               LHfloat dial(Fichier.GetLength());	
               if (dial.DoModal()==IDOK)
               {
                  Fichier.Close();
                  pImageWnd = new CTraitImageWnd(FileData->GetPathName(),dial.m_LARG,dial.m_HAUT,FLOAT32);
               }
               else return;
            }
            pImageWnd->Create(FileData->GetFileTitle(), WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,this);
         }
         else MessageBox("Ce fichier n'existe pas","Message",MB_OK|MB_ICONEXCLAMATION);
      } 
      else MessageBox("Impossible d'ouvrir \n ce type de fichier image","Fichier autre que bmp, mem ou f32",MB_OK|MB_ICONEXCLAMATION);
   }
   delete FileData;
}

afx_msg void CMainFrame::OnOuvrirSequence()           
{     
   BOOL first = TRUE;
   int  num1=1,num2=35;  
   CFile Fichier; // variable fichier BITMAP	 
   int l; 
   CString NomFichier,s;
   char NumExt[20];  
   CFileDialog  *FileData;
   if (premier)FileData=new CFileDialog(TRUE,"*.bmp","c:\\bmp\\*.bmp",OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"BITMAP (*.bmp)|*.bmp||",this);
   else FileData=new CFileDialog(TRUE,"*.bmp",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"BITMAP (*.bmp)|*.bmp||",this);
   if (FileData->DoModal()==IDOK) 
   {
      CDialSeqlimite dial(&num1,&num2); 
      if (dial.DoModal()==IDOK)
      {
         SetCursor(LoadCursor(0,IDC_WAIT));  
         BITMAPFILEHEADER * TheTete = (BITMAPFILEHEADER*) new char[sizeof(BITMAPFILEHEADER)];
         BITMAPINFO * 	TheInfo = (BITMAPINFO *) new unsigned char[sizeof(BITMAPINFOHEADER)+256*sizeof(RGBQUAD)];     
         CTraitImageWnd * pImageWnd;
         l = FileData->GetPathName().GetLength();
         s = FileData->GetPathName().Left(l-7);  
         int NbMax = num2-num1+1;                                   
         for (int i=0;i<NbMax;i++)
         {                             
            sprintf(NumExt,"%.3d.bmp",i+num1); 
            NomFichier = s + NumExt;  
            if (Fichier.Open(NomFichier,CFile::modeRead)!=0)
            {
               Fichier.Read(TheTete,sizeof(BITMAPFILEHEADER)); 
               Fichier.Read(TheInfo,sizeof(BITMAPINFOHEADER)+256*sizeof(RGBQUAD));      
               if (first)
               {	
                  pImageWnd = new CTraitImageWnd(TheInfo->bmiHeader.biWidth,TheInfo->bmiHeader.biHeight,NB,TRUE); 
                  first = FALSE;
                  pImageWnd->Create(FileData->GetFileTitle().Left(FileData->GetFileTitle().GetLength()-3), WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,this);  
                  pImageWnd->AlloueSequence(&NbMax);   
                  Fichier.Read(pImageWnd->Bpile[0],TheInfo->bmiHeader.biWidth*TheInfo->bmiHeader.biHeight); 
                  for(long j=0;j<pImageWnd->Larg*pImageWnd->Haut;j++) pImageWnd->Image[j]=pImageWnd->Bpile[0][j];
                  CString Titreseq;
                  ((CWnd*)pImageWnd)->GetWindowText(Titreseq);
                  Titreseq+=": 1";
                  ((CWnd*)pImageWnd)->SetWindowText(Titreseq);
               }
               else Fichier.Read(pImageWnd->Bpile[i],TheInfo->bmiHeader.biWidth*TheInfo->bmiHeader.biHeight);  
               Fichier.Close(); 

            }	
         } 
         delete TheTete; 
         delete TheInfo;    

         SetCursor(LoadCursor(0,IDC_ARROW));
      }
   }  
   delete FileData;
}
void CMainFrame::OnImagevide()
{
   // TODO: Add your command handler code here
   CLargHaut dial;
   if (dial.DoModal()==IDOK) 
   {
      CTraitImageWnd	* pImageWnd; 
      if (dial.COLOR) pImageWnd = new CTraitImageWnd(dial.m_LARG,dial.m_HAUT,COLOR);
      else pImageWnd = new CTraitImageWnd(dial.m_LARG,dial.m_HAUT,BMP);
      pImageWnd->Create(dial.m_TITRE, WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,this);

   }			
}

void CMainFrame::OnFileMruFile1() 
{
   // TODO: Add your command handler code here
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
      pImageWnd->Create(NomFich.Left(NomFich.GetLength()-4), WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,this);

   }
   else MessageBox("Ce fichier n'existe pas","Message",MB_OK|MB_ICONEXCLAMATION);
}

void CMainFrame::OnFileMruFile2() 
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
      pImageWnd->Create(NomFich.Left(NomFich.GetLength()-4), WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,this);

   }
   else MessageBox("Ce fichier n'existe pas","Message",MB_OK|MB_ICONEXCLAMATION);
}

void CMainFrame::OnFileMruFile3() 
{
   // TODO: Add your command handler code here
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
      pImageWnd->Create(NomFich.Left(NomFich.GetLength()-4), WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,this);

   }
   else MessageBox("Ce fichier n'existe pas","Message",MB_OK|MB_ICONEXCLAMATION);

}

void CMainFrame::OnFileMruFile4() 
{
   // TODO: Add your command handler code here
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
      pImageWnd->Create(NomFich.Left(NomFich.GetLength()-4), WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,this);

   }
   else MessageBox("Ce fichier n'existe pas","Message",MB_OK|MB_ICONEXCLAMATION);

}

void CMainFrame::OnUpdateFileMruFile2(CCmdUI* pCmdUI) 
{
   // TODO: Add your command update UI handler code here
   if (AfxGetApp()->GetProfileString( "Recent File List",  "File1")=="")
      pCmdUI->m_pMenu->DeleteMenu(pCmdUI->m_nID,MF_BYCOMMAND);
}

void CMainFrame::OnUpdateFileMruFile3(CCmdUI* pCmdUI) 
{
   // TODO: Add your command update UI handler code here
   if (AfxGetApp()->GetProfileString( "Recent File List",  "File2")=="")
      pCmdUI->m_pMenu->DeleteMenu(pCmdUI->m_nID,MF_BYCOMMAND);
}

void CMainFrame::OnUpdateFileMruFile4(CCmdUI* pCmdUI) 
{
   // TODO: Add your command update UI handler code here
   if (AfxGetApp()->GetProfileString( "Recent File List",  "File3")=="")
      pCmdUI->m_pMenu->DeleteMenu(pCmdUI->m_nID,MF_BYCOMMAND);
}



int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{

   if (CMDIFrameWnd::OnCreate(lpCreateStruct) == -1)
      return -1;
   // TODO: Add your specialized creation code here
   CSplashWnd::ShowSplashScreen(this);
   //////////////////////

   // TODO: Add your specialized creation code here
   //	CSplashWnd::ShowSplashScreen(this);
   // CG: The following block was inserted by 'Status Bar' component.
   {	
      // Create an array for status bar indicators
      UINT pIndicators[6] = {ID_SEPARATOR, ID_SEPARATOR,ID_INDICATOR_CAPS, ID_INDICATOR_NUM};
      if (!m_wndStatusBar.Create(this, 
         WS_CHILD | WS_VISIBLE | CBRS_BOTTOM, AFX_IDW_STATUS_BAR/*ID_MY_STATUS_BAR*/) ||
         !InitStatusBar(pIndicators, 4, 60))
      {
         TRACE0("Failed to create Status Bar\n");
         return -1;
      }

      UINT nID, nStyle;
      int nWidth;
      m_wndStatusBar.GetPaneInfo(1, nID, nStyle, nWidth);
      m_wndStatusBar.SetPaneInfo(1, nID, nStyle, 150);
      // TODO: If you have a call to EnableDocking() and DockControlBar()
      // before this line, move them to right after this comment.  This
      // insures that the docking area remains above the status bar.
   }
   bar.Create(this);
   bar.m_bAutoDelete=0;
   bar.LoadToolBar(IDR_MAINFRAME);
   // TODO: Remove this if you don't want tool tips or a resizeable toolbar
   bar.SetBarStyle(bar.GetBarStyle() |
      CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC);

   // TODO: Delete these three lines if you don't want the toolbar to
   //  be dockable
   /*	bar.EnableDocking(CBRS_ALIGN_ANY);
   EnableDocking(CBRS_ALIGN_ANY);
   DockControlBar(&bar);
   */	
   //////////////////////
   CToolBarCtrl& pToolBarCtrl=bar.GetToolBarCtrl();
   pToolBarCtrl.HideButton(ID_VISUALISATION_CHOIXDESPARAMETRES,TRUE);
   pToolBarCtrl.HideButton(ID_VISUALISATION_VISUALISERLESPOLYTOPES,TRUE);
   pToolBarCtrl.HideButton(ID_PARZEN,TRUE);

   return 0;
}

void CMainFrame::OnEditToolbar() 
{
   // TODO: Add your command handler code here
   ShowControlBar(&bar,(bar.GetStyle() & WS_VISIBLE)==0, FALSE);
}

void CMainFrame::OnUpdateEditToolbar(CCmdUI* pCmdUI) 
{
   // TODO: Add your command update UI handler code here
   pCmdUI->SetRadio((bar.GetStyle() & WS_VISIBLE)!=0);

}

void CMainFrame::OnUpdateDate(CCmdUI* pCmdUI)
{
   // CG: This function was inserted by 'Status Bar' component.

   // Get current date and format it
   CTime time = CTime::GetCurrentTime();
   CString strDate = time.Format(_T("%A, %B %d, %y "));

   // BLOCK: compute the width of the date string
   CSize size;
   {
      HGDIOBJ hOldFont = NULL;
      HFONT hFont = (HFONT)m_wndStatusBar.SendMessage(WM_GETFONT);
      CClientDC dc(NULL);
      if (hFont != NULL) 
         hOldFont = dc.SelectObject(hFont);
      size = dc.GetTextExtent(strDate);
      if (hOldFont != NULL) 
         dc.SelectObject(hOldFont);
   }

   // Update the pane to reflect the current date
   UINT nID, nStyle;
   int nWidth;
   m_wndStatusBar.GetPaneInfo(m_nDatePaneNo, nID, nStyle, nWidth);
   m_wndStatusBar.SetPaneInfo(m_nDatePaneNo, nID, nStyle, size.cx);
   pCmdUI->SetText(strDate);
   pCmdUI->Enable(TRUE);

}

void CMainFrame::OnUpdateTime(CCmdUI* pCmdUI)
{
   // CG: This function was inserted by 'Status Bar' component.

   // Get current date and format it
   CTime time = CTime::GetCurrentTime();
   CString strTime = time.Format(_T("%X"));

   // BLOCK: compute the width of the date string
   CSize size;
   {
      HGDIOBJ hOldFont = NULL;
      HFONT hFont = (HFONT)m_wndStatusBar.SendMessage(WM_GETFONT);
      CClientDC dc(NULL);
      if (hFont != NULL) 
         hOldFont = dc.SelectObject(hFont);
      size = dc.GetTextExtent(strTime);
      if (hOldFont != NULL) 
         dc.SelectObject(hOldFont);
   }

   // Update the pane to reflect the current time
   UINT nID, nStyle;
   int nWidth;
   m_wndStatusBar.GetPaneInfo(m_nTimePaneNo, nID, nStyle, nWidth);
   m_wndStatusBar.SetPaneInfo(m_nTimePaneNo, nID, nStyle, size.cx);
   pCmdUI->SetText(strTime);
   pCmdUI->Enable(TRUE);

}

BOOL CMainFrame::InitStatusBar(UINT *pIndicators, int nSize, int nSeconds)
{
   // CG: This function was inserted by 'Status Bar' component.

   // Create an index for the DATE pane
   m_nDatePaneNo = nSize++;
   pIndicators[m_nDatePaneNo] = ID_INDICATOR_DATE;
   // Create an index for the TIME pane
   m_nTimePaneNo = nSize++;
   nSeconds = 1;
   pIndicators[m_nTimePaneNo] = ID_INDICATOR_TIME;

   // TODO: Select an appropriate time interval for updating
   // the status bar when idling.
   m_wndStatusBar.SetTimer(0x1000, nSeconds * 1000, NULL);

   return m_wndStatusBar.SetIndicators(pIndicators, nSize);

}

void CMainFrame::OnUpdateGrab(CCmdUI* pCmdUI) 
{
   // TODO: Add your command update UI handler code here
#ifndef ACQUISITION
   pCmdUI->Enable(0);
#endif
   pCmdUI->SetCheck(0);
}

void CMainFrame::OnUpdateSnap(CCmdUI* pCmdUI) 
{
   // TODO: Add your command update UI handler code here
#ifndef ACQUISITION
   pCmdUI->Enable(0);
#endif    
}

void CMainFrame::OnAbout() 
{
   MessageBox("Laboratoire Le2I\nUniversité de Bourgogne\nIUT Le Creusot\n12, Rue de la fonderie, 71200 Le Creusot\ntel. +33 (0)3.85.73.10.00\nfax. +33 (0)3.85.73.10.99","Information",MB_ICONASTERISK & MB_ICONINFORMATION);
}
///////////////////////////////////////////////////////////////
///////// PARZEN 3D ///////////////////////////////////////////

int PGCD(int a, int b)
{
   if(a<b) return PGCD(b,a);
   else if(b==0) return a;
   else return PGCD(b,a%b);
}

void CMainFrame::OnParzen3dCalculdesespacesdeproba() 
{
   MessageBox("Noyau Gaussien");
   CFile fi;
   /*	unsigned __int16  attention*/ long NbMes[256];
   long kmax[256],c;
   BYTE **Va;
   BYTE NbClass,NPAR;
   CString ch;

   unsigned long* espace=new unsigned long[256*256*256];
   //////////////////////////////////////////////////////
   ////////////Lecture des mesures///////////////////////

   CFileDialog * FileData = new CFileDialog(TRUE,"*.mmc",NULL,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT,
      "Mesures mulyi-classes(*.mmc)|*.mmc||",this);

   if (FileData->DoModal()!=IDOK) {delete espace;return;}

   SetCursor(LoadCursor(NULL,IDC_WAIT));

   if(FileData->GetFileExt()!=CString("mmc")) {MessageBox("Impossible d'ouvrir les fichiers autres que mmc ","Warning",MB_ICONEXCLAMATION);delete espace;delete FileData;return;}
   if(fi.Open(FileData->GetPathName(),CFile::modeRead))
   {
      fi.Read(&NbClass,1);
      fi.Read(&NPAR,1);
      if(NPAR!=3/*||NbClass>8*/) {MessageBox("impossible pour des dimensions autres que 3 ou plus de 8 classes");delete espace;delete FileData;return;}

      Va=new BYTE *[NbClass];

      for (c=0;c<NbClass;c++)
      {
         /*attention fi.Read(&(NbMes[c]),2);*/
         fi.Read(&(NbMes[c]),4);

         kmax[c]=(long)NbMes[c]*(long)NPAR;

         Va[c]=new BYTE[kmax[c]];

         fi.Read(Va[c],kmax[c]);


      }//i
      fi.Close();
   }
   else
   {
      MessageBox("Impossible d'ouvrir le fichier "+FileData->GetPathName(),"Warning",MB_ICONEXCLAMATION);
      delete espace;
      delete FileData;
      return;
   }
   delete FileData;
   SetCursor(LoadCursor(NULL,IDC_WAIT));
   ///////////////////////////////////////////////////////
   /////////////remplissage de l'espace//////
   int h;
   ch="Coefficient h";
   CDialVal Dlg(&ch);
   if(Dlg.DoModal()!=IDOK) 
   {
      for(int i=0;i<NbClass;i++) delete Va[i];
      delete Va;
      delete espace;
      return;
   }
   else
   {
      h=atoi(Dlg.m_val);
      if(h<=0||h>30) {MessageBox("h doit être compris entre 1 et 30 !","Warning",MB_OK|MB_ICONEXCLAMATION);for(int i=0;i<NbClass;i++) delete Va[i];delete Va;delete espace;delete FileData;return;}
   }
   ////cacul des poids de Parzen
   int dim, dim2;
   long p,r,t;
   dim=2*h+1;
   dim2=(2*h+1)*(2*h+1);
   unsigned long som=0;
   unsigned long min=0xffffffff,min2=0xffffffff;
   unsigned long * phi=new unsigned long[dim2*dim];
   double h3pi=h*h*h;//acos(-1)*h*h*h;//
   long i;
   for(i=0;i<dim;i++)
      for(int j=0;j<dim;j++)
         for(int k=0;k<dim;k++)
         {
            p=i+j*dim+k*dim2;

            r=((i-h)*(i-h)+(j-h)*(j-h)+(k-h)*(k-h));

            if(r<h*h) 
            {
               //	phi[p]=unsigned long(300000000*double(-(double)(sqrt(r)/h)+1)/h3pi);
               phi[p]=unsigned long(3000000*double(exp(-(4.5*(double)r/(h*h)))/h3pi));
               if(phi[p]==0) {MessageBox("mauvaise dynamique");return;}
               if(phi[p]<=min) min=phi[p];
               else if(phi[p]<min2) min2=phi[p];
            }
            else phi[p]=0;
            som+=phi[p];
         }

         long max=NbMes[0],place;
         for(i=1;i<NbClass;i++) 	if(max<NbMes[i]){max=NbMes[i];place=i;}
         ch.Format("volume=%d\nmin=%d\nmin2=%d\nmax mes=%d",som,min,min2,max);
         MessageBox(ch);

         ///calcul de l'espace prob
         ////calcul coefficient
         long cl,mes=1,val;

         for(cl=0;cl<NbClass;cl++)
         {
            if(cl!=place)
            {
               val=NbMes[cl]/PGCD(NbMes[cl],max);
               if((mes%val)) mes*=val;
            }

         }

         int k,l,m;
         unsigned long pos;
         CFile fich;
         unsigned long compt;
         unsigned long constante,temp;

         for( cl=0;cl<NbClass;cl++)
         {
            for(i=0;i<16777216;i++) espace[i]=0;
            constante=(mes*max)/NbMes[cl];
            compt=0;
            for(i=0;i<NbMes[cl];i++)
            {
               p=Va[cl][i*NPAR]+256*Va[cl][1+i*NPAR]+65536*Va[cl][2+i*NPAR];

               for(k=-h;k<=h;k++)
                  for(l=-h;l<=h;l++)
                     for(m=-h;m<=h;m++)
                     {
                        t=k+l*256+m*65536;
                        temp=constante*phi[k+l*dim+m*dim2+h*(1+dim+dim2)];
                        compt+=temp;
                        if((Va[cl][i*NPAR]+k>=0)&&(Va[cl][i*NPAR]+k<256)&&
                           (Va[cl][1+i*NPAR]+l>=0)&&(Va[cl][1+i*NPAR]+l<256)&&
                           (Va[cl][2+i*NPAR]+m>=0)&&(Va[cl][2+i*NPAR]+m<256)) 
                        {
                           espace[p+t]+=temp;
                        }
                     }
            }	
            ch.Format("tot=%u",compt);
            MessageBox(ch);

            ch.Format("%d-h%d.esp",cl,h);
            if(!fich.Open("d:\\temp\\probclas-"+ch,CFile::modeCreate | CFile::modeWrite)) return;

            for(i=0;i<16777216;i++) 
            {
               pos=i;
               if(espace[pos]!=0) { fich.Write(&pos,4);fich.Write(&espace[pos],4);}
            }
            fich.Close();
         }


         for(i=0;i<NbClass;i++) delete Va[i];
         delete Va;

         delete espace;
         delete phi;	
}

void CMainFrame::OnParzen3dCrerunespacededcision() 
{
   SetCursor(LoadCursor(NULL,IDC_WAIT));
   CFile fich;
   unsigned long tampon;
   unsigned long* espace=new unsigned long[16777216];
   unsigned long pos;
   double compt;
   CProgressDlg Dlg;
   CString ch;
   long i;
   for(i=0;i<16777216;i++) {espace[i]=0;EspClas[i]=0;}

   CFileDialog * FileData = new CFileDialog(TRUE,"*.esp",NULL,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT,
      "Mesures proba(*.esp)|*.esp||",this);
   for(byte k=0;k<256;k++)
   {
      ch.Format("Choix du fichier de proba.:classe %d",k);
      FileData->m_ofn.lpstrTitle=(LPCSTR)ch;

      if (FileData->DoModal()!=IDOK){break;}

      Dlg.Create();
      ch.Format("loading classe %d...",k);
      Dlg.SetWindowText(ch);
      if(!fich.Open(FileData->GetPathName(),CFile::modeRead)) {Dlg.DestroyWindow();delete espace;return;}

      Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
      Dlg.SetStep(1);

      compt=0;
      for(i=0;i<(long)fich.GetLength();i+=8)
      {
         fich.Read(&pos,4);
         fich.Read(&tampon,4);
         if(tampon>espace[pos]) {espace[pos]=tampon;	EspClas[pos]=k+1;}
         compt+=tampon;
         if(!(i%(100000)))
         {
            Dlg.StepIt();
            if(Dlg.CheckCancelButton()) {delete espace;return;}
         }
      }
      fich.Close();
      ch.Format("tot=%1.16e",compt);
      MessageBox(ch);
      Dlg.DestroyWindow();
   }

   /*	if (FileData->DoModal()!=IDOK){delete espace; return;}

   Dlg.Create();
   ch.Format("loading classe 1...");
   Dlg.SetWindowText(ch);

   if(!fich.Open(FileData->GetPathName(),CFile::modeRead)) {Dlg.DestroyWindow();delete espace;return;}
   compt=0;

   Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
   Dlg.SetStep(1);
   for(i=0;i<(long)fich.GetLength();i+=8)
   {
   fich.Read(&pos,4);
   fich.Read(&tampon,4);
   if(tampon>espace[pos]) {espace[pos]=tampon;	EspClas[pos]=2;}
   compt+=tampon;
   if(!(i%(100000)))
   {
   Dlg.StepIt();
   if(Dlg.CheckCancelButton()) {delete espace;return;}
   }
   }
   fich.Close();
   ch.Format("tot=%1.16e",compt);
   MessageBox(ch);
   compt=0;
   Dlg.DestroyWindow();
   Dlg.Create();
   ch.Format("loading classe 2...");
   Dlg.SetWindowText(ch);
   POSITION FilePos;
   if(FileData->GetNextPathName(FilePos)==NULL) {Dlg.DestroyWindow();delete espace;return;}
   if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();delete espace;return;}

   Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
   Dlg.SetStep(1);
   for(i=0;i<(long)fich.GetLength();i+=8)
   {
   fich.Read(&pos,4);
   fich.Read(&tampon,4);
   if(tampon>espace[pos]) {espace[pos]=tampon;	EspClas[pos]=3;}
   compt+=tampon;
   if(!(i%(100000)))
   {
   Dlg.StepIt();
   if(Dlg.CheckCancelButton()) {delete espace;return;}
   }
   }
   fich.Close();
   ch.Format("tot=%1.16e",compt);
   MessageBox(ch);
   compt=0;
   Dlg.DestroyWindow();
   Dlg.Create();
   ch.Format("loading classe 3...");
   Dlg.SetWindowText(ch);

   if(FileData->GetNextPathName(FilePos)==NULL) {Dlg.DestroyWindow();delete espace;return;}
   if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();delete espace;return;}

   Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
   Dlg.SetStep(1);
   for(i=0;i<(long)fich.GetLength();i+=8)
   {
   fich.Read(&pos,4);
   fich.Read(&tampon,4);
   if(tampon>espace[pos]) {espace[pos]=tampon;	EspClas[pos]=4;}
   else if(EspClas[pos]==1) EspClas[pos]=6;
   compt+=tampon;
   if(!(i%(100000)))
   {
   Dlg.StepIt();
   if(Dlg.CheckCancelButton()) {delete espace;return;}
   }
   }
   fich.Close();

   ch.Format("tot=%1.16e",compt);
   MessageBox(ch);
   compt=0;
   Dlg.DestroyWindow();
   Dlg.Create();
   ch.Format("loading classe 4...");
   Dlg.SetWindowText(ch);

   if(FileData->GetNextPathName(FilePos)==NULL) {Dlg.DestroyWindow();delete espace;return;}
   if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();delete espace;return;}

   Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
   Dlg.SetStep(1);

   for(i=0;i<(long)fich.GetLength();i+=8)
   {
   fich.Read(&pos,4);
   fich.Read(&tampon,4);
   if(tampon>espace[pos]) {espace[pos]=tampon;	EspClas[pos]=5;}
   compt+=tampon;
   if(!(i%(100000)))
   {
   Dlg.StepIt();
   if(Dlg.CheckCancelButton()) {delete espace;return;}
   }
   }
   fich.Close();
   ch.Format("tot=%1.16e",compt);
   MessageBox(ch);
   Dlg.DestroyWindow();

   /*	int max=256;
   CTraitImageWnd *toto=new CTraitImageWnd(256,256,FLOAT32);
   toto->Create("projection", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,NULL);
   //	BOOL compt;
   CTraitImageWnd *totoseq=new CTraitImageWnd(256,256,NB,TRUE);
   totoseq->Create("projection en couches", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,NULL);
   totoseq->AlloueSequence(&max);

   for(long j=0;j<max;j++)
   {
   for(i=0;i<65536;i++) 
   {
   //	compt=FALSE;
   pos=i+65536*j;
   if(espace[pos]!=0) toto->ImageFloat[i]=(float)espace[i+65536*j];
   else toto->ImageFloat[i]=0;
   }
   toto->FloatToImage();
   for(i=0;i<65536;i++) totoseq->Bpile[j][i]=toto->Image[i];
   }
   toto->DestroyWindow();	
   */	delete espace;
   delete FileData;
   EspValid=TRUE;
}

void CMainFrame::OnParzen3dOuvrirunespacededcision() 
{
   unsigned long pos;
   CFileDialog * FileData = new CFileDialog(TRUE,"*.dec","",OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT,
      "Espace de décision (*.dec)|*.dec||",this);
   if (FileData->DoModal()!=IDOK) return;

   for(unsigned long i=0;i<16777216;i++)
   {
      EspClas[i]=0;
   }

   CFile fich;
   if(!fich.Open(FileData->GetPathName(),CFile::modeRead)) return;

   CProgressDlg Dlg;
   Dlg.Create();
   CString ch;
   ch.Format("loading "+FileData->GetFileName()+"...");
   Dlg.SetWindowText(ch);
   Dlg.SetRange(0,(long)fich.GetLength()/100000);
   Dlg.SetStep(1);

   for(long j=0;j<(long)fich.GetLength();j+=5)
   {
      fich.Read(&pos,4);
      fich.Read(&EspClas[pos],1);
      if(!(j%(100000)))
      {
         Dlg.StepIt();
         if(Dlg.CheckCancelButton()) {fich.Close();return;}
      }

   }
   fich.Close();
   Dlg.DestroyWindow();
   delete FileData;
   EspValid=TRUE;
}

void CMainFrame::OnParzen3dSauvegarderunespacededcision() 
{
   CFileDialog * FileData = new CFileDialog(FALSE,"*.dec","",OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT,
      "Espace de décision (*.dec)|*.dec||",this);
   if (FileData->DoModal()!=IDOK) return;

   CFile fich;
   if(!fich.Open(FileData->GetPathName(),CFile::modeWrite|CFile::modeCreate)) return;

   CProgressDlg Dlg;
   Dlg.Create();
   CString ch;
   ch.Format("sauvegarde fichier "+FileData->GetFileName()+"...");
   Dlg.SetWindowText(ch);
   Dlg.SetRange(0,167);
   Dlg.SetStep(1);
   for(unsigned long i=0;i<16777216;i++)
   {
      if(EspClas[i]!=0) {fich.Write(&i,4); fich.Write(&EspClas[i],1);}
      if(!(i%(100000)))
      {
         Dlg.StepIt();
         if(Dlg.CheckCancelButton()) {fich.Close();CFile::Remove(FileData->GetPathName());return;}
      }
   }
   fich.Close();
   Dlg.DestroyWindow();
   delete FileData;
}

void CMainFrame::OnUpdateParzen3dSauvegarderunespacededcision(CCmdUI* pCmdUI) 
{
   pCmdUI->Enable(EspValid);
}

void CMainFrame::OnParzen3dVisualiserlespacededcision() 
{
   int max=256;
   CTraitImageWnd *coupedec=new CTraitImageWnd(256,256,NB,TRUE);
   coupedec->Create("projection en couches", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,NULL);
   coupedec->AlloueSequence(&max);

   for(long j=0;j<max;j++)
   {
      for(long i=0;i<65536;i++) if((EspClas[i+65536*j]!=0)&&(EspClas[i+65536*j]!=1)) coupedec->Bpile[j][i]=42*(EspClas[i+65536*j]-1);else coupedec->Bpile[j][i]=20*EspClas[i+65536*j];//toto->Image[i];
   }
}

void CMainFrame::OnUpdateParzen3dVisualiserlespacededcision(CCmdUI* pCmdUI) 
{
   //	pCmdUI->Enable(EspValid);
}
/*
void CMainFrame::OnParzen3dEspaceprobglobal() 
{
SetCursor(LoadCursor(NULL,IDC_WAIT));
CFile fich;
unsigned long tampon;
unsigned __int64 shift;
unsigned long pos;
double compt;
CProgressDlg Dlg;
CString ch;

for(long i=0;i<16777216;i++) {EspProb1[i]=0;EspProb2[i]=0;EspProb3[i]=0;EspProb4[i]=0;}

CFileDialog * FileData = new CFileDialog(TRUE,"*.esp",NULL,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT|OFN_ALLOWMULTISELECT ,
"Mesures proba(*.esp)|*.esp||",this);

char *m_name=new char[1000];
FileData->m_ofn.lpstrFile=m_name;
FileData->m_ofn.nMaxFile=1000;
if (FileData->DoModal()!=IDOK){ return;}
POSITION FilePos;
FilePos=FileData->GetStartPosition();

/*	while(FilePos!=NULL)
{
ch+=FileData->GetNextPathName(FilePos);
ch+="\n";
}
MessageBox(ch);
//

Dlg.Create();
ch.Format("loading classe 0...");
Dlg.SetWindowText(ch);

if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();return;}

Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
Dlg.SetStep(1);

compt=0;
/*	for(i=0;i<(long)fich.GetLength();i+=8)
{
fich.Read(&pos,4);
fich.Read(&tampon,4);
if(tampon>65535) EspProb[pos]=65535;
else EsProb[pos]=tampon;
if(!(i%(100000)))
{
Dlg.StepIt();
if(Dlg.CheckCancelButton()) {return;}
}
}
//	fich.Close();
ch.Format("tot=%1.16e",compt);
//	MessageBox(ch);
Dlg.DestroyWindow();
Dlg.Create();
ch.Format("loading classe 1...");
Dlg.SetWindowText(ch);

if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();return;}
compt=0;

Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
Dlg.SetStep(1);
for( i=0;i<(long)fich.GetLength();i+=8)
{
fich.Read(&pos,4);
fich.Read(&tampon,4);
EspProb1[pos]=tampon;
//	if(tampon>0xFFFF) EspProb1[pos]=0xFFFF;
//	else EspProb1[pos]=(unsigned __int16)tampon;
if(!(i%(100000)))
{
Dlg.StepIt();
if(Dlg.CheckCancelButton()) {return;}
}
}
fich.Close();
ch.Format("tot=%1.16e",compt);
//	MessageBox(ch);
compt=0;
Dlg.DestroyWindow();
Dlg.Create();
ch.Format("loading classe 2...");
Dlg.SetWindowText(ch);

if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();return;}

Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
Dlg.SetStep(1);
for(i=0;i<(long)fich.GetLength();i+=8)
{
fich.Read(&pos,4);
fich.Read(&tampon,4);
EspProb2[pos]=tampon;
//	if(tampon>65535) EspProb2[pos]=0xFFFF;
//	else EspProb2[pos]+=(unsigned __int16)tampon;
if(!(i%(100000)))
{
Dlg.StepIt();
if(Dlg.CheckCancelButton()) {return;}
}
}
fich.Close();
ch.Format("tot=%1.16e",compt);
//	MessageBox(ch);
compt=0;
Dlg.DestroyWindow();
Dlg.Create();
ch.Format("loading classe 3...");
Dlg.SetWindowText(ch);

if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();return;}

Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
Dlg.SetStep(1);
for(i=0;i<(long)fich.GetLength();i+=8)
{
fich.Read(&pos,4);
fich.Read(&tampon,4);
EspProb3[pos]=tampon;
//	shift=tampon;
//	if(tampon>65535) EspProb3[pos]+=0xFFFF;
//	else EspProb3[pos]+=(unsigned __int16)tampon;
if(!(i%(100000)))
{
Dlg.StepIt();
if(Dlg.CheckCancelButton()) {return;}
}
}
fich.Close();

ch.Format("tot=%1.16e",compt);
//	MessageBox(ch);
compt=0;
Dlg.DestroyWindow();
Dlg.Create();
ch.Format("loading classe 4...");
Dlg.SetWindowText(ch);

if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();return;}

Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
Dlg.SetStep(1);

for(i=0;i<(long)fich.GetLength();i+=8)
{
fich.Read(&pos,4);
fich.Read(&tampon,4);
EspProb4[pos]=tampon;
//	shift=tampon;
//	if(tampon>65535) EspProb4[pos]+=0xFFFF;
//	else EspProb4[pos]+=(unsigned __int16)tampon;
if(!(i%(100000)))
{
Dlg.StepIt();
if(Dlg.CheckCancelButton()) {return;}
}
}
fich.Close();
ch.Format("tot=%1.16e",compt);
//	MessageBox(ch);
Dlg.DestroyWindow();

/*	int max=256;
CTraitImageWnd *toto=new CTraitImageWnd(256,256,FLOAT32);
toto->Create("projection", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,NULL);
//	BOOL compt;
CTraitImageWnd *totoseq=new CTraitImageWnd(256,256,NB,TRUE);
totoseq->Create("projection en couches", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,NULL);
totoseq->AlloueSequence(&max);

for(long j=0;j<max;j++)
{
for(i=0;i<65536;i++) 
{
//	compt=FALSE;
pos=i+65536*j;
if(espace[pos]!=0) toto->ImageFloat[i]=(float)espace[i+65536*j];
else toto->ImageFloat[i]=0;
}
toto->FloatToImage();
for(i=0;i<65536;i++) totoseq->Bpile[j][i]=toto->Image[i];
}
toto->DestroyWindow();	
delete espace;*/
/*	delete FileData;
//	EspValid=TRUE;
}*/

void CMainFrame::OnParzen3dMesrurerecouvrementsaveccl0() 
{
   SetCursor(LoadCursor(NULL,IDC_WAIT));
   CFile fich;
   unsigned long tampon;
   unsigned long* espace=new unsigned long[16777216];
   unsigned long pos;
   double compt,compt1,compt2,compt3,compt4,comptot;
   CProgressDlg Dlg;
   CString ch;
   long i;
   for( i=0;i<16777216;i++) {espace[i]=0;EspClas[i]=0;}

   CFileDialog * FileData = new CFileDialog(TRUE,"*.esp",NULL,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT|OFN_ALLOWMULTISELECT ,
      "Mesures proba(*.esp)|*.esp||",this);

   char *m_name=new char[10000];
   FileData->m_ofn.lpstrFile=m_name;
   FileData->m_ofn.nMaxFile=10000;
   if (FileData->DoModal()!=IDOK){delete espace; return;}
   POSITION FilePos;
   FilePos=FileData->GetStartPosition();


   Dlg.Create();
   ch.Format("loading classe 0...");
   Dlg.SetWindowText(ch);

   if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();delete espace;return;}

   Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
   Dlg.SetStep(1);

   compt=0;

   for(i=0;i<(long)fich.GetLength();i+=8)
   {
      fich.Read(&pos,4);
      fich.Read(&espace[pos],4);
      compt+=espace[pos];
      EspClas[pos]=1;
      if(!(i%(100000)))
      {
         Dlg.StepIt();
         if(Dlg.CheckCancelButton()) {delete espace;return;}
      }
   }
   fich.Close();
   Dlg.DestroyWindow();


   Dlg.Create();
   ch.Format("loading classe 1...");
   Dlg.SetWindowText(ch);

   if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();delete espace;return;}

   compt1=0;
   Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
   Dlg.SetStep(1);
   for(i=0;i<(long)fich.GetLength();i+=8)
   {
      fich.Read(&pos,4);
      fich.Read(&tampon,4);
      if(espace[pos]){compt1+=min(tampon,espace[pos]);EspClas[pos]=2;}

      if(!(i%(100000)))
      {
         Dlg.StepIt();
         if(Dlg.CheckCancelButton()) {delete espace;return;}
      }
   }
   fich.Close();
   Dlg.DestroyWindow();

   Dlg.Create();
   ch.Format("loading classe 2...");
   Dlg.SetWindowText(ch);

   if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();delete espace;return;}

   compt2=0;
   comptot=compt1;
   Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
   Dlg.SetStep(1);
   for(i=0;i<(long)fich.GetLength();i+=8)
   {
      fich.Read(&pos,4);
      fich.Read(&tampon,4);
      if(espace[pos]) {tampon=min(tampon,espace[pos]);compt2+=tampon;if(EspClas[pos]==1) comptot+=tampon;EspClas[pos]=3;}

      if(!(i%(100000)))
      {
         Dlg.StepIt();
         if(Dlg.CheckCancelButton()) {delete espace;return;}
      }
   }
   fich.Close();

   Dlg.DestroyWindow();
   Dlg.Create();
   ch.Format("loading classe 3...");
   Dlg.SetWindowText(ch);

   if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();delete espace;return;}

   compt3=0;
   Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
   Dlg.SetStep(1);
   for(i=0;i<(long)fich.GetLength();i+=8)
   {
      fich.Read(&pos,4);
      fich.Read(&tampon,4);
      if(espace[pos]) {tampon=min(tampon,espace[pos]);compt3+=tampon;if(EspClas[pos]==1) comptot+=tampon;EspClas[pos]=4;}

      if(!(i%(100000)))
      {
         Dlg.StepIt();
         if(Dlg.CheckCancelButton()) {delete espace;return;}
      }
   }
   fich.Close();


   Dlg.DestroyWindow();
   Dlg.Create();
   ch.Format("loading classe 4...");
   Dlg.SetWindowText(ch);

   if(!fich.Open(FileData->GetNextPathName(FilePos),CFile::modeRead)) {Dlg.DestroyWindow();delete espace;return;}

   compt4=0;
   Dlg.SetRange(0,1+(int)(fich.GetLength()/100000));
   Dlg.SetStep(1);

   for(i=0;i<(long)fich.GetLength();i+=8)
   {
      fich.Read(&pos,4);
      fich.Read(&tampon,4);
      if(espace[pos]) {tampon=min(tampon,espace[pos]);compt4+=tampon;if(EspClas[pos]==1) comptot+=tampon;EspClas[pos]=5;}
      if(!(i%(100000)))
      {
         Dlg.StepIt();
         if(Dlg.CheckCancelButton()) {delete espace;return;}
      }
   }
   fich.Close();
   Dlg.DestroyWindow();
   ch.Format("total c0 =%1.10e \ntotal c1 =%1.10e\ntotal c2 =%1.10e\ntotal c3 =%1.10e\ntotal c4 =%1.10e\n\nrecouvrements:\n\tc0-c1 = %3.3f %%\n\tc0-c2 = %3.3f %%\n\tc0-c3 = %3.3f %%\n\tc0-c4 = %3.3f %%\n\tc0-ci = %3.3f %%",compt,compt1,compt2,compt3,compt4,(compt1/compt)*100,(compt2/compt)*100,(compt3/compt)*100,(compt4/compt)*100,(comptot/compt)*100);
   MessageBox(ch);

   delete espace;
   delete FileData;		
}

void CMainFrame::OnParzen3dMesuredevolumedesclasses() 
{
   long i;
   long compt=0,volume[256];
   for( i=0;i<256;i++) volume[i]=0;
   for( i=0;i<256*256*256;i++) 
      volume[EspClas[i]]++;
   CString ch,ch0;

   ch="Volume des classes :\n\n";
   for(i=0;i<256;i++) if(volume[i]) {ch0.Format("%d = %d\t%Vtot=%2.3f\n",i,volume[i],(100.0*(float)volume[i])/16777216);ch+=ch0;}
   MessageBox(ch,"volumes");

}

void CMainFrame::OnUpdateParzen3dMesuredevolumedesclasses(CCmdUI* pCmdUI) 
{
   pCmdUI->Enable(EspValid);
}

void CMainFrame::OnTestOuvrir() 
{
   // TODO: Add your command handler code here
   long Larg, Haut,l,c;
   int par,nbmes=0;
   int pos;
   float val;
   CFile Fichier;
   CStdioFile FicData("defautimage.dat",CFile::modeCreate|CFile::modeWrite);
   CString NomRep="C:\\johel\\thomson\\filaments\\defautsbmp\\";
   CString NomFich = NomRep+"Filament_0001_0001_Corps.bmp";

   CFileFind finder;
   BOOL bWorking = finder.FindFile(NomRep+"*.bmp");
   while (bWorking)
   {
      bWorking = finder.FindNextFile();
      nbmes++;
   }
   CString Temp;
   Temp.Format("%d\n",2560);
   FicData.WriteString(Temp);
   Temp.Format("%d\n",nbmes);
   FicData.WriteString(Temp);

   bWorking = finder.FindFile(NomRep+"*.bmp");
   while (bWorking)
   {
      bWorking = finder.FindNextFile();
      NomFich= finder.GetFilePath();

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
         par=0;  
         for (l=32;l<192;l+=4)
         {
            for (c=0;c<Larg;c+=4)
            {
               val=(float)pImageWnd->Image[Larg*l+c];
               val=val/255.0;
               Temp.Format("%f ",val);
               FicData.WriteString(Temp);
               par++;
            }
         }
         FicData.WriteString("\n");
         pImageWnd->Create(NomFich.Left(NomFich.GetLength()-4), WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,this);

         pImageWnd->DestroyWindow(); 
      }
      else MessageBox("Ce fichier n'existe pas","Message",MB_OK|MB_ICONEXCLAMATION);
   }
   FicData.Close();

}
//------------------------------------------
//	Calcul les moments de Zernike à partir 
//	de masques
//  les fichiers masques doivent être dans 
//	'./MasqueExact/...'
//	'./MasqueBinaire/...'
//------------------------------------------
void CMainFrame::OnZernikeBatchmoments() 
{
   int Largeur = 256;
   int Hauteur = 256;
   int NombreMasque = 37;
   int i=0,j=0,k=0,w=0,indchaine=0;;
   int Taille = Largeur*Hauteur;
   CStdioFile FichierMasque;
   CString NomFichier;
   CString ChaineLectureFichier;
   CString ChaineRecopie;
   CString Ouverture ="";
   CString ExtensionExact = ".mse";
   CString ExtensionApprox = ".msd";
   double temp = 0.0f;
   //CString NomDossierExact = "C:\\users\\Cedric Lemaitre\\code curvilineaire\\WINSEG11\\MasqueExact\\";
   CString NomDossierExact = ".\\MasqueExact\\";
   CString NomDossierDither = ".\\MasqueApprox\\";
   long CompTaille = 0; 


   /*double ** TabMasqueBinaire = (double **)malloc(sizeof(double*)*NombreMasque);
   for(i=0;i<NombreMasque;i++)
   {
   TabMasqueBinaire[i]=(double*)malloc(sizeof(double)*Hauteur*Largeur);
   }*/


   double ** TabMasqueExact = (double **)malloc(sizeof(double*)*NombreMasque);
   for(i=0;i<NombreMasque;i++)
   {
      TabMasqueExact[i]=(double*)malloc(sizeof(double)*Hauteur*Largeur);
   }

   for(k=0;k<NombreMasque;k++)
      for(j=0;j<Taille;j++)
      {
         TabMasqueExact[k][j]=0.0f;
      }

      /*for(k=0;k<NombreMasque;k++)
      for(j=0;j<Taille;j++)
      {
      TabMasqueBinaire[k][j]=0.0f;
      }*/


      /*for(k=0;k<NombreMasque;k++)
      {
      NomFichier.Format("%d",k+1);
      Ouverture = NomDossierExact + NomFichier + ExtensionExact;

      if (FichierMasque.Open(Ouverture,CFile::modeRead))
      {
      for(i=0;i<Hauteur;i++)
      {
      FichierMasque.ReadString(ChaineLectureFichier);
      for(j=0;j<Largeur;j++)
      {
      if (j==0)
      {
      for(w=0;w<11;w++)
      {
      ChaineRecopie += ChaineLectureFichier[w];
      }
      temp= atof(ChaineRecopie);
      TabMasqueExact[k][CompTaille] = temp;
      ChaineRecopie = "";
      }
      else
      {
      indchaine=j*12+1;
      for(w=0;w<11;w++)
      {
      ChaineRecopie += ChaineLectureFichier[w];
      }
      temp= atof(ChaineRecopie);
      TabMasqueExact[k][CompTaille] = temp;
      ChaineRecopie = "";
      }
      CompTaille++;
      }
      }
      FichierMasque.Close();
      Ouverture = "";
      }
      }*/


      for(k=0;k<NombreMasque;k++)
      {
         NomFichier.Format("%d",k+1);
         Ouverture = NomDossierExact + NomFichier + ExtensionApprox;

         if (FichierMasque.Open(Ouverture,CFile::modeRead))
         {
            for(i=0;i<Hauteur;i++)
            {
               FichierMasque.ReadString(ChaineLectureFichier);
               for(j=0;j<Largeur;j++)
               {
                  if (j==0)
                  {
                     for(w=0;w<11;w++)
                     {
                        ChaineRecopie += ChaineLectureFichier[w];
                     }
                     temp= atof(ChaineRecopie);
                     TabMasqueExact[k][CompTaille] = temp;
                     ChaineRecopie = "";
                  }
                  else
                  {
                     indchaine=j*12+1;
                     for(w=0;w<11;w++)
                     {
                        ChaineRecopie += ChaineLectureFichier[w];
                     }
                     temp= atof(ChaineRecopie);
                     TabMasqueExact[k][CompTaille] = temp;
                     ChaineRecopie = "";
                  }
                  CompTaille++;
               }
            }
            FichierMasque.Close();
            Ouverture = "";
         }
      }

      //----------------------------------------
      // Ouverture du fichier txt + lecture des
      // images
      // calcul des descripteurs approx
      //----------------------------------------




      //----------------------------------------
      //		Liberation de la memoire
      //----------------------------------------


      for(i=0;i<NombreMasque;i++)
      {
         free(TabMasqueExact[i]);
      }

      free(TabMasqueExact);


}




void CMainFrame::OnEtudecurvi1dGenerationfichierprofilPorte() 
{
   CFile Fichier;
   int i=0;

   int LargeurPorte = 48;
   //int LargeurFourier = 16;
   int LongueurSignal = 256;

   int DemiLargeurPorte = LargeurPorte/2;
   int DemiLongueurSignal = LongueurSignal/2;

   int * SignalPorte = (int *) malloc(sizeof(int)*LongueurSignal);

   int CondGauche = DemiLongueurSignal-DemiLargeurPorte;
   int CondDroite = DemiLongueurSignal+DemiLargeurPorte;
   for (i=0;i<LongueurSignal;i++)
   {
      if ((i>CondGauche)&(i<CondDroite))
      {
         SignalPorte[i]=255;
      }
      else
      {
         SignalPorte[i]=0;
      }
   }

   CStdioFile ExportTxt;
   CString ChaineExport;

   if(Fichier.Open("PorteEtude1DL48.bin",CFile::modeWrite|CFile::modeCreate)&&ExportTxt.Open("ExportXLSPorteEtude1DL48.txt",CFile::modeWrite|CFile::modeCreate))
   {
      for (i=0;i<LongueurSignal;i++)
      {
         ChaineExport.Format("%d\n",SignalPorte[i]);
         ExportTxt.WriteString(ChaineExport);

      }
      Fichier.Write(SignalPorte,sizeof(int)*LongueurSignal);
      Fichier.Close();
      ExportTxt.Close();
   }
   else
   {
      MessageBox("Ouverture du fichier <<PorteEtude1DL48.bin>> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
   }

   free(SignalPorte);	
}

void CMainFrame::OnEtudecurvi1dEtudeavecobjetEtudeporte() 
{
   // Relecture Fichier porte


   CFile Fichier;
   CIma ObjetCIma;
   int i=0,j=0,k=0,l=0;

   int LargeurPorte = 48;
   int LargeurFourier = 16; // doit être une puissance de deux
   int LongueurSignal = 256;

   int DemiLargeurFourier = LargeurFourier/2;

   int RacineLargeurFourier = (int)sqrt((double)LargeurFourier);

   int DemiLargeurPorte = LargeurPorte/2;
   int DemiLongueurSignal = LongueurSignal/2;

   int * SignalPorte = (int *) malloc(sizeof(int));

   int CondGauche = DemiLongueurSignal-DemiLargeurPorte;
   int CondDroite = DemiLongueurSignal+DemiLargeurPorte;

   double * FGaucheRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FDroiteRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreGRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreDRE = (double *) malloc(sizeof(double)*LargeurFourier);


   double * FGaucheIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FDroiteIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreGIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreDIM = (double *) malloc(sizeof(double)*LargeurFourier);



   for(i=0;i<LargeurFourier;i++)
   {

      FGaucheRE[i]=0.0f;
      FDroiteRE[i]=0.0f;
      FCentreGRE[i]=0.0f;
      FCentreDRE[i]=0.0f;

      FGaucheIM[i]=0.0f;
      FDroiteIM[i]=0.0f;
      FCentreGIM[i]=0.0f;
      FCentreDIM[i]=0.0f;

   }
   int IptInit;
   int IptFinal;
   int LargeurGabaritMin;
   int LargeurGabaritMaxi;
   int IptInitMax = LongueurSignal-LargeurFourier-LargeurPorte;
   int IptFinalMax = LongueurSignal-LargeurFourier-1; //attention 
   int IptInitMin = LargeurFourier;
   int IptFinalMin = LargeurFourier + LargeurPorte;

   int cpt=0;
   int cptTabG=0;
   int cptTabD=0;
   int cptSection=0;
   int cptBord=0;

   CStdioFile FichierStdIO;
   CString ChaineDist;

   int TailleTabDist = (IptInitMax-IptInitMin)*(IptFinalMax-IptFinalMin);

   double * TabDistanceGauche = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabDistanceDroite = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabHomogSection = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabBorSim = (double *) malloc(sizeof(double)*TailleTabDist);

   double * MuLTabDis = (double *) malloc(sizeof(double)*TailleTabDist);

   double * TabMaxVisu = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

   //double * TabStoreTemp = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

   double * TabStoreTemp = (double *)malloc(sizeof(double) *(LongueurSignal));

   int * TabLargeur = (int *)malloc(sizeof(int) *(IptInitMax-IptInitMin));



   for (i=0;i<(LongueurSignal);i++)
   {
      TabStoreTemp[i] = 0.0f;
   }

   for (i=0;i<(IptInitMax-IptInitMin);i++)
   {
      TabMaxVisu[i] = 0.0f;
      TabLargeur[i] = 0;
   }

   for(i=0;i<LargeurFourier;i++)
   {

      TabDistanceGauche[i] = 0.0f;
      TabDistanceDroite[i] = 0.0f;


   }

   for(i=0;i<TailleTabDist;i++)
   {
      MuLTabDis[i] = 0.0f;
      TabDistanceGauche[i]=0.0f;
      TabDistanceDroite[i]=0.0f;
      TabHomogSection[i]=0.0f;
      TabBorSim[i]=0.0f;


   }

   if(Fichier.Open("PorteEtude1DL48.bin",CFile::modeRead))
   {

      Fichier.Read(SignalPorte,sizeof(int)*LongueurSignal);
      Fichier.Close();
   }
   else
   {
      MessageBox("Ouverture du fichier <<PorteEtude1DL48.bin>> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
   }



   //------------------------
   //			RAZ 
   //------------------------

   cptTabG=0;
   cptTabD=0;
   cptSection=0;
   cptBord=0;


   //----------------------------------------------------
   // recopie des tab pour fourier + calcul distance
   //----------------------------------------------------
   for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
   {
      for(IptFinal=(IptInit+(LargeurFourier*2)+1);IptFinal<IptFinalMax;IptFinal++) //idem modif 21/11/2006
      {
         //copie Fgauche 

         cpt =0;
         for (i = IptInit-LargeurFourier;i<IptInit;i++)
         {
            FGaucheRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)
         {
            FDroiteRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptInit;i<IptInit+LargeurFourier;i++)
         {
            FCentreGRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptFinal-LargeurFourier+1;i<IptFinal+1;i++)
         {
            FCentreDRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt=0;



         // calcul des FFT
         ObjetCIma.FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);


         TabDistanceGauche[cptTabG] = ObjetCIma.HyperDistanceComplexe(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,DemiLargeurFourier);
         TabDistanceDroite[cptTabD] = ObjetCIma.HyperDistanceComplexe(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,DemiLargeurFourier);
         TabHomogSection[cptSection] = ObjetCIma.HyperDistanceComplexe(FCentreGRE,FCentreGIM,FCentreDRE,FCentreDIM,DemiLargeurFourier);
         TabBorSim[cptBord]=ObjetCIma.HyperDistanceComplexe(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,DemiLargeurFourier);

         cptTabG++;
         cptTabD++;
         cptSection++;
         cptBord++;
      }
   }





   double max = 0.0f;
   int PosMax=0;

   double MaxGauche = 0.0f,MaxDroite = 0.0f,MaxHomoSec = 0.0f, MaxBordSim = 0.0f;



   //-----------------------------------------------------------------------------------------------------
   // recherche des max pour l'inversion des fonctions de probalilité "section homogène" et bord similaire
   //-----------------------------------------------------------------------------------------------------
   for(i=0;i<TailleTabDist;i++)
   {


      if(TabHomogSection[i]>MaxHomoSec)
      {
         MaxHomoSec = TabHomogSection[i];
      }

      if(TabBorSim[i]>MaxBordSim)
      {
         MaxBordSim = TabBorSim[i];
      }

      if(TabDistanceGauche[i]>MaxGauche)
      {
         MaxGauche = TabDistanceGauche[i];
      }

      if(TabDistanceDroite[i]>MaxDroite)
      {
         MaxDroite = TabDistanceDroite[i];
      }

   }




   //--------------------------------------------------------
   //			Fonction export + calcul des couts
   //--------------------------------------------------------

   i=0; // compteur
   j=0;
   int AccuVisu=0;
   double MaxVisu = 0.0f;
   double MaxIter = 0.0f;

   CStdioFile FichierVisu;
   CString ChaineVisu;
   if ((FichierStdIO.Open("ExportBordSimilaire.txt",CFile::modeCreate|CFile::modeWrite)) & (FichierVisu.Open("ExportVisu.txt",CFile::modeCreate|CFile::modeWrite)) )
   {j=0;
   for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
   {

      MaxIter = 0.0f;
      PosMax = IptInit;
      for(IptFinal=(IptInit+(LargeurFourier*2)+1);IptFinal<IptFinalMax;IptFinal++) // modif 21/11/2006 rajout pour limiter les surbalayage
      {

         MuLTabDis[i] = ((MaxHomoSec-TabHomogSection[i])/MaxHomoSec) * (TabDistanceGauche[i]/MaxGauche) * (TabDistanceDroite[i]/MaxDroite)*(MaxBordSim-TabBorSim[i])/MaxBordSim;
         ChaineDist.Format("%f\t",MuLTabDis[i]);
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptInit));
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptFinal));
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptFinal-IptInit));
         FichierStdIO.WriteString(ChaineDist);
         FichierStdIO.WriteString("\n");

         //TabStoreTemp[j] = MuLTabDis[i];

         if(MuLTabDis[i]>MaxIter)
         {
            MaxIter = MuLTabDis[i];
            PosMax = IptFinal;
         }

         i++;
      }//fin IptFinal
      j= (PosMax+IptInit)/2;
      if(TabStoreTemp[j]<MaxIter)TabStoreTemp[j]=MaxIter;// a sauver

      //j++;
      //rechercdu max pour la ligne stockage dans le tableau pour la visu
   }//fin IptInit

   for(j=0;j<(LongueurSignal);j++)
   {
      ChaineVisu.Format("%f\n",TabStoreTemp[j]);
      FichierVisu.WriteString(ChaineVisu);
      /*if(TabStoreTemp[j]>MaxVisu)
      {
      MaxVisu=TabStoreTemp[j];
      TabLargeur[AccuVisu] = ((j+IptFinalMin)-IptInit);
      }*/
   }
   FichierVisu.Close();
   FichierStdIO.Close();
   }


   cpt = 0;
   CString ChaineMessage;
   CString ChaineTemp;
   int PosStart=0,PosEnd=0;
   ChaineMessage = "Config maximale pour IptInit = ";

   for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
   {
      for(IptFinal=IptFinalMin;IptFinal<IptFinalMax;IptFinal++)
      {
         if (cpt == PosMax)
         {
            PosStart=IptInit;
            PosEnd=IptFinal;
            /*ChaineTemp.Format("%d",IptInit);
            ChaineMessage += ChaineTemp;
            ChaineMessage += " et IPtFinal = ";
            ChaineTemp.Format("%d",IptFinal);
            ChaineMessage += ChaineTemp;*/


         }
         cpt++;
      }
   }


   //MessageBox(NULL,ChaineMessage,"Information",MB_OK);



   free(FGaucheRE); 
   free(FDroiteRE);//---> debordement dans les boucles
   free(FCentreGRE); 
   free(FCentreDRE); 


   free(FGaucheIM);
   free(FDroiteIM);
   free(FCentreGIM);
   free(FCentreDIM);

   free(TabDistanceGauche);
   free(TabDistanceDroite);

   free(TabHomogSection);
   free(TabBorSim);


   free(SignalPorte);
}

void CMainFrame::OnEtudecurvi1dGenerationfichierprofilPorte32() 
{
   // TODO: Add your command handler code here
   CFile Fichier;
   int i=0;

   int LargeurPorte = 32;
   //int LargeurFourier = 16;
   int LongueurSignal = 256;

   int DemiLargeurPorte = LargeurPorte/2;
   int DemiLongueurSignal = LongueurSignal/2;

   int * SignalPorte = (int *) malloc(sizeof(int)*LongueurSignal);

   int CondGauche = DemiLongueurSignal-DemiLargeurPorte;
   int CondDroite = DemiLongueurSignal+DemiLargeurPorte;
   for (i=0;i<LongueurSignal;i++)
   {
      if ((i>CondGauche)&(i<CondDroite))
      {
         SignalPorte[i]=255;
      }
      else
      {
         SignalPorte[i]=0;
      }
   }

   CStdioFile ExportTxt;
   CString ChaineExport;
   if(Fichier.Open("PorteEtude1DL32.bin",CFile::modeWrite|CFile::modeCreate)&&ExportTxt.Open("ExportXLSPorteEtude1DL32.txt",CFile::modeWrite|CFile::modeCreate))
   {
      for (i=0;i<LongueurSignal;i++)
      {
         ChaineExport.Format("%d\n",SignalPorte[i]);
         ExportTxt.WriteString(ChaineExport);

      }
      Fichier.Write(SignalPorte,sizeof(int)*LongueurSignal);
      Fichier.Close();
      ExportTxt.Close();
   }
   else
   {
      MessageBox("Ouverture du fichier <<PorteEtude1DL32.bin>> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
   }

   free(SignalPorte);	
}

void CMainFrame::OnEtudecurvi1dGenerationfichierprofilPorte16() 
{
   // TODO: Add your command handler code here
   CFile Fichier;
   int i=0;

   int LargeurPorte = 16;
   //int LargeurFourier = 16;
   int LongueurSignal = 256;

   int DemiLargeurPorte = LargeurPorte/2;
   int DemiLongueurSignal = LongueurSignal/2;

   int * SignalPorte = (int *) malloc(sizeof(int)*LongueurSignal);

   int CondGauche = DemiLongueurSignal-DemiLargeurPorte;
   int CondDroite = DemiLongueurSignal+DemiLargeurPorte;
   for (i=0;i<LongueurSignal;i++)
   {
      if ((i>CondGauche)&(i<CondDroite))
      {
         SignalPorte[i]=255;
      }
      else
      {
         SignalPorte[i]=0;
      }
   }

   CStdioFile ExportTxt;
   CString ChaineExport;

   if(Fichier.Open("PorteEtude1DL16.bin",CFile::modeWrite|CFile::modeCreate)&&ExportTxt.Open("ExportXLSPorteEtude1DL16.txt",CFile::modeWrite|CFile::modeCreate))
   {
      for (i=0;i<LongueurSignal;i++)
      {
         ChaineExport.Format("%d\n",SignalPorte[i]);
         ExportTxt.WriteString(ChaineExport);

      }
      Fichier.Write(SignalPorte,sizeof(int)*LongueurSignal);
      Fichier.Close();
      ExportTxt.Close();
   }
   else
   {
      MessageBox("Ouverture du fichier <<PorteEtude1DL16.bin>> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
   }

   free(SignalPorte);	

}

void CMainFrame::OnEtudecurvi1dGenerationfichierprofilPorte12() 
{
   // TODO: Add your command handler code here
   CFile Fichier;
   int i=0;

   int LargeurPorte = 12;
   //int LargeurFourier = 16;
   int LongueurSignal = 256;

   int DemiLargeurPorte = LargeurPorte/2;
   int DemiLongueurSignal = LongueurSignal/2;

   int * SignalPorte = (int *) malloc(sizeof(int)*LongueurSignal);

   int CondGauche = DemiLongueurSignal-DemiLargeurPorte;
   int CondDroite = DemiLongueurSignal+DemiLargeurPorte;
   for (i=0;i<LongueurSignal;i++)
   {
      if ((i>CondGauche)&(i<CondDroite))
      {
         SignalPorte[i]=255;
      }
      else
      {
         SignalPorte[i]=0;
      }
   }
   CStdioFile ExportTxt;
   CString ChaineExport;
   if(Fichier.Open("PorteEtude1DL12.bin",CFile::modeWrite|CFile::modeCreate)&&ExportTxt.Open("ExportXLSPorteEtude1DL12.txt",CFile::modeWrite|CFile::modeCreate))
   {
      for (i=0;i<LongueurSignal;i++)
      {
         ChaineExport.Format("%d\n",SignalPorte[i]);
         ExportTxt.WriteString(ChaineExport);

      }
      Fichier.Write(SignalPorte,sizeof(int)*LongueurSignal);
      Fichier.Close();
      ExportTxt.Close();
   }
   else
   {
      MessageBox("Ouverture du fichier <<PorteEtude1DL12.bin>> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
   }

   free(SignalPorte);	
}


//-----------------------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------------------
void CMainFrame::OnEtudecurvi1dEtudeavecobjetEtudeporte32() 
{
   // Relecture Fichier porte


   CFile Fichier;
   CIma ObjetCIma;
   int i=0,j=0,k=0,l=0;

   int LargeurPorte = 48;
   int LargeurFourier = 16; // doit être une puissance de deux
   int LongueurSignal = 256;

   int DemiLargeurFourier = LargeurFourier/2;

   int RacineLargeurFourier = (int)sqrt((double)LargeurFourier);

   int DemiLargeurPorte = LargeurPorte/2;
   int DemiLongueurSignal = LongueurSignal/2;

   int * SignalPorte = (int *) malloc(sizeof(int));

   int CondGauche = DemiLongueurSignal-DemiLargeurPorte;
   int CondDroite = DemiLongueurSignal+DemiLargeurPorte;

   double * FGaucheRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FDroiteRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreGRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreDRE = (double *) malloc(sizeof(double)*LargeurFourier);


   double * FGaucheIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FDroiteIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreGIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreDIM = (double *) malloc(sizeof(double)*LargeurFourier);



   for(i=0;i<LargeurFourier;i++)
   {

      FGaucheRE[i]=0.0f;
      FDroiteRE[i]=0.0f;
      FCentreGRE[i]=0.0f;
      FCentreDRE[i]=0.0f;

      FGaucheIM[i]=0.0f;
      FDroiteIM[i]=0.0f;
      FCentreGIM[i]=0.0f;
      FCentreDIM[i]=0.0f;

   }
   int IptInit;
   int IptFinal;
   int LargeurGabaritMin;
   int LargeurGabaritMaxi;
   int IptInitMax = LongueurSignal-LargeurFourier-LargeurPorte;
   int IptFinalMax = LongueurSignal-LargeurFourier-1; //attention 
   int IptInitMin = LargeurFourier;
   int IptFinalMin = LargeurFourier + LargeurPorte;

   int cpt=0;
   int cptTabG=0;
   int cptTabD=0;
   int cptSection=0;
   int cptBord=0;

   CStdioFile FichierStdIO;
   CString ChaineDist;

   int TailleTabDist = (IptInitMax-IptInitMin)*(IptFinalMax-IptFinalMin);

   double * TabDistanceGauche = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabDistanceDroite = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabHomogSection = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabBorSim = (double *) malloc(sizeof(double)*TailleTabDist);

   double * MuLTabDis = (double *) malloc(sizeof(double)*TailleTabDist);

   double * TabMaxVisu = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

   //double * TabStoreTemp = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

   double * TabStoreTemp = (double *)malloc(sizeof(double) *(LongueurSignal));

   int * TabLargeur = (int *)malloc(sizeof(int) *(IptInitMax-IptInitMin));



   for (i=0;i<(LongueurSignal);i++)
   {
      TabStoreTemp[i] = 0.0f;
   }

   for (i=0;i<(IptInitMax-IptInitMin);i++)
   {
      TabMaxVisu[i] = 0.0f;
      TabLargeur[i] = 0;
   }

   for(i=0;i<LargeurFourier;i++)
   {

      TabDistanceGauche[i] = 0.0f;
      TabDistanceDroite[i] = 0.0f;


   }

   for(i=0;i<TailleTabDist;i++)
   {
      MuLTabDis[i] = 0.0f;
      TabDistanceGauche[i]=0.0f;
      TabDistanceDroite[i]=0.0f;
      TabHomogSection[i]=0.0f;
      TabBorSim[i]=0.0f;


   }

   if(Fichier.Open("PorteEtude1DL32.bin",CFile::modeRead))
   {

      Fichier.Read(SignalPorte,sizeof(int)*LongueurSignal);
      Fichier.Close();
   }
   else
   {
      MessageBox("Ouverture du fichier <<PorteEtude1DL32.bin>> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
   }



   //------------------------
   //			RAZ 
   //------------------------

   cptTabG=0;
   cptTabD=0;
   cptSection=0;
   cptBord=0;


   //----------------------------------------------------
   // recopie des tab pour fourier + calcul distance
   //----------------------------------------------------
   for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
   {
      for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) //idem modif 21/11/2006
      {
         //copie Fgauche 

         cpt =0;
         for (i = IptInit-LargeurFourier;i<IptInit;i++)
         {
            FGaucheRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)
         {
            FDroiteRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptInit;i<IptInit+LargeurFourier;i++)
         {
            FCentreGRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptFinal-LargeurFourier+1;i<IptFinal+1;i++)
         {
            FCentreDRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt=0;



         // calcul des FFT
         ObjetCIma.FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);


         TabDistanceGauche[cptTabG] = ObjetCIma.HyperDistanceComplexe(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,DemiLargeurFourier);
         TabDistanceDroite[cptTabD] = ObjetCIma.HyperDistanceComplexe(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,DemiLargeurFourier);
         TabHomogSection[cptSection] = ObjetCIma.HyperDistanceComplexe(FCentreGRE,FCentreGIM,FCentreDRE,FCentreDIM,DemiLargeurFourier);
         TabBorSim[cptBord]=ObjetCIma.HyperDistanceComplexe(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,DemiLargeurFourier);

         cptTabG++;
         cptTabD++;
         cptSection++;
         cptBord++;
      }
   }





   double max = 0.0f;
   int PosMax=0;

   double MaxGauche = 0.0f,MaxDroite = 0.0f,MaxHomoSec = 0.0f, MaxBordSim = 0.0f;



   //-----------------------------------------------------------------------------------------------------
   // recherche des max pour l'inversion des fonctions de probalilité "section homogène" et bord similaire
   //-----------------------------------------------------------------------------------------------------
   for(i=0;i<TailleTabDist;i++)
   {


      if(TabHomogSection[i]>MaxHomoSec)
      {
         MaxHomoSec = TabHomogSection[i];
      }

      if(TabBorSim[i]>MaxBordSim)
      {
         MaxBordSim = TabBorSim[i];
      }

      if(TabDistanceGauche[i]>MaxGauche)
      {
         MaxGauche = TabDistanceGauche[i];
      }

      if(TabDistanceDroite[i]>MaxDroite)
      {
         MaxDroite = TabDistanceDroite[i];
      }

   }




   //--------------------------------------------------------
   //			Fonction export + calcul des couts
   //--------------------------------------------------------

   i=0; // compteur
   j=0;
   int AccuVisu=0;
   double MaxVisu = 0.0f;
   double MaxIter = 0.0f;

   CStdioFile FichierVisu;
   CString ChaineVisu;
   if ((FichierStdIO.Open("ExportBordSimilaireL32.txt",CFile::modeCreate|CFile::modeWrite)) & (FichierVisu.Open("ExportVisuL32.txt",CFile::modeCreate|CFile::modeWrite)) )
   {j=0;
   for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
   {

      MaxIter = 0.0f;
      PosMax = IptInit;
      for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) // modif 21/11/2006 rajout pour limiter les surbalayage
      {

         MuLTabDis[i] = ((MaxHomoSec-TabHomogSection[i])/MaxHomoSec) * (TabDistanceGauche[i]/MaxGauche) * (TabDistanceDroite[i]/MaxDroite)*(MaxBordSim-TabBorSim[i])/MaxBordSim;
         ChaineDist.Format("%f\t",MuLTabDis[i]);
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptInit));
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptFinal));
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptFinal-IptInit));
         FichierStdIO.WriteString(ChaineDist);
         FichierStdIO.WriteString("\n");

         //TabStoreTemp[j] = MuLTabDis[i];

         if(MuLTabDis[i]>MaxIter)
         {
            MaxIter = MuLTabDis[i];
            PosMax = IptFinal;
         }

         i++;
      }//fin IptFinal
      j= (PosMax+IptInit)/2;
      if(TabStoreTemp[j]<MaxIter)TabStoreTemp[j]=MaxIter;// a sauver

      //j++;
      //rechercdu max pour la ligne stockage dans le tableau pour la visu
   }//fin IptInit

   for(j=0;j<(LongueurSignal);j++)
   {
      ChaineVisu.Format("%f\n",TabStoreTemp[j]);
      FichierVisu.WriteString(ChaineVisu);
      /*if(TabStoreTemp[j]>MaxVisu)
      {
      MaxVisu=TabStoreTemp[j];
      TabLargeur[AccuVisu] = ((j+IptFinalMin)-IptInit);
      }*/
   }
   FichierVisu.Close();
   FichierStdIO.Close();
   }






   free(FGaucheRE); 
   free(FDroiteRE);//---> debordement dans les boucles
   free(FCentreGRE); 
   free(FCentreDRE); 


   free(FGaucheIM);
   free(FDroiteIM);
   free(FCentreGIM);
   free(FCentreDIM);

   free(TabDistanceGauche);
   free(TabDistanceDroite);

   free(TabHomogSection);
   free(TabBorSim);


   free(SignalPorte);	
}
//---------------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------------


void CMainFrame::OnEtudecurvi1dEtudeavecobjetEtudeporte16() 
{
   // TODO: Add your command handler code here
   // Relecture Fichier porte


   CFile Fichier;
   CIma ObjetCIma;
   int i=0,j=0,k=0,l=0;

   int LargeurPorte = 48;
   int LargeurFourier = 16; // doit être une puissance de deux
   int LongueurSignal = 256;

   int DemiLargeurFourier = LargeurFourier/2;

   int RacineLargeurFourier = (int)sqrt((double)LargeurFourier);

   int DemiLargeurPorte = LargeurPorte/2;
   int DemiLongueurSignal = LongueurSignal/2;

   int * SignalPorte = (int *) malloc(sizeof(int));

   int CondGauche = DemiLongueurSignal-DemiLargeurPorte;
   int CondDroite = DemiLongueurSignal+DemiLargeurPorte;

   double * FGaucheRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FDroiteRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreGRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreDRE = (double *) malloc(sizeof(double)*LargeurFourier);


   double * FGaucheIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FDroiteIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreGIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreDIM = (double *) malloc(sizeof(double)*LargeurFourier);



   for(i=0;i<LargeurFourier;i++)
   {

      FGaucheRE[i]=0.0f;
      FDroiteRE[i]=0.0f;
      FCentreGRE[i]=0.0f;
      FCentreDRE[i]=0.0f;

      FGaucheIM[i]=0.0f;
      FDroiteIM[i]=0.0f;
      FCentreGIM[i]=0.0f;
      FCentreDIM[i]=0.0f;

   }
   int IptInit;
   int IptFinal;
   int LargeurGabaritMin;
   int LargeurGabaritMaxi;
   int IptInitMax = LongueurSignal-LargeurFourier-LargeurPorte;
   int IptFinalMax = LongueurSignal-LargeurFourier-1; //attention 
   int IptInitMin = LargeurFourier;
   int IptFinalMin = LargeurFourier + LargeurPorte;

   int cpt=0;
   int cptTabG=0;
   int cptTabD=0;
   int cptSection=0;
   int cptBord=0;

   CStdioFile FichierStdIO;
   CString ChaineDist;

   int TailleTabDist = (IptInitMax-IptInitMin)*(IptFinalMax-IptFinalMin);

   double * TabDistanceGauche = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabDistanceDroite = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabHomogSection = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabBorSim = (double *) malloc(sizeof(double)*TailleTabDist);

   double * MuLTabDis = (double *) malloc(sizeof(double)*TailleTabDist);

   double * TabMaxVisu = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

   //double * TabStoreTemp = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

   double * TabStoreTemp = (double *)malloc(sizeof(double) *(LongueurSignal));

   int * TabLargeur = (int *)malloc(sizeof(int) *(IptInitMax-IptInitMin));



   for (i=0;i<(LongueurSignal);i++)
   {
      TabStoreTemp[i] = 0.0f;
   }

   for (i=0;i<(IptInitMax-IptInitMin);i++)
   {
      TabMaxVisu[i] = 0.0f;
      TabLargeur[i] = 0;
   }

   for(i=0;i<LargeurFourier;i++)
   {

      TabDistanceGauche[i] = 0.0f;
      TabDistanceDroite[i] = 0.0f;


   }

   for(i=0;i<TailleTabDist;i++)
   {
      MuLTabDis[i] = 0.0f;
      TabDistanceGauche[i]=0.0f;
      TabDistanceDroite[i]=0.0f;
      TabHomogSection[i]=0.0f;
      TabBorSim[i]=0.0f;


   }

   if(Fichier.Open("PorteEtude1DL16.bin",CFile::modeRead))
   {

      Fichier.Read(SignalPorte,sizeof(int)*LongueurSignal);
      Fichier.Close();
   }
   else
   {
      MessageBox("Ouverture du fichier <<PorteEtude1DL16.bin>> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
   }



   //------------------------
   //			RAZ 
   //------------------------

   cptTabG=0;
   cptTabD=0;
   cptSection=0;
   cptBord=0;


   //----------------------------------------------------
   // recopie des tab pour fourier + calcul distance
   //----------------------------------------------------
   for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
   {
      for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) //idem modif 21/11/2006
      {
         //copie Fgauche 

         cpt =0;
         for (i = IptInit-LargeurFourier;i<IptInit;i++)
         {
            FGaucheRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)
         {
            FDroiteRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptInit;i<IptInit+LargeurFourier;i++)
         {
            FCentreGRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptFinal-LargeurFourier+1;i<IptFinal+1;i++)
         {
            FCentreDRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt=0;



         // calcul des FFT
         ObjetCIma.FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);


         TabDistanceGauche[cptTabG] = ObjetCIma.HyperDistanceComplexe(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,DemiLargeurFourier);
         TabDistanceDroite[cptTabD] = ObjetCIma.HyperDistanceComplexe(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,DemiLargeurFourier);
         TabHomogSection[cptSection] = ObjetCIma.HyperDistanceComplexe(FCentreGRE,FCentreGIM,FCentreDRE,FCentreDIM,DemiLargeurFourier);
         TabBorSim[cptBord]=ObjetCIma.HyperDistanceComplexe(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,DemiLargeurFourier);

         cptTabG++;
         cptTabD++;
         cptSection++;
         cptBord++;
      }
   }





   double max = 0.0f;
   int PosMax=0;

   double MaxGauche = 0.0f,MaxDroite = 0.0f,MaxHomoSec = 0.0f, MaxBordSim = 0.0f;



   //-----------------------------------------------------------------------------------------------------
   // recherche des max pour l'inversion des fonctions de probalilité "section homogène" et bord similaire
   //-----------------------------------------------------------------------------------------------------
   for(i=0;i<TailleTabDist;i++)
   {


      if(TabHomogSection[i]>MaxHomoSec)
      {
         MaxHomoSec = TabHomogSection[i];
      }

      if(TabBorSim[i]>MaxBordSim)
      {
         MaxBordSim = TabBorSim[i];
      }

      if(TabDistanceGauche[i]>MaxGauche)
      {
         MaxGauche = TabDistanceGauche[i];
      }

      if(TabDistanceDroite[i]>MaxDroite)
      {
         MaxDroite = TabDistanceDroite[i];
      }

   }




   //--------------------------------------------------------
   //			Fonction export + calcul des couts
   //--------------------------------------------------------

   i=0; // compteur
   j=0;
   int AccuVisu=0;
   double MaxVisu = 0.0f;
   double MaxIter = 0.0f;

   CStdioFile FichierVisu;
   CString ChaineVisu;
   if ((FichierStdIO.Open("ExportBordSimilaireL16.txt",CFile::modeCreate|CFile::modeWrite)) & (FichierVisu.Open("ExportVisuL16.txt",CFile::modeCreate|CFile::modeWrite)) )
   {j=0;
   for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
   {

      MaxIter = 0.0f;
      PosMax = IptInit;
      for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) // modif 21/11/2006 rajout pour limiter les surbalayage
      {

         MuLTabDis[i] = ((MaxHomoSec-TabHomogSection[i])/MaxHomoSec) * (TabDistanceGauche[i]/MaxGauche) * (TabDistanceDroite[i]/MaxDroite)*(MaxBordSim-TabBorSim[i])/MaxBordSim;
         ChaineDist.Format("%f\t",MuLTabDis[i]);
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptInit));
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptFinal));
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptFinal-IptInit));
         FichierStdIO.WriteString(ChaineDist);
         FichierStdIO.WriteString("\n");

         //TabStoreTemp[j] = MuLTabDis[i];

         if(MuLTabDis[i]>MaxIter)
         {
            MaxIter = MuLTabDis[i];
            PosMax = IptFinal;
         }

         i++;
      }//fin IptFinal
      j= (PosMax+IptInit)/2;
      if(TabStoreTemp[j]<MaxIter)TabStoreTemp[j]=MaxIter;// a sauver

      //j++;
      //rechercdu max pour la ligne stockage dans le tableau pour la visu
   }//fin IptInit

   for(j=0;j<(LongueurSignal);j++)
   {
      ChaineVisu.Format("%f\n",TabStoreTemp[j]);
      FichierVisu.WriteString(ChaineVisu);
      /*if(TabStoreTemp[j]>MaxVisu)
      {
      MaxVisu=TabStoreTemp[j];
      TabLargeur[AccuVisu] = ((j+IptFinalMin)-IptInit);
      }*/
   }
   FichierVisu.Close();
   FichierStdIO.Close();
   }






   free(FGaucheRE); 
   free(FDroiteRE);//---> debordement dans les boucles
   free(FCentreGRE); 
   free(FCentreDRE); 


   free(FGaucheIM);
   free(FDroiteIM);
   free(FCentreGIM);
   free(FCentreDIM);

   free(TabDistanceGauche);
   free(TabDistanceDroite);

   free(TabHomogSection);
   free(TabBorSim);


   free(SignalPorte);

}

void CMainFrame::OnEtudecurvi1dEtudeavecobjetEtudeporte12() 
{
   // TODO: Add your command handler code here

   CFile Fichier;
   CIma ObjetCIma;
   int i=0,j=0,k=0,l=0;

   int LargeurPorte = 48;
   int LargeurFourier = 16; // doit être une puissance de deux
   int LongueurSignal = 256;

   int DemiLargeurFourier = LargeurFourier/2;

   int RacineLargeurFourier = (int)sqrt((double)LargeurFourier);

   int DemiLargeurPorte = LargeurPorte/2;
   int DemiLongueurSignal = LongueurSignal/2;

   int * SignalPorte = (int *) malloc(sizeof(int));

   int CondGauche = DemiLongueurSignal-DemiLargeurPorte;
   int CondDroite = DemiLongueurSignal+DemiLargeurPorte;

   double * FGaucheRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FDroiteRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreGRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreDRE = (double *) malloc(sizeof(double)*LargeurFourier);


   double * FGaucheIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FDroiteIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreGIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreDIM = (double *) malloc(sizeof(double)*LargeurFourier);



   for(i=0;i<LargeurFourier;i++)
   {

      FGaucheRE[i]=0.0f;
      FDroiteRE[i]=0.0f;
      FCentreGRE[i]=0.0f;
      FCentreDRE[i]=0.0f;

      FGaucheIM[i]=0.0f;
      FDroiteIM[i]=0.0f;
      FCentreGIM[i]=0.0f;
      FCentreDIM[i]=0.0f;

   }
   int IptInit;
   int IptFinal;
   int LargeurGabaritMin;
   int LargeurGabaritMaxi;
   int IptInitMax = LongueurSignal-LargeurFourier-LargeurPorte;
   int IptFinalMax = LongueurSignal-LargeurFourier-1; //attention 
   int IptInitMin = LargeurFourier;
   int IptFinalMin = LargeurFourier + LargeurPorte;

   int cpt=0;
   int cptTabG=0;
   int cptTabD=0;
   int cptSection=0;
   int cptBord=0;

   CStdioFile FichierStdIO;
   CString ChaineDist;

   int TailleTabDist = (IptInitMax-IptInitMin)*(IptFinalMax-IptFinalMin);

   double * TabDistanceGauche = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabDistanceDroite = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabHomogSection = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabBorSim = (double *) malloc(sizeof(double)*TailleTabDist);

   double * MuLTabDis = (double *) malloc(sizeof(double)*TailleTabDist);

   double * TabMaxVisu = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

   //double * TabStoreTemp = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

   double * TabStoreTemp = (double *)malloc(sizeof(double) *(LongueurSignal));

   int * TabLargeur = (int *)malloc(sizeof(int) *(IptInitMax-IptInitMin));



   for (i=0;i<(LongueurSignal);i++)
   {
      TabStoreTemp[i] = 0.0f;
   }

   for (i=0;i<(IptInitMax-IptInitMin);i++)
   {
      TabMaxVisu[i] = 0.0f;
      TabLargeur[i] = 0;
   }

   for(i=0;i<LargeurFourier;i++)
   {

      TabDistanceGauche[i] = 0.0f;
      TabDistanceDroite[i] = 0.0f;


   }

   for(i=0;i<TailleTabDist;i++)
   {
      MuLTabDis[i] = 0.0f;
      TabDistanceGauche[i]=0.0f;
      TabDistanceDroite[i]=0.0f;
      TabHomogSection[i]=0.0f;
      TabBorSim[i]=0.0f;


   }

   if(Fichier.Open("PorteEtude1DL12.bin",CFile::modeRead))
   {

      Fichier.Read(SignalPorte,sizeof(int)*LongueurSignal);
      Fichier.Close();
   }
   else
   {
      MessageBox("Ouverture du fichier <<PorteEtude1DL12.bin>> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
   }



   //------------------------
   //			RAZ 
   //------------------------

   cptTabG=0;
   cptTabD=0;
   cptSection=0;
   cptBord=0;


   //----------------------------------------------------
   // recopie des tab pour fourier + calcul distance
   //----------------------------------------------------
   for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
   {
      for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) //idem modif 21/11/2006
      {
         //copie Fgauche 

         cpt =0;
         for (i = IptInit-LargeurFourier;i<IptInit;i++)
         {
            FGaucheRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)
         {
            FDroiteRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptInit;i<IptInit+LargeurFourier;i++)
         {
            FCentreGRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptFinal-LargeurFourier+1;i<IptFinal+1;i++)
         {
            FCentreDRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt=0;



         // calcul des FFT
         ObjetCIma.FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);


         TabDistanceGauche[cptTabG] = ObjetCIma.HyperDistanceComplexe(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,DemiLargeurFourier);
         TabDistanceDroite[cptTabD] = ObjetCIma.HyperDistanceComplexe(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,DemiLargeurFourier);
         TabHomogSection[cptSection] = ObjetCIma.HyperDistanceComplexe(FCentreGRE,FCentreGIM,FCentreDRE,FCentreDIM,DemiLargeurFourier);
         TabBorSim[cptBord]=ObjetCIma.HyperDistanceComplexe(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,DemiLargeurFourier);

         cptTabG++;
         cptTabD++;
         cptSection++;
         cptBord++;
      }
   }





   double max = 0.0f;
   int PosMax=0;

   double MaxGauche = 0.0f,MaxDroite = 0.0f,MaxHomoSec = 0.0f, MaxBordSim = 0.0f;



   //-----------------------------------------------------------------------------------------------------
   // recherche des max pour l'inversion des fonctions de probalilité "section homogène" et bord similaire
   //-----------------------------------------------------------------------------------------------------
   for(i=0;i<TailleTabDist;i++)
   {


      if(TabHomogSection[i]>MaxHomoSec)
      {
         MaxHomoSec = TabHomogSection[i];
      }

      if(TabBorSim[i]>MaxBordSim)
      {
         MaxBordSim = TabBorSim[i];
      }

      if(TabDistanceGauche[i]>MaxGauche)
      {
         MaxGauche = TabDistanceGauche[i];
      }

      if(TabDistanceDroite[i]>MaxDroite)
      {
         MaxDroite = TabDistanceDroite[i];
      }

   }




   //--------------------------------------------------------
   //			Fonction export + calcul des couts
   //--------------------------------------------------------

   i=0; // compteur
   j=0;
   int AccuVisu=0;
   double MaxVisu = 0.0f;
   double MaxIter = 0.0f;

   CStdioFile FichierVisu;
   CString ChaineVisu;
   if ((FichierStdIO.Open("ExportBordSimilaireL12.txt",CFile::modeCreate|CFile::modeWrite)) & (FichierVisu.Open("ExportVisuL12.txt",CFile::modeCreate|CFile::modeWrite)) )
   {j=0;
   for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
   {

      MaxIter = 0.0f;
      PosMax = IptInit;
      for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) // modif 21/11/2006 rajout pour limiter les surbalayage
      {

         MuLTabDis[i] = ((MaxHomoSec-TabHomogSection[i])/MaxHomoSec) * (TabDistanceGauche[i]/MaxGauche) * (TabDistanceDroite[i]/MaxDroite)*(MaxBordSim-TabBorSim[i])/MaxBordSim;
         ChaineDist.Format("%f\t",MuLTabDis[i]);
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptInit));
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptFinal));
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptFinal-IptInit));
         FichierStdIO.WriteString(ChaineDist);
         FichierStdIO.WriteString("\n");

         //TabStoreTemp[j] = MuLTabDis[i];

         if(MuLTabDis[i]>MaxIter)
         {
            MaxIter = MuLTabDis[i];
            PosMax = IptFinal;
         }

         i++;
      }//fin IptFinal
      j= (PosMax+IptInit)/2;
      if(TabStoreTemp[j]<MaxIter)TabStoreTemp[j]=MaxIter;// a sauver

      //j++;
      //rechercdu max pour la ligne stockage dans le tableau pour la visu
   }//fin IptInit

   for(j=0;j<(LongueurSignal);j++)
   {
      ChaineVisu.Format("%f\n",TabStoreTemp[j]);
      FichierVisu.WriteString(ChaineVisu);
      /*if(TabStoreTemp[j]>MaxVisu)
      {
      MaxVisu=TabStoreTemp[j];
      TabLargeur[AccuVisu] = ((j+IptFinalMin)-IptInit);
      }*/
   }
   FichierVisu.Close();
   FichierStdIO.Close();
   }






   free(FGaucheRE); 
   free(FDroiteRE);//---> debordement dans les boucles
   free(FCentreGRE); 
   free(FCentreDRE); 


   free(FGaucheIM);
   free(FDroiteIM);
   free(FCentreGIM);
   free(FCentreDIM);

   free(TabDistanceGauche);
   free(TabDistanceDroite);

   free(TabHomogSection);
   free(TabBorSim);


   free(SignalPorte);
}

void CMainFrame::OnEtudecurvi1dEtudeavecobjetEtudeporte12f4() 
{
   // TODO: Add your command handler code here
   // TODO: Add your command handler code here

   CFile Fichier;
   CIma ObjetCIma;
   int i=0,j=0,k=0,l=0;

   int LargeurPorte = 48;
   int LargeurFourier = 4; // doit être une puissance de deux
   int LongueurSignal = 256;

   int DemiLargeurFourier = LargeurFourier/2;

   int RacineLargeurFourier = (int)sqrt((double)LargeurFourier);

   int DemiLargeurPorte = LargeurPorte/2;
   int DemiLongueurSignal = LongueurSignal/2;

   int * SignalPorte = (int *) malloc(sizeof(int));

   int CondGauche = DemiLongueurSignal-DemiLargeurPorte;
   int CondDroite = DemiLongueurSignal+DemiLargeurPorte;

   double * FGaucheRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FDroiteRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreGRE = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreDRE = (double *) malloc(sizeof(double)*LargeurFourier);


   double * FGaucheIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FDroiteIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreGIM = (double *) malloc(sizeof(double)*LargeurFourier);
   double * FCentreDIM = (double *) malloc(sizeof(double)*LargeurFourier);



   for(i=0;i<LargeurFourier;i++)
   {

      FGaucheRE[i]=0.0f;
      FDroiteRE[i]=0.0f;
      FCentreGRE[i]=0.0f;
      FCentreDRE[i]=0.0f;

      FGaucheIM[i]=0.0f;
      FDroiteIM[i]=0.0f;
      FCentreGIM[i]=0.0f;
      FCentreDIM[i]=0.0f;

   }
   int IptInit;
   int IptFinal;
   int LargeurGabaritMin;
   int LargeurGabaritMaxi;
   int IptInitMax = LongueurSignal-LargeurFourier-LargeurPorte;
   int IptFinalMax = LongueurSignal-LargeurFourier-1; //attention 
   int IptInitMin = LargeurFourier;
   int IptFinalMin = LargeurFourier + LargeurPorte;

   int cpt=0;
   int cptTabG=0;
   int cptTabD=0;
   int cptSection=0;
   int cptBord=0;

   CStdioFile FichierStdIO;
   CString ChaineDist;

   int TailleTabDist = (IptInitMax-IptInitMin)*(IptFinalMax-IptFinalMin);

   double * TabDistanceGauche = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabDistanceDroite = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabHomogSection = (double *) malloc(sizeof(double)*TailleTabDist);
   double * TabBorSim = (double *) malloc(sizeof(double)*TailleTabDist);

   double * MuLTabDis = (double *) malloc(sizeof(double)*TailleTabDist);

   double * TabMaxVisu = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

   //double * TabStoreTemp = (double *)malloc(sizeof(double) *(IptInitMax-IptInitMin));

   double * TabStoreTemp = (double *)malloc(sizeof(double) *(LongueurSignal));

   int * TabLargeur = (int *)malloc(sizeof(int) *(IptInitMax-IptInitMin));



   for (i=0;i<(LongueurSignal);i++)
   {
      TabStoreTemp[i] = 0.0f;
   }

   for (i=0;i<(IptInitMax-IptInitMin);i++)
   {
      TabMaxVisu[i] = 0.0f;
      TabLargeur[i] = 0;
   }

   for(i=0;i<LargeurFourier;i++)
   {

      TabDistanceGauche[i] = 0.0f;
      TabDistanceDroite[i] = 0.0f;


   }

   for(i=0;i<TailleTabDist;i++)
   {
      MuLTabDis[i] = 0.0f;
      TabDistanceGauche[i]=0.0f;
      TabDistanceDroite[i]=0.0f;
      TabHomogSection[i]=0.0f;
      TabBorSim[i]=0.0f;


   }

   if(Fichier.Open("PorteEtude1DL12.bin",CFile::modeRead))
   {

      Fichier.Read(SignalPorte,sizeof(int)*LongueurSignal);
      Fichier.Close();
   }
   else
   {
      MessageBox("Ouverture du fichier <<PorteEtude1DL12.bin>> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
   }



   //------------------------
   //			RAZ 
   //------------------------

   cptTabG=0;
   cptTabD=0;
   cptSection=0;
   cptBord=0;


   //----------------------------------------------------
   // recopie des tab pour fourier + calcul distance
   //----------------------------------------------------
   for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
   {
      for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) //idem modif 21/11/2006
      {
         //copie Fgauche 

         cpt =0;
         for (i = IptInit-LargeurFourier;i<IptInit;i++)
         {
            FGaucheRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)
         {
            FDroiteRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptInit;i<IptInit+LargeurFourier;i++)
         {
            FCentreGRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt =0;
         for (i = IptFinal-LargeurFourier+1;i<IptFinal+1;i++)
         {
            FCentreDRE[cpt]= (double)SignalPorte[i];
            cpt++;
         }

         cpt=0;



         // calcul des FFT
         ObjetCIma.FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);
         ObjetCIma.FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);


         TabDistanceGauche[cptTabG] = ObjetCIma.HyperDistanceComplexe(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,DemiLargeurFourier);
         TabDistanceDroite[cptTabD] = ObjetCIma.HyperDistanceComplexe(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,DemiLargeurFourier);
         TabHomogSection[cptSection] = ObjetCIma.HyperDistanceComplexe(FCentreGRE,FCentreGIM,FCentreDRE,FCentreDIM,DemiLargeurFourier);
         TabBorSim[cptBord]=ObjetCIma.HyperDistanceComplexe(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,DemiLargeurFourier);

         cptTabG++;
         cptTabD++;
         cptSection++;
         cptBord++;
      }
   }





   double max = 0.0f;
   int PosMax=0;

   double MaxGauche = 0.0f,MaxDroite = 0.0f,MaxHomoSec = 0.0f, MaxBordSim = 0.0f;



   //-----------------------------------------------------------------------------------------------------
   // recherche des max pour l'inversion des fonctions de probalilité "section homogène" et bord similaire
   //-----------------------------------------------------------------------------------------------------
   for(i=0;i<TailleTabDist;i++)
   {


      if(TabHomogSection[i]>MaxHomoSec)
      {
         MaxHomoSec = TabHomogSection[i];
      }

      if(TabBorSim[i]>MaxBordSim)
      {
         MaxBordSim = TabBorSim[i];
      }

      if(TabDistanceGauche[i]>MaxGauche)
      {
         MaxGauche = TabDistanceGauche[i];
      }

      if(TabDistanceDroite[i]>MaxDroite)
      {
         MaxDroite = TabDistanceDroite[i];
      }

   }




   //--------------------------------------------------------
   //			Fonction export + calcul des couts
   //--------------------------------------------------------

   i=0; // compteur
   j=0;
   int AccuVisu=0;
   double MaxVisu = 0.0f;
   double MaxIter = 0.0f;

   CStdioFile FichierVisu;
   CString ChaineVisu;
   if ((FichierStdIO.Open("ExportBordSimilaireL12.txt",CFile::modeCreate|CFile::modeWrite)) & (FichierVisu.Open("ExportVisuL12.txt",CFile::modeCreate|CFile::modeWrite)) )
   {j=0;
   for (IptInit=IptInitMin;IptInit<IptInitMax;IptInit++)
   {

      MaxIter = 0.0f;
      PosMax = IptInit;
      for(IptFinal=(IptInit+(LargeurFourier)+1);IptFinal<IptFinalMax;IptFinal++) // modif 21/11/2006 rajout pour limiter les surbalayage
      {

         MuLTabDis[i] = ((MaxHomoSec-TabHomogSection[i])/MaxHomoSec) * (TabDistanceGauche[i]/MaxGauche) * (TabDistanceDroite[i]/MaxDroite)*(MaxBordSim-TabBorSim[i])/MaxBordSim;
         ChaineDist.Format("%f\t",MuLTabDis[i]);
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptInit));
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptFinal));
         FichierStdIO.WriteString(ChaineDist);

         ChaineDist.Format("%d\t",(IptFinal-IptInit));
         FichierStdIO.WriteString(ChaineDist);
         FichierStdIO.WriteString("\n");

         //TabStoreTemp[j] = MuLTabDis[i];

         if(MuLTabDis[i]>MaxIter)
         {
            MaxIter = MuLTabDis[i];
            PosMax = IptFinal;
         }

         i++;
      }//fin IptFinal
      j= (PosMax+IptInit)/2;
      if(TabStoreTemp[j]<MaxIter)TabStoreTemp[j]=MaxIter;// a sauver

      //j++;
      //rechercdu max pour la ligne stockage dans le tableau pour la visu
   }//fin IptInit

   for(j=0;j<(LongueurSignal);j++)
   {
      ChaineVisu.Format("%f\n",TabStoreTemp[j]);
      FichierVisu.WriteString(ChaineVisu);
      /*if(TabStoreTemp[j]>MaxVisu)
      {
      MaxVisu=TabStoreTemp[j];
      TabLargeur[AccuVisu] = ((j+IptFinalMin)-IptInit);
      }*/
   }
   FichierVisu.Close();
   FichierStdIO.Close();
   }






   free(FGaucheRE); 
   free(FDroiteRE);//---> debordement dans les boucles
   free(FCentreGRE); 
   free(FCentreDRE); 


   free(FGaucheIM);
   free(FDroiteIM);
   free(FCentreGIM);
   free(FCentreDIM);

   free(TabDistanceGauche);
   free(TabDistanceDroite);

   free(TabHomogSection);
   free(TabBorSim);


   free(SignalPorte);
}



void CMainFrame::OnEtudecurvi1dGenerationfichierprofilBitexture48() 
{
   // TODO: Add your command handler code here
   CFile Fichier;
   int i=0;

   int LargeurPorte = 48;
   //int LargeurFourier = 16;
   int LongueurSignal = 256;

   int DemiLargeurPorte = LargeurPorte/2;
   int DemiLongueurSignal = LongueurSignal/2;

   int * SignalPorte = (int *) malloc(sizeof(int)*LongueurSignal);

   int CondGauche = DemiLongueurSignal-DemiLargeurPorte;
   int CondDroite = DemiLongueurSignal+DemiLargeurPorte;
   double AmplitudeSinusExt = 10.0f;
   double AmplitudeSinusInt = 10.0f;
   double FrequenceInt = 10.0f;
   double FrequenceExt = 50.0f;
   double ValPi = 3.14159265358979;


   CStdioFile ExportTxt;
   CString ChaineTemp;
   for (i=0;i<LongueurSignal;i++)
   {
      if ((i>CondGauche)&(i<CondDroite))
      {
         SignalPorte[i]=(int)(AmplitudeSinusInt+AmplitudeSinusInt*sin((2*ValPi*FrequenceInt*i)/(CondDroite-CondGauche)));
      }
      else
      {
         SignalPorte[i]=(int)(AmplitudeSinusExt+AmplitudeSinusExt*sin((2*ValPi*FrequenceExt*i)/(CondGauche)));
      }
   }
   if(Fichier.Open("PorteEtude1DL48BiTexture.bin",CFile::modeWrite|CFile::modeCreate)&&ExportTxt.Open("ExportXLS.txt",CFile::modeWrite|CFile::modeCreate))
   {

      for (i=0;i<LongueurSignal;i++)
      {
         ChaineTemp.Format("%d\n",SignalPorte[i]);
         ExportTxt.WriteString(ChaineTemp);
      }
      Fichier.Write(SignalPorte,sizeof(int)*LongueurSignal);
      Fichier.Close();
      ExportTxt.Close();
   }
   else
   {
      MessageBox("Ouverture du fichier <<PorteEtude1DL48BiTexture.bin>> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
   }

   free(SignalPorte);
}

void CMainFrame::OnEtudecurvi1dEtudeavecobjetEtudebitexture48() 
{
#ifdef MP
   CCurviDetect1D Detect1D;

   Detect1D.Quality("PorteEtude1DL48BiTexture.bin","VisuBiTextureL48F16.txt","VisuBiTextureL48F16All.txt",16,256);
#endif
}

void CMainFrame::OnEtudecurvi1dEtudeavecobjetEtudebitexture48f4() 
{
#ifdef MP
   CCurviDetect1D Detect1D;


   Detect1D.Quality("PorteEtude1DL48BiTexture.bin","VisuBiTextureL48F4.txt","VisuBiTextureL48F4All.txt",4,256);
#endif
}

void CMainFrame::OnEtudecurvi1dGenerationfichierprofilPortebordnonsim48() 
{
   // TODO: Add your command handler code here
   CFile Fichier;
   int i=0;

   int LargeurPorte = 48;
   //int LargeurFourier = 16;
   int LongueurSignal = 256;

   int DemiLargeurPorte = LargeurPorte/2;
   int DemiLongueurSignal = LongueurSignal/2;

   int * SignalPorte = (int *) malloc(sizeof(int)*LongueurSignal);

   int CondGauche = DemiLongueurSignal-DemiLargeurPorte;
   int CondDroite = DemiLongueurSignal+DemiLargeurPorte;
   /*double AmplitudeSinusExt = 10.0f;
   double AmplitudeSinusInt = 10.0f;
   double FrequenceInt = 10.0f;
   double FrequenceExt = 50.0f;
   double ValPi = 3.14159265358979;*/


   CStdioFile ExportTxt;
   CString ChaineTemp;
   for (i=0;i<LongueurSignal;i++)
   {
      if ((i>CondGauche)&(i<CondDroite))
      {
         SignalPorte[i]=255;
      }
      if (i<=CondGauche)
      {
         SignalPorte[i]=0;
      }
      if (i>=CondDroite)
      {
         SignalPorte[i]=100;
      }

   }
   if(Fichier.Open("Porte1DL48NonSim.bin",CFile::modeWrite|CFile::modeCreate)&&ExportTxt.Open("ExportXLSPorte1DL48NonSim.txt",CFile::modeWrite|CFile::modeCreate))
   {

      for (i=0;i<LongueurSignal;i++)
      {
         ChaineTemp.Format("%d\n",SignalPorte[i]);
         ExportTxt.WriteString(ChaineTemp);
      }
      Fichier.Write(SignalPorte,sizeof(int)*LongueurSignal);
      Fichier.Close();
      ExportTxt.Close();
   }
   else
   {
      MessageBox("Ouverture du fichier <<PorteEtude1DL48BiTexture.bin>> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
   }

   free(SignalPorte);

}

void CMainFrame::OnEtudecurvi1dEtudeavecobjetEtudeportenonsim48f16() 
{
#ifdef MP
   CCurviDetect1D Detect1D;


   Detect1D.Quality("Porte1DL48NonSim.bin","VisuPorteNonSimL48F16.txt","VisuPorteNomSimL48F16All.txt",16,256);
#endif
}

void CMainFrame::OnEtudecurvi1dGenerationfichierprofilMultisignaux1d() 
{
   // TODO: Add your command handler code here
   CFile Fichier;
   int i=0;

   int LargeurPorte = 48;
   //int LargeurFourier = 16;
   int LongueurSignal = 256;

   int DemiLargeurPorte = LargeurPorte/2;
   int DemiLongueurSignal = LongueurSignal/2;

   int * SignalPorte = (int *) malloc(sizeof(int)*LongueurSignal);

   int CondGauche1 = (LongueurSignal/3)-DemiLargeurPorte;
   int CondDroite1 = (LongueurSignal/3)+DemiLargeurPorte;

   int CondGauche2 = (2*LongueurSignal/3)-DemiLargeurPorte;
   int CondDroite2 = (2*LongueurSignal/3)+DemiLargeurPorte;

   /*double AmplitudeSinusExt = 10.0f;
   double AmplitudeSinusInt = 10.0f;
   double FrequenceInt = 10.0f;
   double FrequenceExt = 50.0f;
   double ValPi = 3.14159265358979;*/


   CStdioFile ExportTxt;
   CString ChaineTemp;
   for (i=0;i<LongueurSignal;i++)
   {
      if ((i<=CondGauche1))
      {
         SignalPorte[i]=0;
      }

      if ((i>CondGauche1)&(i<LongueurSignal/3))
      {
         SignalPorte[i]=255;
      }
      if ((i>=LongueurSignal/3)&(i<CondDroite1))
      {
         SignalPorte[i]=127;
      }

      if ((i>=CondDroite1)&(i<=CondGauche2))
      {
         SignalPorte[i]=0;
      }
      if ((i>CondGauche2)&(i<CondDroite2))
      {
         SignalPorte[i]=255;
      }
      if (i>=CondDroite2)
      {
         SignalPorte[i]=0;
      }

   }
   if(Fichier.Open("SignalMultipleL48.bin",CFile::modeWrite|CFile::modeCreate)&&ExportTxt.Open("ExportXLSSignalMultipleL48.txt",CFile::modeWrite|CFile::modeCreate))
   {

      for (i=0;i<LongueurSignal;i++)
      {
         ChaineTemp.Format("%d\n",SignalPorte[i]);
         ExportTxt.WriteString(ChaineTemp);
      }
      Fichier.Write(SignalPorte,sizeof(int)*LongueurSignal);
      Fichier.Close();
      ExportTxt.Close();
   }
   else
   {
      MessageBox("Ouverture du fichier <<PorteEtude1DL48BiTexture.bin>> impossible","Information",MB_OK|MB_ICONEXCLAMATION);
   }

   free(SignalPorte);


}

void CMainFrame::OnEtudecurvi1dEtudeavecobjetEtudesignauxmultiple48f16() 
{
#ifdef MP
   CCurviDetect1D Detect1D;


   Detect1D.Quality("SignalMultipleL48.bin","VisuSMultiL48F16.txt","VisuSMultiL48F16All.txt",16,256);
#endif
}

void CMainFrame::OnEtudecurvi1dEtudeavecobjetEtudesignauxrels() 
{
   CFile FichierBin;
   CStdioFile FichierXlc,ExportXls;
   int i,pos,Taille;
   int LongueurSignal = 179; ///-------> modif en fonction de la coupe
   int * SignalXlc = (int*)malloc(sizeof(int)*LongueurSignal);
   CString ChaineSplit;
   CString ChaineXlc;
   //converions xlc -> bin et txt pour visu


   if ((FichierXlc.Open("coupe.xlc",CFile::modeRead))&&(FichierBin.Open("CoupeReelle.bin",CFile::modeCreate|CFile::modeWrite))&&(ExportXls.Open("CoupeReelleVisu.txt",CFile::modeCreate|CFile::modeWrite)))
   {
      for(i=0;i<LongueurSignal;i++)
      {
         FichierXlc.ReadString(ChaineXlc);

         pos = ChaineXlc.Find(" ");

         Taille=ChaineXlc.GetLength()-pos-1;

         ChaineSplit = ChaineXlc.Right(Taille);
         SignalXlc[i]=atoi((LPCTSTR)ChaineSplit);
      }


      FichierBin.Write(SignalXlc,LongueurSignal*sizeof(int));


      for(i=0;i<LongueurSignal;i++)
      {
         ChaineXlc.Format("%d\n",SignalXlc[i]);
         ExportXls.WriteString(ChaineXlc);
      }

      ExportXls.Close();
      FichierBin.Close();
      FichierXlc.Close();	
   }
   free(SignalXlc);

   //--------------------------------------
   // Etude Curvi
#ifdef MP
   CCurviDetect1D Detect1D;


   Detect1D.Quality("CoupeReelle.bin","VisuCoupeReelleF16.txt","VisuCoupeReelleF16All.txt",16,LongueurSignal);
#endif
}

void CMainFrame::OnEtudecurvi1dEtudeavecobjetTestnew() 
{
   // TODO: Add your command handler code here
   int LongueurSignal = 256;
#ifdef MP
   CCurviDetect1D Detect1D;


   Detect1D.Quality("PorteEtude1DL48.bin","testnew.txt","testnewF16All.txt",16,LongueurSignal);
#endif
}

void CMainFrame::OnEtudecurvi1dEtudeavecobjetEtudemultifichier() 
{
#ifdef MP
   // TODO: Add your command handler code here

   //-------------------------------------------------------------------------------------
   //					relecture du fichier txt contenant toutes les coupes
   //-------------------------------------------------------------------------------------

   CStdioFile FichierChemin; 
   CStringArray TabChemin;
   int ValeurRtn;
   CString CheminFichier;
   CString ChaineLue;
   BOOL BLu = TRUE;
   int cpt=0;
   int NBLu;
   int i=0,j=0;
   int * TabLongueurCoupe;




   CFile FichierBin;
   CStdioFile FichierXlc,ExportXls;
   int pos,Taille;
   int LongueurSignal ;//= 179; ///-------> modif en fonction de la coupe
   int * SignalXlc ;//= (int*)malloc(sizeof(int)*LongueurSignal);
   CString ChaineSplit;
   CString ChaineXlc;
   CString ChaineExtBin = ".bin";
   CString ChaineCoupe = "coupe";
   CString NomFichierCoupe = "";
   CString ChaineTemp;
   CString NomFichierPourVisu="CoupeReelleVisu_";
   CString ExtensionFichierVisu = ".txt";
   CString NomSaveFichVisu = "";

   CString NomFichierPourVisuAll="CoupeReelleVisu_";
   CString ExtensionFichierVisuAll = "All.txt";
   CString NomSaveFichVisuAll = "";

   CString NomFichierPourVisuAff="CoupeReelleVisua_";
   CString ExtensionFichierVisuAff = "Aff.txt";
   CString NomSaveFichVisuAff = "";

   CString ChaineLectureLongueur = "";

   BOOL ConditionOuverture = FALSE;
   BOOL ConditionLectureLongueur = TRUE;
   CCurviDetect1D Detect1D;


   CFileDialog BoiteSelectFichier(TRUE,NULL,NULL,OFN_OVERWRITEPROMPT,"Fichiers coupes (*.txt)|*.txt||");
   ValeurRtn = BoiteSelectFichier.DoModal();

   CheminFichier = BoiteSelectFichier.GetPathName();


   if(ValeurRtn == IDOK)
   {
      ///MessageBox(CheminFichier);
      if(FichierChemin.Open(CheminFichier,CFile::modeRead))
      {
         //TailleFichier = FichierChemin.GetLength();
         while(BLu == TRUE)
         {
            BLu = FichierChemin.ReadString(ChaineLue);
            TabChemin.Add(ChaineLue);
         }
         FichierChemin.Close();
      }
   }
   else
   {
      MessageBox("Pas de fichiers sélectionnés!");
      return;
   }

   NBLu = TabChemin.GetSize();

   //-------------------------------------------------------------------------------------
   //					Pour toutes les coupes 
   //					1) conversion en fichier *.bin
   //					2) etude curvi sauvegarde des résultats
   //-------------------------------------------------------------------------------------
   TabLongueurCoupe = (int*)malloc(sizeof(int)*NBLu);



   for(i=0;i<NBLu-1;i++)//pour toutes les coupes
   {
      ChaineTemp.Format("%d",i);
      NomFichierCoupe = ChaineCoupe + ChaineTemp + ChaineExtBin;
      NomSaveFichVisu = NomFichierPourVisu + ChaineTemp + ExtensionFichierVisu;
      NomSaveFichVisuAll = NomFichierPourVisuAll + ChaineTemp + ExtensionFichierVisuAll;
      NomSaveFichVisuAff = NomFichierPourVisuAff + ChaineTemp + ExtensionFichierVisuAff;

      ConditionOuverture = (FichierXlc.Open(TabChemin[i],CFile::modeRead))
         &&(FichierBin.Open(NomFichierCoupe,CFile::modeCreate|CFile::modeWrite)
         )&&(ExportXls.Open(NomSaveFichVisu,CFile::modeCreate|CFile::modeWrite));

      if (ConditionOuverture)
      {
         //CString ChaineLectureLongueur = "";

         LongueurSignal=0;
         ConditionLectureLongueur=TRUE;

         while(ConditionLectureLongueur==TRUE)
         {
            ConditionLectureLongueur = FichierXlc.ReadString(ChaineLectureLongueur);
            LongueurSignal++;
         }

         LongueurSignal = LongueurSignal-1;

         FichierXlc.SeekToBegin();

         SignalXlc = (int*)malloc(sizeof(int)*LongueurSignal);

         for(j=0;j<LongueurSignal;j++)
         {
            FichierXlc.ReadString(ChaineXlc);

            pos = ChaineXlc.Find(" ");

            Taille=ChaineXlc.GetLength()-pos-1;

            ChaineSplit = ChaineXlc.Right(Taille);
            SignalXlc[j]=atoi((LPCTSTR)ChaineSplit);
         }


         FichierBin.Write(SignalXlc,LongueurSignal*sizeof(int));


         for(j=0;j<LongueurSignal;j++)
         {
            ChaineXlc.Format("%d\n",SignalXlc[j]);
            ExportXls.WriteString(ChaineXlc);
         }

         ExportXls.Close();
         FichierBin.Close();
         FichierXlc.Close();	
      }

      //calcul proba ici
      Detect1D.Quality(NomFichierCoupe,NomSaveFichVisuAff,NomSaveFichVisuAll,16,LongueurSignal);

      free(SignalXlc);

      NomFichierCoupe ="";
      NomSaveFichVisu="";
      NomSaveFichVisuAll="";
   }

#endif
}

void CMainFrame::OnStegerConvolution1d()
{
   //-------------------------------------------------------------------------------------
   //					relecture du fichier txt contenant toutes les coupes
   //-------------------------------------------------------------------------------------

   CStdioFile FichierChemin; 
   CStringArray TabChemin;
   int ValeurRtn;
   CString CheminFichier;
   CString ChaineLue;
   BOOL BLu = TRUE;
   int cpt=0;
   int NBLu;
   int i=0,j=0;
   int * TabLongueurCoupe;




   CFile FichierBin;
   CStdioFile FichierXlc,ExportXls;
   int pos,Taille,IndTabRslts=0;
   int LongueurSignal=0 ;
   int * SignalXlc ;
   double * SignalDouble;
   double Theta;
   long TailleMasque;
   CString ChaineSplit;
   CString ChaineXlc;
   CString ChaineExtBin = ".bin";
   CString ChaineCoupe = "coupe";
   CString NomFichierCoupe = "";
   CString ChaineTemp;
   CString NomFichierPourVisu="CoupeReelleVisu_";
   CString ExtensionFichierVisu = ".txt";
   CString NomSaveFichVisu = "";

   CString NomFichierPourVisuAll="CoupeReelleVisu_";
   CString ExtensionFichierVisuAll = "_AllResults.txt";
   CString NomSaveFichVisuAll = "";

   CString NomFichierPourVisuAff="CoupeReelleVisua_";
   CString ExtensionFichierVisuAff = "Aff.txt";
   CString NomSaveFichVisuAff = "";

   CString ChaineLectureLongueur = "";

   BOOL ConditionOuverture = FALSE;
   BOOL ConditionLectureLongueur = TRUE;
   CSteger1D Steger1D;

   CString ChaineEportFichier;
   CString ChaineEportFichierTmp;
   CStdioFile FichierExportAllResults;



   CFileDialog BoiteSelectFichier(TRUE,NULL,NULL,OFN_OVERWRITEPROMPT,"Fichiers coupes (*.txt)|*.txt||");
   ValeurRtn = BoiteSelectFichier.DoModal();

   CheminFichier = BoiteSelectFichier.GetPathName();


   if(ValeurRtn == IDOK)
   {
      ///MessageBox(CheminFichier);
      if(FichierChemin.Open(CheminFichier,CFile::modeRead))
      {
         //TailleFichier = FichierChemin.GetLength();
         while(BLu == TRUE)
         {
            BLu = FichierChemin.ReadString(ChaineLue);
            if (BLu == TRUE)
            {
               TabChemin.Add(ChaineLue);
            }
         }
         FichierChemin.Close();
      }
   }
   else
   {
      MessageBox("Pas de fichiers sélectionnés!");
      return;
   }

   NBLu = TabChemin.GetSize();

   //-------------------------------------------------------------------------------------
   //					Pour toutes les coupes 
   //					1) conversion en fichier *.bin
   //					2) etude curvi sauvegarde des résultats
   //-------------------------------------------------------------------------------------
   TabLongueurCoupe = (int*)malloc(sizeof(int)*NBLu);



   for(i=0;i<NBLu;i++)//pour toutes les coupes
   {
      ChaineTemp.Format("%d",i);
      NomFichierCoupe = ChaineCoupe + ChaineTemp + ChaineExtBin;
      NomSaveFichVisu = NomFichierPourVisu + ChaineTemp + ExtensionFichierVisu;
      NomSaveFichVisuAll = NomFichierPourVisuAll + ChaineTemp + ExtensionFichierVisuAll;
      NomSaveFichVisuAff = NomFichierPourVisuAff + ChaineTemp + ExtensionFichierVisuAff;

      ConditionOuverture = (FichierXlc.Open(TabChemin[i],CFile::modeRead))
         &&(FichierBin.Open(NomFichierCoupe,CFile::modeCreate|CFile::modeWrite)
         )&&(ExportXls.Open(NomSaveFichVisu,CFile::modeCreate|CFile::modeWrite));

      if (ConditionOuverture)
      {
         //CString ChaineLectureLongueur = "";

         LongueurSignal=0;
         ConditionLectureLongueur=TRUE;

         while(ConditionLectureLongueur==TRUE)
         {
            ConditionLectureLongueur = FichierXlc.ReadString(ChaineLectureLongueur);
            LongueurSignal++;
         }

         LongueurSignal = LongueurSignal-1;

         FichierXlc.SeekToBegin();

         SignalXlc = (int*)malloc(sizeof(int)*LongueurSignal);
         SignalDouble = (double*)malloc(sizeof(double)*LongueurSignal);

         for(j=0;j<LongueurSignal;j++)
         {
            FichierXlc.ReadString(ChaineXlc);

            pos = ChaineXlc.Find(" ");

            Taille=ChaineXlc.GetLength()-pos-1;

            ChaineSplit = ChaineXlc.Right(Taille);
            SignalXlc[j]=atoi((LPCTSTR)ChaineSplit);
            SignalDouble[j]= (double)SignalXlc[j];
         }


         FichierBin.Write(SignalXlc,LongueurSignal*sizeof(int));

         //Plus nécessaire
         for(j=0;j<LongueurSignal;j++)
         {
            ChaineXlc.Format("%d\n",SignalXlc[j]);
            ExportXls.WriteString(ChaineXlc);
         }

         ExportXls.Close();
         FichierBin.Close();
         FichierXlc.Close();	
      }


      //relecture des coupes où mettre les valeurs dans un Array de double
      TailleMasque=21;
      Theta=1;
      Steger1D.ConvolutionKernels(SignalDouble,LongueurSignal,TailleMasque,Theta);

      if (FichierExportAllResults.Open(NomSaveFichVisuAll,CFile::modeCreate|CFile::modeWrite))
      {
         int CptKernel=0;
         for(IndTabRslts=0;IndTabRslts<LongueurSignal;IndTabRslts++)
         {

            ChaineEportFichierTmp="";
            ChaineEportFichier="";

            ChaineEportFichierTmp.Format("%f",SignalDouble[IndTabRslts]);
            ChaineEportFichier +=ChaineEportFichierTmp;
            ChaineEportFichier += "\t";

            ChaineEportFichierTmp.Format("%f",Steger1D.TabResult[IndTabRslts]);
            ChaineEportFichier +=ChaineEportFichierTmp;
            ChaineEportFichier += "\t";

            ChaineEportFichierTmp.Format("%f",Steger1D.TabResultD1[IndTabRslts]);
            ChaineEportFichier +=ChaineEportFichierTmp;
            ChaineEportFichier += "\t";

            ChaineEportFichierTmp.Format("%f",Steger1D.TabResultD2[IndTabRslts]);
            ChaineEportFichier +=ChaineEportFichierTmp;
            ChaineEportFichier += "\t";
            if(CptKernel<TailleMasque)
            {
               ChaineEportFichierTmp.Format("%f",Steger1D.GaussianKernel[CptKernel]);
               ChaineEportFichier +=ChaineEportFichierTmp;
               ChaineEportFichier += "\t";
               ChaineEportFichierTmp.Format("%f",Steger1D.GaussianKernelD1[CptKernel]);
               ChaineEportFichier +=ChaineEportFichierTmp;
               ChaineEportFichier += "\t";
               ChaineEportFichierTmp.Format("%f",Steger1D.GaussianKernelD2[CptKernel]);
               ChaineEportFichier +=ChaineEportFichierTmp;
               ChaineEportFichier += "\n";
               CptKernel++;
            }
            FichierExportAllResults.WriteString(ChaineEportFichier);
            if(CptKernel==TailleMasque)
            {
               FichierExportAllResults.WriteString("\n");
            }


         }

         FichierExportAllResults.Close();
      }
      Steger1D.EffaceTab();
      free(SignalXlc);
      free(SignalDouble);

      NomFichierCoupe ="";
      NomSaveFichVisu="";
      NomSaveFichVisuAll="";
   }
}

void CMainFrame::OnConversionSignal()
{
   CFileDialog BoiteSelectFichier(TRUE,NULL,NULL,OFN_OVERWRITEPROMPT,"Signal (*.bin)|*.bin||");



   CFile FichierChemin;
   int ValeurLue,j;
   int * SignalBin;
   long LongueurSignal;
   int ConditionLectureLongueur = -1;
   CString ChaineCoupe = "";
   CString ChaineTemp = "";
   CStdioFile FichierCoupe;
   CString NomFichierCoupe ="",CheminFichier="";
   BOOL ValeurRtn;

   ValeurRtn = BoiteSelectFichier.DoModal();

   if(ValeurRtn == IDOK)
   {
      CheminFichier = BoiteSelectFichier.GetPathName();
      NomFichierCoupe = CheminFichier.Left(CheminFichier.GetLength()-4);
      NomFichierCoupe += ".xlc";
      if(FichierChemin.Open(CheminFichier,CFile::modeRead)&&FichierCoupe.Open(NomFichierCoupe,CFile::modeCreate|CFile::modeWrite))
      {

         LongueurSignal = FichierChemin.GetLength()/sizeof(int);

         //FichierChemin.SeekToBegin();
         SignalBin = (int*)malloc(sizeof(int)*LongueurSignal);
         for(j=0;j<LongueurSignal;j++)
            SignalBin[j]=0;

         FichierChemin.Read(SignalBin,LongueurSignal*sizeof(int));
         for(j=0;j<LongueurSignal;j++)
         {
            ChaineCoupe.Format("%d",j+1);
            ChaineCoupe += " ";
            ChaineTemp.Format("%d",SignalBin[j]);
            ChaineCoupe +=ChaineTemp;
            ChaineCoupe +="\n";
            FichierCoupe.WriteString(ChaineCoupe);
         }

         FichierCoupe.Close();
         FichierChemin.Close();
         free(SignalBin);
      }

   }
   else
   {
      MessageBox("Pas de fichiers sélectionnés!");
      return;
   }	

}

void CMainFrame::OnOutilsAjoutbruit()
{
   CFileDialog BoiteSelectFichier(TRUE,NULL,NULL,OFN_OVERWRITEPROMPT,"Signal (*.bin)|*.bin||");



   CFile FichierChemin;
   int ValeurLue,j;
   int * SignalBin;
   long LongueurSignal;
   int ConditionLectureLongueur = -1;
   CString ChaineCoupe = "";
   CString ChaineTemp = "";
   CFile FichierCoupe;
   CString NomFichierCoupe ="",CheminFichier="";
   BOOL ValeurRtn;
   double RandomInteger;
   double PourcentageSignal = 0.1 * 255;

   ValeurRtn = BoiteSelectFichier.DoModal();
   srand((unsigned int)time(0));

   if(ValeurRtn == IDOK)
   {
      CheminFichier = BoiteSelectFichier.GetPathName();
      NomFichierCoupe = CheminFichier.Left(CheminFichier.GetLength()-4);
      NomFichierCoupe += "_bruit.bin";
      if(FichierChemin.Open(CheminFichier,CFile::modeRead)&&FichierCoupe.Open(NomFichierCoupe,CFile::modeCreate|CFile::modeWrite))
      {

         LongueurSignal = FichierChemin.GetLength()/sizeof(int);

         //FichierChemin.SeekToBegin();
         SignalBin = (int*)malloc(sizeof(int)*LongueurSignal);


         FichierChemin.Read(SignalBin,LongueurSignal*sizeof(int));

         for(j=0;j<LongueurSignal;j++)
         {
            RandomInteger = (double)((double)rand()/(double)RAND_MAX);
            RandomInteger = RandomInteger*PourcentageSignal;
            SignalBin[j] = SignalBin[j]	+ (int)RandomInteger;

            /*if (SignalBin[j]>255) 
            {
            SignalBin[j] = SignalBin[j]	- (int)(2*RandomInteger);
            }
            if (SignalBin[j]<0) SignalBin[j]=0;	*/
         }

         FichierCoupe.Write(SignalBin,LongueurSignal*sizeof(int));
         FichierCoupe.Close();
         FichierChemin.Close();
         free(SignalBin);
      }

   }
   else
   {
      MessageBox("Pas de fichiers sélectionnés!");
      return;
   }
}
