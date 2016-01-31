#ifndef __IMAGE_H__     
#define __IMAGE_H__     
                        
                        
#include "stdafx.h"

 
#include "resource.h"

#include "pal.h"                   
#include "constdef.h"   
#include "message.h" 
#include "uselib.h"
#include "DialInfoTrait.h"  
#include "afx.h"
#include "BoiteParamCurviDetect.h"
#include "DialogAffine.h"


#define _CRT_SECURE_NO_DEPRECATE

#include "struct.h"
#include "curvlinear.h"


struct Pixel { unsigned short i;unsigned short j;};
struct point { double X;double Y;double Z;}; 
struct pt    {double X;double Y;};

class CDialLecSeq;
class CProfilWnd;


 /////////////////////////////////////////////////////////////////////////////
class CImageWnd : public CMDIChildWnd
{   
	DECLARE_DYNAMIC(CImageWnd)
		
public:
	CStatusBar  m_wndStatusBarImage;
	//****************************************************
	//************* DONNEES MEMBRES PUBLIQUES ************
	//****************************************************

    // -------------------- le type de l'image ---------------------
    public:
	NBCOLOR TypeColor;   	// TypeColor = NB si N&B
							// TypeColor = COLOR si image couleur
    public:
	TYPEIMAGE TypeImage; // TypeImage = OCTET8 si *.BMP
						// TypeImage = FLOAT32 si Image en float 32 bits
	
	// ************************ BUFFERS image typé *********************
    // 
	public:
	unsigned char  * Image;
	rgb  * ImageRGB;
	// ******************************************************************
    public:
	float  * ImageFloat;
    float mini, maxi;
    //------------------- dimensions de l'image ------------------------
	public:                                          
	
	long Larg;  // largeur
	long Haut;  // hauteur
	DWORD Taille; // larg*haut ou larg*haut*3
    // -----------------------------------------------------------------
	
	// ---- debut et fin de zone sélectionnée au clavier ou souris -----	    
	public:
	long debut_i,fin_i;
	long debut_j,fin_j;   

	// ----------- position de la coupe avec souris -------------------		
	public:
	long Lg_i1,Lg_i2,Lg_j1,Lg_j2; // (i,j)
	long NbPointCoupe; // Nb de points du profil 
	
	
	//*************** variables de séquence d'images ********************

	public:
	BOOL SEQ_EXISTE; //  ouverture d'une séquence d'images du disque :TRUE si séquence
	
	public:
	long TailleSeq;  	// Taille de la séquence 
	int NbImage;     	// Nombre d'images demandé
	int NbMaxImageSeq; 	// Nombre d'images accepté  
	BOOL SEQUENCE_CREE; //    
	unsigned char  * * Bpile;   // Bpile est un tableau d'images
	CDialLecSeq *DialLecSeq;
	long debut_i_seq,fin_i_seq;
	long debut_j_seq,fin_j_seq;
	// ------------------------------------------------------------------
  
	
	// ---------------- position courante de la souris en I,J dans l'image ----------
	public:	
	int I,J;
	// ------------------------------------------------------------------------------
	
	// ----------------- seuil pour binarisation ------------------------------------	
	public:
	unsigned char Seuil,Seuil2;
    // ------------------------------------------------------------------------------
	       
	
	// -------------------- Variable utilisée pour la transformée en Ondelettes       
	public:                                                                    
	int Echelle;
	
	// ----------------- nom du fichier BMP	------------------------------------------       
	public:
	CString FileName ; 	// nom + chemin       
	// ------------------------------------------------------------------------------       
	
	// ----------------- Titre de la fenêtre ----------------------------------------       
	public:
	CString Titre;     
	// -----------------------------------------------------------------------------
	
	// -------------------------- Objet menu chargé	--------------------------------
	public:
	CMenu menu;           
	// -----------------------------------------------------------------------------
    
	public:
	int ZOOM;  // zoom de visualisation
	int m_FacteurZoom;	// Facteur de grossissement
	
//****************************************************
//************* DONNEES MEMBRES PRIVEES **************
//****************************************************
	// BUFFER image	non typé	    
	private:
	void * PImage;       	// Pointeur sur l'image affichée
	void * PImageTemp;	 	// Pointeur sur l'image initiale pour une image fichier    
	SImage Im;				// structure image
	SImageFloat ImFloat;	// structure image float
	SImageRGB ImRGB;
			
	// variables de gestion souris
	private:
	BOOL Appuie; // bouton toujours appuié
	BOOL first;  // premier fois dans MouseMove après l'appui
	BOOL TRACER_COUPE;
	BOOL AfficheCoupe; 
	        
	// pointeur tempon	        
	private:
	unsigned char   * Buf;       
	

	
	private:	 	
	// le rectangle client de la fenetre car les dimensions de la fenetre bouge
	CRect TheRectClient;  
	
	private:
	// première et dernière position de la souris entre premier clic et relaché
	CPoint C1,C2;      
	// premier position pour bouger la zone
	CPoint BC1,BC2;
	//protected:
 
//	HGLOBAL * idpile;    
	  
	protected:
	// affichage d'info sur les traitement en cours
	BOOL				INFO_IMAGE;
	CDialInfoTrait *dialInfo;
	CString s1,s2,s3,s4,s5,s6;
	CString NomTraitement, CommentTraitement;
	CString ConcateneInfo(CString C1="",CString C2="",CString C3="",CString C4="",CString C5="",CString C6="");
	// c'est privé !!
	protected: 
	//pour l'affichage du Bitmap 
	BITMAPINFO	 			*	MyInitInfo; 
	BITMAPFILEHEADER * TheTete;

	//pour affichage de la palette
	BOOL 				PALETTE_AFFICHEE;
	CPaletteWnd 	*	pPaletteWnd;

	// Gestion de la palette courante
	enum TypePal {PALINIT,PALNB,PALCOUL};
	
	LOGPALETTE			* MyLogPaletteInit,*MyLogPaletteNB,*MyLogPaletteCOUL;
	CPalette 			* MyPalette;  
	public:
	TypePal TypePalette;

	void OnRestaureZone();

	protected:
	BOOL LaZone,LaCoupe,zoneSelectionnee,CoupeSelectionnee; 
	BOOL COUPE_VER,COUPE_HOR,COUPE_OBLI;       
	CFile LeFichier; 

	private:
	int HautWndInfo;    
	int LFene,HFene;   
	int PosHFen, PosVFen;
	//position du scrollbar
	int PosVScroll, PosHScroll;

	     
// **********************************************   
// ******* les fonctions membres ****************
// **********************************************
	
	public: 
	// constructeur vide (pour DECLARE_DYNAMIC)                  
	CImageWnd();
	CImageWnd(LONG ,LONG ,NBCOLOR InitColor=NB, BOOL SEQ = FALSE);    // larg, haut, NB ou COLOR 
	// constructeur utilisé pour l'ouverture d'un fichier
	CImageWnd(CString, long, long, TYPEIMAGE);
    // constructeur pour le chargement d'une image FLOAT ou BYTE  
    CImageWnd(long , long , TYPEIMAGE);
			
	// destructeur    
	~CImageWnd();                               
	BOOL Create(LPCSTR szTitle,LONG style = 0,CMDIFrameWnd* pParent = NULL);

	protected: // Gestion de l'impression
	void SetPrintAlign(CDC *pDC, HDC hdcPrn);

	protected:
	DECLARE_MESSAGE_MAP();
	
	protected:	
	//{{AFX_MSG(CImageWnd)
	afx_msg void OnFileMruFile1();
	afx_msg void OnFileMruFile2();
	afx_msg void OnFileMruFile3();
	afx_msg void OnFileMruFile4();
	afx_msg void OnPaint(void);
	afx_msg void OnSauve();  
	afx_msg void OnSelectionneZone(); 
	afx_msg void OnSelectionneZoneClavier();  
	afx_msg void OnRestaureImage();	   
	afx_msg void OnMiroirVert(); 
	afx_msg void OnMiroirH(); 
	afx_msg void OnRoberts();
	afx_msg void OnPrewitt();    
	afx_msg void OnSobel();
    afx_msg void OnMultiEchelle();
	afx_msg void OnFFT();
	afx_msg void OnSeuillerImage1();    
	afx_msg void OnSeuillerImage2();
	afx_msg void OnInversionVideo();    	
	afx_msg void OnCouleur();   
	afx_msg void OnPaletteAffiche(); 
	afx_msg void OnNoirEtBlanc();
	afx_msg void OnInformationGenerale() ; 
	afx_msg void OnSelectionneCoupeV();
	afx_msg void OnSelectionneCoupeOblique();
	afx_msg void OnSelectionneCoupeH();    
	afx_msg void OnHistogramme();  
	afx_msg void OnFonctionRepartition();  
	afx_msg void OnLireSequence();  
	afx_msg void OnSauveSequence();       
	afx_msg void OnDetruireSequence(); 
	afx_msg void OnLButtonDown(UINT nFlags,CPoint Point);  
	afx_msg void OnLButtonUp(UINT nFlags,CPoint Point);
	afx_msg void OnMouseMove(UINT nFlags,CPoint Point);     
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnSize( UINT nType, int cx, int cy ); 
	afx_msg void OnVScroll(UINT nSBCode, UINT nPos, CScrollBar * pScrollBar);   
	afx_msg void OnEtiquettage();
	afx_msg void OnRotationsDroite();
	afx_msg void OnRotationsGauche();
	afx_msg void OnSeparationRGB();
	afx_msg void OnNormalisation();
	afx_msg void OnHaarFloat();
	afx_msg void OnSplineFloat();
	afx_msg void OnGradRoberts();
	afx_msg void OnOpPrewitt();
	afx_msg void OnLaplacien1();
	afx_msg void OnLaplacien2();
	afx_msg void OnEchantLigneCol();
	afx_msg void OnBasse();
	afx_msg void OnCopieEnFloat();
	afx_msg void OnCopieEnByte();
	afx_msg void OnMiseAZero();
	afx_msg void OnLumchro();
	afx_msg void OnZoom100();
	afx_msg void OnZoom25();
	afx_msg void OnZoom50();
	afx_msg void OnHaar();
	afx_msg void OnBSpline();
	afx_msg void OnReconstruct();
	afx_msg void OnCont();
	afx_msg void OnSousEchanCol();
	afx_msg void OnSousEchanLigne();
	afx_msg void OnAffecterEchelle();
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnUpdateAffichePalette(CCmdUI* pCmdUI);
	afx_msg void OnMDIActivate(BOOL bActivate, CWnd* pActivateWnd, CWnd* pDeactivateWnd);
	afx_msg void OnUpdateNb(CCmdUI* pCmdUI);
	afx_msg void OnUpdateCouleur(CCmdUI* pCmdUI);
	afx_msg void OnDestroy();
	afx_msg void OnInfoImage();
	afx_msg void OnUpdateInfoImage(CCmdUI* pCmdUI);
	afx_msg void OnSauveExcel();
	afx_msg void OnUpdateZone(CCmdUI* pCmdUI);
	afx_msg void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnEditCopy();
	afx_msg void OnFilePrint();
	afx_msg void OnEditCut();
	afx_msg void OnEditPaste();
	afx_msg void OnMoyenneurxy121();
	afx_msg void OnMoyenneurxy111();
	afx_msg void OnMoyenneurx111();
	afx_msg void OnMoyenneury111();
	afx_msg void OnZoomX();
	afx_msg void OnExtraPaul();
	afx_msg void OnUpdateHorizontale(CCmdUI* pCmdUI);
	afx_msg void OnUpdateOblique(CCmdUI* pCmdUI);
	afx_msg void OnUpdateVerticale(CCmdUI* pCmdUI);
	afx_msg void OnSauverSous();
	afx_msg void OnClose();
	afx_msg void OnUpdateLireSequence(CCmdUI* pCmdUI);
	afx_msg void OnUpdateDetruireSequence(CCmdUI* pCmdUI);
	afx_msg void OnUpdateSauveSequence(CCmdUI* pCmdUI);
	afx_msg void OnGrab();
	afx_msg void OnUpdateGrab(CCmdUI* pCmdUI);
	afx_msg void OnSnap();
	afx_msg void OnUpdateSnap(CCmdUI* pCmdUI);
	afx_msg void OnShowWindow(BOOL bShow, UINT nStatus);
	afx_msg void OnTraitementMorphologieErosion();
	afx_msg void OnTraitementMorphologieDilatation();
	afx_msg void OnTraitementMorphologieOuverture();
	afx_msg void OnTraitementMorphologieFermeture();
	afx_msg void OnOutilsPseudorelief();
	afx_msg void OnTransformationConversionenbspline();
	afx_msg void OnTraitementNormalisationdescomposantes();
	afx_msg void OnTraitementPetrou();
	afx_msg void OnTraitementContourpetrou();
	afx_msg void OnTraitementConvolutiongaussien();
	afx_msg void OnTraitementConvolutiongaussienGaussien1dx();
	afx_msg void OnTraitementConvolutiongaussienGaussien1dy();
	afx_msg void OnTraitementConvolutiongaussienCurvilineardetection();
	afx_msg void OnTraitementCurvilinear();
	afx_msg void OnColorCurvi();
	afx_msg void OnRButtonDblClk(UINT nFlags, CPoint point);
	afx_msg void OnGradY();
	afx_msg void OnTraitementEstimationtempsfft();
	afx_msg void OnTraitementHypergradientParbloc();
	afx_msg void OnTraitementHypergradientFenetreglissante();
	afx_msg void OnTraitementHypergradientFenetreglissantezeropadding();
	afx_msg void OnTraitementOutilsEtudedistancemaxi();
	//}}AFX_MSG
	public:
	int CloseFen();
	// ******* FIN DE REPONSE DES MESSAGES ***************
	public:	   
	SImage  GetSImage();
	SImageFloat  GetSImageFloat();
	SImageRGB  GetSImageRGB();

	protected:
	void * AlloueMemoire(long,long);   
	float  * AlloueMemoireFloat(long , long);
	
	
	// permet de passer de l'image ImageFloat à Image avec une normalisation entre 0 et 255
	public:
	void  FloatToImage();
	void  ImageToFloat();                                         
	void Affichage();      
		
	protected:
	void Normalisation(long ,long ,long ,long ); // normalisation dans une zone
	    
	// affichage du BITAMP 
   
   // routine de mémorisation de l'image avant traitement
	void Memorise();
	BOOL IMAGE_MODIFIEE;
	// routine de traitement 
	void Histo(long TabHisto[256]);  
	void FonctionRepartition(long TabFRepar[256],long TabHisto[256]);  
		
	// fonctions du timer
	virtual void MettreLeContinu();
	virtual void ArreteLeContinu();    
	
	virtual void ChargeMenu();    
	// transformation écran -> image et image -> écran
	
	protected:
	CPoint Get_ij(CPoint point);   
	CPoint Get_ij(int x,int y);
	CPoint Inv_Get_ij(int i,int j);
	
	protected:
	
	virtual void DessineRectangle();
	void DessineRectangle(int , int , int ,int , COLORREF);
	void DessineLigne(long , long , long ,long , COLORREF );   
	void DessineLigne(long i1, long j1, long i2,long j2);
	
	protected:
	void MoveWindow( int , int , int , int ,BOOL);
					
	LONG OnSetCoupeVariable(UINT uiParam,LONG lparam );
	CProfilWnd * pCoupeVariable;
	int NbProfil;
	void MettreAJourLaCoupeActive();

    private:
	void SetSImage();    
	void SetSImageFloat();   
	void SetSImageRGB();
	
	
	// C'est privé !!
	private:
	void SourisRelachePourZone();   
	void SourisRelachePourBougeZone();
	void SourisBougePourZone(CPoint Point);  
	void SourisBougePourBougeZone(CPoint Point);
	void SourisRelachePourCoupe();    
	void SourisBougePourCoupe(CPoint Point);

	public:
	void AfficheSequence(int ,int );    
	void AlloueSequence(int *);
		
	private:
	void SetImage(unsigned char  *);     
 	void SetImage(rgb  *); 
 	
 	public:
	float GetPixelNB(long);

	private:
	void InitialisationPourBMP();
	void InitialisationPourFLOAT();
	void InitialisationPourMEM();
	void InitGlobal(long, long);

	private:			
	void TestMinMax(long& i1,long& j1,long& i2,long& j2); 
	void Troncature();
	void Mapping_Palette(LOGPALETTE*);                                 
	void Palette_Couleurs(); 
	void Palette_NB();  
	protected:	
	void AfficheNiveau(CPoint Point); 
	private:			
	void AfficheDimensions();              
	void AffichePositionCoupe();       
	  	
	void RemplitTableauIJCoupe(Pixel TabIJ[1000]);
	
// Fonctions utilisees par traitements standards
	void CalculLissage(float [3][3]); // Moyenneur

	protected:
	void Filtre_2D(SMasque masq,CString Chaine);
	void Filtre_2D_GX_ET_GY(SMasque masqGX,SMasque masqGY,CString Chaine);
	void Filtre_2D_GX_PUIS_GY(SMasque masqGX,SMasque masqGY,CString Chaine);	
	void Extraction(SFiltre h,SFiltre g,SImage SSource,SImage SDest,CString chaine);
	void Extraction(SFiltre h,SFiltre g,SImageFloat SSource,SImageFloat SDest,CString chaine);

public:
	afx_msg void OnOutilsSuividecontours();
	afx_msg void OnOutilsParamcurvidetect();
	int m_LargeurEtude;
	int m_LargeurFourier;
	int m_LongueurMiniForme;
	double m_SeuilReponse;
	int m_LargeurMini;
	//CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> TableauListesPointsFiltrees;
	//CList <PosDetectCurvi,PosDetectCurvi> * ListePointsFiltres;
	//CList <PosDetectCurvi,PosDetectCurvi> * ListePointsFiltresComp;

	//CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> TableauListesPointsFiltreesC2D;
	//CList <PosDetectCurvi,PosDetectCurvi> * ListePointsFiltresC2D;
	double VLocaleLargeur;
	afx_msg void OnMatchingSelection3points();
	BOOL MATCHINGCURVI;
	BOOL MATCHINGCURVI4;
	CStdioFile FichierMatching;
	CStdioFile FichierMatching4;
	int CptManipMatching;
	afx_msg void OnTransformationimageParamtransformation();
	double ParamAffine[7];
	double ParamHomogra[9];
	afx_msg void OnTransformationimageRun();
	afx_msg void OnMatchingCalculmatricepassage();
	CDialogAffine BoiteParamAffine;
	BoiteParamCurviDetect BoiteReglageParam;
	afx_msg void OnMatchingCalculmatricepassage504();
	afx_msg void OnMatchingSelection4points();
	afx_msg void OnTransformationimageRun506();
	void Convolution2DCedricCurvi(BYTE * ImageDst, BYTE * ImageSrc, double SigmaMasque,int Largeur, int Hauteur);
	void Filtre_2D(BYTE * ImageDst, BYTE * ImageSrc, SMasque masq,int Largeur, int Hauteur);
	long Xmin;
	long Xmax;
	long Ymin;
	long Ymax;
	CRgn * TraceRegionCurvi;
	rgb * ImageRep1;
	rgb * ImageRep2;
	void TireCouleurAlea(BYTE * r, BYTE * g, BYTE * b);
	void IncrementeCouleur(BYTE * r, BYTE * g, BYTE * b,long cpt);
	afx_msg void OnMatchingRepetabilit507();
	int IsAnEdge(int *Image,int Larg,CPoint A);
   void VisualiseCurvi(bool DebugEdge, bool DebugTracking, int Larg, int Haut, utls::BAry *edge, TableOfCurvilinearRegions &TableauListesPointsFiltrees, TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D);

/*
	void BuildImageFromList(int *Plan1,CList<PosDetectCurvi,PosDetectCurvi>* ListeCurvi);
	void MergeConnectedRegions(CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> & TableauListesPointsFiltreesC2D,CImageWnd * pImage);
	void TryToMerge(CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> & TableauListesPointsFiltreesC2D,CArray <CList<int,int>*,CList<int,int>*> &TabNumListesConnex,int *TabListeOk,CImageWnd * pImage);
	void TwoDimConstraints(CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> & TableauListesPointsFiltreesC2D,CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> & TableauListesPointsFiltrees,int m_LongueurMiniRegionFiltrage,int SeuilDepass);
	void SectionOptimization(CImageWnd * pImageTracking,CArray <CList<PosDetectCurvi,PosDetectCurvi>*,CList<PosDetectCurvi,PosDetectCurvi>*> & TableauListesPointsFiltrees);
	double Covariance(CArray <Sfpoint,Sfpoint> & Tab, Sfpoint Moy);
	double Covariance(Sfpoint *pTab, long NbPt, Sfpoint Moy);
	double Variance(CArray <Sfpoint,Sfpoint> & Tab, double Moy);
	double Variance(Sfpoint *pTab, long NbPt, double Moy);
	Sfpoint Projection(Sfpoint A, Sfpoint B,Sfpoint C,double *pr);
	double Scalaire(Sfpoint A, Sfpoint B,Sfpoint C);
	static int CompareDouble(const void *a,const void *b);
	void InvertTable(STabIndice *Tabi, STabIndice *Tabf,int Taille);
	double Distance(Sfpoint p1,Sfpoint p2);
	void PetrouCouleur(BYTE *pDest,rgb *pSource,int SeuilPetrou, int SigmaGaussien,ImagePetrou *I_x,ImagePetrou *I_y);
*/
public:
	afx_msg void OnTraitementPetroucouleur();
public:
	afx_msg void OnTraitementCurvicouleur();
public:
	afx_msg void OnTraitementGaussiencouleur();
};

#endif // __IMAGE_H__
