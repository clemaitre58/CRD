#ifndef __TIMAGE_H__     
#define __TIMAGE_H__  
#include "image.h"  
//#include "imagecla.h"
#include "GraphWnd.h"
#include "dialpara.h"

typedef unsigned char byte;
typedef unsigned short int word;
typedef unsigned int  dword;
// *****************************************************************          
// classe qui sert pour inclure des routines qui ne se rapportent
// pas à l'acquisition des images
// *****************************************************************
// cette classe DERIVE de CImageWnd                                
// ***************************************************************** 

class CVue3d;
class CParam;
class CImageClas;
class CDialLegend;

class CTraitImageWnd : public CImageWnd 
{  
	DECLARE_DYNAMIC(CTraitImageWnd)

	// constructeur vide (pour DECLARE_DYNAMIC)                  
	CTraitImageWnd();
    public:
	//  constructeur pour créer une image à partir de Larg et Haut et
	// du type : NB ou COLOR (par défaut NB)
	// SEQ = TRUE lorsque c'est une séquence
	CTraitImageWnd(LONG InitLarg,LONG InitHaut,NBCOLOR InitColor=NB, BOOL SEQ = FALSE);
	
    //  constructeur utilisé pour l'ouverture d'un fichier
	CTraitImageWnd(CString ,long , long , TYPEIMAGE );
	CTraitImageWnd(long , long , TYPEIMAGE );
		
	
	// destructeur
	~CTraitImageWnd();  


	protected:
	DECLARE_MESSAGE_MAP();
	
	public:  
		void Rotation(float);
		
		int PGCD(int,int);
		BOOL MasqueCGQuartmoy();
		BOOL MasquCGQuartier();
		BOOL MasqueCG();

		void MontrerLesVisuPar();
	//{{AFX_MSG(CTraitImageWnd)
	afx_msg void OnTraitementLuminancemoyenneseuil();
	afx_msg void OnTraitementLuminancemoyenne();
	afx_msg void OnTraitementSobel();
	afx_msg void OnTraitementFiltresett();
	afx_msg void OnTraitementFiltressettFiltresettnonnormalis();
	afx_msg void OnTraitementFiltressettHautbas();
	afx_msg void OnTraitementCorrelation();
	afx_msg void OnDst();
	afx_msg void OnContoursContourhystrsis1();
	afx_msg void OnContoursContourhystrsis2();
	afx_msg void OnContoursContourhystrsisfen88();
	afx_msg void OnTraitementReductionsur16niveaux();
	afx_msg void OnTraitementEgalisationdhistogramme();
	afx_msg void OnTraitementModificationdhistogrammeExpansiondynamique();
	afx_msg void OnEdition255();
	afx_msg void OnTraitementMaxmin();
	afx_msg void OnMasque5paramtres();
	afx_msg void OnParamtrageParamtres();
	afx_msg void OnUpdateParamtrageParamtres(CCmdUI* pCmdUI);
	afx_msg void OnParamtrageVisualisationparamtres();
	afx_msg void OnUpdateParamtrageVisualisationparamtres(CCmdUI* pCmdUI);
	afx_msg void OnParamtrageVisualisationmesoupol();
	afx_msg void OnUpdateParamtrageVisualisationmesoupol(CCmdUI* pCmdUI);
	afx_msg void OnMasqueRechercher();
	afx_msg void OnMasqueLamain();
	afx_msg void OnMasqueAnnuler();
	afx_msg void OnUpdateMasqueAnnuler(CCmdUI* pCmdUI);
	afx_msg void OnTraitementSupprimer();
	afx_msg void OnUpdateTraitementSupprimer(CCmdUI* pCmdUI);
	afx_msg void OnMasqueContraste();
	afx_msg void OnUpdateMasqueContraste(CCmdUI* pCmdUI);
	afx_msg void OnMasqueMaxmin();
	afx_msg void OnUpdateMasqueMaxmin(CCmdUI* pCmdUI);
	afx_msg void OnMasqueEntropie();
	afx_msg void OnUpdateMasqueEntropie(CCmdUI* pCmdUI);
	afx_msg void OnMasqueLuminancemoyenne();
	afx_msg void OnUpdateMasqueLuminancemoyenne(CCmdUI* pCmdUI);
	afx_msg void OnMasqueGradientmoyen();
	afx_msg void OnUpdateMasqueGradientmoyen(CCmdUI* pCmdUI);
	afx_msg void OnMasqueDimensions();
	afx_msg void OnUpdateMasqueDimensions(CCmdUI* pCmdUI);
	afx_msg void OnMasqueTous();
	afx_msg void OnUpdateMasqueTous(CCmdUI* pCmdUI);
	afx_msg void OnMasqueNbpixels();
	afx_msg void OnUpdateMasqueNbpixels(CCmdUI* pCmdUI);
	afx_msg void OnMasqueHistogramme();
	afx_msg void OnUpdateMasqueHistogramme(CCmdUI* pCmdUI);
	afx_msg void OnMasqueEssai();
	afx_msg void OnUpdateMasqueEssai(CCmdUI* pCmdUI);
	afx_msg void OnMasqueEssai2();
	afx_msg void OnUpdateMasqueEssai2(CCmdUI* pCmdUI);
	afx_msg void OnMasqueEssai3();
	afx_msg void OnUpdateMasqueEssai3(CCmdUI* pCmdUI);
	afx_msg void OnMasqueEssai4();
	afx_msg void OnMasqueEssai5();
	afx_msg void OnMasqueEssai6();
	afx_msg void OnMasqueDcoupage();
	afx_msg void OnUpdateMasqueDcoupage(CCmdUI* pCmdUI);
	afx_msg void OnCalculTexture();
	afx_msg void OnUpdateCalculTexture(CCmdUI* pCmdUI);
	afx_msg void OnDestroy();
	afx_msg void OnMasqueCalculslocauxLuminancemoyenne();
	afx_msg void OnUpdateMasqueCalculslocauxLuminancemoyenne(CCmdUI* pCmdUI);
	afx_msg void OnMasqueLocauxcontraste();
	afx_msg void OnUpdateMasqueLocauxcontraste(CCmdUI* pCmdUI);
	afx_msg void OnMasqueCalculslocauxMaxmin();
	afx_msg void OnUpdateMasqueCalculslocauxMaxmin(CCmdUI* pCmdUI);
	afx_msg void OnMasqueCalculslocauxEntropie();
	afx_msg void OnUpdateMasqueCalculslocauxEntropie(CCmdUI* pCmdUI);
	afx_msg void OnMasqueCalculslocauxGradientmoyen();
	afx_msg void OnUpdateMasqueCalculslocauxGradientmoyen(CCmdUI* pCmdUI);
	afx_msg void OnMasqueCalculslocauxTous();
	afx_msg void OnUpdateMasqueCalculslocauxTous(CCmdUI* pCmdUI);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnTraitementGradientderob();
	afx_msg void OnApprentissageGnrerunprojet();
	afx_msg void OnApprentissageDfinirlesparamtres();
	afx_msg void OnUpdateApprentissageDfinirlesparamtres(CCmdUI* pCmdUI);
	afx_msg void OnApprentissageGnrerlefichierdapprentissage();
	afx_msg void OnUpdateApprentissageGnrerlefichierdapprentissage(CCmdUI* pCmdUI);
	afx_msg void OnApprentissageConstruirelefichiermulticlasses();
	afx_msg void OnApprentissageGnrerunfichierpolytopesmulticlasses();
	afx_msg void OnDcisionConfigurationdeladcision();
	afx_msg void OnDcisionDcision();
	afx_msg void OnUpdateDcisionDcision(CCmdUI* pCmdUI);
	afx_msg void OnDcisionDcisionavecmasque();
	afx_msg void OnUpdateDcisionDcisionavecmasque(CCmdUI* pCmdUI);
	afx_msg void OnCaracteclair();
	afx_msg void OnMasqueCaractrisationglobal();
	afx_msg void OnDcisionDcision2masques();
	afx_msg void OnUpdateDcisionDcision2masques(CCmdUI* pCmdUI);
	afx_msg void OnTraitementSoustraction();
	afx_msg void OnMasqueMasquefixe();
	afx_msg void OnTraitementMax121232b();
	afx_msg void OnTraitementTestcontratse();
	afx_msg void OnTraitement5paramtrescommun();
	afx_msg void OnTraitement5paramtresetessais5paramtrescommunmasque();
	afx_msg void OnTraitement5paramtresetessaisEssaitot();
	afx_msg void OnTraitementRecalage();
	afx_msg void OnCalculellipse();
	afx_msg void OnDcisionDcision3masques();
	afx_msg void OnUpdateDcisionDcision3masques(CCmdUI* pCmdUI);
	afx_msg void OnDcisionDecision3masquesenregist();
	afx_msg void OnUpdateDcisionDecision3masquesenregist(CCmdUI* pCmdUI);
	afx_msg void OnDcisionDcision3masquesadapt();
	afx_msg void OnUpdateDcisionDcision3masquesadapt(CCmdUI* pCmdUI);
	afx_msg void OnTraitementCentre();
	afx_msg void OnDcisionDcision3masquescg();
	afx_msg void OnUpdateDcisionDcision3masquescg(CCmdUI* pCmdUI);
	afx_msg void OnDcisionDcision3masquessansrech();
	afx_msg void OnUpdateDcisionDcision3masquessansrech(CCmdUI* pCmdUI);
	afx_msg void OnTraitementDessinmesures();
	afx_msg void OnDcisionDecisionparzencircularit();
	afx_msg void OnUpdateDcisionDecisionparzencircularit(CCmdUI* pCmdUI);
	afx_msg void OnDcisionParzenDcisionenregistre();
	afx_msg void OnUpdateDcisionParzenDcisionenregistre(CCmdUI* pCmdUI);
	afx_msg void OnMasqueInversion();
	afx_msg void OnUpdateMasqueInversion(CCmdUI* pCmdUI);
	afx_msg void OnDcisionTest();
	afx_msg void OnUpdateDcisionTest(CCmdUI* pCmdUI);
	afx_msg void OnDcisionRemplissagedelespacedesparam();
	afx_msg void OnUpdateDcisionRemplissagedelespacedesparam(CCmdUI* pCmdUI);
	afx_msg void OnDcisionRemplissagedelespdesparamseq();
	afx_msg void OnUpdateDcisionRemplissagedelespdesparamseq(CCmdUI* pCmdUI);
	afx_msg void OnTraitementRductiondeniveaux();
	afx_msg void OnDcisionParzenImagecomplte();
	afx_msg void OnUpdateDcisionParzenImagecomplte(CCmdUI* pCmdUI);
	afx_msg void OnTraitementCorrelationCorrelation2();
	afx_msg void OnTraitementCodageextrema();
	afx_msg void OnTraitementDensitdanszone();
	afx_msg void OnTraitementRotation();
	afx_msg void OnTraitementOprateursunairesSeuillagedewen();
	afx_msg void OnTraitementMorphologieBarbulage();
	afx_msg void OnTraitementMorphologieAmincissement();
	afx_msg void OnTraitementMorphologieEpaississement();
	afx_msg void OnDcisionDecisionsurfloat();
	afx_msg void OnDcisionDecidesurfloatdebug();
	afx_msg void OnDcisionDecidesvm();
	afx_msg void OnDcisionGenerefichierdat();
	afx_msg void OnTraitementComparaisonimbin();
	afx_msg void OnTraitementAjustement();
	afx_msg void OnTraitementAnalysecontour();
	afx_msg void OnTraitementEquilibrage();
	afx_msg void OnDcisionGeneredata();
	//}}AFX_MSG
	// autre traitement   

//*****************************************************//
//		FONCTIONS MEMBRES de la Classe CTRAITIMAGEWND	   //
//*****************************************************//
	protected:
		void ajusty(int l1,int l2,int c,byte *point,byte *bram);
		void ajustx(int c1,int c2,int l,byte *point,byte *bram,BOOL);
		void InitialiseObjetCTrait();
		void MaximumGlissant1212(unsigned char *);
		void MinimumGlissant1212(unsigned char *);
		void Ouverture1212();
		void Fermeture1212();
		void GradientGlissant1212();
		void ContrasteGlissant1212();
		void LuminanceGlissante1212();
		void FenetreGlissanteRapide2(unsigned char*,unsigned __int16*,
								unsigned __int16*,int ,int ,int ,int);
		void GradientRapide2(unsigned __int32 *,unsigned __int32*,int,int);
		void GradientRapide3(unsigned __int32 *lram,unsigned __int32 *bram);
		void GradientRapide4(unsigned __int32 *lram,unsigned __int32 *bram,int* debutl,int*finl,int jl);
		void GradientRapide5(unsigned __int32 *lram,unsigned __int32 *bram);

		void FenetreGlissanteRapide(unsigned char*,unsigned __int16*,
								unsigned __int16*,int ,int ,int ,int);
		void GradientRapide(unsigned __int32 *,int ,int );
		long FenetreGlissanteRapideMoy(unsigned char* ,unsigned __int16* ,
											unsigned __int16* ,int ,int ,int ,int);
		long FenetreGlissanteRapideMoy2(unsigned char* ,unsigned __int16* ,
											unsigned __int16* ,int ,int ,int ,int );
		void fenetregliscombi(unsigned char* aram,unsigned char* bram,unsigned __int16* cufen,
											unsigned __int16* lufen,unsigned __int16* gcufen,
											unsigned __int16* glufen);
		void fenetregliscombi2(unsigned char* aram,unsigned char* bram,unsigned __int16* cufen,
											unsigned __int16* lufen,unsigned __int16* gcufen,
											unsigned __int16* glufen,int*debutl,int*finl,int*debuth,int*finh,int jl,int ih);
		void fenetregliscombi3(unsigned char* aram,unsigned char* bram,unsigned __int16* cufen,
											unsigned __int16* lufen,unsigned __int16* acufen,
											unsigned __int16* alufen);


	/////////////polytopes multi-classes
	protected:
		float MasqEucl(int&);
		BOOL MasqCGMin();
		BOOL OnCalculellipseAutom();
		void decoupageZone(int m_seuil=40);
		void CacherLesVisuPar();
		BOOL SauveParam(CString);
		BOOL ChargementDesParam(CString);
		BOOL CalculParam(int);
		void Calc(long p);
		void AlloueParam();
		LONG OnFinClasse(UINT uiParam,LONG lparam );
		LONG OnFinVisuMes(UINT uiParam,LONG lparam );

	/////////////masque
		void AffectRegionMasque();
		void OnMasqueRechercheAutom(int);
		float OnMasqueLuminancemoyenneAutom();
		virtual void DessineRectangle();

//*****************************************************//
//		DONNEES MEMBRES de la Classe CTRAITIMAGEWND	   //
//*****************************************************//
	////contour///////////
	int DEBUTL[303],FINL[303],DEBUTH[291],FINH[291],JL,IH;
	unsigned char Seuilcont1,Seuilcont2;
	long SeuillongContour, LongmaxContour;

	///masque
	BOOL Masque,Decoupage;

	BYTE Masquergn[296*308];//[272*272];
	////menu masque
	int PositionMenuMasque;
 
	///vue3d
	CVue3d *PVue3d;
	BOOL VUE3D;
	CRgn RegionEllipse;
	
	//////////polytopes multi-classes/////

	CString NomFichProj,NomFichAppr,NomFichParam,NomImageAppr;	
	CString NomFichMMC,NomFichCMC,NomFichPMC;
	CImageWnd *VisuPar[256];
    CParam *PParametre;
    CDialPara DialPara;
	CAfficheMesPol *AffMesPol;

	BOOL ParamCalcul;
    BOOL Visu;
	BOOL ValideMenu;
	BOOL AffectNonClas;
	BOOL Visumes;
	BOOL ValideDecision,ValideVisuMes;

	int CoeffPoly;

	int ClassAffect;

	byte Buffer[256];//*Buffer;
    BYTE NbParam,NBPMAX;

	CImageClas* ImageAppr;
	BOOL DefImagAppr;
	CDialLegend *Legende;
	/////////////////////
	//////Parzen 3D//////
	/////////////////////
//	byte EspClas[256*256*256];

};  
          

#endif // __TIMAGE_H__
