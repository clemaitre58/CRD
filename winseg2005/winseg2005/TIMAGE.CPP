#include "stdafx.h"
#include "timage.h"
#include "imagecla.h"
//#include <mmsystem.h> inclure winmm.lib 

#include "AfficheMesPol.h"
#include "Vue3d.h"

#include "param.h"
#include "polmc.h"
#include "polmcfloat.h"


#include "dialogue.h"
#include "DialSetT.h"
#include "DialContour.h"
#include "ProgDlg.h"
#include "DialProjAppr.h"
#include "DialConfigPoly.h"
#include "DialNbPointAppr.h"
#include "DialPmc.h"
#include "DialConfigDec.h"
#include "DialLegend.h"

#pragma warning( disable : 4996 )

//////////////////////////////////////////////
/////////variables globales////////////////////
typedef unsigned char byte;

using namespace utls;

const UINT uiMsg=::RegisterWindowMessage("WM_FINCLASSE");
const UINT uiMsg2=::RegisterWindowMessage("WM_FINVISUMES");
//const UINT uiMsg4=::RegisterWindowMessage("WM_FINVUE3D");


unsigned char TABABSO[0x10000];
unsigned char CONTRAST[0x10000];
unsigned char *Inittampon;//[272*272];
unsigned __int16 Initcufen[272*272],Initllufen[272*272],Initlgradfen[272*272];
unsigned __int16 ESSAI[400*400];


int AJUST=8;
int AJUSTH=12;
int AJUSTD=0;
int AJUSTG=0;
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!
// déclaration des messages!
// !!!!!!!!!!!!!!!!!!!!!!!!!
BEGIN_MESSAGE_MAP(CTraitImageWnd, CImageWnd)
	ON_WM_CREATE()
	// fonctions déclarées dans CImage
	//{{AFX_MSG_MAP(CTraitImageWnd)    
	ON_COMMAND(ID_TRAITEMENT_LUMINANCEMOYENNESEUIL, OnTraitementLuminancemoyenneseuil)
	ON_COMMAND(ID_TRAITEMENT_LUMINANCEMOYENNE, OnTraitementLuminancemoyenne)
	ON_COMMAND(ID_TRAITEMENT_SOBEL, OnTraitementSobel)
	ON_COMMAND(ID_TRAITEMENT_FILTRESETT, OnTraitementFiltresett)
	ON_COMMAND(ID_TRAITEMENT_FILTRESSETT_FILTRESETTNONNORMALIS, OnTraitementFiltressettFiltresettnonnormalis)
	ON_COMMAND(ID_TRAITEMENT_FILTRESSETT_HAUTBAS, OnTraitementFiltressettHautbas)
	ON_COMMAND(ID_TRAITEMENT_CORRELATION, OnTraitementCorrelation)
	ON_COMMAND(ID_DST, OnDst)
	ON_COMMAND(ID_CONTOURS_CONTOURHYSTRSIS1, OnContoursContourhystrsis1)
	ON_COMMAND(ID_CONTOURS_CONTOURHYSTRSIS2, OnContoursContourhystrsis2)
	ON_COMMAND(ID_CONTOURS_CONTOURHYSTRSISFEN88, OnContoursContourhystrsisfen88)
	ON_COMMAND(ID_TRAITEMENT_REDUCTIONSUR16NIVEAUX, OnTraitementReductionsur16niveaux)
	ON_COMMAND(ID_TRAITEMENT_MODIFICATIONDHISTOGRAMME_EGALISATIONDHISTOGRAMME, OnTraitementEgalisationdhistogramme)
	ON_COMMAND(ID_TRAITEMENT_MODIFICATIONDHISTOGRAMME_EXPANSIONDYNAMIQUE, OnTraitementModificationdhistogrammeExpansiondynamique)
	ON_COMMAND(ID_EDITION_255, OnEdition255)
	ON_COMMAND(ID_TRAITEMENT_MAXMIN, OnTraitementMaxmin)
	ON_COMMAND(ID_MASQUE_5PARAMTRES, OnMasque5paramtres)
	ON_COMMAND(ID_PARAMTRAGE_PARAMTRES, OnParamtrageParamtres)
	ON_UPDATE_COMMAND_UI(ID_PARAMTRAGE_PARAMTRES, OnUpdateParamtrageParamtres)
	ON_COMMAND(ID_PARAMTRAGE_VISUALISATIONPARAMTRES, OnParamtrageVisualisationparamtres)
	ON_UPDATE_COMMAND_UI(ID_PARAMTRAGE_VISUALISATIONPARAMTRES, OnUpdateParamtrageVisualisationparamtres)
	ON_COMMAND(ID_PARAMTRAGE_VISUALISATIONMESOUPOL, OnParamtrageVisualisationmesoupol)
	ON_UPDATE_COMMAND_UI(ID_PARAMTRAGE_VISUALISATIONMESOUPOL, OnUpdateParamtrageVisualisationmesoupol)
	ON_COMMAND(ID_MASQUE_RECHERCHER, OnMasqueRechercher)
	ON_COMMAND(ID_MASQUE_LAMAIN, OnMasqueLamain)
	ON_COMMAND(ID_MASQUE_ANNULER, OnMasqueAnnuler)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_ANNULER, OnUpdateMasqueAnnuler)
	ON_COMMAND(ID_TRAITEMENT_SUPPRIMER, OnTraitementSupprimer)
	ON_UPDATE_COMMAND_UI(ID_TRAITEMENT_SUPPRIMER, OnUpdateTraitementSupprimer)
	ON_COMMAND(ID_MASQUE_CONTRASTE, OnMasqueContraste)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_CONTRASTE, OnUpdateMasqueContraste)
	ON_COMMAND(ID_MASQUE_MAXMIN, OnMasqueMaxmin)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_MAXMIN, OnUpdateMasqueMaxmin)
	ON_COMMAND(ID_MASQUE_ENTROPIE, OnMasqueEntropie)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_ENTROPIE, OnUpdateMasqueEntropie)
	ON_COMMAND(ID_MASQUE_LUMINANCEMOYENNE, OnMasqueLuminancemoyenne)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_LUMINANCEMOYENNE, OnUpdateMasqueLuminancemoyenne)
	ON_COMMAND(ID_MASQUE_GRADIENTMOYEN, OnMasqueGradientmoyen)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_GRADIENTMOYEN, OnUpdateMasqueGradientmoyen)
	ON_COMMAND(ID_MASQUE_DIMENSIONS, OnMasqueDimensions)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_DIMENSIONS, OnUpdateMasqueDimensions)
	ON_COMMAND(ID_MASQUE_TOUS, OnMasqueTous)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_TOUS, OnUpdateMasqueTous)
	ON_COMMAND(ID_MASQUE_NBPIXELS, OnMasqueNbpixels)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_NBPIXELS, OnUpdateMasqueNbpixels)
	ON_COMMAND(ID_MASQUE_HISTOGRAMME, OnMasqueHistogramme)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_HISTOGRAMME, OnUpdateMasqueHistogramme)
	ON_COMMAND(ID_MASQUE_ESSAI, OnMasqueEssai)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_ESSAI, OnUpdateMasqueEssai)
	ON_COMMAND(ID_MASQUE_ESSAI2, OnMasqueEssai2)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_ESSAI2, OnUpdateMasqueEssai2)
	ON_COMMAND(ID_MASQUE_ESSAI3, OnMasqueEssai3)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_ESSAI3, OnUpdateMasqueEssai3)
	ON_COMMAND(ID_MASQUE_ESSAI4, OnMasqueEssai4)
	ON_COMMAND(ID_MASQUE_ESSAI5, OnMasqueEssai5)
	ON_COMMAND(ID_MASQUE_DCOUPAGE, OnMasqueDcoupage)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_DCOUPAGE, OnUpdateMasqueDcoupage)
	ON_COMMAND(ID_CALCUL_TEXTURE, OnCalculTexture)
	ON_UPDATE_COMMAND_UI(ID_CALCUL_TEXTURE, OnUpdateCalculTexture)
	ON_WM_DESTROY()
	ON_COMMAND(ID_MASQUE_CALCULSLOCAUX_LUMINANCEMOYENNE, OnMasqueCalculslocauxLuminancemoyenne)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_CALCULSLOCAUX_LUMINANCEMOYENNE, OnUpdateMasqueCalculslocauxLuminancemoyenne)
	ON_COMMAND(ID_MASQUE_LOCAUXCONTRASTE, OnMasqueLocauxcontraste)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_LOCAUXCONTRASTE, OnUpdateMasqueLocauxcontraste)
	ON_COMMAND(ID_MASQUE_CALCULSLOCAUX_MAXMIN, OnMasqueCalculslocauxMaxmin)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_CALCULSLOCAUX_MAXMIN, OnUpdateMasqueCalculslocauxMaxmin)
	ON_COMMAND(ID_MASQUE_CALCULSLOCAUX_ENTROPIE, OnMasqueCalculslocauxEntropie)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_CALCULSLOCAUX_ENTROPIE, OnUpdateMasqueCalculslocauxEntropie)
	ON_COMMAND(ID_MASQUE_CALCULSLOCAUX_GRADIENTMOYEN, OnMasqueCalculslocauxGradientmoyen)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_CALCULSLOCAUX_GRADIENTMOYEN, OnUpdateMasqueCalculslocauxGradientmoyen)
	ON_COMMAND(ID_MASQUE_CALCULSLOCAUX_TOUS, OnMasqueCalculslocauxTous)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_CALCULSLOCAUX_TOUS, OnUpdateMasqueCalculslocauxTous)
	ON_WM_CREATE()
	ON_COMMAND(ID_TRAITEMENT_GRADIENTDEROB, OnTraitementGradientderob)
	ON_COMMAND(ID_APPRENTISSAGE_GNRERUNPROJET, OnApprentissageGnrerunprojet)
	ON_COMMAND(ID_APPRENTISSAGE_DFINIRLESPARAMTRES, OnApprentissageDfinirlesparamtres)
	ON_UPDATE_COMMAND_UI(ID_APPRENTISSAGE_DFINIRLESPARAMTRES, OnUpdateApprentissageDfinirlesparamtres)
	ON_COMMAND(ID_APPRENTISSAGE_GNRERLEFICHIERDAPPRENTISSAGE, OnApprentissageGnrerlefichierdapprentissage)
	ON_UPDATE_COMMAND_UI(ID_APPRENTISSAGE_GNRERLEFICHIERDAPPRENTISSAGE, OnUpdateApprentissageGnrerlefichierdapprentissage)
	ON_COMMAND(ID_APPRENTISSAGE_CONSTRUIRELEFICHIERMULTICLASSES, OnApprentissageConstruirelefichiermulticlasses)
	ON_COMMAND(ID_APPRENTISSAGE_GNRERUNFICHIERPOLYTOPESMULTICLASSES, OnApprentissageGnrerunfichierpolytopesmulticlasses)
	ON_COMMAND(ID_DCISION_CONFIGURATIONDELADCISION, OnDcisionConfigurationdeladcision)
	ON_COMMAND(ID_DCISION_DCISION, OnDcisionDcision)
	ON_UPDATE_COMMAND_UI(ID_DCISION_DCISION, OnUpdateDcisionDcision)
	ON_COMMAND(ID_DCISION_DCISIONAVECMASQUE, OnDcisionDcisionavecmasque)
	ON_UPDATE_COMMAND_UI(ID_DCISION_DCISIONAVECMASQUE, OnUpdateDcisionDcisionavecmasque)
	ON_COMMAND(ID_CARACTECLAIR, OnCaracteclair)
	ON_COMMAND(ID_MASQUE_CARACTRISATIONGLOBAL, OnMasqueCaractrisationglobal)
	ON_COMMAND(ID_DCISION_DCISION2MASQUES, OnDcisionDcision2masques)
	ON_UPDATE_COMMAND_UI(ID_DCISION_DCISION2MASQUES, OnUpdateDcisionDcision2masques)
	ON_COMMAND(ID_TRAITEMENT_SOUSTRACTION, OnTraitementSoustraction)
	ON_COMMAND(ID_MASQUE_MASQUEFIXE, OnMasqueMasquefixe)
	ON_COMMAND(ID_TRAITEMENT_MAX121232B, OnTraitementMax121232b)
	ON_COMMAND(ID_TRAITEMENT_TESTCONTRATSE, OnTraitementTestcontratse)
	ON_COMMAND(ID_TRAITEMENT_5PARAMTRESCOMMUN, OnTraitement5paramtrescommun)
	ON_COMMAND(ID_TRAITEMENT_5PARAMTRESETESSAIS_5PARAMTRESCOMMUNMASQUE, OnTraitement5paramtresetessais5paramtrescommunmasque)
	ON_COMMAND(ID_TRAITEMENT_5PARAMTRESETESSAIS_ESSAITOT, OnTraitement5paramtresetessaisEssaitot)
	ON_COMMAND(ID_TRAITEMENT_RECALAGE, OnTraitementRecalage)
	ON_COMMAND(ID_CALCULELLIPSE, OnCalculellipse)
	ON_COMMAND(ID_DCISION_DCISION3MASQUES, OnDcisionDcision3masques)
	ON_UPDATE_COMMAND_UI(ID_DCISION_DCISION3MASQUES, OnUpdateDcisionDcision3masques)
	ON_COMMAND(ID_DCISION_DECISION3MASQUESENREGIST, OnDcisionDecision3masquesenregist)
	ON_UPDATE_COMMAND_UI(ID_DCISION_DECISION3MASQUESENREGIST, OnUpdateDcisionDecision3masquesenregist)
	ON_COMMAND(ID_DCISION_DCISION3MASQUESADAPT, OnDcisionDcision3masquesadapt)
	ON_UPDATE_COMMAND_UI(ID_DCISION_DCISION3MASQUESADAPT, OnUpdateDcisionDcision3masquesadapt)
	ON_COMMAND(ID_TRAITEMENT_CENTRE, OnTraitementCentre)
	ON_COMMAND(ID_DCISION_DCISION3MASQUESCG, OnDcisionDcision3masquescg)
	ON_UPDATE_COMMAND_UI(ID_DCISION_DCISION3MASQUESCG, OnUpdateDcisionDcision3masquescg)
	ON_COMMAND(ID_DCISION_DCISION3MASQUESSANSRECH, OnDcisionDcision3masquessansrech)
	ON_UPDATE_COMMAND_UI(ID_DCISION_DCISION3MASQUESSANSRECH, OnUpdateDcisionDcision3masquessansrech)
	ON_COMMAND(ID_TRAITEMENT_DESSINMESURES, OnTraitementDessinmesures)
	ON_COMMAND(ID_DCISION_DECISIONPARZENCIRCULARIT, OnDcisionDecisionparzencircularit)
	ON_UPDATE_COMMAND_UI(ID_DCISION_DECISIONPARZENCIRCULARIT, OnUpdateDcisionDecisionparzencircularit)
	ON_COMMAND(ID_DCISION_PARZEN_DCISIONENREGISTRE, OnDcisionParzenDcisionenregistre)
	ON_UPDATE_COMMAND_UI(ID_DCISION_PARZEN_DCISIONENREGISTRE, OnUpdateDcisionParzenDcisionenregistre)
	ON_COMMAND(ID_MASQUE_INVERSION, OnMasqueInversion)
	ON_UPDATE_COMMAND_UI(ID_MASQUE_INVERSION, OnUpdateMasqueInversion)
	ON_COMMAND(ID_DCISION_TEST, OnDcisionTest)
	ON_UPDATE_COMMAND_UI(ID_DCISION_TEST, OnUpdateDcisionTest)
	ON_COMMAND(ID_DCISION_REMPLISSAGEDELESPACEDESPARAM, OnDcisionRemplissagedelespacedesparam)
	ON_UPDATE_COMMAND_UI(ID_DCISION_REMPLISSAGEDELESPACEDESPARAM, OnUpdateDcisionRemplissagedelespacedesparam)
	ON_COMMAND(ID_DCISION_REMPLISSAGEDELESPDESPARAMSEQ, OnDcisionRemplissagedelespdesparamseq)
	ON_UPDATE_COMMAND_UI(ID_DCISION_REMPLISSAGEDELESPDESPARAMSEQ, OnUpdateDcisionRemplissagedelespdesparamseq)
	ON_COMMAND(ID_TRAITEMENT_RDUCTIONDENIVEAUX, OnTraitementRductiondeniveaux)
	ON_COMMAND(ID_DCISION_PARZEN_IMAGECOMPLTE, OnDcisionParzenImagecomplte)
	ON_UPDATE_COMMAND_UI(ID_DCISION_PARZEN_IMAGECOMPLTE, OnUpdateDcisionParzenImagecomplte)
	ON_COMMAND(ID_TRAITEMENT_CORRELATION_CORRELATION2, OnTraitementCorrelationCorrelation2)
	ON_COMMAND(ID_TRAITEMENT_CODAGEEXTREMA, OnTraitementCodageextrema)
	ON_COMMAND(ID_TRAITEMENT_DENSITDANSZONE, OnTraitementDensitdanszone)
	ON_COMMAND(ID_TRAITEMENT_ROTATION, OnTraitementRotation)
	ON_COMMAND(ID_TRAITEMENT_OPRATEURSUNAIRES_SEUILLAGEDEWEN, OnTraitementOprateursunairesSeuillagedewen)
	ON_COMMAND(ID_TRAITEMENT_MORPHOLOGIE_BARBULAGE, OnTraitementMorphologieBarbulage)
	ON_COMMAND(ID_TRAITEMENT_MORPHOLOGIE_AMINCISSEMENT, OnTraitementMorphologieAmincissement)
	ON_COMMAND(ID_TRAITEMENT_MORPHOLOGIE_EPAISSISSEMENT, OnTraitementMorphologieEpaississement)
	ON_COMMAND(ID_DCISION_DECISIONSURFLOAT, OnDcisionDecisionsurfloat)
	ON_COMMAND(ID_DCISION_DECIDESURFLOATDEBUG, OnDcisionDecidesurfloatdebug)
	ON_COMMAND(ID_DCISION_DECIDESVM, OnDcisionDecidesvm)
	ON_COMMAND(ID_DCISION_GENEREFICHIERDAT, OnDcisionGenerefichierdat)
	ON_COMMAND(ID_TRAITEMENT_COMPARAISONIMBIN, OnTraitementComparaisonimbin)
	ON_COMMAND(ID_TRAITEMENT_AJUSTEMENT, OnTraitementAjustement)
	ON_COMMAND(ID_TRAITEMENT_ANALYSECONTOUR, OnTraitementAnalysecontour)
	ON_COMMAND(ID_TRAITEMENT_EQUILIBRAGE, OnTraitementEquilibrage)
	ON_REGISTERED_MESSAGE(uiMsg, OnFinClasse)
	ON_REGISTERED_MESSAGE(uiMsg2, OnFinVisuMes)
	ON_COMMAND(ID_DCISION_GENEREDATA, OnDcisionGeneredata)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()                            

IMPLEMENT_DYNAMIC(CTraitImageWnd, CImageWnd)
                   
// contructeur vide                   
CTraitImageWnd::CTraitImageWnd():CImageWnd()
{
	MessageBox("Vous avez ouvert un objet VIDE","ERREUR", MB_ICONEXCLAMATION);
}
// acquisition
CTraitImageWnd::CTraitImageWnd(LONG InitLarg,LONG InitHaut,NBCOLOR InitColor,BOOL SEQ)
					:CImageWnd( InitLarg, InitHaut,InitColor,SEQ)
{
	InitialiseObjetCTrait();
}

// constructeur appelé pour ouvrir une image à partir d'un fichier.                             
CTraitImageWnd::CTraitImageWnd(CString InitFileName,long InitLarg, long InitHaut, TYPEIMAGE InitLeType)
		:CImageWnd( InitFileName, InitLarg,  InitHaut,InitLeType)
{
	InitialiseObjetCTrait();
}


// constructeur appelé pour ouvrir une image résultat d'un traitement.
CTraitImageWnd::CTraitImageWnd(long InitLarg, long InitHaut, TYPEIMAGE InitLeType)
		:CImageWnd( InitLarg,  InitHaut,InitLeType)
{
	InitialiseObjetCTrait();
}

  
CTraitImageWnd::~CTraitImageWnd()
{                                
// rajoute ici
}

// ****   Initialisation des donnes membres de la classe
void CTraitImageWnd::InitialiseObjetCTrait()
{
	SetCursor(LoadCursor(NULL,IDC_WAIT));
	/////////contours
	Seuilcont1=128;
	Seuilcont2=128;
	LongmaxContour=1000;
	SeuillongContour=50;
////////////////////////////////caract eclairage
/*	CPoint h,sup[3],inf[3];
	CSize difference;
	CRect rect,rect2,rect3;
	CRgn trianglesup,triangleinf,rgntemp,rgntemp1,rgntemp2,rgntemp3,ellipse1,rgn20,rgn21,
		rgn30,rgn31,rgn32,rgn33,rgn40,rgn41,rgn42,rgn43,rgn44,rgn45,rgn46,rgn47;


//////////////////************////////////////**************/////////////
	///////////		recherche et affectation des différentes zones
//////////////////************////////////////**************/////////////
/*	debut_i=12;
	debut_j=24;
	fin_i=269;
	fin_j=281;
	AffectRegionMasque();//trouve la grande ellipse

	RegionEllipse.GetRgnBox(rect);
	sup[0]=rect.CenterPoint();
	sup[1]=rect.TopLeft();
	sup[2].y=rect.TopLeft().y;
	sup[2].x=rect.BottomRight().x;
	inf[0]=rect.CenterPoint();
	inf[1]=rect.BottomRight();
	inf[2].x=rect.TopLeft().x;
	inf[2].y=rect.BottomRight().y;
	trianglesup.CreatePolygonRgn( sup, 3, ALTERNATE );
	triangleinf.CreatePolygonRgn( inf, 3, ALTERNATE );


	rgn40.CreateEllipticRgnIndirect(rect);
	rgn41.CreateEllipticRgnIndirect(rect);
	rgn42.CreateEllipticRgnIndirect(rect);
	rgn43.CreateEllipticRgnIndirect(rect);
	rgn44.CreateEllipticRgnIndirect(rect);
	rgn45.CreateEllipticRgnIndirect(rect);
	rgn46.CreateEllipticRgnIndirect(rect);
	rgn47.CreateEllipticRgnIndirect(rect);

	//calcul la 3eme ellipse
	rect2=rect;
	difference.cx=(long)(rect.Width()*0.15843497);
	difference.cy=(long)(rect.Height()*0.15843497);
	rect2.DeflateRect(difference);
	rgn30.CreateEllipticRgnIndirect(rect2);
	rgn31.CreateEllipticRgnIndirect(rect2);
	rgn32.CreateEllipticRgnIndirect(rect2);
	rgn33.CreateEllipticRgnIndirect(rect2);

	//calcul la 2eme ellipse
	rect2=rect;
	difference.cx=(long)(rect.Width()*0.2763932);
	difference.cy=(long)(rect.Height()*0.2763932);
	rect2.DeflateRect(difference);
	rgn20.CreateEllipticRgnIndirect(rect2);
	rgn21.CreateEllipticRgnIndirect(rect2);

	//calcul l'ellipse centrale
	rect2=rect;
	difference.cx=(long)(rect.Width()*0.3709005);
	difference.cy=(long)(rect.Height()*0.3709005);
	rect2.DeflateRect(difference);
	ellipse1.CreateEllipticRgnIndirect(rect2);

	//calcul les disques
	rgn40.CombineRgn(&rgn40,&rgn30,RGN_DIFF);
	rgn41.CombineRgn(&rgn41,&rgn30,RGN_DIFF);
	rgn42.CombineRgn(&rgn42,&rgn30,RGN_DIFF);
	rgn43.CombineRgn(&rgn43,&rgn30,RGN_DIFF);
	rgn44.CombineRgn(&rgn44,&rgn30,RGN_DIFF);
	rgn45.CombineRgn(&rgn45,&rgn30,RGN_DIFF);
	rgn46.CombineRgn(&rgn46,&rgn30,RGN_DIFF);
	rgn47.CombineRgn(&rgn47,&rgn30,RGN_DIFF);

	rgn30.CombineRgn( &rgn30, &rgn20, RGN_DIFF );
	rgn31.CombineRgn( &rgn31, &rgn20, RGN_DIFF );
	rgn32.CombineRgn( &rgn32, &rgn20, RGN_DIFF );
	rgn33.CombineRgn( &rgn33, &rgn20, RGN_DIFF );
	rgn20.CombineRgn( &rgn20, &ellipse1, RGN_DIFF );
	rgn21.CombineRgn( &rgn21, &ellipse1, RGN_DIFF );;

	//calcul les parties sup et inf
	rect3=rect2=rect;
	rect2.OffsetRect(0,int(rect.Height()/2-0.5));
	rect3.OffsetRect(0,-rect.Height()/2);
	rgntemp.CreateRectRgnIndirect(rect2);
	rgntemp1.CreateRectRgnIndirect(rect3);

	rgn21.CombineRgn(&rgn21,&rgntemp,RGN_DIFF);
	rgn32.CombineRgn(&rgn32,&rgntemp,RGN_DIFF);
	rgn33.CombineRgn(&rgn33,&rgntemp,RGN_DIFF);
	rgn44.CombineRgn(&rgn44,&rgntemp,RGN_DIFF);
	rgn45.CombineRgn(&rgn45,&rgntemp,RGN_DIFF);
	rgn46.CombineRgn(&rgn46,&rgntemp,RGN_DIFF);
	rgn47.CombineRgn(&rgn47,&rgntemp,RGN_DIFF);


	rgn20.CombineRgn(&rgn20,&rgntemp1,RGN_DIFF);
	rgn30.CombineRgn(&rgn30,&rgntemp1,RGN_DIFF);
	rgn31.CombineRgn(&rgn31,&rgntemp1,RGN_DIFF);
	rgn40.CombineRgn(&rgn40,&rgntemp1,RGN_DIFF);
	rgn41.CombineRgn(&rgn41,&rgntemp1,RGN_DIFF);
	rgn42.CombineRgn(&rgn42,&rgntemp1,RGN_DIFF);
	rgn43.CombineRgn(&rgn43,&rgntemp1,RGN_DIFF);
	
	//calcul les parties droites et gauches
	rect3=rect2=rect;
	rect2.OffsetRect(int(rect.Width()/2-0.5),0);
	rect3.OffsetRect(-rect.Width()/2,0);
	rgntemp2.CreateRectRgnIndirect(rect2);
	rgntemp3.CreateRectRgnIndirect(rect3);

	rgn31.CombineRgn(&rgn31,&rgntemp2,RGN_DIFF);
	rgn32.CombineRgn(&rgn32,&rgntemp2,RGN_DIFF);
	rgn42.CombineRgn(&rgn42,&rgntemp2,RGN_DIFF);
	rgn43.CombineRgn(&rgn43,&rgntemp2,RGN_DIFF);
	rgn44.CombineRgn(&rgn44,&rgntemp2,RGN_DIFF);
	rgn45.CombineRgn(&rgn45,&rgntemp2,RGN_DIFF);
	
	rgn30.CombineRgn(&rgn30,&rgntemp3,RGN_DIFF);
	rgn33.CombineRgn(&rgn33,&rgntemp3,RGN_DIFF);
	rgn40.CombineRgn(&rgn40,&rgntemp3,RGN_DIFF);
	rgn41.CombineRgn(&rgn41,&rgntemp3,RGN_DIFF);
	rgn46.CombineRgn(&rgn46,&rgntemp3,RGN_DIFF);
	rgn47.CombineRgn(&rgn47,&rgntemp3,RGN_DIFF);

	//soustraction des triangles
	rgn40.CombineRgn(&rgn40,&triangleinf,RGN_DIFF);
	rgn42.CombineRgn(&rgn42,&triangleinf,RGN_AND);
	rgn41.CombineRgn(&rgn41,&triangleinf,RGN_AND);
	rgn43.CombineRgn(&rgn43,&triangleinf,RGN_DIFF);

	rgn44.CombineRgn(&rgn44,&trianglesup,RGN_DIFF);
	rgn46.CombineRgn(&rgn46,&trianglesup,RGN_AND);
	rgn45.CombineRgn(&rgn45,&trianglesup,RGN_AND);
	rgn47.CombineRgn(&rgn47,&trianglesup,RGN_DIFF);

	trianglesup. DeleteObject( );
	triangleinf. DeleteObject( );
//////////////////************////////////////**************/////////////
	///////////	calcul dans les zones
//////////////////************////////////////**************/////////////
	//balayage rect
/*	for(int j=0;j<308;j++)
		for(int i=0;i<296;i++)
		{
			h=Inv_Get_ij(i,j);
			int l=i+j*296;
			if(rgn47.PtInRegion(h))	{Masquergn[l]=15;}
			else if(rgn46.PtInRegion(h))	{Masquergn[l]=14;}
			else if(rgn45.PtInRegion(h))	{Masquergn[l]=13;}
			else if(rgn44.PtInRegion(h))	{Masquergn[l]=12;}
			else if(rgn43.PtInRegion(h))	{Masquergn[l]=11;}
			else if(rgn42.PtInRegion(h))	{Masquergn[l]=10;}
			else if(rgn41.PtInRegion(h))	{Masquergn[l]=9;}
			else if(rgn40.PtInRegion(h))	{Masquergn[l]=8;}
			else if(rgn33.PtInRegion(h))	{Masquergn[l]=7;}
			else if(rgn32.PtInRegion(h))	{Masquergn[l]=6;}
			else if(rgn31.PtInRegion(h))	{Masquergn[l]=5;}
			else if(rgn30.PtInRegion(h))	{Masquergn[l]=4;}
			else if(rgn21.PtInRegion(h))	{Masquergn[l]=3;}
			else if(rgn20.PtInRegion(h))	{Masquergn[l]=2;}
			else if(ellipse1.PtInRegion(h))	{Masquergn[l]=1;}
			else Masquergn[l]=0;
		}*/
debut_i=debut_j=0;
fin_i=Larg-1;
fin_j=Haut-1;
//////////////////////////////////////////////////////////////////////////////////////
/*	BOOL prem;
	CPoint h;
	debut_i=debut_j=0;
	fin_i=293;
	fin_j=305;
	AffectRegionMasque();
	int cpt=0,i1,j1;
	for(int j=0;j<Haut;j++)
	{
		prem=TRUE;
		for(int i=0;i<Larg;i++)
		{		
			h=Inv_Get_ij(i,j);
			if(prem && RegionEllipse.PtInRegion(h))
			{
				i1=i;
				prem=FALSE;
			}
			if(!prem &&!RegionEllipse.PtInRegion(h))
			{
				if((i-i1)>13) 
				{
					if(!cpt) JL=j;
					DEBUTL[cpt]=i1;
					FINL[cpt++]=i-1;
				}
				prem=TRUE;
			}
		}
	}
	cpt=0;
	for(int i=0;i<Larg;i++)
	{
		prem=TRUE;
		for(int j=0;j<Haut;j++)
		{		
			h=Inv_Get_ij(i,j);
			if(prem && RegionEllipse.PtInRegion(h))
			{
				j1=j;
				prem=FALSE;
			}
			if(!prem &&!RegionEllipse.PtInRegion(h))
			{
				if((j-j1)>13)
				{
					if(!cpt) IH=i;
					DEBUTH[cpt]=j1;
					FINH[cpt++]=j-1;		
				}
				prem=TRUE;
			}
		}
	}*/
debut_i=debut_j=0;
fin_i=Larg-1;
fin_j=Haut-1;


	//////////masque//////////////////
	Masque=FALSE;
	Decoupage=FALSE;
	////////tabulation de abs//////////
	int l,c,vt;
	unsigned long pt=0;
	unsigned char contras;
	for(pt=0;pt<=0xFFFF;pt++)
	{	l=pt/0x100;c=pt&0xFF;
		if (l>c) 
		{
			vt=l-c;
			contras=(unsigned char)((float)255*(l-c)/(l+c));
		}
		else
		{
			vt=c-l;
			if(l==c) contras=0;
			else contras=(unsigned char)((float)255*(c-l)/(l+c));
		}
		TABABSO[pt]=unsigned char(vt/2);
		CONTRAST[pt]=contras;
	}



	/////////masque fixe////////////////////
	//*********************
/*	CPoint	h1 = Inv_Get_ij((int)10,(int)12); 
	CPoint	h2 = Inv_Get_ij((int)271-13,(int)271-13);
	CRgn	grandellipse;
			grandellipse.CreateEllipticRgnIndirect(CRect(h1,h2));
			h1 = Inv_Get_ij((int)30,(int)32); 
			h2 = Inv_Get_ij((int)271-33,(int)271-33);
	CRgn	petitellipse;
			petitellipse.CreateEllipticRgnIndirect(CRect(h1,h2));
			
	int j;
	CPoint h;
	for(j=0;j<272;j++)
		for(int i=0;i<272;i++)
		{
			h=Inv_Get_ij(i,j);
			if(petitellipse.PtInRegion(h))
			{
				Masquergn[i+272*j]=1;
			}
			else if(grandellipse.PtInRegion(h))
			{
				Masquergn[i+272*j]=2;
			}
			else Masquergn[i+272*j]=0;
	}
*/	
	//////////MASQUE POUR TEXTURE 
/*	CPoint	h1 = CPoint((int)23,(int)23); 
	CPoint	h2 = CPoint((int)235,(int)235);
	CRgn	grandellipse;
			grandellipse.CreateEllipticRgnIndirect(CRect(h1,h2));			
	
//	CPoint h;
	for(int j=0;j<272;j++)
		for(int i=0;i<272;i++)
		{
			h=CPoint(i,j);
			if(grandellipse.PtInRegion(h))
			{
				Masquergn[i+272*j]=1;
			}
			else Masquergn[i+272*j]=0;
	}
*/
///////vue3d////////////////////////////
//	PVue3d=new CVue3d(Image, Larg, Haut,(CWnd*)this);
//	VUE3D=FALSE;


	//////////////////////////////////////////////////
	///////polytopes multi-classes////////////////////
	//////////////////////////////////////////////////
	NomFichProj=NomFichAppr=NomFichParam=NomImageAppr="";
	
	NomFichMMC=NomFichCMC=NomFichPMC="";

	ImageAppr=NULL;
	DefImagAppr=FALSE;
	
	ValideMenu=FALSE;
	Visu=FALSE; 

	PParametre=NULL;
	ParamCalcul=FALSE;

	NbParam=0;
	NBPMAX=0;

	AffectNonClas=FALSE;
	ClassAffect=0;

	ValideDecision=FALSE;

	ValideVisuMes=FALSE;

	Visumes=FALSE;

	AffMesPol=NULL;

	CoeffPoly=7;

	for (int i=0;i<256;i++)  VisuPar[i]=NULL;
	
	Legende=new CDialLegend(this);

	SetCursor(LoadCursor(NULL,IDC_ARROW));

}

////////////////////////////////////////////////////////////
/////////////Traitements//////////////////////////////////
//////////////////////////////////////////////////////////
void CTraitImageWnd::OnTraitementLuminancemoyenneseuil() 
{
	long l,i,j,k,m;
	long temp;
	BOOL flag;
	
	ArreteLeContinu();
	Memorise();
	CDialSeuil dial (&Seuil,&flag); 
	if (dial.DoModal()==IDOK)
	{
		SetCursor(LoadCursor(0,IDC_WAIT));
	
	MessageBox("Calcul luminance moyenne fenêtre fixe 10*10+ seuil ");
	for( j=debut_j;j<=fin_j-10;j+=10)
		for( i=debut_i;i<=fin_i-10;i+=10)
			{   
				l=j*Larg+i;
				temp=0;
			 	for( k=0;k<10;k++)
					for( m=0;m<10;m++)
					{
						temp+=Image[l+k+m*Larg];
					}
				
				if(temp<100*Seuil) temp=0;
				else temp=255;
			for( k=0;k<10;k++)
					for( m=0;m<10;m++)
					{
						Image[l+k+m*Larg]=(unsigned char)temp;
					}
			}
	Invalidate(FALSE);	
	SetCursor(LoadCursor(0,IDC_ARROW));
	}//if dial.	
}

void CTraitImageWnd::OnTraitementLuminancemoyenne() 
{
	long l;
	long temp;
	
	
	ArreteLeContinu();
	Memorise();

	SetCursor(LoadCursor(0,IDC_WAIT));
	
//	MessageBox("Calcul luminance moyenne fenêtre glissante 12*12 ");
	for( long j=debut_j;j<=fin_j-12;j++)
		for(long i=debut_i;i<=fin_i-12;i++)
			{   
				l=j*Larg+i;
				temp=0;
			 	for(int k=0;k<12;k++)
					for(int m=0;m<12;m++)
					{
						temp+=Image[l+k+m*Larg];
					}
							
				Image[l]=(unsigned char)(temp/144);
			}
	Invalidate(FALSE);	
	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::OnTraitementGradientderob() 
{
	long l;
	ArreteLeContinu();
	Memorise();
	for(long j=debut_j;j<fin_j;j++)
		{for(long i=debut_i;i<fin_i;i++)
			{   
				l=j*Larg+i;
				Image[l]=(abs(Image[l]-Image[l+1+Larg])/2+abs(Image[l+Larg]-Image[l+1])/2); 
			}
		 
		}
	Invalidate(FALSE);
}

void CTraitImageWnd::OnTraitementSobel() 
{
	unsigned char * temp;
	temp=new unsigned char[Larg*Haut];
	ArreteLeContinu();
	Memorise();
	long i;
	for(i=0;i<Larg*Haut;i++) temp[i]=Image[i];
	for(i=Larg;i<(Larg-1)*Haut;i++)
	{
	Image[i]= (unsigned char)(labs(2*(long)temp[i-1]+(long)temp[i+Larg-1]+(long)temp[i-Larg-1]-2*(long)temp[i+1]-(long)temp[i+Larg+1]-(long)temp[i-Larg+1])+
  					 labs((long)temp[i+Larg-1]+2*(long)temp[i+Larg]+(long)temp[i+Larg+1]-(long)temp[i-Larg-1]-2*(long)temp[i-Larg]-(long)temp[i-Larg+1]));
	}
	Invalidate(FALSE);	
	delete temp;
}
///////////////filtres S et T///////////////////////////////
void CTraitImageWnd::OnTraitementFiltresett() 
{
	long valS,valT,temp,temp2;
	valS=16;
	valT=16;

	CDialSetT dlg(this);
	ArreteLeContinu();
	Memorise();
	
	if(dlg.DoModal()==IDOK)
	{
	  valS=dlg.m_S;
	  valT=dlg.m_T;
	} 
	else return;

	SetCursor(LoadCursor(0,IDC_WAIT));  

	for(int j=0;j<Haut;j++)
	  for(int i=valT;i<Larg-valS;i++)
	  {
		temp=0;
		temp2=1;
		for (int q=0;q<valS;q++) 
		{
			temp+=labs(Image[i-valT+q+Larg*j]-Image[i+q+j*Larg]);
//			temp2+=Image[i-valT/2+(int) pow(-1.0,q)*(1.0+q)/2.0+j*Larg];
			// supprimé car devenu inutile lol
		}

		//Image[i-valT+j*Larg]=(unsigned char) (temp2*temp/(255*valS*valS));
		Image[i-valT+j*Larg]=(unsigned char) (255*(double)temp/temp2);
	  }
	
	Invalidate(FALSE);
	SetCursor(LoadCursor(0,IDC_ARROW));	
}

void CTraitImageWnd::OnTraitementFiltressettFiltresettnonnormalis() 
{
	long valS,valT,temp;
	valS=16;
	valT=16;
	CDialSetT dlg(this);
	ArreteLeContinu();
	Memorise();

	SetCursor(LoadCursor(0,IDC_WAIT));  

	if(dlg.DoModal()==IDOK)
	{
	  valS=dlg.m_S;
	  valT=dlg.m_T;
	} 
	else return;



	for(int j=0;j<Haut;j++)
	for(int i=valT;i<Larg-valS;i++)
	{
	temp=0;
	for (int q=0;q<valS;q++) 
	{
		temp+=labs(Image[i-valT+q+Larg*j]-Image[i+q+j*Larg]);
	}

	Image[i-valT+j*Larg]=(unsigned char)(temp/valS);
	}
	Invalidate(FALSE);
	SetCursor(LoadCursor(0,IDC_ARROW)); 
}

void CTraitImageWnd::OnTraitementFiltressettHautbas() 
{
	long valS,valT,temp;
	unsigned char *tampon;
	valS=16;
	valT=16;
	tampon=new unsigned char [Larg*Haut];
	CDialSetT dlg(this);
	ArreteLeContinu();
	Memorise();
	

	if(dlg.DoModal()==IDOK)
	{
	  valS=dlg.m_S;
	  valT=dlg.m_T;
	} 
	else return;
	SetCursor(LoadCursor(0,IDC_WAIT));


	for(int j=0;j<Haut;j++)
	for(int i=valT;i<Larg-valS;i++)
	{
	temp=0;
	for (int q=0;q<valS;q++) 
	{
		temp+=labs(Image[i-valT+q+Larg*j]-Image[i+q+j*Larg]);
	}

	tampon[i+j*Larg]=(unsigned char)(temp/valS);
	}

	for(int i=0;i<Larg;i++)
	for(int j=valT;j<Haut-valS;j++)
	{
	temp=0;
	for (int q=0;q<valS;q++) 
	{
		temp+=labs(Image[i+Larg*(j+q-valT)]-Image[i+(j+q)*Larg]);
	}

	Image[i+(j-valT)*Larg]=(unsigned char)((temp/valS)+tampon[i+j*Larg])/2;
	}


	Invalidate(FALSE);	
	SetCursor(LoadCursor(0,IDC_ARROW));
	delete tampon;
}
//////////////////correlation/////////////////////

void CTraitImageWnd::OnTraitementCorrelation() 
{
	unsigned char *temp;
	int large,haute;
	long cor,moyim,moymas=0;
	large=fin_i-debut_i+1;
	haute=fin_j-debut_j+1;

	temp=new unsigned char [large*haute];
	ArreteLeContinu();
	Memorise();
	
	SetCursor(LoadCursor(0,IDC_WAIT));
	long i,j;
	for( i=0;i<large;i++)
		for( j=0;j<haute;j++)
	{
		temp[i+j*large]=Image[i+debut_i+(j+debut_j)*Larg];
		moymas+=Image[i+debut_i+(j+debut_j)*Larg];
	}


	moymas/=large*haute;
	long q,p;
	for(i=0;i<=Larg-large;i++)
		for(j=0;j<=Haut-haute;j++)
	{
		cor=0;
		moyim=0;
		for(q=0;q<haute;q++)
			for(p=0;p<large;p++)
		{
			moyim=Image[i+p+(j+q)*Larg];
		}
		moyim/=large*haute;


		for( q=0;q<haute;q++)
			for(p=0;p<large;p++)
		{
			if(((temp[p+q*large]<moymas)&&(Image[i+p+(j+q)*Larg]<moyim))||
				((temp[p+q*large]>moymas)&&(Image[i+p+(j+q)*Larg]>moyim)))
				cor+=temp[p+q*large]*Image[i+p+(j+q)*Larg];
			else
				cor-=temp[p+q*large]*Image[i+p+(j+q)*Larg];
		}
		if (cor<0) Image[i+j*Larg]=0;
		else Image[i+j*Larg]=cor/(64*large*haute);

	}
	Invalidate(FALSE);
	SetCursor(LoadCursor(0,IDC_ARROW));	
	delete temp;
}

void CTraitImageWnd::OnTraitementCorrelationCorrelation2() 
{
	float *temp,cor;
	int large,haute;
	long moymas=0;
	large=fin_i-debut_i+1;
	haute=fin_j-debut_j+1;

	if(large==0||haute==0) return;
	CTraitImageWnd *resultat=new CTraitImageWnd(Larg-large+1,Haut-haute+1,FLOAT32);	
	resultat->Create(Titre+": corrélation",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());

	temp=new float [large*haute];
	ArreteLeContinu();
	Memorise();
	
	SetCursor(LoadCursor(0,IDC_WAIT));


	long i,j;

	for( i=0;i<large;i++)
		for( j=0;j<haute;j++)
	{
		temp[i+j*large]=Image[i+debut_i+(j+debut_j)*Larg];
		moymas+=Image[i+debut_i+(j+debut_j)*Larg];
	}
	moymas/=large*haute;
	for( i=0;i<large;i++)
		for(j=0;j<haute;j++)
	{
		temp[i+j*large]-=moymas;
	}
moymas=0;
	for( i=0;i<Larg;i++)
		for( j=0;j<Haut;j++)
	{
		moymas+=Image[i+j*Larg];
	}

moymas/=Larg*Haut;
	for( i=0;i<resultat->Larg;i++)
		for( j=0;j<resultat->Haut;j++)
	{
		cor=0;
/*	moyim=0;
		for(long q=0;q<haute;q++)
			for(long p=0;p<large;p++)
		{
			moyim=Image[i+p+(j+q)*Larg];
		}
		moyim/=large*haute;
*/
		for(long q=0;q<haute;q++)
			for(long p=0;p<large;p++)
		{
			//if(((temp[p+q*large]<moymas)&&(Image[i+p+(j+q)*Larg]<moyim))||
			//	((temp[p+q*large]>moymas)&&(Image[i+p+(j+q)*Larg]>moyim)))
				cor+=(temp[p+q*large]*(Image[i+p+(j+q)*Larg]-moymas));
		//	else
			//	cor-=temp[p+q*large]*Image[i+p+(j+q)*Larg];
		}
		//if (cor<0) Image[i+j*Larg]=0;
	//	else Image[i+j*Larg]=(unsigned char)(255.0*((double)cor/(moycar)));
		//	cor/=large*haute;
			resultat->ImageFloat[i+j*resultat->Larg]=cor;

	}
	resultat->FloatToImage();
	Invalidate(FALSE);
	SetCursor(LoadCursor(0,IDC_ARROW));	
	delete temp;	
}
/////////////////////dst/////////////////////
void CTraitImageWnd::OnDst() 
{
	double rac2;
	double *tampon[9];
	long l,p,q;
	ArreteLeContinu();
	Memorise();
	
	SetCursor(LoadCursor(0,IDC_WAIT));
	rac2=sqrt(2.0);
	double temp[9][9]={	0.25,rac2/4,0.25,rac2/4,0.5,rac2/4,0.25,rac2/4,0.25,
			rac2/4,0,-rac2/4,0.5,-rac2/2,0.5,rac2/4,-0.5,rac2/4,
			0.25,-rac2/4,0.25,rac2/4,-0.5,rac2/4,0.25,-rac2/4,0.25,
			rac2/4,-0.5,rac2/4,0,0,0,rac2/4,0.5,rac2/4,
			0.5,0,-0.5,0,0,0,0.5,0,-0.5,
			0.5,-rac2/4,0.5,0,0,0,-0.5,rac2/4,-0.5,
			0.25,rac2/4,0.25,-rac2/4,-0.5,-rac2/4,0.25,rac2/4,0.25,
			rac2/4,0,-rac2/4,-0.5,0,0.5,rac2/4,0,-rac2/4,
			0.25,-rac2/4,0.25,-rac2/4,0.5,-rac2/4,0.25,-rac2/4,0.25};

	for(long i=0;i<9;i++) tampon[i]=new double [Larg*Haut];


	for(long j=0;j<Haut-2;j++)
	{
		
		for(long i=0;i<Larg-2;i++)
		{
			l=i+j*Larg;
			for(int t=0;t<9;t++) tampon[t][l]=0;
			for(p=0;p<9;p++)
			{
				
				for(q=0;q<9;q++)
				{
					tampon[p][l]+=temp[p][q]*Image[i+(q%3)+(j+2-q/3)*Larg];
				}
				if(p==0) tampon[p][l]/=(6+4*rac2)/3;
				//if(p==1) tampon[p][l]*=2;
				//if(p==2) tampon[p][l]*=2;
				if(p==3) tampon[p][l]/=rac2;
				if(p==4) tampon[p][l]*=2;
				if(p==5) tampon[p][l]/=rac2/4;
			//if(p==6) tampon[p][l]/=0.5;
				if(p==7) tampon[p][l]/=rac2/4;
				if(p==8) tampon[p][l]/=(6-4*rac2);
				tampon[p][l]=fabs(tampon[p][l]);
			}
			
		}
	}
	int i;
	CImageWnd *un,*deux,*trois,*quatre,*cinq,*six,*sept,*huit,*neuf;
	un=new CImageWnd(Larg,Haut);
	deux=new CImageWnd(Larg,Haut);
	trois=new CImageWnd(Larg,Haut);
	quatre=new CImageWnd(Larg,Haut);
	cinq=new CImageWnd(Larg,Haut);
	six=new CImageWnd(Larg,Haut);
	sept=new CImageWnd(Larg,Haut);
	huit=new CImageWnd(Larg,Haut);
	neuf=new CImageWnd(Larg,Haut);
	for(i=0;i<Haut*Larg;i++)
	{
		un->Image[i]=(unsigned char)tampon[0][i];
		deux->Image[i]=(unsigned char)tampon[1][i];
		trois->Image[i]=(unsigned char)tampon[2][i];
		quatre->Image[i]=(unsigned char)tampon[3][i];
		cinq->Image[i]=(unsigned char)tampon[4][i];
		six->Image[i]=(unsigned char)tampon[5][i];
		sept->Image[i]=(unsigned char)tampon[6][i];
		huit->Image[i]=(unsigned char)tampon[7][i];
		neuf->Image[i]=(unsigned char)tampon[8][i];
	}

	un->Create("premier voisin", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,NULL);
	deux->Create("deuxième voisin", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,NULL);
	trois->Create("troisième voisin", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,NULL);
	quatre->Create("quatrième voisin", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,NULL);
	cinq->Create("cinquième voisin",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,NULL);
	six->Create("sixième voisin", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,NULL);
	sept->Create("septième voisin", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,NULL);
	huit->Create("huitième voisin",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,NULL);
	neuf->Create("neuvième voisin",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,NULL);
	un->ShowWindow(SW_SHOW);
	deux->ShowWindow(SW_SHOW);
	trois->ShowWindow(SW_SHOW);
	quatre->ShowWindow(SW_SHOW);
	cinq->ShowWindow(SW_SHOW);
	six->ShowWindow(SW_SHOW);
	sept->ShowWindow(SW_SHOW);
	huit->ShowWindow(SW_SHOW);
	neuf->ShowWindow(SW_SHOW);
	Invalidate(FALSE);
	SetCursor(LoadCursor(0,IDC_ARROW));	
	for( i=0;i<9;i++) delete tampon[i];

}

//*****************************************************
////////////contours//////////////////////////////////
//****************************************************

void CTraitImageWnd::OnContoursContourhystrsis1() 
{
	long i,l,p,k=0,temp,courant,dir,cont;
	long * sommet;
	unsigned char * contour;
	unsigned char *ImageSommet;
	long * compteur;
	long memocontour[20000];
	char ch[90];
	BOOL check;

	sommet=new long [Larg*Haut];

	compteur=new long[Larg*Haut/4 +1];


	CTraitImageWnd * pImageSom =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,NB);   
	ImageSommet=pImageSom->Image;
	
	CTraitImageWnd * pImageContour =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,NB);   
	contour=pImageContour->Image;

	CTraitImageWnd * pImageContoursupseuil =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,NB);   

	compteur[0]=0;

	ArreteLeContinu();
	Memorise();
	
	
	////boite de dialogue des seuils
	CDialContour Dilg(Seuilcont1,Seuilcont2,LongmaxContour,SeuillongContour,this);
	if (Dilg.DoModal()==IDOK)
	{
		Seuilcont1=(unsigned char) Dilg.m_seuilhaut;
		Seuilcont2=(unsigned char)Dilg.m_seuilbas;
		LongmaxContour=Dilg.m_longmax;
		SeuillongContour=Dilg.m_seuilcont;
		check=Dilg.m_check;
	
	}
	else return;

	SetCursor(LoadCursor(0,IDC_WAIT));

	////tour de l'image à 0
	for(i=0;i<Larg;i++) 
	{
		Image[i]=0;
		Image[(Haut-1)*Larg+i]=0;
	}
	for(i=0;i<Haut;i++)
	{
		Image[i*Larg]=0;
		Image[(i+1)*Larg-1]=0;
	}
	/////initialisation images

	for ( i=0; i< Larg*Haut; i++) 
	{
		contour[i]=0;
		ImageSommet[i]=0;
	}
	////recherche sommets
	l=0;
	for (i=0;i< Haut*Larg;i++)
		{
			if(Image[i]>=Seuilcont1)
				{
					sommet[l]=i;
					ImageSommet[i]=255;
					l++;
				}
		}

	
	///recherche contours
	for (i=0;i<l;i++)
		{
			courant=sommet[i];
			
			while(Image[courant-1]>=Seuilcont2) 
				{courant--;};
			temp=courant;

			if (contour[courant]!=255)
			{
				contour[courant]=255;
				compteur[0]++;
				compteur[compteur[0]]=1;
				dir=0;
				cont=1;
				p=0;
				memocontour[p++]=courant;////pour le seuillage sur la longueur
				

				do
				{
					
					switch(dir+cont)
					{
						case 0:
							courant-=Larg;
							dir=6;
							break;
						case 1:
							courant+=Larg;
							dir=2;
							break;
						case 2:
							courant++;
							dir=0;
							break;
						case 3:
							courant--;
							dir=4;
							break;
						case 4:
							courant+=Larg;
							dir=2;
							break;
						case 5:
							courant-=Larg;
							dir=6;
							break;
						case 6:
							courant--;
							dir=4;
							break;
						case 7:
							courant++;
							dir=0;
							break;
					}
					
					if (Image[courant]>=Seuilcont2)
						{
							cont=1;
							if(contour[courant]!=255)
							{
								compteur[compteur[0]]++;
								contour[courant]=255;
								memocontour[p++]=courant;
								if(compteur[compteur[0]]>=LongmaxContour){MessageBox("contour trop long!");return;}
								
							}
						}
					else cont=0;
				}
				while(courant!=temp);

				if(!check) 
				{
					if(p>=SeuillongContour)
					{	k++;
						for(long j=0;j<p;j++) pImageContoursupseuil->Image[memocontour[j]]=255;
					}
				}
				else
				{
					if(p<SeuillongContour)
					{	k++;
						for(long j=0;j<p;j++) pImageContoursupseuil->Image[memocontour[j]]=255;
					}
				}
			}////if image=sommet
			
		}//for l
	pImageSom->Create(Titre+" : Sommets",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
	pImageContour->Create(Titre+" : contour",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
	
	sprintf(ch,"%d",SeuillongContour);
	pImageContoursupseuil->Create(Titre+" : contours supérieurs à "+ch,WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
	GetMDIFrame()->MDITile(MDITILE_VERTICAL);
	this->SetFocus();
	SetCursor(LoadCursor(0,IDC_ARROW)); 
	sprintf(ch,"nombre de sommets:%d\nnbre total de contours: %d\nnbre de contour de long>=seuil: %d",l,compteur[0],k);
	MessageBox(ch,"contours");
	delete sommet;
	delete compteur;
}

void CTraitImageWnd::OnContoursContourhystrsis2() 
{
	long i,l,p,k=0,courant,suivant,dir,dirinit,p0,p1;
	BOOL fini;
	long * sommet;
	long memocontour[4000];
	unsigned char * contour;
	unsigned char *ImageSommet;
	long * compteur;
	char ch[90];

	sommet=new long [Larg*Haut];
	compteur=new long[Larg*Haut/4 +1];
	CTraitImageWnd * pImageSom =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,NB);   
	ImageSommet=pImageSom->Image;
	CTraitImageWnd * pImageContour =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,NB);   
	contour=pImageContour->Image;
	CTraitImageWnd * pImageContoursupseuil =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,NB);   
	
	compteur[0]=0;

	ArreteLeContinu();
	Memorise();

	////boite de dialogue des seuils
	CDialContour Dilg(Seuilcont1,Seuilcont2,LongmaxContour,SeuillongContour,this);
	if (Dilg.DoModal()==IDOK)
	{
		Seuilcont1=(unsigned char)Dilg.m_seuilhaut;
		Seuilcont2=(unsigned char)Dilg.m_seuilbas;
		LongmaxContour=Dilg.m_longmax;
		SeuillongContour=Dilg.m_seuilcont;
	
	}
	else return;

	SetCursor(LoadCursor(0,IDC_WAIT));
	

	////tour de l'image à 0
	for(i=0;i<Larg;i++) 
	{
		Image[i]=0;
		Image[(Haut-1)*Larg+i]=0;
	}
	for(i=0;i<Haut;i++)
	{
		Image[i*Larg]=0;
		Image[(i+1)*Larg-1]=0;
	}

	//////initialisation des images sommet et contour et contour supseuil
	
	for ( i=0; i< Larg*Haut; i++) 
	{
		contour[i]=0;
		ImageSommet[i]=0;
		pImageContoursupseuil->Image[i]=0;
	}
	
	//////recherche des sommets avec affichage dans imagesommet et mémorisation dans sommet
	l=0;
	for (i=0;i< Haut*Larg;i++)
		{
			if(Image[i]>= Seuilcont1)
				{
					sommet[l]=i;
					ImageSommet[i]=255;
					l++;
				}
		}
		
	////détections des contours des sommets et comptage
	p=0;
	for(i=0;i<l;i++)
	{
		courant=sommet[i];///courant sur le sommet n°i
	
		//////détection d'un point de contour à partir du sommet
		while(Image[courant-1]>=Seuilcont2) 
			{courant--;};
		
		/////si le point de contour n'a pas déjà été trouvé rechercher le contour 
		if (contour[courant]!=255)
		{
			fini=FALSE;
			dir=4;
			p0=courant;
			contour[courant]=255;
			compteur[0]++;
			compteur[compteur[0]]=1;
			p=0;
			memocontour[p++]=courant;////pour le seuillage sur la longueur

			if(!((Image[p0-1]<Seuilcont2)&&(Image[p0+1]<Seuilcont2)&&(Image[p0+Larg]<Seuilcont2)&&
				(Image[p0+Larg-1]<Seuilcont2)&&(Image[p0+Larg+1]<Seuilcont2)&&
				(Image[p0-Larg]<Seuilcont2)&&(Image[p0-Larg+1]<Seuilcont2)&&
				(Image[p0-Larg-1]<Seuilcont2)))
			{
				while(!fini)
				{
					do
					{
						dir=(dir+1)%8;

						switch (dir)
						{
							case 0:
								suivant=courant+1;
								break;
							case 1:
								suivant=courant+Larg+1;
								break;
							case 2:
								suivant=courant+Larg;
								break;
							case 3:
								suivant=courant+Larg-1;
								break;
							case 4:
								suivant=courant-1;
								break;
							case 5:
								suivant=courant-Larg-1;
								break;
							case 6:
								suivant=courant-Larg;
								break;
							case 7:
								suivant=courant-Larg+1;
								break;
						}
						dirinit=(dir+5)%8;
					}
					while(Image[suivant]<Seuilcont2);

					if(compteur[compteur[0]]==1) p1=suivant;
					else
					{
						if((courant==p0)&&(suivant==p1)) fini=TRUE;
						if(compteur[compteur[0]]>=LongmaxContour){MessageBox("contour trop long!");return;}
					}
					courant=suivant;
					if(contour [courant]!=255)
					{
						contour[courant]=255;
						compteur[compteur[0]]++;
						memocontour[p++]=courant;
					}
					dir=dirinit;		
				}//while fini
			}///if !....

			/////pour le seuillage sur la longueur
			if(p>=SeuillongContour) 
			{	
				k++;
				for(long j=0;j<p;j++) pImageContoursupseuil->Image[memocontour[j]]=255;
			}
		}//if contour
	}//for i=0->l

	pImageSom->Create(Titre+" : Sommets",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
	pImageContour->Create(Titre+" : contours",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
	sprintf(ch,"%d",SeuillongContour);
	pImageContoursupseuil->Create(Titre+" : contours supérieurs à "+ch,WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
	GetMDIFrame()->MDITile(MDITILE_VERTICAL);
	this->SetFocus();
	SetCursor(LoadCursor(0,IDC_ARROW)); 

	sprintf(ch,"nombre de sommets:%d\nnbre total de contours: %d\nnbre de contour de long>=seuil: %d",l,compteur[0],k);
	MessageBox(ch,"contours");
	delete sommet;
	delete compteur;
}

void CTraitImageWnd::OnContoursContourhystrsisfen88() 
{
	long i,j,li,co,l,k=0,courant,suivant,dir,dirinit,p0,p1;
	BOOL fini;
	long sommet[64][2];
	int  compteur,max;
	int contour[100];
	unsigned char imageTampon[100];


	compteur=0;

	ArreteLeContinu();
	Memorise();
	
	////boite de dialogue des seuils
	CDialContour Dilg(Seuilcont1,Seuilcont2,LongmaxContour,SeuillongContour,this);
	if (Dilg.DoModal()==IDOK)
	{
		Seuilcont1=(unsigned char)Dilg.m_seuilhaut;
		Seuilcont2=(unsigned char)Dilg.m_seuilbas;
		LongmaxContour=Dilg.m_longmax;
		SeuillongContour=Dilg.m_seuilcont;
	
	}
	else return;

	SetCursor(LoadCursor(0,IDC_WAIT));
	
	
	for(i=0;i<Larg-8;i+=8)
		for(j=0;j<Haut-8;j+=8)
		{
			///recherche des sommets dans fen 8*8 et intialisation du buffer contour et image tampon
			
			
			for(k=0;k<100;k++) 
			{
				imageTampon[k]=0;
				contour[k]=0;
			}

			l=0;
			for(co=0;co<8;co++)
				for(li=0;li<8;li++)
				{
					imageTampon[1+co+10*(li+1)]=Image[i+co+(j+li)*Larg];
					if(Image[i+co+(j+li)*Larg]>=Seuilcont1) {sommet[l][0]=co;sommet[l][1]=li;l++;}
				}//for co for li
			
			///recherche des contours
			compteur=0;
			max=0;
			for(k=0;k<l;k++)
			{
				courant=1+sommet[k][0]+10*(1+sommet[k][1]);

				while(imageTampon[courant-1]>=Seuilcont2)
					{courant--;};
				
				if(contour[courant]==0)
				{
					fini=FALSE;
					dir=4;
					p0=courant;
					contour[courant]=1;
					compteur=1;
					
					if(!((imageTampon[p0-1]<Seuilcont2)&&(imageTampon[p0+1]<Seuilcont2)&&(imageTampon[p0+10]<Seuilcont2)&&
						(imageTampon[p0+10-1]<Seuilcont2)&&(imageTampon[p0+10+1]<Seuilcont2)&&
						(imageTampon[p0-10]<Seuilcont2)&&(imageTampon[p0-10+1]<Seuilcont2)&&
						(imageTampon[p0-10-1]<Seuilcont2)))
					{
						while(!fini)
						{
							do
							{
								dir=(dir+1)%8;

								switch (dir)
								{
									case 0:
										suivant=courant+1;
										break;
									case 1:
										suivant=courant+10+1;
										break;
									case 2:
										suivant=courant+10;
										break;
									case 3:
										suivant=courant+10-1;
										break;
									case 4:
										suivant=courant-1;
										break;
									case 5:
										suivant=courant-10-1;
										break;
									case 6:
										suivant=courant-10;
										break;
									case 7:
										suivant=courant-10+1;
										break;
								}
								dirinit=(dir+5)%8;
							}
							while(imageTampon[suivant]<Seuilcont2);

							if(compteur==1) p1=suivant;
							else
							{
								if((courant==p0)&&(suivant==p1)) fini=TRUE;
								if(compteur>=LongmaxContour){MessageBox("contour trop long!");return;}
							}
							courant=suivant;
							if(contour [courant]==0)
							{
								contour[courant]=1;
								compteur++;
							}
							dir=dirinit;		
						}//while fini
					}///if !....

				}///if !contour
				if(compteur>max) max=compteur;
			}//for k

			for(co=0;co<8;co++)
				for(li=0;li<8;li++)
				{
					if(8*max<256) Image[i+co+(j+li)*Larg]=8*max;
					else Image[i+co+(j+li)*Larg]=255;
				}
			
		}//for i for j
Invalidate(FALSE);
SetCursor(LoadCursor(0,IDC_ARROW));  
	
}
/////////////////////////////////////////////////////////////
//////////////////////modification d'histogramme///////////
/////////////////////////////////////////////////////////

void CTraitImageWnd::OnTraitementReductionsur16niveaux() 
{
	ArreteLeContinu();
	Memorise();
	
	for(long i=0;i<Larg*Haut;i++)
	{
		Image[i]=(unsigned char) ldiv(Image[i],16).quot;
	}
	Invalidate(FALSE);
}

void CTraitImageWnd::OnTraitementEgalisationdhistogramme() 
{
	long histo[256];
	long i,j;
	ArreteLeContinu();
	Memorise();

	for( j=0;j<256;j++) histo[j]=0;
	for( j=debut_j;j<=fin_j;j++)
	for( i=debut_i;i<=fin_i;i++)
		{
			
			histo[Image[i+j*Larg]]++;
			
		}
	for(j=1;j<256;j++) histo[j]+=histo[j-1];
	for( j=debut_j;j<=fin_j;j++)
	for( i=debut_i;i<=fin_i;i++)
		{
						
			Image[i+j*Larg]=(unsigned char)((255*histo[Image[i+j*Larg]])/((fin_i-debut_i)*(fin_j-debut_j)));
			
		}
Invalidate(FALSE);
	
}

void CTraitImageWnd::OnTraitementModificationdhistogrammeExpansiondynamique() 
{
	int min=255,max=0;
	long l;
	double a,b;
	ArreteLeContinu();
	Memorise();
	long i,j;
	for( j=debut_j;j<=fin_j;j++)
	for( i=debut_i;i<=fin_i;i++)
		{
			l=i+j*Larg;
			if(Image[l]<min) min=Image[l];
			if(Image[l]>max) max=Image[l];
		}

	if(max==min) {MessageBox("Expansion dynamique impossible");return;}
	if((max==255)&&(min==0)) {MessageBox("Dynamique maximum!");return;}
	a=(double)(-255.0*min/(max-min));
	b=(double)(255.0/(max-min));
	for( j=debut_j;j<=fin_j;j++)
	for( i=debut_i;i<=fin_i;i++)
		{
			l=i+j*Larg;			
			Image[l]=(unsigned char) (b*Image[l]+a);
			
		}
Invalidate(FALSE);
	
}


//*****************************************************
////////////////divers///////////////////////
//*****************************************************

void CTraitImageWnd::OnEdition255() 
{
	ArreteLeContinu();
	Memorise();

	for(long i=0;i<Larg*Haut;i++) if(Image[i]==255) Image[i]=254;
	Invalidate(FALSE);	
}
/////////////////////////////////////////
//////////////////////Max-min//////////////////
/////////////////////////////////////////

void CTraitImageWnd::OnTraitementMaxmin() 
{
	char ch[20];
	int q,l,i,j;
	unsigned char classe[256];
	unsigned char Max,Maxc,courant,ancien;

	unsigned char *maximuml;
	unsigned char *maximumc;
	
	maximuml=new unsigned char[Larg*Haut];
	maximumc=new unsigned char[Larg*Haut];

	ArreteLeContinu();
	Memorise();

	for(int i=0;i<64;i++) classe[i]=classe[255-i]=classe[127-i]=classe[128+i]=0;
	//traitement ligne
	
	DWORD Time=GetCurrentTime();
	int t;
for( t=0;t<100;t++)
	{
		q=l=0;

		//////////////////////////////lignes
		for( j=0;j<Haut;j++)
		{
			Max=0;
			for( i=0;i<12;i++)
			{
				courant=Image[q];
				classe[courant]++;
				(((courant) > (Max)) ? (Max=courant) : (0));
				//Max=__max(courant,Max);
			//	maximuml[q++]=Max;
				q++;
				
			}//for i 0->11
			
			for(i=12;i<Larg;i++)
			{
				courant=Image[q];
				ancien=Image[l++];
				classe[courant]++;
				classe[ancien]--;
				if(courant>Max) Max=courant;
				else if(ancien==Max) 
				{
					while(!classe[Max])
					{
						Max--;
					}
				}
				maximuml[q++]=Max;
			}//for i 12->Larg
		for(i=0;i<12;i++) classe[Image[l++]]--;
		}///for j+


		///////////////////////////////colonnes
		for(i=0;i<Larg;i++)
		{
			q=l=i;Maxc=0;
			for( j=0;j<12;j++)
			{
				courant=maximuml[q];
				classe[courant]++;
				(((courant) > (Maxc)) ? (Maxc=courant) : (0));
			//	maximumc[q]=Maxc;
				q+=Larg;
			}

			for(j=12;j<Haut;j++)
			{
				courant=maximuml[q];
				ancien=maximuml[l];
				classe[courant]++;
				classe[ancien]--;
				if(courant>Maxc) Maxc=courant;
				else if(ancien==Maxc) 
				{
					while(!classe[Maxc])
					{
						Maxc--;
					}
				}
				maximumc[q]=Maxc;
				q+=Larg;
				l+=Larg;
			}
			for(j=0;j<12;j++) {classe[maximuml[l]]--;l+=Larg;}
		}//for i


	}/// for t


	Time=GetCurrentTime()-Time;

	sprintf(ch,"time=%d",Time);
	MessageBox(ch);

	for( i=0;i<Larg*Haut;i++) Image[i]=maximumc[i];

	delete	maximuml;
	delete maximumc;
	Invalidate(FALSE);	
}
/////////////////////////////////////////
////////soustraction
////////////////////////////////////////
void CTraitImageWnd::OnTraitementSoustraction() 
{
	int  num1=1,num2=35;
	CString s,NomFichier;
	CFile Fichier;
	long large,haute;
	DWORD DWComp;
	WORD BitCount;
	BITMAPFILEHEADER Tete;
	CString NomFichier2;
	unsigned char *pBuffer;

	if(SEQUENCE_CREE) 
		if(MessageBox("Différence de deux séquences?","Soustraction",MB_YESNO)==IDYES)
		{
			CFileDialog * FileData = new CFileDialog
			(TRUE,"*.bmp",NULL,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT,"Images Bitmap(*.bmp)|*.bmp||",this);

			if (FileData->DoModal()==IDOK)
			{
				if (FileData->GetFileExt() == "bmp")
				{
					if (Fichier.Open(FileData->GetPathName(), CFile::modeRead))
					{
						Fichier.Read(&Tete,sizeof(BITMAPFILEHEADER));
						Fichier.Seek(sizeof(DWORD),CFile::current);
						Fichier.Read(&large,sizeof(LONG)); 
						Fichier.Read(&haute,sizeof(LONG));
						Fichier.Seek(sizeof(WORD),CFile::current);
						Fichier.Read(&BitCount,sizeof(WORD));
						Fichier.Read(&DWComp,sizeof(DWORD));
						if ((DWComp != BI_RGB)||(BitCount!=8 && BitCount!=24))
						{
							MessageBox("Impossible de lire les fichiers BMP compressés","Erreur",MB_OK|MB_ICONSTOP);
							Fichier.Close();
							return;
						}
  	 					else if ((large!=Larg)||(haute!=Haut))
						{
							MessageBox("Taille incompatible!","Erreur",MB_OK|MB_ICONSTOP);
							Fichier.Close();
						}
						else
						{
							Fichier.Close();
							CDialSeqlimite dial(&num1,&num2); 
							if (dial.DoModal()!=IDOK) {delete FileData;return;}
							if(num2-num1+1>NbMaxImageSeq) { MessageBox("Différence impossible");delete FileData;return;}

							s = FileData->GetPathName().Left(FileData->GetPathName().GetLength()-7); 
							SetCursor(LoadCursor(0,IDC_WAIT));   
							pBuffer = (unsigned char *) malloc(Larg*Haut);
							for (int j=0;j<=num2-num1;j++)
							{   
								NomFichier.Format(s+"%.3d.bmp",j+num1);
								if (!Fichier.Open(NomFichier, CFile::modeRead)) {delete FileData; return;}

								Fichier.Seek(Tete.bfOffBits,CFile::begin);
							
								Fichier.Read(pBuffer,Larg*Haut);

								for(long i=0;i<Larg*Haut;i++)
								 Bpile[j][i] = abs(Bpile[j][i]-pBuffer[i]);
								Fichier.Close();
							}
							free(pBuffer);
							AfficheSequence(0,1);
						}
					
					}
				}
				else MessageBox("Différence d'images entre des fichiers BMP uniquement");
			}
			delete FileData;
			return;
		}


	CFileDialog * FileData = new CFileDialog
	(TRUE,"*.bmp",NULL,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT,"Images Bitmap(*.bmp)|*.bmp||",this);

	if (FileData->DoModal()==IDOK)
	{
		if (FileData->GetFileExt() == "bmp")
		{
			if (Fichier.Open(FileData->GetPathName(), CFile::modeRead))
			{
				Fichier.Read(&Tete,sizeof(BITMAPFILEHEADER));
				Fichier.Seek(sizeof(DWORD),CFile::current);
				Fichier.Read(&large,sizeof(LONG)); 
				Fichier.Read(&haute,sizeof(LONG));
				Fichier.Seek(sizeof(WORD),CFile::current);
				Fichier.Read(&BitCount,sizeof(WORD));
				Fichier.Read(&DWComp,sizeof(DWORD));
				if ((DWComp != BI_RGB)||(BitCount!=8 && BitCount!=24))
				{
					MessageBox("Impossible de lire les fichiers BMP compressés","Erreur",MB_OK|MB_ICONSTOP);
					Fichier.Close();
					return;
				}
  	 			else if ((large!=Larg)||(haute!=Haut))
					MessageBox("Taille incompatible!","Erreur",MB_OK|MB_ICONSTOP);
				else
				{
					Fichier.Seek(Tete.bfOffBits,CFile::begin);
					pBuffer = (unsigned char *) malloc(Larg*Haut);
					Fichier.Read(pBuffer,Larg*Haut);

					for(long i=0;i<Larg*Haut;i++)
					 Image[i] = abs(Image[i]-pBuffer[i]);

					Invalidate(FALSE);
					free(pBuffer);
				}
				Fichier.Close();
			}
		}
		else MessageBox("Différence d'images entre des fichiers BMP uniquement");
	}
	delete FileData;
}
/////////////////////////////////////////
//// 5 paramètres
/////////////////////////////////////////

void CTraitImageWnd::OnMasque5paramtres() 
{
	ArreteLeContinu();
//	SetCursor(LoadCursor(0,IDC_WAIT));


	Fermeture1212();

	GradientGlissant1212();

	Ouverture1212();

	LuminanceGlissante1212();

	ContrasteGlissant1212();

	this->SetFocus();
	GetMDIFrame()->MDITile(MDITILE_VERTICAL);	
//	SetCursor(LoadCursor(0,IDC_ARROW));
}


void CTraitImageWnd::Fermeture1212()
{
//	SetCursor(LoadCursor(0,IDC_WAIT));
	int stop,l,i;
	unsigned char *resultat;
	CTraitImageWnd *fermeture;

	resultat=new unsigned char[Larg*Haut];
	fermeture=new CTraitImageWnd(Larg,Haut,BMP);
//	DWORD time=GetCurrentTime();
//	for( i=0;i<100;i++)
	{	
		MaximumGlissant1212(fermeture->Image);
		fermeture->MinimumGlissant1212(resultat);
	}
/*	time=GetCurrentTime()-time;
	CString text;
	text.Format("temps: %d",time);
	MessageBox(text,"Fermeture");
*/	if((Larg<12)||(Haut<12)){MessageBox("image trop petite","Traitement impossible");return;}
	stop=12+Larg*12;
	l=stop;
	for( i=0;i<Larg*Haut-stop;i++) fermeture->Image[i]=resultat[l++];
	for(i=Larg*Haut-stop;i<Larg*Haut;i++) fermeture->Image[i]=0;
	fermeture->Create("fermeture glis 12*12 de "+Titre, WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());
	delete resultat;
//	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::Ouverture1212()
{
//	SetCursor(LoadCursor(0,IDC_WAIT));
	int stop,l,i;
	unsigned char *resultat;
	CTraitImageWnd *ouverture;

	resultat=new unsigned char[Larg*Haut];
	ouverture=new CTraitImageWnd(Larg,Haut,BMP);
//	DWORD time=GetCurrentTime();
//	for( i=0;i<100;i++)
	{
		MinimumGlissant1212(ouverture->Image);
		ouverture->MaximumGlissant1212(resultat);
	}
/*	time=GetCurrentTime()-time;
	CString text;
	text.Format("temps: %d",time);
	MessageBox(text,"Ouverture");
*/
	if((Larg<12)||(Haut<12)){MessageBox("image trop petite","Traitement impossible");return;}
	stop=12+Larg*12;
	l=stop;
	for( i=0;i<Larg*Haut-stop;i++) ouverture->Image[i]=resultat[l++];	
	for(i=Larg*Haut-stop;i<Larg*Haut;i++) ouverture->Image[i]=0;
	ouverture->Create("ouverture glis 12*12 de "+Titre, WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());
	delete resultat;
//	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::MaximumGlissant1212(unsigned char *maximumc)
{
//	SetCursor(LoadCursor(0,IDC_WAIT));
	int i,j,q,l;
	unsigned char classe[256];
	unsigned char Max,Maxc,courant,ancien;

	unsigned char *maximuml;

	
	maximuml=new unsigned char[Larg*Haut];

	for(int i=0;i<64;i++) classe[i]=classe[255-i]=classe[127-i]=classe[128+i]=0;
	//traitement ligne
	q=l=0;

	//////////////////////////////lignes
	for( j=0;j<Haut;j++)
	{
		Max=0;
		for( i=0;i<12;i++)
		{
			courant=Image[q];
			classe[courant]++;
			Max=(((courant) > (Max)) ? (courant) : (Max));
			//Max=__max(courant,Max);
			maximuml[q++]=Max;
			
		}//for i 0->11
		
		for(i=12;i<Larg;i++)
		{
			courant=Image[q];
			ancien=Image[l++];
			classe[courant]++;
			Max=(((courant) > (Max)) ? (courant) : (Max));
			classe[ancien]--;
			if(ancien==Max) 
			{
				while(!classe[Max])
				{
					Max--;
				}
			}
			maximuml[q++]=Max;
		}//for i 12->Larg
	for(i=0;i<12;i++) classe[Image[l++]]--;
	}///for j+


	///////////////////////////////colonnes
	for(i=0;i<Larg;i++)
	{
		q=l=i;Maxc=0;
		for( j=0;j<12;j++)
		{
			courant=maximuml[q];
			classe[courant]++;
			Maxc=(((courant) > (Maxc)) ? (courant) : (Maxc));
			maximumc[q]=Maxc;
			q+=Larg;
		}

		for(j=12;j<Haut;j++)
		{
			courant=maximuml[q];
			ancien=maximuml[l];
			classe[courant]++;
			Maxc=(((courant) > (Maxc)) ? (courant) : (Maxc));
			classe[ancien]--;
			if(ancien==Maxc) 
			{
				while(!classe[Maxc])
				{
					Maxc--;
				}
			}
			maximumc[q]=Maxc;
			q+=Larg;
			l+=Larg;
		}
		for(j=0;j<12;j++) {classe[maximuml[l]]--;l+=Larg;}
	}//for i
	delete maximuml;
//	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::MinimumGlissant1212(unsigned char *Minimumc)
{
//	SetCursor(LoadCursor(0,IDC_WAIT));
	int q,l;
	unsigned char classe[256];
	unsigned char Min,Minc,courant,ancien;

	unsigned char *Minimuml;
	int i,j;
	
	Minimuml=new unsigned char[Larg*Haut];

	for(int i=0;i<64;i++) classe[i]=classe[255-i]=classe[127-i]=classe[128+i]=0;
	//traitement ligne
	q=l=0;

	//////////////////////////////lignes
	for( j=0;j<Haut;j++)
	{
		Min=255;
		for( i=0;i<12;i++)
		{
			courant=Image[q];
			classe[courant]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			//Min=__Min(courant,Min);
			Minimuml[q++]=Min;
			
		}//for i 0->11
		
		for(i=12;i<Larg;i++)
		{
			courant=Image[q];
			ancien=Image[l++];
			classe[courant]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			classe[ancien]--;
			if(ancien==Min) 
			{
				while(!classe[Min])
				{
					Min++;
				}
			}
			Minimuml[q++]=Min;
		}//for i 12->Larg
	for(i=0;i<12;i++) classe[Image[l++]]--;
	}///for j+


	///////////////////////////////colonnes
	for(i=0;i<Larg;i++)
	{
		q=l=i;Minc=255;
		for( j=0;j<12;j++)
		{
			courant=Minimuml[q];
			classe[courant]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Minimumc[q]=Minc;
			q+=Larg;
		}

		for(j=12;j<Haut;j++)
		{
			courant=Minimuml[q];
			ancien=Minimuml[l];
			classe[courant]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			classe[ancien]--;
			if(ancien==Minc) 
			{
				while(!classe[Minc])
				{
					Minc++;
				}
			}
			Minimumc[q]=Minc;
			q+=Larg;
			l+=Larg;
		}
		for(j=0;j<12;j++) {classe[Minimuml[l]]--;l+=Larg;}
	}//for i
	delete Minimuml;
//	SetCursor(LoadCursor(0,IDC_ARROW));
}


void CTraitImageWnd::GradientGlissant1212()
{
//	SetCursor(LoadCursor(0,IDC_WAIT));

	int i;
	unsigned char *tamp;
	unsigned __int16 *templ,*tempc;
	int large=Larg-Larg%4;
	int haute=Haut-Haut%4;
	CTraitImageWnd *gradMoy;
	long stop,l;
	templ=new unsigned __int16[large*haute];
	tempc=new unsigned __int16[large*haute];
	tamp=new unsigned char[large*haute];
	gradMoy=new CTraitImageWnd(large,haute,BMP);
//	DWORD time=GetCurrentTime();
//	for(i=0;i<100;i++)
	{
		GradientRapide2((unsigned __int32*)Image,(unsigned __int32*)tamp,large,haute);
		FenetreGlissanteRapide(tamp,templ,tempc,large,haute,12,12);
	}
/*	time=GetCurrentTime()-time;
	CString text;
	text.Format("temps: %d",time);
	MessageBox(text,"Gradient");
*/	stop=6+6*Larg;;
	l=stop;
	for(i=0;i<large*haute-stop;i++) gradMoy->Image[i]=tempc[l++]/144;
	for(i=large*haute-stop;i<large*haute;i++) gradMoy->Image[i]=0;
	gradMoy->Create("Gradient Moy glis 12*12 de "+Titre, WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());
	delete templ;
	delete tempc;
	delete tamp;
//	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::LuminanceGlissante1212()
{
//	SetCursor(LoadCursor(0,IDC_WAIT));
	
	int i;
	unsigned __int16 *templ,*tempc;
	int large=Larg-Larg%4;
	int haute=Haut-Haut%4;
	CTraitImageWnd * lumMoy;
	long stop ,l;

	lumMoy=new CTraitImageWnd(large,haute,BMP);
	templ=new unsigned __int16[large*haute];
	tempc=new unsigned __int16[large*haute];

//	DWORD time=GetCurrentTime();
//	for(i=0;i<100;i++)
	FenetreGlissanteRapide2(Image,templ,tempc,large,haute,12,12);
/*	time=GetCurrentTime()-time;
	CString text;
	text.Format("temp: %d",time);
	MessageBox(text,"Luminance");
*/	stop=6+6*Larg;
	l=stop;
	for(i=0;i<large*haute-stop;i++) lumMoy->Image[i]=tempc[l++]/144;
	for(i=large*haute-stop;i<large*haute;i++) lumMoy->Image[i]=0;
	lumMoy->Create("Luminance Moy glis 12*12 de "+Titre, WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());
	delete templ;
	delete tempc;
//	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::ContrasteGlissant1212()
{
//	SetCursor(LoadCursor(0,IDC_WAIT));
	long stop,l;
	unsigned char *minimum,*maximum,min,max;
	CTraitImageWnd *contraste;

	minimum=new unsigned char[Larg*Haut];
	maximum=new unsigned char[Larg*Haut];
	contraste=new CTraitImageWnd(Larg,Haut,BMP);
//	DWORD time=GetCurrentTime();
//	for(int j=0;j<100;j++)
	{
		MaximumGlissant1212(maximum);
		MinimumGlissant1212(minimum);
		stop=6+6*Larg;
		l=stop;
		for(int i=0;i<Larg*Haut-stop;i++)
		{
			min=minimum[l];
			max=maximum[l++];

			if((min==0)&&(max==0)) contraste->Image[i]=0;
			else contraste->Image[i]=(unsigned char)((float)255*(max-min)/(max+min));
			
		}
	}
/*	time=GetCurrentTime()-time;
	CString text;
	text.Format("temps: %d",time);
	MessageBox(text,"Contraste");
*/
	for(int i=Larg*Haut-stop;i<Larg*Haut;i++) contraste->Image[i]=0;
	delete minimum;
	delete maximum;
	contraste->Create("contraste glis 12*12 de "+Titre, WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());
//	SetCursor(LoadCursor(0,IDC_ARROW));
}


/////////////////////////////////////////////////////////////
////////////////calculs "rapides"///////////////////////////

void CTraitImageWnd::FenetreGlissanteRapide2(unsigned char* bram,unsigned __int16* cufen,
											unsigned __int16* lufen,int large,int haute,int LF,int HF)
{
	unsigned __int32 p,pth,ptl,pb;
	unsigned __int32 ph,pl;
	unsigned __int16 som;
	unsigned __int32 *lram;
	unsigned __int32 som2,*dcufen,*dlufen;
	int l,c,c4,lf4,c2,c4tot;//44 msc

	lram=(unsigned __int32 *)bram;
	////////colonnes
	c4=large/4;lf4=LF/4;c4tot=Larg/4;
	lram=(unsigned __int32 *)bram;
	pb=0;
	for(l=0;l<haute;l++)
	{	som=0;p=l*Larg;
		ptl=l*c4tot;pth=ptl+lf4;
		for(c=0;c<LF;c++) 
		{som=som+bram[p++];cufen[pb++]=som;}
		for(c=lf4;c<c4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			som=som+unsigned __int16((ph&0xff)-(pl&0xff));
			cufen[pb++]=som;
			som=som+unsigned __int16(((ph&0xff00)>>8)-((pl&0xff00)>>8));
			cufen[pb++]=som;
			som=som+unsigned __int16(((ph&0xff0000)>>16)-((pl&0xff0000)>>16));
			cufen[pb++]=som;
			som=som+unsigned __int16((ph>>24)-(pl>>24));
			cufen[pb++]=som;
			pth++;ptl++;
		/*
			_asm
			{	
				;mov esi,pb 
				;shl esi,1	;multiplie *2 l'index pb
				mov ebx,ph
				mov ecx,pl
				mov dx,som
				;***********
				movzx ax,bl
				add dx,ax
				movzx ax,cl
				sub dx,ax
				mov word ptr cufen[esi],dx
				;***********
				movzx ax,bh
				add dx,ax
				movzx ax,ch
				sub dx,ax
				mov word ptr cufen[esi+2],dx
				;***********
				bswap ebx
				bswap ecx
				;***********
				movzx ax,bh
				add dx,ax
				movzx ax,ch
				sub dx,ax
				mov word ptr cufen[esi+4],dx
				;***********
				movzx ax,bl
				add dx,ax
				movzx ax,cl
				sub dx,ax
				mov word ptr cufen[esi+6],dx
				;***********
				mov som,dx
			}
			pth++;ptl++;
			pb+=4;
			*/
		}//for c
	}//for l

//////////////////////ligne
	dcufen=(unsigned __int32 *)cufen;
	dlufen=(unsigned __int32 *)lufen;
	c2=large/2;
	for(c=0;c<c2;c++)
	{	pth=c;som2=0;ptl=c;
		for(l=0;l<HF;l++) 
		{	som2=som2+dcufen[pth];//attention pas de carry
			dlufen[pth]=som2;
			pth+=c2;
		}//for l
		
		for(l=HF;l<haute;l++)
		{	som2=som2+(dcufen[pth]-dcufen[ptl]);
			dlufen[pth]=som2;
			pth+=c2;ptl+=c2;
		}//for l
	}//for c
}

void CTraitImageWnd::GradientRapide2(unsigned __int32 *lram,unsigned __int32 *bram,int large,int haute)
{
	int l,c,c4,c4tot;
	unsigned long pth,ptl,pt;
	unsigned __int32 ph,pl,gradph;
	unsigned char newph,newpl,oldph=0,oldpl=0;
	c4=large/4;
	c4tot=Larg/4;
	pt=0;
	for(l=0;l<(haute-1);l++)
	{
		pth=l*c4tot;ptl=pth+c4tot;
		for(c=0;c<c4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			_asm
			{
				mov ebx,ph
				mov ecx,pl
				bswap ecx
				bswap ebx
				mov newph,bl
				mov newpl,cl
				xor eax,eax
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov ah,oldph
				mov bl,TABABSO[eax]
				mov ax,cx
				mov ah,oldpl
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				;	ror ecx,8
				;***************
				mov al,newph
				mov oldph,al
				mov ah,newpl
				mov oldpl,ah
				bswap ebx
				mov gradph,ebx
				;xchg ah,al	;12->21
				;bswap eax	;1234->4321
			}
			bram[pt++]=gradph;
			pth++;ptl++;
		}//for c
	}//for l
	
}
void CTraitImageWnd::GradientRapide3(unsigned __int32 *lram,unsigned __int32 *bram)
{
	int l,c,c4;
	unsigned long pth,ptl;
	unsigned __int32 ph,pl,gradph;
	unsigned char newph,newpl,oldph=0,oldpl=0;
	c4=Larg/4;
	pth=0;ptl=pth+c4;
	for(l=0;l<(Haut-1);l++)
	{ for(c=0;c<c4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			_asm
			{
				mov ebx,ph
				mov ecx,pl
				bswap ecx
				bswap ebx
				mov newph,bl
				mov newpl,cl
				xor eax,eax
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov ah,oldph
				mov bl,TABABSO[eax]
				mov ax,cx
				mov ah,oldpl
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				;	ror ecx,8
				;***************
				mov al,newph
				mov oldph,al
				mov ah,newpl
				mov oldpl,ah
				bswap ebx
				mov gradph,ebx
				;xchg ah,al	;12->21
				;bswap eax	;1234->4321
			}
			bram[pth]=gradph;
			pth++;ptl++;
		}//for c
	}//for l
}

void CTraitImageWnd::FenetreGlissanteRapide(unsigned char* bram,unsigned __int16* cufen,
											unsigned __int16* lufen,int large,int haute,int LF,int HF)
{
	unsigned __int32 p,pth,ptl,pb;
	unsigned __int32 ph,pl;
	unsigned __int16 som;
	unsigned __int32 *lram;
	unsigned __int32 som2,*dcufen,*dlufen;
	int l,c,c4,lf4,c2;//44 msc

	lram=(unsigned __int32 *)bram;
	////////colonnes
	c4=large/4;lf4=LF/4;
	lram=(unsigned __int32 *)bram;
	pb=0;ptl=0;pth=ptl+lf4;
	for(l=0;l<haute;l++)
	{	som=0;p=l*large;
		for(c=0;c<LF;c++) 
		{som=som+bram[p++];cufen/*ESSAI*/[pb++]=som;}
		for(c=lf4;c<c4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			som=som+unsigned __int16((ph&0xff)-(pl&0xff));
			cufen[pb++]=som;
			som=som+unsigned __int16(((ph&0xff00)>>8)-((pl&0xff00)>>8));
			cufen[pb++]=som;
			som=som+unsigned __int16(((ph&0xff0000)>>16)-((pl&0xff0000)>>16));
			cufen[pb++]=som;
			som=som+unsigned __int16((ph>>24)-(pl>>24));
			cufen[pb++]=som;
			pth++;ptl++;
		/*
			_asm
			{	
				mov esi,pb 
				shl esi,1	;multiplie *2 l'index pb
				mov ebx,ph
				mov ecx,pl
				mov dx,som
				;***********
				movzx ax,bl
				add dx,ax
				movzx ax,cl
				sub dx,ax
				mov word ptr ESSAI[esi],dx
				;***********
				movzx ax,bh
				add dx,ax
				movzx ax,ch
				sub dx,ax
				mov word ptr ESSAI[esi+2],dx
				;***********
				bswap ebx
				bswap ecx
				;***********
				movzx ax,bh
				add dx,ax
				movzx ax,ch
				sub dx,ax
				mov word ptr ESSAI[esi+4],dx
				;***********
				movzx ax,bl
				add dx,ax
				movzx ax,cl
				sub dx,ax
				mov word ptr ESSAI[esi+6],dx
				;***********
				mov som,dx
			}
			pth++;ptl++;
			pb+=4;
		*/	
		}//for c
		pth+=lf4;ptl+=lf4;
	}//for l

//////////////////////ligne
	dcufen=(unsigned __int32 *)cufen/*ESSAI*/;
	dlufen=(unsigned __int32 *)lufen;
	c2=large/2;
	for(c=0;c<c2;c++)
	{	pth=c;som2=0;ptl=c;
		for(l=0;l<HF;l++) 
		{	som2=som2+dcufen[pth];//attention pas de carry
			dlufen[pth]=som2;
			pth+=c2;
		}//for l
		
		for(l=HF;l<haute;l++)
		{	som2=som2+(dcufen[pth]-dcufen[ptl]);
			dlufen[pth]=som2;
			pth+=c2;ptl+=c2;
		}//for l
	}//for c
}

void CTraitImageWnd::GradientRapide(unsigned __int32 *lram,int large,int haute)
{
	int l,c,c4;
	unsigned long pth,ptl;
	unsigned __int32 ph,pl,gradph;
	unsigned char newph,newpl,oldph=0,oldpl=0;
	c4=large/4;
	pth=0;ptl=pth+c4;
	for(l=0;l<(haute-1);l++)
	{ for(c=0;c<c4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			_asm
			{
				mov ebx,ph
				mov ecx,pl
				bswap ecx
				bswap ebx
				mov newph,bl
				mov newpl,cl
				xor eax,eax
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov ah,oldph
				mov bl,TABABSO[eax]
				mov ax,cx
				mov ah,oldpl
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				;	ror ecx,8
				;***************
				mov al,newph
				mov oldph,al
				mov ah,newpl
				mov oldpl,ah
				bswap ebx
				mov gradph,ebx
				;xchg ah,al	;12->21
				;bswap eax	;1234->4321
			}
			lram[pth]=gradph;
			pth++;ptl++;
		}//for c
	}//for l
	
}

long CTraitImageWnd::FenetreGlissanteRapideMoy(unsigned char* bram,unsigned __int16* cufen,
											unsigned __int16* lufen,int large,int haute,int LF,int HF)
{
	unsigned __int32 p,pth,ptl,pb,moy;
	unsigned __int32 ph,pl;
	unsigned __int16 som;
	unsigned __int32 *lram;
	unsigned __int32 som2,*dcufen,*dlufen;
	int l,c,c4,lf4,c2;//44 msc
	
	lram=(unsigned __int32 *)bram;
	////////colonnes
	c4=large/4;lf4=LF/4;
	lram=(unsigned __int32 *)bram;
	pb=0;ptl=0;pth=ptl+lf4;
	for(l=0;l<haute;l++)
	{	som=0;p=l*large;
		for(c=0;c<LF;c++) 
		{som=som+bram[p++];cufen[pb++]=som;}
		for(c=lf4;c<c4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			som=som+unsigned __int16((ph&0xff)-(pl&0xff));
			cufen[pb++]=som;
			som=som+unsigned __int16(((ph&0xff00)>>8)-((pl&0xff00)>>8));
			cufen[pb++]=som;
			som=som+unsigned __int16(((ph&0xff0000)>>16)-((pl&0xff0000)>>16));
			cufen[pb++]=som;
			som=som+unsigned __int16((ph>>24)-(pl>>24));
			cufen[pb++]=som;
			pth++;ptl++;
		}//for c
		pth+=lf4;ptl+=lf4;
	}//for l

//////////////////////ligne
	dcufen=(unsigned __int32 *)cufen;
	dlufen=(unsigned __int32 *)lufen;
	c2=large/2;moy=0;
	for(c=0;c<c2;c++)
	{	pth=c;som2=0;ptl=c;
		for(l=0;l<HF;l++) 
		{	som2=som2+dcufen[pth];//attention pas de carry
			dlufen[pth]=som2;
			///
			moy+=lufen[2*pth]+lufen[2*pth+1];
			///
			pth+=c2;
			
		}//for l
		
		for(l=HF;l<haute;l++)
		{	som2=som2+(dcufen[pth]-dcufen[ptl]);
			dlufen[pth]=som2;
			///
			moy+=lufen[2*pth]+lufen[2*pth+1];
			///
			pth+=c2;ptl+=c2;
		}//for l
	}//for c
	moy/=large*haute;
	return moy;
}

long CTraitImageWnd::FenetreGlissanteRapideMoy2(unsigned char* bram,unsigned __int16* cufen,
											unsigned __int16* lufen,int large,int haute,int LF,int HF)
{
	unsigned __int32 p,pth,ptl,pb,moy;
	unsigned __int32 ph,pl;
	unsigned __int16 som;
	unsigned __int32 *lram;
	unsigned __int32 som2,*dcufen,*dlufen;
	int l,c,c4,lf4,c2,c4tot;;//44 msc
	
	lram=(unsigned __int32 *)bram;
	////////colonnes
	c4=large/4;lf4=LF/4;c4tot=Larg/4;
	lram=(unsigned __int32 *)bram;
	pb=0;
	for(l=0;l<haute;l++)
	{	som=0;p=l*Larg;
		ptl=l*c4tot;pth=ptl+lf4;
		for(c=0;c<LF;c++) 
		{som=som+bram[p++];cufen[pb++]=som;}
		for(c=lf4;c<c4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			som=som+unsigned __int16((ph&0xff)-(pl&0xff));
			cufen[pb++]=som;
			som=som+unsigned __int16(((ph&0xff00)>>8)-((pl&0xff00)>>8));
			cufen[pb++]=som;
			som=som+unsigned __int16(((ph&0xff0000)>>16)-((pl&0xff0000)>>16));
			cufen[pb++]=som;
			som=som+unsigned __int16((ph>>24)-(pl>>24));
			cufen[pb++]=som;
			pth++;ptl++;
			/*som=som+unsigned __int16((ph>>24)-(pl>>24));
			cufen[pb++]=som;
			som=som+unsigned __int16(((ph&0xff0000)>>16)-((pl&0xff0000)>>16));
			cufen[pb++]=som;
			som=som+unsigned __int16(((ph&0xff00)>>8)-((pl&0xff00)>>8));
			cufen[pb++]=som;
			som=som+unsigned __int16((ph&0xff)-(pl&0xff));
			cufen[pb++]=som;
			pth++;ptl++;*/
		
		/*	__asm
			{	
				mov esi,pb 
				shl esi,1	;multiplie *2 l'index pb
				mov ebx,ph
				mov ecx,pl
				mov dx,som
				;***********
				movzx ax,bl
				add dx,ax
				movzx ax,cl
				sub dx,ax
				mov word ptr cufen[esi],dx
				;***********
				movzx ax,bh
				add dx,ax
				movzx ax,ch
				sub dx,ax
				mov word ptr cufen[esi+2],dx
				;***********
				bswap ebx
				bswap ecx
				;***********
				movzx ax,bh
				add dx,ax
				movzx ax,ch
				sub dx,ax
				mov word ptr cufen[esi+4],dx
				;***********
				movzx ax,bl
				add dx,ax
				movzx ax,cl
				sub dx,ax
				mov word ptr cufen[esi+6],dx
				;***********
				mov som,dx
			}
			pb=pb+4;
			pth++;ptl++;*/
		}//for c
	
	}//for l

//////////////////////ligne
	dcufen=(unsigned __int32 *)cufen;
	dlufen=(unsigned __int32 *)lufen;
	c2=large/2;moy=0;
	for(c=0;c<c2;c++)
	{	pth=c;som2=0;ptl=c;
		for(l=0;l<HF;l++) 
		{	som2=som2+dcufen[pth];//attention pas de carry
			dlufen[pth]=som2;
			///
			moy+=lufen[2*pth]+lufen[2*pth+1];
			///
			pth+=c2;
			
		}//for l
		
		for(l=HF;l<haute;l++)
		{	som2=som2+(dcufen[pth]-dcufen[ptl]);
			dlufen[pth]=som2;
			///
			moy+=lufen[2*pth]+lufen[2*pth+1];
			///
			pth+=c2;ptl+=c2;
		}//for l
	}//for c
	moy/=large*haute;
	return moy;
}

////////////////////////////////////////////////////////////////////////////
//////////////////////// masque//////////////////////
////////////////////////////////////////////////////////////////////////////
void CTraitImageWnd::DessineRectangle()
{
	if(Masque)
		{
			CClientDC DC(this);
			CPoint h1,h2; 
			
			CPen stylo(PS_SOLID,1,RGB(0,0,255));
			DC.SelectObject(&stylo); 
			DC.SelectStockObject(NULL_BRUSH);
			
			h1 = Inv_Get_ij((int)debut_i,(int)debut_j); 
			h2 = Inv_Get_ij((int)fin_i,(int)fin_j);
			CRect rect(h1,h2);
			DC.Ellipse(rect);
			if(Decoupage)
			{
				CRect rect2;
				CSize difference;
				rect2=rect;
				difference.cx=(long)(rect.Width()*0.15843497);
				difference.cy=(long)(rect.Height()*0.15843497);
				rect2.DeflateRect(difference);
				DC.Ellipse(rect2);
				
				rect2=rect;
				difference.cx=(long)(rect.Width()*0.2763932);
				difference.cy=(long)(rect.Height()*0.2763932);
				rect2.DeflateRect(difference);
				DC.Ellipse(rect2);

				rect2=rect;
				difference.cx=(long)(rect.Width()*0.3709005);
				difference.cy=(long)(rect.Height()*0.3709005);
				rect2.DeflateRect(difference);
				DC.Ellipse(rect2);
			}
			ReleaseDC(&DC); 
		}
	else CImageWnd::DessineRectangle();
}

void CTraitImageWnd::OnMasqueRechercher() 
{
	long l,i1,i2,j1,j2,i3,j3,i,j;
	int Nbpoint=0;
	BOOL flag;

	ArreteLeContinu();
	OnMasqueAnnuler();
	
	Seuil=40;
	CDialSeuil dial (&Seuil,&flag); 
	if (dial.DoModal()==IDOK)
	{
		SetCursor(LoadCursor(0,IDC_WAIT));
		

		for( j=debut_j;j<=fin_j;j++)
		{for( i=debut_i;i<=fin_i;i++)
			{   
				l=j*Larg+i;
				if (((Image[l]>Seuil)&&flag)||((Image[l]<Seuil)&&(!flag))) {i1=i;j1=j;Nbpoint=1;break;}
				
			}
		if (Nbpoint==1) break;
		 }//détection du premier point

		if(Nbpoint==0) {MessageBox("Impossible de trouver le masque");return;}

		for(j=fin_j;j>j1;j--)
		{   
				l=j*Larg+i1;
				if(((Image[l]>Seuil)&&flag)||((Image[l]<Seuil)&&(!flag))) {j2=j;Nbpoint=2;break;}
				
		}//détection du second point en face

		if(Nbpoint==1) {MessageBox("Impossible de trouver le masque");return;}

		for( i=debut_i;i<=fin_i;i++)
		{for( j=debut_j;j<=fin_j;j++)
			{   
				l=j*Larg+i;
				if(((Image[l]>Seuil)&&flag)||((Image[l]<Seuil)&&(!flag))) {i2=i;j3=j;Nbpoint=3;break;}
				
			}
			if (Nbpoint==3) break;
		 }//détection du troisième point

		if(Nbpoint==2) {MessageBox("Impossible de trouver le masque");return;}

		for(i=fin_i;i>i2;i--)
		{   
				l=j3*Larg+i;
				if(((Image[l]>Seuil)&&flag)||((Image[l]<Seuil)&&(!flag))) {i3=i;Nbpoint=4;break;}
				
		}//détection du quatrième point en face du troisième

		if(Nbpoint==3) {MessageBox("Impossible de trouver le masque");return;}
		Masque=TRUE;
		menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(6,MF_ENABLED|MF_BYPOSITION);
		menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(7,MF_ENABLED|MF_BYPOSITION);
		debut_i	= i2;
		fin_i	= i3;
		debut_j	= j1;
		fin_j	= j2;
		zoneSelectionnee=TRUE;
		SetCursor(LoadCursor(0,IDC_ARROW));
		Invalidate(FALSE);
		Seuil=128;
	}
}

void CTraitImageWnd::OnMasqueLamain() 
{
	Masque=TRUE;
	menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(6,MF_ENABLED|MF_BYPOSITION);
	menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(7,MF_ENABLED|MF_BYPOSITION);

	ArreteLeContinu();
	Memorise();
	CImageWnd::OnSelectionneZone();
//	CImageWnd::OnSelectionneZone();
}

void CTraitImageWnd::OnMasqueRechercheAutom(int m_seuil) 
{
	long l,i1,i2,j1,j2,i3,j3,i,j;
	int Nbpoint=0;
	

	OnMasqueAnnuler();
	ArreteLeContinu();


	SetCursor(LoadCursor(0,IDC_WAIT));
	

	for( j=debut_j;j<=fin_j;j++)
	{for( i=debut_i;i<=fin_i;i++)
		{   
			l=j*Larg+i;
			if (Image[l]>m_seuil) {i1=i;j1=j;Nbpoint=1;break;};
			
		}
	if (Nbpoint==1) break;
	 }//détection du premier point

	if(Nbpoint==0) {MessageBox("Impossible de trouver le masque");return;}

	for(j=fin_j;j>j1;j--)
	{   
			l=j*Larg+i1;
			if(Image[l]>m_seuil) {j2=j;Nbpoint=2;break;};
			
	}//détection du second point en face

	if(Nbpoint==1) {MessageBox("Impossible de trouver le masque");return;}

	for( i=debut_i;i<=fin_i;i++)
	{for( j=debut_j;j<=fin_j;j++)
		{   
			l=j*Larg+i;
			if(Image[l]>m_seuil) {i2=i;j3=j;Nbpoint=3;break;};
			
		}
		if (Nbpoint==3) break;
	 }//détection du troisième point

	if(Nbpoint==2) {MessageBox("Impossible de trouver le masque");return;}

	for(i=fin_i;i>i2;i--)
	{   
			l=j3*Larg+i;
			if(Image[l]>m_seuil){i3=i;Nbpoint=4;break;};
			
	}//détection du quatrième point en face du troisième

	if(Nbpoint==3) {MessageBox("Impossible de trouver le masque");return;}
	Masque=TRUE;
	menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(6,MF_ENABLED|MF_BYPOSITION);
	menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(7,MF_ENABLED|MF_BYPOSITION);

	debut_i	= i2;
	fin_i	= i3;
	debut_j	= j1;
	fin_j	= j2;
	
	SetCursor(LoadCursor(0,IDC_ARROW));
	Invalidate(FALSE);	
}
void CTraitImageWnd::OnMasqueAnnuler() 
{
	Masque=FALSE;
	menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(6,MF_GRAYED|MF_BYPOSITION);
	menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(7,MF_GRAYED|MF_BYPOSITION);

	Decoupage=FALSE;
	RegionEllipse.DeleteObject();
	OnRestaureZone();
}

void CTraitImageWnd::OnUpdateMasqueAnnuler(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);
}

void CTraitImageWnd::OnTraitementSupprimer() 
{
	CPoint h;
	OnEditCopy();
	if(Masque)
	{
		AffectRegionMasque();
		ArreteLeContinu();
		for(int j=debut_j;j<=fin_j;j++)
		for(int i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			if(RegionEllipse.PtInRegion(h))
			{	
				Image[i+Larg*j]=255;
			}
		}
	}
	else
	{
		for (long i=debut_i;i<=fin_i;i++) 	
		for(long j=debut_j;j<=fin_j;j++) Image[i+j*Larg]=255;
	}
	RegionEllipse.DeleteObject();
	Masque=FALSE;
	menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(6,MF_GRAYED|MF_BYPOSITION);
	menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(7,MF_GRAYED|MF_BYPOSITION);

	OnRestaureZone();
}

void CTraitImageWnd::OnUpdateTraitementSupprimer(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);
}

void CTraitImageWnd::AffectRegionMasque()
{ 
	CPoint h1,h2;
	RegionEllipse.DeleteObject();
	h1 = Inv_Get_ij((int)debut_i,(int)debut_j-1); 
	h2 = Inv_Get_ij((int)fin_i+1,(int)fin_j);
	RegionEllipse.CreateEllipticRgnIndirect(CRect(h1,h2));	
}

float CTraitImageWnd::OnMasqueLuminancemoyenneAutom() 
{
	long moy=0,nbpoint=0;
	float result;
	CPoint h;
	char ch[50];

	SetCursor(LoadCursor(0,IDC_WAIT));  
	AffectRegionMasque();
	Memorise();
	for(int j=debut_j;j<=fin_j;j++)
		for(int i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			if(RegionEllipse.PtInRegion(h))
			{	
				
				moy+=Image[i+j*Larg];
				nbpoint++;
				Image[i+Larg*j]=255-Image[i+Larg*j];
				
			}
		}
	Invalidate(FALSE);
	result=(float) moy/nbpoint;
	sprintf(ch,"luminance moyenne=%3.2f",result);
	SetCursor(LoadCursor(0,IDC_ARROW)); 
//	MessageBox(ch,"Masque");
	OnRestaureImage();
	return result;	 
}
void CTraitImageWnd::OnMasqueInversion() 
{
	CPoint h;
	AffectRegionMasque();
	Memorise();
	SetCursor(LoadCursor(0,IDC_WAIT));  
	for(int j=debut_j;j<=fin_j;j++)
		for(int i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			if(RegionEllipse.PtInRegion(h))
			{	
				Image[i+Larg*j]=255-Image[i+Larg*j];
			}
		}
	Invalidate(FALSE);
}

void CTraitImageWnd::OnUpdateMasqueInversion(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);
}

void CTraitImageWnd::OnMasqueContraste() 
{
	unsigned char maxi,mini,temp;
	double contraste;
	CPoint h;
	char ch[50];
	maxi=0;
	mini=255;

	AffectRegionMasque();
	Memorise();
	SetCursor(LoadCursor(0,IDC_WAIT));  
	for(int j=debut_j;j<=fin_j;j++)
		for(int i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			if(RegionEllipse.PtInRegion(h))
			{	
				
				temp=Image[i+j*Larg];
				if(temp>maxi) maxi=temp;
				if (temp<mini) mini=temp;
				Image[i+Larg*j]=255-Image[i+Larg*j];
				
			}
		}
	Invalidate(FALSE);
	if((maxi+mini)==0) contraste=0;
	else contraste= ((double) (maxi-mini)/(maxi+mini))*100;
	sprintf(ch,"contraste= %3.2f%%",contraste);
	SetCursor(LoadCursor(0,IDC_ARROW));
	MessageBox(ch,"Masque");
	OnRestaureImage();
}

void CTraitImageWnd::OnUpdateMasqueContraste(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);
}

void CTraitImageWnd::OnMasqueMaxmin() 
{
	unsigned char maxi,mini,temp;
	double contraste;
	CPoint h;
	char ch[50];
	maxi=0;
	mini=255;

	SetCursor(LoadCursor(0,IDC_WAIT));  
	Memorise();
	AffectRegionMasque();

	
	for(int j=debut_j;j<=fin_j;j++)
		for(int i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			if(RegionEllipse.PtInRegion(h))
			{	
				
				temp=Image[i+j*Larg];
				if(temp>maxi) maxi=temp;
				if (temp<mini) mini=temp;
				Image[i+Larg*j]=255-Image[i+Larg*j];
				
			}
		}
	Invalidate(FALSE);
	if((maxi+mini)==0) contraste=0;
	else contraste= ((double) (maxi-mini)/255.0)*100;
	sprintf(ch,"max-min= %3.2f%%",contraste);
	SetCursor(LoadCursor(0,IDC_ARROW));
	MessageBox(ch,"Masque");
	OnRestaureImage();
}

void CTraitImageWnd::OnUpdateMasqueMaxmin(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);	
}

void CTraitImageWnd::OnMasqueEntropie() 
{
	long p,ligne,colo,nbpoint=0,nbniveau;
	double d;
	float *pilnpi;
	long  *cpt;
	double som;
	char ch[20];
	int i,j;
	CPoint h;
    CDialVal DialVal;

	SetCursor(LoadCursor(0,IDC_WAIT));  
	AffectRegionMasque();
	Memorise();

	if(DialVal.DoModal()==IDOK) nbniveau=atoi(DialVal.m_val);
	else return;
	if(nbniveau<=0) {MessageBox("erreur dans les parametres","Entropie",MB_OK|MB_ICONSTOP);return;} 
	
	for( j=debut_j;j<=fin_j;j++)
		for( i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			if(RegionEllipse.PtInRegion(h))
			{	
				nbpoint++;
			}
		}

	pilnpi=new float[nbpoint+1];
	cpt=new long[nbniveau];


	som=0;
	pilnpi[0]=(float)0.0;

	for ( i=1;i<=nbpoint;i++)  
	{
		d=exp(1.0)*((double) i/nbpoint)*log((double) (i)/(nbpoint));
		pilnpi[i]=(float)fabs(d);
	}	

	for (i=0;i<nbniveau;i++) cpt[i]=0; 

	for (ligne=debut_j;ligne<=fin_j;ligne++)
	{   
		for (colo=debut_i;colo<=fin_i;colo++)
		{
			p=Larg*ligne+colo;
			h=Inv_Get_ij(colo,ligne);
			if(RegionEllipse.PtInRegion(h))
			{	
				cpt[(Image[p]*nbniveau)/256]++;
				Image[p]=(Image[p]*nbniveau)/256;
			}
			
		} // colo
	} // ligne	
	
	for (i=0;i<nbniveau;i++)
	{
		som=som+pilnpi[cpt[i]];
	}
	Invalidate(FALSE);
	som=(100*som)/nbniveau;
	sprintf(ch,"Entropie moyenne=%2.1f%%",som);
	SetCursor(LoadCursor(0,IDC_ARROW));
	MessageBox(ch,"Masque");
	delete pilnpi;
	delete cpt;
	OnRestaureImage();
}

void CTraitImageWnd::OnUpdateMasqueEntropie(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);	
}

void CTraitImageWnd::OnMasqueLuminancemoyenne() 
{
	long moy=0,nbpoint=0;
	double result;
	CPoint h;
	char ch[50];

	SetCursor(LoadCursor(0,IDC_WAIT));  
	AffectRegionMasque();
	Memorise();
	for(int j=debut_j;j<=fin_j;j++)
		for(int i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			if(RegionEllipse.PtInRegion(h))
			{	
				
				moy+=Image[i+j*Larg];
				nbpoint++;
				Image[i+Larg*j]=255-Image[i+Larg*j];
				
			}
		}
	Invalidate(FALSE);
	result=(double) moy/nbpoint;
	sprintf(ch,"luminance moyenne=%3.2f",result);
	SetCursor(LoadCursor(0,IDC_ARROW)); 
	MessageBox(ch,"Masque");
	OnRestaureImage();
}

void CTraitImageWnd::OnUpdateMasqueLuminancemoyenne(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);	
}

void CTraitImageWnd::OnMasqueCalculslocauxLuminancemoyenne() 
{
	double * luminance,moyenne,ecarttype,pourcent;
	int sygma;
	BOOL existe;
	CPoint h;
	char ch[80];
	int compteur=0,comp=0,lumi;
	BOOL l,k=FALSE;
	
	luminance=new double[(debut_i-fin_i+1)*(debut_j-fin_j+1)/40+1];

	AffectRegionMasque();
	Memorise();

	for(int j=debut_j;j<=fin_j;j+=20)
	{
		l=k;
		for(int i=debut_i;i<=fin_i;i+=20)
		{	
			lumi=0;
			comp=0;
			existe=FALSE;
			for(int li=0;li<20;li++)
			for(int co=0;co<20;co++)
			{
				if (((i+co)<=fin_i)||((li+j)<=fin_j))
				{
					h=Inv_Get_ij(i+co,j+li);

					if(RegionEllipse.PtInRegion(h))
					{	
						
						lumi+=Image[i+co+(li+j)*Larg];
						comp++;
						if(l) Image[i+co+Larg*(j+li)]=255-Image[i+co+Larg*(j+li)];
						existe=TRUE;
					}
				}
				
			}
			
			if(l) l=FALSE;
				else l=TRUE;
			if (existe) 
			{
				luminance[compteur++]=((double) lumi/comp);
				
			}
			
		}
		if(k)k=FALSE;
		else k=TRUE;
	}
	Invalidate(FALSE);

	moyenne=0;
	int i;
	for ( i=0;i<compteur;i++) moyenne+=luminance[i];
	moyenne/=compteur;

	ecarttype=0;
	for(i=0;i<compteur;i++) ecarttype+=pow((luminance[i]-moyenne),2);
	ecarttype/=compteur;
	ecarttype=sqrt(ecarttype);

	sygma=0;
	for(i=0;i<compteur;i++) if((luminance[i]<(moyenne+ecarttype/2))&&
		(luminance[i]>(moyenne+1-ecarttype/2))||(luminance[i]==moyenne)) sygma++;

	pourcent=(100*sygma)/compteur;

	sprintf(ch,"moyenne luminance=%3.2f\nécart type=%3.2f\npourcentage à 1 sygma =%3.2f"
		,moyenne,ecarttype,pourcent);
	MessageBox(ch,"Masque luminance");
	delete luminance;
	OnRestaureImage();
}

void CTraitImageWnd::OnUpdateMasqueCalculslocauxLuminancemoyenne(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);
}
void CTraitImageWnd::OnMasqueLocauxcontraste() 
{
	unsigned char maxi,mini,temp;
	double * contraste,moyenne,ecarttype,pourcent;
	int sygma;
	BOOL existe;
	CPoint h;
	char ch[80];
	int compteur=0;
	BOOL l,k=FALSE;

	contraste=new double[(debut_i-fin_i+1)*(debut_j-fin_j+1)/40+1];

	AffectRegionMasque();
	Memorise();
	int i,j;

	for( j=debut_j;j<=fin_j;j+=20)
	{
		l=k;
		for( i=debut_i;i<=fin_i;i+=20)
		{	
			maxi=0;
			mini=255;
			existe=FALSE;
			for(int li=0;li<20;li++)
			for(int co=0;co<20;co++)
			{
				if (((i+co)<=fin_i)||((li+j)<=fin_j))
				{
					h=Inv_Get_ij(i+co,j+li);

					if(RegionEllipse.PtInRegion(h))
					{	
						
						temp=Image[i+co+(li+j)*Larg];
						if(temp>maxi) maxi=temp;
						if (temp<mini) mini=temp;
						if(l) Image[i+co+Larg*(j+li)]=255-Image[i+co+Larg*(j+li)];
						existe=TRUE;
					}
				}
				
			}
			if(l) l=FALSE;
			else l=TRUE;
			
			if (existe) 
			{
				if((maxi+mini)==0) contraste[compteur++]=0;
				else contraste[compteur++]=((double) (maxi-mini)/(maxi+mini))*100;
			}
			
		}
		if(k)k=FALSE;
		else k=TRUE;
	}
	Invalidate(FALSE);

	moyenne=0;
	for (int i=0;i<compteur;i++) moyenne+=contraste[i];
	moyenne/=compteur;

	ecarttype=0;
	for(i=0;i<compteur;i++) ecarttype+=pow((contraste[i]-moyenne),2);
	ecarttype/=compteur;
	ecarttype=sqrt(ecarttype);

	sygma=0;
	for(i=0;i<compteur;i++) if((contraste[i]<(moyenne+ecarttype/2))&&
		(contraste[i]>(moyenne+1-ecarttype/2))||(contraste[i]==moyenne)) sygma++;

	pourcent=(100*sygma)/compteur;

	sprintf(ch,"moyenne du contraste=%3.2f\nécart type=%3.2f\npourcentage à 1 sygma =%3.2f"
		,moyenne,ecarttype,pourcent);
	MessageBox(ch,"Masque contraste");
	delete contraste;
	OnRestaureImage();
}

void CTraitImageWnd::OnUpdateMasqueLocauxcontraste(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);	
}

void CTraitImageWnd::OnMasqueCalculslocauxMaxmin() 
{
	unsigned char maxi,mini,temp;
	double * contraste,moyenne,ecarttype,pourcent;
	int sygma;
	BOOL existe;
	CPoint h;
	char ch[80];
	int compteur=0;
	BOOL l,k=FALSE;

	contraste=new double[(debut_i-fin_i+1)*(debut_j-fin_j+1)/40+1];

	AffectRegionMasque();
	Memorise();

	for(int j=debut_j;j<=fin_j;j+=20)
	{
		l=k;
		for(int i=debut_i;i<=fin_i;i+=20)
		{	
			maxi=0;
			mini=255;
			existe=FALSE;
			for(int li=0;li<20;li++)
			for(int co=0;co<20;co++)
			{
				if (((i+co)<=fin_i)||((li+j)<=fin_j))
				{
					h=Inv_Get_ij(i+co,j+li);

					if(RegionEllipse.PtInRegion(h))
					{	
						
						temp=Image[i+co+(li+j)*Larg];
						if(temp>maxi) maxi=temp;
						if (temp<mini) mini=temp;
						if(l) Image[i+co+Larg*(j+li)]=255-Image[i+co+Larg*(j+li)];
						existe=TRUE;
					}
				}
				
			}
			if(l) l=FALSE;
			else l=TRUE;
			
			if (existe) 
			{
				contraste[compteur++]=((double) (maxi-mini)/255)*100;
			}
			
		}
		if(k)k=FALSE;
		else k=TRUE;
	}
	Invalidate(FALSE);

	moyenne=0;
	int i;
	for ( i=0;i<compteur;i++) moyenne+=contraste[i];
	moyenne/=compteur;

	ecarttype=0;
	for(i=0;i<compteur;i++) ecarttype+=pow((contraste[i]-moyenne),2);
	ecarttype/=compteur;
	ecarttype=sqrt(ecarttype);

	sygma=0;
	for(i=0;i<compteur;i++) if((contraste[i]<(moyenne+ecarttype/2))&&
		(contraste[i]>(moyenne+1-ecarttype/2))||(contraste[i]==moyenne)) sygma++;

	pourcent=(100*sygma)/compteur;

	sprintf(ch,"moyenne du max-min=%3.2f\nécart type=%3.2f\npourcentage à 1 sygma =%3.2f"
		,moyenne,ecarttype,pourcent);
	MessageBox(ch,"Masque max-min");
	delete contraste;
	OnRestaureImage();	
}

void CTraitImageWnd::OnUpdateMasqueCalculslocauxMaxmin(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);
}

void CTraitImageWnd::OnMasqueCalculslocauxEntropie() 
{
	double *pilnpi,*entropie,moyenne,ecarttype,pourcent,d,som;
	int sygma,nbniveau;
	long *cpt;
	CPoint h;
	char ch[500];
	int compteur=0,nbpoint=0,nb=0;
	BOOL l,k=FALSE;
	CDialVal DialVal;

	
	entropie=new double[(fin_i-debut_i+1)*(fin_j-debut_j+1)/40];

	AffectRegionMasque();
	Memorise();
	
	if(DialVal.DoModal()==IDOK) nbniveau=atoi(DialVal.m_val);
	else return;
	if(nbniveau<=0) {MessageBox("erreur dans les parametres","Entropie",MB_OK|MB_ICONSTOP);return;} 
	
	cpt=new long[nbniveau];

	for(int j=debut_j;j<=fin_j;j+=20)
	{
		l=k;
		for(int i=debut_i;i<=fin_i;i+=20)
		{	
			nbpoint=0;
			nb=0;
			for(int li=0;li<20;li++)
			for(int co=0;co<20;co++)
			{
				if (((i+co)<=fin_i)||((li+j)<=fin_j))
				{
					h=Inv_Get_ij(i+co,j+li);

					if(RegionEllipse.PtInRegion(h))
					{	
						nbpoint++;
					}
				}			
			}//for li for co
			
			int q,li,co;
			if (nbpoint>0) 
			{
				pilnpi=new double[nbpoint+1];
				
				pilnpi[0]=0.0;

				for (q=1;q<=nbpoint;q++)  
				{
					d=exp(1.0)*((double) q/nbpoint)*log((double) (q)/(nbpoint));
					pilnpi[q]=fabs(d);
				}	

				for (q=0;q<nbniveau;q++) cpt[q]=0; 

				for(li=0;li<20;li++)
				for(co=0;co<20;co++)
				{
					if (((i+co)<=fin_i)||((li+j)<=fin_j))
					{
						h=Inv_Get_ij(i+co,j+li);

						if(RegionEllipse.PtInRegion(h))
						{	
							cpt[(Image[i+co+Larg*(j+li)]*nbniveau)/256]++;
							nb++;
							if(l) Image[i+co+Larg*(j+li)]=255;
						}
					}			
				}//for li for co
				
				som=0;
				for(q=0;q<nbniveau;q++) som+=pilnpi[cpt[q]];
				Invalidate(FALSE);
				
				entropie[compteur++]=100*som/nbniveau;
				delete pilnpi;

					
			}//if existe
			if(l) l=FALSE;
				else l=TRUE;
		}//for int i
		if(k)k=FALSE;
		else k=TRUE;
	}//for int j
	Invalidate(FALSE);

	moyenne=0;
	int i;
	for ( i=0;i<compteur;i++) moyenne+=entropie[i];
	if (compteur!=0) moyenne/=compteur;

	ecarttype=0;
	for(i=0;i<compteur;i++) ecarttype+=pow((entropie[i]-moyenne),2);
	if (compteur!=0) 
	{
		ecarttype/=compteur;
		ecarttype=sqrt(ecarttype);
	}

	sygma=0;
	for(i=0;i<compteur;i++) if((entropie[i]<(moyenne+ecarttype/2))&&
		(entropie[i]>(moyenne+1-ecarttype/2))||(entropie[i]==moyenne)) sygma++;

	if (compteur!=0) pourcent=(100*sygma)/compteur;
	
	
	sprintf(ch,"nb carre=%d\nentropie=%3.2f\nécart type=%3.2f\npourcentage à 1 sygma =%3.2f"
		,compteur,moyenne,ecarttype,pourcent);
	MessageBox(ch,"Masque entropie");
	delete entropie;
	delete cpt;
	OnRestaureImage();
}

void CTraitImageWnd::OnUpdateMasqueCalculslocauxEntropie(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);	
}

void CTraitImageWnd::OnMasqueCalculslocauxGradientmoyen() 
{
	double * luminance,moyenne,ecarttype,pourcent;
	int sygma;
	BOOL existe;
	CPoint h;
	char ch[80];
	int compteur=0,comp=0,lumi;
	BOOL l,k=FALSE;
	
	luminance=new double[(debut_i-fin_i+1)*(debut_j-fin_j+1)/40+1];

	AffectRegionMasque();
	Memorise();
	OnTraitementGradientderob();

	for(int j=debut_j;j<=fin_j;j+=20)
	{
		l=k;
		for(int i=debut_i;i<=fin_i;i+=20)
		{	
			lumi=0;
			comp=0;
			existe=FALSE;
			for(int li=0;li<20;li++)
			for(int co=0;co<20;co++)
			{
				if (((i+co)<=fin_i)||((li+j)<=fin_j))
				{
					h=Inv_Get_ij(i+co,j+li);

					if(RegionEllipse.PtInRegion(h))
					{	
						
						lumi+=Image[i+co+(li+j)*Larg];
						comp++;
						if(l) Image[i+co+Larg*(j+li)]=255-Image[i+co+Larg*(j+li)];
						existe=TRUE;
					}
				}
				
			}
			
			if(l) l=FALSE;
				else l=TRUE;
			if (existe) 
			{
				luminance[compteur++]=((double) lumi/comp);
				
			}
			
		}
		if(k)k=FALSE;
		else k=TRUE;
	}
	Invalidate(FALSE);

	moyenne=0;
	int i;
	for ( i=0;i<compteur;i++) moyenne+=luminance[i];
	moyenne/=compteur;

	ecarttype=0;
	for(i=0;i<compteur;i++) ecarttype+=pow((luminance[i]-moyenne),2);
	ecarttype/=compteur;
	ecarttype=sqrt(ecarttype);

	sygma=0;
	for(i=0;i<compteur;i++) if((luminance[i]<(moyenne+ecarttype/2))&&
		(luminance[i]>(moyenne+1-ecarttype/2))||(luminance[i]==moyenne)) sygma++;

	pourcent=(100*sygma)/compteur;

	sprintf(ch,"moyenne gradient=%3.2f\nécart type=%3.2f\npourcentage à 1 sygma =%3.2f"
		,moyenne,ecarttype,pourcent);
	MessageBox(ch,"Masque gradient");
	delete luminance;
	OnRestaureImage();	
}

void CTraitImageWnd::OnUpdateMasqueCalculslocauxGradientmoyen(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);
}

void CTraitImageWnd::OnMasqueCalculslocauxTous() 
{
	OnMasqueLocauxcontraste();
	OnMasqueCalculslocauxMaxmin();
	OnMasqueCalculslocauxLuminancemoyenne();
	OnMasqueCalculslocauxEntropie();
	OnMasqueCalculslocauxGradientmoyen();
}

void CTraitImageWnd::OnUpdateMasqueCalculslocauxTous(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);	
}

void CTraitImageWnd::OnMasqueGradientmoyen() 
{
	long grad=0,nbpoint=0,temp;
	double result;
	CPoint h1,h2,h3,h4;
	char ch[50];

	SetCursor(LoadCursor(0,IDC_WAIT));  
	AffectRegionMasque();
	Memorise();
	for(int j=debut_j;j<=fin_j;j++)
		for(int i=debut_i;i<=fin_i;i++)
		{
			h1=Inv_Get_ij(i,j);
			h2=Inv_Get_ij(i+1,j);
			h3=Inv_Get_ij(i,j+1);
			h4=Inv_Get_ij(i+1,j+1);
			if((RegionEllipse.PtInRegion(h1))&&(RegionEllipse.PtInRegion(h2))&&
				(RegionEllipse.PtInRegion(h3))&&(RegionEllipse.PtInRegion(h4)))
			{	
				
				temp=labs(Image[i+j*Larg]-Image[i+1+(j+1)*Larg])+
					labs(Image[i+1+j*Larg]-Image[i+(j+1)*Larg]);
				grad+=temp;
				nbpoint++;
				Image[i+Larg*j]=(unsigned char)(temp/2.0);
				
			}
		}
	Invalidate(FALSE);
	result=(double) grad/(2*nbpoint);
	sprintf(ch,"gradient moyen=%3.2f",result);
	SetCursor(LoadCursor(0,IDC_ARROW));
	MessageBox(ch,"Masque");
	OnRestaureImage();
}

void CTraitImageWnd::OnUpdateMasqueGradientmoyen(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);	
}

void CTraitImageWnd::OnMasqueDimensions() 
{
	char ch[60];
	int largeur, hauteur;
	float rapport;

	AffectRegionMasque();
	hauteur=fin_j-debut_j+1;
	largeur=fin_i-debut_i+1;
	if (largeur!=0) rapport=(float) hauteur/largeur;
	else rapport=(float)0.0;
	sprintf(ch,"hauteur=%d \nlargeur=%d \nrapport=%1.2f",hauteur,largeur,rapport);
	MessageBox(ch,"Masque");
}

void CTraitImageWnd::OnUpdateMasqueDimensions(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);	
}

void CTraitImageWnd::OnMasqueTous() 
{
	OnMasqueContraste();
	OnMasqueMaxmin();
	OnMasqueEntropie();
	OnMasqueLuminancemoyenne();
	OnMasqueGradientmoyen();
	OnMasqueDimensions();
}

void CTraitImageWnd::OnUpdateMasqueTous(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);	
}

void CTraitImageWnd::OnMasqueNbpixels() 
{
	char ch[60];
	int l=0;
	CPoint h;

	SetCursor(LoadCursor(0,IDC_WAIT));  
	AffectRegionMasque();
	for(int j=debut_j;j<=fin_j;j++)
		for(int i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			if(RegionEllipse.PtInRegion(h))
			{
				l++;
			}
		}

	sprintf(ch,"nombre de points dans le masque:%d",l);
	SetCursor(LoadCursor(0,IDC_ARROW));  
	MessageBox(ch,"Masque");
}

void CTraitImageWnd::OnUpdateMasqueNbpixels(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);	
}


void CTraitImageWnd::OnMasqueHistogramme() 
{
	long TabHisto[256];
	CPoint h;
	
	SetCursor(LoadCursor(0,IDC_WAIT));
	AffectRegionMasque();
	Memorise();
	
	for (long i=0;i<256;i++) TabHisto[i] = 0;

	for(int j=debut_j;j<=fin_j;j++)
		for(int i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			if(RegionEllipse.PtInRegion(h))
			{
				TabHisto[Image[i+Larg*j]]++;
			}
		}

	CGraphWnd * pHisto = (CGraphWnd*) new CGraphWnd(TabHisto,256,"Nombre","Niveau");  
	if (!pHisto->Create("Histogramme Masque: "+Titre,GetMDIFrame())) return;
	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::OnUpdateMasqueHistogramme(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);
}

void CTraitImageWnd::OnMasqueEssai() 
{
	CPoint h;
	BOOL existe=FALSE;
	long temp,comp=0,comp2=0;
	debut_i+=23;
	debut_j+=23;
	fin_i-=23;
	fin_j-=23;
	
	Invalidate(FALSE);
	AffectRegionMasque();
	Memorise();
	SetCursor(LoadCursor(0,IDC_WAIT)); 
	OnTraitementGradientderob();
	for(int j=debut_j;j<=fin_j;j++)
		for(int i=debut_i;i<=fin_i;i++)
		{
			temp=0;
			comp=0;
			existe=FALSE;
			for(int li=0;li<8;li++)
			for(int co=0;co<8;co++)
			{
				h=Inv_Get_ij(i,j);

				if(RegionEllipse.PtInRegion(h))
				
				{
					if (((i+co)<=fin_i)||((li+j)<=fin_j))
					{
						h=Inv_Get_ij(i+co,j+li);

						if(RegionEllipse.PtInRegion(h))
						{	
							
							temp+=Image[i+co+Larg*(j+li)];
							existe=TRUE;
							comp++;
						}
					}
				}
			}
			if((existe)&&(temp/comp<36)) 
			{
				
				Image[i+j*Larg]=0;
				comp2++;
			}
			else Image[i+j*Larg]=255;
		}
	Invalidate(FALSE);
	if(comp2>200) MessageBox("défaut");
	else MessageBox("ok");
	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::OnUpdateMasqueEssai(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);
}

void CTraitImageWnd::OnMasqueEssai2() 
{
	CPoint h;
	BOOL existe=FALSE;
	long temp,comp=0,comp2=0;
	debut_i+=23;
	debut_j+=23;
	fin_i-=23;
	fin_j-=23;
	
	Invalidate(FALSE);
	AffectRegionMasque();
	Memorise();
	SetCursor(LoadCursor(0,IDC_WAIT)); 
	OnTraitementGradientderob();
	for(int j=debut_j;j<=fin_j;j++)
		for(int i=debut_i;i<=fin_i;i++)
		{
			temp=0;
			comp=0;
			existe=FALSE;
			for(int li=0;li<8;li++)
			for(int co=0;co<8;co++)
			{
				
				if (((i+co)<=fin_i)||((li+j)<=fin_j))
					{
						h=Inv_Get_ij(i,j);

						if(RegionEllipse.PtInRegion(h))
						{	
							temp+=Image[i+co+Larg*(j+li)];
							existe=TRUE;
							comp++;
						}
					}
			}
			if((existe)&&(temp/comp<36)) 
			{
				
				Image[i+j*Larg]=0;
				comp2++;
			}
			else Image[i+j*Larg]=255;
		}
	Invalidate(FALSE);
	if(comp2>200) MessageBox("défaut");
	else MessageBox("ok");
	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::OnUpdateMasqueEssai2(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);	
}

void CTraitImageWnd::OnMasqueEssai3() 
{
	long i,j,l,c;
	CPoint h;
	BOOL existe=FALSE;
	long complisse=0,compbourgeon=0,temp,comp=0,moy=0,nbpoint=0;
	char ch[50];
	unsigned char *tampon;

	tampon=new unsigned char [Haut*Larg];
	SetCursor(LoadCursor(0,IDC_WAIT)); 
	for(i=0;i<Haut*Larg;i++)tampon[i]=Image[i];

	for(j=0;j<=Haut-8;j++)
		for(i=0;i<=Larg-8;i++)
	{
		temp=0;
		for(l=0;l<8;l++)
			for(c=0;c<8;c++)
		{
			temp+=tampon[i+c+(j+l)*Larg];
		}
		tampon[i+j*Larg]=(unsigned char)(temp/64);

	}

	AffectRegionMasque();
	Memorise();
	
	OnTraitementGradientderob();
	OnTraitementLuminancemoyenne();
	debut_i+=19;
	debut_j+=19;
	fin_i-=27;
	fin_j-=27;
	AffectRegionMasque();

	for( j=debut_j;j<=fin_j;j++)
		for( i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			if(RegionEllipse.PtInRegion(h)) 
			{				
				moy+=tampon[i+j*Larg];
				nbpoint++;
			}
		}
	moy/=nbpoint;

	for( j=debut_j;j<=fin_j;j++)
		for( i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			l=i+j*Larg;
			if(RegionEllipse.PtInRegion(h)) 
			{
				if((tampon[l]>(moy+30))&&(Image[l]<40)) 
				{
					Image[l]=128;
					compbourgeon++;
				}
				else
				{ if(Image[l]<34)
					{
						Image[l]=0;
						complisse++;
					}
					else Image[l]=255;
				}
			}
			
		}
	Invalidate(FALSE);
	sprintf(ch,"lisses=%d\nbourgeon=%d",complisse,compbourgeon);
	MessageBox(ch);
	delete tampon;
	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::OnUpdateMasqueEssai3(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(Masque);
}

void CTraitImageWnd::OnMasqueEssai4() 
{
	long i,j,l,c;
	CPoint h,h1,h2;
	BOOL existe=FALSE;
	long complisse=0,compbourgeon=0,compmanque=0,temp,comp=0,moy=0,nbpoint=0;
	char ch[50];
	unsigned char *tampon;
	CRgn petitellipse;

	ArreteLeContinu();
	OnMasqueRechercher();

	tampon=new unsigned char [Haut*Larg];
	SetCursor(LoadCursor(0,IDC_WAIT)); 
	for(i=0;i<Haut*Larg;i++)tampon[i]=Image[i];//stockage de l'image dans tampon


	//calcul de la luminance moyenne dans tampon
	for(j=0;j<=Haut-8;j++)
		for(i=0;i<=Larg-8;i++)
	{
		temp=0;
		for(l=0;l<8;l++)
			for(c=0;c<8;c++)
		{
			temp+=tampon[i+c+(j+l)*Larg];
		}
		tampon[i+j*Larg]=(unsigned char)(temp/64);

	}

	debut_i-=4;
	debut_j-=4;
	fin_i+=4;
	fin_j+=4;
	AffectRegionMasque();
	Memorise();
	
	OnTraitementGradientderob();//gradient dans image
	OnTraitementLuminancemoyenne();//luminance moyenne dans image
	
	h1 = Inv_Get_ij((int)debut_i+25,(int)debut_j+22); 
	h2 = Inv_Get_ij((int)fin_i-32,(int)fin_j-31);
	petitellipse.CreateEllipticRgnIndirect(CRect(h1,h2));//petite ellipse
	debut_i+=7;
	debut_j+=5;
	fin_i-=15;
	fin_j-=13;

	AffectRegionMasque();

	//moyenne de luminance dans tampon
	for( j=debut_j;j<=fin_j;j++)
		for( i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			if((petitellipse.PtInRegion(h))) 
			{				
				moy+=tampon[i+j*Larg];
				nbpoint++;
			}
		}
	moy/=nbpoint;
	
	//debut de détection
	for( j=debut_j;j<=fin_j;j++)
		for( i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			l=i+j*Larg;
			if((petitellipse.PtInRegion(h))) //dans la petite ellipse
			{
				if((tampon[l]<moy-52)/*&&(Image[l]>33)*/)//manque dans petite ellipse
				{
					compmanque++;
					Image[l]=160;
				}
				else
					if((tampon[l]>(moy+30))&&(Image[l]<35))//bourgeon dans petite ellipse 
				{
					Image[l]=255;
					compbourgeon++;
				}
				else if(Image[l]<30)//lisse dans petite ellipse
					{
						Image[l]=100;
						complisse++;
					}
				 
					else Image[l]=0;
				
			}
			else if(RegionEllipse.PtInRegion(h))//dans couronne
			{
				if((tampon[l]>(moy+32))&&(Image[l]<38))//bourgeon dans couronne
				{
					Image[l]=255;
					compbourgeon++;
				}
				else if((tampon[l]<moy-58)||(Image[l]>80))//manque dans couronne
				{
					compmanque++;
					Image[l]=160;
				}
				else Image[l]=60;
			}

			
		}
	Invalidate(FALSE);
	sprintf(ch,"lisses=%d\nbourgeon=%d\nmanque=%d",complisse,compbourgeon,compmanque);
	MessageBox(ch);
	delete tampon;
	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::OnMasqueEssai5() 
{
	long j,i,l;
	char ch[80];
	unsigned __int32 moy;
	int large,haute;
	CPoint h,h1,h2;
	BOOL existe=FALSE;
	long complisse=0,compbourgeon=0,compmanque=0,comp=0,nbpoint=0;
	unsigned char *tampon;//,*tampongrad;
	unsigned __int16 *cufen,*llufen,*lgradfen,lumi,grad;
	CImageWnd * essai,*essai2;
	CRgn petitellipse,grandellipse;

	ArreteLeContinu();	
	OnMasqueRechercher();
	SetCursor(LoadCursor(0,IDC_WAIT)); 
	debut_i+=2;
	debut_j-=4;
	fin_i=debut_i+271;//4*((fin_i-debut_i)/4)+debut_i+11;
	fin_j=debut_j+271;//4*((fin_j-debut_j)/4)+debut_j+11;
	

	
	h1 = Inv_Get_ij((int)debut_i+10,(int)debut_j+12); 
	h2 = Inv_Get_ij((int)fin_i-13,(int)fin_j-13);
	grandellipse.CreateEllipticRgnIndirect(CRect(h1,h2));
	h1 = Inv_Get_ij((int)debut_i+30,(int)debut_j+32); 
	h2 = Inv_Get_ij((int)fin_i-33,(int)fin_j-33);
	petitellipse.CreateEllipticRgnIndirect(CRect(h1,h2));
	large=fin_i-debut_i+1;
	haute=fin_j-debut_j+1;
	
	tampon=new unsigned char [haute*large];
	//tampongrad=new unsigned char  [haute*large];
	
	cufen=new unsigned __int16[large*haute];
	llufen=new unsigned __int16[large*haute];
	lgradfen=new unsigned __int16[large*haute];
	essai=new CImageWnd(large,haute);
	essai2=new CImageWnd(large,haute);
	for(i=0;i<large*haute;i++) 
	{
		tampon[i]/*=tampongrad[i]*/=Image[i%large+debut_i+(debut_j+i/large)*Larg];
	}

	//MessageBox("debut");
	//for(i=0;i<1000;i++)
	{
		moy=FenetreGlissanteRapideMoy(tampon,cufen,llufen,large,haute,8,8);
		GradientRapide((unsigned __int32*)tampon/*tampongrad*/,large,haute);
		FenetreGlissanteRapide(tampon/*tampongrad*/,cufen,lgradfen,large,haute,8,8);
		
		
	}
	//MessageBox("fin");
	//moy=10560;
	
	for( j=0;j<haute;j++)
		for( i=0;i<large;i++)
		{
			h=Inv_Get_ij(i+debut_i,j+debut_j);
			l=i+j*large;
			lumi=llufen[l];
			grad=lgradfen[l];
			if((petitellipse.PtInRegion(h))) //dans la petite ellipse
			{
				if((lumi<moy)/*&&(grad>1050)*/)//manque dans petite ellipse
				{
					compmanque++;
					Image[i+debut_i+Larg*(debut_j+j)]=160;
				}
				else
					if((lumi>(moy+5625))&&(grad<1075))//bourgeon dans petite ellipse 
				{
					Image[i+debut_i+Larg*(debut_j+j)]=255;
					compbourgeon++;
				}
				else if(grad<1000)//lisse dans petite ellipse
					{
						Image[i+debut_i+Larg*(debut_j+j)]=100;
						complisse++;
					}
				 
					else Image[i+debut_i+Larg*(debut_j+j)]=0;
				
			}
			else if(grandellipse.PtInRegion(h))//dans couronne
			{
				if((lumi>(moy+5760))&&(grad<1216))//bourgeon dans couronne
				{
					Image[i+debut_i+Larg*(debut_j+j)]=255;
					compbourgeon++;
				}
				else if((lumi<moy-150)||(grad>2700))//manque dans couronne
				{
					compmanque++;
					Image[i+debut_i+Larg*(debut_j+j)]=160;
				}
				else Image[i+debut_i+Larg*(debut_j+j)]=60;
			}
			else Image[i+debut_i+Larg*(debut_j+j)]=0;

			
		}
	Invalidate(FALSE);
	sprintf(ch,"lisses=%d\nbourgeon=%d\nmanque=%d\n\nlargeur=%d\nhauteur=%d\nmoy=%d",complisse,compbourgeon,compmanque,large,haute,moy);
	MessageBox(ch);

	for(i=0;i<large*haute;i++) 
	{
		essai->Image[i]=(unsigned char) ((unsigned __int16) (llufen[i]/64));
		essai2->Image[i]=(unsigned char) ((unsigned __int16) (lgradfen[i]/64));
	}

	essai->Create("lumi moy",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,NULL);
	essai2->Create("grad moy",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,NULL);
	essai->ShowWindow(SW_SHOW);
	essai2->ShowWindow(SW_SHOW);
	delete tampon;
	delete cufen;
	delete llufen;
	delete lgradfen;
	SetCursor(LoadCursor(0,IDC_ARROW));
}
/*
void CTraitImageWnd::OnMasqueEssai6() 
{
	long j,i,k,l;
	char ch[80];
	unsigned __int32 moy;
	long complisse,compbourgeon,compmanque;
	unsigned __int16 lumi,grad;
	BOOL detect1,detect2;
	unsigned char *tampon2;
	tampon2=new unsigned char[272*272];
	ArreteLeContinu();	
	SetCursor(LoadCursor(0,IDC_WAIT)); 
	
	MessageBox("debut");
//	for(int q=0;q<500;q++)
	{
	i=j=k=0;
	detect1=detect2=FALSE;
	while(!detect1)
	{
		if (Image[i/*%248+136+Larg*(i/248)*//*]>40) 
		{
			debut_j=i/Larg;///248;
			detect1=TRUE;
		}
		i++;
	}
	i=debut_j+120;
	k=0;
	while(!detect2)
	{
		if (Image[j/30+(i+j%30)*Larg]>40) 
		{
			debut_i=j/30;
			detect2=TRUE;
		}
		j++;
	
	}

	debut_i+=2;
	debut_j-=4;

	Inittampon=Image+debut_i+(Larg*debut_j);
	
	moy=FenetreGlissanteRapideMoy2(Inittampon,Initcufen,Initllufen,272,272,12,12);
	GradientRapide2((unsigned __int32*)Inittampon,(unsigned __int32*)tampon2,272,272);
	FenetreGlissanteRapide(tampon2,Initcufen,Initlgradfen,272,272,12,12);

	complisse=compbourgeon=compmanque=0;
//	for( j=0;j<272;j++)
		for( l=0;l<272*272;l++)
		{
			//l=i+j*272;
			lumi=Initllufen[l];
			grad=Initlgradfen[l];
			if(Masquergn[l]==1) //dans la petite ellipse
			{
				if((lumi<moy))//&&(grad>1050))//manque dans petite ellipse
				{
					compmanque++;
					//Image[i+debut_i+Larg*(debut_j+j)]=160;
				}
				else
					if((lumi>(moy+5625))&&(grad<1075))//bourgeon dans petite ellipse 
				{
					//Image[i+debut_i+Larg*(debut_j+j)]=255;
					compbourgeon++;
				}
				else if(grad<1000)//lisse dans petite ellipse
					{
						//Image[i+debut_i+Larg*(debut_j+j)]=100;
						complisse++;
					}
				 
					//else Image[i+debut_i+Larg*(debut_j+j)]=0;
				
			}
			else if(Masquergn[l]==2)//dans couronne
			{
				if((lumi>(moy+5760))&&(grad<1216))//bourgeon dans couronne
				{
					//Image[i+debut_i+Larg*(debut_j+j)]=255;
					compbourgeon++;
				}
				else if((lumi<moy-150)||(grad>2700))//manque dans couronne
				{
					compmanque++;
					//Image[i+debut_i+Larg*(debut_j+j)]=160;
				}
				//else Image[i+debut_i+Larg*(debut_j+j)]=60;
			}
			//else Image[i+debut_i+Larg*(debut_j+j)]=0;

			
	}
}
	Invalidate(FALSE);
	sprintf(ch,"lisses=%d\nbourgeon=%d\nmanque=%d\n\nmoy=%d",complisse,compbourgeon,compmanque,moy);
	MessageBox(ch);	
	CImageWnd *essai,*essai2;
	essai=new CImageWnd(272,272);
	essai2=new CImageWnd(272,272);
	for(i=0;i<272*272;i++) 
	{
		essai->Image[i]=(unsigned char) ((unsigned __int16) (Initllufen[i]/144));
		essai2->Image[i]=(unsigned char) ((unsigned __int16) (Initlgradfen[i]/144));
	}

	essai->Create("lumi moy",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,NULL);
	essai2->Create("grad moy",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,NULL);
	essai->ShowWindow(SW_SHOW);
	essai2->ShowWindow(SW_SHOW);
	delete tampon2;
	SetCursor(LoadCursor(0,IDC_ARROW));
}
*/
void CTraitImageWnd::OnMasqueDcoupage() 
{
	long n;
	double ecarttype[16],sygma[16],moy;
	unsigned char temp;
	CString ch0;
	FILE *Fichier;
	CPoint h,sup[3],inf[3];
	long nbPoint[15],lumPoint[15],l;
	double lumMoy[15];
	char ch [300];
	CSize difference;
	CRect rect,rect2,rect3;
	CRgn trianglesup,triangleinf,rgntemp,rgntemp1,rgntemp2,rgntemp3,ellipse1,rgn20,rgn21,
		rgn30,rgn31,rgn32,rgn33,rgn40,rgn41,rgn42,rgn43,rgn44,rgn45,rgn46,rgn47;


//////////////////************////////////////**************/////////////
	///////////		recherche et affectation des différentes zones
//////////////////************////////////////**************/////////////


	SetCursor(LoadCursor(0,IDC_WAIT));
	AffectRegionMasque();//trouve la grande ellipse

	RegionEllipse.GetRgnBox(rect);
	sup[0]=rect.CenterPoint();
	sup[1]=rect.TopLeft();
	sup[2].y=rect.TopLeft().y;
	sup[2].x=rect.BottomRight().x;
	inf[0]=rect.CenterPoint();
	inf[1]=rect.BottomRight();
	inf[2].x=rect.TopLeft().x;
	inf[2].y=rect.BottomRight().y;
	trianglesup.CreatePolygonRgn( sup, 3, ALTERNATE );
	triangleinf.CreatePolygonRgn( inf, 3, ALTERNATE );


	rgn40.CreateEllipticRgnIndirect(rect);
	rgn41.CreateEllipticRgnIndirect(rect);
	rgn42.CreateEllipticRgnIndirect(rect);
	rgn43.CreateEllipticRgnIndirect(rect);
	rgn44.CreateEllipticRgnIndirect(rect);
	rgn45.CreateEllipticRgnIndirect(rect);
	rgn46.CreateEllipticRgnIndirect(rect);
	rgn47.CreateEllipticRgnIndirect(rect);

	//calcul la 3eme ellipse
	rect2=rect;
	difference.cx=(long)(rect.Width()*0.15843497);
	difference.cy=(long)(rect.Height()*0.15843497);
	rect2.DeflateRect(difference);
	rgn30.CreateEllipticRgnIndirect(rect2);
	rgn31.CreateEllipticRgnIndirect(rect2);
	rgn32.CreateEllipticRgnIndirect(rect2);
	rgn33.CreateEllipticRgnIndirect(rect2);

	//calcul la 2eme ellipse
	rect2=rect;
	difference.cx=(long)(rect.Width()*0.2763932);
	difference.cy=(long)(rect.Height()*0.2763932);
	rect2.DeflateRect(difference);
	rgn20.CreateEllipticRgnIndirect(rect2);
	rgn21.CreateEllipticRgnIndirect(rect2);

	//calcul l'ellipse centrale
	rect2=rect;
	difference.cx=(long)(rect.Width()*0.3709005);
	difference.cy=(long)(rect.Height()*0.3709005);
	rect2.DeflateRect(difference);
	ellipse1.CreateEllipticRgnIndirect(rect2);

	//calcul les disques
	rgn40.CombineRgn(&rgn40,&rgn30,RGN_DIFF);
	rgn41.CombineRgn(&rgn41,&rgn30,RGN_DIFF);
	rgn42.CombineRgn(&rgn42,&rgn30,RGN_DIFF);
	rgn43.CombineRgn(&rgn43,&rgn30,RGN_DIFF);
	rgn44.CombineRgn(&rgn44,&rgn30,RGN_DIFF);
	rgn45.CombineRgn(&rgn45,&rgn30,RGN_DIFF);
	rgn46.CombineRgn(&rgn46,&rgn30,RGN_DIFF);
	rgn47.CombineRgn(&rgn47,&rgn30,RGN_DIFF);

	rgn30.CombineRgn( &rgn30, &rgn20, RGN_DIFF );
	rgn31.CombineRgn( &rgn31, &rgn20, RGN_DIFF );
	rgn32.CombineRgn( &rgn32, &rgn20, RGN_DIFF );
	rgn33.CombineRgn( &rgn33, &rgn20, RGN_DIFF );
	rgn20.CombineRgn( &rgn20, &ellipse1, RGN_DIFF );
	rgn21.CombineRgn( &rgn21, &ellipse1, RGN_DIFF );;

	//calcul les parties sup et inf
	rect3=rect2=rect;
	rect2.OffsetRect(0,int(rect.Height()/2-0.5));
	rect3.OffsetRect(0,-rect.Height()/2);
	rgntemp.CreateRectRgnIndirect(rect2);
	rgntemp1.CreateRectRgnIndirect(rect3);

	rgn21.CombineRgn(&rgn21,&rgntemp,RGN_DIFF);
	rgn32.CombineRgn(&rgn32,&rgntemp,RGN_DIFF);
	rgn33.CombineRgn(&rgn33,&rgntemp,RGN_DIFF);
	rgn44.CombineRgn(&rgn44,&rgntemp,RGN_DIFF);
	rgn45.CombineRgn(&rgn45,&rgntemp,RGN_DIFF);
	rgn46.CombineRgn(&rgn46,&rgntemp,RGN_DIFF);
	rgn47.CombineRgn(&rgn47,&rgntemp,RGN_DIFF);


	rgn20.CombineRgn(&rgn20,&rgntemp1,RGN_DIFF);
	rgn30.CombineRgn(&rgn30,&rgntemp1,RGN_DIFF);
	rgn31.CombineRgn(&rgn31,&rgntemp1,RGN_DIFF);
	rgn40.CombineRgn(&rgn40,&rgntemp1,RGN_DIFF);
	rgn41.CombineRgn(&rgn41,&rgntemp1,RGN_DIFF);
	rgn42.CombineRgn(&rgn42,&rgntemp1,RGN_DIFF);
	rgn43.CombineRgn(&rgn43,&rgntemp1,RGN_DIFF);
	
	//calcul les parties droites et gauches
	rect3=rect2=rect;
	rect2.OffsetRect(int(rect.Width()/2-0.5),0);
	rect3.OffsetRect(-rect.Width()/2,0);
	rgntemp2.CreateRectRgnIndirect(rect2);
	rgntemp3.CreateRectRgnIndirect(rect3);

	rgn31.CombineRgn(&rgn31,&rgntemp2,RGN_DIFF);
	rgn32.CombineRgn(&rgn32,&rgntemp2,RGN_DIFF);
	rgn42.CombineRgn(&rgn42,&rgntemp2,RGN_DIFF);
	rgn43.CombineRgn(&rgn43,&rgntemp2,RGN_DIFF);
	rgn44.CombineRgn(&rgn44,&rgntemp2,RGN_DIFF);
	rgn45.CombineRgn(&rgn45,&rgntemp2,RGN_DIFF);
	
	rgn30.CombineRgn(&rgn30,&rgntemp3,RGN_DIFF);
	rgn33.CombineRgn(&rgn33,&rgntemp3,RGN_DIFF);
	rgn40.CombineRgn(&rgn40,&rgntemp3,RGN_DIFF);
	rgn41.CombineRgn(&rgn41,&rgntemp3,RGN_DIFF);
	rgn46.CombineRgn(&rgn46,&rgntemp3,RGN_DIFF);
	rgn47.CombineRgn(&rgn47,&rgntemp3,RGN_DIFF);

	//soustraction des triangles
	rgn40.CombineRgn(&rgn40,&triangleinf,RGN_DIFF);
	rgn42.CombineRgn(&rgn42,&triangleinf,RGN_AND);
	rgn41.CombineRgn(&rgn41,&triangleinf,RGN_AND);
	rgn43.CombineRgn(&rgn43,&triangleinf,RGN_DIFF);

	rgn44.CombineRgn(&rgn44,&trianglesup,RGN_DIFF);
	rgn46.CombineRgn(&rgn46,&trianglesup,RGN_AND);
	rgn45.CombineRgn(&rgn45,&trianglesup,RGN_AND);
	rgn47.CombineRgn(&rgn47,&trianglesup,RGN_DIFF);

	trianglesup. DeleteObject( );
	triangleinf. DeleteObject( );
//////////////////************////////////////**************/////////////
	///////////	calcul dans les zones
//////////////////************////////////////**************/////////////
	int i,j;
	for( i=0;i<15;i++)
		{
			nbPoint[i]=lumPoint[i]=0;
			ecarttype[i]=sygma[i]=0.0;
		}
	//balayage rect
	for( j=debut_j;j<=fin_j;j++)
		for( i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			l=i+j*Larg;
			temp=Image[l];
			if(rgn47.PtInRegion(h))	{nbPoint[14]++;lumPoint[14]+=temp;Image[l]=255;}
			if(rgn46.PtInRegion(h))	{nbPoint[13]++;lumPoint[13]+=temp;Image[l]=200;}
			if(rgn45.PtInRegion(h))	{nbPoint[12]++;lumPoint[12]+=temp;Image[l]=175;}
			if(rgn44.PtInRegion(h))	{nbPoint[11]++;lumPoint[11]+=temp;Image[l]=150;}
			if(rgn43.PtInRegion(h))	{nbPoint[10]++;lumPoint[10]+=temp;Image[l]=125;}
			if(rgn42.PtInRegion(h))	{nbPoint[9]++;lumPoint[9]+=temp;Image[l]=100;}
			if(rgn41.PtInRegion(h))	{nbPoint[8]++;lumPoint[8]+=temp;Image[l]=75;}
			if(rgn40.PtInRegion(h))	{nbPoint[7]++;lumPoint[7]+=temp;Image[l]=0;}
			if(rgn33.PtInRegion(h))	{nbPoint[6]++;lumPoint[6]+=temp;Image[l]=100;}
			if(rgn32.PtInRegion(h))	{nbPoint[5]++;lumPoint[5]+=temp;Image[l]=0;}
			if(rgn31.PtInRegion(h))	{nbPoint[4]++;lumPoint[4]+=temp;Image[l]=175;}
			if(rgn30.PtInRegion(h))	{nbPoint[3]++;lumPoint[3]+=temp;Image[l]=200;}
			if(rgn21.PtInRegion(h))	{nbPoint[2]++;lumPoint[2]+=temp;Image[l]=200;}
			if(rgn20.PtInRegion(h))	{nbPoint[1]++;lumPoint[1]+=temp;Image[l]=100;}
			if(ellipse1.PtInRegion(h))	{nbPoint[0]++;lumPoint[0]+=temp;Image[l]=0;}
		}
	Invalidate(FALSE);

	sprintf(ch,
			"centre=%d\ndisque20=%d\ndisque21=%d\ndisque30=%d\ndisque31=%d\ndisque32=%d\ndisque33=%d\ndisque40=%d\ndisque41=%d\ndisque42=%d\ndisque43=%d\ndisque44=%d\ndisque45=%d\ndisque46=%d\ndisque47=%d",
			nbPoint[0],nbPoint[1],nbPoint[2],nbPoint[3],nbPoint[4],nbPoint[5],nbPoint[6],nbPoint[7],nbPoint[8],
			nbPoint[9],nbPoint[10],nbPoint[11],nbPoint[12],nbPoint[13],nbPoint[14]);
	MessageBox(ch,"Decoupage");

	OnRestaureImage();

	for(i=0;i<15;i++) lumMoy[i]=(double) lumPoint[i]/nbPoint[i];

	for( j=debut_j;j<=fin_j;j++)
		for( i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			l=i+j*Larg;
			temp=Image[l];
			if(rgn47.PtInRegion(h))	{ecarttype[14]+=pow((temp-lumMoy[14]),2);}
			if(rgn46.PtInRegion(h))	{ecarttype[13]+=pow((temp-lumMoy[13]),2);}
			if(rgn45.PtInRegion(h))	{ecarttype[12]+=pow((temp-lumMoy[12]),2);}
			if(rgn44.PtInRegion(h))	{ecarttype[11]+=pow((temp-lumMoy[11]),2);}
			if(rgn43.PtInRegion(h))	{ecarttype[10]+=pow((temp-lumMoy[10]),2);}
			if(rgn42.PtInRegion(h))	{ecarttype[9]+=pow((temp-lumMoy[9]),2);}
			if(rgn41.PtInRegion(h))	{ecarttype[8]+=pow((temp-lumMoy[8]),2);}
			if(rgn40.PtInRegion(h))	{ecarttype[7]+=pow((temp-lumMoy[7]),2);}
			
			if(rgn33.PtInRegion(h))	{ecarttype[6]+=pow((temp-lumMoy[6]),2);}
			if(rgn32.PtInRegion(h))	{ecarttype[5]+=pow((temp-lumMoy[5]),2);}
			if(rgn31.PtInRegion(h))	{ecarttype[4]+=pow((temp-lumMoy[4]),2);}
			if(rgn30.PtInRegion(h))	{ecarttype[3]+=pow((temp-lumMoy[3]),2);}
			if(rgn21.PtInRegion(h))	{ecarttype[2]+=pow((temp-lumMoy[2]),2);}
			if(rgn20.PtInRegion(h))	{ecarttype[1]+=pow((temp-lumMoy[1]),2);}
			if(ellipse1.PtInRegion(h))	{ecarttype[0]+=pow((temp-lumMoy[0]),2);}
		}

	///////////////////////////////////ecarttype des zones
	for(i=0;i<15;i++)
	{
		ecarttype[i]=ecarttype[i]/nbPoint[i];
		ecarttype[i]=sqrt(ecarttype[i]);
	}

	for( j=debut_j;j<=fin_j;j++)
		for( i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			l=i+j*Larg;
			temp=Image[l];
			if(rgn47.PtInRegion(h))	
			{
				if((temp<(lumMoy[14]+ecarttype[14]/2))&&(temp>(lumMoy[14]-ecarttype[14]/2)))
					sygma[14]++;
				Image[l]=(unsigned char)lumMoy[14];
			}
			if(rgn46.PtInRegion(h))	
			{
				if((temp<(lumMoy[13]+ecarttype[13]/2))&&(temp>(lumMoy[13]-ecarttype[13]/2)))
					sygma[13]++;
				Image[l]=(unsigned char)lumMoy[13];
			}
			if(rgn45.PtInRegion(h))	
			{
				if((temp<(lumMoy[12]+ecarttype[12]/2))&&(temp>(lumMoy[12]-ecarttype[12]/2)))
					sygma[12]++;
				Image[l]=(unsigned char)lumMoy[12];
			}
			if(rgn44.PtInRegion(h))	
			{
				if((temp<(lumMoy[11]+ecarttype[11]/2))&&(temp>(lumMoy[11]-ecarttype[11]/2)))
					sygma[11]++;
				Image[l]=(unsigned char)lumMoy[11];
			}
			if(rgn43.PtInRegion(h))	
			{
				if((temp<(lumMoy[10]+ecarttype[10]/2))&&(temp>(lumMoy[10]-ecarttype[10]/2)))
					sygma[10]++;
				Image[l]=(unsigned char)lumMoy[10];
			}
			if(rgn42.PtInRegion(h))
			{
				if((temp<(lumMoy[9]+ecarttype[9]/2))&&(temp>(lumMoy[9]-ecarttype[9]/2)))
					sygma[9]++;
				Image[l]=(unsigned char)lumMoy[9];
			}
			if(rgn41.PtInRegion(h))	
			{
				if((temp<(lumMoy[8]+ecarttype[8]/2))&&(temp>(lumMoy[8]-ecarttype[8]/2)))
					sygma[8]++;
				Image[l]=(unsigned char)lumMoy[8];
			}
			if(rgn40.PtInRegion(h))
			{
				if((temp<(lumMoy[7]+ecarttype[7]/2))&&(temp>(lumMoy[7]-ecarttype[7]/2)))
					sygma[7]++;
				Image[l]=(unsigned char)lumMoy[7];
			}
			
			if(rgn33.PtInRegion(h))
			{
				if((temp<(lumMoy[6]+ecarttype[6]/2))&&(temp>(lumMoy[6]-ecarttype[6]/2)))
					sygma[6]++;
				Image[l]=(unsigned char)lumMoy[6];
			}
			if(rgn32.PtInRegion(h))
			{
				if((temp<(lumMoy[5]+ecarttype[5]/2))&&(temp>(lumMoy[5]-ecarttype[5]/2)))
					sygma[5]++;
				Image[l]=(unsigned char)lumMoy[5];
			}
			if(rgn31.PtInRegion(h))
			{
				if((temp<(lumMoy[4]+ecarttype[4]/2))&&(temp>(lumMoy[4]-ecarttype[4]/2)))
					sygma[4]++;
				Image[l]=(unsigned char)lumMoy[4];
			}
			if(rgn30.PtInRegion(h))
			{
				if((temp<(lumMoy[3]+ecarttype[3]/2))&&(temp>(lumMoy[3]-ecarttype[3]/2)))
					sygma[3]++;
				Image[l]=(unsigned char)lumMoy[3];
			}
			if(rgn21.PtInRegion(h))
			{
				if((temp<(lumMoy[2]+ecarttype[2]/2))&&(temp>(lumMoy[2]-ecarttype[2]/2)))
					sygma[2]++;
				Image[l]=(unsigned char)lumMoy[2];
			}
			if(rgn20.PtInRegion(h))
			{
				if((temp<(lumMoy[1]+ecarttype[1]/2))&&(temp>(lumMoy[1]-ecarttype[1]/2)))
					sygma[1]++;
				Image[l]=(unsigned char)lumMoy[1];
			}
			if(ellipse1.PtInRegion(h))
			{
				if((temp<(lumMoy[0]+ecarttype[0]/2))&&(temp>(lumMoy[0]-ecarttype[0]/2)))
					sygma[0]++;
				Image[l]=(unsigned char)lumMoy[0];
			}
		}


	OnCouleur();
	Invalidate(FALSE);
	

	for(i=0;i<15;i++)	sygma[i]=(100*sygma[i])/nbPoint[i];
	
	
	///////////////////////////////////moy ecarttype et % à 1sygma sur total
	moy=0;
	for(i=0;i<15;i++) moy+=lumMoy[i];
	moy/=15;

	ecarttype[15]=0;
	for(i=0;i<15;i++) ecarttype[15]+=pow((lumMoy[i]-moy),2);
	ecarttype[15]/=15;
	ecarttype[15]=sqrt(ecarttype[15]);

	sygma[15]=0;
	for(i=0;i<15;i++) if((lumMoy[i]<(moy+ecarttype[15]/2))&&
		(lumMoy[i]>(moy-ecarttype[15]/2))) sygma[15]++;

	sygma[15]=(100*sygma[15])/15;

	///////////////////////////////////////////////////////////////////////////
	sprintf(ch,
			"centre=%3.2f\ndisque20=%3.2f\ndisque21=%3.2f\ndisque30=%3.2f\ndisque31=%3.2f\ndisque32=%3.2f\ndisque33=%3.2f\ndisque40=%3.2f\ndisque41=%3.2f\ndisque42=%3.2f\ndisque43=%3.2f\ndisque44=%3.2f\ndisque45=%3.2f\ndisque46=%3.2f\ndisque47=%3.2f\n\nmoyenne=%3.2f\nécart-type=%2.2f,\n%% à1 sygma=%2.2f",
			lumMoy[0],lumMoy[1],lumMoy[2],lumMoy[3],lumMoy[4],lumMoy[5],lumMoy[6],lumMoy[7],lumMoy[8],
			lumMoy[9],lumMoy[10],lumMoy[11],lumMoy[12],lumMoy[13],lumMoy[14],moy,ecarttype[15],sygma[15]);
	MessageBox(ch,"Decoupage");

	SetCursor(LoadCursor(0,IDC_ARROW));
	
	CFileDialog FileData(FALSE,"*.xls",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"EXCEL (*.xls)|*.xls||",this);
	if (FileData.DoModal()==IDOK)
	{
		if (Fichier=fopen(FileData.GetPathName(),"w"))
		{
			SetCursor(LoadCursor(0,IDC_WAIT));
			fprintf(Fichier,Titre);
			fprintf(Fichier,"\nCoefficient d'homoégnéité\n");
			fprintf(Fichier,"zone;moyenne;ecart-type;%% à 1 sygma\n");

			fprintf(Fichier,"Totale;");
			sprintf(ch,"%3.2f",moy);
			ch0=ch;
			ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
			fprintf(Fichier,ch0+";");
			sprintf(ch,"%2.2f",ecarttype[15]);
			ch0=ch;
			ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
			fprintf(Fichier,ch0+";");

			sprintf(ch,"%2.2f",sygma[15]);
			ch0=ch;
			ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
			fprintf(Fichier,ch0+";\n");
			n=0;l=1;
			for(i=1;i<5;i++)
			{
				
				for(j=0;j<l;j++)
				{
					fprintf(Fichier,"%d",i);
					fprintf(Fichier,"%d;",j);
					sprintf(ch,"%3.2f",lumMoy[n]);
					ch0=ch;
					ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
					fprintf(Fichier,ch0+";");
					sprintf(ch,"%2.2f",ecarttype[n]);
					ch0=ch;
					ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
					fprintf(Fichier,ch0+";");

					sprintf(ch,"%2.2f",sygma[n++]);
					ch0=ch;
					ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
					fprintf(Fichier,ch0+";\n");
				}
				l*=2;


			}
			fclose(Fichier);
			SetCursor(LoadCursor(0,IDC_ARROW));
		}
	} 
}/////////decoupage

void CTraitImageWnd::OnCaracteclair() 
{
	decoupageZone(40);
}/////////OnCaracteclair

void CTraitImageWnd::decoupageZone(int m_seuil)
{
	long n;
	double ecarttype,sygma,moy;
	unsigned char temp;
	CString ch0;
	FILE *Fichier;
	CPoint h,sup[3],inf[3];
	long nbPoint[15],lumPoint[15],l;
	double lumMoy[15];
	char ch [300];
	CSize difference;
	CRect rect,rect2,rect3;
	CRgn trianglesup,triangleinf,rgntemp,rgntemp1,rgntemp2,rgntemp3,ellipse1,rgn20,rgn21,
		rgn30,rgn31,rgn32,rgn33,rgn40,rgn41,rgn42,rgn43,rgn44,rgn45,rgn46,rgn47;


//////////////////************////////////////**************/////////////
	///////////		recherche et affectation des différentes zones
//////////////////************////////////////**************/////////////
	OnMasqueRechercheAutom(m_seuil);


	SetCursor(LoadCursor(0,IDC_WAIT));
	AffectRegionMasque();//trouve la grande ellipse

	RegionEllipse.GetRgnBox(rect);
	sup[0]=rect.CenterPoint();
	sup[1]=rect.TopLeft();
	sup[2].y=rect.TopLeft().y;
	sup[2].x=rect.BottomRight().x;
	inf[0]=rect.CenterPoint();
	inf[1]=rect.BottomRight();
	inf[2].x=rect.TopLeft().x;
	inf[2].y=rect.BottomRight().y;
	trianglesup.CreatePolygonRgn( sup, 3, ALTERNATE );
	triangleinf.CreatePolygonRgn( inf, 3, ALTERNATE );


	rgn40.CreateEllipticRgnIndirect(rect);
	rgn41.CreateEllipticRgnIndirect(rect);
	rgn42.CreateEllipticRgnIndirect(rect);
	rgn43.CreateEllipticRgnIndirect(rect);
	rgn44.CreateEllipticRgnIndirect(rect);
	rgn45.CreateEllipticRgnIndirect(rect);
	rgn46.CreateEllipticRgnIndirect(rect);
	rgn47.CreateEllipticRgnIndirect(rect);

	//calcul la 3eme ellipse
	rect2=rect;
	difference.cx=(long)(rect.Width()*0.15843497);
	difference.cy=(long)(rect.Height()*0.15843497);
	rect2.DeflateRect(difference);
	rgn30.CreateEllipticRgnIndirect(rect2);
	rgn31.CreateEllipticRgnIndirect(rect2);
	rgn32.CreateEllipticRgnIndirect(rect2);
	rgn33.CreateEllipticRgnIndirect(rect2);

	//calcul la 2eme ellipse
	rect2=rect;
	difference.cx=(long)(rect.Width()*0.2763932);
	difference.cy=(long)(rect.Height()*0.2763932);
	rect2.DeflateRect(difference);
	rgn20.CreateEllipticRgnIndirect(rect2);
	rgn21.CreateEllipticRgnIndirect(rect2);

	//calcul l'ellipse centrale
	rect2=rect;
	difference.cx=(long)(rect.Width()*0.3709005);
	difference.cy=(long)(rect.Height()*0.3709005);
	rect2.DeflateRect(difference);
	ellipse1.CreateEllipticRgnIndirect(rect2);

	//calcul les disques
	rgn40.CombineRgn(&rgn40,&rgn30,RGN_DIFF);
	rgn41.CombineRgn(&rgn41,&rgn30,RGN_DIFF);
	rgn42.CombineRgn(&rgn42,&rgn30,RGN_DIFF);
	rgn43.CombineRgn(&rgn43,&rgn30,RGN_DIFF);
	rgn44.CombineRgn(&rgn44,&rgn30,RGN_DIFF);
	rgn45.CombineRgn(&rgn45,&rgn30,RGN_DIFF);
	rgn46.CombineRgn(&rgn46,&rgn30,RGN_DIFF);
	rgn47.CombineRgn(&rgn47,&rgn30,RGN_DIFF);

	rgn30.CombineRgn( &rgn30, &rgn20, RGN_DIFF );
	rgn31.CombineRgn( &rgn31, &rgn20, RGN_DIFF );
	rgn32.CombineRgn( &rgn32, &rgn20, RGN_DIFF );
	rgn33.CombineRgn( &rgn33, &rgn20, RGN_DIFF );
	rgn20.CombineRgn( &rgn20, &ellipse1, RGN_DIFF );
	rgn21.CombineRgn( &rgn21, &ellipse1, RGN_DIFF );;

	//calcul les parties sup et inf
	rect3=rect2=rect;
	rect2.OffsetRect(0,int(rect.Height()/2-0.5));
	rect3.OffsetRect(0,-rect.Height()/2);
	rgntemp.CreateRectRgnIndirect(rect2);
	rgntemp1.CreateRectRgnIndirect(rect3);

	rgn21.CombineRgn(&rgn21,&rgntemp,RGN_DIFF);
	rgn32.CombineRgn(&rgn32,&rgntemp,RGN_DIFF);
	rgn33.CombineRgn(&rgn33,&rgntemp,RGN_DIFF);
	rgn44.CombineRgn(&rgn44,&rgntemp,RGN_DIFF);
	rgn45.CombineRgn(&rgn45,&rgntemp,RGN_DIFF);
	rgn46.CombineRgn(&rgn46,&rgntemp,RGN_DIFF);
	rgn47.CombineRgn(&rgn47,&rgntemp,RGN_DIFF);


	rgn20.CombineRgn(&rgn20,&rgntemp1,RGN_DIFF);
	rgn30.CombineRgn(&rgn30,&rgntemp1,RGN_DIFF);
	rgn31.CombineRgn(&rgn31,&rgntemp1,RGN_DIFF);
	rgn40.CombineRgn(&rgn40,&rgntemp1,RGN_DIFF);
	rgn41.CombineRgn(&rgn41,&rgntemp1,RGN_DIFF);
	rgn42.CombineRgn(&rgn42,&rgntemp1,RGN_DIFF);
	rgn43.CombineRgn(&rgn43,&rgntemp1,RGN_DIFF);
	
	//calcul les parties droites et gauches
	rect3=rect2=rect;
	rect2.OffsetRect(int(rect.Width()/2-0.5),0);
	rect3.OffsetRect(-rect.Width()/2,0);
	rgntemp2.CreateRectRgnIndirect(rect2);
	rgntemp3.CreateRectRgnIndirect(rect3);

	rgn31.CombineRgn(&rgn31,&rgntemp2,RGN_DIFF);
	rgn32.CombineRgn(&rgn32,&rgntemp2,RGN_DIFF);
	rgn42.CombineRgn(&rgn42,&rgntemp2,RGN_DIFF);
	rgn43.CombineRgn(&rgn43,&rgntemp2,RGN_DIFF);
	rgn44.CombineRgn(&rgn44,&rgntemp2,RGN_DIFF);
	rgn45.CombineRgn(&rgn45,&rgntemp2,RGN_DIFF);
	
	rgn30.CombineRgn(&rgn30,&rgntemp3,RGN_DIFF);
	rgn33.CombineRgn(&rgn33,&rgntemp3,RGN_DIFF);
	rgn40.CombineRgn(&rgn40,&rgntemp3,RGN_DIFF);
	rgn41.CombineRgn(&rgn41,&rgntemp3,RGN_DIFF);
	rgn46.CombineRgn(&rgn46,&rgntemp3,RGN_DIFF);
	rgn47.CombineRgn(&rgn47,&rgntemp3,RGN_DIFF);

	//soustraction des triangles
	rgn40.CombineRgn(&rgn40,&triangleinf,RGN_DIFF);
	rgn42.CombineRgn(&rgn42,&triangleinf,RGN_AND);
	rgn41.CombineRgn(&rgn41,&triangleinf,RGN_AND);
	rgn43.CombineRgn(&rgn43,&triangleinf,RGN_DIFF);

	rgn44.CombineRgn(&rgn44,&trianglesup,RGN_DIFF);
	rgn46.CombineRgn(&rgn46,&trianglesup,RGN_AND);
	rgn45.CombineRgn(&rgn45,&trianglesup,RGN_AND);
	rgn47.CombineRgn(&rgn47,&trianglesup,RGN_DIFF);

	trianglesup. DeleteObject( );
	triangleinf. DeleteObject( );
//////////////////************////////////////**************/////////////
	///////////	calcul dans les zones
//////////////////************////////////////**************/////////////
	int i,j;
	for( i=0;i<15;i++)
		{
			nbPoint[i]=lumPoint[i]=0;
		}
	//balayage rect
	for( j=debut_j;j<=fin_j;j++)
		for( i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			l=i+j*Larg;
			temp=Image[l];
			if(rgn47.PtInRegion(h))	{nbPoint[14]++;lumPoint[14]+=temp;}
			if(rgn46.PtInRegion(h))	{nbPoint[13]++;lumPoint[13]+=temp;}
			if(rgn45.PtInRegion(h))	{nbPoint[12]++;lumPoint[12]+=temp;}
			if(rgn44.PtInRegion(h))	{nbPoint[11]++;lumPoint[11]+=temp;}
			if(rgn43.PtInRegion(h))	{nbPoint[10]++;lumPoint[10]+=temp;}
			if(rgn42.PtInRegion(h))	{nbPoint[9]++;lumPoint[9]+=temp;}
			if(rgn41.PtInRegion(h))	{nbPoint[8]++;lumPoint[8]+=temp;}
			if(rgn40.PtInRegion(h))	{nbPoint[7]++;lumPoint[7]+=temp;}
			if(rgn33.PtInRegion(h))	{nbPoint[6]++;lumPoint[6]+=temp;}
			if(rgn32.PtInRegion(h))	{nbPoint[5]++;lumPoint[5]+=temp;}
			if(rgn31.PtInRegion(h))	{nbPoint[4]++;lumPoint[4]+=temp;}
			if(rgn30.PtInRegion(h))	{nbPoint[3]++;lumPoint[3]+=temp;}
			if(rgn21.PtInRegion(h))	{nbPoint[2]++;lumPoint[2]+=temp;}
			if(rgn20.PtInRegion(h))	{nbPoint[1]++;lumPoint[1]+=temp;}
			if(ellipse1.PtInRegion(h))	{nbPoint[0]++;lumPoint[0]+=temp;}
		}
	for(i=0;i<15;i++) lumMoy[i]=(double) lumPoint[i]/nbPoint[i];
	for( j=debut_j;j<=fin_j;j++)
		for( i=debut_i;i<=fin_i;i++)
		{
			h=Inv_Get_ij(i,j);
			l=i+j*Larg;
			temp=Image[l];
			if(rgn47.PtInRegion(h))	
			{
				Image[l]=(unsigned char)lumMoy[14];
			}
			if(rgn46.PtInRegion(h))	
			{
				Image[l]=(unsigned char)lumMoy[13];
			}
			if(rgn45.PtInRegion(h))	
			{
				Image[l]=(unsigned char)lumMoy[12];
			}
			if(rgn44.PtInRegion(h))	
			{
				Image[l]=(unsigned char)lumMoy[11];
			}
			if(rgn43.PtInRegion(h))	
			{
				Image[l]=(unsigned char)lumMoy[10];
			}
			if(rgn42.PtInRegion(h))
			{
				Image[l]=(unsigned char)lumMoy[9];
			}
			if(rgn41.PtInRegion(h))	
			{
				Image[l]=(unsigned char)lumMoy[8];
			}
			if(rgn40.PtInRegion(h))
			{
				Image[l]=(unsigned char)lumMoy[7];
			}
			
			if(rgn33.PtInRegion(h))
			{
				Image[l]=(unsigned char)lumMoy[6];
			}
			if(rgn32.PtInRegion(h))
			{
				Image[l]=(unsigned char)lumMoy[5];
			}
			if(rgn31.PtInRegion(h))
			{
				Image[l]=(unsigned char)lumMoy[4];
			}
			if(rgn30.PtInRegion(h))
			{
				Image[l]=(unsigned char)lumMoy[3];
			}
			if(rgn21.PtInRegion(h))
			{
				Image[l]=(unsigned char)lumMoy[2];
			}
			if(rgn20.PtInRegion(h))
			{
				Image[l]=(unsigned char)lumMoy[1];
			}
			if(ellipse1.PtInRegion(h))
			{
				Image[l]=(unsigned char)lumMoy[0];
			}
		}


	if(TypePalette!=PALCOUL) OnCouleur();
	Invalidate(FALSE);
	

	
	
	///////////////////////////////////moy ecarttype et % à 1sygma sur total
	moy=0;
	for(i=0;i<15;i++) moy+=lumMoy[i];
	moy/=15;

	ecarttype=0;
	for(i=0;i<15;i++) ecarttype+=pow((lumMoy[i]-moy),2);
	ecarttype/=15;
	ecarttype=sqrt(ecarttype);

	sygma=0;
	for(i=0;i<15;i++) if((lumMoy[i]<(moy+ecarttype/2))&&
		(lumMoy[i]>(moy-ecarttype/2))) sygma++;

	sygma=(100*sygma)/15;

	sprintf(ch,"Moyenne Totale=%3.2f\nEcart Type Total=%3.2f\n%% à 1 ygma Total=%3.2f\n\nMoy centre=%3.2f\nMoy 20=%3.2f\nMoy 21=%3.2f\nMoy 30=%3.2f\nMoy 31=%3.2f\nMoy 32=%3.2f\nMoy 33=%3.2f\nMoy 40=%3.2f\nMoy 41=%3.2f\nMoy 42=%3.2f\nMoy 43=%3.2f\nMoy 44=%3.2f\nMoy 45=%3.2f\nMoy 46=%3.2f\nMoy 47=%3.2f",
		moy,ecarttype,sygma,lumMoy[0],lumMoy[1],lumMoy[2],lumMoy[3],lumMoy[4],lumMoy[5],lumMoy[6],lumMoy[7],lumMoy[8],
			lumMoy[9],lumMoy[10],lumMoy[11],lumMoy[12],lumMoy[13],lumMoy[14]);
	MessageBox(ch,"Decoupage");

	SetCursor(LoadCursor(0,IDC_ARROW));

	CFileDialog FileData(FALSE,"*.xls",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"EXCEL (*.xls)|*.xls||",this);
	if (FileData.DoModal()==IDOK)
	{
		if (Fichier=fopen(FileData.GetPathName(),"w"))
		{
			SetCursor(LoadCursor(0,IDC_WAIT));
			fprintf(Fichier,Titre);
			fprintf(Fichier,"\nCoefficient d'homoégnéité\n");
			fprintf(Fichier,"zone;moyenne;ecart-type;%% à 1 sygma\n");

			fprintf(Fichier,"Totale;");
			sprintf(ch,"%3.2f",moy);
			ch0=ch;
			ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
			fprintf(Fichier,ch0+";");
			sprintf(ch,"%2.2f",ecarttype);
			ch0=ch;
			ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
			fprintf(Fichier,ch0+";");

			sprintf(ch,"%2.2f",sygma);
			ch0=ch;
			ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
			fprintf(Fichier,ch0+";\n");
			n=0;l=1;
			for(i=1;i<5;i++)
			{
				
				for(j=0;j<l;j++)
				{
					fprintf(Fichier,"%d",i);
					fprintf(Fichier,"%d;",j);
					sprintf(ch,"%3.2f",lumMoy[n++]);
					ch0=ch;
					ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
					fprintf(Fichier,ch0+";\n");
				}
				l*=2;
			}
			fclose(Fichier);
			SetCursor(LoadCursor(0,IDC_ARROW));
		}
	} 

}

void CTraitImageWnd::OnMasqueCaractrisationglobal() 
{
	ArreteLeContinu();
	int i;
	unsigned char *tamp;
	unsigned __int16 *templ,*tempc;
	int large=Larg-Larg%4;
	int haute=Haut-Haut%4;
	CTraitImageWnd *gradMoy;
	CTraitImageWnd * lumMoy;
	long stop,l;
	templ=new unsigned __int16[large*haute];
	tempc=new unsigned __int16[large*haute];
	tamp=new unsigned char[large*haute];

	gradMoy=new CTraitImageWnd(large,haute);
	lumMoy=new CTraitImageWnd(large,haute);
	SetCursor(LoadCursor(0,IDC_WAIT));
	stop=6+6*Larg;
	/////////lum fen gliss
	FenetreGlissanteRapide2(Image,templ,tempc,large,haute,12,12);
	l=stop;
	for(i=0;i<large*haute-stop;i++) lumMoy->Image[i]=tempc[l++]/144;
	for(i=large*haute-stop;i<large*haute;i++) lumMoy->Image[i]=0;
	lumMoy->Create("Luminance Moy glis 12*12 de "+Titre, WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());
	lumMoy->decoupageZone(40);
	/////////grad fen gliss
	GradientRapide2((unsigned __int32*)Image,(unsigned __int32*)tamp,large,haute);
	FenetreGlissanteRapide(tamp,templ,tempc,large,haute,12,12);

	l=stop;
	for(i=0;i<large*(haute-1)-stop;i++) gradMoy->Image[i]=tempc[l++]/144;
	for(i=large*(haute-1)-stop;i<large*haute;i++) gradMoy->Image[i]=0;
	gradMoy->Create("Gradient Moy glis 12*12 de "+Titre, WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());
	gradMoy->decoupageZone(10);
	delete templ;
	delete tempc;
	delete tamp;
	
}

void CTraitImageWnd::OnUpdateMasqueDcoupage(CCmdUI* pCmdUI) 
{
//	pCmdUI->Enable(Masque);	
}

void CTraitImageWnd::OnCalculTexture() 
{
	CString ch0;
//	char ch[50];
	FILE * Fichier;
	float temp;
//	CProgressDlg Dlg;

	unsigned char *tamp;
	unsigned __int16 *templ,*tempc;
	int large=Larg-Larg%4;
	int haute=Haut-Haut%4;
	long stop,l,j,k,p;
	templ=new unsigned __int16[large*haute];
	tempc=new unsigned __int16[large*haute];
	tamp=new unsigned char[large*haute];
	stop=6+6*large;

	CFileDialog FileData(FALSE,"*.xls",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"EXCEL (*.xls)|*.xls||",this);
	if (FileData.DoModal()==IDOK)
	{
		if (Fichier=fopen(FileData.GetPathName(),"w"))
		{
	/*		Dlg.Create();
			sprintf(ch,"Texture");
			Dlg.SetWindowText(ch);
			Dlg.SetRange(1,NbMaxImageSeq);
			Dlg.SetStep(1);
	*/		SetCursor(LoadCursor(0,IDC_WAIT));
			fprintf(Fichier,Titre+" moyenne du gradient dans fen gliss 12*12 sur la surface de la cathode\n");
			fprintf(Fichier,"seuil:10\n");
			for (int i=0;i<NbMaxImageSeq;i++)
			{

				AfficheSequence(i,i+1);
				GradientRapide2((unsigned __int32*)Image,(unsigned __int32*)tamp,large,haute);
				FenetreGlissanteRapide((unsigned char*)tamp,templ,tempc,large,haute,12,12);
				l=stop;
				p=0;
				for(j=0;j<haute-6;j++)
				{
					p=j*Larg;
					for(k=0;k<large;k++)
					{
						Image[p]=tempc[l++]/144;
						p++;
					}
				}

				OnMasqueRechercheAutom(10);
				debut_i+=24;
				debut_j+=24;
				fin_i-=24;
				fin_j-=24;
				temp=OnMasqueLuminancemoyenneAutom();
				OnMasqueAnnuler();
		//		sprintf(ch,"Texture image n°%d =%3.2f",i+1,temp);
		//		Dlg.SetWindowText(ch);
			//	sprintf(ch,"%3.2f",temp);
				ch0.Format("%3.2f",temp);
				ch0.SetAt(ch0.Find('.'),',');
				fprintf(Fichier,"image %d\t",i);
				fprintf(Fichier,ch0+"\n");
			//	Dlg.StepIt();
			}
			fclose(Fichier);
		//	Dlg.DestroyWindow();
			SetCursor(LoadCursor(0,IDC_ARROW));
		}
	}   
	delete templ;
	delete tempc;
	delete tamp;   
}

void CTraitImageWnd::OnUpdateCalculTexture(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(SEQUENCE_CREE);
}
//////////////////////////////////////////////////////
///////////////////vue3d/////////////////////////////
/////////////////////////////////////////////////////
/*
void CTraitImageWnd::OnOutilsVue3d() 
{
	ArreteLeContinu();
	if(!IsWindow(PVue3d->GetSafeHwnd()))
	{
		PVue3d->Create("VUE 3D");
	}
	VUE3D=!VUE3D;
	if(VUE3D) PVue3d->ShowWindow(SW_SHOW);
	else PVue3d->ShowWindow(SW_HIDE);
	SetFocus();
	
}
*/
/*
void CTraitImageWnd::OnUpdateOutilsVue3d(CCmdUI* pCmdUI) 
{
	pCmdUI->SetRadio(VUE3D);
}
*/
//////////////////////////////////////////////////
//////////////////////////////////////////////////
void CTraitImageWnd::OnDestroy() 
{
	CImageWnd::OnDestroy();
	
//	if (DialPara!=NULL) delete DialPara;
//	if ((Nbapp!=NULL)) delete Nbapp;
	if (PParametre!=NULL) delete PParametre;
	//delete Buffer;
	for (int i=0;i<(int)NBPMAX;i++) 
	{ 
		if(IsWindow(VisuPar[i]->GetSafeHwnd()))//&&(VisuPar[i]!=NULL))//(VisuPar[i]->m_hWnd!=0xdddddddd&&(VisuPar[i]!=NULL)) 
		
			VisuPar[i]->DestroyWindow();
	}
	Visu=FALSE;
	if(DefImagAppr) ImageAppr->DestroyWindow();
	if(Visumes)  AffMesPol->DestroyWindow();
//	GetMDIFrame()->SendMessageToDescendants( uiMsg4,0,(LPARAM)this,TRUE,FALSE);
	delete Legende;
	
}


int CTraitImageWnd::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CImageWnd::OnCreate(lpCreateStruct) == -1)
		return -1;
	PositionMenuMasque=menu.GetMenuItemCount()-2;///////////attention au nombre de menus
	// TODO: Add your specialized creation code here
	if(TypeColor==NB&&TypeImage!=FLOAT32)
	{
		menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(6,MF_GRAYED|MF_BYPOSITION);
		menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(7,MF_GRAYED|MF_BYPOSITION);
		menu.GetSubMenu(PositionMenuMasque)->EnableMenuItem(12,MF_GRAYED|MF_BYPOSITION);
		menu.GetSubMenu(PositionMenuMasque-1)->EnableMenuItem(1,MF_GRAYED|MF_BYPOSITION);
	}
		return 0;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void CTraitImageWnd::OnApprentissageGnrerunprojet() 
{
	ArreteLeContinu();
	//cache les visupar
	CacherLesVisuPar();
	Visu=FALSE;
	//détruit l'image d'apprentissage s'il faut
	if(DefImagAppr) ImageAppr->DestroyWindow();

	//choisi un nom de projet par défaut
	if (NomFichProj=="") NomFichProj=FileName.Left(FileName.ReverseFind('.'))+".prj";

	//construction de la boite de dialogue et mise en place de l'apprentissage
	CDialProjAppr MonProjet(this,NomFichProj);
	if(MonProjet.DoModal()==IDOK)
	{
		SetCursor(LoadCursor(NULL,IDC_WAIT));

		NomFichProj=MonProjet.FileProj;
		NomImageAppr=MonProjet.ImageAppr;
		NomFichAppr=MonProjet.FileAppr;
		NomFichParam=MonProjet.FileParam;

		ImageAppr = new CImageClas(this->m_hWnd,NomImageAppr,Larg,Haut,BMP); 
		ImageAppr->Create(NomImageAppr.Right(NomImageAppr.GetLength()-NomImageAppr.ReverseFind('\\')-1)+"  Zones d'apprentissage", WS_MINIMIZE|WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,NULL);
		DefImagAppr=TRUE;
		this->SetFocus();

		DialPara.nbpara=0;
		ChargementDesParam(NomFichParam);
		AlloueParam();
		ParamCalcul=FALSE;
		ValideDecision=FALSE;
		ValideMenu=TRUE;
		menu.GetSubMenu(PositionMenuMasque-1)->EnableMenuItem(1,MF_ENABLED|MF_BYPOSITION);
		SetCursor(LoadCursor(0,IDC_ARROW));

	}
	else 
	{
		ValideMenu=FALSE;
		menu.GetSubMenu(PositionMenuMasque-1)->EnableMenuItem(1,MF_GRAYED|MF_BYPOSITION);
	}
	ArreteLeContinu();
}

LONG CTraitImageWnd::OnFinClasse(UINT uiParam,LONG lparam )
{
	if ((CImageClas *)lparam==ImageAppr) DefImagAppr=FALSE;
	return 0;
}

BOOL CTraitImageWnd::ChargementDesParam(CString initFile)
{
	CFile Fichier;
	BYTE NPAR,nva,nbn;
	int Pos,nvint;
	CString Cval,CTemp;

	if(Fichier.Open(initFile, CFile::modeRead))
	{
		Fichier.Read(&NPAR,sizeof(BYTE));

		for (int i=0;i<(int)NPAR;i++)
		{

			Fichier.Read(&nva,sizeof(BYTE));
			DialPara.IChoisi[i]=nva;// indice du param
			Pos=DialPara.m_dispo[nva].Find(":");
			if (Pos==-1) 
			{
				MessageBox("Erreur syntaxe paramètre (:)","Lecture .param",MB_OK|MB_ICONSTOP);
				return 0; 
			}
			DialPara.m_choisi[i]=DialPara.m_dispo[nva].Left(Pos+1);
			Fichier.Read(&nbn,sizeof(BYTE));
			Cval="";
			for (int j=0;j<nbn;j++)
			{
				Fichier.Read(&nvint,sizeof(int));
				CTemp.Format("%d",nvint);
				if (nvint>=0) Cval+="+";
				Cval=Cval+CTemp;
				if (j!=(nbn-1)) Cval+=",";
			}
			DialPara.m_choisi[i]+=Cval;
		}//i 
		DialPara.nbpara=NPAR;
	}
	else return FALSE;

	return TRUE;
}

void CTraitImageWnd::AlloueParam() 
{
	int i;
	NbParam=DialPara.nbpara;
	SetCursor(LoadCursor(NULL,IDC_WAIT));

	if (NbParam==0)
	{
		MessageBox("Définissez les parametres!","Attention",MB_OK);
		return;
	}
	if (NbParam>=100)
	{
		MessageBox("Nombre de paramètres très important!\nRisque de dépassement mémoire!","Attention",MB_OK);
		return;
	}
	if(NbParam>NBPMAX) NBPMAX=NbParam;
	for (i=0;i<NbParam;i++)
	{
		if (VisuPar[i]==NULL)  
		{	
			VisuPar[i]=new CImageWnd(Larg,Haut);
  			if (VisuPar[i]->Image==NULL) {MessageBox("Erreur memoire paramètre","Message",MB_OK);return;}
  			VisuPar[i]->Create((LPCSTR)DialPara.m_choisi[i], WS_CHILD | WS_CAPTION | WS_BORDER | WS_THICKFRAME ,NULL);
			VisuPar[i]->ShowWindow(SW_HIDE); 
		}
	}  // i
	SetCursor(LoadCursor(0,IDC_ARROW));

/*	if (Visu==TRUE) 
	{
		MontrerLesVisuPar();
	} // visu*/
} // alloue 

void CTraitImageWnd::OnApprentissageDfinirlesparamtres() 
{
	ArreteLeContinu();
	

	for (int i=0;i<NbParam;i++)
	{
		if (VisuPar[i]!=NULL) VisuPar[i]->ShowWindow(SW_HIDE);	

		Visu=FALSE;
	}

	if (DialPara.DoModal()==IDOK) 
	{
		AlloueParam();
		SauveParam(NomFichParam);
		ParamCalcul=FALSE;
	}
}

void CTraitImageWnd::OnUpdateApprentissageDfinirlesparamtres(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideMenu);
}

BOOL CTraitImageWnd::SauveParam(CString paramFile)
{
	CFile Fichier;
	BYTE nbn;
	int Pos,Decale;
	CString Cval,CTemp;

	if(Fichier.Open(paramFile,CFile::modeCreate | CFile::modeWrite))
	{
		//ecriture du nb de param
		Fichier.Write(&DialPara.nbpara,sizeof(BYTE));

		for (int i=0;i<(long)DialPara.nbpara;i++)
		{
			Fichier.Write(&DialPara.IChoisi[i],sizeof(BYTE)); // indice du param

			//recherche du nombre de n dans le paramètre
			Pos=DialPara.m_choisi[i].Find(":");
			Cval=DialPara.m_choisi[i].Right(DialPara.m_choisi[i].GetLength()-Pos-1);

			nbn=1;
			Pos=Cval.Find(",");
			while(Pos!=-1)
			{
				nbn++;
				Cval=Cval.Right(Cval.GetLength()-Pos-1);
				Pos=Cval.Find(",");
			} 

			//ecriture du nb de n
			Fichier.Write(&nbn,sizeof(BYTE));
			
			//recherche et ecriture des n
			Cval=DialPara.m_choisi[i];
			Pos=Cval.Find(":");
			do
			{	
				Cval=Cval.Right(Cval.GetLength()-Pos-1);
				Pos=Cval.Find(",");

				CTemp=Cval.SpanExcluding(",");
				if (CTemp==Cval) Decale=atoi(Cval);
				else 	Decale=atoi(Cval.SpanExcluding(","));

				Fichier.Write(&Decale,sizeof(int));
			} while(CTemp!=Cval);
		}//i
	}
	else return FALSE;

	return TRUE;
}

void CTraitImageWnd::OnParamtrageVisualisationparamtres() 
{
	ArreteLeContinu();
	if (!ParamCalcul)
		for (int i=0;i<NbParam;i++) 
		{ 
			if(!CalculParam(i)) 
			{
				ParamCalcul=FALSE;
				return;
			}
		}

	if(Visu==TRUE) 
	{
		CacherLesVisuPar();
		Visu=FALSE;
	}
	else 
	{
		MontrerLesVisuPar();
		Visu=TRUE;
		GetMDIFrame()->MDITile(MDITILE_VERTICAL);
	}
}

void CTraitImageWnd::OnUpdateParamtrageVisualisationparamtres(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideMenu &&(NbParam!=0));
	pCmdUI->SetCheck(Visu);
}


void CTraitImageWnd::CacherLesVisuPar()
{
	for (int i=0;i<NbParam;i++)
  		{           
	  		if (VisuPar[i]!=NULL) 
			{
				VisuPar[i]->ShowWindow(SW_HIDE);
				VisuPar[i]->SetWindowText(DialPara.m_choisi[i]);
			}
  		}
}

void CTraitImageWnd::MontrerLesVisuPar()
{
	for (int i=0;i<NbParam;i++)
  	{           
  		if (VisuPar[i]!=NULL) 
		{
			VisuPar[i]->SetWindowText(DialPara.m_choisi[i]);
  			VisuPar[i]->ShowWindow(SW_SHOW); 
		}
  	} 
}

BOOL CTraitImageWnd::CalculParam(int i)
{
	SetCursor(LoadCursor(0,IDC_WAIT));
	ParamCalcul=TRUE;
	if (PParametre==NULL) PParametre=new CParam();
	switch (DialPara.IChoisi[i])
 	{          
 		case 0 :return PParametre->LumPix(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut); break;
 		case 1 :return PParametre->LumMoyGlis(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
 		case 2 :return PParametre->LumMoyFix(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
 		case 3 :return PParametre->GradPix(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
 		case 4 :return PParametre->GradMoyGlis(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
 		case 5 :return PParametre->GradMoyFix(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 6 :return PParametre->LumMoyDecLigFix(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 7 :return PParametre->Entropie(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 8 :return PParametre->EntropieGlis(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 9 :return PParametre->ContrastLoc(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 10:return PParametre->Max_Min(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 11:return PParametre->ContrastGlis(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 12:return PParametre->Max_MinGlis(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 13:return PParametre->MedianGlissant(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 14:return PParametre->MedianFix(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 15:return PParametre->ImageEnregistre(VisuPar[i]->Image,Larg,Haut);break;
		case 16:return PParametre->SobelFix(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 17:return PParametre->SobelGlissant(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 18:return PParametre->RehaussementFix(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 19:return PParametre->RehaussementGlis(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 20:return PParametre->Ouverture(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 21:return PParametre->Fermeture(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 22:return PParametre->EcartmedFen(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 23:return PParametre->EcartmedFix(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 24:return PParametre->MoyAdap(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 25:return PParametre->MoyAdapAuto(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 26:return PParametre->MoyTronqueglis(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 27:return PParametre->MoyTronqueFix(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 28:return PParametre->Maxrepglis(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 29:return PParametre->NbMaxrepglis(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 30:return PParametre->GradMoyGlis2(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 31:return PParametre->FermetureRapide(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 32:return PParametre->OuvertureRapide(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 33:return PParametre->ContrastrapideGlis(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 34:return PParametre->Gradientrapide(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 35:return PParametre->Lumirapide(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 36:return PParametre->GradientrapideNonNorm(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 37:return PParametre->EcartTypeFenGlis(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 38:return PParametre->MaxFenGlis(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 39:return PParametre->DensMax33(DialPara.m_choisi[i],VisuPar[i]->Image,Image,Larg,Haut);break;
		case 0xff :return FALSE; break;
 	} // switch                     
return FALSE;
} // calcul

void CTraitImageWnd::OnApprentissageGnrerlefichierdapprentissage() 
{
	ArreteLeContinu();

	int l,c,t;
	double taille,rd;
	long   *hram,*tempo,p=0;
	long n;
	CFile Fichier;

	CDialNbPointAppr Dlg(this);
	if(Dlg.DoModal()==IDOK)
	{
		SetCursor(LoadCursor(0,IDC_WAIT));
		n=Dlg.m_nbpoint;

		hram=new long[n];
		if(hram==NULL) {MessageBox("Erreur hram","Message",MB_OK);return;}

		tempo=new long[Haut*Larg];
		if (tempo==NULL) {MessageBox("Erreur tempo","Message",MB_OK);return;}

		if (!ParamCalcul) 
			for (int i=0;i<NbParam;i++) 
			{
				if(!CalculParam(i)) 
				{
					ParamCalcul=FALSE;
					return;// FALSE;
				}
			}

		p=0;

		for (l=0;l<Haut;l++)
			for (c=0;c<Larg;c++)
			{
			   if (ImageAppr->Image[Larg*l+c]==255) {tempo[p++]=Larg*l+c;} 
			}

		taille=double (p);
		if(taille>RAND_MAX) 
		{
			MessageBox("Zone d'apprentissage trop grande","Erreur",MB_ICONSTOP);
			return;// FALSE;
		}

		if(taille==0)
		{
			MessageBox("Zone d'apprentissage non définie","Erreur",MB_ICONSTOP);
			return;// FALSE;
		}

		srand((unsigned)time(NULL));
		p=0;	  
		do 
		{
			rd=(double)rand()/(double)(RAND_MAX);
			t=(long)(rd*taille);
			if (t<0 ) { MessageBox("Erreur aléatoire");return;}
			//Image[tempo[t]]=255;
			hram[p]=tempo[t];
			p++;
		} while(p<n);

		if(Fichier.Open(NomFichAppr,CFile::modeCreate|CFile::modeWrite))
		{
			Fichier.Write((LPCSTR)NomFichParam,1024);
			Fichier.Write(&NbParam,sizeof(BYTE));
			Fichier.Write(&n,sizeof(long));

			for(t=0;t<n;t++)
			{
				p=hram[t];
				Calc(p);
				Fichier.Write(Buffer,NbParam);
			}
			Fichier.Close();
		}
		else	MessageBox("Nom du fichier .appr incorrect!","Erreur",MB_ICONSTOP);
		CStdioFile MonFic;
		CString ChaineTemp;
		float fval;
		CString NomFicDat;
		NomFicDat=NomFichAppr.Left(NomFichAppr.GetLength()-5)+".dat";
		if(MonFic.Open(NomFicDat,CFile::modeCreate|CFile::modeWrite))
		{
			ChaineTemp.Format("%d\n",NbParam);
			MonFic.WriteString(ChaineTemp);
			ChaineTemp.Format("%d\n",n);
			MonFic.WriteString(ChaineTemp);
			for(t=0;t<n;t++)
			{
				p=hram[t];
				Calc(p);
				for (int ppar=0;ppar<NbParam;ppar++)
				{
					fval=(float)Buffer[ppar]/255.0;
					ChaineTemp.Format("%f ",fval);
					MonFic.WriteString(ChaineTemp);
				}

					MonFic.WriteString("\n");
			}
			MonFic.Close();
		}
		else	MessageBox("Nom du fichier .appr incorrect!","Erreur",MB_ICONSTOP);


		delete tempo;
		delete hram;
		SetCursor(LoadCursor(0,IDC_ARROW));
	}//domodal
}

void CTraitImageWnd::OnUpdateApprentissageGnrerlefichierdapprentissage(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideMenu&&DefImagAppr&&(NbParam!=0));
}

void CTraitImageWnd::OnApprentissageConstruirelefichiermulticlasses() 
{
	ArreteLeContinu();
	CacherLesVisuPar();
	Visu=FALSE;
	CDialConfigPoly Dlg(this,NomFichCMC);
	if(Dlg.DoModal()==IDOK)
	{
		NomFichCMC=Dlg.NomFich;
		NomFichMMC=NomFichCMC.Left(NomFichCMC.ReverseFind('.'))+".mmc";
		if(Dlg.m_tree.FichParam!=NomFichParam)
		{
			NomFichParam=Dlg.m_tree.FichParam;
			DialPara.nbpara=0;
			ChargementDesParam(NomFichParam);
			AlloueParam();
		}
		ParamCalcul=FALSE;
		ValideDecision=FALSE;
		ValideMenu=FALSE;
		menu.GetSubMenu(PositionMenuMasque-1)->EnableMenuItem(1,MF_GRAYED|MF_BYPOSITION);
		ValideVisuMes=TRUE;
	}
}

void CTraitImageWnd::OnApprentissageGnrerunfichierpolytopesmulticlasses() 
{
	CDialPmc Dlg(this,NomFichMMC,NomFichPMC,CoeffPoly);
	if(Dlg.DoModal()==IDOK)
	{
		if(Dlg.NomFichMMC!=NomFichMMC)
		{
			ValideVisuMes=FALSE;
		}
		NomFichMMC=Dlg.NomFichMMC;
		NomFichPMC=Dlg.NomFichPMC;
		CoeffPoly=Dlg.m_coefR;

		CPolyMc PolyMc;

		SetCursor(LoadCursor(NULL,IDC_WAIT));
		
		if(PolyMc.Init(NomFichMMC))
			{
				if (PolyMc.Verif_redondance()) MessageBox("Mesures redondantes de classes différentes !",MB_OK);
				PolyMc.coef=CoeffPoly;
				PolyMc.polytope();
			/*	PolyMc.Verifrecouvrtemp();
				CString ch;
				ch.Format("nb recouvrement = %d",PolyMc.Verif_recouvr());
				MessageBox(ch,"recouvrement");	
			*/	PolyMc.sauve_pol(NomFichPMC); 
				PolyMc.FinAlloc();
				PolyMc.AffPol();
			} else MessageBox("Erreur","",MB_OK);
		SetCursor(LoadCursor(0,IDC_ARROW));		
	}
}



void CTraitImageWnd::OnDcisionConfigurationdeladcision() 
{
	CacherLesVisuPar();
	Visu=FALSE;
	if(Visumes) delete AffMesPol;
	CDialConfigDec Dlg(this,NomFichPMC,NomFichParam,AffectNonClas,ClassAffect);	

	if(Dlg.DoModal()==IDOK)
	{	
		NomFichPMC=Dlg.NomFichPMC;
		if(NomFichParam!=Dlg.NomFichParam) ParamCalcul=FALSE;
		NomFichParam=Dlg.NomFichParam;
		AffectNonClas=Dlg.m_check;
		ClassAffect=Dlg.m_classeaffect;
		ChargementDesParam(NomFichParam);
		AlloueParam();
		ValideDecision=TRUE;
		ValideMenu=FALSE;
		ValideVisuMes=FALSE;
		menu.GetSubMenu(PositionMenuMasque-1)->EnableMenuItem(1,MF_GRAYED|MF_BYPOSITION);
	}
	else ValideDecision=FALSE;
}

void CTraitImageWnd::OnDcisionDcision() 
{
	CPolyMc Polymc;
	long p,kpol,clas,k;
	int l,c,par,i,nbclas;
	BYTE ext,decide,aa,nonclas;


	long *zone;
	unsigned char *classecoul;

	ArreteLeContinu();
	
	Polymc.LirePol(NomFichPMC);	
	Polymc.AffPol();

	nbclas=Polymc.NbClass;

	if(!AffectNonClas) 
	{
		nbclas++;
	}	

	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
		zone[i]=0;
	}
	  for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}

	SetCursor(LoadCursor(NULL,IDC_WAIT));
	
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			Calc(p);
			decide=0;
			nonclas=1;

			for( clas=0;clas<Polymc.NbClass;clas++)
			{
			   for(kpol=0;kpol<Polymc.Npol[clas];kpol++)
  				{ ext=0;
  				 for(par=0;par<Polymc.NPAR;par++)
 				   { 
  					   k=kpol*Polymc.NPAR+par;
 					   aa=Buffer[par];
 					   if (aa<Polymc.poli[clas][k]) {ext=1;break;}
 					   else   if (aa>Polymc.pols[clas][k]) {ext=1;break;}
   				   }//for par
				   if (!ext) {decide=classecoul[clas];zone[clas]++;nonclas=0;break;}
				}//for pol
			}//for clas
			if(nonclas)
			{
				if(AffectNonClas)
				{
					zone[ClassAffect]++;
					decide=classecoul[ClassAffect];
				}
				else 
				{
					zone[nbclas-1]++;
					decide=classecoul[nbclas-1];
				}
			}
			Image[p]=decide;
			p++;
		} // c
	} // l
	if(TypePalette!=PALCOUL) OnCouleur();
	Invalidate(FALSE);

	CString ch="",ch0;;
	for(i=0;i<Polymc.NbClass;i++)
	{
		ch0.Format("Classe %d= %d pixels\n\n",i,zone[i]);
		ch+=ch0;
	}
	if(!AffectNonClas)
	{
		ch0.Format("Non classés= %d pixels\n\n",zone[i]);
		ch+=ch0;
	}
	if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(Polymc.NbClass,!AffectNonClas);Legende->Invalidate(TRUE);}
	else
	{
		Legende->Init(Polymc.NbClass,!AffectNonClas);
		Legende->ShowWindow(SW_SHOW);
		Legende->Invalidate(TRUE);
	}
	MessageBox(ch);
	Polymc.Fin();
	delete classecoul;
	delete zone;
	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::OnUpdateDcisionDcision(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision);
}

void CTraitImageWnd::OnDcisionDcisionavecmasque() 
{
	CPolyMc Polymc;
	long p,kpol,clas,k;
	int l,c,par,i,nbclas;
	BYTE ext,decide,aa,nonclas;
	CPoint h;

	long *zone;
	unsigned char *classecoul;

	ArreteLeContinu();
	
	Polymc.LirePol(NomFichPMC);	
	Polymc.AffPol();

	nbclas=Polymc.NbClass;

	if(!AffectNonClas) 
	{
		nbclas++;
	}	

	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
		zone[i]=0;
	}
//	AlloueParam();
	// if (!ParamCalcul) 
	  for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}
	/////////////////////////////////////////
	OnMasqueRechercheAutom(40);
	debut_i+=12;//8
	debut_j+=12;//8
	fin_j-=12;//8
	fin_i-=12;//8
	AffectRegionMasque();
	SetCursor(LoadCursor(NULL,IDC_WAIT));
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			h=Inv_Get_ij(c,l);
			if(RegionEllipse.PtInRegion(h))
			{	
				Calc(p);
				decide=0;
				nonclas=1;
				for( clas=0;clas<Polymc.NbClass;clas++)
				{
					
					for(kpol=0;kpol<Polymc.Npol[clas];kpol++)
  					{ 
						ext=0;

						for(par=0;par<Polymc.NPAR;par++)
						{ 
						   k=kpol*Polymc.NPAR+par;
						   aa=Buffer[par];
						   if ((aa<Polymc.poli[clas][k])||(aa>Polymc.pols[clas][k])) {ext=1;break;}
						}//for par
						if (!ext){decide=classecoul[clas];zone[clas]++;nonclas=0;break;}				 
					}//for pol					
					if(!ext) break;
					
				}//for clas
				if(nonclas)
				{
					if(AffectNonClas)
					{
						zone[ClassAffect]++;
						decide=classecoul[ClassAffect];
					}
					else 
					{
						zone[nbclas-1]++;
						decide=classecoul[nbclas-1];
					}
				}
				Image[p]=decide;
			}//if region
			p++;
		} // c
	} // l
	if(TypePalette!=PALCOUL) OnCouleur();
	Invalidate(FALSE);

	CString ch="",ch0;
	for(i=0;i<Polymc.NbClass;i++)
	{
		ch0.Format("Classe %d= %d pixels\n\n",i,zone[i]);
		ch+=ch0;
	}
	if(!AffectNonClas)
	{
		ch0.Format("Non classés= %d pixels\n\n",zone[i]);
		ch+=ch0;
	}
	if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(Polymc.NbClass,!AffectNonClas);Legende->Invalidate(TRUE);}
	else
	{
		Legende->Init(Polymc.NbClass,!AffectNonClas);
		Legende->ShowWindow(SW_SHOW);
		Legende->Invalidate(TRUE);
	}
	MessageBox(ch);
//	OnMasqueAnnuler();
	Polymc.Fin();
	SetCursor(LoadCursor(0,IDC_ARROW));
	delete classecoul;
	delete zone;
}

void CTraitImageWnd::OnUpdateDcisionDcisionavecmasque(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision);
}
void CTraitImageWnd::OnDcisionDcision2masques() 
{
	CPolyMc Polymc;
	long p,kpol,clas,k;
	int l,c,par,i,nbclas;
	BYTE ext,decide,aa,nonclas;
	CPoint h;
	CRgn petiteellipse;
	long *zone;
	unsigned char *classecoul;

	ArreteLeContinu();
	
	Polymc.LirePol(NomFichPMC);	
	Polymc.AffPol();

	nbclas=Polymc.NbClass;

	if(!AffectNonClas) 
	{
		nbclas++;
	}	

	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
		zone[i]=0;
	}
//	AlloueParam();
	// if (!ParamCalcul) 
	  for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}
	/////////////////////////////////////////
	OnCalculellipse();
	OnMasqueRechercheAutom(10);
/*	à remettre
	debut_i+=12;//8
	debut_j+=12;//8
	fin_j-=12;//8
	fin_i-=12;//8

*/	
	debut_i+=7;//8
	debut_j+=7;//8
	fin_j=debut_j+245;//8
	fin_i=debut_i+245;//8
	AffectRegionMasque();
	petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+20,debut_j+19),Inv_Get_ij(fin_i-19,fin_j-20)));
	SetCursor(LoadCursor(NULL,IDC_WAIT));
	
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			h=Inv_Get_ij(c,l);
			if(RegionEllipse.PtInRegion(h))
			{	
				Calc(p);
				decide=0;
				nonclas=1;

				for( clas=0;clas<Polymc.NbClass;clas++)
				{
					if(((clas!=2)&&(clas!=1/*à enlever*/))||(petiteellipse.PtInRegion(h)))
					{
						for(kpol=0;kpol<Polymc.Npol[clas];kpol++)
  						{ 
							ext=0;

							for(par=0;par<Polymc.NPAR;par++)
							{ 
							   k=kpol*Polymc.NPAR+par;
							   aa=Buffer[par];
							   if ((aa<Polymc.poli[clas][k])||(aa>Polymc.pols[clas][k])) {ext=1;break;}
							}//for par
							if (!ext) {decide=classecoul[clas];zone[clas]++;nonclas=0;break;}				 
						}//for pol					
					if(!ext) break;
					}
				}//for clas
				if(nonclas)
				{
					if(AffectNonClas)
					{
						zone[ClassAffect]++;
						decide=classecoul[ClassAffect];
					}
					else 
					{
						zone[nbclas-1]++;
						decide=classecoul[nbclas-1];
					}
				}
				Image[p]=decide;
			}//if region
			p++;
		} // c
	} // l
	if(TypePalette!=PALCOUL) OnCouleur();
	
	Invalidate(FALSE);

	CString ch="",ch0;
	for(i=0;i<Polymc.NbClass;i++)
	{
		ch0.Format("Classe %d= %d pixels\n\n",i,zone[i]);
		ch+=ch0;
	}
	if(!AffectNonClas)
	{
		ch0.Format("Non classés= %d pixels\n\n",zone[i]);
		ch+=ch0;
	}
	if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(Polymc.NbClass,!AffectNonClas);Legende->Invalidate(TRUE);}
	else
	{
		Legende->Init(Polymc.NbClass,!AffectNonClas);
		Legende->ShowWindow(SW_SHOW);
		Legende->Invalidate(TRUE);
	}
	if(zone[1]>37) 
		ch+="MANQUE\n";//80
	if(zone[2]>=50) 
		ch+="LISSE\n";
	if(zone[3]>=14) 
		ch+="BOURGEON\n";
	if (zone[4]>=42)
		ch+="ARRACHEMENT\n";	
	if((zone[1]<=37)&&(zone[2]<50)&&(zone[3]<14)&&(zone[4]<42)) 
		ch+="BONNE";
	MessageBox(ch);
//	OnMasqueAnnuler();
	Polymc.Fin();
	SetCursor(LoadCursor(0,IDC_ARROW));
	delete classecoul;
	delete zone;
	//sndPlaySound("c:\\Users\\Pierreg\\bonne.wav",SND_SYNC);
}

void CTraitImageWnd::OnUpdateDcisionDcision2masques(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision);
}
void CTraitImageWnd::OnDcisionDcision3masques() 
{
	CPolyMc Polymc;
	long p,kpol,clas,k;
	int l,c,par,i,nbclas;
	BYTE ext,decide,aa,nonclas;
	CPoint h;
	CRgn petiteellipse,moyenneellipse;
	long *zone;
	unsigned char *classecoul;

	ArreteLeContinu();
	
	Polymc.LirePol(NomFichPMC);	
	Polymc.AffPol();

	nbclas=Polymc.NbClass;

	if(!AffectNonClas) 
	{
		nbclas++;
	}	

	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
		zone[i]=0;
	}
//	AlloueParam();
	// if (!ParamCalcul) 
	  for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}
	  OnCalculellipse();
	/////////////////////////////////////////
	OnMasqueRechercheAutom(10);
/*	à remettre
	debut_i+=12;//8
	debut_j+=12;//8
	fin_j-=12;//8
	fin_i-=12;//8*/
	///////////////////////////////////
	////////2 lignes à remettre pour le masque fixe
	////////////////////////////////
	fin_j=debut_j+258;//8
	fin_i=debut_i+258;//8

	AffectRegionMasque();

	moyenneellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+11,debut_j+10),Inv_Get_ij(fin_i-10,fin_j-11)));
	petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+27,debut_j+26),Inv_Get_ij(fin_i-26,fin_j-27)));
	SetCursor(LoadCursor(NULL,IDC_WAIT));
	
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			h=Inv_Get_ij(c,l);
			if(RegionEllipse.PtInRegion(h))
			{	
				Calc(p);
				decide=0;
				nonclas=1;

				for( clas=0;clas<Polymc.NbClass;clas++)
				{
					if((clas==3)||(moyenneellipse.PtInRegion(h)&&clas==4)||(petiteellipse.PtInRegion(h)))
					{
						for(kpol=0;kpol<Polymc.Npol[clas];kpol++)
  						{ 
							ext=0;

							for(par=0;par<Polymc.NPAR;par++)
							{ 
							   k=kpol*Polymc.NPAR+par;
							   aa=Buffer[par];
							   if ((aa<Polymc.poli[clas][k])||(aa>Polymc.pols[clas][k])) {ext=1;break;}
							}//for par
							if (!ext) {decide=classecoul[clas];zone[clas]++;nonclas=0;break;}				 
						}//for pol					
					if(!ext) break;
					}
				}//for clas
				if(nonclas)
				{
					if(AffectNonClas)
					{
						zone[ClassAffect]++;
						decide=classecoul[ClassAffect];
					}
					else 
					{
						zone[nbclas-1]++;
						decide=classecoul[nbclas-1];
					}
				}
				Image[p]=decide;
			}//if region
			else Image[p]=0;
			p++;
		} // c
	} // l
	if(TypePalette!=PALCOUL) OnCouleur();
//	OnCalculellipse();
	Invalidate(FALSE);

	CString ch="",ch0;
	for(i=0;i<Polymc.NbClass;i++)
	{
		ch0.Format("Classe %d= %d pixels\n\n",i,zone[i]);
		ch+=ch0;
	}
	if(!AffectNonClas)
	{
		ch0.Format("Non classés= %d pixels\n\n",zone[i]);
		ch+=ch0;
	}
	if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(Polymc.NbClass,!AffectNonClas);Legende->Invalidate(TRUE);}
	else
	{
		Legende->Init(Polymc.NbClass,!AffectNonClas);
		Legende->ShowWindow(SW_SHOW);
		Legende->Invalidate(TRUE);
	}
	if(zone[1]>37) 
		ch+="MANQUE\n";//80
	if(zone[2]>=50) 
		ch+="LISSE\n";
	if(zone[3]>=14) 
		ch+="BOURGEON\n";
	if (zone[4]>=42)
		ch+="ARRACHEMENT\n";	
	if((zone[1]<=37)&&(zone[2]<50)&&(zone[3]<14)&&(zone[4]<42)) 
		ch+="BONNE";
	MessageBox(ch);
//	OnMasqueAnnuler();
	Polymc.Fin();
	SetCursor(LoadCursor(0,IDC_ARROW));
	delete classecoul;
	delete zone;
}

void CTraitImageWnd::OnUpdateDcisionDcision3masques(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision);
}

void CTraitImageWnd::OnDcisionDcision3masquesadapt() 
{
	CPolyMc Polymc;
	long p,kpol,clas,k;
	int l,c,par,i,nbclas;
	BYTE ext,decide,aa,nonclas;
	CPoint h;
	CRgn petiteellipse,moyenneellipse;
	long *zone;
	unsigned char *classecoul;

	ArreteLeContinu();
	
	Polymc.LirePol(NomFichPMC);	
	Polymc.AffPol();

	nbclas=Polymc.NbClass;

	if(!AffectNonClas) 
	{
		nbclas++;
	}	

	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
		zone[i]=0;
	}
//	AlloueParam();
	// if (!ParamCalcul) 
	  for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}
	  OnCalculellipse();
	/////////////////////////////////////////
	OnMasqueRechercheAutom(10);
/*	à remettre
	debut_i+=12;//8
	debut_j+=12;//8
	fin_j-=12;//8
	fin_i-=12;//8/*
	///////////////////////////////////
	////////2 lignes à remettre pour le masque fixe
	////////////////////////////////
	fin_j=debut_j+258;//8
	fin_i=debut_i+258;//8
*/
	AffectRegionMasque();

	moyenneellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+11,debut_j+10),Inv_Get_ij(fin_i-10,fin_j-11)));
	petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+27,debut_j+26),Inv_Get_ij(fin_i-26,fin_j-27)));
	SetCursor(LoadCursor(NULL,IDC_WAIT));
	
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			h=Inv_Get_ij(c,l);
			if(RegionEllipse.PtInRegion(h))
			{	
				Calc(p);
				decide=0;
				nonclas=1;

				for( clas=0;clas<Polymc.NbClass;clas++)
				{
					if((clas==3)||(moyenneellipse.PtInRegion(h)&&clas==4)||(petiteellipse.PtInRegion(h)))
					{
						for(kpol=0;kpol<Polymc.Npol[clas];kpol++)
  						{ 
							ext=0;

							for(par=0;par<Polymc.NPAR;par++)
							{ 
							   k=kpol*Polymc.NPAR+par;
							   aa=Buffer[par];
							   if ((aa<Polymc.poli[clas][k])||(aa>Polymc.pols[clas][k])) {ext=1;break;}
							}//for par
							if (!ext) {decide=classecoul[clas];zone[clas]++;nonclas=0;break;}				 
						}//for pol					
					if(!ext) break;
					}
				}//for clas
				if(nonclas)
				{
					if(AffectNonClas)
					{
						zone[ClassAffect]++;
						decide=classecoul[ClassAffect];
					}
					else 
					{
						zone[nbclas-1]++;
						decide=classecoul[nbclas-1];
					}
				}
				Image[p]=decide;
			}//if region
			else Image[p]=0;
			p++;
		} // c
	} // l
	if(TypePalette!=PALCOUL) OnCouleur();
//	OnCalculellipse();
	Invalidate(FALSE);

	CString ch="",ch0;
	for(i=0;i<Polymc.NbClass;i++)
	{
		ch0.Format("Classe %d= %d pixels\n\n",i,zone[i]);
		ch+=ch0;
	}
	if(!AffectNonClas)
	{
		ch0.Format("Non classés= %d pixels\n\n",zone[i]);
		ch+=ch0;
	}
	if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(Polymc.NbClass,!AffectNonClas);Legende->Invalidate(TRUE);}
	else
	{
		Legende->Init(Polymc.NbClass,!AffectNonClas);
		Legende->ShowWindow(SW_SHOW);
		Legende->Invalidate(TRUE);
	}
	if(zone[1]>37) 
		ch+="MANQUE\n";//80
	if(zone[2]>=50) 
		ch+="LISSE\n";
	if(zone[3]>=14) 
		ch+="BOURGEON\n";
	if (zone[4]>=42)
		ch+="ARRACHEMENT\n";	
	if((zone[1]<=37)&&(zone[2]<50)&&(zone[3]<14)&&(zone[4]<42)) 
		ch+="BONNE";
	MessageBox(ch);
//	OnMasqueAnnuler();
	Polymc.Fin();
	SetCursor(LoadCursor(0,IDC_ARROW));
	delete classecoul;
	delete zone;
}

void CTraitImageWnd::OnUpdateDcisionDcision3masquesadapt(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision);
}

void CTraitImageWnd::OnDcisionDcision3masquescg() 
{
	CPolyMc Polymc;
	long p,kpol,clas,k;
	int l,c,par,i,nbclas;
	BYTE ext,decide,aa,nonclas;
	CPoint h;
	CRgn petiteellipse,moyenneellipse;
	long *zone;
	unsigned char *classecoul;

	ArreteLeContinu();
	
	Polymc.LirePol(NomFichPMC);	
	Polymc.AffPol();

	nbclas=Polymc.NbClass;

	if(!AffectNonClas) 
	{
		nbclas++;
	}	

	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
		zone[i]=0;
	}
//	AlloueParam();
	// if (!ParamCalcul) 
	  for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}
	  if(!MasqueCG()) {MessageBox("problèmes de bords","décision",MB_OK); }
/*	à remettre
	debut_i+=12;//8
	debut_j+=12;//8
	fin_j-=12;//8
	fin_i-=12;//8/*
	///////////////////////////////////
	////////2 lignes à remettre pour le masque fixe
	////////////////////////////////
	fin_j=debut_j+258;//8
	fin_i=debut_i+258;//8
*/
	AffectRegionMasque();
/*
	moyenneellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+11,debut_j+10),Inv_Get_ij(fin_i-10,fin_j-11)));
	petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+27,debut_j+26),Inv_Get_ij(fin_i-26,fin_j-27)));
*/
	int deci=(int) (((float)(fin_i-debut_i))/23.45);
	int deciprim=(int) (((float)(fin_i-debut_i))/25.8);
	int decj=(int) (((float)(fin_j-debut_j))/25.8);
	int decjprim=(int) (((float)(fin_j-debut_j))/23.45);
	moyenneellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+deci,debut_j+decj),Inv_Get_ij(fin_i-deciprim,fin_j-decjprim)));
	deci=(int) (((float)(fin_i-debut_i))/9.5);
	decj=(int) (((float)(fin_j-debut_j))/9.92);
	deciprim=(int) (((float)(fin_i-debut_i))/9.92);
	decjprim=(int) (((float)(fin_j-debut_j))/9.5);
	petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+deci,debut_j+decj),Inv_Get_ij(fin_i-deciprim,fin_j-decjprim)));

	SetCursor(LoadCursor(NULL,IDC_WAIT));
	
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			h=Inv_Get_ij(c,l);
			if(RegionEllipse.PtInRegion(h))
			{	
				Calc(p);
				decide=0;
				nonclas=1;

				for( clas=0;clas<Polymc.NbClass;clas++)
				{
					if((clas==3)||(moyenneellipse.PtInRegion(h)&&clas==4)||(petiteellipse.PtInRegion(h)))
					{
						for(kpol=0;kpol<Polymc.Npol[clas];kpol++)
  						{ 
							ext=0;

							for(par=0;par<Polymc.NPAR;par++)
							{ 
							   k=kpol*Polymc.NPAR+par;
							   aa=Buffer[par];
							   if ((aa<Polymc.poli[clas][k])||(aa>Polymc.pols[clas][k])) {ext=1;break;}
							}//for par
							if (!ext) {decide=classecoul[clas];zone[clas]++;nonclas=0;break;}				 
						}//for pol					
					if(!ext) break;
					}
				}//for clas
				if(nonclas)
				{
					if(AffectNonClas)
					{
						zone[ClassAffect]++;
						decide=classecoul[ClassAffect];
					}
					else 
					{
						zone[nbclas-1]++;
						decide=classecoul[nbclas-1];
					}
				}
				Image[p]=decide;
			}//if region
			else Image[p]=0;
			p++;
		} // c
	} // l
	if(TypePalette!=PALCOUL) OnCouleur();
//	OnCalculellipse();
	Invalidate(FALSE);

	CString ch="",ch0;
	for(i=0;i<Polymc.NbClass;i++)
	{
		ch0.Format("Classe %d= %d pixels\n\n",i,zone[i]);
		ch+=ch0;
	}
	if(!AffectNonClas)
	{
		ch0.Format("Non classés= %d pixels\n\n",zone[i]);
		ch+=ch0;
	}
	if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(Polymc.NbClass,!AffectNonClas);Legende->Invalidate(TRUE);}
	else
	{
		Legende->Init(Polymc.NbClass,!AffectNonClas);
		Legende->ShowWindow(SW_SHOW);
		Legende->Invalidate(TRUE);
	}
	if(zone[1]>37) 
		ch+="MANQUE\n";//80
	if(zone[2]>=50) 
		ch+="LISSE\n";
	if(zone[3]>=14) 
		ch+="BOURGEON\n";
	if (zone[4]>=42)
		ch+="ARRACHEMENT\n";	
	if((zone[1]<=37)&&(zone[2]<50)&&(zone[3]<14)&&(zone[4]<42)) 
		ch+="BONNE";
	MessageBox(ch);
//	OnMasqueAnnuler();
	Polymc.Fin();
	SetCursor(LoadCursor(0,IDC_ARROW));
	delete classecoul;
	delete zone;
}

void CTraitImageWnd::OnUpdateDcisionDcision3masquescg(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision);
}
void CTraitImageWnd::OnDcisionDcision3masquessansrech() 
{
	CPolyMc Polymc;
	long p,kpol,clas,k;
	int l,c,par,i,nbclas;
	BYTE ext,decide,aa,nonclas;
	CPoint h;
	CRgn petiteellipse;//,moyenneellipse;
	long *zone;
	unsigned char *classecoul;

	ArreteLeContinu();
	
	Polymc.LirePol(NomFichPMC);	
	Polymc.AffPol();

	nbclas=Polymc.NbClass;

	if(!AffectNonClas) 
	{
		nbclas++;
	}	

	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];//////////////////

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
		zone[i]=0;
	}
//	AlloueParam();
	// if (!ParamCalcul) 
	  for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}

	AffectRegionMasque();
	int deci=(int) (((float)(fin_i-debut_i))/12.38);// / 16.5);// /23.45);// 
	int deciprim=(int) (((float)(fin_i-debut_i))/13);// /17.33);// /25.8);// 
	int decj=(int) (((float)(fin_j-debut_j))/13);// /17.33);// /25.8);// 
	int decjprim=(int) (((float)(fin_j-debut_j))/12.38);// /16.5);// /23.45);// 
/*	moyenneellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+deci,debut_j+decj),Inv_Get_ij(fin_i-deciprim,fin_j-decjprim)));
	deci=(int) (((float)(fin_i-debut_i))/12.38);// /9.5);// 
	decj=(int) (((float)(fin_j-debut_j))/13);// /9.92);// 
	deciprim=(int) (((float)(fin_i-debut_i))/13);// /9.92);// 
	decjprim=(int) (((float)(fin_j-debut_j))/12.38);// /9.5);// 
*/	petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+deci,debut_j+decj),Inv_Get_ij(fin_i-deciprim,fin_j-decjprim)));

	SetCursor(LoadCursor(NULL,IDC_WAIT));
	
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			h=Inv_Get_ij(c,l);
			if(RegionEllipse.PtInRegion(h))
			{	
				Calc(p);
				decide=0;
				nonclas=1;

				for( clas=0;clas<Polymc.NbClass;clas++)
				{
					//if((clas==3)||(moyenneellipse.PtInRegion(h)&&clas==4)||(petiteellipse.PtInRegion(h)))
					//if((clas==3)||(clas==4)||(moyenneellipse.PtInRegion(h)&&clas==1)||(petiteellipse.PtInRegion(h)))
/*	attention*/		if((clas==3)||(clas==4)||(petiteellipse.PtInRegion(h)))
					{
						for(kpol=0;kpol<Polymc.Npol[clas];kpol++)
  						{ 
							ext=0;

							for(par=0;par<Polymc.NPAR;par++)
							{ 
							   k=kpol*Polymc.NPAR+par;
							   aa=Buffer[par];
							   if ((aa<Polymc.poli[clas][k])||(aa>Polymc.pols[clas][k])) {ext=1;break;}
							}//for par
							if (!ext) {/*if(petiteellipse.PtInRegion(h)&&clas==1){zone[5]++;decide=classecoul[clas]+20;}else */{zone[clas]++;decide=classecoul[clas];}nonclas=0;break;}				 
						}//for pol					
					if(!ext) break;
					}
				}//for clas
				if(nonclas)
				{
					if(AffectNonClas)
					{
						zone[ClassAffect]++;
						decide=classecoul[ClassAffect];
					}
					else 
					{
						zone[nbclas-1]++;
						decide=classecoul[nbclas-1];
					}
				}
				Image[p]=decide;
			}//if region
			else Image[p]=0;
			p++;
		} // c
	} // l
	if(TypePalette!=PALCOUL) OnCouleur();
//	OnCalculellipse();
	Invalidate(FALSE);

	CString ch="",ch0;
	for(i=0;i<Polymc.NbClass;i++)
	{
		ch0.Format("Classe %d= %d pixels\n\n",i,zone[i]);
		ch+=ch0;
	}
	if(!AffectNonClas)
	{
		ch0.Format("Non classés= %d pixels\n\n",zone[i]);
		ch+=ch0;
	}
	if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(Polymc.NbClass,!AffectNonClas);Legende->Invalidate(TRUE);}
	else
	{
		Legende->Init(Polymc.NbClass,!AffectNonClas);
		Legende->ShowWindow(SW_SHOW);
		Legende->Invalidate(TRUE);
	}
	if(zone[1]>37) 
		ch+="MANQUE\n";//80
	if(zone[2]>=50) 
		ch+="LISSE\n";
	if(zone[3]>12) 
		ch+="BOURGEON\n";
	if (zone[4]>=20)
		ch+="ARRACHEMENT\n";	
	if((zone[1]<=37)&&(zone[2]<50)&&(zone[3]<=12)&&(zone[4]<20)) 
		ch+="BONNE";
	MessageBox(ch);
//	OnMasqueAnnuler();
	Polymc.Fin();
	SetCursor(LoadCursor(0,IDC_ARROW));
	delete classecoul;
	delete zone;
}

void CTraitImageWnd::OnUpdateDcisionDcision3masquessansrech(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision && Masque);	
}
void CTraitImageWnd::OnDcisionParzenDcisionenregistre() 
{
	long p;
	byte decide;
	int l,c,i,nbclas=6,max2;
	CPoint h;
	CRgn petiteellipse;
	long *zone;
//	unsigned char *classecoul;
	float contamin;
	int deci,deciprim,decj,decjprim;
	ArreteLeContinu();

	byte* copEspClas;
	copEspClas=((CMainFrame*)GetMDIFrame())->EspClas;


//	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];

/*	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
	//	zone[i]=0;
	}
*/	CFileDialog FileData(FALSE,"*.xls",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"EXCEL (*.xls)|*.xls||",this);
	if (FileData.DoModal()!=IDOK) return;
	CStdioFile Fichier;
	CString ch;

	if (Fichier.Open(FileData.GetPathName(),CFile::modeCreate|CFile::modeWrite))
	{
		Fichier.WriteString(Titre);
		Fichier.WriteString("\n");
		ch.Format("cathode \t circularité \t homogénéité \t manque int \t lisse \t bourgeon \t arrachement \t risque \t TextureOK \n");
		Fichier.WriteString(ch);

		for(int r=0;r<NbMaxImageSeq;r++)
		{
			for(i=0;i<nbclas;i++) 
			{
				zone[i]=0;
			}
			
			

			AfficheSequence(r,r+1);
			
			contamin=MasqEucl(max2);

			deci=(int) (((float)(fin_i-debut_i))/12.38);// /9.5);
			decj=(int) (((float)(fin_j-debut_j))/13);// /9.92);
			deciprim=(int) (((float)(fin_i-debut_i))/13);// /9.92);
			decjprim=(int) (((float)(fin_j-debut_j))/12.38);// /9.5);
			petiteellipse.DeleteObject();
			petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+deci,debut_j+decj),Inv_Get_ij(fin_i-deciprim,fin_j-decjprim)));
			////luminance moyenne
		/*	long moy=0,nbpoint=0;
			for (l=debut_j;l<=fin_j;l++)
			{
				p=Larg*l+debut_i;
			
				for (c=debut_i;c<=fin_i;c++)
				{
					h=Inv_Get_ij(c,l);
					if(petiteellipse.PtInRegion(h))
					{	
						moy+=Image[p];
						nbpoint++;
					}//if region
					p++;
				} // c
			} // l
			moy/=nbpoint;
			if((moy>143)&&(moy<150))
			{
				moy-=143;
				for (l=debut_j;l<=fin_j;l++)
				{
					p=Larg*l+debut_i;
				
					for (c=debut_i;c<=fin_i;c++)
					{
						h=Inv_Get_ij(c,l);
						if(RegionEllipse.PtInRegion(h))
						{	
							if(Image[p]>=(unsigned char)moy) Image[p]-=(unsigned char)moy;
						}//if region
						p++;
					} // c
				} // l
			}
		*/	//////////////////////////
			for (i=0;i<NbParam;i++) 
			{
				CalculParam(i);
			}
			SetCursor(LoadCursor(NULL,IDC_WAIT));
			
			for (l=debut_j;l<=fin_j;l++)
			{
				p=Larg*l+debut_i;
			
				for (c=debut_i;c<=fin_i;c++)
				{
					h=Inv_Get_ij(c,l);
					if(RegionEllipse.PtInRegion(h))
					{	
						Calc(p);
						decide=0;
						decide=copEspClas[Buffer[0]+256*Buffer[1]+65536*Buffer[2]];
						if(decide==0||(!petiteellipse.PtInRegion(h)&&(decide==2||decide==3))) decide=1;
						decide--;
					//	Image[p]=classecoul[decide];
						zone[decide]++;
					}//if region
					else Image[p]=0;


					p++;
				} // c
			} // l
			if(TypePalette!=PALCOUL) OnCouleur();
		//	OnCalculellipse();
			Invalidate(FALSE);
			ch.Format("%d \t %2.2f \t %d \t %d \t %d \t %d \t %d \t %d\t %d\n",r+1,contamin,max2,zone[1],zone[2],zone[3],zone[4],zone[5],zone[0]);
			Fichier.WriteString(ch);
			if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(nbclas,FALSE);Legende->Invalidate(TRUE);}
			else
			{
				Legende->Init(nbclas,FALSE);
				Legende->ShowWindow(SW_SHOW);
				Legende->Invalidate(TRUE);
			}
			
			OnMasqueAnnuler();
		
	}
		Fichier.Close();
	}
	SetCursor(LoadCursor(0,IDC_ARROW));
//	delete classecoul;	
	delete zone;
}

void CTraitImageWnd::OnUpdateDcisionParzenDcisionenregistre(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision&&SEQ_EXISTE);
}

void CTraitImageWnd::OnDcisionDecisionparzencircularit() 
{
	SetCursor(LoadCursor(NULL,IDC_WAIT));
	long p;
	int l,c,i,nbclas=6,max2;
	CPoint h;
	CRgn petiteellipse;
	long *zone;
	unsigned char *classecoul;
	float contamin;
	byte decide;
	ArreteLeContinu();
	byte* copEspClas;
	copEspClas=((CMainFrame*)GetMDIFrame())->EspClas;
	if(!AffectNonClas) nbclas++;
	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];//////////////////

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
		zone[i]=0;
	}
	 contamin=MasqEucl(max2); 
	 for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}
	
	AffectRegionMasque();
	int deci=(int) (((float)(fin_i-debut_i))/12.38);
	int deciprim=(int) (((float)(fin_i-debut_i))/13);
	int decj=(int) (((float)(fin_j-debut_j))/13);
	int decjprim=(int) (((float)(fin_j-debut_j))/12.38);
	
	petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+deci,debut_j+decj),Inv_Get_ij(fin_i-deciprim,fin_j-decjprim)));

	SetCursor(LoadCursor(NULL,IDC_WAIT));
	
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			h=Inv_Get_ij(c,l);
			if(RegionEllipse.PtInRegion(h))
			{	
				Calc(p);
				decide=copEspClas[Buffer[0]+256*Buffer[1]+65536*Buffer[2]];
			//	if(Image[p]==0) Image[p]=100;//////////////ATTENTION
				if(!petiteellipse.PtInRegion(h)&&(decide==2||decide==3||decide==0)) decide=1;
				if(decide==0)
				{
					if(AffectNonClas) decide=1;
					else decide=nbclas;
				}
				decide--;
				Image[p]=classecoul[decide];
				zone[decide]++;
			}//if region
			else Image[p]=0;

			p++;
		} // c
	} // l
	if(TypePalette!=PALCOUL) OnCouleur();
//	OnCalculellipse();
	Invalidate(FALSE);

	CString ch="",ch0;
	for(i=0;i<nbclas;i++)
	{
		ch0.Format("Classe %d= %d pixels\n\n",i,zone[i]);
		ch+=ch0;
	}
	if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(nbclas,!AffectNonClas);Legende->Invalidate(TRUE);}
	else
	{
		Legende->Init(nbclas,!AffectNonClas);
		Legende->ShowWindow(SW_SHOW);
		Legende->Invalidate(TRUE);
	}
	if(zone[1]>37) 
		ch+="MANQUE\n";//80
	if(zone[2]>=50) 
		ch+="LISSE\n";
	if(zone[3]>12) 
		ch+="BOURGEON\n";
	if (zone[4]>=20)
		ch+="ARRACHEMENT\n";	
	if((zone[1]<=37)&&(zone[2]<50)&&(zone[3]<=12)&&(zone[4]<20)) 
		ch+="BONNE";
	MessageBox(ch);

	SetCursor(LoadCursor(0,IDC_ARROW));
	delete classecoul;
	delete zone;	
}

void CTraitImageWnd::OnUpdateDcisionDecisionparzencircularit(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision);
}

void CTraitImageWnd::OnDcisionDecision3masquesenregist() 
{
	CPolyMc Polymc;
	long p,kpol,clas,k;
	int l,c,par,i,nbclas,max2;
	BYTE ext,decide,aa,nonclas;
	CPoint h;
	CRgn petiteellipse,moyenneellipse;
	long *zone;
	unsigned char *classecoul;
	float contamin;
	int deci,deciprim,decj,decjprim;
	ArreteLeContinu();
	if(!AffectNonClas||ClassAffect!=0) 
	{
	//	nbclas++;
		MessageBox("Les non-classés doivent être affecté à la classe 0","Stop",MB_OK);
		return;
	}	
	Polymc.LirePol(NomFichPMC);	
	Polymc.AffPol();

	nbclas=Polymc.NbClass;

	

	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
	//	zone[i]=0;
	}
	CFileDialog FileData(FALSE,"*.xls",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"EXCEL (*.xls)|*.xls||",this);
	if (FileData.DoModal()!=IDOK) return;
	CStdioFile Fichier;
	CString ch;

	if (Fichier.Open(FileData.GetPathName(),CFile::modeCreate|CFile::modeWrite))
	{
		Fichier.WriteString(Titre);
		Fichier.WriteString("\n");
		ch.Format("cathode \t circularité \t manque int \t lisse \t bourgeon \t arrachement \t TextureOK \n");
		Fichier.WriteString(ch);

		for(int r=0;r<NbMaxImageSeq;r++)
		{
			for(i=0;i<nbclas;i++) 
			{
				zone[i]=0;
			}
			
			

			AfficheSequence(r,r+1);
			
			contamin=MasqEucl(max2);//!MasquCGQuartier();//OnCalculellipseAutom();	
			
			for (i=0;i<NbParam;i++) 
			{
				CalculParam(i);
			}
			
			/////////////////////////////////////////
			//OnMasqueRechercheAutom(10);
		/*	fin_j=debut_j+258;//8  //////////////////////à remettre pour le masque fixe
			fin_i=debut_i+258;//8*/
//			AffectRegionMasque();



	/*		deci=(int) (((float)(fin_i-debut_i))/ 16.5);// /23.45);
			deciprim=(int) (((float)(fin_i-debut_i))/17.33);// /25.8);
			decj=(int) (((float)(fin_j-debut_j))/17.33);// /25.8);
			decjprim=(int) (((float)(fin_j-debut_j))/16.5);// /23.45);

			moyenneellipse.DeleteObject();
			moyenneellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+deci,debut_j+decj),Inv_Get_ij(fin_i-deciprim,fin_j-decjprim)));
		*/	deci=(int) (((float)(fin_i-debut_i))/12.38);// /9.5);
			decj=(int) (((float)(fin_j-debut_j))/13);// /9.92);
			deciprim=(int) (((float)(fin_i-debut_i))/13);// /9.92);
			decjprim=(int) (((float)(fin_j-debut_j))/12.38);// /9.5);
			petiteellipse.DeleteObject();
			petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+deci,debut_j+decj),Inv_Get_ij(fin_i-deciprim,fin_j-decjprim)));

	/*		moyenneellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+11,debut_j+10),Inv_Get_ij(fin_i-10,fin_j-11)));
			petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+27,debut_j+26),Inv_Get_ij(fin_i-26,fin_j-27)));
		*/	SetCursor(LoadCursor(NULL,IDC_WAIT));
			
			for (l=debut_j;l<=fin_j;l++)
			{
				p=Larg*l+debut_i;
			
				for (c=debut_i;c<=fin_i;c++)
				{
					h=Inv_Get_ij(c,l);
					if(RegionEllipse.PtInRegion(h))
					{	
						Calc(p);
						decide=0;
						nonclas=1;

						for( clas=0;clas<Polymc.NbClass;clas++)
						{
						//	if((clas==3)||(moyenneellipse.PtInRegion(h)&&clas==4)||(petiteellipse.PtInRegion(h)))
						//	if((clas==3))||(clas==4)||(moyenneellipse.PtInRegion(h)&&clas==1)||(petiteellipse.PtInRegion(h)))
/*	attention		*/		if((clas==3)||(clas==4)||(petiteellipse.PtInRegion(h)))
							{
								for(kpol=0;kpol<Polymc.Npol[clas];kpol++)
  								{ 
									ext=0;

									for(par=0;par<Polymc.NPAR;par++)
									{ 
									   k=kpol*Polymc.NPAR+par;
									   aa=Buffer[par];
									   if ((aa<Polymc.poli[clas][k])||(aa>Polymc.pols[clas][k])) {ext=1;break;}
									}//for par
									if (!ext) {zone[clas]++;decide=classecoul[clas];nonclas=0;break;}				 
								}//for pol					
							if(!ext) break;
							}
						}//for clas
						if(nonclas)
						{
							if(AffectNonClas)
							{
								zone[ClassAffect]++;
								decide=classecoul[ClassAffect];
							}
							else 
							{
								zone[nbclas-1]++;
								decide=classecoul[nbclas-1];
							}
						}
						Image[p]=decide;
					}//if region
					else Image[p]=0;
					p++;
				} // c
			} // l
			if(TypePalette!=PALCOUL) OnCouleur();
		//	OnCalculellipse();
			Invalidate(FALSE);
			ch.Format("%d \t %2.2f \t %d \t %d \t %d \t %d \t %d\n",r+1,contamin,zone[1],zone[2],zone[3],zone[4],zone[0]);
			Fichier.WriteString(ch);
			if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(Polymc.NbClass,!AffectNonClas);Legende->Invalidate(TRUE);}
			else
			{
				Legende->Init(Polymc.NbClass,!AffectNonClas);
				Legende->ShowWindow(SW_SHOW);
				Legende->Invalidate(TRUE);
			}
			
			OnMasqueAnnuler();
		
	}
		Fichier.Close();
	}
	Polymc.Fin();
	SetCursor(LoadCursor(0,IDC_ARROW));
	delete classecoul;	
	delete zone;
}

void CTraitImageWnd::OnUpdateDcisionDecision3masquesenregist(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision&&SEQ_EXISTE);
}
inline void CTraitImageWnd::Calc(long t)
{
int i;
if(t<0){ MessageBox("erreur de dépassement");return;}
for (i=0;i<NbParam;i++)
	{
	Buffer[i]=VisuPar[i]->Image[t];
	}
}  // clac

LONG CTraitImageWnd::OnFinVisuMes(UINT uiParam,LONG lparam )
{
	Visumes=FALSE;
	return 0;
}

void CTraitImageWnd::OnParamtrageParamtres() 
{
 ArreteLeContinu();
 ParamCalcul=FALSE;
 int i; 

 for (i=0;i<NbParam;i++)
  {
   if (VisuPar[i]!=NULL) VisuPar[i]->ShowWindow(SW_HIDE);
   Visu=FALSE;
  }

 if (DialPara.DoModal()==IDOK)  AlloueParam();
}

void CTraitImageWnd::OnUpdateParamtrageParamtres(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideMenu);
}

void CTraitImageWnd::OnParamtrageVisualisationmesoupol() 
{
	ArreteLeContinu();
	if(Visumes) delete AffMesPol;
	else
	{
		AffMesPol=new CAfficheMesPol();
		if (AffMesPol->Init(this,NomFichMMC,/*NomFichPMC*/"",&DialPara)) 
		{
			AffMesPol->Create("mesures de "+
			NomFichMMC.Right(NomFichMMC.GetLength()-NomFichMMC.ReverseFind('\\')-1));
			AffMesPol->SetFocus();
			Visumes=TRUE;
		}
		else delete AffMesPol;
	}
}

void CTraitImageWnd::OnUpdateParamtrageVisualisationmesoupol(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideVisuMes);
	pCmdUI->SetCheck(Visumes);
}

/*
void CTraitImageWnd::OnTest() 
{
	DWORD Time;
	long tab[14];
	long s=654321,i;
	Time=GetCurrentTime();

	for(int r=0;r<100;r++)
	{
	for(i=0;i<14;i++) tab[i]=0xFF159FFF;

	for(i=0;i<5*52340;i++)
	{
		s&=tab[0];
		s&=tab[1];
		s&=tab[2];
		s&=tab[3];
		s&=tab[4];
		s&=tab[5];
		s&=tab[6];
		s&=tab[7];
		s&=tab[8];
		s&=tab[9];
		s&=tab[10];
		s&=tab[11];
		s&=tab[12];
		s&=tab[13];



	}
	}
	Time=GetCurrentTime()-Time;

	CString ch;
	ch.Format("%d",Time);
	MessageBox(ch);
}*/

void CTraitImageWnd::OnMasqueMasquefixe() 
{
	if(Larg<300 || Haut<300) {MessageBox("Traitement Impossible sur ce format d'image","Attention!",MB_OK|MB_ICONEXCLAMATION);return;}
	ArreteLeContinu();
	unsigned char *tamp;
	unsigned __int16 *templ,*tempc,*grdtempl,*grdtempc;
	long p,q;
	int i1,j1;
	long nbPoint[15],lumPoint[15];
	double lumMoy[15];
	double ecarttype,sygma,moy;
	char ch [300];


	templ=new unsigned __int16[296*308];
	tempc=new unsigned __int16[296*308];
	grdtempl=new unsigned __int16[296*308];
	grdtempc=new unsigned __int16[296*308];
	tamp=new unsigned char[Larg*Haut];

	int j,i;
	for( j=0;j<Haut;j++)
	{
		p=Larg*j;
		for( i=0;i<Larg;i++)
		{
			if(Image[p]>40) break;
			p++;
		}
		if(Image[p]>40) {j1=j;i1=i;break;}
	}
	
	for( i=0;i<i1;i++)
	{
		q=i+(j1+103)*Larg;
		for( j=0;j<50;j++)
		{
			if(Image[q]>40) break;
			q+=Larg;
		}
		if (Image[q]>40) {i1=i;break;}
	}

	i1-=12;
	j1-=24;
	debut_i=i1;
	debut_j=j1;
	fin_i=i1+295;
	fin_j=j1+307;

	unsigned __int32* debutImage32=(unsigned __int32*)(Image+debut_i+debut_j*Larg);
	unsigned __int32* debutTampon32=(unsigned __int32*)(tamp+debut_i+debut_j*Larg);
	unsigned char* debutImage=Image+debut_i+debut_j*Larg;
	unsigned char* debutTampon=tamp+debut_i+debut_j*Larg;

	GradientRapide5(debutImage32,debutTampon32);
	fenetregliscombi3(debutTampon,debutImage,templ,tempc,grdtempl,grdtempc);
p=6+6*296;
	for(j=debut_j;j<=fin_j-7;j++)
	for(i=debut_i;i<=fin_i;i++) Image[i+j*Larg]=tempc[p++]/144;
	
	Invalidate();
	MessageBox("luminance");
	BYTE l;
	for( i=0;i<15;i++) {nbPoint[i]=0;lumPoint[i]=0;}

	for( i=0;i<296;i++)
	for( j=0;j<308;j++)
	{
		l=Masquergn[i+j*296];
		if(l) {nbPoint[l-1]++;lumPoint[l-1]+=Image[i+debut_i+(j+debut_j)*Larg];}
	}

	for(i=0;i<15;i++) lumMoy[i]=(double) lumPoint[i]/nbPoint[i];

	for( i=0;i<296;i++)
	for( j=0;j<308;j++)
	{
		l=Masquergn[i+j*296];
		if(l) Image[debut_i+i+(j+debut_j)*Larg]=(unsigned char)lumMoy[l-1];
	}

	moy=0;
	for(i=0;i<15;i++) moy+=lumMoy[i];
	moy/=15;

	ecarttype=0;
	for(i=0;i<15;i++) ecarttype+=pow((lumMoy[i]-moy),2);
	ecarttype/=15;
	ecarttype=sqrt(ecarttype);

	sygma=0;
	for(i=0;i<15;i++) if((lumMoy[i]<(moy+ecarttype/2))&&
		(lumMoy[i]>(moy-ecarttype/2))) sygma++;

	sygma=(100*sygma)/15;

	sprintf(ch,"Moyenne Totale=%3.2f\nEcart Type Total=%3.2f\n%% à 1 ygma Total=%3.2f\n\nMoy centre=%3.2f\nMoy 20=%3.2f\nMoy 21=%3.2f\nMoy 30=%3.2f\nMoy 31=%3.2f\nMoy 32=%3.2f\nMoy 33=%3.2f\nMoy 40=%3.2f\nMoy 41=%3.2f\nMoy 42=%3.2f\nMoy 43=%3.2f\nMoy 44=%3.2f\nMoy 45=%3.2f\nMoy 46=%3.2f\nMoy 47=%3.2f",
		moy,ecarttype,sygma,lumMoy[0],lumMoy[1],lumMoy[2],lumMoy[3],lumMoy[4],lumMoy[5],lumMoy[6],lumMoy[7],lumMoy[8],
			lumMoy[9],lumMoy[10],lumMoy[11],lumMoy[12],lumMoy[13],lumMoy[14]);
	Invalidate(FALSE);
	MessageBox(ch,"Decoupage luminance");

	FILE *Fichier;
	CString ch0;
	int n;
	CFileDialog FileData(FALSE,"*.xls",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"EXCEL (*.xls)|*.xls||",this);
	if (FileData.DoModal()==IDOK)
	{
		if (Fichier=fopen(FileData.GetPathName(),"w"))
		{
			SetCursor(LoadCursor(0,IDC_WAIT));
			fprintf(Fichier,Titre);
			fprintf(Fichier,"\nCoefficient d'homoégnéité luminance\n");
			fprintf(Fichier,"zone\tmoyenne\tecart-type\t%% à 1 sygma\n");

			fprintf(Fichier,"Totale\t");
			sprintf(ch,"%3.2f",moy);
			ch0=ch;
			ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
			fprintf(Fichier,ch0+"\t");
			sprintf(ch,"%2.2f",ecarttype);
			ch0=ch;
			ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
			fprintf(Fichier,ch0+"\t");

			sprintf(ch,"%2.2f",sygma);
			ch0=ch;
			ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
			fprintf(Fichier,ch0+"\t\n");
			n=0;l=1;
			for(i=1;i<5;i++)
			{
				
				for(j=0;j<l;j++)
				{
					fprintf(Fichier,"%d",i);
					fprintf(Fichier,"%d\t",j);
					sprintf(ch,"%3.2f",lumMoy[n++]);
					ch0=ch;
					ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
					fprintf(Fichier,ch0+"\t\n");
				}
				l*=2;
			}
			fclose(Fichier);
			SetCursor(LoadCursor(0,IDC_ARROW));
		}
	} 

	//////////////////////////////////////////////////////

p=6+6*296;
	for(j=debut_j;j<=fin_j-7;j++)
	for(i=debut_i;i<=fin_i;i++) Image[i+j*Larg]=grdtempc[p++]/144;
	Invalidate();
	MessageBox("gradient");


	for( i=0;i<15;i++) {nbPoint[i]=0;lumPoint[i]=0;}

	for( i=0;i<296;i++)
	for( j=0;j<308;j++)
	{
		l=Masquergn[i+j*296];
		if(l) {nbPoint[l-1]++;lumPoint[l-1]+=Image[i+debut_i+(j+debut_j)*Larg];}
	}

	for(i=0;i<15;i++) lumMoy[i]=(double) lumPoint[i]/nbPoint[i];

	for( i=0;i<296;i++)
	for( j=0;j<308;j++)
	{
		l=Masquergn[i+j*296];
		if(l) Image[debut_i+i+(j+debut_j)*Larg]=(unsigned char)lumMoy[l-1];
	}

	moy=0;
	for(i=0;i<15;i++) moy+=lumMoy[i];
	moy/=15;

	ecarttype=0;
	for(i=0;i<15;i++) ecarttype+=pow((lumMoy[i]-moy),2);
	ecarttype/=15;
	ecarttype=sqrt(ecarttype);

	sygma=0;
	for(i=0;i<15;i++) if((lumMoy[i]<(moy+ecarttype/2))&&
		(lumMoy[i]>(moy-ecarttype/2))) sygma++;

	sygma=(100*sygma)/15;

	sprintf(ch,"Moyenne Totale=%3.2f\nEcart Type Total=%3.2f\n%% à 1 ygma Total=%3.2f\n\nMoy centre=%3.2f\nMoy 20=%3.2f\nMoy 21=%3.2f\nMoy 30=%3.2f\nMoy 31=%3.2f\nMoy 32=%3.2f\nMoy 33=%3.2f\nMoy 40=%3.2f\nMoy 41=%3.2f\nMoy 42=%3.2f\nMoy 43=%3.2f\nMoy 44=%3.2f\nMoy 45=%3.2f\nMoy 46=%3.2f\nMoy 47=%3.2f",
		moy,ecarttype,sygma,lumMoy[0],lumMoy[1],lumMoy[2],lumMoy[3],lumMoy[4],lumMoy[5],lumMoy[6],lumMoy[7],lumMoy[8],
			lumMoy[9],lumMoy[10],lumMoy[11],lumMoy[12],lumMoy[13],lumMoy[14]);
	Invalidate(FALSE);
	MessageBox(ch,"Decoupage gradient");
	FILE *Fichier2;
	CFileDialog FileData2(FALSE,"*.xls",NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"EXCEL (*.xls)|*.xls||",this);
	if (FileData2.DoModal()==IDOK)
	{
		if (Fichier2=fopen(FileData2.GetPathName(),"w"))
		{
			SetCursor(LoadCursor(0,IDC_WAIT));
			fprintf(Fichier2,Titre);
			fprintf(Fichier2,"\nCoefficient d'homoégnéité gradient\n");
			fprintf(Fichier2,"zone\tmoyenne\tecart-type\t%% à 1 sygma\n");

			fprintf(Fichier2,"Totale\t");
			sprintf(ch,"%3.2f",moy);
			ch0=ch;
			ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
			fprintf(Fichier2,ch0+"\t");
			sprintf(ch,"%2.2f",ecarttype);
			ch0=ch;
			ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
			fprintf(Fichier2,ch0+"\t");

			sprintf(ch,"%2.2f",sygma);
			ch0=ch;
			ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
			fprintf(Fichier2,ch0+"\t\n");
			n=0;l=1;
			for(i=1;i<5;i++)
			{
				
				for(j=0;j<l;j++)
				{
					fprintf(Fichier2,"%d",i);
					fprintf(Fichier2,"%d\t",j);
					sprintf(ch,"%3.2f",lumMoy[n++]);
					ch0=ch;
					ch0=ch0.Left(ch0.Find("."))+","+ch0.Right(ch0.GetLength()-ch0.Find(".")-1);
					fprintf(Fichier2,ch0+"\t\n");
				}
				l*=2;
			}
			fclose(Fichier2);
			SetCursor(LoadCursor(0,IDC_ARROW));
		}
	} 
}

void CTraitImageWnd::OnTraitementMax121232b() 
{
//	char ch[20];
	int q,l;
	int i,j;
	unsigned char classe[256],classe1[256],classe2[256],classe3[256];
	unsigned char Max,Maxc,Maxc1,Maxc2,Maxc3,courant,courant1,courant2,courant3,ancien,ancien1,ancien2,ancien3;

	unsigned char *maximuml;
	unsigned char *maximumc;
	
	unsigned __int32 *Image32,*maximuml32;
	unsigned __int32 courant32,ancien32;
	maximuml=new unsigned char[Larg*Haut];
	maximumc=new unsigned char[Larg*Haut];

	ArreteLeContinu();
	Memorise();

	for( i=0;i<64;i++) 
	{
		classe[i]=classe[255-i]=classe[127-i]=classe[128+i]=0;
		classe2[i]=classe2[255-i]=classe2[127-i]=classe2[128+i]=0;
		classe3[i]=classe3[255-i]=classe3[127-i]=classe3[128+i]=0;
		classe1[i]=classe1[255-i]=classe1[127-i]=classe1[128+i]=0;
	}
	//traitement ligne
	
	DWORD Time=GetCurrentTime();
	Image32=(unsigned __int32*) Image;
//	for(int t=0;t<1000;t++)
	{
		q=l=0;
		//////////////////////////////lignes
		for( j=0;j<Haut;j++)
		{
			Max=0;
			for( i=0;i<3;i++)
			{
				courant32=Image32[q/4];
				courant=(unsigned char)(courant32&0xff);
				classe[courant]++;
				Max=(((courant) > (Max)) ? (courant) : (Max));
				//Max=__max(courant,Max);
				maximuml[q]=Max;
				
				courant=(unsigned char)((courant32&0xff00)>>8);
				classe[courant]++;
				Max=(((courant) > (Max)) ? (courant) : (Max));
				//Max=__max(courant,Max);
				maximuml[q+1]=Max;
				
				courant=(unsigned char)((courant32&0xff0000)>>16);
				classe[courant]++;
				Max=(((courant) > (Max)) ? (courant) : (Max));
				//Max=__max(courant,Max);
				maximuml[q+2]=Max;

				courant=(unsigned char)(courant32>>24);
				classe[courant]++;
				Max=(((courant) > (Max)) ? (courant) : (Max));
				//Max=__max(courant,Max);
				maximuml[q+3]=Max;
				q+=4;
				
			}//for i 0->11
			
			for(i=12;i<Larg/4;i++)
			{
				courant32=Image32[q/4];
				ancien32=Image32[l/4];
				courant=(unsigned char)(courant32&0xff);
				ancien=(unsigned char)(ancien32&0xff);
				classe[courant]++;
				Max=(((courant) > (Max)) ? (courant) : (Max));
				//if(courant>Max) Max=courant;
				classe[ancien]--;
				if(ancien==Max) 
				{
					while(!classe[Max])
					{
						Max--;
					}
				}
				maximuml[q]=Max;
				
				courant=(unsigned char)((courant32&0xff00)>>8);
				ancien=(unsigned char)((ancien32&0xff00)>>8);
				classe[courant]++;
				Max=(((courant) > (Max)) ? (courant) : (Max));
				//if(courant>Max) Max=courant;
				classe[ancien]--;
				if(ancien==Max) 
				{
					while(!classe[Max])
					{
						Max--;
					}
				}
				maximuml[q+1]=Max;

				courant=(unsigned char)((courant32&0xff0000)>>16);
				ancien=(unsigned char)((ancien32&0xff0000)>>16);
				classe[courant]++;
				Max=(((courant) > (Max)) ? (courant) : (Max));
				//if(courant>Max) Max=courant;
				classe[ancien]--;
				if(ancien==Max) 
				{
					while(!classe[Max])
					{
						Max--;
					}
				}
				maximuml[q+2]=Max;

				courant=(unsigned char)(courant32>>24);
				ancien=(unsigned char)(ancien32>>24);
				classe[courant]++;
				Max=(((courant) > (Max)) ? (courant) : (Max));
				//if(courant>Max) Max=courant;
				classe[ancien]--;
				if(ancien==Max) 
				{
					while(!classe[Max])
					{
						Max--;
					}
				}
				maximuml[q+3]=Max;
				q+=4;
				l+=4;
			}//for i 12->Larg
		for(i=0;i<12;i++) classe[Image[l++]]--;
		}///for j+


		///////////////////////////////colonnes
		maximuml32=(unsigned __int32*)maximuml;
		for(i=0;i<Larg/4;i++)
		{
			q=l=i;Maxc=Maxc1=Maxc2=Maxc3=0;
			for( j=0;j<12;j++)
			{
				courant32=maximuml32[q];
				courant=(unsigned char)(courant32&0xff);
				courant1=(unsigned char)((courant32&0xff00)>>8);
				courant2=(unsigned char)((courant32&0xff0000)>>16);
				courant3=(unsigned char)(courant32>>24);
				classe[courant]++;
				classe1[courant1]++;
				classe2[courant2]++;
				classe3[courant3]++;
				Maxc=(((courant) > (Maxc)) ? (courant) : (Maxc));
				Maxc1=(((courant1) > (Maxc1)) ? (courant1) : (Maxc1));
				Maxc2=(((courant2) > (Maxc2)) ? (courant2) : (Maxc2));
				Maxc=(((courant3) > (Maxc3)) ? (courant3) : (Maxc3));
				maximumc[4*q]=Maxc;
				maximumc[4*q+1]=Maxc1;
				maximumc[4*q+2]=Maxc2;
				maximumc[4*q+3]=Maxc3;
				q+=Larg/4;
			}

			for(j=12;j<Haut;j++)
			{
				courant32=maximuml32[q];
				ancien32=maximuml32[l];
				courant=(unsigned char)(courant32&0xff);
				courant1=(unsigned char)((courant32&0xff00)>>8);
				courant2=(unsigned char)((courant32&0xff0000)>>16);
				courant3=(unsigned char)(courant32>>24);
				ancien=(unsigned char)(ancien32&0xff);
				ancien1=(unsigned char)((ancien32&0xff00)>>8);
				ancien2=(unsigned char)((ancien32&0xff0000)>>16);
				ancien3=(unsigned char)(ancien32>>24);
				classe[courant]++;
				classe1[courant1]++;
				classe2[courant2]++;
				classe3[courant3]++;
				Maxc=(((courant) > (Maxc)) ? (courant) : (Maxc));
				Maxc1=(((courant1) > (Maxc1)) ? (courant1) : (Maxc1));
				Maxc2=(((courant2) > (Maxc2)) ? (courant2) : (Maxc2));
				Maxc3=(((courant3) > (Maxc3)) ? (courant3) : (Maxc3));
				//if(courant>Maxc) Maxc=courant;
				classe[ancien]--;
				classe1[ancien1]--;
				classe2[ancien2]--;
				classe3[ancien3]--;
				if(ancien==Maxc) 
				{
					while(!classe[Maxc])
					{
						Maxc--;
					}
				}
				if(ancien1==Maxc1) 
				{
					while(!classe1[Maxc1])
					{
						Maxc1--;
					}
				}
				if(ancien2==Maxc2) 
				{
					while(!classe2[Maxc2])
					{
						Maxc2--;
					}
				}
				if(ancien3==Maxc3) 
				{
					while(!classe3[Maxc3])
					{
						Maxc3--;
					}
				}
				maximumc[4*q]=Maxc;
				maximumc[4*q+1]=Maxc1;
				maximumc[4*q+2]=Maxc2;
				maximumc[4*q+3]=Maxc3;
				q+=Larg/4;
				l+=Larg/4;
			}
			for(j=0;j<12;j++) 
			{
				classe[maximuml32[l]&0xff]--;
				classe1[(maximuml32[l]&0xff00)>>8]--;
				classe2[(maximuml32[l]&0xff0000)>>16]--;
				classe3[maximuml32[l]>>24]--;
				l+=Larg/4;
			}
		}//for i


	}/// for t


	Time=GetCurrentTime()-Time;
	CString ch;
	ch.Format("time=%d",Time);
	MessageBox(ch);

	for( i=0;i<Larg*Haut;i++) Image[i]=maximumc[i];

	delete	maximuml;
	delete maximumc;
	Invalidate(FALSE);
}

void CTraitImageWnd::OnTraitementTestcontratse() 
{
	unsigned char cont[0x10000];
	int l,c,contras,q;
	unsigned long pt=0;
	unsigned char classe[256],classe2[256];
	unsigned char Min,Minc,courant,courant2,ancien2,ancien,Max,Maxc;

	unsigned char *Minimuml,*Maximuml,*result;
//////////////////////remplissage du tablo
	for( pt=0;pt<=0xFFFF;pt++)
	{	l=pt/0x100;c=pt&0xFF;
		if (l>c) 
		{
			contras=(unsigned char)((float)255*(l-c)/(l+c));
				
		}
		else
		{
			if(l==c) contras=0;
			else contras=(unsigned char)((float)255*(c-l)/(l+c));
		}
		cont[pt]=contras;
	}
	

	
	Minimuml=new unsigned char[Larg*Haut];
	Maximuml=new unsigned char[Larg*Haut];
	result=new unsigned char[Larg*Haut];

	for(int i=0;i<64;i++) 
	{
		classe[i]=classe[255-i]=classe[127-i]=classe[128+i]=0;
		classe2[i]=classe2[255-i]=classe2[127-i]=classe2[128+i]=0;
	}
	//traitement ligne

DWORD time=GetCurrentTime();
	//////////////////////////////lignes
//for(int r=0;r<100;r++)
{
	q=l=0;
	int i,j;
	for( j=0;j<Haut;j++)
	{
		Min=255;Max=0;
		for( i=0;i<12;i++)
		{
			courant=Image[q];
			classe[courant]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant)	> (Max)) ? (courant) : (Max));
			Minimuml[q]=Min;
			Maximuml[q++]=Max;
		}//for i 0->11
		
		for(i=12;i<Larg;i++)
		{
			courant=Image[q];
			ancien=Image[l++];
			classe[courant]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant)	> (Max)) ? (courant) : (Max));
			classe[ancien]--;
			if(ancien==Min) 
			{
				while(!classe[Min])
				{
					Min++;
				}
			}
			if(ancien==Max) 
			{
				while(!classe[Max])
				{
					Max--;
				}
			}
			Maximuml[q]=Max;
			Minimuml[q++]=Min;
		}//for i LF->Largl
	for(i=0;i<12;i++) classe[Image[l++]]--;
	}///for j+

	///////////////////////////////colonnes
	for(i=0;i<Larg;i++)
	{
		q=l=i;Minc=255;Maxc=0;
		for( j=0;j<12;j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
			result[q]=cont[Maxc*0x100+Minc];
			q+=Larg;
		}

		for(j=12;j<Haut;j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			ancien=Minimuml[l];
			ancien2=Maximuml[l];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
			classe[ancien]--;
			classe2[ancien2]--;
			if(ancien==Minc) 
			{
				while(!classe[Minc])
				{
					Minc++;
				}
			}
			if(ancien2==Maxc) 
			{
				while(!classe2[Maxc])
				{
					Maxc--;
				}
			}
			result[q]=cont[Maxc*0x100+Minc];
			q+=Larg;
			l+=Larg;
		}
		for(j=0;j<12;j++) {classe[Minimuml[l]]--;classe2[Maximuml[l]]--;l+=Larg;}
	}//for i
}
	time=GetCurrentTime()-time;
	CString ch;
	ch.Format("time: %d",time);
	MessageBox(ch,"temps");
	int i;
	for (i=0;i<Larg*Haut;i++) Image[i]=result[i];

	Invalidate();
}

void CTraitImageWnd::OnTraitement5paramtrescommun() 
{
	////////////////Attention aux multiples de 4 !!!!!!!!!!!!!!!!
	/////////////////////////luminance 12*12 et gradient 12*12
	int i,q,l,j;
	unsigned char *tamp;
	unsigned __int16 *templ,*tempc,*grdtempl,*grdtempc;

	unsigned char classe[256],classe2[256];
	unsigned char Min,Minc,courant,courant2,ancien2,ancien,Max,Maxc;

	unsigned char *Minimuml,*Maximuml,*result,*ouverture,*fermeture;
	unsigned __int16 *MaxMin,temp;

	templ=new unsigned __int16[Larg*Haut];
	tempc=new unsigned __int16[Larg*Haut];
	grdtempl=new unsigned __int16[Larg*Haut];
	grdtempc=new unsigned __int16[Larg*Haut];
	tamp=new unsigned char[Larg*Haut];

	Minimuml=new unsigned char[Larg*Haut];
	Maximuml=new unsigned char[Larg*Haut];
	MaxMin=new unsigned __int16[Larg*Haut];
	result=new unsigned char[Larg*Haut];
	ouverture=new unsigned char[Larg*Haut];
	fermeture=new unsigned char[Larg*Haut];

	for( i=0;i<64;i++) 
	{
		classe[i]=classe[255-i]=classe[127-i]=classe[128+i]=0;
		classe2[i]=classe2[255-i]=classe2[127-i]=classe2[128+i]=0;
	}

DWORD time=GetCurrentTime();
//for(int r=0;r<100;r++)
{
	GradientRapide3((unsigned __int32*)Image,(unsigned __int32*)tamp);
	fenetregliscombi(tamp,Image,templ,tempc,grdtempl,grdtempc);

//////////////////////////////lignes contrast

	q=l=0;
	for(j=0;j<Haut;j++)
	{
		Min=255;Max=0;
		for( i=0;i<12;i++)
		{
			courant=Image[q];
			classe[courant]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant)	> (Max)) ? (courant) : (Max));
			//Minimuml[q]=Min;
			//Maximuml[q++]=Max;
			q++;
		}//for i 0->11
		
		for(i=12;i<Larg;i++)
		{
			courant=Image[q];
			ancien=Image[l++];
			classe[courant]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant)	> (Max)) ? (courant) : (Max));
			classe[ancien]--;
			if(ancien==Min) 
			{
				while(!classe[Min])
				{
					Min++;
				}
			}
			if(ancien==Max) 
			{
				while(!classe[Max])
				{
					Max--;
				}
			}
			Maximuml[q]=Max;
			Minimuml[q++]=Min;
		}//for i LF->Largl
	for(i=0;i<12;i++) classe[Image[l++]]--;
	}///for j+

	int j;
	///////////////////////////////colonnes
	for(i=0;i<Larg;i++)
	{
		q=l=i;Minc=255;Maxc=0;
		for(j=0;j<12;j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
		//	temp=(Minc+Maxc*0x100);
		//	result[q]=cont[temp];
		//	MaxMin[q]=temp;
			q+=Larg;
		}

		for(j=12;j<Haut;j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			ancien=Minimuml[l];
			ancien2=Maximuml[l];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
			classe[ancien]--;
			classe2[ancien2]--;
			if(ancien==Minc) 
			{
				while(!classe[Minc])
				{
					Minc++;
				}
			}
			if(ancien2==Maxc) 
			{
				while(!classe2[Maxc])
				{
					Maxc--;
				}
			}
			temp=(Minc+Maxc*0x100);
			result[q]=CONTRAST[temp];
			MaxMin[q]=temp;
			q+=Larg;
			l+=Larg;
		}
		for(j=0;j<12;j++) {classe[Minimuml[l]]--;classe2[Maximuml[l]]--;l+=Larg;}
	}//for i

	for(i=0;i<Larg*Haut;i++) Image[i]=(unsigned char) (MaxMin[i]&0xff);
	Invalidate();
	MessageBox("min");
/*	for(i=0;i<Larg*Haut;i++) Image[i]=(unsigned char) ((MaxMin[i]&0xff00)>>8);
	Invalidate();
	MessageBox("max");
*/
	/////////// 1/2 ouverture fermeture
	q=l=0;
	for( j=0;j<Haut;j++)
	{
		Min=255;Max=0;
		for( i=0;i<12;i++)
		{
			courant=(unsigned char)((MaxMin[q]&0xff00)>>8);
			courant2=(unsigned char)(MaxMin[q]&0xff);
			classe[courant]++;
			classe2[courant2]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant2)> (Max)) ? (courant2) : (Max));
		//	Minimuml[q]=Min;
		//	Maximuml[q++]=Max;
			q++;
		}//for i 0->11
		
		for(i=12;i<Larg;i++)
		{
			courant=(unsigned char)((MaxMin[q]&0xff00)>>8);
			courant2=(unsigned char)(MaxMin[q]&0xff);
			ancien=(unsigned char)((MaxMin[l]&0xff00)>>8);
			ancien2=(unsigned char)(MaxMin[l++]&0xff);
			classe[courant]++;
			classe2[courant2]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant2)	> (Max)) ? (courant2) : (Max));
			classe[ancien]--;
			classe2[ancien2]--;
			if(ancien==Min) 
			{
				while(!classe[Min])
				{
					Min++;
				}
			}
			if(ancien2==Max) 
			{
				while(!classe2[Max])
				{
					Max--;
				}
			}
			Maximuml[q]=Max;
			Minimuml[q++]=Min;
		}//for i LF->Largl
		for(i=0;i<12;i++) 
		{
			classe[(unsigned char)((MaxMin[l]&0xff00)>>8)]--;
			classe2[(unsigned char)(MaxMin[l++]&0xff)]--;
		}
	}///for j+


	///////////////////////////////colonnes
	for(i=0;i<Larg;i++)
	{
		q=l=i;Minc=255;Maxc=0;
		for(j=0;j<12;j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
		//	ouvert[q]=Minc;
		//	ferme[q]=Maxc;
			q+=Larg;
		}

		for(j=12;j<Haut;j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			ancien=Minimuml[l];
			ancien2=Maximuml[l];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
			classe[ancien]--;
			classe2[ancien2]--;
			if(ancien==Minc) 
			{
				while(!classe[Minc])
				{
					Minc++;
				}
			}
			if(ancien2==Maxc) 
			{
				while(!classe2[Maxc])
				{
					Maxc--;
				}
			}
			fermeture[q]=Minc;
			ouverture[q]=Maxc;
			q+=Larg;
			l+=Larg;
		}
		for(j=0;j<12;j++) {classe[Minimuml[l]]--;classe2[Maximuml[l]]--;l+=Larg;}
	}//for i

}
	time=GetCurrentTime()-time;
	CString ch;
	ch.Format("time: %d",time);
	MessageBox(ch,"tous");
/*	CImageWnd *Im=new CImageWnd(Larg,Haut,FLOAT32);
	for(i=0;i<Larg*Haut;i++) Im->ImageFloat[i]=tempc[i];
	Im->Create("resultat float", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());
	Invalidate();
	MessageBox("luminance");
*/
/*	for(i=0;i<Larg*Haut;i++) Image[i]=grdtempc[i]/144;
	Invalidate();
	MessageBox("gradient");

	for (i=0;i<Larg*Haut;i++) Image[i]=result[i];
	Invalidate();
	MessageBox("contraste");

	for (i=0;i<Larg*Haut;i++) Image[i]=ouverture[i];
	Invalidate();
	MessageBox("ouverture");

	for (i=0;i<Larg*Haut;i++) Image[i]=fermeture[i];
	Invalidate();
	MessageBox("fermeture");

	OnRestaureImage();*/
}


void CTraitImageWnd::fenetregliscombi(unsigned char* aram,unsigned char* bram,unsigned __int16* cufen,
											unsigned __int16* lufen,unsigned __int16* acufen,
											unsigned __int16* alufen)
{
	unsigned __int32 p,pth,ptl,pb;
	unsigned __int32 ph,pl,aph,apl;
	unsigned __int16 som,asom;
	unsigned __int32 *lram,*alram;
	unsigned __int32 som2,*dcufen,*dlufen,asom2,*adcufen,*adlufen;
	int l,c,c4,lf4,c2;//44 msc

	lram=(unsigned __int32 *)bram;
	alram=(unsigned __int32 *)aram;
	////////colonnes
	c4=Larg/4;lf4=3;//12/4;
	pb=0;ptl=0;pth=ptl+lf4;
	for(l=0;l<Haut;l++)
	{	som=asom=0;p=l*Larg;
		for(c=0;c<12;c++) 
		{
			som=som+bram[p];//cufen[pb]=som;
			asom=asom+aram[p++];//acufen[pb++]=som;
			pb++;
		}
		for(c=lf4;c<c4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			aph=alram[pth];apl=alram[ptl];

			som=som+unsigned __int16((ph&0xff)-(pl&0xff));
			asom=asom+unsigned __int16((aph&0xff)-(apl&0xff));
			cufen[pb]=som;
			acufen[pb++]=asom;

			som=som+unsigned __int16(((ph&0xff00)>>8)-((pl&0xff00)>>8));
			asom=asom+unsigned __int16(((aph&0xff00)>>8)-((apl&0xff00)>>8));
			cufen[pb]=som;
			acufen[pb++]=asom;

			som=som+unsigned __int16(((ph&0xff0000)>>16)-((pl&0xff0000)>>16));
			asom=asom+unsigned __int16(((aph&0xff0000)>>16)-((apl&0xff0000)>>16));
			cufen[pb]=som;
			acufen[pb++]=asom;

			som=som+unsigned __int16((ph>>24)-(pl>>24));
			asom=asom+unsigned __int16((aph>>24)-(apl>>24));
			cufen[pb]=som;
			acufen[pb++]=asom;

			pth++;ptl++;
		}//for c
		pth+=lf4;ptl+=lf4;
	}//for l

//////////////////////ligne
	dcufen=(unsigned __int32 *)cufen;
	dlufen=(unsigned __int32 *)lufen;
	adcufen=(unsigned __int32 *)acufen;
	adlufen=(unsigned __int32 *)alufen;

	c2=Larg/2;
	for(c=0;c<c2;c++)
	{	pth=c;som2=asom2=0;ptl=c;
		for(l=0;l<12;l++) 
		{	som2=som2+dcufen[pth];//attention pas de carry
			asom2=asom2+adcufen[pth];//attention pas de carry
			//dlufen[pth]=som2;
			//adlufen[pth]=asom2;
			pth+=c2;
		}//for l
		
		for(l=12;l<Haut;l++)
		{	som2=som2+(dcufen[pth]-dcufen[ptl]);
			asom2=asom2+(adcufen[pth]-adcufen[ptl]);
			dlufen[pth]=som2;
			adlufen[pth]=asom2;
			pth+=c2;ptl+=c2;
		}//for l
	}//for c
}


void CTraitImageWnd::OnTraitement5paramtresetessais5paramtrescommunmasque() 
{
	OnMasqueRechercheAutom(40);	
	debut_i-=12;
	debut_j-=24;
	fin_i+=24;
	fin_j+=24;
	AffectRegionMasque();

	int jl,ih,debutl[303],finl[303],debuth[291],finh[291];

	BOOL prem;
	CPoint h;
	int cpt=0,i1,j1;
	for(int j=0;j<Haut;j++)
	{
		prem=TRUE;
		for(int i=0;i<Larg;i++)
		{		
			h=Inv_Get_ij(i,j);
			if(prem && RegionEllipse.PtInRegion(h))
			{
				i1=i;
				prem=FALSE;
			}
			if(!prem &&!RegionEllipse.PtInRegion(h))
			{
				if((i-i1)>13) 
				{
					if(!cpt) jl=j;
				//	cpt++;
				//	Image[i1+j*Larg]=255;
				//	Image[i-1+j*Larg]=128;
					debutl[cpt]=i1;
					finl[cpt++]=i-1;
				}
				prem=TRUE;
			}
		}
	}
/*CString c;
c.Format("%d",cpt);
MessageBox(c);*/
	cpt=0;
	for(int i=0;i<Larg;i++)
	{
		prem=TRUE;
		for(int j=0;j<Haut;j++)
		{		
			h=Inv_Get_ij(i,j);
			if(prem && RegionEllipse.PtInRegion(h))
			{
				j1=j;
				prem=FALSE;
			}
			if(!prem &&!RegionEllipse.PtInRegion(h))
			{
				if((j-j1)>13)
				{
					if(!cpt) ih=i;
				//	cpt++;
				//	Image[i+j1*Larg]=255;
				//	Image[i+(j-1)*Larg]=128;
					debuth[cpt]=j1;
					finh[cpt++]=j-1;		
				}
				prem=TRUE;
			}
		}
	}
/*	c.Format("%d",cpt);
MessageBox(c);
*/
/*	int deb=jl;
	for (j=0;j<280;j++)
	{
		for(int i=debutl[j];i<Larg;i++)	Image[i+deb*Larg]=255;
		deb++;
	}
	Invalidate(FALSE);
	MessageBox("ligne");
	OnRestaureImage();
	deb=ih;
	for(i=0;i<282;i++)
	{
		for(int j=debuth[i];j<Haut;j++) Image[deb+j*Larg]=255;
		deb++;
	}
	Invalidate();
	MessageBox("colonne");
	OnRestaureImage();

*/
////////////////Attention aux multiples de 4 !!!!!!!!!!!!!!!!
	/////////////////////////luminance 12*12 et gradient 12*12
	int q,l;
	unsigned char *tamp;
	unsigned __int16 *templ,*tempc,*grdtempl,*grdtempc;

	unsigned char classe[256],classe2[256];
	unsigned char Min,Minc,courant,courant2,ancien2,ancien,Max,Maxc;

	unsigned char *Minimuml,*Maximuml,*result,*ouverture,*fermeture;
	unsigned __int16 *MaxMin,temp;

	templ=new unsigned __int16[Larg*Haut];
	tempc=new unsigned __int16[Larg*Haut];
	grdtempl=new unsigned __int16[Larg*Haut];
	grdtempc=new unsigned __int16[Larg*Haut];
	tamp=new unsigned char[Larg*Haut];

	Minimuml=new unsigned char[Larg*Haut];
	Maximuml=new unsigned char[Larg*Haut];
	MaxMin=new unsigned __int16[Larg*Haut];
	result=new unsigned char[Larg*Haut];
	ouverture=new unsigned char[Larg*Haut];
	fermeture=new unsigned char[Larg*Haut];
	int i;
	for( i=0;i<64;i++) 
	{
		classe[i]=classe[255-i]=classe[127-i]=classe[128+i]=0;
		classe2[i]=classe2[255-i]=classe2[127-i]=classe2[128+i]=0;
	}

DWORD time=GetCurrentTime();
for(int r=0;r<100;r++)
{
	GradientRapide3((unsigned __int32*)Image,(unsigned __int32*)tamp);
	fenetregliscombi(tamp,Image,templ,tempc,grdtempl,grdtempc);
//	GradientRapide4((unsigned __int32*)Image,(unsigned __int32*)tamp,debutl,finl,jl);
//	fenetregliscombi2(tamp,Image,templ,tempc,grdtempl,grdtempc,debutl,finl,debuth,finh,jl,ih);

//////////////////////////////lignes contrast

	int j;
	int	deb=jl;
	for(j=0;j<303;j++)
	{
		Min=255;Max=0;
		l=q=debutl[j]+Larg*(deb++);
		for( i=0;i<12;i++)
		{
			courant=Image[q];
			classe[courant]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant)	> (Max)) ? (courant) : (Max));
			//Minimuml[q]=Min;
			//Maximuml[q++]=Max;
			q++;
		}//for i 0->11
		
		for(i=debutl[j]+12;i<=finl[j];i++)
		{
			courant=Image[q];
			ancien=Image[l++];
			classe[courant]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant)	> (Max)) ? (courant) : (Max));
			classe[ancien]--;
			if(ancien==Min) 
			{
				while(!classe[Min])
				{
					Min++;
				}
			}
			if(ancien==Max) 
			{
				while(!classe[Max])
				{
					Max--;
				}
			}
			Maximuml[q]=Max;
			Minimuml[q++]=Min;
		}//for i LF->Largl
	for(i=0;i<12;i++) classe[Image[l++]]--;
	}///for j+

/*	for(i=0;i<Larg*Haut;i++) Image[i]=Maximuml[i];
	Invalidate();
	MessageBox("maxl");
	for(i=0;i<Larg*Haut;i++) Image[i]=Minimuml[i];
	Invalidate();
	MessageBox("minl");
*/
	///////////////////////////////colonnes
	deb=ih;
	for(i=0;i<291;i++)
	{
		q=l=Larg*debuth[i]+deb++;Minc=255;Maxc=0;
		for( j=0;j<12;j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
		//	temp=(Minc+Maxc*0x100);
		//	result[q]=cont[temp];
		//	MaxMin[q]=temp;
			q+=Larg;
		}

		for(j=debuth[i]+12;j<=finh[i];j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			ancien=Minimuml[l];
			ancien2=Maximuml[l];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
			classe[ancien]--;
			classe2[ancien2]--;
			if(ancien==Minc) 
			{
				while(!classe[Minc])
				{
					Minc++;
				}
			}
			if(ancien2==Maxc) 
			{
				while(!classe2[Maxc])
				{
					Maxc--;
				}
			}
			temp=Maxc;
			temp=(temp<<8);
			temp+=Minc;
			result[q]=CONTRAST[temp];
			MaxMin[q]=temp;
			q+=Larg;
			l+=Larg;
		}
		for(j=0;j<12;j++) {classe[Minimuml[l]]--;classe2[Maximuml[l]]--;l+=Larg;}
	}//for i
/*	for(i=0;i<Larg*Haut;i++) Image[i]=(unsigned char) (MaxMin[i]&0xff);
	Invalidate();
	MessageBox("min");
	for(i=0;i<Larg*Haut;i++) Image[i]=(unsigned char) ((MaxMin[i]&0xff00)>>8);
	Invalidate();
	MessageBox("max");
*/	/////////// 1/2 ouverture fermeture
//	q=l=0;
	deb=jl;
	for( j=0;j<303;j++)
	{
		Min=255;Max=0;
		l=q=debutl[j]+Larg*(deb++);
		for( i=0;i<12;i++)
		{
			courant=(unsigned char)((MaxMin[q]&0xff00)>>8);
			courant2=(unsigned char)(MaxMin[q]&0xff);
			classe[courant]++;
			classe2[courant2]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant2)> (Max)) ? (courant2) : (Max));
		//	Minimuml[q]=Min;
		//	Maximuml[q++]=Max;
			q++;
		}//for i 0->11
		
		for(i=debutl[j]+12;i<=finl[j];i++)
		{
			courant=(unsigned char)((MaxMin[q]&0xff00)>>8);
			courant2=(unsigned char)(MaxMin[q]&0xff);
			ancien=(unsigned char)((MaxMin[l]&0xff00)>>8);
			ancien2=(unsigned char)(MaxMin[l++]&0xff);
			classe[courant]++;
			classe2[courant2]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant2)	> (Max)) ? (courant2) : (Max));
			classe[ancien]--;
			classe2[ancien2]--;
			if(ancien==Min) 
			{
				while(!classe[Min])
				{
					Min++;
				}
			}
			if(ancien2==Max) 
			{
				while(!classe2[Max])
				{
					Max--;
				}
			}
			Maximuml[q]=Max;
			Minimuml[q++]=Min;
		}//for i LF->Largl
		for(i=0;i<12;i++) 
		{
			classe[(unsigned char)((MaxMin[l]&0xff00)>>8)]--;
			classe2[(unsigned char)(MaxMin[l++]&0xff)]--;
		}
	}///for j+


	///////////////////////////////colonnes
	deb=ih;
	for(i=0;i<291;i++)
	{
		q=l=Larg*debuth[i]+deb++;Minc=255;Maxc=0;
		for(j=0;j<12;j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
		//	ouvert[q]=Minc;
		//	ferme[q]=Maxc;
			q+=Larg;
		}

		for(j=debuth[i]+12;j<=finh[i];j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			ancien=Minimuml[l];
			ancien2=Maximuml[l];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
			classe[ancien]--;
			classe2[ancien2]--;
			if(ancien==Minc) 
			{
				while(!classe[Minc])
				{
					Minc++;
				}
			}
			if(ancien2==Maxc) 
			{
				while(!classe2[Maxc])
				{
					Maxc--;
				}
			}
			fermeture[q]=Minc;
			ouverture[q]=Maxc;
			q+=Larg;
			l+=Larg;
		}
		for(j=0;j<12;j++) {classe[Minimuml[l]]--;classe2[Maximuml[l]]--;l+=Larg;}
	}//for i

}
	time=GetCurrentTime()-time;
	CString ch;
	ch.Format("time: %d",time);
	MessageBox(ch,"tous");

	for(i=0;i<Larg*Haut;i++) Image[i]=tempc[i]/144;
	Invalidate();
	MessageBox("luminance");

	for(i=0;i<Larg*Haut;i++) Image[i]=grdtempc[i]/144;
	Invalidate();
	MessageBox("gradient");

	for (i=0;i<Larg*Haut;i++) Image[i]=result[i];
	Invalidate();
	MessageBox("contraste");

	for (i=0;i<Larg*Haut;i++) Image[i]=ouverture[i];
	Invalidate();
	MessageBox("ouverture");

	for (i=0;i<Larg*Haut;i++) Image[i]=fermeture[i];
	Invalidate();
	MessageBox("fermeture");

	OnRestaureImage();
}

void CTraitImageWnd::GradientRapide4(unsigned __int32 *lram,unsigned __int32 *bram,int* debutl,int*finl,int jl)
{
	int l,c,c4;
	unsigned long pth,ptl;
	unsigned __int32 ph,pl,gradph;
	unsigned char newph,newpl,oldph=0,oldpl=0;
	c4=Larg/4;
//	pth=0;ptl=pth+c4;
	for(l=0;l<304;l++)
	{ 
		pth=debutl[l]/4+(jl++)*Larg/4;
		ptl=pth+c4;
		for(c=debutl[l]/4;c<finl[l]/4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			_asm
			{
				mov ebx,ph
				mov ecx,pl
				bswap ecx
				bswap ebx
				mov newph,bl
				mov newpl,cl
				xor eax,eax
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov ah,oldph
				mov bl,TABABSO[eax]
				mov ax,cx
				mov ah,oldpl
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				;	ror ecx,8
				;***************
				mov al,newph
				mov oldph,al
				mov ah,newpl
				mov oldpl,ah
				bswap ebx
				mov gradph,ebx
				;xchg ah,al	;12->21
				;bswap eax	;1234->4321
			}
			bram[pth]=gradph;
			pth++;ptl++;
		}//for c
	}//for l
}

void CTraitImageWnd::fenetregliscombi2(unsigned char* aram,unsigned char* bram,unsigned __int16* cufen,
											unsigned __int16* lufen,unsigned __int16* acufen,
											unsigned __int16* alufen,int*debutl,int*finl,int*debuth,int*finh,int jl,int ih)
{
	unsigned __int32 p,pth,ptl,pb;
	unsigned __int32 ph,pl,aph,apl;
	unsigned __int16 som,asom;
	unsigned __int32 *lram,*alram;
	unsigned __int32 som2,*dcufen,*dlufen,asom2,*adcufen,*adlufen;
	int l,c,c4,lf4,c2;//44 msc

	lram=(unsigned __int32 *)bram;
	alram=(unsigned __int32 *)aram;
	////////colonnes
	c4=Larg/4;lf4=3;//12/4;
	
	for(l=0;l<304;l++)
	{	
		ptl=debutl[l]/4+jl*Larg/4;
		pth=ptl+lf4;
		pb=debutl[l]+jl*Larg;
		som=asom=0;p=debutl[l]+(jl++)*Larg;

		for(c=debutl[l]/4;c<debutl[l]/4+12;c++) 
		{
			som=som+bram[p];//cufen[pb]=som;
			asom=asom+aram[p++];//acufen[pb++]=som;
			pb++;
		}
		for(c=lf4+debutl[l]/4;c<finl[l]/4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			aph=alram[pth];apl=alram[ptl];

			som=som+unsigned __int16((ph&0xff)-(pl&0xff));
			asom=asom+unsigned __int16((aph&0xff)-(apl&0xff));
			cufen[pb]=som;
			acufen[pb++]=asom;

			som=som+unsigned __int16(((ph&0xff00)>>8)-((pl&0xff00)>>8));
			asom=asom+unsigned __int16(((aph&0xff00)>>8)-((apl&0xff00)>>8));
			cufen[pb]=som;
			acufen[pb++]=asom;

			som=som+unsigned __int16(((ph&0xff0000)>>16)-((pl&0xff0000)>>16));
			asom=asom+unsigned __int16(((aph&0xff0000)>>16)-((apl&0xff0000)>>16));
			cufen[pb]=som;
			acufen[pb++]=asom;

			som=som+unsigned __int16((ph>>24)-(pl>>24));
			asom=asom+unsigned __int16((aph>>24)-(apl>>24));
			cufen[pb]=som;
			acufen[pb++]=asom;

			pth++;ptl++;
		}//for c
	//	pth+=lf4;ptl+=lf4;
	}//for l

//////////////////////ligne
	dcufen=(unsigned __int32 *)cufen;
	dlufen=(unsigned __int32 *)lufen;
	adcufen=(unsigned __int32 *)acufen;
	adlufen=(unsigned __int32 *)alufen;

	c2=Larg/2;
	for(c=0;c<142;c++)
	{	ptl=pth=(ih++)+debuth[c*2]*c2;
		som2=asom2=0;
		for(l=0;l<12;l++) 
		{	som2=som2+dcufen[pth];//attention pas de carry
			asom2=asom2+adcufen[pth];//attention pas de carry
			//dlufen[pth]=som2;
			//adlufen[pth]=asom2;
			pth+=c2;
		}//for l
		
		for(l=12+debuth[c*2];l<finh[c*2];l++)
		{	som2=som2+(dcufen[pth]-dcufen[ptl]);
			asom2=asom2+(adcufen[pth]-adcufen[ptl]);
			dlufen[pth]=som2;
			adlufen[pth]=asom2;
			pth+=c2;ptl+=c2;
		}//for l
	}//for c

}


/*

void CTraitImageWnd::OnTest2() 
{
	int p,q,j1,i1;
	DWORD time=GetCurrentTime();
//	for(int r=0;r<100;r++)
	{
		for(int j=0;j<572;j++)
		{
			p=130+Larg*j;
			for(int i=130;i<382;i++)
			{
				if(Image[p]>40) break;
		//		else Image[p]=255;
				p++;
			}
			if(Image[p]>40) {j1=j;i1=i;break;}
		}
		
		for(int i=i1-200;i<i1-50;i++)
		{
			q=i+(j1+103)*Larg;
			for(int j=0;j<50;j++)
			{
				if(Image[q]>40) break;
			//	else Image[q]=255;
				q+=Larg;
			}
			if (Image[q]>40) {i1=i;break;}
		}
	}

	i1-=12;
	j1-=24;
	debut_i=i1;
	debut_j=j1;
	fin_i=i1+295;
	fin_j=j1+307;
	time=GetCurrentTime()-time;
	CString c;
	c.Format("%d",time);
//	Image[i1+j1*Larg]=255;
	Invalidate();
	MessageBox(c);
}
*/
void CTraitImageWnd::OnTraitement5paramtresetessaisEssaitot() 
{
	ArreteLeContinu();
	int p,j1,i1;
	int q,l;
	unsigned char *tamp;
	unsigned __int16 *templ,*tempc,*grdtempl,*grdtempc;

	unsigned char classe[256],classe2[256];
	unsigned char Min,Minc,courant,courant2,ancien2,ancien,Max,Maxc;

	unsigned char *Minimuml,*Maximuml,*result,*ouverture,*fermeture;
	unsigned __int16 *MaxMin,temp;

	templ=new unsigned __int16[296*308];
	tempc=new unsigned __int16[296*308];
	grdtempl=new unsigned __int16[296*308];
	grdtempc=new unsigned __int16[296*308];
	tamp=new unsigned char[Larg*Haut];

	Minimuml=new unsigned char[Larg*Haut];
	Maximuml=new unsigned char[Larg*Haut];
	MaxMin=new unsigned __int16[Larg*Haut];
	result=new unsigned char[Larg*Haut];
	ouverture=new unsigned char[Larg*Haut];
	fermeture=new unsigned char[Larg*Haut];

	for(int i=0;i<64;i++) 
	{
		classe[i]=classe[255-i]=classe[127-i]=classe[128+i]=0;
		classe2[i]=classe2[255-i]=classe2[127-i]=classe2[128+i]=0;
	}
	int recopdebl[304],recopfinl[304],recopdebh[292],recopfinh[292];
	int m;
	for(m=0;m<303;m++){recopdebl[m]=DEBUTL[m];recopfinl[m]=FINL[m];}
	for(m=0;m<291;m++){recopdebh[m]=DEBUTH[m];recopfinh[m]=FINH[m];}
	int i,j;
		for( j=0;j<572;j++)
		{
			p=130+Larg*j;
			for( i=130;i<382;i++)
			{
				if(Image[p]>40) break;
			//	else Image[p]=255;
				p++;
			}
			if(Image[p]>40) {j1=j;i1=i;break;}
		}
		
		for( i=i1-200;i<i1-50;i++)
		{
			q=i+(j1+103)*Larg;
			for( j=0;j<50;j++)
			{
				if(Image[q]>40) break;
		//		else Image[q]=255;
				q+=Larg;
			}
			if (Image[q]>40) {i1=i;break;}
		}
	

	i1-=12;
	j1-=24;
	debut_i=i1;
	debut_j=j1;
	fin_i=i1+295;
	fin_j=j1+307;

////////////////Attention aux multiples de 4 !!!!!!!!!!!!!!!!
	/////////////////////////luminance 12*12 et gradient 12*12
DWORD time=GetCurrentTime();
//for(int r=0;r<100;r++)
{

unsigned __int32* debutImage32=(unsigned __int32*)(Image+debut_i+debut_j*Larg);
unsigned __int32* debutTampon32=(unsigned __int32*)(tamp+debut_i+debut_j*Larg);
unsigned char* debutImage=Image+debut_i+debut_j*Larg;
unsigned char* debutTampon=tamp+debut_i+debut_j*Larg;

	GradientRapide5(debutImage32,debutTampon32);
	fenetregliscombi3(debutTampon,debutImage,templ,tempc,grdtempl,grdtempc);
//	GradientRapide4((unsigned __int32*)Image,(unsigned __int32*)tamp,debutl,finl,jl);
//	fenetregliscombi2(tamp,Image,templ,tempc,grdtempl,grdtempc,debutl,finl,debuth,finh,jl,ih);

//////////////////////////////lignes contrast
	int	deb=JL+debut_j;
	for(j=0;j<303;j++)
	{
		Min=255;Max=0;
		DEBUTL[j]=recopdebl[j]+debut_i;
		FINL[j]=recopfinl[j]+debut_i;

		l=q=DEBUTL[j]+Larg*(deb++);
		for( i=0;i<12;i++)
		{
			courant=Image[q];
			classe[courant]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant)	> (Max)) ? (courant) : (Max));
			//Minimuml[q]=Min;
			//Maximuml[q++]=Max;
			q++;
		}//for i 0->11
		
		for(i=DEBUTL[j]+12;i<=FINL[j];i++)
		{
			courant=Image[q];
			ancien=Image[l++];
			classe[courant]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant)	> (Max)) ? (courant) : (Max));
			classe[ancien]--;
			if(ancien==Min) 
			{
				while(!classe[Min])
				{
					Min++;
				}
			}
			if(ancien==Max) 
			{
				while(!classe[Max])
				{
					Max--;
				}
			}
			Maximuml[q]=Max;
			Minimuml[q++]=Min;
		}//for i LF->Largl
	for(i=0;i<12;i++) classe[Image[l++]]--;
	}///for j+

/*	for(i=0;i<Larg*Haut;i++) Image[i]=Maximuml[i];
	Invalidate();
	MessageBox("maxl");
	for(i=0;i<Larg*Haut;i++) Image[i]=Minimuml[i];
	Invalidate();
	MessageBox("minl");

*/	///////////////////////////////colonnes
	deb=IH+debut_i;
	for(i=0;i<291;i++)
	{
		Minc=255;Maxc=0;
		DEBUTH[i]=recopdebh[i]+debut_j;
		FINH[i]=recopfinh[i]+debut_j;

		q=l=Larg*DEBUTH[i]+deb++;
		for( j=0;j<12;j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
		//	temp=(Minc+Maxc*0x100);
		//	result[q]=cont[temp];
		//	MaxMin[q]=temp;
			q+=Larg;
		}

		for(j=DEBUTH[i]+12;j<=FINH[i];j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			ancien=Minimuml[l];
			ancien2=Maximuml[l];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
			classe[ancien]--;
			classe2[ancien2]--;
			if(ancien==Minc) 
			{
				while(!classe[Minc])
				{
					Minc++;
				}
			}
			if(ancien2==Maxc) 
			{
				while(!classe2[Maxc])
				{
					Maxc--;
				}
			}
			temp=(Minc+Maxc*0x100);
			result[q]=CONTRAST[temp];
			MaxMin[q]=temp;
			q+=Larg;
			l+=Larg;
		}
		for(j=0;j<12;j++) {classe[Minimuml[l]]--;classe2[Maximuml[l]]--;l+=Larg;}
	}//for i
/*	for(i=0;i<Larg*Haut;i++) Image[i]=(unsigned char) (MaxMin[i]&0xff);
	Invalidate();
	MessageBox("min");
	for(i=0;i<Larg*Haut;i++) Image[i]=(unsigned char) ((MaxMin[i]&0xff00)>>8);
	Invalidate();
	MessageBox("max");
*/	/////////// 1/2 ouverture fermeture
	q=l=0;
	deb=JL+debut_j;
	for( j=0;j<303;j++)
	{
		Min=255;Max=0;
		l=q=DEBUTL[j]+Larg*(deb++);
		for( i=0;i<12;i++)
		{
			courant=(unsigned char)((MaxMin[q]&0xff00)>>8);
			courant2=(unsigned char)(MaxMin[q]&0xff);
			classe[courant]++;
			classe2[courant2]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant2)> (Max)) ? (courant2) : (Max));
		//	Minimuml[q]=Min;
		//	Maximuml[q++]=Max;
			q++;
		}//for i 0->11
		
		for(i=DEBUTL[j]+12;i<=FINL[j];i++)
		{
			courant=(unsigned char)((MaxMin[q]&0xff00)>>8);
			courant2=(unsigned char)(MaxMin[q]&0xff);
			ancien=(unsigned char)((MaxMin[l]&0xff00)>>8);
			ancien2=(unsigned char)(MaxMin[l++]&0xff);
			classe[courant]++;
			classe2[courant2]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			Max=(((courant2)	> (Max)) ? (courant2) : (Max));
			classe[ancien]--;
			classe2[ancien2]--;
			if(ancien==Min) 
			{
				while(!classe[Min])
				{
					Min++;
				}
			}
			if(ancien2==Max) 
			{
				while(!classe2[Max])
				{
					Max--;
				}
			}
			Maximuml[q]=Max;
			Minimuml[q++]=Min;
		}//for i LF->Largl
		for(i=0;i<12;i++) 
		{
			classe[(unsigned char)((MaxMin[l]&0xff00)>>8)]--;
			classe2[(unsigned char)(MaxMin[l++]&0xff)]--;
		}
	}///for j+


	///////////////////////////////colonnes
	deb=IH+debut_i;
	for(i=0;i<291;i++)
	{
		q=l=Larg*DEBUTH[i]+deb++;Minc=255;Maxc=0;
		for(j=0;j<12;j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
		//	ouvert[q]=Minc;
		//	ferme[q]=Maxc;
			q+=Larg;
		}

		for(j=DEBUTH[i]+12;j<=FINH[i];j++)
		{
			courant=Minimuml[q];
			courant2=Maximuml[q];
			ancien=Minimuml[l];
			ancien2=Maximuml[l];
			classe[courant]++;
			classe2[courant2]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Maxc=(((courant2) > (Maxc)) ? (courant2) : (Maxc));
			classe[ancien]--;
			classe2[ancien2]--;
			if(ancien==Minc) 
			{
				while(!classe[Minc])
				{
					Minc++;
				}
			}
			if(ancien2==Maxc) 
			{
				while(!classe2[Maxc])
				{
					Maxc--;
				}
			}
			fermeture[q]=Minc;
			ouverture[q]=Maxc;
			q+=Larg;
			l+=Larg;
		}
		for(j=0;j<12;j++) {classe[Minimuml[l]]--;classe2[Maximuml[l]]--;l+=Larg;}
	}//for i

}
	time=GetCurrentTime()-time;
	CString ch;
	ch.Format("time: %d",time);
	MessageBox(ch,"tous");
p=0;
	for( j=debut_j;j<=fin_j;j++)
	for(i=debut_i;i<=fin_i;i++) Image[i+j*Larg]=tempc[p++]/144;
	
	Invalidate();
	MessageBox("luminance");

p=0;
	for( j=debut_j;j<=fin_j;j++)
	for(i=debut_i;i<=fin_i;i++) Image[i+j*Larg]=grdtempc[p++]/144;
	Invalidate();
	MessageBox("gradient");

	for (i=0;i<Larg*Haut;i++) Image[i]=result[i];
	Invalidate();
	MessageBox("contraste");

	for (i=0;i<Larg*Haut;i++) Image[i]=ouverture[i];
	Invalidate();
	MessageBox("ouverture");

	for (i=0;i<Larg*Haut;i++) Image[i]=fermeture[i];
	Invalidate();
	MessageBox("fermeture");

	OnRestaureImage();
	for( m=0;m<303;m++){DEBUTL[m]=recopdebl[m];FINL[m]=recopfinl[m];}
	for(m=0;m<291;m++){DEBUTH[m]=recopdebh[m];FINH[m]=recopfinh[m];}

}
void CTraitImageWnd::GradientRapide5(unsigned __int32 *lram,unsigned __int32 *bram)
{
	int l,c,c4,c4tot;
	unsigned long pth,ptl,pt=0;
	unsigned __int32 ph,pl,gradph;
	unsigned char newph,newpl,oldph=0,oldpl=0;
	c4tot=Larg/4;
	c4=74;
	for(l=0;l<307;l++)
	{
		pth=l*c4tot;ptl=pth+c4tot;
		for(c=0;c<c4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			_asm
			{
				mov ebx,ph
				mov ecx,pl
				bswap ecx
				bswap ebx
				mov newph,bl
				mov newpl,cl
				xor eax,eax
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,TABABSO[eax]
				mov ax,cx
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov ah,oldph
				mov bl,TABABSO[eax]
				mov ax,cx
				mov ah,oldpl
				mov cl,TABABSO[eax]
				add bl,cl
				ror ebx,8
				;	ror ecx,8
				;***************
				mov al,newph
				mov oldph,al
				mov ah,newpl
				mov oldpl,ah
				bswap ebx
				mov gradph,ebx
				;xchg ah,al	;12->21
				;bswap eax	;1234->4321
			}
			bram[pth]=gradph;
			pth++;ptl++;
		}//for c
	}//for l
}
void CTraitImageWnd::fenetregliscombi3(unsigned char* aram,unsigned char* bram,unsigned __int16* cufen,
											unsigned __int16* lufen,unsigned __int16* acufen,
											unsigned __int16* alufen)
{
	unsigned __int32 p,pth,ptl,pb;
	unsigned __int32 ph,pl,aph,apl;
	unsigned __int16 som,asom;
	unsigned __int32 *lram,*alram;
	unsigned __int32 som2,*dcufen,*dlufen,asom2,*adcufen,*adlufen;
	int l,c,c4,c4tot,lf4,c2;//44 msc

	lram=(unsigned __int32 *)bram;
	alram=(unsigned __int32 *)aram;
	////////colonnes
	c4tot=Larg/4;lf4=3;//12/4;
	c4=74;
	pb=0;
	for(l=0;l<308;l++)
	{	
		ptl=l*c4tot;pth=ptl+lf4;
		som=asom=0;p=l*Larg;
		for(c=0;c<12;c++) 
		{
			som=som+bram[p];//cufen[pb]=som;
			asom=asom+aram[p++];//acufen[pb++]=som;
			pb++;
		}
		for(c=lf4;c<c4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			aph=alram[pth];apl=alram[ptl];

			som=som+unsigned __int16((ph&0xff)-(pl&0xff));
			asom=asom+unsigned __int16((aph&0xff)-(apl&0xff));
			cufen[pb]=som;
			acufen[pb++]=asom;

			som=som+unsigned __int16(((ph&0xff00)>>8)-((pl&0xff00)>>8));
			asom=asom+unsigned __int16(((aph&0xff00)>>8)-((apl&0xff00)>>8));
			cufen[pb]=som;
			acufen[pb++]=asom;

			som=som+unsigned __int16(((ph&0xff0000)>>16)-((pl&0xff0000)>>16));
			asom=asom+unsigned __int16(((aph&0xff0000)>>16)-((apl&0xff0000)>>16));
			cufen[pb]=som;
			acufen[pb++]=asom;

			som=som+unsigned __int16((ph>>24)-(pl>>24));
			asom=asom+unsigned __int16((aph>>24)-(apl>>24));
			cufen[pb]=som;
			acufen[pb++]=asom;

			pth++;ptl++;
		}//for c
	//	pth+=lf4;ptl+=lf4;
	}//for l

//////////////////////ligne
	dcufen=(unsigned __int32 *)cufen;
	dlufen=(unsigned __int32 *)lufen;
	adcufen=(unsigned __int32 *)acufen;
	adlufen=(unsigned __int32 *)alufen;

	c2=148;
	for(c=0;c<c2;c++)
	{	pth=c;som2=asom2=0;ptl=c;
		for(l=0;l<12;l++) 
		{	som2=som2+dcufen[pth];//attention pas de carry
			asom2=asom2+adcufen[pth];//attention pas de carry
			//dlufen[pth]=som2;
			//adlufen[pth]=asom2;
			pth+=c2;
		}//for l
		
		for(l=12;l<308;l++)
		{	som2=som2+(dcufen[pth]-dcufen[ptl]);
			asom2=asom2+(adcufen[pth]-adcufen[ptl]);
			dlufen[pth]=som2;
			adlufen[pth]=asom2;
			pth+=c2;ptl+=c2;
		}//for l
	}//for c
}

void CTraitImageWnd::OnTraitementRecalage() 
{
	CTraitImageWnd *recal=new CTraitImageWnd(296,308,BMP);
	for(int i =0;i<296*308;i++)recal->Image[i]=Image[i];
	recal->Create("recalage de  "+Titre, WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());

}

void CTraitImageWnd::OnCalculellipse() 
{
	Memorise();
	int k,centre_i,centre_j,rayon,nbPoints,x,y,x1,y1;
	double Pi=acos(-1.0);
	OnMasqueRechercheAutom(10);
	rayon=135;
	debut_i-=5;
	debut_j-=5;
	fin_i=debut_i+270;
	fin_j=debut_j+270;
	centre_i=debut_i+rayon;
	centre_j=debut_j+rayon;
	nbPoints=(int)(2*Pi*rayon)/10;
	for(k=0;k<nbPoints/2;k++)
	{
		x1=centre_i+(int)(rayon*cos((2*Pi*(k+1))/nbPoints));
		y1=(centre_j+(int)(rayon*sin((2*Pi*(k+1))/nbPoints)));
		x=centre_i+(int)(rayon*cos(-(2*Pi*k)/nbPoints));
		y=(centre_j+(int)(rayon*sin(-(2*Pi*k)/nbPoints)));
		if(Image[x+Larg*y]>5||Image[x1+Larg*y1]>5) 
		{
			MessageBox("Contamination !","Défaut !",MB_ICONEXCLAMATION|MB_OK);
			break;
		}
		Image[x+Larg*y]=255;
		Image[x1+Larg*y1]=255;
	}
	OnMasqueAnnuler();
	Invalidate();
	MessageBox("ok");

	OnRestaureImage();
	
}
BOOL CTraitImageWnd::OnCalculellipseAutom()
{
	Memorise();
	int k,centre_i,centre_j,rayon,nbPoints,x,y,x1,y1;
	double Pi=acos(-1.0);
	OnMasqueRechercheAutom(10);
	rayon=135;
	debut_i-=5;
	debut_j-=5;
	fin_i=debut_i+270;
	fin_j=debut_j+270;
	centre_i=debut_i+rayon;
	centre_j=debut_j+rayon;
	nbPoints=(int)(2*Pi*rayon)/10;
	for(k=0;k<nbPoints/2;k++)
	{
		x1=centre_i+(int)(rayon*cos((2*Pi*(k+1))/nbPoints));
		y1=(centre_j+(int)(rayon*sin((2*Pi*(k+1))/nbPoints)));
		x=centre_i+(int)(rayon*cos(-(2*Pi*k)/nbPoints));
		y=(centre_j+(int)(rayon*sin(-(2*Pi*k)/nbPoints)));
		if(Image[x+Larg*y]>5||Image[x1+Larg*y1]>5) 
		{
			//MessageBox("Contamination !","Défaut !",MB_ICONEXCLAMATION|MB_OK);
			OnMasqueAnnuler();
			OnRestaureImage();
			return TRUE;
		}
		Image[x+Larg*y]=255;
		Image[x1+Larg*y1]=255;
	}
	OnMasqueAnnuler();
	OnRestaureImage();
	return FALSE;
}
/*
void CTraitImageWnd::OnTestminmax() 
{
	////////////////attention au nb de menus
	int l,q,k;
	long p;
	unsigned char * TabMax,Max,Maxg,Maxb,*ImAdd;
	int *PmaxX,*PmaxY,pgx,pgy,pbx,pby;
	TabMax=new unsigned char [Larg*Haut];
	PmaxX=new int[Larg*Haut];
	PmaxY=new int[Larg*Haut];
	DWORD time;

	q=0;
	for(int i=0;i<Larg;i++)
	{
		TabMax[q]=0;
		PmaxX[q]=11;
		PmaxY[q++]=11;
		
	}
	l=Larg;
	for(int j=1;j<Haut;j++)
	{
		TabMax[l]=0;
		PmaxX[l]=11;
		PmaxY[l]=11;
		l+=Larg;
	}

time = GetCurrentTime();
for(int r=0;r<100;r++)
{
	q=Larg;
	l=0;
	for( j=1;j<Haut;j++)
	{
		Max=TabMax[q++];
		//l=1+Larg*(j-1);q=l+Larg;
		l++;
		for(int i=1;i<Larg;i++)
		{
			
			pgx=PmaxX[q-1];
			pgy=PmaxY[q-1];
			pbx=PmaxX[l];
			pby=PmaxY[l];
			Maxg=Max;
			Maxb=TabMax[l];
			Max=Image[q];
		//	Max=TabMax+q;
			if(pgx&&pby)//si aucun sur le bord
			{
				if(Maxg>=Maxb)//si gauche >= Haut
				{
						/*if(ImTemp>=Maxg) //si courant >= gauche
						{
							Max=ImTemp;
							PmaxX[q]=11;
							PmaxY[q]=11;
						}
						else
						{
							Max=Maxg;
							PmaxX[q]=pgx-1;
							PmaxY[q]=pgy;
						}*/

/*
					((Max>=Maxg)?(PmaxX[q]=11,PmaxY[q]=11):(Max=Maxg,PmaxX[q]=pgx-1,PmaxY[q]=pgy));
				}
				else
					((Max>=Maxb)?(PmaxX[q]=11,PmaxY[q]=11):(Max=Maxb,PmaxX[q]=pbx,PmaxY[q]=pby-1));
				/*if(ImTemp>=Maxb)//si courant >= haut avec haut>gauche
				{
					Max=ImTemp;
					PmaxX[q]=11;
					PmaxY[q]=11;
				}
				else
				{
					Max=Maxb;
					PmaxX[q]=pbx;
					PmaxY[q]=pby-1;
				}*/
/*			}
			else if(!pgx && pby) // si gauche sur le bord
			{
				PmaxX[q]=11;
				ImAdd=Image+q-1;
				for(k=10;k>=0;k--)
				{
					((*ImAdd>Max)?(Max=*ImAdd,PmaxX[q]=k):(0));
					/*if(*ImAdd>Max)
					{
						Max=*ImAdd;
						PmaxX[q]=k;
					}*/
/*					ImAdd--;
				}

				((Max>=Maxb)?(PmaxY[q]=11):(Max=Maxb,PmaxX[q]=pbx,PmaxY[q]=pby-1));
			/*	if(Max>=Maxb)
				{
					PmaxY[q]=11;
				}
				else
				{
					Max=Maxb;
					PmaxX[q]=pbx;
					PmaxY[q]=pby-1;
				}*/

/*			}
			else if(pgx&&!pby)//haut sur le bord
			{
				PmaxY[q]=11;
				ImAdd=Image+q-Larg;
				for(k=10;k>=0;k--)
				{
					((*ImAdd>Max)?(Max=*ImAdd,PmaxY[q]=k):(0));
					/*if(*ImAdd>Max)
					{
						Max=*ImAdd;
						PmaxY[q]=k;
					}*/
/*					ImAdd-=Larg;
				}
				((Max>=Maxg)?(PmaxX[q]=11):(Max=Maxg,PmaxX[q]=pgx-1,PmaxY[q]=pgy));
			/*	if(Max>=Maxg)
				{
					PmaxX[q]=11;
				}
				else
				{
					Max=Maxg;
					PmaxX[q]=pgx-1;
					PmaxY[q]=pgy;
				}*/
/*			}//si seul sur le bord
			else
			{
				//Max=0;
				ImAdd=Image+q-11-11*Larg;
				for(k=0;k<12;k++)
				{
					for(p=0;p<12;p++)
					{
					/*	if(*(ImAdd+p)>=Max)
						{
							Max=*(ImAdd+p);
							PmaxY[q]=k;
							PmaxX[q]=p;
						}*/
/*						((*(ImAdd+p)>=Max) ? (Max=*(ImAdd+p),PmaxY[q]=k,PmaxX[q]=p):(0));
					}
					ImAdd+=Larg;
					
				}

			}	
			TabMax[q]=Max;
			q++;
			l++;
		}//for i
	}//for j
}		
	time = GetCurrentTime() - time;
	CString ch;
	ch.Format("%d",time);
	MessageBox(ch);
	for(j=0;j<Larg*Haut;j++)
	Image[j]=TabMax[j];
	Invalidate();
	delete TabMax;
	delete PmaxX;
	delete PmaxY;

}//test*/
///////////////////////////////////////////////////////
/////test polytopes
//////////////////////////////////////////////////////
typedef struct{
	DWORD p1;
	DWORD p2;
	DWORD p3; 
	DWORD p4;
	DWORD p5;
	DWORD p6;
/*	DWORD p7;
	DWORD p8;
	DWORD p9;
	DWORD p10;
	DWORD p11;
	DWORD p12;
	DWORD p13;
	DWORD p14;
	DWORD p15;
	DWORD p16;
	DWORD p17;
	DWORD p18;
	DWORD p19;
	DWORD p20;
	DWORD p21;
	DWORD p22;
	DWORD p23;
	DWORD p24;*/
} ELEM_POLY;

#define NBR_PIX 52337
#define TAILLE_AXE1 256
#define TAILLE_AXE2 256
#define TAILLE_AXE3 256
//#define TAILLE_AXE4 256
//#define TAILLE_AXE5 256
#define NBR_POLYTOPE 6*32

	ELEM_POLY axe1[TAILLE_AXE1];
	ELEM_POLY axe2[TAILLE_AXE2];
	ELEM_POLY axe3[TAILLE_AXE3];
//	ELEM_POLY axe4[TAILLE_AXE4];
//	ELEM_POLY axe5[TAILLE_AXE5];

	DWORD res1,res2,res3,res4,res5,res6;//,res7,res8,res9,res10,res11,res12,res13,res14,res15,res16,res17,res18,res19,res20,res21,res22,res23,res24;
	BYTE tab_classe[6]; //resultat de la classif pour tous les pixels.
	BYTE tab_polytope[NBR_POLYTOPE]; 


	void init_polytope(void)
{
	int i=0,j=0;

	srand( (unsigned)time( NULL ) );

	for( ; i< TAILLE_AXE1 ; i++ )
	{
		axe1[i].p1=0;//rand();
		axe1[i].p2=0;//rand();
		axe1[i].p3=0;//rand();
		axe1[i].p4=0;//rand();
		axe1[i].p5=0;//rand();
		axe1[i].p6=0;//rand();
	/*	axe1[i].p7=rand();
		axe1[i].p8=rand();
		axe1[i].p9=rand();
		axe1[i].p10=rand();
		axe1[i].p11=rand();
		axe1[i].p12=rand();
		axe1[i].p13=rand();
		axe1[i].p14=rand();
		axe1[i].p15=rand();
		axe1[i].p16=rand();
		axe1[i].p17=rand();
		axe1[i].p18=rand();
		axe1[i].p19=rand();
		axe1[i].p20=rand();
		axe1[i].p21=rand();
		axe1[i].p22=rand();
		axe1[i].p23=rand();
		axe1[i].p24=rand();*/
	}

	for( i=0 ; i< TAILLE_AXE2 ; i++ )
	{
		axe2[i].p1=0;//rand();
		axe2[i].p2=0;//rand();
		axe2[i].p3=0;//rand();
		axe2[i].p4=0;//rand();
		axe2[i].p5=0;//rand();
		axe2[i].p6=0;//rand();
/*		axe2[i].p7=rand();
		axe2[i].p8=rand();
		axe2[i].p9=rand();
		axe2[i].p10=rand();
		axe2[i].p11=rand();
		axe2[i].p12=rand();
		axe2[i].p13=rand();
		axe2[i].p14=rand();
		axe2[i].p15=rand();
		axe2[i].p16=rand();
		axe2[i].p17=rand();
		axe2[i].p18=rand();
		axe2[i].p19=rand();
		axe2[i].p20=rand();
		axe2[i].p21=rand();
		axe2[i].p22=rand();
		axe2[i].p23=rand();
		axe2[i].p24=rand();*/
	}

	for( i=0 ; i< TAILLE_AXE3 ; i++ )
	{
		axe3[i].p1=0;//rand();
		axe3[i].p2=0;//rand();
		axe3[i].p3=0;//rand();
		axe3[i].p4=0;//rand();
		axe3[i].p5=0;//rand();
		axe3[i].p6=0;//rand();
	/*	axe3[i].p7=rand();
		axe3[i].p8=rand();
		axe3[i].p9=rand();
		axe3[i].p10=rand();
		axe3[i].p11=rand();
		axe3[i].p12=rand();
		axe3[i].p13=rand();
		axe3[i].p14=rand();
		axe3[i].p15=rand();
		axe3[i].p16=rand();
		axe3[i].p17=rand();
		axe3[i].p18=rand();
		axe3[i].p19=rand();
		axe3[i].p20=rand();
		axe3[i].p21=rand();
		axe3[i].p22=rand();
		axe3[i].p23=rand();
		axe3[i].p24=rand();*/
	}

	/*	for( i=0 ; i< TAILLE_AXE4 ; i++ )
	{
		axe4[i].p1=rand();
		axe4[i].p2=rand();
		axe4[i].p3=rand();
		axe4[i].p4=rand();
		axe4[i].p5=rand();
		axe4[i].p6=rand();
	axe4[i].p7=rand();
		axe4[i].p8=rand();
		axe4[i].p9=rand();
		axe4[i].p10=rand();
		axe4[i].p11=rand();
		axe4[i].p12=rand();
		axe4[i].p13=rand();
		axe4[i].p14=rand();
		axe4[i].p15=rand();
		axe4[i].p16=rand();
		axe4[i].p17=rand();
		axe4[i].p18=rand();
		axe4[i].p19=rand();
		axe4[i].p20=rand();
		axe4[i].p21=rand();
		axe4[i].p22=rand();
		axe4[i].p23=rand();
		axe4[i].p24=rand();
	}

	for( i=0 ; i< TAILLE_AXE5 ; i++ )
	{
		axe5[i].p1=rand();
		axe5[i].p2=rand();
		axe5[i].p3=rand();
		axe5[i].p4=rand();
		axe5[i].p5=rand();
		axe5[i].p6=rand();
		axe5[i].p7=rand();
		axe5[i].p8=rand();
		axe5[i].p9=rand();
		axe5[i].p10=rand();
		axe5[i].p11=rand();
		axe5[i].p12=rand();
		axe5[i].p13=rand();
		axe5[i].p14=rand();
		axe5[i].p15=rand();
		axe5[i].p16=rand();
		axe5[i].p17=rand();
		axe5[i].p18=rand();
		axe5[i].p19=rand();
		axe5[i].p20=rand();
		axe5[i].p21=rand();
		axe5[i].p22=rand();
		axe5[i].p23=rand();
		axe5[i].p24=rand();
	}
*/
	for( i=0,j=0 ; i < NBR_POLYTOPE ; i++)
		tab_polytope[i]=(++j%6); // de 0..5: bonnes, 4 mauvaises, non classées.
}


void simu_poly1(void)
{


	int i=0;
	//unsigned __int16 tmp1,tmp2;
/*	unsigned char tmp1,tmp2;
	unsigned char tmp3,tmp4,tmp5;
*/
	for( ; i<NBR_PIX ; i++)
	{
/*		tmp1=tempc2[i];
		tmp2=grdtempc2[i];
		//tmp1=result1[i];	//	param sur 8 bits...
		//tmp2=tamp1[i];
		tmp3=result2[i];
		tmp4=ouverture1[i];
		tmp5=fermeture1[i];
		res1=axe1[tmp1].p1 & axe2[tmp2].p1 & axe3[tmp3].p1 & axe4[tmp4].p1 & axe5[tmp5].p1 ;
		res2=axe1[tmp1].p2 & axe2[tmp2].p2 & axe3[tmp3].p2 & axe4[tmp4].p2 & axe5[tmp5].p2 ;
		res3=axe1[tmp1].p3 & axe2[tmp2].p3 & axe3[tmp3].p3 & axe4[tmp4].p3 & axe5[tmp5].p3 ;
		res4=axe1[tmp1].p4 & axe2[tmp2].p4 & axe3[tmp3].p4 & axe4[tmp4].p4 & axe5[tmp5].p4 ;
		res5=axe1[tmp1].p5 & axe2[tmp2].p5 & axe3[tmp3].p5 & axe4[tmp4].p5 & axe5[tmp5].p5 ;
		res6=axe1[tmp1].p6 & axe2[tmp2].p6 & axe3[tmp3].p6 & axe4[tmp4].p6 & axe5[tmp5].p6 ;
		res7=axe1[tmp1].p7 & axe2[tmp2].p7 & axe3[tmp3].p7 & axe4[tmp4].p7 & axe5[tmp5].p7 ;
		res8=axe1[tmp1].p8 & axe2[tmp2].p8 & axe3[tmp3].p8 & axe4[tmp4].p8 & axe5[tmp5].p8 ;
		res9=axe1[tmp1].p9 & axe2[tmp2].p9 & axe3[tmp3].p9 & axe4[tmp4].p9 & axe5[tmp5].p9 ;
		res10=axe1[tmp1].p10 & axe2[tmp2].p10 & axe3[tmp3].p10 & axe4[tmp4].p10 & axe5[tmp5].p10 ;
		res11=axe1[tmp1].p11 & axe2[tmp2].p11 & axe3[tmp3].p11 & axe4[tmp4].p11 & axe5[tmp5].p11 ;
		res12=axe1[tmp1].p12 & axe2[tmp2].p12 & axe3[tmp3].p12 & axe4[tmp4].p12 & axe5[tmp5].p12 ;
		res13=axe1[tmp1].p13 & axe2[tmp2].p13 & axe3[tmp3].p13 & axe4[tmp4].p13 & axe5[tmp5].p13 ;
		res14=axe1[tmp1].p14 & axe2[tmp2].p14 & axe3[tmp3].p14 & axe4[tmp4].p14 & axe5[tmp5].p14 ;
		res15=axe1[tmp1].p15 & axe2[tmp2].p15 & axe3[tmp3].p15 & axe4[tmp4].p15 & axe5[tmp5].p15 ;
		res16=axe1[tmp1].p16 & axe2[tmp2].p16 & axe3[tmp3].p16 & axe4[tmp4].p16 & axe5[tmp5].p16 ;
		res17=axe1[tmp1].p17 & axe2[tmp2].p17 & axe3[tmp3].p17 & axe4[tmp4].p17 & axe5[tmp5].p17 ;
		res18=axe1[tmp1].p18 & axe2[tmp2].p18 & axe3[tmp3].p18 & axe4[tmp4].p18 & axe5[tmp5].p18 ;
		res19=axe1[tmp1].p19 & axe2[tmp2].p19 & axe3[tmp3].p19 & axe4[tmp4].p19 & axe5[tmp5].p19 ;
		res20=axe1[tmp1].p20 & axe2[tmp2].p20 & axe3[tmp3].p20 & axe4[tmp4].p20 & axe5[tmp5].p20 ;
		res21=axe1[tmp1].p21 & axe2[tmp2].p21 & axe3[tmp3].p21 & axe4[tmp4].p21 & axe5[tmp5].p21 ;
		res22=axe1[tmp1].p22 & axe2[tmp2].p22 & axe3[tmp3].p22 & axe4[tmp4].p22 & axe5[tmp5].p22 ;
		res23=axe1[tmp1].p23 & axe2[tmp2].p23 & axe3[tmp3].p23 & axe4[tmp4].p23 & axe5[tmp5].p23 ;
		res24=axe1[tmp1].p24 & axe2[tmp2].p24 & axe3[tmp3].p24 & axe4[tmp4].p24 & axe5[tmp5].p24 ;
*/	}

}


/*

void CTraitImageWnd::OnTestpolytop() 
{

	MessageBox("fichier d=6 3 param brillbord...");
		////////////////attention au nb de menus
	init_polytope();
		CPolyMc Polymc;
	long p,kpol,clas,k;
	int l,c,par,i,nbclas;
	BYTE decide,nonclas;
	CPoint h;
	CRgn petiteellipse,moyenneellipse;
	long *zone;
	unsigned char *classecoul;

	ArreteLeContinu();
	
	Polymc.LirePol(NomFichPMC);	
	Polymc.AffPol();
	nbclas=Polymc.NbClass;

	if(!AffectNonClas) 
	{
		nbclas++;
	}	

	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
		zone[i]=0;
	}
	

	ELEM_POLY * axe_courant;
	int numeropoly=0;
	int IndiceP,IndiceRang;
	for(clas=0;clas<Polymc.NbClass;clas++) 
		for(kpol=0;kpol<Polymc.Npol[clas];kpol++)
  		{ 	
			IndiceP=numeropoly/32;
			IndiceRang=numeropoly%32;
			for(par=0;par<Polymc.NPAR;par++)
				{ 
				
					k=kpol*Polymc.NPAR+par;
					if(par==0) axe_courant=axe1;
					else if(par==1) axe_courant=axe2;
					else if(par==2)axe_courant=axe3;
					switch (IndiceP)
					{
					case 0:
							for (i=Polymc.poli[clas][k];i<=Polymc.pols[clas][k];i++)
							{
								axe_courant[i].p1+=1<<IndiceRang;
							}
						break;
					case 1:
					
							for (i=Polymc.poli[clas][k];i<=Polymc.pols[clas][k];i++)
							{
								axe_courant[i].p2+=1<<IndiceRang;
							}
						break;
					case 2:
							for (i=Polymc.poli[clas][k];i<=Polymc.pols[clas][k];i++)
							{
								axe_courant[i].p3+=1<<IndiceRang;
							}
						break;
					case 3:
							for (i=Polymc.poli[clas][k];i<=Polymc.pols[clas][k];i++)
							{
								axe_courant[i].p4+=1<<IndiceRang;
							}
						break;
					case 4:
							for (i=Polymc.poli[clas][k];i<=Polymc.pols[clas][k];i++)
							{
								axe_courant[i].p5+=1<<IndiceRang;
							}
						break;
					case 5:
							for (i=Polymc.poli[clas][k];i<=Polymc.pols[clas][k];i++)
							{
								axe_courant[i].p6+=1<<IndiceRang;
							}
						break;
					}
			}
				numeropoly++;
		}
		
	CString ch;
	ch.Format("%d",numeropoly);
//	MessageBox(ch,"nb polytopes");
		////////////////////////////////////
	  for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}
	 // OnCalculellipse();
	/////////////////////////////////////////
	OnMasqueRechercheAutom(10);
	fin_j=debut_j+258;//8
	fin_i=debut_i+258;//8
	AffectRegionMasque();

	moyenneellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+11,debut_j+10),Inv_Get_ij(fin_i-10,fin_j-11)));
	petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+27,debut_j+26),Inv_Get_ij(fin_i-26,fin_j-27)));
	SetCursor(LoadCursor(NULL,IDC_WAIT));

	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			h=Inv_Get_ij(c,l);
			if(RegionEllipse.PtInRegion(h))
			{	
				Calc(p);
				decide=0;
				nonclas=1;

				res1=axe1[Buffer[0]].p1 & axe2[Buffer[1]].p1 & axe3[Buffer[2]].p1 ;
				res2=axe1[Buffer[0]].p2 & axe2[Buffer[1]].p2 & axe3[Buffer[2]].p2 ;
				res3=axe1[Buffer[0]].p3 & axe2[Buffer[1]].p3 & axe3[Buffer[2]].p3;
				res4=axe1[Buffer[0]].p4 & axe2[Buffer[1]].p4 & axe3[Buffer[2]].p4 ;
				res5=axe1[Buffer[0]].p5 & axe2[Buffer[1]].p5 & axe3[Buffer[2]].p5 ;
				res6=axe1[Buffer[0]].p6 & axe2[Buffer[1]].p6 & axe3[Buffer[2]].p6 ;
	/*
				if(moyenneellipse.PtInRegion(h)&&((res4&0xfffff800)||res5))//||res6)
				{
					zone[4]++;
					decide=classecoul[4];
				}
				else
				if((res3&0xfc000000)||(res4&0x7ff))
				{
						zone[3]++; decide=classecoul[3];
				}
				else
				if ((res3&0x2000000)&&petiteellipse.PtInRegion(h))
				{
					zone[2]++;
					decide=classecoul[2];
				}
				else
				if(((res2&0xffffe000)||(res3&0x1ffffff))&&petiteellipse.PtInRegion(h))
				{
					zone[1]++;
					decide=classecoul[1];
				}
				else
				if((res1 || (res2&0x1fff))&&petiteellipse.PtInRegion(h))
				{
					zone[0]++;
					decide=classecoul[0];
				}
				else 
				if(AffectNonClas)
				{
					zone[ClassAffect]++;
					decide=classecoul[ClassAffect];
				}
				else 
				{
					zone[nbclas-1]++;
					decide=classecoul[nbclas-1];
				}*/
				
	
				
/*				if((res1 || (res2&0x1fff))&&petiteellipse.PtInRegion(h))
				{
					zone[0]++;
					decide=classecoul[0];
				}
				else 
				if(((res2&0xffffe000)||(res3&0x1ffffff))&&petiteellipse.PtInRegion(h))
				{
					zone[1]++;
					decide=classecoul[1];
				}
				else
				if ((res3&0x2000000)&&petiteellipse.PtInRegion(h))
				{
					zone[2]++;
					decide=classecoul[2];
				}
				else
				if((res3&0xfc000000)||(res4&0x7fffff))
				{
						zone[3]++; decide=classecoul[3];
				}
				else
				if(moyenneellipse.PtInRegion(h)&&((res4&0xff800000)||res5||res6))
				{
					zone[4]++;
					decide=classecoul[4];
				}
				else 
				if(AffectNonClas)
				{
					zone[ClassAffect]++;
					decide=classecoul[ClassAffect];
				}
				else 
				{
					zone[nbclas-1]++;
					decide=classecoul[nbclas-1];
				}
				Image[p]=decide;
			}//if region
			else Image[p]=0;
			p++;
		} // c
	} // l
	if(TypePalette!=PALCOUL) OnCouleur();
//	OnCalculellipse();
	Invalidate(FALSE);

	CString ch0;ch="";
	for(i=0;i<Polymc.NbClass;i++)
	{
		ch0.Format("Classe %d= %d pixels\n\n",i,zone[i]);
		ch+=ch0;
	}
	if(!AffectNonClas)
	{
		ch0.Format("Non classés= %d pixels\n\n",zone[i]);
		ch+=ch0;
	}
	if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(Polymc.NbClass,!AffectNonClas);Legende->Invalidate(TRUE);}
	else
	{
		Legende->Init(Polymc.NbClass,!AffectNonClas);
		Legende->ShowWindow(SW_SHOW);
		Legende->Invalidate(TRUE);
	}
	if(zone[1]>150) 
		ch+="MANQUE\n";//80
	if(zone[2]>150) 
		ch+="LISSE\n";
	if(zone[3]>10) 
		ch+="BOURGEON\n";
	if (zone[4]>10)
		ch+="ARRACHEMENT\n";	
	if((zone[1]<=150)&&(zone[2]<=150)&&(zone[3]<=10)&&(zone[4]<=10)) 
		ch+="BONNE";
	MessageBox(ch);
//	OnMasqueAnnuler();
	Polymc.Fin();
	SetCursor(LoadCursor(0,IDC_ARROW));
	delete classecoul;
	delete zone;	
}*/

////////////////////////////test//////////

void CTraitImageWnd::OnTraitementCentre() 
{
/*	CString ch;
	long*Ysom,*Xsom,Nbp=0,p=0;
	long Xg=0,Yg=0;
	int rayon,nbPoints,x,y; 
	double Pi=acos(-1.0);


	Ysom=(long*) malloc( Haut*sizeof(long));
	Xsom=(long*) malloc( Larg*sizeof(long));;

	for(long j=0;j<Haut;j++)
		for(long i=0;i<Larg;i++)
		{
			if(Image[p]>5) 
			{
				Xsom[i]++;
				Ysom[j]++;
				Nbp++;
			}
			p++;
		}
	
	for(long i=1;i<Larg;i++) Xg+=i*Xsom[i];

	for(i=1;i<Haut;i++) Yg+=i*Ysom[i];
	Xg/=Nbp;
	Yg/=Nbp;
	
	rayon=140;

	nbPoints=48;//(int)(2*Pi*rayon)/10;
	double pente;
	long po =Xg+Larg*Yg;
	int* longueur=new int[nbPoints],dist;
	BOOL *stop=new BOOL[nbPoints];
	
	for(i=0;i<nbPoints;i++) {longueur[i]=0;stop[i]=FALSE;}
	int o=0;
	for(int k=-(nbPoints)/8;k<nbPoints/8;k++)
	{
		pente=tan((2*Pi*k)/nbPoints);
		dist=(int)(rayon*cos((2*Pi*k)/nbPoints));
		switch(o)
		{
			case 0:
				dist+=3;
				break;
			case 1:
				dist++;
				break;
			case 2:
				dist++;
				break;
			case 4:
				dist--;
				break;
			case 6:
				dist--;
				break;
			case 8:
				dist--;
				break;
			
			case 10:
				dist++;
				break;
			case 11:
				dist++;
				break;

		}
		o++;
		for(int l = dist;l>50;l--)
		{
			x=l;
			y=(int)(pente*l);

			p=po+x+Larg*y;
			if(!stop[k+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/8]++;}
				else stop[k+nbPoints/8]=TRUE;

			p=po-y+Larg*x;
			if(!stop[k+nbPoints/4+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/4+nbPoints/8]++;}
				else stop[k+nbPoints/4+nbPoints/8]=TRUE;

			p=po-x-Larg*y;
			if(!stop[k+nbPoints/2+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/2+nbPoints/8]++;}
				else stop[k+nbPoints/2+nbPoints/8]=TRUE;

			p=po+y-Larg*x;
			if(!stop[k+(nbPoints/4)*3+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+(nbPoints/4)*3+nbPoints/8]++;}
				else stop[k+(nbPoints/4)*3+nbPoints/8]=TRUE;
			//else break;*/

//		}
	/*	for(l = -(int)(rayon*cos((2*Pi*k)/nbPoints));l<-1;l++)
		{
			p=l+Xg+Larg*(Yg+(int)(tan((2*Pi*k)/nbPoints)*l));
			if(Image[p]<5) Image[p]=255;
			else break;
		}
		for(l = (int)(rayon*sin((2*Pi*k)/nbPoints+Pi/2));l>1;l--)
		{
			p=Xg+(int) (l*1.0/tan((2*Pi*k)/nbPoints+Pi/2))+Larg*(Yg+l);
			if(Image[p]<5) Image[p]=255;
			else break;
		}
		for(l = -(int)(rayon*sin((2*Pi*k)/nbPoints+Pi/2));l<-1;l++)
		{
			p=Xg+(int) (l*1.0/tan((2*Pi*k)/nbPoints+Pi/2))+Larg*(Yg+l);
			if(Image[p]<5) Image[p]=255;
			else break;
		}*/
/*	}


	Invalidate();

	ch.Format("Xg = %d\nYg= %d\n\n",Xg,Yg);
	CString cho;
	int moy=0;
	for(i=0;i<nbPoints-1;i++) {moy+=longueur[i];cho.Format("Longueur[%d]=%d\n",i,longueur[i]);ch+=cho;if(labs(longueur[i]-longueur[i+1])>3) {cho.Format("          delta Longueur[%d]=%d\n",i,labs(longueur[i]-longueur[i+1]));ch+=cho;}}
	moy+=longueur[i];
	moy/=nbPoints;
	cho.Format("moyenne =%d",moy);
	ch+=cho;
	MessageBox(ch,"Centre de gravité");
	debut_i=Xg-140+moy;
	debut_j=Yg-140+moy;
	fin_i=Xg+140-moy;
	fin_j=Yg+140-moy;
	Masque=TRUE;
	free(Ysom);
	free(Xsom);*/
//MasqueCG();
/*MasqCGMin();
MasqueCGQuartmoy();*/
int a;
MasqEucl(a);
//MasquCGQuartier();
}

BOOL CTraitImageWnd::MasqueCG()
{
	long*Ysom,*Xsom,Nbp=0,p=0;
	long Xg=0,Yg=0;
	int rayon,nbPoints,x,y; 
	double Pi=acos(-1.0);

	Memorise();
	Ysom=(long*) malloc( Haut*sizeof(long));
	Xsom=(long*) malloc( Larg*sizeof(long));;
	long i,j;
	for( j=0;j<Haut;j++)
		for( i=0;i<Larg;i++)
		{
			if(Image[p]>5) 
			{
				Xsom[i]++;
				Ysom[j]++;
				Nbp++;
			}
			p++;
		}
	
	for( i=1;i<Larg;i++) Xg+=i*Xsom[i];

	if(Nbp==0) return FALSE;
	for(i=1;i<Haut;i++) Yg+=i*Ysom[i];
	Xg/=Nbp;
	Yg/=Nbp;
	
	rayon=140;

	nbPoints=48;//(int)(2*Pi*rayon)/10;
	double pente;
	long po =Xg+Larg*Yg;
	int* longueur=new int[nbPoints],dist;
	BOOL *stop=new BOOL[nbPoints];
	
	for(i=0;i<nbPoints;i++) {longueur[i]=0;stop[i]=FALSE;}
	int o=0;
	for(int k=-(nbPoints)/8;k<nbPoints/8;k++)
	{
		pente=tan((2*Pi*k)/nbPoints);
		dist=(int)(rayon*cos((2*Pi*k)/nbPoints));
		switch(o)
		{
			case 0:
				dist+=3;
				break;
			case 1:
				dist++;
				break;
			case 2:
				dist++;
				break;
			case 4:
				dist--;
				break;
			case 6:
				dist--;
				break;
			case 8:
				dist--;
				break;
			
			case 10:
				dist++;
				break;
			case 11:
				dist++;
				break;

		}
		o++;
		for(int l = dist;l>50;l--)
		{
			x=l;
			y=(int)(pente*l);

			p=po+x+Larg*y;
			if(!stop[k+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/8]++;}
				else stop[k+nbPoints/8]=TRUE;

			p=po-y+Larg*x;
			if(!stop[k+nbPoints/4+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/4+nbPoints/8]++;}
				else stop[k+nbPoints/4+nbPoints/8]=TRUE;

			p=po-x-Larg*y;
			if(!stop[k+nbPoints/2+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/2+nbPoints/8]++;}
				else stop[k+nbPoints/2+nbPoints/8]=TRUE;

			p=po+y-Larg*x;
			if(!stop[k+(nbPoints/4)*3+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+(nbPoints/4)*3+nbPoints/8]++;}
				else stop[k+(nbPoints/4)*3+nbPoints/8]=TRUE;
		}
	}


	Invalidate();
	BOOL sortie=TRUE;
	int moy=0;
	int max=longueur[0],min=longueur[0];
	CString c,c0;
	c="";
	for(i=0;i<nbPoints-1;i++) {c0.Format("longueur[%d]=%d \n",i,longueur[i]);c+=c0;moy+=longueur[i];if(longueur[i]>max) max=longueur[i];if(longueur[i]<min) min=longueur[i];if(labs(longueur[i]-longueur[i+1])>3) {sortie=FALSE;}}
	if(longueur[i]>max) max=longueur[i];
	if(longueur[i]<min) min=longueur[i];
	c0.Format("longueur[%d]=%d\n\n",i,longueur[i]);
	c+=c0;
	moy+=longueur[i];
	moy/=nbPoints;
	c0.Format("moyenne=%d\n\nmax=%d\nmin=%d\n dif=%d",moy,max,min,max-min);
	c+=c0;
	debut_i=Xg-140+moy;
	debut_j=Yg-140+moy;
	fin_i=Xg+140-moy;
	fin_j=Yg+140-moy;

	Masque=TRUE;
	////////////à retirer //////////////////////
	 if(!sortie) MessageBox(c0,"problèmes de bords");
	 else MessageBox(c,"cg MOY ok");
	//////////////////////////////
	OnRestaureImage();
	free(Ysom);
	free(Xsom);
	return sortie;
}


BOOL CTraitImageWnd::MasqCGMin()
{
	long*Ysom,*Xsom,Nbp=0,p=0;
	long Xg=0,Yg=0;
	int rayon,nbPoints,x,y; 
	double Pi=acos(-1.0);

	Memorise();
	Ysom=(long*) malloc( Haut*sizeof(long));
	Xsom=(long*) malloc( Larg*sizeof(long));;

	for(long j=0;j<Haut;j++)
		for(long i=0;i<Larg;i++)
		{
			if(Image[p]>5) 
			{
				Xsom[i]++;
				Ysom[j]++;
				Nbp++;
			}
			p++;
		}
	long i;
	for( i=1;i<Larg;i++) Xg+=i*Xsom[i];

	if(Nbp==0) return FALSE;
	for(i=1;i<Haut;i++) Yg+=i*Ysom[i];
	Xg/=Nbp;
	Yg/=Nbp;
	
	rayon=140;

	nbPoints=48;//(int)(2*Pi*rayon)/10;
	double pente;
	long po =Xg+Larg*Yg;
	int* longueur=new int[nbPoints],dist;
	BOOL *stop=new BOOL[nbPoints];
	
	for(i=0;i<nbPoints;i++) {longueur[i]=0;stop[i]=FALSE;}
	int o=0;
	for(int k=-(nbPoints)/8;k<nbPoints/8;k++)
	{
		pente=tan((2*Pi*k)/nbPoints);
		dist=(int)(rayon*cos((2*Pi*k)/nbPoints));
		switch(o)
		{
			case 0:
				dist+=3;
				break;
			case 1:
				dist++;
				break;
			case 2:
				dist++;
				break;
			case 4:
				dist--;
				break;
			case 6:
				dist--;
				break;
			case 8:
				dist--;
				break;
			
			case 10:
				dist++;
				break;
			case 11:
				dist++;
				break;

		}
		o++;
		for(int l = dist;l>50;l--)
		{
			x=l;
			y=(int)(pente*l);

			p=po+x+Larg*y;
			if(!stop[k+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/8]++;}
				else stop[k+nbPoints/8]=TRUE;

			p=po-y+Larg*x;
			if(!stop[k+nbPoints/4+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/4+nbPoints/8]++;}
				else stop[k+nbPoints/4+nbPoints/8]=TRUE;

			p=po-x-Larg*y;
			if(!stop[k+nbPoints/2+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/2+nbPoints/8]++;}
				else stop[k+nbPoints/2+nbPoints/8]=TRUE;

			p=po+y-Larg*x;
			if(!stop[k+(nbPoints/4)*3+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+(nbPoints/4)*3+nbPoints/8]++;}
				else stop[k+(nbPoints/4)*3+nbPoints/8]=TRUE;
		}
	}


	Invalidate();
	BOOL sortie=TRUE;
	int max=longueur[0];
	for(i=0;i<nbPoints-1;i++) {if(longueur[i]>max) max=longueur[i];if(labs(longueur[i]-longueur[i+1])>3) {sortie=FALSE;}}
	if(longueur[i]>max) max=longueur[i];
	debut_i=Xg-140+max;
	debut_j=Yg-140+max;
	fin_i=Xg+140-max;
	fin_j=Yg+140-max;
	Masque=TRUE;
	////////////à retirer //////////////////////
	 if(!sortie) MessageBox("problèmes de bords");
	 else MessageBox("cg MIN ok");
	 //////////////////////////////
	OnRestaureImage();
	free(Ysom);
	free(Xsom);
	return sortie;
}



BOOL CTraitImageWnd::MasquCGQuartier()
{
	long*Ysom,*Xsom,Nbp=0,p=0;
	long Xg=0,Yg=0;
	int rayon,nbPoints,x,y; 
	double Pi=acos(-1.0);

	Memorise();
	Ysom=(long*) malloc( Haut*sizeof(long));
	Xsom=(long*) malloc( Larg*sizeof(long));;

	for(long j=0;j<Haut;j++)
		for(long i=0;i<Larg;i++)
		{
			if(Image[p]>5) 
			{
				Xsom[i]++;
				Ysom[j]++;
				Nbp++;
			}
			p++;
		}
	long i;
	for( i=1;i<Larg;i++) Xg+=i*Xsom[i];

	if(Nbp==0) return FALSE;
	for(i=1;i<Haut;i++) Yg+=i*Ysom[i];
	Xg/=Nbp;
	Yg/=Nbp;
	
	rayon=140;

	nbPoints=48;//(int)(2*Pi*rayon)/10;
	double pente;
	long po =Xg+Larg*Yg;
	int* longueur=new int[nbPoints],dist;
	BOOL *stop=new BOOL[nbPoints];
	
	for(i=0;i<nbPoints;i++) {longueur[i]=0;stop[i]=FALSE;}
	int o=0;
	for(int k=-(nbPoints)/8;k<nbPoints/8;k++)
	{
		pente=tan((2*Pi*k)/nbPoints);
		dist=(int)(rayon*cos((2*Pi*k)/nbPoints));
		switch(o)
		{
			case 0:
				dist+=3;
				break;
			case 1:
				dist++;
				break;
			case 2:
				dist++;
				break;
			case 4:
				dist--;
				break;
			case 6:
				dist--;
				break;
			case 8:
				dist--;
				break;
			
			case 10:
				dist++;
				break;
			case 11:
				dist++;
				break;

		}
		o++;
		for(int l = dist;l>50;l--)
		{
			x=l;
			y=(int)(pente*l);

			p=po+x+Larg*y;
			if(!stop[k+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/8]++;}
				else stop[k+nbPoints/8]=TRUE;

			p=po-y+Larg*x;
			if(!stop[k+nbPoints/4+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/4+nbPoints/8]++;}
				else stop[k+nbPoints/4+nbPoints/8]=TRUE;

			p=po-x-Larg*y;
			if(!stop[k+nbPoints/2+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/2+nbPoints/8]++;}
				else stop[k+nbPoints/2+nbPoints/8]=TRUE;

			p=po+y-Larg*x;
			if(!stop[k+(nbPoints/4)*3+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+(nbPoints/4)*3+nbPoints/8]++;}
				else stop[k+(nbPoints/4)*3+nbPoints/8]=TRUE;
		}
	}


	Invalidate();
	BOOL sortie=TRUE;
	int max,MAX,MIN;
	MAX=MIN=max=longueur[0];
	for(i=0;i<12;i++) 
	{
		if(longueur[i]>max) max=longueur[i];
		if(longueur[i]<MIN) MIN=longueur[i];
		if(labs(longueur[i]-longueur[i+1])>3) {sortie=FALSE;}
	}
	MAX=max;

	fin_i=Xg+140-max;
	max=longueur[i];
	for(i=12;i<24;i++) 
	{
		if(longueur[i]>max) max=longueur[i];
		if(longueur[i]<MIN) MIN=longueur[i];
		if(labs(longueur[i]-longueur[i+1])>3) {sortie=FALSE;}
	}
	if(max>MAX) MAX=max;

	fin_j=Yg+140-max;
	max=longueur[i];
	for(i=24;i<36;i++) 
	{
		if(longueur[i]>max) max=longueur[i];
		if(longueur[i]<MIN) MIN=longueur[i];
		if(labs(longueur[i]-longueur[i+1])>3) {sortie=FALSE;}
	}
	if(max>MAX) MAX=max;

	debut_i=Xg-140+max;
	max=longueur[i];
	for(i=36;i<47;i++) 
	{
		if(longueur[i]>max) max=longueur[i];
		if(longueur[i]<MIN) MIN=longueur[i];
		if(labs(longueur[i]-longueur[i+1])>3) {sortie=FALSE;}
	}
	if(max>MAX) MAX=max;

	debut_j=Yg-140+max;
	
	Masque=TRUE;
	if(MAX-MIN>5)
	{
		if(debut_i>Xg-140+MIN+5) debut_i=Xg-140+MIN+2;
		if(debut_j>Yg-140+MIN+5) debut_j=Yg-140+MIN+2;
		if(fin_i<Xg+140-MIN-5) fin_i=Xg+140-MIN-2;
		if(fin_j<Yg+140-MIN-5) fin_j=Yg+140-MIN-2;
	}
	////////////à retirer //////////////////////
	 if(!sortie) MessageBox("problèmes de bords");
	 else if(MAX-MIN>5) MessageBox("Masque Max-Min");
		 else MessageBox("cg QUARTIER MIN ok");
	 //////////////////////////////
	OnRestaureImage();
	free(Ysom);
	free(Xsom);
	return sortie;

}

BOOL CTraitImageWnd::MasqueCGQuartmoy()
{
	long*Ysom,*Xsom,Nbp=0,p=0;
	long Xg=0,Yg=0;
	int rayon,nbPoints,x,y; 
	double Pi=acos(-1.0);

	Memorise();
	Ysom=(long*) malloc( Haut*sizeof(long));
	Xsom=(long*) malloc( Larg*sizeof(long));;

	for(long j=0;j<Haut;j++)
		for(long i=0;i<Larg;i++)
		{
			if(Image[p]>5) 
			{
				Xsom[i]++;
				Ysom[j]++;
				Nbp++;
			}
			p++;
		}
	long i;
	for( i=1;i<Larg;i++) Xg+=i*Xsom[i];

	if(Nbp==0) return FALSE;
	for(i=1;i<Haut;i++) Yg+=i*Ysom[i];
	Xg/=Nbp;
	Yg/=Nbp;
	
	rayon=140;

	nbPoints=48;//(int)(2*Pi*rayon)/10;
	double pente;
	long po =Xg+Larg*Yg;
	int* longueur=new int[nbPoints],dist;
	BOOL *stop=new BOOL[nbPoints];
	
	for(i=0;i<nbPoints;i++) {longueur[i]=0;stop[i]=FALSE;}
	int o=0;
	for(int k=-(nbPoints)/8;k<nbPoints/8;k++)
	{
		pente=tan((2*Pi*k)/nbPoints);
		dist=(int)(rayon*cos((2*Pi*k)/nbPoints));
		switch(o)
		{
			case 0:
				dist+=3;
				break;
			case 1:
				dist++;
				break;
			case 2:
				dist++;
				break;
			case 4:
				dist--;
				break;
			case 6:
				dist--;
				break;
			case 8:
				dist--;
				break;
			
			case 10:
				dist++;
				break;
			case 11:
				dist++;
				break;

		}
		o++;
		for(int l = dist;l>50;l--)
		{
			x=l;
			y=(int)(pente*l);

			p=po+x+Larg*y;
			if(!stop[k+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/8]++;}
				else stop[k+nbPoints/8]=TRUE;

			p=po-y+Larg*x;
			if(!stop[k+nbPoints/4+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/4+nbPoints/8]++;}
				else stop[k+nbPoints/4+nbPoints/8]=TRUE;

			p=po-x-Larg*y;
			if(!stop[k+nbPoints/2+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+nbPoints/2+nbPoints/8]++;}
				else stop[k+nbPoints/2+nbPoints/8]=TRUE;

			p=po+y-Larg*x;
			if(!stop[k+(nbPoints/4)*3+nbPoints/8])
				if(Image[p]<5) {Image[p]=255;longueur[k+(nbPoints/4)*3+nbPoints/8]++;}
				else stop[k+(nbPoints/4)*3+nbPoints/8]=TRUE;
		}
	}


	Invalidate();
	BOOL sortie=TRUE;
	int moy=0;

	for(i=0;i<12;i++) 
	{
		moy+=longueur[i];
		if(labs(longueur[i]-longueur[i+1])>3) {sortie=FALSE;}
	}
	moy/=12;
	fin_i=Xg+140-moy;
	moy=0;
	for(i=12;i<24;i++) 
	{
		moy+=longueur[i];
		if(labs(longueur[i]-longueur[i+1])>3) {sortie=FALSE;}
	}
	moy/=12;
	fin_j=Yg+140-moy;
	moy=0;
	for(i=24;i<36;i++) 
	{
		moy+=longueur[i];
		if(labs(longueur[i]-longueur[i+1])>3) {sortie=FALSE;}
	}
	moy/=12;
	debut_i=Xg-140+moy;
	moy=0;
	for(i=36;i<47;i++) 
	{
		moy+=longueur[i];
		if(labs(longueur[i]-longueur[i+1])>3) {sortie=FALSE;}
	}
	moy/=12;
	debut_j=Yg-140+moy;

	
	
	
	Masque=TRUE;
	////////////à retirer //////////////////////
	 if(!sortie) MessageBox("problèmes de bords");
	 else MessageBox("cg QUARTIER MOY ok");
	 //////////////////////////////
	OnRestaureImage();
	free(Ysom);
	free(Xsom);
	return sortie;

}

float CTraitImageWnd::MasqEucl(int &max2)
{
	int *gradcirc;


	long*Ysom,*Xsom,Nbp=0,p=0;//,*longueur,*tabmoy;
	float *longueur,*tabmoy;
	long Xg=0,Yg=0,po;//,moyloc;
	float moyloc;
	int rayon,nbPoints,loc=25;
	int/*float*/ x,y; 
	double Pi=acos(-1.0);
	BOOL *stop;

	ArreteLeContinu();

	rayon=140;
	nbPoints=200;
	moyloc=0;
	longueur=new float[nbPoints];//long[nbPoints];
	tabmoy=new float[nbPoints];//long[nbPoints];
	stop=new BOOL[nbPoints];

	gradcirc=new int[nbPoints];

	Memorise();
	Ysom=(long*) malloc( Haut*sizeof(long));
	Xsom=(long*) malloc( Larg*sizeof(long));;

	for(long j=0;j<Haut;j++)
		for(long i=0;i<Larg;i++)
		{
			if(Image[p]>5) 
			{
				Xsom[i]++;
				Ysom[j]++;
				Nbp++;
			}
			p++;
		}
	long i;
	for( i=1;i<Larg;i++) Xg+=i*Xsom[i];

	if(Nbp==0) return 10.0;
	for(i=1;i<Haut;i++) Yg+=i*Ysom[i];
	Xg/=Nbp;
	Yg/=Nbp;
	
	po =Xg+Larg*Yg;
	
	for(i=0;i<nbPoints;i++) {longueur[i]=0;stop[i]=FALSE;}
	
	float pente;
	int dist,xcour;
	long p1,p2,p3,p4,p11,p22,p33,p44;
	for(int k=-(nbPoints)/8;k<nbPoints/8;k++)
	{
		///////////////////////////////////
		pente=(float)tan((2*Pi*k)/nbPoints);
		dist=(int)(rayon*cos((2*Pi*k)/nbPoints));
		///////////////////////////////////////////
		if(k==-(nbPoints)/8) dist++;
		if(k==0) dist--;
	
		for(int l = dist/*rayon*/;l>(int)((rayon-40)*cos((2*Pi*k)/nbPoints))/*105*/;l--)
		{
			//x=(float)(l*cos((2*Pi*k)/nbPoints));
			//y=(float)(l*sin((2*Pi*k)/nbPoints));
			x=l;
			y=(int)(pente*l);
			p=po+(int)x+Larg*(int)y;
			if(l==dist) {p1=p;p11=p;xcour=x;}
			if(!stop[k+nbPoints/8])
			{
				if(((Image[p]+Image[p1]+Image[p11])>=160)&&(Image[p]<240)) //
				{
					stop[k+nbPoints/8]=TRUE;
					longueur[k+nbPoints/8]=(float) sqrt((dist-xcour)*(dist-xcour)*(1+pente*pente));//rayon-l;

					gradcirc[k+nbPoints/8]=(int)Image[po+l-3+Larg*(int)(pente*(l-3))]+(int)Image[po+l-2+Larg*(int)(pente*(l-2))]+(int)Image[po+l-4+Larg*(int)(pente*(l-4))];
//					Image[p]=255;
				}
				p11=p1;
				p1=p;
			}
			else if(Image[p]<=40||Image[p]>240)
			{
				longueur[k+nbPoints/8]=(float) sqrt((dist-x)*(dist-x)*(1+pente*pente));//rayon-l;
			
	//			gradcirc[k+nbPoints/8]=(int)Image[po+l-3+Larg*(int)(pente*(l-3))];
			}
			
//
			p=po-(int)y+Larg*(int)x;
			if(l==dist){p2=p;p22=p;}
			if(!stop[k+nbPoints/4+nbPoints/8])
			{
				if(((Image[p]+Image[p2]+Image[p22])>=160)&&(Image[p]<240))
				{
					stop[k+nbPoints/4+nbPoints/8]=TRUE;
					longueur[k+nbPoints/4+nbPoints/8]=(float) sqrt((dist-xcour)*(dist-xcour)*(1+pente*pente));//rayon-l;
				
					gradcirc[k+nbPoints/4+nbPoints/8]=(int)Image[po-(int)(pente*(l-3))+Larg*(l-3)]+(int)Image[po-(int)(pente*(l-2))+Larg*(l-2)]+(int)Image[po-(int)(pente*(l-4))+Larg*(l-4)];
//					Image[p]=255;
				}
				p22=p2;
				p2=p;
			}
			else if(Image[p]<=40||Image[p]>240) 
			{
				longueur[k+nbPoints/4+nbPoints/8]=(float) sqrt((dist-x)*(dist-x)*(1+pente*pente));//rayon-l;
				
//				gradcirc[k+nbPoints/4+nbPoints/8]=(int)Image[po-(int)(pente*(l-3))+Larg*(l-3)];
			}
			
//
			p=po-(int)x-Larg*(int)y;
			if(l==dist) {p3=p;p33=p;}
			if(!stop[k+nbPoints/2+nbPoints/8])
			{
				if(((Image[p]+Image[p3]+Image[p33])>=160)&&(Image[p]<240)) 
				{
					stop[k+nbPoints/2+nbPoints/8]=TRUE;
					longueur[k+nbPoints/2+nbPoints/8]=(float) sqrt((dist-xcour)*(dist-xcour)*(1+pente*pente));//rayon-l;
				
					gradcirc[k+nbPoints/2+nbPoints/8]=(int)Image[po-l+3-Larg*(int)(pente*(l-3))]+(int)Image[po-l+2-Larg*(int)(pente*(l-2))]+(int)Image[po-l+4-Larg*(int)(pente*(l-4))];
//					Image[p]=255;
				}
				p33=p3;
				p3=p;
			}
			else if(Image[p]<=40||Image[p]>240)
			{
				longueur[k+nbPoints/2+nbPoints/8]=(float) sqrt((dist-x)*(dist-x)*(1+pente*pente));//rayon-l;
				
	//			gradcirc[k+nbPoints/2+nbPoints/8]=(int)Image[po-l+3-Larg*(int)(pente*(l-3))];
			}
			
//			
			p=po+(int)y-Larg*(int)x;
			if(l==dist) {p4=p;p44=p;}
			if(!stop[k+(nbPoints/4)*3+nbPoints/8])
			{
				if(((Image[p]+Image[p4]+Image[p44])>=160)&&(Image[p]<240)) 
				{
					stop[k+(nbPoints/4)*3+nbPoints/8]=TRUE;
					longueur[k+(nbPoints/4)*3+nbPoints/8]=(float) sqrt((dist-xcour)*(dist-xcour)*(1+pente*pente));//rayon-l;
				
					gradcirc[k+(nbPoints/4)*3+nbPoints/8]=(int)Image[po+(int)(pente*(l-3))-Larg*(l-3)]+(int)Image[po+(int)(pente*(l-2))-Larg*(l-2)]+(int)Image[po+(int)(pente*(l-4))-Larg*(l-4)];
//					Image[p]=255;
				}
				p44=p4;
				p4=p;
			}
			else if(Image[p]<=40||Image[p]>240)
			{
				longueur[k+(nbPoints/4)*3+nbPoints/8]=(float) sqrt((dist-x)*(dist-x)*(1+pente*pente));//rayon-l;
		
	//			gradcirc[k+(nbPoints/4)*3+nbPoints/8]=(int)Image[po+(int)(pente*(l-3))-Larg*(l-3)];
			}
			
//			
			xcour=x;
		}//pour l
	}//pour k
	BOOL sortie=TRUE;
	
	moyloc=longueur[0];
	for(i=1;i<=loc;i++) 
	{
		moyloc+=longueur[i];
		moyloc+=longueur[nbPoints-i];
	}

float /*long*/ max=0;
	for(i=0;i<nbPoints/4;i++)
	{
		tabmoy[i]=(float)(longueur[i]-moyloc/(2*loc+1));
		if(max<longueur[i]) max=longueur[i];
		moyloc+=longueur[(i+loc+1)%nbPoints];
		moyloc-=longueur[(nbPoints-loc+i)%nbPoints];
	}
	fin_i=Xg+rayon-(int)max-1;

	max=0;
	for(i=nbPoints/4;i<nbPoints/2;i++)
	{
		tabmoy[i]=(float)(longueur[i]-moyloc/(2*loc+1));
		if(max<longueur[i]) max=longueur[i];
		moyloc+=longueur[(i+loc+1)%nbPoints];
		moyloc-=longueur[(nbPoints-loc+i)%nbPoints];
	}
	fin_j=Yg+rayon-(int)max;

	max=0;
	for(i=nbPoints/2;i<3*nbPoints/4;i++)
	{
		tabmoy[i]=(float)(longueur[i]-moyloc/(2*loc+1));
		if(max<longueur[i]) max=longueur[i];
		moyloc+=longueur[(i+loc+1)%nbPoints];
		moyloc-=longueur[(nbPoints-loc+i)%nbPoints];
	}
	debut_i=Xg-rayon+(int)max;

	max=0;
	for(i=3*nbPoints/4;i<nbPoints;i++)
	{
		tabmoy[i]=(float)(longueur[i]-moyloc/(2*loc+1));
		if(max<longueur[i]) max=longueur[i];
		moyloc+=longueur[(i+loc+1)%nbPoints];
		moyloc-=longueur[(nbPoints-loc+i)%nbPoints];
	}
	debut_j=Yg-rayon+(int)max+1;
Masque=TRUE;
	debut_i+=4;
	debut_j+=4;
	fin_i-=4;
	fin_j-=4;
AffectRegionMasque();
max=0;
float tampon;
for(i=0;i<nbPoints;i++) 
{
	tampon=(tabmoy[i]+tabmoy[(i+1)%nbPoints]+tabmoy[(nbPoints+i-1)%nbPoints])/3;
	if(max<tampon) max=tampon;
}
float max3=0;
float DerivLong[200];
for(i=0;i<nbPoints;i++)
{
	//beta=20 alpha=0.44 -12 à +12
	DerivLong[i]=	(float)-(longueur[(nbPoints-12+i)%nbPoints]+2*longueur[(nbPoints-11+i)%nbPoints]+
							2*longueur[(nbPoints-10+i)%nbPoints]+3*longueur[(nbPoints-9+i)%nbPoints]+
							5*longueur[(nbPoints-8+i)%nbPoints]+6*longueur[(nbPoints-7+i)%nbPoints]+
							9*longueur[(nbPoints-6+i)%nbPoints]+11*longueur[(nbPoints-5+i)%nbPoints]+
							14*longueur[(nbPoints-4+i)%nbPoints]+16*longueur[(nbPoints-3+i)%nbPoints]+
							17*longueur[(nbPoints-2+i)%nbPoints]+13*longueur[(nbPoints-1+i)%nbPoints])
					+(float)(longueur[(12+i)%nbPoints]+2*longueur[(11+i)%nbPoints]+
							2*longueur[(10+i)%nbPoints]+3*longueur[(9+i)%nbPoints]+
							5*longueur[(8+i)%nbPoints]+6*longueur[(7+i)%nbPoints]+
							9*longueur[(6+i)%nbPoints]+11*longueur[(5+i)%nbPoints]+14*longueur[(4+i)%nbPoints]+
							16*longueur[(3+i)%nbPoints]+17*longueur[(2+i)%nbPoints]+13*longueur[(1+i)%nbPoints]);
	DerivLong[i]/=99;
	if(fabs(DerivLong[i])>fabs(max3)) max3=DerivLong[i];
}

//////////////////////////////////////////////
//////////////////////////////////////////////
/*
CFile Fich;
	long large,haute;
	DWORD DWComp;
	WORD BitCount;
	BITMAPFILEHEADER Tete;
	CString NomFichier2;
	unsigned char *pBuffer;

	CFileDialog * FileData = new CFileDialog
	(TRUE,"*.bmp",NULL,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT,"Images Bitmap(*.bmp)|*.bmp||",this);

	if (FileData->DoModal()==IDOK)
	{
		if (FileData->GetFileExt() == "bmp")
		{
			if (Fich.Open(FileData->GetPathName(), CFile::modeRead))
			{
				Fich.Read(&Tete,sizeof(BITMAPFILEHEADER));
				Fich.Seek(sizeof(DWORD),CFile::current);
				Fich.Read(&large,sizeof(LONG)); 
				Fich.Read(&haute,sizeof(LONG));
				Fich.Seek(sizeof(WORD),CFile::current);
				Fich.Read(&BitCount,sizeof(WORD));
				Fich.Read(&DWComp,sizeof(DWORD));
				if ((DWComp != BI_RGB)||(BitCount!=8 && BitCount!=24))
				{
					MessageBox("Impossible de lire les fichiers BMP compressés","Erreur",MB_OK|MB_ICONSTOP);
					Fich.Close();
					return 10.0;
				}
  	 			else if ((large!=Larg)||(haute!=Haut))
					MessageBox("Taille incompatible!","Erreur",MB_OK|MB_ICONSTOP);
				else
				{
					Fich.Seek(Tete.bfOffBits,CFile::begin);
					pBuffer = new unsigned char[Larg*Haut];
					Fich.Read(pBuffer,Larg*Haut);
				}
				Fich.Close();
			}
		}
		else MessageBox("Fichiers BMP uniquement");
	}
	else return 10.0;
	delete FileData;*/
//////////////////////////////////////

//////////////////////////////
//gradient circulaire
	int largde=10;
	float deriv=0.0;
	float derivtab[200];
	float derivtab2[200];
	float moycirc=(float)gradcirc[0];
	float moycirc2=(float) gradcirc[0];
	int *res=new int[nbPoints];
	for(i=1;i<11;i++) {moycirc+=gradcirc[i];moycirc+=gradcirc[nbPoints-i];}
	for(i=1;i<4;i++) {moycirc2+=gradcirc[i];moycirc2+=gradcirc[nbPoints-i];}

	for(i=1;i<=largde;i++) {deriv+=gradcirc[i];deriv-=gradcirc[nbPoints-i];}

	float max4=0.0;
	for(i=0;i<nbPoints;i++)
	{
		derivtab[i]=(float)(deriv/20.0);
		res[i]=abs(gradcirc[i]-(int)(moycirc/21.0));
		moycirc+=gradcirc[(i+10+1)%nbPoints];
		moycirc-=gradcirc[(nbPoints-10+i)%nbPoints];
		moycirc2+=gradcirc[(i+3+1)%nbPoints];
		moycirc2-=gradcirc[(nbPoints-3+i)%nbPoints];

		deriv+=gradcirc[(i+largde+1)%nbPoints]+gradcirc[(nbPoints-largde+i)%nbPoints];
		deriv-=gradcirc[i]+gradcirc[(i+1)%nbPoints];
//beta =20 alpha=0.56 -10 à +10
		derivtab2[i]=	(float)-(gradcirc[(nbPoints-10+i)%nbPoints]+gradcirc[(nbPoints-9+i)%nbPoints]+
								2*gradcirc[(nbPoints-8+i)%nbPoints]+3*gradcirc[(nbPoints-7+i)%nbPoints]+
								4*gradcirc[(nbPoints-6+i)%nbPoints]+6*gradcirc[(nbPoints-5+i)%nbPoints]+
								9*gradcirc[(nbPoints-4+i)%nbPoints]+11*gradcirc[(nbPoints-3+i)%nbPoints]+
								13*gradcirc[(nbPoints-2+i)%nbPoints]+11*gradcirc[(nbPoints-1+i)%nbPoints])+
						(gradcirc[(10+i)%nbPoints]+gradcirc[(9+i)%nbPoints]+2*gradcirc[(8+i)%nbPoints]+
						3*gradcirc[(7+i)%nbPoints]+4*gradcirc[(6+i)%nbPoints]+6*gradcirc[(5+i)%nbPoints]+
						9*gradcirc[(4+i)%nbPoints]+11*gradcirc[(3+i)%nbPoints]+13*gradcirc[(2+i)%nbPoints]+
						11*gradcirc[(1+i)%nbPoints]);

		derivtab2[i]/=48;
		if(fabs(derivtab2[i])>fabs(max4)) max4=derivtab2[i];
	}

	int prem=res[0],deux=res[1];
	float premoy=tabmoy[0],deuxmoy=tabmoy[1];
	max2=0;
	for(i=0;i<nbPoints-2;i++) {res[i]=(res[i]+res[i+1]+res[i+2])/3;tabmoy[i]=(tabmoy[i]+tabmoy[i+1]+tabmoy[i+2])/3;if(res[i]>max2) max2=res[i];}
	res[i]=(res[i]+res[i+1]+prem);
	tabmoy[i]=(tabmoy[i]+tabmoy[i+1]+premoy);
	if(res[i]>max2) max2=res[i];
	res[++i]=(res[i]+prem+deux);
	tabmoy[i]=(tabmoy[i]+premoy+deuxmoy)/3;
	if(res[i]>max2) max2=res[i];

////////////////////////////;
	BOOL METAL=FALSE;
	CPoint h;
	int j;
	for(i=Xg-rayon;i<Xg+rayon;i++)
		for(j=Yg-rayon;j<Yg+rayon;j++)
		{
			h=Inv_Get_ij(i,j);
			p=i+Larg*j;
			if(!RegionEllipse.PtInRegion(h)) 
			{
				if(Image[p]<139)Image[p]=147;
				if(Image[p]>=240) METAL=TRUE;
			}//pBuffer[p];//

		}
	Invalidate();
	if(METAL) max=10;
	CString ch;

	CGraphWnd *graphique0=(CGraphWnd*) new CGraphWnd(longueur,nbPoints,"Nombre","Niveau");
	graphique0->Create("longueur de rayons "+Titre,GetMDIFrame());
//	ch.Format("2x%d",largde);
//	CGraphWnd *graphique5=(CGraphWnd*) new CGraphWnd(derivtab,nbPoints,"Nombre","Niveau");
//	graphique5->Create("derivée luminance sur " +ch+" rayons de "+Titre,GetMDIFrame());

	ch.Format(" circul max=%2.2f\tmax=%2.2f\n lumi max=%d\tmax=%2.2f\n",max,max3,max2,max4);
	if(METAL)ch+="METAL";
	CGraphWnd *graphique1=(CGraphWnd*) new CGraphWnd(tabmoy,nbPoints,"Nombre","Niveau");
	graphique1->Create("Courbe de circularité sur "+Titre,GetMDIFrame());

	CGraphWnd *graphique10=(CGraphWnd*) new CGraphWnd(DerivLong,nbPoints,"Nombre","Niveau");
	graphique10->Create("Courbe de circularité deriche "+Titre,GetMDIFrame());

	CGraphWnd *graphique6=(CGraphWnd*) new CGraphWnd(derivtab2,nbPoints,"Nombre","Niveau");
	graphique6->Create("derivée deriche luminance sur 10 rayons de "+Titre,GetMDIFrame());

	//////////gradient circulaire
	CGraphWnd *graphique2=(CGraphWnd*) new CGraphWnd(res,nbPoints,"Nombre","Niveau");
	graphique2->Create("gradient circulaire "+Titre,GetMDIFrame());
	if(max>=2.6) MessageBox(ch,"PB de bords");
	else MessageBox(ch,"ok");
//	delete pBuffer;
	delete longueur;
	delete stop;
	
	delete gradcirc;
	delete res;

	free(Ysom);
	free(Xsom);
	return max;
}

void CTraitImageWnd::OnTraitementDessinmesures() 
{
	long*Ysom,*Xsom,Nbp=0,p=0;
	long Xg=0,Yg=0,po;
	int rayon,nbPoints;
	int/*float*/ x,y; 
	double Pi=acos(-1.0);
	
	rayon=140;
	nbPoints=200;
	
	Memorise();
	Ysom=(long*) malloc( Haut*sizeof(long));
	Xsom=(long*) malloc( Larg*sizeof(long));;

	for(long j=0;j<Haut;j++)
		for(long i=0;i<Larg;i++)
		{
			if(Image[p]>5) 
			{
				Xsom[i]++;
				Ysom[j]++;
				Nbp++;
			}
			p++;
		}
	long i;
	for( i=1;i<Larg;i++) Xg+=i*Xsom[i];

	if(Nbp==0) return ;
	for(i=1;i<Haut;i++) Yg+=i*Ysom[i];
	Xg/=Nbp;
	Yg/=Nbp;
	
	po =Xg+Larg*Yg;
	
	float pente;
	int dist;
	
	for(int k=-(nbPoints)/8;k<nbPoints/8;k++)
	{
		 
		///////////////////////////////////
		pente=(float)tan((2*Pi*k)/nbPoints);
		dist=(int)(rayon*cos((2*Pi*k)/nbPoints));
		///////////////////////////////////////////
		if(k==-(nbPoints)/8) dist++;
		if(k==0) dist--;
		for(int l = dist/*rayon*/;l>(int)((rayon-40)*cos((2*Pi*k)/nbPoints))/*105*/;l--)
		{
			//x=(float)(l*cos((2*Pi*k)/nbPoints));
			//y=(float)(l*sin((2*Pi*k)/nbPoints));
			x=l;
			y=(int)(pente*l);

			p=po+(int)x+Larg*(int)y;
			Image[p]=255;

			p=po-(int)y+Larg*(int)x;
			Image[p]=255;
			
			p=po-(int)x-Larg*(int)y;
			Image[p]=255;

			p=po+(int)y-Larg*(int)x;
			Image[p]=255;

		}//pour l
	}//pour k
Invalidate();	
}



 //Affichage coupe polytopes
void CTraitImageWnd::OnDcisionTest() 
{
	CPolyMc Polymc;
	long kpol,clas,k;
	int i,l,j,nbclas;
	unsigned char *classecoul;

	ArreteLeContinu();
	
	Polymc.LirePol(NomFichPMC);	
	Polymc.AffPol();
	int max=255;
	CTraitImageWnd *coupedec=new CTraitImageWnd(256,256,NB,TRUE);
	coupedec->Create("polytopes :projection en couches", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,NULL);
	coupedec->AlloueSequence(&max);

	nbclas=Polymc.NbClass;

	classecoul=new unsigned char[nbclas];

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
	}
	classecoul[0]=20;
	SetCursor(LoadCursor(NULL,IDC_WAIT));
	
	for( clas=0;clas<Polymc.NbClass;clas++)
	{
		for(kpol=0;kpol<Polymc.Npol[clas];kpol++)
  		{ 
			k=kpol*Polymc.NPAR;
			for(l=Polymc.poli[clas][k];l<=Polymc.pols[clas][k];l++)
				for(i=Polymc.poli[clas][k+1];i<=Polymc.pols[clas][k+1];i++)
					for(j=Polymc.poli[clas][k+2];(j<=Polymc.pols[clas][k+2])&&(j<255);j++)
					{
						if((coupedec->Bpile[j][l+256*i]!=classecoul[clas])&&(coupedec->Bpile[j][l+256*i]!=0))	coupedec->Bpile[j][l+256*i]=255;
						else coupedec->Bpile[j][l+256*i]=classecoul[clas];
					}
		}
	}
	Invalidate(FALSE);
	Polymc.Fin();
	SetCursor(LoadCursor(0,IDC_ARROW));
	delete classecoul;
}

void CTraitImageWnd::OnUpdateDcisionTest(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision);	
}



void CTraitImageWnd::OnDcisionRemplissagedelespacedesparam() 
{
	ArreteLeContinu();
	SetCursor(LoadCursor(NULL,IDC_WAIT));
	if(MessageBox("Voulez vous  détruire l'espace de décsision charger!","Attention",MB_YESNO)==IDNO) return;
	long p;
	int l,c,i;
	CPoint h;
	CRgn petiteellipse;
	
	byte* copEspClas;
	copEspClas=((CMainFrame*)GetMDIFrame())->EspClas;

	for(unsigned long k=0;k<16777216;k++)
	{
		copEspClas[k]=0;
	}
	int a;
	MasqEucl(a); 
	for (i=0;i<NbParam;i++) 
	{
		CalculParam(i);
	}
	
	AffectRegionMasque();
	int deci=(int) (((float)(fin_i-debut_i))/12.38);
	int deciprim=(int) (((float)(fin_i-debut_i))/13);
	int decj=(int) (((float)(fin_j-debut_j))/13);
	int decjprim=(int) (((float)(fin_j-debut_j))/12.38);
	
	petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+deci,debut_j+decj),Inv_Get_ij(fin_i-deciprim,fin_j-decjprim)));

	SetCursor(LoadCursor(NULL,IDC_WAIT));
	
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			h=Inv_Get_ij(c,l);
			if(RegionEllipse.PtInRegion(h))
			{	
				Calc(p);
				if(petiteellipse.PtInRegion(h))copEspClas[Buffer[0]+256*Buffer[1]+65536*Buffer[2]]=4;
				else copEspClas[Buffer[0]+256*Buffer[1]+65536*Buffer[2]]=3;
			
			}//if region
			p++;
		} // c
	} // l

	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::OnUpdateDcisionRemplissagedelespacedesparam(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision);	
}

void CTraitImageWnd::OnDcisionRemplissagedelespdesparamseq() 
{
	ArreteLeContinu();

	long p;
	int l,c,i,a;
	CPoint h;
	CRgn petiteellipse;
	int deci,deciprim,decj,decjprim;

	byte* copEspClas;
	copEspClas=((CMainFrame*)GetMDIFrame())->EspClas;



	for(int r=0;r<NbMaxImageSeq;r++)
	{
		AfficheSequence(r,r+1);
		
		MasqEucl(a);
		
		for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}
		
		deci=(int) (((float)(fin_i-debut_i))/12.38);// /9.5);
		decj=(int) (((float)(fin_j-debut_j))/13);// /9.92);
		deciprim=(int) (((float)(fin_i-debut_i))/13);// /9.92);
		decjprim=(int) (((float)(fin_j-debut_j))/12.38);// /9.5);
		petiteellipse.DeleteObject();
		petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+deci,debut_j+decj),Inv_Get_ij(fin_i-deciprim,fin_j-decjprim)));

		SetCursor(LoadCursor(NULL,IDC_WAIT));
		
		for (l=debut_j;l<=fin_j;l++)
		{
			p=Larg*l+debut_i;
		
			for (c=debut_i;c<=fin_i;c++)
			{
				h=Inv_Get_ij(c,l);
				if(RegionEllipse.PtInRegion(h))
				{	
					Calc(p);
					if(petiteellipse.PtInRegion(h))copEspClas[Buffer[0]+256*Buffer[1]+65536*Buffer[2]]=4;
					else copEspClas[Buffer[0]+256*Buffer[1]+65536*Buffer[2]]=3;
				
				}//if region
				p++;			
			} // c
		} // l
		OnMasqueAnnuler();

	}
	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::OnUpdateDcisionRemplissagedelespdesparamseq(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision&&SEQ_EXISTE);
}
/*
void CTraitImageWnd::OnDcisionTestprob() 
{
	SetCursor(LoadCursor(NULL,IDC_WAIT));
	long p;
	int l,c,i,nbclas=5;
	CPoint h;
	CRgn petiteellipse;
	double *zone;
	unsigned long /*zone,temp,max,tampon;
	unsigned char *classecoul;
	unsigned __int16 /*temp,/*max,nbpix[4]={0,0,0,0};
	float contamin;
	byte decide;
	ArreteLeContinu();
	unsigned __int32* copEspProb1,* copEspProb2,* copEspProb3,* copEspProb4;
	copEspProb1=((CMainFrame*)GetMDIFrame())->EspProb1;
	copEspProb2=((CMainFrame*)GetMDIFrame())->EspProb2;
	copEspProb3=((CMainFrame*)GetMDIFrame())->EspProb3;
	copEspProb4=((CMainFrame*)GetMDIFrame())->EspProb4;
	if(!AffectNonClas) nbclas++;
	classecoul=new unsigned char[nbclas];
	zone=new /*unsigned long double[nbclas];//////////////////

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
		zone[i]=0;
	}
	 contamin=MasqEucl(); 
	 for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}
	
	AffectRegionMasque();
	int deci=(int) (((float)(fin_i-debut_i))/12.38);
	int deciprim=(int) (((float)(fin_i-debut_i))/13);
	int decj=(int) (((float)(fin_j-debut_j))/13);
	int decjprim=(int) (((float)(fin_j-debut_j))/12.38);
	
	petiteellipse.CreateEllipticRgnIndirect(CRect(Inv_Get_ij(debut_i+deci,debut_j+decj),Inv_Get_ij(fin_i-deciprim,fin_j-decjprim)));

	SetCursor(LoadCursor(NULL,IDC_WAIT));
	
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			h=Inv_Get_ij(c,l);
			if(RegionEllipse.PtInRegion(h))
			{	
				Calc(p);
				decide =0;
				tampon=Buffer[0]+256*Buffer[1]+65536*Buffer[2];
				if(petiteellipse.PtInRegion(h))
				{
					temp=copEspProb1[tampon];
					max=temp;
					if(temp!=0){zone[1]+=temp;nbpix[0]++;decide=1;}
					temp=copEspProb2[tampon];
					if(temp!=0){zone[2]+=temp;nbpix[1]++;if(temp>max){decide=2;max=temp;}}
					temp=copEspProb3[tampon];
					if(temp!=0){zone[3]+=temp;nbpix[2]++;if(temp>max){decide=3;max=temp;}}
					temp=copEspProb4[tampon];
					if(temp!=0){zone[4]+=temp;nbpix[3]++;if(temp>max){decide=4;max=temp;}}
				}
				else
				{
					temp=copEspProb3[tampon];
					max=temp;
					if(temp!=0){zone[3]+=temp;nbpix[2]++;decide=3;}
					temp=copEspProb4[tampon];
					if(temp!=0){zone[4]+=temp;nbpix[3]++;if(temp>max){decide=4;max=temp;}}
				}
				Image[p]=classecoul[decide];
			}//if region
			else Image[p]=0;

			p++;
		} // c
	} // l
	if(TypePalette!=PALCOUL) OnCouleur();
//	OnCalculellipse();
	Invalidate(FALSE);

	CString ch="",ch0;
	double prob;
	for(i=1;i<nbclas;i++)
	{
		prob=((double)zone[i])/3000000;
		ch0.Format("Classe %d= %5.2f proba nb=%d\n\n",i,prob,nbpix[i-1]);
		ch+=ch0;
	}
	if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(5,!AffectNonClas);Legende->Invalidate(TRUE);}
	else
	{
		Legende->Init(5,!AffectNonClas);
		Legende->ShowWindow(SW_SHOW);
		Legende->Invalidate(TRUE);
	}
/*	if(zone[1]>37) 
		ch+="MANQUE\n";//80
	if(zone[2]>=50) 
		ch+="LISSE\n";
	if(zone[3]>12) 
		ch+="BOURGEON\n";
	if (zone[4]>=20)
		ch+="ARRACHEMENT\n";	
	if((zone[1]<=37)&&(zone[2]<50)&&(zone[3]<=12)&&(zone[4]<20)) 
		ch+="BONNE";
	MessageBox(ch);

	SetCursor(LoadCursor(0,IDC_ARROW));
	delete classecoul;
	delete zone;	
}

void CTraitImageWnd::OnUpdateDcisionTestprob(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision);
}
*/

void CTraitImageWnd::OnTraitementRductiondeniveaux() 
{
	Memorise();
	int dec;
	CDialVal Dlg;	
	if(Dlg.DoModal()==IDOK)
	{
		dec=atoi(Dlg.m_val);
		for(long i=0;i<Larg*Haut;i++) {if(Image[i]+dec>=0) Image[i]+=dec;else Image[i]=0;}
		Invalidate();
	}

}

void CTraitImageWnd::OnDcisionParzenImagecomplte() 
{
	SetCursor(LoadCursor(NULL,IDC_WAIT));
	long p;
	int l,c,i,nbclas=6;
	CPoint h;
	CRgn petiteellipse;
	long *zone;
	unsigned char *classecoul;
	byte decide;
	ArreteLeContinu();
	byte* copEspClas;
	copEspClas=((CMainFrame*)GetMDIFrame())->EspClas;
	if(!AffectNonClas) nbclas++;
	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];//////////////////

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
		zone[i]=0;
	}
	 for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}
	
	SetCursor(LoadCursor(NULL,IDC_WAIT));
	
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			{	
				Calc(p);
				decide=copEspClas[Buffer[0]+256*Buffer[1]+65536*Buffer[2]];
				if(decide==0)
				{
					if(AffectNonClas) decide=1;
					else decide=nbclas;
				}
				decide--;
				Image[p]=classecoul[decide];
				zone[decide]++;
			}//if region
			

			p++;
		} // c
	} // l
	if(TypePalette!=PALCOUL) OnCouleur();
//	OnCalculellipse();
	Invalidate(FALSE);

	CString ch="",ch0;
	for(i=0;i<nbclas;i++)
	{
		ch0.Format("Classe %d= %d pixels\n\n",i,zone[i]);
		ch+=ch0;
	}
	if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(nbclas,!AffectNonClas);Legende->Invalidate(TRUE);}
	else
	{
		Legende->Init(nbclas,!AffectNonClas);
		Legende->ShowWindow(SW_SHOW);
		Legende->Invalidate(TRUE);
	}
	if(zone[1]>37) 
		ch+="MANQUE\n";//80
	if(zone[2]>=50) 
		ch+="LISSE\n";
	if(zone[3]>12) 
		ch+="BOURGEON\n";
	if (zone[4]>=20)
		ch+="ARRACHEMENT\n";	
	if((zone[1]<=37)&&(zone[2]<50)&&(zone[3]<=12)&&(zone[4]<20)) 
		ch+="BONNE";
	MessageBox(ch);

	SetCursor(LoadCursor(0,IDC_ARROW));
	delete classecoul;
	delete zone;		
}

void CTraitImageWnd::OnUpdateDcisionParzenImagecomplte(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(ValideDecision);
}
void CTraitImageWnd::OnTraitementCodageextrema() 
{
	int c,l;
	unsigned char cod,cour;

	CTraitImageWnd * pImagecod =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,NB);   

	for(long i=Larg;i<Larg*Haut-Larg;i++)
	{ 
		cod=0;
		cour=Image[i];
		for(c=-1;c<2;c++)
			for(l=-1;l<2;l++)
			{
				if((c!=0 || l!=0) && cour>=Image[i+c+l*Larg]) cod++;
			}
	 pImagecod->Image[i]=31*cod;
	}
	pImagecod ->Create(Titre+"codage ",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
}


void CTraitImageWnd::OnTraitementDensitdanszone() 
{
	if(debut_i==0&&debut_j==0&&fin_i==Larg-1&&fin_j==Haut-1) {MessageBox("Définir une zone !");return;}
	CString Tit="Niveau de gris";
	CDialVal dlg(&Tit);
	if(dlg.DoModal()!=IDOK) return;

	int co=(fin_i-debut_i+1)/2;
	int lo=(fin_j-debut_j+1)/2;
	int nbpix=(2*lo+1)*(2*co+1);
	int compt,l,c;
	unsigned char niv,* tampon=new unsigned char[Larg*Haut];
	long i,fin=Larg-2*co-1;

	niv=abs(atoi(dlg.m_val));

	for(i=0;i<Larg*Haut;i++) tampon[i]=0;
	for(long j=lo;j<Haut-lo;j++)
	{
		for(i=co;i<Larg-co;i++)
		{ 
			compt=0;
			for(c=-co;c<=co;c++)
				for(l=-lo;l<=lo;l++)
				{
					if(Image[i+c+(l+j)*Larg]==niv) compt++;
				}
			tampon[i+j*Larg]=(unsigned char) (255*((float) compt/nbpix));
		}
	}
	for(i=0;i<Larg*Haut;i++) Image[i]=tampon[i];
	Invalidate();
	delete tampon;	
}



void CTraitImageWnd::Rotation(float radians)
{



	// Compute the cosine and sine only once
	float cosine = (float)cos(radians);
	float sine = (float)sin(radians);

	// Compute dimensions of the resulting bitmap
	// First get the coordinates of the 3 corners other than origin
	int x1 = (int)(-Haut * sine);
	int y1 = (int)(Haut * cosine);
	int x2 = (int)(Larg * cosine - Haut * sine);
	int y2 = (int)(Haut * cosine + Larg * sine);
	int x3 = (int)(Larg * cosine);
	int y3 = (int)(Larg * sine);

	int minx = min(0,min(x1, min(x2,x3)));
	int miny = min(0,min(y1, min(y2,y3)));
	int maxx = max(x1, max(x2,x3));
	int maxy = max(y1, max(y2,y3));

	int w = maxx - minx;
	int h = maxy - miny;

    CTraitImageWnd *Rota=new CTraitImageWnd(w,h,BMP);
//	 CTraitImageWnd *Rota2=new CTraitImageWnd(w,h,BMP);

	w=Rota->Larg;
	// Now do the actual rotating - a pixel at a time
	// Computing the destination point for each source point
	// will leave a few pixels that do not get covered
	// So we use a reverse transform - e.i. compute the source point
	// for each destination point
	double r1,r2,r3,r4;
	for( int y = 0; y < h; y++ )
	{
		for( int x = 0; x < w; x++ )
		{
			
			double subx=((x+minx)*cosine + (y+miny)*sine);
			double suby=((y+miny)*cosine - (x+minx)*sine);
			int sourcex = (int)subx;
			subx-=sourcex;
			int sourcey = (int)suby;
			suby-=sourcey;

			r1=sqrt(subx*subx+suby*suby);
			r2=sqrt((1-subx)*(1-subx)+suby*suby);
			r3=sqrt(subx*subx+(1-suby)*(1-suby));
			r4=sqrt((1-subx)*(1-subx)+(1-suby)*(1-suby));
			if( sourcex >= 0 && sourcex < Larg && sourcey >= 0 
				&& sourcey < Haut )
			{
					BYTE pixel;//,pixel2 ;
					if(subx!=0 || suby!=0)
					{
					/*	pixel2 =(unsigned char)(((long)2*Image[sourcex+ Larg*sourcey]+
						(Image[sourcex+1+Larg*sourcey]-Image[sourcex+ Larg*sourcey])*(subx)+
						(Image[sourcex+Larg*(sourcey+1)]-Image[sourcex+ Larg*sourcey])*(suby))/2);

					*/	pixel =	(unsigned char)((double)(r2*r3*r4*Image[sourcex+Larg*sourcey]
												+r1*r3*r4*Image[sourcex+1+Larg*sourcey]
												+r1*r2*r4*Image[sourcex+Larg*(sourcey+1)]
												+r1*r2*r3*Image[sourcex+1+Larg*(sourcey+1)])/(r1*r2*r3+r1*r2*r4+r1*r3*r4+r2*r3*r4));
					}							
					else {pixel=Image[sourcex+Larg*sourcey];}//pixel2=pixel;}
					Rota->Image[x + w*y]= pixel;
				//	Rota2->Image[x+w*y]=pixel2;
			}
			else 
			{
				// Draw the background color. The background color
				// has already been drawn for 8 bits per pixel and less
				Rota->Image[x + w*y]= 0;
			//	Rota2->Image[x + w*y]= 0;
			}
		}
	}
	Rota->Create("Rotation  "+Titre, WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());
//	Rota2->Create("Rotation simple "+Titre, WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());
}

void CTraitImageWnd::OnTraitementRotation() 
{
	CString Tit="angle de rotation";
	CDialVal fenedial(&Tit);
	if(fenedial.DoModal()!=IDOK) return;
  	 	
  
   	float nn=(float)atof(fenedial.m_val); 
	Rotation((float)(nn*acos(-1.0)/180));
}

void CTraitImageWnd::OnTraitementOprateursunairesSeuillagedewen() 
{
  //	Seuillage a n niveaux (n<10) par Wen d'une image float 
	long j,i,ll;      
   	int k,n,nn,nm,lar,ni;
    float histo[255];
    byte min,max;
	byte p;
    double kk,zz;
    double m[20];
    double mat[110];
    double c[10],z[10],q[10];
    double a[10],b[10],w[10];
    float seuil[10];
    double rep,qt;
    double s,ds;
	                               
   	// traitement sur l'image de PImageFene
   	// Calcul de l'histogramme
	CString Tit="Nb de niveaux";
	CDialVal fenedial(&Tit);
	if(fenedial.DoModal()!=IDOK) return;
  	 	
  
   	nn=atoi(fenedial.m_val); // nn est le nombre de niveaux, compris entre 2 et 9
	if(nn>10||nn<2) {MessageBox("nb de niveaux entre 2 et 10!");return;}
	CTraitImageWnd *PImageFene = new CTraitImageWnd(fin_i-debut_i+1, fin_j-debut_j+1, NB);   

   	min=Image[0];
   	max=Image[0];
   	
   	for (j=debut_j;j<=fin_j;j++) 
   	{	                   
   		for (i=debut_i;i<=fin_i;i++) 
   		{	 
			 p=Image[j*Larg+i];
			 if(p>max)max=p;
			 if(p<min)min=p;
   		}
   	}         
   	for (k=0;k<256;k++)histo[k]=0;

   	ll=0;
   	for (j=debut_j;j<=fin_j;j++) 
   	{	                   
   		for (i=debut_i;i<=fin_i;i++) 
   		{	 
			 p=Image[j*Larg+i];
			 histo[p]++;
			 ll++;
   		}
   	}      
   	for (k=0;k<256;k++)histo[k]=histo[k]/ll;
   	
   	
	// Calcul des moments
	for (n=1;n<2*nn;n++)m[n]=0;
	m[0]=1;
	for (k=0;k<256;k++)
	{   
		kk=1;
		for (n=1;n<2*nn;n++)
		{   
			kk=kk*k;
			m[n]=m[n]+histo[k]*kk;
		}
	}                    
	
	// calcul des coefficients C[]
	lar=nn+1;	
	for (j=0;j<nn;j++)
	{
		for (i=0;i<nn;i++) mat[j*lar+i]=m[2*nn-2-i-j];
	 	mat[j*lar+nn]=-m[2*nn-j-1];
	}	                   
	for (k=0;k<nn;k++)
	{
	 	for (i=nn;i>k;i--)mat[k*lar+i]=mat[k*lar+i]/mat[k*lar+k];
	 	mat[k*lar+k]=1;
	 	for (j=0;j<k;j++)
	 	{
	 		for (i=nn;i>k;i--)mat[j*lar+i]=mat[j*lar+i]-mat[j*lar+k]*mat[k*lar+i];
	 		mat[j*lar+k]=0;
	 	}
	 	for (j=k+1;j<nn;j++)
	 	{
	 		for (i=nn;i>k;i--)mat[j*lar+i]=mat[j*lar+i]-mat[j*lar+k]*mat[k*lar+i];
	 		mat[j*lar+k]=0;
	 	}
	}	               
	for (i=0;i<nn;i++)c[i]=mat[(nn-1-i)*lar+nn];		

	// Calcul des niveaux equivalents
	for (n=1;n<=nn;n++)
	{
		a[n]=c[nn-n];
		z[n]=0;
	}           
	a[0]=1;
	for (nm=nn;nm>=1;nm--)
	{
		s=0;ni=0;
		do
		{   
			ni++;
			b[0]=1;w[0]=0;
			for (n=1;n<=nm;n++)
			{
				b[n]=a[n]+s*b[n-1];
				w[n]=b[n-1]+s*w[n-1];
			}
			ds=-b[nm]/w[nm];
			s=s+ds;	
		}
		while((fabs(b[nm])>1E-9)&&(ni<30));
		z[nn-nm+1]=s;
		for (n=1;n<nm;n++)a[n]=b[n];
	}
	
	// Calcul des effectifs
	lar=nn+1;	
	for (i=0;i<nn;i++)
	{
		zz=1;
		for(j=nn-1;j>=0;j--)
		{
			mat[j*lar+i]=zz;
			zz=zz*z[nn-i];
		}
	}
	for (j=0;j<nn;j++)mat[j*lar+nn]=m[nn-j-1];
	for (k=0;k<nn;k++)
	{
	 	for (i=nn;i>k;i--)mat[k*lar+i]=mat[k*lar+i]/mat[k*lar+k];
	 	mat[k*lar+k]=1;
	 	for (j=0;j<k;j++)
	 	{
	 		for (i=nn;i>k;i--)mat[j*lar+i]=mat[j*lar+i]-mat[j*lar+k]*mat[k*lar+i];
	 		mat[j*lar+k]=0;
	 	}
	 	for (j=k+1;j<nn;j++)
	 	{
	 		for (i=nn;i>k;i--)mat[j*lar+i]=mat[j*lar+i]-mat[j*lar+k]*mat[k*lar+i];
	 		mat[j*lar+k]=0;
	 	}
	}	               
	for (i=0;i<nn;i++)q[i+1]=mat[(nn-1-i)*lar+nn];		
	
	// Calcul des seuils
	rep=0;k=0;    
	seuil[0]=min;qt=0;
	for (n=1;n<nn;n++)
	{                 
		qt=qt+q[n];
		while (rep<qt)
		{
			rep=rep+histo[k];k++;
		}                        
		seuil[n]=(float)((k-1))+min;
	}
	seuil[nn]=max;
   	
   	// seuillage
   	ll=0;
   	for (j=debut_j;j<=fin_j;j++) 
   	{	                   
   		for (i=debut_i;i<PImageFene->Larg+debut_i;i++) 
   		{	 
			p=Image[j*Larg+i];
			n=1;
			while(p>seuil[n])n++;
			PImageFene->Image[ll++]=int(z[n]);
		}
   	}   
   	// fin du traitement
   	PImageFene->Create(Titre+"  Image seuillée multiniveau", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());	
}	

void CTraitImageWnd::OnTraitementMorphologieAmincissement() 
{
	Memorise();
	unsigned char * BufMorphoAmin;
	DWORD TailleBuf = Larg*Haut;
	int  tmp;
	unsigned char mini,maxi;
	long p,i,j,l=0, compteur, boucle=0;
	int Tab[8]={-1-Larg,-1,-1+Larg,Larg,1+Larg,1,1-Larg,-Larg};
	// Coefficients sur pixels du voisinage 8-connexe
	//		2	3	4
	// tab	1		5
	//		0	7	6
	BufMorphoAmin = (unsigned char*)malloc(TailleBuf);            
	long r;
	for(r=0;r<(int)TailleBuf;r++)BufMorphoAmin[r]=Image[r];
	SetCursor(LoadCursor(0,IDC_WAIT));

	do	{
		compteur=0;
			for(l=0;l<8;l++)
			{
				for(j=1;j<Haut-1;j++)
				{
					p=1+j*Larg;
					for(i=1;i<Larg-1;i++)
					{
						maxi=max(Image[p+Tab[0]],max(Image[p+Tab[7]],Image[p+Tab[6]]));
						mini=min(Image[p+Tab[2]],min(Image[p+Tab[3]],Image[p+Tab[4]]));
						if(Image[p]>maxi &&Image[p]<=mini){BufMorphoAmin[p]=maxi;compteur++;}
						p++;
					}
				}
				for(r=0;r<(int)TailleBuf;r++)Image[r]=BufMorphoAmin[r];
				tmp = Tab[0];
				for (i=1;i<8;i++)
					Tab[i-1] = Tab[i];
				Tab[7] = tmp;
			}

	} while (compteur);

  	free(BufMorphoAmin);
	SetCursor(LoadCursor(0,IDC_ARROW));
	Invalidate(TRUE);
}

void CTraitImageWnd::OnTraitementMorphologieEpaississement() 
{
	Memorise();
	unsigned char * BufMorphoAmin;
	DWORD TailleBuf = Larg*Haut;
	int  tmp;
	unsigned char mini,maxi;
	long p,i,j,l=0, compteur, boucle=0;
	int Tab[8]={-1-Larg,-1,-1+Larg,Larg,1+Larg,1,1-Larg,-Larg};
	// Coefficients sur pixels du voisinage 8-connexe
	//		2	3	4
	// tab	1		5
	//		0	7	6
	BufMorphoAmin = (unsigned char*)malloc(TailleBuf);            
	long r;
	for(r=0;r<(int)TailleBuf;r++)BufMorphoAmin[r]=Image[r];
	SetCursor(LoadCursor(0,IDC_WAIT));

	do	{
			compteur=0;
			for(l=0;l<8;l++)
			{
				for(j=1;j<Haut-1;j++)
				{
					p=1+j*Larg;
					for(i=1;i<Larg-1;i++)
					{
						maxi=max(Image[p+Tab[0]],max(Image[p+Tab[7]],Image[p+Tab[6]]));
						mini=min(Image[p+Tab[2]],min(Image[p+Tab[3]],Image[p+Tab[4]]));
						if(Image[p]>=maxi &&Image[p]<mini){BufMorphoAmin[p]=mini;compteur++;}
						p++;
					}
				}
				for(r=0;r<(int)TailleBuf;r++)Image[r]=BufMorphoAmin[r];
				tmp = Tab[0];
				for (i=1;i<8;i++)
					Tab[i-1] = Tab[i];
				Tab[7] = tmp;
			}		
	} while (compteur);//&& !boucle);

  	free(BufMorphoAmin);
	SetCursor(LoadCursor(0,IDC_ARROW));
	Invalidate(TRUE);	
}

void CTraitImageWnd::OnTraitementMorphologieBarbulage() 
{
	CString Tit="nombre d'itération";
	CDialVal fenedial(&Tit);
	if(fenedial.DoModal()!=IDOK) return;
	Memorise();
	unsigned char * BufMorphoAmin;
	DWORD TailleBuf = Larg*Haut;
	int  tmp;
	unsigned char maxi;
	long p,i,j,l=0, compteur, boucle=0;
	int Tab[8]={-1-Larg,-1,-1+Larg,Larg,1+Larg,1,1-Larg,-Larg};
	// Coefficients sur pixels du voisinage 8-connexe
	//		2	3	4
	// tab	1		5
	//		0	7	6
	BufMorphoAmin = (unsigned char*)malloc(TailleBuf);            
	long r;
	for( r=0;r<(int)TailleBuf;r++)BufMorphoAmin[r]=Image[r];
	SetCursor(LoadCursor(0,IDC_WAIT));
	compteur=atoi(fenedial.m_val);

	do	{
			for(l=0;l<8;l++)
			{
				for(j=1;j<Haut-1;j++)
				{
					p=1+j*Larg;
					for(i=1;i<Larg-1;i++)
					{
						maxi=max(Image[p+Tab[0]],max(Image[p+Tab[1]],max(Image[p+Tab[7]],max(Image[p+Tab[6]],Image[p+Tab[5]]))));
						if(Image[p]>maxi){BufMorphoAmin[p]=maxi;}
						p++;
					}
				}
				for(r=0;r<(int)TailleBuf;r++)Image[r]=BufMorphoAmin[r];
				tmp = Tab[0];
				for (i=1;i<8;i++)
					Tab[i-1] = Tab[i];
				Tab[7] = tmp;
			}
			compteur--;
	} while (compteur);

  	free(BufMorphoAmin);
	SetCursor(LoadCursor(0,IDC_ARROW));
	Invalidate(TRUE);
}




void CTraitImageWnd::OnDcisionDecisionsurfloat() 
{
	// TODO: Add your command handler code here

}

void CTraitImageWnd::OnDcisionDecidesurfloatdebug() 
{
	// TODO: Add your command handler code here
			CPolyMcFloat Polymc;
	long p,kpol,clas,k;
	int l,c,par,i,nbclas;
	BYTE ext,decide,aa,nonclas;


		float aafloat,coefnorm,Tmin[256],Tmax[256];
	BYTE NPAR;
	long *zone;
	unsigned char *classecoul;

	CFile fi;
	if(fi.Open("d:\\svm\\reapp.lim",CFile::modeRead))
	{
		fi.Read(&NPAR,sizeof(NPAR));
		fi.Read(Tmin,NPAR*sizeof(float));
		fi.Read(Tmax,NPAR*sizeof(float));
		fi.Close();
	}
	else MessageBox("fichier lim absent");

	ArreteLeContinu();
	
	Polymc.LirePol("d:\\svm\\reapp.pmf");	
	Polymc.AffPol();

	nbclas=Polymc.NbClass;

	if(!AffectNonClas) 
	{
		nbclas++;
	}	

	classecoul=new unsigned char[nbclas];
	zone=new long[nbclas];

	for(i=0;i<nbclas;i++) 
	{
		classecoul[i]=(unsigned char) (i*((double)255/nbclas));
		zone[i]=0;
	}
	  for (i=0;i<NbParam;i++) 
		{
			CalculParam(i);
		}

	SetCursor(LoadCursor(NULL,IDC_WAIT));
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			Calc(p);
			decide=0;
			nonclas=1;

			for( clas=0;clas<Polymc.NbClass;clas++)
			{
			   for(kpol=0;kpol<Polymc.Npol[clas];kpol++)
  				{ ext=0;
  				 for(par=0;par<Polymc.NPAR;par++)
 				   { 
  					   	coefnorm=(Tmax[par]-Tmin[par])/(0.8);
						k=kpol*Polymc.NPAR+par;
 					   aa=Buffer[par];
					   if (aa>(BYTE)Tmax[par]) aa=(BYTE)Tmax[par];
					   if (aa<(BYTE)Tmin[par]) aa=(BYTE)Tmin[par];

					   aafloat=0.1+((float)(aa)-Tmin[par])/coefnorm;

 					   if (aafloat<Polymc.poli[clas][k]) {ext=1;break;}
 					   else   if (aafloat>Polymc.pols[clas][k]) {ext=1;break;}
   				   }//for par
				   if (!ext) {decide=classecoul[clas];zone[clas]++;nonclas=0;break;}
				}//for pol
			}//for clas
			if(nonclas)
			{
				if(AffectNonClas)
				{
					zone[ClassAffect]++;
					decide=classecoul[ClassAffect];
				}
				else 
				{
					zone[nbclas-1]++;
					decide=classecoul[nbclas-1];
				}
			}
			Image[p]=decide;
			p++;
		} // c
	} // l
	if(TypePalette!=PALCOUL) OnCouleur();
	Invalidate(FALSE);

	CString ch="",ch0;;
	for(i=0;i<Polymc.NbClass;i++)
	{
		ch0.Format("Classe %d= %d pixels\n\n",i,zone[i]);
		ch+=ch0;
	}
	if(!AffectNonClas)
	{
		ch0.Format("Non classés= %d pixels\n\n",zone[i]);
		ch+=ch0;
	}
	if(!IsWindow(Legende->GetSafeHwnd())){Legende->Create(Polymc.NbClass,!AffectNonClas);Legende->Invalidate(TRUE);}
	else
	{
		Legende->Init(Polymc.NbClass,!AffectNonClas);
		Legende->ShowWindow(SW_SHOW);
		Legende->Invalidate(TRUE);
	}
	MessageBox(ch);
	Polymc.Fin();
	delete classecoul;
	delete zone;
	SetCursor(LoadCursor(0,IDC_ARROW));
}

void CTraitImageWnd::OnDcisionDecidesvm() 
{
	// TODO: Add your command handler code here
	long p,kpol,clas,k;
	int l,c,par,i,nbclasinit,nbclas;
	BYTE ext,decide,aa,nonclas;
	CStdioFile FichierResult;
	CString Chaine;



//	WinExec("d:\\svm\\svmpredict.exe esfil_svm.dat esfil_model.model result.txt",SW_SHOW);
	FichierResult.Open("d:\\svm\\result.txt",CFile::modeRead);
	
	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{

			FichierResult.ReadString(Chaine);
			if (Chaine=="1") decide=0;
			else decide=255;

			Image[p]=decide;
			p++;
		} // c
	} // l

	FichierResult.Close();

	if(TypePalette!=PALCOUL) OnCouleur();
	Invalidate(FALSE);


}

void CTraitImageWnd::OnDcisionGenerefichierdat() // pour svm
{
	// TODO: Add your command handler code here
	long p,kpol,clas,k;
	int l,c,par,i,nbclasinit,nbclas;
	CStdioFile FichierSVM;
	float fa;
	CString Chaine,Temp="";
	ArreteLeContinu();
	
	FichierSVM.Open("essai_svm.dat",CFile::modeWrite|CFile::modeCreate);

	nbclasinit=2; // debug
 
	for (i=0;i<NbParam;i++) 
	{
		CalculParam(i);
	}

	for (l=debut_j;l<=fin_j;l++)
	{
		p=Larg*l+debut_i;
	
		for (c=debut_i;c<=fin_i;c++)
		{
			Temp="1 ";
			Calc(p); // rempli buffer 
 			 for(par=0;par<NbParam;par++)
 			   { 
 				   fa=(float)Buffer[par];
				   Chaine.Format("%d:%3.6f ",par+1,fa);
				   Temp+=Chaine;
   			   }//for par
			Temp+=" \n";
			FichierSVM.WriteString(Temp);
 	

			p++;
		} // c
	} // l
	FichierSVM.Close();
	if(TypePalette!=PALCOUL) OnCouleur();
	Invalidate(FALSE);	
}

void CTraitImageWnd::OnTraitementComparaisonimbin() 
{
	// TODO: Add your command handler code here
	
		int  num1=1,num2=35;
	CString s,NomFichier;
	CFile Fichier;
	long large,haute;
	DWORD DWComp;
	WORD BitCount;
	BITMAPFILEHEADER Tete;
	CString NomFichier2;
	unsigned char *pBuffer;


	long e0=0,e1=0;
	long noir=0,blanc=0;

	CFileDialog * FileData = new CFileDialog
	(TRUE,"*.bmp",NULL,OFN_HIDEREADONLY|OFN_OVERWRITEPROMPT,"Images Bitmap(*.bmp)|*.bmp||",this);

	if (FileData->DoModal()==IDOK)
	{
		if (FileData->GetFileExt() == "bmp")
		{
			if (Fichier.Open(FileData->GetPathName(), CFile::modeRead))
			{
				Fichier.Read(&Tete,sizeof(BITMAPFILEHEADER));
				Fichier.Seek(sizeof(DWORD),CFile::current);
				Fichier.Read(&large,sizeof(LONG)); 
				Fichier.Read(&haute,sizeof(LONG));
				Fichier.Seek(sizeof(WORD),CFile::current);
				Fichier.Read(&BitCount,sizeof(WORD));
				Fichier.Read(&DWComp,sizeof(DWORD));
				if ((DWComp != BI_RGB)||(BitCount!=8 && BitCount!=24))
				{
					MessageBox("Impossible de lire les fichiers BMP compressés","Erreur",MB_OK|MB_ICONSTOP);
					Fichier.Close();
					return;
				}
  	 			else if ((large!=Larg)||(haute!=Haut))
					MessageBox("Taille incompatible!","Erreur",MB_OK|MB_ICONSTOP);
				else
				{
					Fichier.Seek(Tete.bfOffBits,CFile::begin);
					pBuffer = (unsigned char *) malloc(Larg*Haut);
					Fichier.Read(pBuffer,Larg*Haut);

					for(long i=0;i<Larg*Haut;i++)
					{
						if (Image[i]==0 && pBuffer[i]!=0) e0++;
						if (Image[i]!=0 && pBuffer[i]==0) e1++;
						if (Image[i]==0) noir++; else blanc++;

					}

					Invalidate(FALSE);
					free(pBuffer);
				}
				Fichier.Close();
			}
		}
		else MessageBox("Différence d'images entre des fichiers BMP uniquement");
	}
	delete FileData;
	CString Ch;
	float f0,f1,ft;
	f0=float(e0)/float(noir);
	f0=f0*100;
	f1=float(e1)/float(blanc);
	f1=f1*100;
	ft=(float(e0)+float(e1))/(float(noir+blanc));
	ft=ft*100;
	Ch.Format(" f0= %3.5f    f1=%3.5f\n ft=%3.5f",f0,f1,ft);
	MessageBox(Ch);




	

}

void CTraitImageWnd::OnTraitementAjustement() 
{
	// TODO: Add your command handler code here
	unsigned char * temp;
	unsigned char * Sobel;
	unsigned char * Masque;
	long Taille=Larg*Haut,i,l,p,c1,p1,c2,p2,c,l1,l2;
	Sobel=new unsigned char[Taille];
	temp=new unsigned char[Taille];
	Masque=new unsigned char[Taille];

	ArreteLeContinu();
	Memorise();

	int LignesPasBlanche=0;
	AJUSTD=0;
	AJUSTG=0;

//	Erode(GetSImage(),1);
	OnMoyenneurxy111();
	OnMoyenneurxy111();

	for(i=0;i<Taille;i++) temp[i]=Image[i];
	for(i=Larg;i<(Larg-1)*Haut;i++)
	{
		Sobel[i]= (unsigned char)(labs(2*(long)temp[i-1]+(long)temp[i+Larg-1]+(long)temp[i-Larg-1]-2*(long)temp[i+1]-(long)temp[i+Larg+1]-(long)temp[i-Larg+1])+
  					 labs((long)temp[i+Larg-1]+2*(long)temp[i+Larg]+(long)temp[i+Larg+1]-(long)temp[i-Larg-1]-2*(long)temp[i-Larg]-(long)temp[i-Larg+1]));
	}

	for(i=0;i<Taille;i++)
	{
		if (Image[i]>83) Masque[i]= 255;
		else Masque [i]=0;
	}
	
	
	for(l=1;l<Haut;l++)
	{
		 p=Larg*l;
	//	 do
	//	 {
			do
			{
				p++;
			} while(Masque[p]==0 && (p%Larg)<(Larg-AJUST));
			if ((p%Larg)<(Larg-AJUST))
			{
				c1=p%Larg;
				p1=p;
				do
				{
					p++;
				} while(Masque[p]==255 && (p%Larg)<(Larg-AJUST));

				c2=(p%Larg);
				p2=p;
				ajustx(c1,c2,l,Sobel,Masque,FALSE);
				LignesPasBlanche++;
			} // si pas ligne blanche
	  // } while ((p%Larg)<(Larg-AJUST));
	 
	} // for l
	AJUSTD/=(LignesPasBlanche);
	AJUSTG/=(LignesPasBlanche);
	
	CString Stemp;
	Stemp.Format("AjustD=%d		AjustG=%d",AJUSTD,AJUSTG);
	MessageBox(Stemp);

	for(l=1;l<Haut;l++)
	{
	 p=Larg*l;
//	 do
//	 {
        do
        {
			p++;
        } while(Masque[p]==0 && (p%Larg)<(Larg-AJUST));
        if ((p%Larg)<(Larg-AJUST))
        {
			 c1=p%Larg;
			 p1=p;

		 do
			{
				p++;
			} while(Masque[p]==255 && (p%Larg)<(Larg-AJUST));
			c2=(p%Larg);
			p2=p;
			ajustx(c1,c2,l,Sobel,Masque,TRUE);
        } // si pas ligne blanche
  // } while ((p%Larg)<(Larg-AJUST));
}
	
	





for(c=1;c<Larg;c++)
{
   p=c;
        do
        {
        p+=Larg;
        } while(Masque[p]==0 && (p/Larg)<(Haut-AJUSTH));

        if ((p/Larg)<(Haut-AJUSTH))
        {
			l1=p/Larg;
			p1=p;

			do
			{
			p+=Larg;
			} while(Masque[p]==255 && (p/Larg)<(Haut-AJUSTH));
			l2=(p/Larg);
			p2=p;
        
			ajusty(l1,l2,c,Sobel,Masque);
		} // si pas colonne blanche

} // pour les colonnes





	for(i=0;i<Taille;i++) Image[i]=Masque[i];

	Invalidate(FALSE);	
	delete temp;
	delete Sobel;
	delete Masque;
}



void CTraitImageWnd::ajustx(int c1,int c2,int l,byte *point,byte *bram,BOOL Flag)
{
	byte gradmax=0;
	word pgrad,p,p1,p2;
	int c,ci;
		
	p2=Larg*l+c2;
	p1=Larg*l+c1;

	/* Gauche */
	pgrad=p1;	
	//        for (c=c1-AJUST;(c<=c1+AJUST && c<c2);c++)
	if (!Flag)
	{
		for (c=c1-AJUST;c<=c1;c++)
		{
			p=Larg*l+c;
			if (point[p]>gradmax && c>0 && c<Larg) {gradmax=point[p];pgrad=p;}
		}
		
		AJUSTG+=p1-pgrad;

	}
	else
	{
		for (p=p1-AJUSTG;p<p1;p++) bram[p]=255;
	}

/* Droite */
    gradmax=0;
    pgrad=p2;
    if ((c2-AJUST)<c1) ci=c1; else ci=c2-AJUST;
//        for (c=ci;c<=c2+AJUST;c++)
    if (!Flag)
	{
		for (c=c2;c<=c2+AJUST;c++)
		{
			p=Larg*l+c;
			if (point[p]>gradmax && c>0 && c<Larg) {gradmax=point[p];pgrad=p;}
		}
		AJUSTD+=pgrad-p2;
	}
	else
	{
		for (p=p2;p<=p2+AJUSTD;p++) bram[p]=255;
	}
//	CString Stemp;
//	Stemp.Format("AjustD=%d		AjustG=%d",AJUSTD,AJUSTG);
//	MessageBox(Stemp);

}

void CTraitImageWnd::ajusty(int l1,int l2,int c,byte *point,byte *bram)
{

		byte gradmax=0;
		word pgrad,p,p1,p2;
		int l,li;
        p1=l1*Larg+c;
        p2=l2*Larg+c;
        pgrad=p1;
		
		/*bas*/

		for (l=l1-AJUSTH;(l<=l1 && l<l2);l++)
        {
        	p=Larg*l+c;
        	if (point[p]>gradmax && l>0 && l<Haut) {gradmax=point[p];pgrad=p;}
        }
		if (pgrad<p1)
        {
        	for (p=pgrad;p<p1;p+=Larg) bram[p]=255;
        }      

        
		/*haut*/
		pgrad=p2;
        gradmax=0;
		
        //if ((l2-AJUST)<l1) li=l1; else li=l2-AJUST;

		for (l=l2;l<=l2+AJUSTH;l++)
        {
        	p=Larg*l+c;
        	if (point[p]>gradmax && l>0 && l<Haut) {gradmax=point[p];pgrad=p;}
        }

        if (p2<pgrad)
        {
        	for (p=p2;p<=pgrad;p+=Larg) bram[p]=255;
        }      

}// ajusty




void CTraitImageWnd::OnTraitementAnalysecontour() 
{
	// TODO: Add your command handler code here
		long c,i,l,p,k=0,temp,courant,dir,cont;
	long * sommet;
	unsigned char * contour;
	unsigned char *ImageSommet;
	long * compteur;
	long memocontour[10000];
	char ch[90];
	BOOL check;

	sommet=new long [Larg*Haut];

	compteur=new long[Larg*Haut/4 +1];



	BOOL Flag=FALSE;
	for (c=0;c<Larg;c++) 
	{	
		p=0;
		for (l=0;l<Haut;l++) if (Image[l*Larg+c]>20) p++; // compte le nombre de lignes dont lum>seuil
		for (l=0;l<Haut;l++) 
		{
			if (p<(Haut/3))	Image[l*Larg+c]=0; // si on a moins d'un tiers de la hauteur =mise à zéro des pattes (à voir)
			else Flag=TRUE; // on s'arrete au premier
		}
		if (Flag) break; // terminé
	}

/*
	CGraphWnd * pMoy = (CGraphWnd*) new CGraphWnd(Moy,Larg,"Nombre","Niveau");  
	if (!pMoy->Create("Moyenne: "+Titre,GetMDIFrame())) return;
*/

	CTraitImageWnd * pImageSom =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,NB);   
	ImageSommet=pImageSom->Image;
	
	CTraitImageWnd * pImageContour =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,NB);   
	
	contour=pImageContour->Image;

	CTraitImageWnd * pImageContoursupseuil =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,NB);   

	compteur[0]=0;

	ArreteLeContinu();
	Memorise();

	OnMoyenneurxy111();
	OnMoyenneurxy111();
	
	
	////boite de dialogue des seuils
/*	CDialContour Dilg(Seuilcont1,Seuilcont2,LongmaxContour,SeuillongContour,this);
	if (Dilg.DoModal()==IDOK)
	{
		Seuilcont1=(unsigned char) Dilg.m_seuilhaut;
		Seuilcont2=(unsigned char)Dilg.m_seuilbas;
		LongmaxContour=Dilg.m_longmax;
		SeuillongContour=Dilg.m_seuilcont;
		check=Dilg.m_check;
	
	}
	else return;
*/
		Seuilcont1=50;//(unsigned char) Dilg.m_seuilhaut;
		Seuilcont2=17;//(unsigned char)Dilg.m_seuilbas;
		LongmaxContour=1500;//Dilg.m_longmax;
		SeuillongContour=500;//Dilg.m_seuilcont;
		check=FALSE;
	

	SetCursor(LoadCursor(0,IDC_WAIT));

	////tour de l'image à 0
	for(i=0;i<Larg;i++) 
	{
		Image[i]=0;
		Image[(Haut-1)*Larg+i]=0;
	}
	for(i=0;i<Haut;i++)
	{
		Image[i*Larg]=0;
		Image[(i+1)*Larg-1]=0;
	}
	/////initialisation images

	for ( i=0; i< Larg*Haut; i++) 
	{
		contour[i]=0;
		ImageSommet[i]=0;
	}
	////recherche sommets
	l=0;
	for (i=0;i< Haut*Larg;i++)
		{
			if(Image[i]>=Seuilcont1)
				{
					sommet[l]=i;
					ImageSommet[i]=255;
					l++;
				}
		}

	
	///recherche contours
	for (i=0;i<l;i++)
		{
			courant=sommet[i];
			
			while(Image[courant-1]>=Seuilcont2) 
				{courant--;};
			temp=courant;

			if (contour[courant]!=255)
			{
				contour[courant]=255;
				compteur[0]++;
				compteur[compteur[0]]=1;
				dir=0;
				cont=1;
				p=0;
				memocontour[p++]=courant;////pour le seuillage sur la longueur
				

				do
				{
					
					switch(dir+cont)
					{
						case 0:
							courant-=Larg;
							dir=6;
							break;
						case 1:
							courant+=Larg;
							dir=2;
							break;
						case 2:
							courant++;
							dir=0;
							break;
						case 3:
							courant--;
							dir=4;
							break;
						case 4:
							courant+=Larg;
							dir=2;
							break;
						case 5:
							courant-=Larg;
							dir=6;
							break;
						case 6:
							courant--;
							dir=4;
							break;
						case 7:
							courant++;
							dir=0;
							break;
					}
					
					if (Image[courant]>=Seuilcont2)
						{
							cont=1;
							if(contour[courant]!=255)
							{
								compteur[compteur[0]]++;
								contour[courant]=255;
								memocontour[p++]=courant;
								if(compteur[compteur[0]]>=LongmaxContour){MessageBox("contour trop long!");return;}
								
							}
						}
					else cont=0;
				}
				while(courant!=temp);

				if(!check) 
				{
					if(p>=SeuillongContour)
					{	k++;
//						for(long j=0;j<p;j++) pImageContoursupseuil->Image[memocontour[j]]=255;
					}
					
				}
				else
				{
					if(p<SeuillongContour)
					{	k++;
						for(long j=0;j<p;j++) pImageContoursupseuil->Image[memocontour[j]]=255;
					}
				}
				if (k==1) break;
			}////if image=sommet
			if (k==1) break;
			
		}//for l
	
	long *pTabLigne=new long[p];
	long j;
	for( j=0;j<p;j++) 
	{
		pImageContoursupseuil->Image[memocontour[j]]=255;
		pTabLigne[j]=memocontour[j]/Larg;
	}

	for(j=0;j<p-1;j++) 
	{
			pTabLigne[j]=(pTabLigne[j]+pTabLigne[j+1])/2;
	}
	
	/*detection max*/
	long TabMax[1500];

	long Dif,Adif=0,max=0;
	long xpos,xneg;

	for(j=0;j<p-1;j++) 
	{		
			Dif=pTabLigne[j]-pTabLigne[j+1];
			if (Dif>0) xpos=j+1;
			if (Dif<0) xneg=j+1;
			if (Dif<0 && Adif>0) TabMax[max++]=(xpos+xneg)/2;
			if (Dif!=0) Adif=Dif;
	}
	for (j=0;j<max;j++)
	{
		Image[memocontour[TabMax[j]]]=255;
	}


		/*detection min*/
	long TabMin[1500];
	long min=0;
	for(j=0;j<p-1;j++) 
	{		
			Dif=pTabLigne[j]-pTabLigne[j+1];
			if (Dif>0) xpos=j+1;
			if (Dif<0) xneg=j+1;

			if (Dif>0 && Adif<0) TabMin[min++]=(xpos+xneg)/2;
			if (Dif!=0) Adif=Dif;
	}
	for (j=0;j<min;j++)
	{
		Image[memocontour[TabMin[j]]]=255;
	}

	CGraphWnd * pTour = (CGraphWnd*) new CGraphWnd(pTabLigne,p,"Nombre","Niveau");  
	if (!pTour->Create("Histogramme Masque: "+Titre,GetMDIFrame())) return;


//	pImageSom->Create(Titre+" : Sommets",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
//	pImageContour->Create(Titre+" : contour",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
	
	sprintf(ch,"%d",SeuillongContour);
	pImageContoursupseuil->Create(Titre+" : contours supérieurs à "+ch,WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
	GetMDIFrame()->MDITile(MDITILE_VERTICAL);
	this->SetFocus();
	SetCursor(LoadCursor(0,IDC_ARROW)); 
	sprintf(ch,"nombre de sommets:%d\nnbre total de contours: %d\nnbre de contour de long>=seuil: %d",l,compteur[0],k);
	MessageBox(ch,"contours");
	delete sommet;
	delete compteur;
	delete pTabLigne;
	Invalidate(FALSE);
}

void CTraitImageWnd::OnTraitementEquilibrage() 
{
	// TODO: Add your command handler code here
		CTraitImageWnd * pPente =  (CTraitImageWnd * ) new CTraitImageWnd(Larg,Haut,NB);
		
		int l,c,p,NbPix,m,PremierNonNul;
		long *TabMoy;
		float *TabTemp=new float[Larg*Haut];
		TabMoy=new long[Haut];
		for (l=0;l<Haut;l++) TabMoy[l]=0;

		OnMoyenneurxy111();
		OnMoyenneurxy111();

		BOOL Flag=FALSE;
		for (c=0;c<Larg;c++) 
		{	
			p=0;
			for (l=0;l<Haut;l++) if (Image[l*Larg+c]>20) p++; // compte le nombre de lignes dont lum>seuil
			for (l=0;l<Haut;l++) 
			{
				if (p<(Haut/3))	Image[l*Larg+c]=0; // si on a moins d'un tiers de la hauteur =mise à zéro des pattes (à voir)
				else Flag=TRUE; // on s'arrete au premier
			}
			if (Flag) break; // terminé
		}

		Flag=FALSE;

		for (l=0;l<Haut;l++)
		{	
			NbPix=0;
			for (c=0;c<Larg;c++)
			{
				p=Larg*l+c;
				if (Image[p]>15) 
				{	
					NbPix++;
					TabMoy[l]+=Image[p];
				}
			}
			if (NbPix>50) 
			{
				if (!Flag) {PremierNonNul=l;Flag=TRUE;}
				TabMoy[l]/=NbPix;
			}
			else if (Flag) TabMoy[l]=TabMoy[l-1];
		}

		float min=255,Val;
		float max=0;
	
		for (l=0;l<PremierNonNul;l++)
		{
			TabMoy[l]=TabMoy[PremierNonNul];
		}

		for (l=0;l<Haut-2;l++)  
		{
			TabMoy[l]=(TabMoy[l]+TabMoy[l+1]+TabMoy[l+2])/3;
		}
		

		for (l=0;l<Haut;l++)
		{	
			for (c=0;c<Larg;c++)
			{
				p=Larg*l+c;
				if (Image[p]>15) TabTemp[p]=(float)Image[p]-1.5*(float)TabMoy[l];
				else Image[p]=0;
				if (TabTemp[p]<min) min=TabTemp[p];

			}
		}

		for (l=0;l<Haut;l++)
		{	
			for (c=0;c<Larg;c++)
			{
				p=Larg*l+c;
				Val=TabTemp[p]-min;
				pPente->Image[p]=(BYTE)Val;

			}
		}




		CGraphWnd * pMoy = (CGraphWnd*) new CGraphWnd(TabMoy,Haut,"Nombre","Niveau");  
		if (!pMoy->Create("Moyenne lignes: "+Titre,GetMDIFrame())) return;

		pPente->Create(Titre+" : Pente",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW ,GetMDIFrame());
		delete TabTemp;
}

void CTraitImageWnd::OnDcisionGeneredata() 
{
	// TODO: Add your command handler code here
	ArreteLeContinu();

	int l,c,t;
	double taille,rd;
	long n=Larg*Haut;
	CFile Fichier;

		SetCursor(LoadCursor(0,IDC_WAIT));


		if (!ParamCalcul) 
			for (int i=0;i<NbParam;i++) 
			{
				if(!CalculParam(i)) 
				{
					ParamCalcul=FALSE;
					return;// FALSE;
				}
			}

		CStdioFile MonFic;
		CString ChaineTemp;
		float fval;
		CString NomFicDat;
		NomFicDat=NomFichAppr.Left(NomFichAppr.GetLength()-5)+".dat";
		if(MonFic.Open(NomFicDat,CFile::modeCreate|CFile::modeWrite))
		{
			ChaineTemp.Format("%d\n",NbParam);
			MonFic.WriteString(ChaineTemp);
			ChaineTemp.Format("%d\n",n);
			MonFic.WriteString(ChaineTemp);
			for(t=0;t<n;t++)
			{
				Calc(t);
				for (int ppar=0;ppar<NbParam;ppar++)
				{
					fval=(float)Buffer[ppar]/255.0;
					ChaineTemp.Format("%f ",fval);
					MonFic.WriteString(ChaineTemp);
				}

					MonFic.WriteString("\n");
			}
			MonFic.Close();
		}
		else	MessageBox("Nom du fichier .dat incorrect!","Erreur",MB_ICONSTOP);


		SetCursor(LoadCursor(0,IDC_ARROW));
}
