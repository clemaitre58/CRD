#include "struct.h"
// enum

// struct       

              
//  ********** ROUTINES TRAITEMENTS ********************              
extern 	void   MiroirV		(SImage);  
extern   void  MiroirH		(SImage); 
extern   void  Inversion	(SImage); 
extern   void  Seuillage1	(SImage ,	unsigned char,	BOOL);
extern   void  Seuillage2	(SImage , 	unsigned char, 	unsigned char);

extern   void  Roberts		(SImage ,	unsigned char,	unsigned char);
extern   void  Sobel		(SImage ,	unsigned char,	unsigned char);
extern   void  Prewitt		(SImage ,	unsigned char,	unsigned char);

extern   void  FFT2D		(SImage , 	SImage, SImage);   
extern   long  puissance2		(long );   
extern   long  TestPuissance2	(long ); 

extern   long  Etiquett(SImage ,SImage ,unsigned char ); 
extern   void  CentresGravites(SImage ,SImage , long );

// bricolage
extern	void  RotationD(SImage , SImage );
extern	void  RotationG(SImage , SImage );
extern	void  RotationD(SImageFloat , SImageFloat );
extern	void  RotationG(SImageFloat , SImageFloat );

extern 	void SousEchantillonnageLigneColonne(SImage ,SImage);
extern	void SousEchantillonnageLigne(SImage ,SImage );
extern 	void SousEchantillonnageColonne(SImage ,SImage );  

// Transformée en ondelette
extern void TOD (SFiltre h,SFiltre g,SImage source,SImage result,SImage reconst,SImage cont,BOOL ok_reconst,BOOL ok_cont,BOOL ok_basse,BOOL ok_extra,int echelle);
extern void TOD (SFiltre h,SFiltre g,SImageFloat source,SImageFloat result,SImageFloat reconst,SImageFloat cont,BOOL ok_reconst,BOOL ok_cont,BOOL ok_basse,int echelle);

extern void TOD_EX (SFiltre h,SFiltre g,SImage source,SImage ,SImage result,float *  *tab_v_j,float *  *tab_d_h,float *  *tab_d_v,float * *tab_d_o,int echelle,int cpt);
extern void TOD_EX (SFiltre h,SFiltre g,SImageFloat source,SImageFloat ,SImageFloat result,float *  *tab_v_j,float *  *tab_d_h,float *  *tab_d_v,float * *tab_d_o,int echelle,int cpt);

extern void recop_image(unsigned char  *im,float  *dest, long taille);
extern void recop_image(float  *im,float  *dest, long taille);

extern void min_max(float  *f,float  * min,float * max,long taille);
extern void r_affichage(float  *f,unsigned char  *im,float maxi,float mini,long haut,long larg,int echelle);
extern void r_affichage(float  *f,float  *im,long haut,long larg,int echelle);

extern void masque_convolution_2D(float  *source,float  *dest,SImage info,SMasque masq);


// couleur
extern void LuminanceChrominance(SImageRGB,SImage,SImage,SImage);
extern void SeparationRGB(SImageRGB , SImage , SImage , SImage );	
               
//  ********** ROUTINES ACQUISITION ******************** 
/*  
extern   NBCOLOR    InitialisationCarte();
extern   void   FermeCarte(); 
extern   void   Acquerir(int mode,unsigned char  * Ima, int i1, int j1, int i2, int j2);     
extern   long   ReturnLarg();
extern   long   ReturnHaut();    
extern   void   Snap(); 
extern   void   Freeze(); 
extern   void   Grab();                           
extern   DWORD   SequenceCarte(int NbMaxImageSeq,unsigned char  **Bpile,int Temps, int i1,int j1, int i2,int j2);
extern   void	SelectCamera(short);
*/
//************** Morphologie Mathématique **************
extern void Erode(SImage,int);
extern void Dilate(SImage,int);
//************** Outils BSpline ************************
extern void OnConvBspline(SImageFloat im,SImageFloat dest);
extern void OnConvBsplineFloat(SImage im,SImageFloat dest);
extern void OnGradientBspline(SImageFloat im,SImageFloat dest);
extern void OnDerDirSec(SImageFloat im,SImageFloat dest);
extern void OnLigneZero(SImageFloat im,SImage dest);
extern void OnDecompo(SImageFloat im,SImageFloat dest,int prof);
extern void OnRecomp(SImageFloat im,SImageFloat dest, int prof);
extern void OnTraitOnd(SImageFloat im,SImageFloat dest,int prof,float *m_a);

//--------------------------------------------------------------------------------
//				FFT et outils
//--------------------------------------------------------------------------------

extern   void FFT1D(long q, SComplexe  *TFFT1D_S);
extern   long  TestPuissance2(long valeur);
