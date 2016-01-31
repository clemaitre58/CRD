
#include "struct.h"     
#include <windows.h>    
#include <windowsx.h>
// ***************************************************************** 

void masque_convolution_2D(float  *source,float  *dest,SImage info,SMasque masq);	
void TOD_EX (SFiltre h,SFiltre g,SImage source,SImage dest,SImage extrac,float * *tab_v_j,float * *tab_d_h,float * *tab_d_v,float * *tab_d_o,int echelle,int cpt);
void TOD_EX (SFiltre h,SFiltre g,SImageFloat source,SImageFloat dest,SImageFloat extrac,float * *tab_v_j,float * *tab_d_h,float * *tab_d_v,float * *tab_d_o,int echelle,int cpt);

void reset(float  *f, long taille);  
void min_max(float  *f,float  * min,float * max,long taille);

void recop_image(unsigned char  *im,float  *dest, long taille);
void recop_image(float  *im,float  *dest, long taille);
	
void ss_echant_col(float  *im,float   *f,long haut,long larg,int echelle);
void ss_echant_lig(float  *im,float  *f,long haut,long larg,int echelle);      
void sur_echant_col(float  *im,float   *f,long haut,long larg,int echelle);
void sur_echant_lig(float  *im,float   *f,long haut,long larg,int echelle);	
void convolution_ligne(float  *f, float  *f_conv,float  *filtre,int fin, int centre,long haut,long larg,int echelle);
void convolution_colonne(float  *f, float  *f_conv,float  *filtre,int fin, int centre, long haut,long larg,int echelle);
//void masque_convolution_2D(float  *source,float  *dest,SImage info,float *filtre,int centre_l, int centre_c,int f_larg,int f_haut,int echelle);
void affichage(float  *f,unsigned char   *im,long offset,float maxi, float mini,long haut,long larg,int echelle);
void affichage(float  *f,float    *im,long offset,long haut,long larg,int echelle);

void r_affichage(float  *f,unsigned char  *im,float maxi,float mini,long haut,long larg,int echelle);
void r_affichage(float  *f,float  *im,long haut,long larg,int echelle);

void reconstruction(float  *h, int indice_h_centre,int indice_h_fin,float  *g,int indice_g_centre,int indice_g_fin,float  *v_j,float  *d_h,
 float  *d_v,float  *d_o, float  *f_conv,long larg, long haut,int echelle);
void TOD (SFiltre h,SFiltre g,SImage source,SImage result,SImage reconst,SImage cont,BOOL ok_reconst,BOOL ok_cont,BOOL ok_basse,BOOL ok_extra,int echelle);
void TOD (SFiltre h,SFiltre g,SImageFloat source,SImageFloat result,SImageFloat reconst,SImageFloat cont,BOOL ok_reconst,BOOL ok_cont,BOOL ok_basse,int echelle);

