// Petrou.cpp: implementation of the CPetrou class.
//
//////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <cmath>
#include "struct.h"
#include "petrou.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
/* Voici differents masque de lissage utilise par J.Kittler et M.Petrou   */

float masquepetrou1[11][11] =
{
   {0.0f     , 0.0f     , 0.00049f , 0.02135f , 0.06300f , 0.08423f , 0.06300f , 0.02135f , 0.00049f , 0.0f     , 0.0f  },
   {0.0f     , 0.00394f , 0.08423f , 0.26622f , 0.45502f , 0.53485f , 0.45502f , 0.26622f , 0.08423f , 0.00394f , 0.0f  },
   {0.00049f , 0.08423f , 0.38409f , 0.83549f , 1.24773f , 1.14549f , 1.24773f , 0.83549f , 0.38409f , 0.08423f , 0.00049f },
   {0.02135f , 0.26622f , 0.83549f , 1.60099f , 2.28535f , 2.56577f , 2.28535f , 1.60099f , 0.83549f , 0.26622f , 0.02135f },
   {0.06300f , 0.45502f , 1.24773f , 2.28535f , 3.23228f , 3.63454f , 3.23228f , 2.28535f , 1.24773f , 0.45502f , 0.06300f },
   {0.08423f , 0.53485f , 1.14549f , 2.56577f , 3.63454f , 4.10412f , 3.63454f , 2.56577f , 1.14549f , 0.53485f , 0.08423f },
   {0.06300f , 0.45502f , 1.24773f , 2.28535f , 3.23228f , 3.63454f , 3.23228f , 2.28535f , 1.24773f , 0.45502f , 0.06300f },
   {0.02135f , 0.26622f , 0.83549f , 1.60099f , 2.28535f , 2.56577f , 2.28535f , 1.60099f , 0.83549f , 0.26622f , 0.02135f },
   {0.00049f , 0.08423f , 0.38409f , 0.83549f , 1.24773f , 1.14549f , 1.24773f , 0.83549f , 0.38409f , 0.08423f , 0.00049f },
   {0.0f     , 0.00394f , 0.08423f , 0.26622f , 0.45502f , 0.53485f , 0.45502f , 0.26622f , 0.08423f , 0.00394f , 0.0f },
   {0.0f     , 0.0f     , 0.00049f , 0.02135f , 0.06300f , 0.08423f , 0.06300f , 0.02135f , 0.00049f , 0.0f     , 0.0f }
};


float masquepetrou3[11][11] = 
{ 
   { 0.0f     , 0.0f     , 0.00053f , 0.02298f , 0.06760f , 0.09028f , 0.06760f , 0.02298f , 0.00053f , 0.0f     , 0.0f},
   { 0.0f     , 0.00425f , 0.09028f , 0.28319f , 0.48086f , 0.56374f , 0.48086f , 0.28319f , 0.09028f , 0.00425f , 0.0f},
   { 0.00053f , 0.09028f , 0.40687f , 0.87229f , 1.28632f , 1.45179f , 1.28632f , 0.87229f , 0.40687f , 0.09028f , 0.00053f},
   { 0.02298f , 0.28319f , 0.87229f , 1.63271f , 2.28095f , 2.53746f , 2.28095f , 1.63271f , 0.87229f , 0.28319f , 0.02298f},
   { 0.06760f , 0.48086f , 1.28632f , 2.28095f , 3.12403f , 3.46058f , 3.12403f , 2.28095f , 1.28632f , 0.48086f , 0.06760f},
   { 0.09028f , 0.56374f , 1.45179f , 2.53746f , 3.46058f , 3.83311f , 3.46058f , 2.53746f , 1.45179f , 0.56374f , 0.09028f},
   { 0.06760f , 0.48086f , 1.28632f , 2.28095f , 3.12403f , 3.46058f , 3.12403f , 2.28095f , 1.28632f , 0.48086f , 0.06760f},
   { 0.02298f , 0.28319f , 0.87229f , 1.63271f , 2.28095f , 2.53746f , 2.28095f , 1.63271f , 0.87229f , 0.28319f , 0.02298f},
   { 0.00053f , 0.09028f , 0.40687f , 0.87229f , 1.28632f , 1.45179f , 1.28632f , 0.87229f , 0.40687f , 0.09028f , 0.00053f},
   { 0.0f     , 0.00425f , 0.09028f , 0.28319f , 0.48086f , 0.56374f , 0.48086f , 0.28319f , 0.09028f , 0.00425f , 0.0f},
   { 0.0f     , 0.0f     , 0.00053f , 0.02298f , 0.06760f , 0.09028f , 0.06760f , 0.02298f , 0.00053f , 0.0f     , 0.0f}
};


float masquepetrou4[11][11] =
{ 
   {0.0f     , 0.0f     , 0.00053f , 0.02324f , 0.06837f , 0.09129f , 0.06837f , 0.02324f , 0.00053f , 0.0f     , 0.0f},
   {0.0f     , 0.00430f , 0.09129f , 0.28605f , 0.48251f , 0.56859f , 0.48251f , 0.28605f , 0.09129f , 0.00430f , 0.0f},
   {0.00053f , 0.09129f , 0.41071f , 0.87841f , 1.29256f , 1.45757f , 1.29256f , 0.87841f , 0.41071f , 0.09129f , 0.00053f},
   {0.02324f , 0.28605f , 0.87841f , 1.63761f , 2.27956f , 2.53211f , 2.27956f , 1.63761f , 0.87841f , 0.28605f , 0.02324f},
   {0.06837f , 0.48251f , 1.29256f , 2.27956f , 3.10617f , 3.43309f , 3.10617f , 2.27956f , 1.29256f , 0.48251f , 0.06837f},
   {0.09129f , 0.56859f , 1.45757f , 2.53211f , 3.43309f , 3.79245f , 3.43309f , 2.53211f , 1.45757f , 0.56859f , 0.09129f}, 
   {0.06837f , 0.48251f , 1.29256f , 2.27956f , 3.10617f , 3.43309f , 3.10617f , 2.27956f , 1.29256f , 0.48251f , 0.06837f},
   {0.02324f , 0.28605f , 0.87841f , 1.63761f , 2.27956f , 2.53211f , 2.27956f , 1.63761f , 0.87841f , 0.28605f , 0.02324f},
   {0.00053f , 0.09129f , 0.41071f , 0.87841f , 1.29256f , 1.45757f , 1.29256f , 0.87841f , 0.41071f , 0.09129f , 0.00053f},
   {0.0f     , 0.00430f , 0.09129f , 0.28605f , 0.48251f , 0.56859f , 0.48251f , 0.28605f , 0.09129f , 0.00430f , 0.0f},
   {0.0f     , 0.0f     , 0.00053f , 0.02324f , 0.06837f , 0.09129f , 0.06837f , 0.02324f , 0.00053f , 0.0f     , 0.0f},
};


template <typename ImageType>
void lissage(const ImageType* A, ImageType* B, int w, int d)
{
   int i, j, p,q ;
   d=1;  
   /* la taille du masquepetrou est de (2w-1)*(2w-1) */
   /* Pour contourner les problemes de bord nous n appliquons le lissage sur les bords */
   /* c est pourquoi on commence a w-1 Pixelps de bord */ 
   for (i=w-1; i<A->rows()-w+1; i++)
   {
      for (j=w-1; j<A->cols()-w+1; j++)
      {
         if (i<w-1)	B->el[i][j]=A->el[i][j];
         else
            if (i>=A->rows()-w+1)	  B->el[i][j]=A->el[i][j];
            else 
               if (j<w-1)	   B->el[i][j]=A->el[i][j];
               else
                  if (j>=A->cols()-w+1)	 B->el[i][j]=A->el[i][j];
                  else
                  {
                     B->el[i][j]=0;
                     for (p=-w+1; p<w; p++)
                     {
                        for (q=-w+1; q<w; q++)
                        {

                           switch (d)
                           {
                           case 1 :
                              B->el[i][j]+=A->el[i+p][j+q]*(masquepetrou1[w+p-1][w+q-1]/100);
                              break;
                           case 3 :
                              B->el[i][j]+=A->el[i+p][j+q]*(masquepetrou3[w+p-1][w+q-1]/100); 
                              break;
                           case 4 :
                              B->el[i][j]+=A->el[i+p][j+q]*(masquepetrou4[w+p-1][w+q-1]/100); 
                              break; 
                           }/*produit de convolution*/
                        }
                     }
                  }
      }
   }
}


template <typename ImageType>
void derive_Y(const ImageType* A, ImageType* I)
{
   int i, j;
   for (i=1; i<A->rows()-1; i++)
      for (j=0; j<A->cols(); j++)
         I->el[i][j]=A->el[i-1][j]-A->el[i+1][j];  /*derive en y au point (i,j)*/

}

template <typename ImageType>
void derive_X(const ImageType* A, ImageType* I)
{
   int i, j;
   for (i=0; i<A->rows(); i++)
      for (j=1; j<A->cols()-1; j++)
         I->el[i][j]=A->el[i][j-1]-A->el[i][j+1];  /*derive en x au point (i,j)*/
}


CPetrou::CPetrou(const ImagePetrou* A, int w, int seuil, int d)
{
   Output = new ImagePetrou(A->rows(), A->cols());
   I_x    = new ImagePetrou(A->rows(), A->cols());
   I_y    = new ImagePetrou(A->rows(), A->cols());
   compute(A, w, seuil, d);
}

CPetrou::CPetrou(const ImagePetrouRGB* A, int w, int seuil, int d)
{
   Output = new ImagePetrou(A->rows(), A->cols());
   I_x    = new ImagePetrou(A->rows(), A->cols());
   I_y    = new ImagePetrou(A->rows(), A->cols());
   compute(A, w, seuil, d);
}

CPetrou::~CPetrou()
{
   delete I_x;
   delete I_y;
   delete Output;
}

void CPetrou::compute(const ImagePetrou* A, int w, int seuil, int d)
{
   lissage(A, Output, w, d);
   derive_X(Output, I_x);
   derive_Y(Output, I_y);
   extrema(Output, I_x, I_y, seuil, d, w);
}

void CPetrou::compute(const ImagePetrouRGB* RGB, int w, int seuil, int d)
{
   // color version picking up maximum of gradients
   size_t len = RGB->cols()*RGB->rows();
   ImagePetrouRGB *tmp_x  = new ImagePetrouRGB(RGB->rows(), RGB->cols());
   ImagePetrouRGB *tmp_y  = new ImagePetrouRGB(RGB->rows(), RGB->cols());
   ImagePetrouRGB *tmpRGB = new ImagePetrouRGB(RGB->rows(), RGB->cols());

   lissage(RGB, tmpRGB, w, d);
   derive_X(tmpRGB, tmp_x); derive_Y(tmpRGB, tmp_y);

   ImagePetrou::pointer    ox, oy;
   ImagePetrouRGB::pointer ix, iy;

   ix = tmp_x->first();
   iy = tmp_y->first();

   ox = I_x->first();
   oy = I_y->first();

   for (size_t i=0; i<len; i++)
   {
      ImagePetrouRGB::value v;
      // gradient magnitude
      v = utls::sqr(*ix)+utls::sqr(*iy);
      if (v.r >= v.g && v.r >= v.b)
      {
         // red is strongest
         *ox = ix->r; *oy = iy->r;
      } else {
         // red is ruled out
         if (v.g >= v.b)
         { 
            // green is strongest
            *ox = ix->g; *oy = iy->g;
         } else {
            // otherwise blue is strongest
            *ox = ix->b; *oy = iy->b;
         }
      }
      ix++; iy++; ox++; oy++;
   }

   delete tmp_x; delete tmp_y; delete tmpRGB;
   extrema(Output, I_x, I_y, seuil, d, w);
}

/*Recherche des extremas locaux dans la direction du gradient */
void CPetrou::extrema(ImagePetrou* B, ImagePetrou* I_x,ImagePetrou* I_y, int seuil , int d, int w)
{
   int i,j;
   int m=0,n=0;
   int q;
   float tangente;
   float a=0.414f;
   float b=2.414f;
   d=1;
   //int compteur=0;
   ImagePetrou* Inter = new ImagePetrou(B->rows(), B->cols());

   /*Calcul de la norme du gradient*/ 
   B->clear();

   for (i=1;i<Inter->rows()-1;i++)
   {
      for (j=1;j<Inter->cols()-1;j++)
      {
         Inter->el[i][j]=sqrt(I_x->el[i][j]*I_x->el[i][j]+I_y->el[i][j]*I_y->el[i][j]);
      }
   }

   for (i=w+1;i<B->rows()-w-1;i++)
   {
      for (j=w+1;j<B->cols()-w-1;j++)
      {
         B->el[i][j]=sqrt(I_x->el[i][j]*I_x->el[i][j]+I_y->el[i][j]*I_y->el[i][j]);
         /*Calcul de la tangente*/
         if (I_x->el[i][j]==0)
         {
            if (I_y->el[i][j]>=0) {m=1;n=0;} else {m=-1;n=0;}
         }
         else 
         {
            tangente=I_y->el[i][j]/I_x->el[i][j];
            if (tangente>-a && tangente<a) {m=0; if (I_x->el[i][j]>=0) n=1; else n=-1;}
            if (tangente>a && tangente<b && I_x->el[i][j]>=0 && I_y->el[i][j]>=0) {m=1;n=1;}
            if (tangente>a && tangente<b && I_x->el[i][j]<0 && I_y->el[i][j]<0) {m=-1;n=-1;}
            if (fabs(tangente)>b) {n=0; if(I_y->el[i][j]>=0) m=1;else m=-1;}
            if (tangente>-b && tangente<-a && I_x->el[i][j]<=0 && I_y->el[i][j]>=0) {m=1;n=-1;}
            if (tangente>-b && tangente<-a && I_x->el[i][j]>=0 && I_y->el[i][j]<=0) {m=-1;n=1;}
         }

         for (q=1;q<d+1;q++)
         {  
            if (Inter->el[i][j]<Inter->el[i+q*m][j+q*n] 
            || Inter->el[i][j]<Inter->el[i-q*m][j-q*n])
               B->el[i][j]=0; 
         }

         if (B->el[i][j]>seuil) B->el[i][j]=255; else B->el[i][j]=0;

      } // fin for
   }
   delete Inter;
}

void CPetrou::dynamique(ImagePetrou *pImagePetrou,int H,int L)
{
   int i,j;
   float valmin,valmax;
   ImagePetrou *ptr;

   ptr=pImagePetrou;

   valmin=pImagePetrou->el[0][0];
   valmax=pImagePetrou->el[0][0];

   for (i=0;i<H;i++)
      for (j=0;j<L;j++)
      {
         if (pImagePetrou->el[i][j]<valmin) valmin=pImagePetrou->el[i][j];
         if (pImagePetrou->el[i][j]>valmax) valmax=pImagePetrou->el[i][j];
      }

      for (i=0;i<H;i++)
         for (j=0;j<L;j++)
            ptr->el[i][j]=ImagePetrou::value((pImagePetrou->el[i][j]-valmin)*255.0)/(valmax-valmin);
}