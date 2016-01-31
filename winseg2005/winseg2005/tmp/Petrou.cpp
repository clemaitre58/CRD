// Petrou.cpp: implementation of the CPetrou class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "STRUCT.h"
#include "winseg.h"
#include "Petrou.h"
#include "regress.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
/* Voici differents masque de lissage utilise par J.Kittler et M.Petrou   */

float masquepetrou1[11][11] =
{{0       , 0       , 0.00049 , 0.02135 , 0.06300 , 0.08423 , 0.06300 , 0.02135 , 0.00049 ,  0      , 0  },
{0       , 0.00394 , 0.08423 , 0.26622 , 0.45502 , 0.53485 , 0.45502 , 0.26622 , 0.08423 , 0.00394 , 0  }
,{0.00049 , 0.08423 , 0.38409 , 0.83549 , 1.24773 , 1.14549 , 1.24773 , 0.83549 , 0.38409 , 0.08423 , 0.00049 }
,{0.02135 , 0.26622 , 0.83549 , 1.60099 , 2.28535 , 2.56577 , 2.28535 , 1.60099 , 0.83549 , 0.26622 , 0.02135 }
,{0.06300 , 0.45502 , 1.24773 , 2.28535 , 3.23228 , 3.63454 , 3.23228 , 2.28535 , 1.24773 , 0.45502 , 0.06300 }
,{0.08423 , 0.53485 , 1.14549 , 2.56577 , 3.63454 , 4.10412 , 3.63454 , 2.56577 , 1.14549 , 0.53485 , 0.08423 }
,{0.06300 , 0.45502 , 1.24773 , 2.28535 , 3.23228 , 3.63454 , 3.23228 , 2.28535 , 1.24773 , 0.45502 , 0.06300 }
,{0.02135 , 0.26622 , 0.83549 , 1.60099 , 2.28535 , 2.56577 , 2.28535 , 1.60099 , 0.83549 , 0.26622 , 0.02135 }
,{0.00049 , 0.08423 , 0.38409 , 0.83549 , 1.24773 , 1.14549 , 1.24773 , 0.83549 , 0.38409 , 0.08423 , 0.00049 }
,{0       , 0.00394 , 0.08423 , 0.26622 , 0.45502 , 0.53485 , 0.45502 , 0.26622 , 0.08423 , 0.00394 , 0 }
,{0       , 0       , 0.00049 , 0.02135 , 0.06300 , 0.08423 , 0.06300 , 0.02135 , 0.00049 ,  0      , 0 }};


float masquepetrou3[11][11] = 
{ 
  { 0       , 0       , 0.00053 , 0.02298 , 0.06760 , 0.09028 , 0.06760 , 0.02298 , 0.00053 , 0       , 0},
  { 0       , 0.00425 , 0.09028 , 0.28319 , 0.48086 , 0.56374 , 0.48086 , 0.28319 , 0.09028 , 0.00425 , 0},
  { 0.00053 , 0.09028 , 0.40687 , 0.87229 , 1.28632 , 1.45179 , 1.28632 , 0.87229 , 0.40687 , 0.09028 , 0.00053},
  { 0.02298 , 0.28319 , 0.87229 , 1.63271 , 2.28095 , 2.53746 , 2.28095 , 1.63271 , 0.87229 , 0.28319 , 0.02298},
  { 0.06760 , 0.48086 , 1.28632 , 2.28095 , 3.12403 , 3.46058 , 3.12403 , 2.28095 , 1.28632 , 0.48086 , 0.06760},
  { 0.09028 , 0.56374 , 1.45179 , 2.53746 , 3.46058 , 3.83311 , 3.46058 , 2.53746 , 1.45179 , 0.56374 , 0.09028},
  { 0.06760 , 0.48086 , 1.28632 , 2.28095 , 3.12403 , 3.46058 , 3.12403 , 2.28095 , 1.28632 , 0.48086 , 0.06760},
  { 0.02298 , 0.28319 , 0.87229 , 1.63271 , 2.28095 , 2.53746 , 2.28095 , 1.63271 , 0.87229 , 0.28319 , 0.02298},
  { 0.00053 , 0.09028 , 0.40687 , 0.87229 , 1.28632 , 1.45179 , 1.28632 , 0.87229 , 0.40687 , 0.09028 , 0.00053},
  { 0       , 0.00425 , 0.09028 , 0.28319 , 0.48086 , 0.56374 , 0.48086 , 0.28319 , 0.09028 , 0.00425 , 0},
  { 0       , 0       , 0.00053 , 0.02298 , 0.06760 , 0.09028 , 0.06760 , 0.02298 , 0.00053 , 0       , 0}
};


float masquepetrou4[11][11] =
{ 
  {0       , 0       , 0.00053 , 0.02324 , 0.06837 , 0.09129 , 0.06837 , 0.02324 , 0.00053 , 0       , 0},
  {0       , 0.00430 , 0.09129 , 0.28605 , 0.48251 , 0.56859 , 0.48251 , 0.28605 , 0.09129 , 0.00430 , 0},
  {0.00053 , 0.09129 , 0.41071 , 0.87841 , 1.29256 , 1.45757 , 1.29256 , 0.87841 , 0.41071 , 0.09129 , 0.00053},
  {0.02324 , 0.28605 , 0.87841 , 1.63761 , 2.27956 , 2.53211 , 2.27956 , 1.63761 , 0.87841 , 0.28605 , 0.02324},
  {0.06837 , 0.48251 , 1.29256 , 2.27956 , 3.10617 , 3.43309 , 3.10617 , 2.27956 , 1.29256 , 0.48251 , 0.06837},
  {0.09129 , 0.56859 , 1.45757 , 2.53211 , 3.43309 , 3.79245 , 3.43309 , 2.53211 , 1.45757 , 0.56859 , 0.09129}, 
  {0.06837 , 0.48251 , 1.29256 , 2.27956 , 3.10617 , 3.43309 , 3.10617 , 2.27956 , 1.29256 , 0.48251 , 0.06837},
  {0.02324 , 0.28605 , 0.87841 , 1.63761 , 2.27956 , 2.53211 , 2.27956 , 1.63761 , 0.87841 , 0.28605 , 0.02324},
  {0.00053 , 0.09129 , 0.41071 , 0.87841 , 1.29256 , 1.45757 , 1.29256 , 0.87841 , 0.41071 , 0.09129 , 0.00053},
  {0       , 0.00430 , 0.09129 , 0.28605 , 0.48251 , 0.56859 , 0.48251 , 0.28605 , 0.09129 , 0.00430 , 0},
  {0       , 0       , 0.00053 , 0.02324 , 0.06837 , 0.09129 , 0.06837 , 0.02324 , 0.00053 , 0       , 0},
};

CPetrou::CPetrou()
{
	m_TE = 3;
	m_MaxLength = 3000;
	m_MaxWidth = 6;
	m_MaxDw = 6;
	m_CheckDw = FALSE;
	m_Check_Reg=FALSE;
	m_MinLength = 1;
	m_MaxReg=0.8;
}

CPetrou::~CPetrou()
{
	/*int i;
	for (i=0;i<A.dim_y;i++)
	{
		free(A.ima[i]);
		free(B.ima[i]);
		free(I_x.ima[i]);
		free(I_y.ima[i]);
	}
	free(A.ima);
	free(B.ima);
	free(I_x.ima);
	free(I_y.ima);*/
}

void CPetrou::lissage(const ImagePetrou* A, ImagePetrou* B, int w, int d)
{
 
  int i, j,p,q;
d=1;  
  /* la taille du masquepetrou est de (2w-1)*(2w-1) */
  /* Pour contourner les problemes de bord nous n appliquons le lissage sur les bords */
  /* c est pourquoi on commence a w-1 Pixelps de bord */ 


	for (i=w-1; i<A->dim_y-w+1; i++)
    {
     for (j=w-1; j<A->dim_x-w+1; j++)
       {
		 if (i<w-1)	B->ima[i][j]=A->ima[i][j];
		 else
		   if (i>=A->dim_y-w+1)	  B->ima[i][j]=A->ima[i][j];
		   else 
			 if (j<w-1)	   B->ima[i][j]=A->ima[i][j];
			 else
				 if (j>=A->dim_x-w+1)	 B->ima[i][j]=A->ima[i][j];
		       else
			 {
			   B->ima[i][j]=0;
			   for (p=-w+1; p<w; p++)
			   {
				  for (q=-w+1; q<w; q++)
				 {
			     
					  switch (d)
					 {
					 case 1 :
					   B->ima[i][j]+=A->ima[i+p][j+q]*(masquepetrou1[w+p-1][w+q-1]/100);
					   break;
					 case 3 :
					   B->ima[i][j]+=A->ima[i+p][j+q]*(masquepetrou3[w+p-1][w+q-1]/100); 
					   break;
					 case 4 :
					   B->ima[i][j]+=A->ima[i+p][j+q]*(masquepetrou4[w+p-1][w+q-1]/100); 
					   break; 
					   }/*produit de convolution*/
				}
		     }
		 }
       }
    }

/*   for (i=w-1; i<A->dim_y-w+1; i++)
     for (j=w-1; j<A->dim_x-w+1; j++)
			B->ima[i][j]=A->ima[i][j];*/

}



void CPetrou::derive_Y(const ImagePetrou* A, ImagePetrou* I)
{
  
  int i, j;
 
 // for (j=0; j<A->dim_x; j++)  {B->ima[0][j]=0;B->ima[B->dim_y-1][j]=0;}

   for (i=1; i<A->dim_y-1; i++)
    for (j=0; j<A->dim_x; j++)
       I->ima[i][j]=A->ima[i-1][j]-A->ima[i+1][j];  /*derive en y au point (i,j)*/
    
}
 
  
void CPetrou::derive_X(const ImagePetrou* A, ImagePetrou* I)
{
  
  int i, j;
 
  //for (i=0; i<A->dim_y; i++)   {B->ima[i][0]=0;B->ima[i][B->dim_x-1]=0;}

   for (i=0; i<A->dim_y; i++)
    for (j=1; j<A->dim_x-1; j++)
      I->ima[i][j]=A->ima[i][j-1]-A->ima[i][j+1];  /*derive en x au point (i,j)*/
}







/*Recherche des extremas locaux dans la direction du gradient */
void CPetrou::extrema(ImagePetrou* B, ImagePetrou* I_x,ImagePetrou* I_y, int seuil , int d, int w)
{
  int i,j;
  int m=0,n=0;
  int q;
  float tangente;
  float a=0.414;
  float b=2.414;
d=1;
  //int compteur=0;
  ImagePetrou* Inter;


  Inter = (ImagePetrou *) malloc(sizeof(ImagePetrou));
  Inter->dim_x = B->dim_x;
  Inter->dim_y = B->dim_y;
  Inter->ima=(Pixelp **) malloc(sizeof(Pixelp*)*Inter->dim_y);
  for (i=0;i<Inter->dim_y;i++)
    Inter->ima[i]=(Pixelp *) malloc(Inter->dim_x*sizeof(Pixelp));

  /*Calcul de la norme du gradient*/ 
  for (i=0;i<Inter->dim_y;i++)
    for (j=0;j<Inter->dim_x;j++)
		B->ima[i][j]=0;
   
  for (i=1;i<Inter->dim_y-1;i++)
    for (j=1;j<Inter->dim_x-1;j++)
      {
		Inter->ima[i][j]=sqrt(I_x->ima[i][j]*I_x->ima[i][j]+I_y->ima[i][j]*I_y->ima[i][j]);
      }




  for (i=w+1;i<B->dim_y-w-1;i++)
    for (j=w+1;j<B->dim_x-w-1;j++)
      {
		B->ima[i][j]=sqrt(I_x->ima[i][j]*I_x->ima[i][j]+I_y->ima[i][j]*I_y->ima[i][j]);
	/*Calcul de la tangente*/

	
    	if (I_x->ima[i][j]==0) 
		{
			if (I_y->ima[i][j]>=0) {m=1;n=0;} else {m=-1;n=0;}
		}
		else 
		{
			tangente=I_y->ima[i][j]/I_x->ima[i][j];
			if (tangente>-a && tangente<a) {m=0; if (I_x->ima[i][j]>=0) n=1; else n=-1;}
			if (tangente>a && tangente<b && I_x->ima[i][j]>=0 && I_y->ima[i][j]>=0) {m=1;n=1;}
			if (tangente>a && tangente<b && I_x->ima[i][j]<0 && I_y->ima[i][j]<0) {m=-1;n=-1;}
			if (fabs(tangente)>b) {n=0; if(I_y->ima[i][j]>=0) m=1;else m=-1;}
			if (tangente>-b && tangente<-a && I_x->ima[i][j]<=0 && I_y->ima[i][j]>=0) {m=1;n=-1;}
			if (tangente>-b && tangente<-a && I_x->ima[i][j]>=0 && I_y->ima[i][j]<=0) {m=-1;n=1;}

		
		}

	
	for (q=1;q<d+1;q++)
	  {  
	    if (Inter->ima[i][j]<Inter->ima[i+q*m][j+q*n] 
			|| Inter->ima[i][j]<Inter->ima[i-q*m][j-q*n])
	      B->ima[i][j]=0; 
	  }

	if (B->ima[i][j]>seuil) 
	  {
	    B->ima[i][j]=255;
	  }
	  else B->ima[i][j]=0;



	} // fin for

  for (i=0;i<Inter->dim_y;i++)
    free(Inter->ima[i]);
  free(Inter->ima);
  free(Inter);

}



void CPetrou::dynamique(ImagePetrou *pImagePetrou,int H,int L)
{
  int i,j;
  float valmin,valmax;
  ImagePetrou *ptr;

  ptr=pImagePetrou;
/*
  ptr=(ImagePetrou *) malloc(sizeof(ImagePetrou));

  ptr->dim_x = L;
  ptr->dim_y = H;
  ptr->ima=(Pixelp **) malloc(sizeof(Pixelp*)*ptr->dim_y);
  for (i=0;i<ptr->dim_y;i++)
    ptr->ima[i]=(Pixelp *) malloc(ptr->dim_x*sizeof(Pixelp));
*/
  
  valmin=pImagePetrou->ima[0][0];
  valmax=pImagePetrou->ima[0][0];

  for (i=0;i<H;i++)
    for (j=0;j<L;j++)
      {
		if (pImagePetrou->ima[i][j]<valmin) valmin=pImagePetrou->ima[i][j];
		if (pImagePetrou->ima[i][j]>valmax) valmax=pImagePetrou->ima[i][j];
	  }
  
  //printf("%f %f\n",valmin,valmax);
  for (i=0;i<H;i++)
    for (j=0;j<L;j++)
      ptr->ima[i][j]=((pImagePetrou->ima[i][j]-valmin)*255.0)/(valmax-valmin);
//      ptr->ima[i][j]=255/(valmax-valmin)*pImagePetrou->ima[i][j]-255*valmin/(valmax-valmin);
  
 // return(ptr);
} 
 
void CPetrou::filtre_petrou(ImagePetrou* A, ImagePetrou* I, int w, int seuil, int d)
{
//  ImagePetrou I_lissee;
  int i, j;

 // I_lissee= (ImagePetrou *) malloc(sizeof(ImagePetrou));
/*
  I_lissee.dim_x = A->dim_x;
  I_lissee.dim_y = A->dim_y;
  I_lissee.ima=(Pixelp **) malloc(sizeof(Pixelp*)*I_lissee.dim_y);
  for (i=0;i<I_lissee.dim_y;i++)
    I_lissee.ima[i]=(Pixelp *) malloc(I_lissee.dim_x*sizeof(Pixelp));

 */
 
  //printf("Lissage \n");
  lissage(A,I, w, d);
  //printf("Derivation en x\n");
  derive_X(I,I_x);
  //printf("Derivation en y\n");
  derive_Y(I,I_y);
/*
   for (i=0;i<I_lissee.dim_y;i++) free(I_lissee.ima[i]);
   free(I_lissee.ima);*/

  /* Extremas locaux du gradient dans la direction du gradient */
 // printf("recherche des extremas\n");
  extrema(I,I_x,I_y,seuil,d,w);
  
}

void CPetrou::filtre_petrou_couleur(ImagePetrou* R,ImagePetrou* G,ImagePetrou* B, ImagePetrou* I, int w, int seuil, int d)
{
  int i, j;
  float grc1,grc2;
  ImagePetrou I_Temp_X;
  I_Temp_X.dim_x = R->dim_x;
  I_Temp_X.dim_y = R->dim_y;
  I_Temp_X.ima=(Pixelp **) malloc(sizeof(Pixelp*)*I_Temp_X.dim_y);
  ImagePetrou I_Temp_Y;
  I_Temp_Y.dim_x = R->dim_x;
  I_Temp_Y.dim_y = R->dim_y;
  I_Temp_Y.ima=(Pixelp **) malloc(sizeof(Pixelp*)*I_Temp_Y.dim_y);

  for (i=0;i<I_Temp_X.dim_y;i++)
  {
    I_Temp_X.ima[i]=(Pixelp *) malloc(I_Temp_X.dim_x*sizeof(Pixelp));
    I_Temp_Y.ima[i]=(Pixelp *) malloc(I_Temp_Y.dim_x*sizeof(Pixelp));
  }
 

  lissage(R,I, w, d);
  derive_X(I,I_x);
  derive_Y(I,I_y);
 // FILE *fp;
 // fp=fopen("c:\\test.raw","wb");
 // for (i=0;i<I_Temp_X.dim_y;i++) 
	//for (j=0;j<I_Temp_X.dim_x;j++)
	//{
	//	BYTE val=(BYTE)I_x->ima[i][j];
	//	fwrite(&val,1,1,fp);
	//}
 // fclose(fp);

  lissage(G,I, w, d);
  derive_X(I,&I_Temp_X);
  derive_Y(I,&I_Temp_Y);

   for (i=0;i<I_Temp_X.dim_y;i++) 
	for (j=0;j<I_Temp_X.dim_x;j++)
	{
		grc1=I_Temp_X.ima[i][j]*I_Temp_X.ima[i][j]+I_Temp_Y.ima[i][j]*I_Temp_Y.ima[i][j];
		grc2=I_x->ima[i][j]*I_x->ima[i][j]+I_y->ima[i][j]*I_x->ima[i][j];
		if (grc1>grc2)
		{
			I_x->ima[i][j]=I_Temp_X.ima[i][j];
			I_y->ima[i][j]=I_Temp_Y.ima[i][j];
		}
	}
  lissage(B,I, w, d);
  derive_X(I,&I_Temp_X);
  derive_Y(I,&I_Temp_Y);

    for (i=0;i<I_Temp_X.dim_y;i++) 
	for (j=0;j<I_Temp_X.dim_x;j++)
	{
		grc1=I_Temp_X.ima[i][j]*I_Temp_X.ima[i][j]+I_Temp_Y.ima[i][j]*I_Temp_Y.ima[i][j];
		grc2=I_x->ima[i][j]*I_x->ima[i][j]+I_y->ima[i][j]*I_x->ima[i][j];
		if (grc1>grc2)
		{
			I_x->ima[i][j]=I_Temp_X.ima[i][j];
			I_y->ima[i][j]=I_Temp_Y.ima[i][j];
		}
	}


  for (i=0;i<I_Temp_X.dim_y;i++) {free(I_Temp_X.ima[i]);free(I_Temp_Y.ima[i]);}

  free(I_Temp_X.ima);
  free(I_Temp_Y.ima);

  /* Extremas locaux du gradient dans la direction du gradient */
 // printf("recherche des extremas\n");
  extrema(I,I_x,I_y,seuil,d,w);
  
}


void CPetrou::Go(BYTE *pImagePetrou, BYTE *dest,int Hauteur,int Largeur,int seuil,int d)
{
  int i,j;
  int w ;


  /*Chargement de l ImagePetrou dans nos variables*/
  A->dim_x=Largeur;
  A->dim_y=Hauteur;
  B->dim_x=Largeur;
  B->dim_y=Hauteur;

  I_x->dim_x = A->dim_x;
  I_x->dim_y = A->dim_y;

  I_y->dim_x = A->dim_x;
  I_y->dim_y = A->dim_y;

  A->ima=(Pixelp**) malloc(sizeof(Pixelp*)*A->dim_y);
  B->ima=(Pixelp **) malloc(sizeof(Pixelp)*B->dim_y);
  I_x->ima=(Pixelp **) malloc(sizeof(Pixelp*)*I_x->dim_y);
  I_y->ima=(Pixelp **) malloc(sizeof(Pixelp*)*I_y->dim_y);

  for (i=0;i<A->dim_y;i++)
  {
    A->ima[i]=(Pixelp *) malloc(A->dim_x*sizeof(Pixelp));
    B->ima[i]=(Pixelp *) malloc(B->dim_x*sizeof(Pixelp));
    I_x->ima[i]=(Pixelp *) malloc(I_x->dim_x*sizeof(Pixelp));
    I_y->ima[i]=(Pixelp *) malloc(I_y->dim_x*sizeof(Pixelp));
 }


  for (i=0;i<Hauteur;i++)
    for (j=0;j<Largeur;j++)
	{
		A->ima[i][j]=(float)pImagePetrou[j+i*Largeur];
		B->ima[i][j]=0;
		I_x->ima[i][j]=0;
		I_y->ima[i][j]=0;
	}

  /* pour le moment s=0 car il ne sert a rien et w=6 */
  w=6;
  filtre_petrou(A,B,w,seuil, d );

 if (d==1) ProfilOblique();
  
  dynamique(B,Hauteur,Largeur);

  for(i=0;i<Hauteur;i++)
    for(j=0;j<Largeur;j++)
		dest[j+i*Largeur]=(BYTE)B->ima[i][j];
  
  
 

}


void CPetrou::RemplitTableauIJCoupe(POINT *TabIJ,int x,int y,float dx,float dy,int w,int *Centre)
{
	//float Pente=0.0F, OrdonneeOrigine=0.0F;
	//char equation='0';
	//long signe;
	//long i,j,l;
	//int wx;
	//int wy;
	//*Centre=w;
	//
	//int NbPointCoupe =2*w+1;
	//
 //
	//if (dx==0) equation = '1';   //la perp est verticale
	//else if(fabs(dy)<1) equation = '2';   // la perp est horizontale
	// else equation = '3';   //equation de la forme y=ax+b
 //   
	// switch (equation)
	//{  
	//    case '1':
	//		*Centre=w;
	//		j=y-w;  // last change
	//		i=x;
 //           for (l=0; l<NbPointCoupe; l++)
 //           {   
	//			if (i>=A->dim_x) i=A->dim_x-1;
	//			if (j>=A->dim_y) j=A->dim_y-1;
	//			if (i<0) i=0;
	//			if (j<0) j=0;
	//			TabIJ[l].x = (unsigned short) i;
 //               TabIJ[l].y = (unsigned short) j;
	//		//	if (TabIJ[l].x==x && TabIJ[l].y==y) *Centre=l;
	//			j ++;
	//		}
	//		break;
	//	case '2':
	//		*Centre=w;
	//		i=x-w;
	//		j=y;
	//		for (l=0; l<NbPointCoupe; l++)
	//		{   
	//			if (i>=A.dim_x) i=A.dim_x-1;
	//			if (j>=A.dim_y) j=A.dim_y-1;
	//			if (i<0) i=0;
	//			if (j<0) j=0;
 //				TabIJ[l].x = (unsigned short) i;
 //               TabIJ[l].y = (unsigned short) j;
	//		//	if (TabIJ[l].x==x && TabIJ[l].y==y) *Centre=l;
	//			i ++;
	//		}
	//		break;
	//	case '3':
	//		Pente = dy/dx; //initialise la pente a
	//		OrdonneeOrigine = y - Pente*x;  //initialise l'ordonnee b de y=a*x+b
 //			if (Pente>0.0)
	//		{

	//			if (Pente>1.0)	
	//			{

	//				j = y-w;
	//				i =  (long)  (((float)j-OrdonneeOrigine)/Pente);
	//				for (l=0; l<NbPointCoupe; l++) 
	//				{   
	//					if (i>=A.dim_x) i=A.dim_x-1;
	//					if (j>=A.dim_y) j=A.dim_y-1;
	//					if (i<0) i=0;
	//					if (j<0) j=0;
	//					TabIJ[l].x = (unsigned short) i;
 //               		TabIJ[l].y = (unsigned short) j;
	//					if (TabIJ[l].x==x && TabIJ[l].y==y) *Centre=l;
	//					j++;   
	//					i = (long)  (((float)j-OrdonneeOrigine)/Pente);   

	//				}	
	//			}
	//			else
	//			{	
	//				i =  x-w;
	//				j = Pente*(float)i+OrdonneeOrigine;
	//				for (l=0; l<NbPointCoupe; l++)         
	//				{
	//					if (i>=A.dim_x) i=A.dim_x-1;
	//					if (j>=A.dim_y) j=A.dim_y-1;
	//					if (i<0) i=0;
	//					if (j<0) j=0;
	//					TabIJ[l].x = (unsigned short) i;
 //              			TabIJ[l].y = (unsigned short) j;
	//					if (TabIJ[l].x==x && TabIJ[l].y==y) *Centre=l;
	//				 	i ++ ;
	//				 	j = (long)(Pente*(float)i+OrdonneeOrigine);
	//				}                                   
	//			}	
	//		}
	//		else
	//		{   
	//			if (Pente<-1.0)
	//			{
	//				j = y+w;
	//				i =  (long)  (((float)j-OrdonneeOrigine)/Pente);
	//				for (l=0; l<NbPointCoupe; l++)   
	//				{  
	//					if (i>=A.dim_x) i=A.dim_x-1;
	//					if (j>=A.dim_y) j=A.dim_y-1;
	//					if (i<0) i=0;
	//					if (j<0) j=0;
	//					TabIJ[l].x = (unsigned short) i;
 //               		TabIJ[l].y = (unsigned short) j;
	//					if (TabIJ[l].x==x && TabIJ[l].y==y) *Centre=l;
	//					j--;   
	//					i = (long) (((float)j-OrdonneeOrigine)/Pente);   
	//				}	
	//			}
	//			else
	//			{
	//				i =  x-w;
	//				j = Pente*(float)i+OrdonneeOrigine;
	//				for (l=0; l<NbPointCoupe; l++)     
	//				{
	//					if (i>=A.dim_x) i=A.dim_x-1;
	//					if (j>=A.dim_y) j=A.dim_y-1;
	//					if (i<0) i=0;
	//					if (j<0) j=0;
	//					TabIJ[l].x = (unsigned short) i;
	//					TabIJ[l].y = (unsigned short) j;
	//					if (TabIJ[l].x==x && TabIJ[l].y==y) *Centre=l;
	//				 	i++ ;
	//				 	j = (long)(Pente*(float)i+OrdonneeOrigine);
	//				}	
	//			}
	//		
	//		}	
	//		break;
	//		
	//	default:
	//		break;
	//}
}

void CPetrou::ProfilOblique()
{

//float pis2=(3.14159/2.0);
////pis2=pis2*0.1;
//int li,co,ii,jj,k,m,n,p,i,j;
//float d,dmax;
////int w=20;
//dmax=20.0;	
//int dw=10;
//int OldWidth,DifWidth;
//
//int wmax=(A.dim_x*m_MaxWidth)/100;
//if (wmax>=A.dim_x) wmax=(A.dim_x*10)/100;
//int LengthMax=m_MaxLength;
//int NbPts=2*(wmax+dw)+1;;
//int NbPtsMax=2*(wmax)+1;
//POINT ***TabProfil;
//
//POINT *Tab;
//int l,c,u,Oldl,Oldc;
//bool Contour;
//float OldSigneDx,SigneDx;
//float Angle, AngleInitial;
//int erreur=0;
//int compteur=0;
//CPoint Axis;
//POINT *tAxis;
//int Width,WidthLeft,WidthRight,Centre;
//BOOL Oppose;
////pis2=pis2-(pis2*0.1);
///*
//float *ProfilInitial;
//ProfilInitial=new float[NbPtsMax];
//
//float *ProfilBase;
//ProfilBase=new float[NbPts];
//float *ProfilCourant;
//ProfilCourant=new float[NbPts];
//*/
//int ImageSize=B.dim_x*B.dim_y;
//Tab= (POINT*) malloc(NbPts*sizeof(POINT));
//TabProfil=(POINT ***) malloc(ImageSize*sizeof(POINT **)); // alloc maxi !!
//int **TabWidth = new int*[ImageSize];
//int *TabLength= new int [ImageSize];
//int Milieu ;
//int Dir=+1;
//BYTE ValGrad;
//int Lim_x=A.dim_x-1;
//int Lim_y=A.dim_y-1;
//
//
//ImagePetrou Inter;
//Inter.dim_x = B.dim_x;
//Inter.dim_y = B.dim_y;
//Inter.ima=(Pixelp **) malloc(sizeof(Pixelp*)*Inter.dim_y);
//for (i=0;i<Inter.dim_y;i++)
//Inter.ima[i]=(Pixelp *) malloc(Inter.dim_x*sizeof(Pixelp));
//
//for(i=0;i<B.dim_y;i++)
//for(j=0;j<B.dim_x;j++)
//	Inter.ima[i][j]=0;//B.ima[i][j];
//CString Temp;
//	float DeltaAngle,DeltaAngleMax=0.5;
//	POINTF GradFin,GradInit;
//
////for (k=0;k<NbPtsMax;k++) ProfilInitial[k]=0;
//
//
//for (li=dw;li<B.dim_y-dw;li++)
//{
//for (co=dw;co<B.dim_x-dw;co++)
// {
//	l=li;c=co;
//	Width=0;
//	if(B.ima[l][c]==255.0) 
//	{
//		B.ima[l][c]=127.0;
//		Contour=SuiviInitial(&l,&c,1);
//		if (Contour)
//		{
//			if (I_x.ima[l][c]>=0)  SigneDx=1 ;
//			else SigneDx=-1;
//			OldSigneDx=SigneDx;
//
//			Width=wmax;
//			RemplitTableauIJCoupe(Tab,c,l,I_x.ima[l][c],I_y.ima[l][c],wmax,&Centre);
//	//		for (k=0;k<2*wmax+1;k++)	 Inter.ima[Tab[k].y][Tab[k].x]=50+2*k;
///*			for (k=0;k<NbPtsMax;k++)
//			{
//				if (Tab[k].y>=0 && Tab[k].y<A.dim_y && Tab[k].x>=0 && Tab[k].x<A.dim_x)
//				{
//					ProfilInitial[k]=B.ima[Tab[k].y][Tab[k].x];
//				//	Inter.ima[Tab[k].y][Tab[k].x]=200;
//				}
//				else ProfilInitial[k]=0;
//			}
//*/
//			// determine orientation gradient initial
//			//AngleInitial=AngleGradient(l,c);
//			// recherche à gauche	
//			//Angle=AngleInitial;
//			GradInit.x=I_x.ima[l][c];
//			GradInit.y=I_y.ima[l][c];
//			Oppose=FALSE;
//			k=Centre-2;
//			while(k>1 && Oppose==FALSE)
//			{
//				k--;
//				ValGrad=IsEdge(Tab[k].y,Tab[k].x,-1);
//				if (ValGrad==255) //ProfilInitial[k]!=0)
//				{
//				//	Angle=AngleGradient(Tab[k].y,Tab[k].x);
//				//	if (fabs(Angle-AngleInitial)>pis2) Oppose=TRUE;
//
//					GradFin.x=I_x.ima[Tab[k].y][Tab[k].x];
//					GradFin.y=I_y.ima[Tab[k].y][Tab[k].x];
//					DeltaAngle=DiffAngle(GradInit,GradFin);
//					if (DeltaAngle<0.5) Oppose=TRUE;
//
//
//
//				}
//			} 
//			if (!Oppose)	WidthLeft=wmax;// pas de contour opposé dans cette direction
//			else WidthLeft=abs(Centre-k);
//			// recherche à droite
//			//Angle=AngleInitial;
//			Oppose=FALSE;
//			k=Centre+2;
//			while(k<NbPtsMax-1 && Oppose==FALSE)
//			{
//				k++;
//				ValGrad=IsEdge(Tab[k].y,Tab[k].x,1);
//				if (ValGrad==255) 
//				{	
//					/*Angle=AngleGradient(Tab[k].y,Tab[k].x);
//					if (fabs(Angle-AngleInitial)>pis2) Oppose=TRUE;*/
//
//					GradFin.x=I_x.ima[Tab[k].y][Tab[k].x];
//					GradFin.y=I_y.ima[Tab[k].y][Tab[k].x];
//					DeltaAngle=DiffAngle(GradInit,GradFin);
//					if (DeltaAngle<0.5) Oppose=TRUE;
//
//				}
//			} 
//			if (Oppose==FALSE)	WidthRight=wmax;// pas de contour opposé dans cette direction
//			else WidthRight=abs(k-Centre);
//
//			
//			if (WidthLeft>WidthRight) 
//			{ 
//				Width=WidthRight ;
//			//	Milieu=Centre+(Width/2);
//				Dir=+1;
//			}
//			else 
//			{
//				Width=WidthLeft;
//			//	Milieu=Centre-(Width/2);
//				Dir=-1;
//			}
//			
//			OldWidth=Width;
//
//			TabLength[compteur]=0;
//
//			if (Width<=wmax) //si bien un contour
//			{
//			u=0;
//			do
//			{
////CStdioFile MonFic("test.txt",CFile::modeWrite|CFile::modeCreate);
//				Oppose=FALSE;
///*
//				if (Milieu>0 && Milieu<NbPtsMax)
//				{
//					Axis.x=Tab[Milieu].x;
//					Axis.y=Tab[Milieu].y;
//				}
//				else {Axis.x=0; Axis.y=0; break;}
//			//	if (Axis.x<Inter.dim_x && Axis.x>0 && Axis.y>0 && Axis.y<Inter.dim_y) Inter.ima[Axis.y][Axis.x]=255;
//
//				NbPts=2*Width+1;
//				RemplitTableauIJCoupe(Tab,Axis.x,Axis.y,I_x.ima[l][c],I_y.ima[l][c],Width,&Centre);
//			//	for (k=0;k<NbPts;k++)	 Inter.ima[Tab[k].y][Tab[k].x]=10+3*k;
//*/				
//				//Oldl=l;
//				//Oldc=c;
//				Contour=SuiviInitial(&l,&c,1);
//				if (Contour)
//				{
//
//					B.ima[l][c]=127.0;
//					//Nettoie(Oldl,Oldc);
//					NbPts=2*(Width+dw)+1;
//					RemplitTableauIJCoupe(Tab,c,l,I_x.ima[l][c],I_y.ima[l][c],Width+dw,&Centre);
//			//	for (k=0;k<NbPts;k++)	 Inter.ima[Tab[k].y][Tab[k].x]=255-10+3*u;
//				
//				  //	for (k=0;k<NbPts;k++)
//				//	{
//				//		if (Tab[k].y>=0 && Tab[k].y<A.dim_y && Tab[k].x>=0 && Tab[k].x<A.dim_x)
//				//			ProfilInitial[k]=B.ima[Tab[k].y][Tab[k].x];
//				//		else ProfilInitial[k]=0;
//			//	Temp.Format("%3.2f\n",ProfilInitial[k]);
//			//	MonFic.WriteString(Temp);
//			//		}
//				//	AngleInitial=AngleGradient(l,c);
//				//	Angle=AngleInitial;
//					GradInit.x=I_x.ima[l][c];
//					GradInit.y=I_y.ima[l][c];
//					Oppose=FALSE;
//					if (I_x.ima[l][c]>=0) SigneDx=1 ;else SigneDx=-1;
//					if (SigneDx!=OldSigneDx) 
//					{
//						Dir=-Dir;
//					}
//					OldSigneDx=SigneDx;
//
//					
//					k=Centre+2*Dir;
//
//					while(Oppose==FALSE && (k>1 && k<(NbPts-1)))
//					{
//						k+=Dir; // on part 2 pixels + loin
//						ValGrad=IsEdge(Tab[k].y,Tab[k].x,Dir);
//						if (ValGrad==255) 
//						{
//							/*Angle=AngleGradient(Tab[k].y,Tab[k].x);
//							if (fabs(Angle-AngleInitial)>pis2) Oppose=TRUE;*/
//
//							GradFin.x=I_x.ima[Tab[k].y][Tab[k].x];
//							GradFin.y=I_y.ima[Tab[k].y][Tab[k].x];
//							DeltaAngle=DiffAngle(GradInit,GradFin);
//							if (DeltaAngle<0.5) Oppose=TRUE;
//
//						}
//
//					} 
//					if (Oppose) 
//					{
//						Width=(k-Centre)*Dir;
//						DifWidth=abs(Width-OldWidth);
//						OldWidth=Width;
//						if (!(m_CheckDw && DifWidth>m_MaxDw))
//						{
//							if (u==0) // si premier point de la region
//							{
//								TabProfil[compteur]=(POINT **) malloc(LengthMax*sizeof(POINT*));
//								if (TabProfil[compteur]==NULL) 
//								{Temp.Format("cpt=%d",compteur);MessageBox(NULL,"erreur",Temp,MB_OK);}
//								TabWidth[compteur]= new int[LengthMax];
//							//	Temp.Format("cpt=%d",compteur);MessageBox(NULL,Temp,Temp,MB_OK);
//								compteur++;
//							}
//							TabProfil[compteur-1][TabLength[compteur-1]]=(POINT*)malloc((Width+1)*sizeof(POINT));
//							if (TabProfil[compteur-1][TabLength[compteur-1]]==NULL) MessageBox(NULL,"erreur","",MB_OK);
//							TabWidth[compteur-1][TabLength[compteur-1]]=Width;
//							int itab=0;
//							if (Dir==1)
//							{	itab=0;
//								for (int t=Centre+1;t<=Centre+Width;t++)
//								{
//								//	Inter.ima[Tab[t].y][Tab[t].x]=10+16*compteur+TabLength[compteur-1];
//									TabProfil[compteur-1][TabLength[compteur-1]][itab].x=Tab[t].x;
//									TabProfil[compteur-1][TabLength[compteur-1]][itab++].y=Tab[t].y;
//								}
//							}
//							else 
//							{	itab=0;
//								for (int t=Centre;t>=Centre-Width;t--)
//								{
//								//Inter.ima[Tab[t].y][Tab[t].x]=255-10+16*compteur+TabLength[compteur-1];
//									TabProfil[compteur-1][TabLength[compteur-1]][itab].x=Tab[t].x;
//									TabProfil[compteur-1][TabLength[compteur-1]][itab++].y=Tab[t].y;
//								}
//							}
//							TabLength[compteur-1]++;
//							if (TabLength[compteur-1]>LengthMax-1) Width=wmax;
//
//							// for regress testing
//							if (TabLength[compteur-1]>=m_MinLength && m_Check_Reg)
//							{
//								tAxis=(POINT * )malloc(sizeof(POINT)*TabLength[compteur-1]);
//								for (int lu=0;lu<TabLength[compteur-1];lu++)
//								{
//									tAxis[lu].x=(TabProfil[compteur-1][lu][0].x+TabProfil[compteur-1][lu][TabWidth[compteur-1][lu]-1].x)/2;
//									tAxis[lu].y=(TabProfil[compteur-1][lu][0].y+TabProfil[compteur-1][lu][TabWidth[compteur-1][lu]-1].y)/2;
//								} // for lu
//								CRegress Regress(tAxis,TabLength[compteur-1]);
//								free(tAxis);
//								if (Regress.r<m_MaxReg-0.3)
//								{
//									Oppose=FALSE; // si pas droit on arrete
//								}
//							// end regress testing
//							}
//
//						}
//						else Oppose=FALSE; // si dw/dt >dwmax on arrete
//					} // if opposite point
//						else 
//						{int kkkk=0;} // debug
//				} // if contour
////				MonFic.Close();
//				u++;
//			} while(u<LengthMax && Width<wmax && Oppose==TRUE);
////				Width=wmax;
//			 
//
//		   } // si forme detectee
//
//
//			 // Temp.Format("  long prof=%d",u);
//			//if (u!=0)  MessageBox(NULL,Temp,"",MB_OK);
//
//
//		} // si point de contour apres suivi initial de 5
//
//	Nettoie();
//	} // si B
////	if (compteur==1) break;
// } // for col
////	if (compteur==1) break;
//} // for ligne
//
//int cptFinal=0;
//srand(compteur);
//for (i=0;i<compteur;i++)
//{
//  int r=(rand()*compteur)/RAND_MAX;
//  if (/*TabLength[i]>TabWidth[i][0] &&*/ TabLength[i]>=m_MinLength)
//  {
//	  DifWidth=abs(TabWidth[i][0]-TabWidth[i][TabLength[i]-1]);
//	  if (!(m_CheckDw && DifWidth>m_MaxDw))
//	  {
//		  tAxis=(POINT * )malloc(sizeof(POINT)*TabLength[i]);
//			
//		  for (u=0;u<TabLength[i];u++)
//			{
//				tAxis[u].x=(TabProfil[i][u][0].x+TabProfil[i][u][TabWidth[i][u]-1].x)/2;
//				tAxis[u].y=(TabProfil[i][u][0].y+TabProfil[i][u][TabWidth[i][u]-1].y)/2;
//		/*	
//				Inter.ima[Axis.y][Axis.x]=compteur;//16+16*compteur;
//				Inter.ima[TabProfil[i][u][0].y][TabProfil[i][u][0].x]=r;//16+16*i;
//				Inter.ima[TabProfil[i][u][TabWidth[i][u]-1].y][TabProfil[i][u][TabWidth[i][u]-1].x]=r;//16+16*i;
//		*/
//		  //		for (j=0;j<TabWidth[i][u];j++)
//		//		{
//		//			Inter.ima[TabProfil[i][u][j].y][TabProfil[i][u][j].x]=16+16*i;
//		//		} 
//		  } // for u
//		 CRegress Regress(tAxis,TabLength[i]);
//		 if (Regress.r>=m_MaxReg)
//		 {
//			 cptFinal++;
//			 for (u=0;u<TabLength[i];u++)
//			{
//				Inter.ima[tAxis[u].y][tAxis[u].x]=compteur;//16+16*compteur;
//				Inter.ima[TabProfil[i][u][0].y][TabProfil[i][u][0].x]=r;//16+16*i;
//				Inter.ima[TabProfil[i][u][TabWidth[i][u]-1].y][TabProfil[i][u][TabWidth[i][u]-1].x]=r;//16+16*i;
//			} // for u
//
//		 } // if regress
//		 free(tAxis);
//	  } // if not triangle
//  } // if length
//}
//
//Temp.Format("Initial Nb Curvilinear regions=%d\nFinal Nb Curvilinear regions=%d",compteur,cptFinal);
//MessageBox(NULL,Temp,"End",MB_OK);
//
///*
//RemplitTableauIJCoupe(Tab,j,i,I_x.ima[i][j],I_y.ima[i][j], w,&centre);
//for (k=0;k<NbPts;k++)	 ProfilBase[k]=A.ima[Tab[k].y][Tab[k].x];
//
//for (ii=i+1;ii>=i-1;ii--)
//{
//	for (jj=j-1;jj<=j+1;jj++)
//	{
//		if (B.ima[ii][jj]==255) { i=ii;j=jj; break;}
//	}
//	if (B.ima[i][j]==255) break;
//}
//
//for (p=0;p<750;p++)
//{
//	RemplitTableauIJCoupe(Tab,j,i,I_x.ima[i][j],I_y.ima[i][j], w,&centre);
//	for (k=0;k<NbPts;k++)	 Inter.ima[Tab[k].y][Tab[k].x]=127;
//	for (k=0;k<NbPts;k++)	 ProfilCourant[k]=A.ima[Tab[k].y][Tab[k].x];
//	Inter.ima[i][j]=B.ima[i][j]=127;
//
//	d=Distance(ProfilBase,ProfilCourant,NbPts);
//	if (d>dmax) break;
//
//
//	for (ii=i+1;ii>=i-1;ii--)
//	{
//		for (jj=j-1;jj<=j+1;jj++)
//		{
//			if (B.ima[ii][jj]==255) { i=ii;j=jj; break;}
//		}
//		if (B.ima[i][j]==255) break;
//	}
//	if (B.ima[i][j]==0) 
//	{
//		erreur++;
//	}
//}
//*/
////  CString Temp;
//
//  for(i=0;i<B.dim_y;i++)
//    for(j=0;j<B.dim_x;j++)
////		B.ima[i][j]=I_x.ima[i][j];
//		B.ima[i][j]=Inter.ima[i][j];
////		B.ima[i][j]=AngleGradient(i,j);
//
//
//	for (i=0;i<Inter.dim_y;i++)    free(Inter.ima[i]);
//  free(Inter.ima);
//
//
//for (i=0;i<compteur;i++)
// for (u=0;u<TabLength[i];u++)	free(TabProfil[i][u]);
//
//for (i=0;i<compteur;i++)
//	free(TabProfil[i]);
//free(TabProfil);
//
// for (i=0;i<compteur;i++)
//	delete TabWidth[i];
//
// 
//delete TabWidth;
//delete TabLength;
//free( Tab);
//
///*delete ProfilInitial;
//delete ProfilBase;
//delete ProfilCourant;*/
// 
}

float CPetrou::Distance(float *ta, float *tb, int n)
{
 float d=0;
 for (int i=0;i<n;i++)
 {
	d=d+abs(ta[i]-tb[i]);
 }
 d=d/n;
 return d;

}

bool CPetrou::SuiviInitial(int *pl, int *pc, int lmax)
{
	int ii,jj,i,j;
	bool Contour;
	i=*pl;
	j=*pc;
	for (int k=0;k<lmax;k++) // on cherche les lmax suivants
	{
		Contour=FALSE;
		for (ii=i+1;ii>=i-1;ii--)
		{
			for (jj=j-1;jj<=j+1;jj++)
			{
				if (B->ima[ii][jj]==255) 
				{
					B->ima[ii][jj]=127.0;
					if (!Contour) {Contour=TRUE;i=ii;j=jj;break;}					
				}
			}
			if (Contour) break;
		}
		if (!Contour) break;
	}
	*pl=i;
	*pc=j;
	return Contour;
	
}

float CPetrou::AngleGradient(int l, int c)
{

	float angle=0;
	float dx=I_x->ima[l][c];
	float dy=I_y->ima[l][c];
	return atan2((double)dy,(double)dx);

}

float CPetrou::DiffAngle(POINTF p1, POINTF p2)
{

	float dotp=p1.x*p2.x+p1.y*p2.y;
	float normp=sqrt(p1.x*p1.x+p1.y*p1.y)*sqrt(p2.x*p2.x+p2.y*p2.y);
	float d;
	if (normp!=0.0) d=(1.0+dotp/normp)/2.0;
	else d=1.0;

	return d;
}

void CPetrou::Nettoie()
{
	int ii,jj,i,j,mi,mj;
	int l,c;
	int cpt=0;
	for (l=0;l<A->dim_y;l++)
	 for (c=0;c<A->dim_x;c++)
	 {
		// if (B.ima[l][c]==127.0) B.ima[l][c]=0;

		 if (B->ima[l][c]==255.0)
		 {	cpt=0;
			for (ii=l+1;ii>=l-1;ii--)
			{
				for (jj=c-1;jj<=c+1;jj++)
				{
						if (jj>0 && ii>0 && jj<B->dim_x && ii<B->dim_y) 
							if (B->ima[ii][jj]==255.0) {mi=ii;mj=jj;cpt++;}// B.ima[ii][jj]=64;//0;
				}
			}
			if (cpt==1)  B->ima[mi][mj]=127;
		 }
	 }
	
}

int CPetrou::IsEdge(int l, int c,int dir)
{
	if (l>0 && c>0 && l<B->dim_y-1 && c<B->dim_x-1)
	{
	if (B->ima[l][c]!=0
	    || B->ima[l+dir][c+dir]!=0 
		|| B->ima[l+dir][c]!=0
		|| B->ima[l][c+dir]!=0)
		return (255);
	else return 0;
	}
	else return 0;

}

CPetrou::CPetrou(ImagePetrou * ImageA, ImagePetrou * ImageB,ImagePetrou * ImageGradX, ImagePetrou * ImageGradY)
{
	int i=0,j=0;
	m_TE = 3;
	m_MaxLength = 3000;
	m_MaxWidth = 6;
	m_MaxDw = 6;
	m_CheckDw = FALSE;
	m_Check_Reg=FALSE;
	m_MinLength = 1;
	m_MaxReg=0.8;

	 I_x =  ImageGradX;
	 I_y =  ImageGradY;
	 A =  ImageA;
	 B =  ImageB;


}
