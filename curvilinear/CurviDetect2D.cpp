#include <cmath>
#include <ary.h>
#include "struct.h"
#include "petrou.h"
#include "curvidetect2d.h"

CCurviDetect2D::CCurviDetect2D(void)
{
   TabListePoint.clear();
}

CCurviDetect2D::~CCurviDetect2D(void)
{	
   for (TablePointList::iterator it=TabListePoint.begin(); it!=TabListePoint.end();it++)
   {
      delete *it;
   }
   TabListePoint.clear();
} 

void CCurviDetect2D::BordAZero(byte *pImage, int Largeur, int Hauteur)
{
   // --- périmètre de sécurité --------------------------------------------- //
   int k;
   long Taille = Largeur*Hauteur;
   for (k=0;k<Largeur;k++) pImage[k]=0;
   for (k=Taille-Largeur+1;k<Taille;k++) pImage[k]=0;
   for (k=0;k<Hauteur;k++) pImage[k*Largeur]=0;
   for (k=0;k<Hauteur;k++) pImage[k*Largeur+(Largeur-1)]=0;
}

#define MARKED_PIXEL 128

bool CCurviDetect2D::SuiviFreeman(utls::BAry *Image)
{
   // init tableau
   byte FreemanHoraire[] = {0, 2, 4, 6, 1, 3, 5, 7};

   // par sécurité
   BordAZero(Image->el[0],Image->cols(),Image->rows());

   /////////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////////recherche du premier point/////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////
   utls::BAry::pointer ptr=Image->first(), size=Image->last();
   IPoint2D FirstPoint(0,0);
   do 
   {  
      // find first unvisited point in the image
      for (; ptr<size && (*ptr)!=255; ptr++) ;

      // still in the image?
      if (ptr<size)
      {
         // translate offset to x,y
         FirstPoint.y = (ptr - Image->first()         ) / Image->cols();
         FirstPoint.x = (ptr - Image->el[FirstPoint.y]);
         /////////////////////////////////////////////////////////////////////////////////////
         //////////////////////////Etude de la connexité du premier point/////////////////////
         /////////////////////////////////////////////////////////////////////////////////////

         IPoint2D CurrentPoint = FirstPoint;
         bool Rebouclage = false, ConditionArret=true;   

         //allocation d'un pointeur pour ajouter une nouvelle liste à l'array
         ConnexityDegreeAFreemanCode ConnexityFirstPoint;
         ConnexityFirstPoint.ConnexityDegree = 0;
         for(int i=0;i<8;i++)
         {
            IPoint2D NewPoint = CurrentPoint;
            switch (i)
            {
            case 0 : NewPoint.x++; break;
            case 1 :
               {
                  NewPoint.x++;
                  NewPoint.y++;
               } 
               break;		
            case 2 : NewPoint.y++; break;
            case 3 : 
               {
                  NewPoint.x--;
                  NewPoint.y++;
               }
               break;
            case 4 : NewPoint.x--; break;
            case 5 : 
               {
                  NewPoint.x--;
                  NewPoint.y--; 
               }
               break;
            case 6 : NewPoint.y--; break;
            case 7 : 
               {
                  NewPoint.x++;
                  NewPoint.y--;
               }
               break;
            } // switch
            if (Image->el[NewPoint.y][NewPoint.x]==255)
            {
               ConnexityFirstPoint.ConnexityDegree++;
               ConnexityFirstPoint.FreemanCode[i]=true;
            }
            else
               ConnexityFirstPoint.FreemanCode[i]=false;
         } // for i - les 8 directions 

         /////////////////////////////////////////////////////////////////////////////////////
         //////////////////////////////////Parcours pour le suivi Freeman/////////////////////
         /////////////////////////////////////////////////////////////////////////////////////
         if(ConnexityFirstPoint.ConnexityDegree!=0)
         {
            PointList *ListePointContour = new PointList;
            ListePointContour->push_back(FirstPoint);
            CurrentPoint = FirstPoint;
            do
            {		
               int  CodeFreeman = 0, NbZero = 0, NbMoinsUn = 0;
               bool bFinFreeman = false, BoolRetour=false, CandidatTrouve = false;
               ConditionArret=true;
               do	
               {
                  IPoint2D NewPoint = CurrentPoint;
                  switch (FreemanHoraire[CodeFreeman])
                  {
                  case 0 : NewPoint.x++; break;
                  case 1 :
                     {
                        NewPoint.x++;
                        NewPoint.y++;
                     } 
                     break;		
                  case 2 : NewPoint.y++; break;
                  case 3 : 
                     {
                        NewPoint.x--;
                        NewPoint.y++;
                     }
                     break;
                  case 4 : NewPoint.x--; break;
                  case 5 : 
                     {
                        NewPoint.x--;
                        NewPoint.y--; 
                     }
                     break;
                  case 6 : NewPoint.y--; break;
                  case 7 : 
                     {
                        NewPoint.x++;
                        NewPoint.y--;
                     }
                     break;
                  } // fin switch
                  utls::BAry::value v = Image->el[NewPoint.y][NewPoint.x];
                  if (v == 0)            NbZero++;
                  if (v == MARKED_PIXEL) NbMoinsUn++;
                  if (v == 255)
                  {
                     ListePointContour->push_back(NewPoint);
                     CurrentPoint = NewPoint;
                     CandidatTrouve = true;
                     Image->el[NewPoint.y][NewPoint.x] = MARKED_PIXEL;
                     if (CurrentPoint == FirstPoint)
                     {
                        ConditionArret = false; //on sort si la tete = la queue
                        Rebouclage = true; // pas la peine de traiter la connexité
                     }
                     bFinFreeman=true; // pas la peine de continuer				
                  }
                  else// si on a pas trouvé de point blanc
                  {
                     if (CodeFreeman == 7)// et qu'on a balayé tous les freemans alors on sort
                     {
                        if (BoolRetour)
                           ConditionArret = false;
                        else
                        {
                           if(ListePointContour->size() > 1)
                           {                                 
                              BoolRetour=true;
                              Image->el[CurrentPoint.y][CurrentPoint.x] = MARKED_PIXEL;
                              PointList::iterator PositionQueue = ListePointContour->end();// ListePointContour->GetTailPosition();
                              PositionQueue--;
                              PositionQueue--;
                              CurrentPoint = *PositionQueue;
                              CodeFreeman = -1;
                              NbZero=0;
                              NbMoinsUn=0;
                           }
                           else
                           {	
                              ConditionArret = false;
                           }
                        } // fin else
                     }
                  }
                  CodeFreeman++;
                  if (CodeFreeman > 7) bFinFreeman=true;			
               } while (!bFinFreeman);

               if ((NbZero==7) && (CandidatTrouve==false))//extremité
               {
                  Image->el[CurrentPoint.y][CurrentPoint.x] = MARKED_PIXEL;
                  ConditionArret=false;
               }

               if ((NbZero<=7)&&(NbMoinsUn<=7)&&(CandidatTrouve==false))//cas d'un point emprisonné
               {
                  Image->el[CurrentPoint.y][CurrentPoint.x] = MARKED_PIXEL;
                  ConditionArret=false;
               }

               if ((NbZero==8)&&(CandidatTrouve==false))//cas d'un point isolé
               {
                  Image->el[CurrentPoint.y][CurrentPoint.x] = MARKED_PIXEL;
                  ConditionArret=false;
               }
            } while (ConditionArret); //condition d'arrêt

            if (ConnexityFirstPoint.ConnexityDegree>1 && !Rebouclage)
            {
               CurrentPoint = FirstPoint;
               do
               {
                  int  CodeFreeman = 0, NbZero = 0, NbMoinsUn = 0;
                  bool bFinFreeman = false, BoolRetour=false, CandidatTrouve = false;
                  ConditionArret=true;
                  do	
                  {
                     IPoint2D NewPoint = CurrentPoint;
                     switch (FreemanHoraire[CodeFreeman])
                     {
                     case 0 : NewPoint.x++; break;
                     case 1 :
                        {
                           NewPoint.x++;
                           NewPoint.y++;
                        } 
                        break;		
                     case 2 : NewPoint.y++; break;
                     case 3 : 
                        {
                           NewPoint.x--;
                           NewPoint.y++;
                        }
                        break;
                     case 4 : NewPoint.x--; break;
                     case 5 : 
                        {
                           NewPoint.x--;
                           NewPoint.y--; 
                        }
                        break;
                     case 6 : NewPoint.y--; break;
                     case 7 : 
                        {
                           NewPoint.x++;
                           NewPoint.y--;
                        }
                        break;
                     } // fin switch
                     utls::BAry::value v = Image->el[NewPoint.y][NewPoint.x];
                     if (v == 0)            NbZero++;
                     if (v == MARKED_PIXEL) NbMoinsUn++;
                     if (v == 255)
                     {
                        ListePointContour->push_front(NewPoint);
                        CurrentPoint = NewPoint;
                        CandidatTrouve = true;
                        Image->el[NewPoint.y][NewPoint.x] = MARKED_PIXEL;
                        if(CurrentPoint == FirstPoint)
                        {
                           Image->el[FirstPoint.y][FirstPoint.x] = MARKED_PIXEL;
                           ConditionArret = false; //on sort si la tete = la queue
                           Rebouclage = true; // pas la peine de traiter la connexité
                        }
                        bFinFreeman=true; // pas la peine de continuer
                     }
                     else// si on a pas trouvé de point blanc
                     {
                        if(CodeFreeman == 7)// et qu'on a balayé tous les freemans alors on sort
                        {
                           if(BoolRetour)
                              ConditionArret = false;
                           else
                           {
                              if(ListePointContour->size() > 1)
                              {
                                 BoolRetour=true;
                                 Image->el[CurrentPoint.y][CurrentPoint.x] = MARKED_PIXEL;
                                 PointList::iterator PositionQueue = ListePointContour->begin();   
                                 PositionQueue++;
                                 CurrentPoint = *PositionQueue;
                                 CodeFreeman = -1;
                                 NbZero=0;
                                 NbMoinsUn=0;
                              }
                              else
                                 ConditionArret = false;
                           }
                        }
                     }
                     CodeFreeman++;
                     if (CodeFreeman >7) bFinFreeman=true;
                  } while (!bFinFreeman);

                  if ((NbZero==7) && (CandidatTrouve==false))//extremité -- end of the edge
                  {
                     Image->el[CurrentPoint.y][CurrentPoint.x] = MARKED_PIXEL;
                     ConditionArret=false;
                  }

                  if ((NbZero<=7)&&(NbMoinsUn<=7)&&(CandidatTrouve==false))//cas d'un point emprisonné  -- problem blocked poimt
                  {
                     Image->el[CurrentPoint.y][CurrentPoint.x] = MARKED_PIXEL;
                     ConditionArret=false;
                  }

                  if ((NbZero==8)&&(CandidatTrouve==false))//cas d'un point isolé -- just one point
                  {
                     Image->el[CurrentPoint.y][CurrentPoint.x] = MARKED_PIXEL;
                     ConditionArret=false;
                  }
               }while (ConditionArret); //condition d'arrêt; //condition d'arrêt

            }//fin traitement cas U

            TabListePoint.push_back(ListePointContour);
         }// fin if point isolé
         else
            Image->el[FirstPoint.y][FirstPoint.x] = MARKED_PIXEL;
      }// fin if (ptr<size)
   } while (ptr<size);
   return true;
}

void CCurviDetect2D::RemplitTableauIJCoupeBiDir(IPoint2D *TabIJ,int x,int y,float dx,float dy,int w,int dw,int Dir,int *IptInit, long Hauteur, long Largeur)
{
   float Pente=0.0F, OrdonneeOrigine=0.0F;
   float Pente2=0.0F, OrdonneeOrigine2=0.0F;
   char equation='0';
   long i,j,l,ip,jp;

   //*IptInit=dw;

   *IptInit=w;

   int NbPointCoupe =(2*(dw+w))+1; // attention modif nov 2007
   Pente=Pente2=0;

   if (dx!=0.0) // pas vertical
   {
      Pente = dy/dx;
      OrdonneeOrigine = (float)y - Pente*(float)x; //initialise l'ordonnee b de y=a*x+b
      equation = '3';
   }
   else equation = '1';

   if (dy!=0.0) // pas horizontal
   {
      Pente2=dx/dy;
      OrdonneeOrigine2=(float)x - Pente2*(float)y;
   }
   else equation = '2';
   if(Pente>100.0) equation = '1';
   if (Pente2>100.0) equation ='2';

   // si dx et dy nuls en meme temps = soucis

   switch (equation)
   {

   case '1':
      /* cas vertical */
      //if(dy>=0.0) Dir=1; else Dir=-1;
      *IptInit=w;
      j=y-w*Dir;
      // last change //-------------------->> A voir
      i=x;
      for (l=0; l<NbPointCoupe; l++)
      {
         ip=i;jp=j;
         if (ip>=Largeur) ip=Largeur-1;
         if (jp>=Hauteur) jp=Hauteur-1;
         if (ip<0) ip=0;
         if (jp<0) jp=0;
         TabIJ[l].x = (unsigned short) ip;
         TabIJ[l].y = (unsigned short) jp;
         j +=Dir;
      }
      break;
   case '2': /*cas horizontal*/
      *IptInit=w;
      i=x-w*Dir;
      j=y;
      for (l=0; l<NbPointCoupe; l++)
      {
         ip=i;jp=j;
         if (i>=Largeur) ip=Largeur-1;
         if (j>=Hauteur) jp=Hauteur-1;
         if (ip<0) ip=0;
         if (jp<0) jp=0;
         TabIJ[l].x = (unsigned short) ip;
         TabIJ[l].y = (unsigned short) jp;
         i +=Dir;

      }
      break;
   case '3':
      if (Pente>0.0)
      {
         if (Pente>1.0)
         {
            j = y-w*Dir;
            i=(long)(Pente2*(float)j+OrdonneeOrigine2);
            for (l=0; l<NbPointCoupe; l++)
            {
               ip=i;jp=j;
               if (j>=Hauteur) {jp=Hauteur-1;ip=(long)(Pente2*(float)jp+OrdonneeOrigine2);}
               if (ip>=Largeur) {ip=Largeur-1;jp=(long)(((float)ip-OrdonneeOrigine2)/Pente2);}
               if (jp<0) {jp=0;ip=(long)OrdonneeOrigine2;}
               if (ip<0) {ip=0;jp=(long)(-OrdonneeOrigine2/Pente2);}
               if (ip<0 || ip>=Largeur || jp<0 || jp>=Hauteur) 
               {ip=0;jp=0;}
               TabIJ[l].x = (unsigned short) ip;
               TabIJ[l].y = (unsigned short) jp;
               if (TabIJ[l].x==x && TabIJ[l].y==y) *IptInit=l;
               j+=Dir;
               i=(long)(Pente2*(float)j+OrdonneeOrigine2);
            }
            // for l
         }
         else
         {
            i = x-w*Dir;
            j =(long)(Pente*(float)i+OrdonneeOrigine);
            for (l=0; l<NbPointCoupe; l++)
            {
               ip=i;jp=j;
               if (i>=Largeur) {ip=Largeur-1;jp=(long)(Pente*(float)ip+OrdonneeOrigine);}
               if (jp>=Hauteur) {jp=Hauteur-1;ip=(long)(((float)jp-OrdonneeOrigine)/Pente);}
               if (ip<0) {ip=0;jp=(long)OrdonneeOrigine;}
               if (jp<0) {jp=0;ip=(long)(-OrdonneeOrigine/Pente);}
               if (ip<0 || ip>=Largeur || jp<0 || jp>=Hauteur) 
               {ip=0;jp=0;}	
               TabIJ[l].x = (unsigned short) ip;
               TabIJ[l].y = (unsigned short) jp;
               if (TabIJ[l].x==x && TabIJ[l].y==y) *IptInit=l;
               i +=Dir ;
               j = (long)(Pente*(float)i+OrdonneeOrigine);
            }
         }
      }
      else
      {
         if (Pente<-1.0)
         {
            j = y+w*Dir;
            i=(long)(Pente2*(float)j+OrdonneeOrigine2);
            for (l=0; l<NbPointCoupe; l++)
            {
               ip=i;jp=j;
               if (j>=Hauteur) {jp=Hauteur-1;ip=(long)(Pente2*(float)jp+OrdonneeOrigine2);}
               if (ip>=Largeur) {ip=Largeur-1;jp=(long)(((float)ip-OrdonneeOrigine2)/Pente2);}
               if (jp<0) {jp=0;ip=(long)OrdonneeOrigine2;}
               if (ip<0) {ip=0;jp=(long)(-OrdonneeOrigine2/Pente2);}
               if (ip<0 || ip>=Largeur || jp<0 || jp>=Hauteur) 
               {ip=0;jp=0;}
               TabIJ[l].x = (unsigned short) ip;
               TabIJ[l].y = (unsigned short) jp;
               if (TabIJ[l].x==x && TabIJ[l].y==y) *IptInit=l;
               j-=Dir;
               i=(long)(Pente2*(float)j+OrdonneeOrigine2);
            }
            // for l
         }
         else
         {
            i = x-w*Dir;
            j = (long)(Pente*(float)i+OrdonneeOrigine);
            for (l=0; l<NbPointCoupe; l++)
            {
               ip=i;jp=j;
               if (i>=Largeur) {ip=Largeur-1;jp=(long)(Pente*(float)ip+OrdonneeOrigine);}
               if (jp>=Hauteur) {jp=Hauteur-1;ip=(long)(((float)jp-OrdonneeOrigine)/Pente);}
               if (ip<0) {ip=0;jp=(long)OrdonneeOrigine;}
               if (jp<0) {jp=0;ip=(long)(-OrdonneeOrigine/Pente);}
               if (ip<0 || ip>=Largeur || jp<0 || jp>=Hauteur) 
               {ip=0;jp=0;}
               TabIJ[l].x = (unsigned short) ip;
               TabIJ[l].y = (unsigned short) jp;
               if (TabIJ[l].x==x && TabIJ[l].y==y) *IptInit=l;
               i+=Dir ;
               j = (long)(Pente*(float)i+OrdonneeOrigine);
            }
         }
      }
      break;
   default:
      break;
   }
}

void CCurviDetect2D::RemplitTableauIJCoupeBiDirOld(
   IPoint2D *TabIJ, int x,int y, float dx,float dy, int w,int dw, int Dir, int *IptInit, long Hauteur, long Largeur)
{
   double Pente=0.0F;
   double Pente2=0.0F;
   char equation='3';
   int i,j,l,d;
   *IptInit=w;

   int NbPointCoupe = 2*w+1;
   Pente=Pente2=0;

   if (dx!=0.0) // pas vertical
      Pente = dy/dx;
   else 
      equation = '1';

   if (dy!=0.0) // pas horizontal
      Pente2=dx/dy;
   else 
      equation = '2'; 

   // si dx et dy nuls en meme temps = soucis
   switch (equation)
   {

   case '1':
      /* cas vertical */
      j=y-w*Dir;
      // last change //-------------------->> A voir
      i=x;
      for (l=0; l<NbPointCoupe; l++)
      {
         if (i>=Largeur) i=Largeur-1;
         if (j>=Hauteur) j=Hauteur-1;
         if (i<0) i=0;
         if (j<0) j=0;
         TabIJ[l].x = i;
         TabIJ[l].y = j;
         j += Dir;
      }
      break;
   case '2': /*cas horizontal*/
      i=x-w*Dir;
      j=y;
      for (l=0; l<NbPointCoupe; l++)
      {
         if (i>=Largeur) i=Largeur-1;
         if (j>=Hauteur) j=Hauteur-1;
         if (i<0) i=0;
         if (j<0) j=0;
         TabIJ[l].x = i;
         TabIJ[l].y = j;
         i += Dir;

      }
      break;
   case '3':
      d = -(w+1)*Dir;
      if (Pente>0.0)
      {
         if (Pente>1.0)
         {
            for (l=0; l<NbPointCoupe; l++)
            {
               j = y+d;
               i = x+int(double(d)*Pente2);
               d += Dir;
               if (i>=Largeur) i=Largeur-1;
               if (j>=Hauteur) j=Hauteur-1;
               if (i<0) i=0;
               if (j<0) j=0;
               TabIJ[l].x = i;
               TabIJ[l].y = j;
            }
         }
         else
         {
            for (l=0; l<NbPointCoupe; l++)
            {
               i = x+d;
               j = y+int(double(d)*Pente);
               d += Dir;
               if (i>=Largeur) i=Largeur-1;
               if (j>=Hauteur) j=Hauteur-1;
               if (i<0) i=0;
               if (j<0) j=0;
               TabIJ[l].x = i;
               TabIJ[l].y = j;
            }
         }
      }
      else
      {
         if (Pente<-1.0)
         {
            for (l=0; l<NbPointCoupe; l++)
            {
               j = y-d;
               i = x-int(double(d)*Pente2);               
               d += Dir;
               if (i>=Largeur) i=Largeur-1;
               if (j>=Hauteur) j=Hauteur-1;
               if (i<0) i=0;
               if (j<0) j=0;
               TabIJ[l].x = i;
               TabIJ[l].y = j;
            }
         }
         else
         {
            for (l=0; l<NbPointCoupe; l++)
            {
               i = x+d;
               j = y+int(double(d)*Pente);
               d += Dir;
               if (i>=Largeur) i=Largeur-1;
               if (j>=Hauteur) j=Hauteur-1;
               if (i<0) i=0;
               if (j<0) j=0;
               TabIJ[l].x = i;
               TabIJ[l].y = j;
            }
         }
      }
      break;
   default:
      break;
   }

} // fin fonction

bool CCurviDetect2D::TestCroisement(IPoint2D x1, IPoint2D x2, IPoint2D xp1, IPoint2D xp2)
{
   bool Cond1 = false, Cond2 = false;

   //test en x
   if(((xp1.x-x1.x)>0)&&((xp2.x-x2.x)>0))
   {
      Cond1=true;
   }
   //test en x
   if(((x1.y-x2.y)>0)&&((xp1.y-xp2.y)>0))
   {
      Cond2=true;
   }
   return (Cond1 && Cond2);
}

void CCurviDetect2D::InvertionPoint(IPoint2D * xp1, IPoint2D * xp2)
{
   IPoint2D Temp;
   Temp.x = xp1->x;
   Temp.y = xp1->y;
   xp1=xp2;
   xp2->x= Temp.x;
   xp2->y= Temp.y;
}