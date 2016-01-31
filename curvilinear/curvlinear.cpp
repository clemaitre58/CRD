#include <list>
#include <assert.h>
#ifdef WIN32
#include <windows.h>
#endif
#include <timeutls.h>
#include <raster.h>
#include <ary.h>
#include "curvlinear.h"

using namespace utls;

double PointDistance(const IPoint2D &pt1, const IPoint2D &pt2)
{
   return sqrt(double((pt1.x-pt2.x)*(pt1.x-pt2.x)+(pt1.y-pt2.y)*(pt1.y-pt2.y)));
}

void BuildImageFromList(IAry *plan1, CurvilinearRegion *ListeCurvi, IPoint2D &bbmin, IPoint2D &bbmax)
{
   CurvilinearRegion::iterator pos = ListeCurvi->begin(), pos1;
   if (pos != ListeCurvi->end())
   {
      // at least one point?
      pos1 = pos; pos1++;
   } else
      return;
   bbmin.x = min(pos->PosInit.x, pos->PosFinal.x);
   bbmax.x = max(pos->PosInit.x, pos->PosFinal.x);
   bbmin.y = min(pos->PosInit.y, pos->PosFinal.y);
   bbmax.y = max(pos->PosInit.y, pos->PosFinal.y);
   while (pos1 != ListeCurvi->end())
	{
      bbmin.x = min(bbmin.x, pos1->PosInit.x, pos1->PosFinal.x);
      bbmax.x = max(bbmax.x, pos1->PosInit.x, pos1->PosFinal.x);
      bbmin.y = min(bbmin.y, pos1->PosInit.y, pos1->PosFinal.y);
      bbmax.y = max(bbmax.y, pos1->PosInit.y, pos1->PosFinal.y);

      if (triangleArea2D(pos->PosInit, pos1->PosInit, pos1->PosFinal) < 0)
         fill_small_triangle(plan1, pos->PosInit, pos1->PosInit, pos1->PosFinal, IAry::value(1));
      else
         fill_small_triangle(plan1, pos->PosInit, pos1->PosFinal, pos1->PosInit, IAry::value(1));

      if (triangleArea2D(pos1->PosFinal, pos->PosFinal, pos->PosInit) < 0)
         fill_small_triangle(plan1, pos1->PosFinal, pos->PosFinal, pos->PosInit, IAry::value(1));
      else
         fill_small_triangle(plan1, pos1->PosFinal, pos->PosInit, pos->PosFinal, IAry::value(1));

      pos++; pos1++;
	}
   // Check image boundaries
   bbmin.x = min(max(plan1->lb2, bbmin.x), plan1->ub2+1);
   bbmax.x = max(min(bbmax.x, plan1->ub2+1), plan1->lb2);
   bbmin.y = min(max(plan1->lb1, bbmin.y), plan1->ub1+1);
   bbmax.y = max(min(bbmax.y, plan1->ub1+1), plan1->lb1);
}

template <class PointType, class PointType2>
void PointVectorStatsAndBoundingBox(
                        const std::vector<PointType> &tab, 
                        PointType  &bbmin,
                        PointType  &bbmax,
                        PointType2 &Mean, 
                        PointType2 &Variance, 
                        double     &Covariance)
{
   typename std::vector<PointType>::const_iterator it = tab.begin();
   Mean.x = Mean.y = 0;
   Variance.x = Variance.y = 0;
   Covariance = 0;
   double cnt = 0;

   if (it != tab.end())
   {
      bbmin.x = bbmax.x = it->x;
      bbmin.y = bbmax.y = it->y;
   }
   for (it = tab.begin(); it != tab.end(); it++)
   {
      // accumulate stats
      Mean.x     += it->x; Mean.y += it->y;
      Variance.x += it->x * it->x; Variance.y += it->y * it->y;
      Covariance += it->x * it->y;
      
      // update bounding box
      if (bbmin.x > it->x) bbmin.x = it->x;
      if (bbmin.y > it->y) bbmin.y = it->y;
      if (bbmax.x < it->x) bbmax.x = it->x;
      if (bbmax.y < it->y) bbmax.x = it->y;

      cnt ++;
   }
   // var X = E X*X - (E X) * (E X);
   Mean.x /= cnt; Mean.y /= cnt; Variance.x /= cnt; Variance.y /= cnt; Covariance /= cnt;
   Variance.x = Variance.x - Mean.x*Mean.x;
   Variance.y = Variance.y - Mean.y*Mean.y;
   Covariance = Covariance - Mean.y*Mean.x;
}

void MergeConnectedRegions(
                           TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D,
                           int Larg, int Haut)
{
   //-----------------------------------------------------------------------------------
   // merge connected regions
   //-----------------------------------------------------------------------------------
   IAry *plan1, *plan2;
   plan1 = new IAry(Haut, Larg); plan1->clear();
   plan2 = new IAry(Haut, Larg); plan2->clear();
   int *Plan1 = plan1->first(),
       *Plan2 = plan2->first();
   int uTail, uHead;
   int AddSens, ParcoursSens;

   int NbOverlap;
   int NbOverlapB;
   double SeuilDistGrav = 8.0;
   double SeuilLargMoy = 2.0;
   double LargeurA, LargeurB; LargeurA = LargeurB = 0;
   size_t SeuilOverlapp = 3;
   int Taille=Larg*Haut;
   int u;

   double LargMoyA,LargMoyB;

   TableOfCurvilinearRegions::iterator it, jt;

   for (size_t IndiceReg = 0; IndiceReg < TableauListesPointsFiltreesC2D.size(); IndiceReg ++) // for all regions
   {
      IPoint2D bb1min, bb1max;
      CurvilinearRegion *r1 = TableauListesPointsFiltreesC2D[IndiceReg];
      if (r1->merged || r1->OverlappingRegions.size() == 0)
         // if already merged or not overlapping regions, continue with next region
         continue;
      BuildImageFromList(plan1, r1, bb1min, bb1max);

      for (size_t IndiceRegB = 0; IndiceRegB < r1->OverlappingRegions.size(); IndiceRegB ++) // for all box-connected regions
      {
         int i2 = r1->OverlappingRegions[IndiceRegB];
         CurvilinearRegion *r2 = TableauListesPointsFiltreesC2D[i2];
         if (r2->merged)
            // if r2 was already merged, continue with next region
            continue;

         IPoint2D bb2min, bb2max; 
         plan2->clear(); BuildImageFromList(plan2, r2, bb2min, bb2max);

         // compute the gravity center of each axis points in the overlapping region
         DPoint2D Grav1, Grav2;

         Grav1.x=Grav1.y=0.0;
         Grav2.x=Grav2.y=0.0;
         NbOverlap=0;
         NbOverlapB=0;
         ParcoursSens=-1;
         LargMoyA=0;
         LargMoyB=0;

         PosDetectCurvi PosN, PosNHead, PosNTail;
         CurvilinearRegion::iterator pos;

         // find the position of the overlapp : head or tail or neither ?
         pos   = r1->end(); pos--;
         PosNTail  = *pos;
         uTail = Larg*(int)PosNTail.PosAxe.y+(int)PosNTail.PosAxe.x;

         pos   = r1->begin();
         PosNHead  = *pos;
         uHead = Larg*(int)PosNHead.PosAxe.y+(int)PosNHead.PosAxe.x;
         
         AddSens=-1;

         if (Plan2[uTail]!=0 && Plan2[uHead]!=0) // A is in B, no merging with this part
         {
            AddSens=2;
         }
         if (Plan2[uTail]!=0 && Plan2[uHead]==0) // l'axe tail de A est dans B, we need add tail during merging
         {
            AddSens=1; LargeurA = PosNTail.Largeur;
         }
         if (Plan2[uTail]==0 && Plan2[uHead]!=0) // l'axe head de A est dans B, we need add head during merging
         {
            AddSens=0; LargeurA = PosNHead.Largeur;
         }

         // if (AddSens!=-1)
         if (AddSens == 0 || AddSens == 1)
         {
            std::vector<DPoint2D> ListeAxeA;
            while (pos != r1->end())
            {
               PosN = *pos;
               u = Larg*(int)PosN.PosAxe.y+(int)PosN.PosAxe.x;
               LargMoyA += PosN.Largeur;

               if (Plan1[u] != 0 && Plan2[u] != 0) 
                  ListeAxeA.push_back(PosN.PosAxe);

               pos++;
            }
            LargMoyA /= r1->size();

            if (ListeAxeA.size() > SeuilOverlapp) //  if  overlap, try to merge
            {
               pos = r2->end(); pos --;
               PosNTail = *pos;
               uTail = Larg*(int)PosN.PosAxe.y+(int)PosNTail.PosAxe.x;

               pos = r2->begin();
               PosNHead = *pos;
               uHead = Larg*(int)PosN.PosAxe.y+(int)PosNHead.PosAxe.x;	
               
               ParcoursSens=-1;
               
               if (Plan1[uTail]!=0 && Plan1[uHead]!=0) // B is in A, no merging with this part
               {
                  ParcoursSens=2;
               }
               
               if (Plan1[uTail]!=0 && Plan1[uHead]==0) // l'axe tail de B est dans A, we need PREV during merging
               {
                  ParcoursSens=1;
                  LargeurB = PosNTail.Largeur; 
               }
               
               if (Plan1[uTail]==0 && Plan1[uHead]!=0) // l'axe head de B est dans 1, we need NEXT during merging
               {
                  ParcoursSens=0;
                  LargeurB = PosNHead.Largeur;
               }

               // if (ParcoursSens != -1)
               if (ParcoursSens == 0 || ParcoursSens == 1)
               {
                  std::vector<DPoint2D> ListeAxeB;
                  int numpt = 0;
                  pos = r2->begin();
                  while (pos != r2->end())
                  {
                     PosN = *pos;
                     LargMoyB+=PosN.Largeur;
                     u=Larg*(int)PosN.PosAxe.y+(int)PosN.PosAxe.x;
                     if (Plan1[u]!=0 && Plan2[u]!=0) 
                     {
                        ListeAxeB.push_back(PosN.PosAxe);
                        if (ParcoursSens == 0) // on prend la derniere largeur comme LargeurB
                           LargeurB=PosN.Largeur;
                        if (ParcoursSens == 1) // on prend la premiere largeur comme LargeurB
                           if (numpt==0) { LargeurB=PosN.Largeur; numpt++;}
                     }
                     pos++; // on passe à la section suivante
                  }
                  LargMoyB /= r2->size();

                  if (ListeAxeB.size() != 0)
                  {
                     // compute stats of overlapping part of region A
                     DPoint2D GravA, CentreA, VarianceA, CentreMinA, CentreMaxA;
                     double CovarianceA;

                     PointVectorStatsAndBoundingBox(ListeAxeA, CentreMinA, CentreMaxA, GravA, VarianceA, CovarianceA);

                     CentreA.x = (CentreMinA.x + CentreMaxA.x) / 2.0;
                     CentreA.y = (CentreMinA.y + CentreMaxA.y) / 2.0;

                     // compute stats of overlapping part of region B
                     DPoint2D GravB, CentreB, VarianceB, CentreMinB, CentreMaxB;
                     double CovarianceB;
                     
                     PointVectorStatsAndBoundingBox(ListeAxeB, CentreMinB, CentreMaxB, GravB, VarianceB, CovarianceB);

                     CentreB.x = (CentreMinB.x + CentreMaxB.x) / 2.0;
                     CentreB.y = (CentreMinB.y + CentreMaxB.y) / 2.0;

                     double 
                        AlphaA = atan2(CovarianceA, VarianceA.x),
                        AlphaB = atan2(CovarianceB, VarianceB.x),
                        SumDis = fabs(AlphaA - AlphaB) * 180.0 / 3.14159;

                     // if the distance between the gravity centers is low, merge the two lists
                     double DistGrav = utls::l2norm2D(CentreA.x - CentreB.x, CentreA.y - CentreB.y);
                     // double DistLargMoy = fabs(LargMoyB - LargMoyA)/log(utls::max(LargMoyA, LargMoyB));
                     double DistLargMoy = fabs(LargeurB - LargeurA) / log(utls::max(LargeurA, LargeurB));
                     if (
                        DistGrav    < utls::max(LargMoyA, LargMoyB)/4 && 
                        DistLargMoy < SeuilLargMoy &&
                        SumDis      < 45.0 ) //  merge
                     {
                        if (ParcoursSens==0)
                        {
                           // add r2 from head to r1
                           pos = r2->begin();
                           while (pos != r2->end())
                           {
                              PosN = *pos;
                              u = Larg * (int)PosN.PosAxe.y + (int)PosN.PosAxe.x;
                              if (Plan1[u] == 0)  // if the axis is not already in A region
                              {
                                 if (AddSens == 1) r1->push_back(PosN);
                                 if (AddSens == 0) r1->push_front(PosN);
                              }
                              pos++;
                           } // whike pos
                        }

                        if (ParcoursSens==1)
                        {
                           // add r2 from tail to r1
                           CurvilinearRegion::reverse_iterator pos;
                           pos = r2->rbegin();
                           while (pos != r2->rend())
                           {
                              PosN = *pos;
                              u = Larg * (int)PosN.PosAxe.y + (int)PosN.PosAxe.x;
                              if (Plan1[u] == 0)  // if the axis is not already in the A region
                              {
                                 if (AddSens == 1) r1->push_back(PosN);
                                 if (AddSens == 0) r1->push_front(PosN);
                              }
                              pos++;
                           }
                        }

                        // update the plan1 and plan2 to reflect merged r1 & r2
                        plan1->clear(bb1min.y, bb1max.y, bb1min.x, bb1max.x);
                        plan2->clear(bb2min.y, bb2max.y, bb2min.x, bb2max.x);
                        BuildImageFromList(plan1, r1, bb1min, bb1max);
                                                
                        // update overlapping region tables, start from next region
                        for (size_t IndiceReg3 = IndiceReg + 1; IndiceReg3 < TableauListesPointsFiltreesC2D.size(); IndiceReg3++)
                        {
                           CurvilinearRegion *r3 = TableauListesPointsFiltreesC2D[IndiceReg3];
                           // if a reference to r2 is found in some region r3
                           // remove it from r3 and add r3 to r1 to preserve forward ordering of overlapping regions
                           for (IndexVector::iterator jt = r3->OverlappingRegions.begin(); jt != r3->OverlappingRegions.end(); jt++)
                           {
                              if (*jt == i2)
                              {
                                 IndexVector::iterator kt;
                                 // insert r3 to appropriate location in r1->overlapping...
                                 for (kt = r1->OverlappingRegions.begin() + IndiceRegB + 1; kt != r1->OverlappingRegions.end(); kt++)
                                 {
                                    if (*kt >= IndiceReg3)
                                       break;
                                 }
                                 if (kt == r1->OverlappingRegions.end() || *kt != IndiceReg3)
                                    // if there is no reference to r3 already, add it
                                    r1->OverlappingRegions.insert(kt, IndiceReg3);
                              }
                           }
                        }

                        // mark r2 as merged and clear it's data
                        r2->merged = true; r2->clear();
                     }
                  } // end if merging is possible
               } // end if B is in A
            } // end if no overlap
         } // end if no merging
         plan2->clear(bb2min.y, bb2max.y, bb2min.x, bb2max.x);
      }// fin for IndiceRegB
      plan1->clear(bb1min.y, bb1max.y, bb1min.x, bb1max.x);
   }//fin for IndiceReg
   delete Plan1;
   delete Plan2;
}

// Filter of circles and squared shapes
void FilterCirclesAndSquares(TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D, CCurvLinearParams &p)
{  
   int m_LongueurMiniRegionFiltrage = p.LongueurMiniForme/2;
   int LongueurTableauListesPointsFiltreesC2D = TableauListesPointsFiltreesC2D.size();

   for (int IndiceReg=0; IndiceReg<LongueurTableauListesPointsFiltreesC2D; IndiceReg++)//pour toutes les régions
   {
      CurvilinearRegion *r = TableauListesPointsFiltreesC2D[IndiceReg];
      if (r->size()> 0)// si la liste n'est pas vide
      {
         // compute bounding box bas
         r->UpdateAxisBoundingBox();
         if (r->bbsize < r->awidth*1.5)
         {
            TableauListesPointsFiltreesC2D.erase(TableauListesPointsFiltreesC2D.begin()+IndiceReg);
            IndiceReg--;
            LongueurTableauListesPointsFiltreesC2D--;
         }
      } // si longueur ok
   } // for indice
}

// Find overlapping regions
void FindOverlappingRegions(TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D, CCurvLinearParams &p)
{
   TableOfCurvilinearRegions::iterator it, jt;
   for (it = TableauListesPointsFiltreesC2D.begin(); it != TableauListesPointsFiltreesC2D.end(); it++)
   {
      CurvilinearRegion *r1 = *it;
      for (jt = it + 1; jt != TableauListesPointsFiltreesC2D.end(); jt++)
      {
         CurvilinearRegion *r2 = *jt;
         if (!(
            r1->bbmin.x > r2->bbmax.x || r1->bbmax.x < r2->bbmin.x || // x coordinate of r1 is outside r2
            r1->bbmin.y > r2->bbmax.y || r1->bbmax.y < r2->bbmin.y))   // y coordinate of r1 is outside r2
         {
            r1->OverlappingRegions.push_back(int(jt - TableauListesPointsFiltreesC2D.begin()));
            // do only forward references...
            // r2->OverlappingRegions.push_back(r1);
         }
      }
   }
}

void TwoDimConstraints(TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D,
                       TableOfCurvilinearRegions &TableauListesPointsFiltrees,
                       CCurvLinearParams &p)
{
   //------------------------------------------
   // Contraintes 2D
   //------------------------------------------
   int IndiceReg,longueur,cptDepassw=0;
   double LargTemp;
   double ddw;
   double AngleCourb;
   double SeuilCourb=1.0;
   CurvilinearRegion *ListePointsFiltresC2D;
   bool Init = false;
   CurvilinearRegion::iterator pos;
   PosDetectCurvi PosNM1,PosN,PosNP1;


   //PosDetectCurvi *TabPointTemp=(PosDetectCurvi *)malloc(sizeof (PosDetectCurvi)*SeuilDepass);
   int longueurTableauListesPointsFiltrees=TableauListesPointsFiltrees.size();
   for(IndiceReg=0;IndiceReg<longueurTableauListesPointsFiltrees;IndiceReg++)//pour toutes les régions
   {
      pos = TableauListesPointsFiltrees[IndiceReg]->begin();
      Init = false;
      if (int(TableauListesPointsFiltrees[IndiceReg]->size()) > p.TailleRegionSortie)// si la liste n'est pas vide
      {	
         ListePointsFiltresC2D = new CurvilinearRegion;//on alloue une nouvelle liste
         while (pos != TableauListesPointsFiltrees[IndiceReg]->end())
         {
            if (!Init)//initialisation pour avoir n-1
            {
               PosNM1 = *pos; // on récupère ce qui sera n-1
               pos++;         // on pointure le suivant
               ListePointsFiltresC2D->push_front(PosNM1);
               Init = true;
               cptDepassw=0;
               LargTemp=(double)PosNM1.Largeur;

            }
            else // phase de traitement sur la largeur
            {
               PosN = *pos;// On récupère la section n 
               AngleCourb=-10;
               if (pos != TableauListesPointsFiltrees[IndiceReg]->end())// si le pointeur de position ok
               {
                  //PosN=TableauListesPointsFiltrees[IndiceReg]->GetAt(pos);

                  ddw=fabs((double)PosN.Largeur-LargTemp);
                  ddw=ddw/log(LargTemp);

                  pos++;                  
                  if (pos != TableauListesPointsFiltrees[IndiceReg]->end()) 
                  {
                     PosNP1=*pos;
                     pos--;
                     if (ddw>(p.VLocaleLargeur*4.0) )
                     {
                        cptDepassw++;
                        if (cptDepassw >= p.SeuilDepass)
                        {
                           if (int(ListePointsFiltresC2D->size())<p.TailleRegionSortie) 
                              delete ListePointsFiltresC2D;
                           else 
                              TableauListesPointsFiltreesC2D.push_back(ListePointsFiltresC2D);
                           ListePointsFiltresC2D = new CurvilinearRegion; //on alloue une nouvelle liste parce que la section en cours ne respecte pas la contrainte
                           for (int tt=0; tt < p.SeuilDepass; tt++) 
                              pos--;                             
                           Init=false; // retour au debut
                           cptDepassw=0;
                           LargTemp=(double)PosN.Largeur;
                        }
                     } else {
                        // no more need of local curvature neither axe distance
                        ListePointsFiltresC2D->push_front(PosN);//on ajoute la section dans la liste en cours
                        cptDepassw=0;
                        LargTemp=(double)PosN.Largeur;
                     }
                     PosNM1 = PosN; // n-1=n;
                     pos++;
                     //if(pos)PosN = TableauListesPointsFiltrees[IndiceReg]->GetAt(pos);

                  }// si la position est valable
               } // si pos
            }// else si init ok
         }// fin while tant qu'on est dans liste

         longueur=ListePointsFiltresC2D->size();
         //if (longueur>m_LongueurMiniForme)	
         TableauListesPointsFiltreesC2D.push_back(ListePointsFiltresC2D);//on ajoute la liste en cours dans le tableau de liste
      }// si les listes ne sont pas vide
   }//fin for IndiceReg
}

int MeasureOverlap(const int *Ima1,
                   const int *Ima2,
                   int *HistoSurfIm1,
                   int *HistoInter,
                   int Larg,
                   int Haut,
                   int MaxIndex, 
                   double SeuilRecouv, 
                   int &MaxOverlapIndex)
{
   int p, w;
   int HistoSurfIm2;
   int Taille = Larg*Haut;
   for (w = 0; w <= MaxIndex; w++) HistoInter[w] = HistoSurfIm1[w] = 0; 

   HistoSurfIm2 = 0;
   for (p = 0; p < Taille; p++)
   {
      // for all non zero pixels of region2
      if (Ima2[p] != 0) 
      { 
         // compute intersection histogram with regions already in plan1
         HistoInter[Ima1[p]]++;
         // and number of such pixels
         HistoSurfIm2++;
      }
      // also compute number of pixels of regions in plan1
      HistoSurfIm1[Ima1[p]]++;
   }

   int max_overlap = HistoInter[1];
   MaxOverlapIndex = 1;
   for (w = 1; w <= MaxIndex; w++)
   {
      if (MaxOverlapIndex < HistoInter[w])
      {
         max_overlap = HistoInter[w];
         MaxOverlapIndex = w;
      }
   }

   double Measure = (double)max_overlap / (double)(utls::min(HistoSurfIm2, HistoSurfIm1[MaxOverlapIndex]));

   if (Measure > SeuilRecouv && HistoSurfIm1[MaxOverlapIndex] > HistoSurfIm2)
      return 12; // on garde la 1 et on vire la 2

   if (Measure > SeuilRecouv && HistoSurfIm1[MaxOverlapIndex] <= HistoSurfIm2)
      return 21; //on garde la 2 et on vire la 1

   return 0;
}

void FilterOverlappingAndShortRegions( 
                              TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D,
                              CCurvLinearParams &p,
                              int Larg, int Haut)
{
   IAry *plan1, *plan2;
   plan1 = new IAry(Haut, Larg); plan1->clear();
   plan2 = new IAry(Haut, Larg); plan2->clear();
   int *Plan1 = plan1->first(),
       *Plan2 = plan2->first();

   size_t m_LongueurMiniRegionFiltrage = p.LongueurMiniForme/2;
   size_t LongueurTableauListesPointsFiltreesC2D;
   int Taille = Larg*Haut;

   int Intersect=0,Surf1=0,Surf2=0;
   double SeuilRecou=0.7;
   int NbReg1=1;

   
   LongueurTableauListesPointsFiltreesC2D=TableauListesPointsFiltreesC2D.size();
  
   int *HistoSurfIm1 = new int[LongueurTableauListesPointsFiltreesC2D+1];
   int *HistoInter   = new int[LongueurTableauListesPointsFiltreesC2D+1];

   for (size_t IndiceReg = 0; IndiceReg<LongueurTableauListesPointsFiltreesC2D; IndiceReg++)//pour toutes les régions
   {
      CurvilinearRegion *r = TableauListesPointsFiltreesC2D[IndiceReg];
      if (r->size() < m_LongueurMiniRegionFiltrage)
      {
         // kill too short regions (that were not merged)
         r->merged = true; r->clear();
      }
   }

   IPoint2D bb1min, bb1max;
   BuildImageFromList(plan1, TableauListesPointsFiltreesC2D[0], bb1min, bb1max);

   for (size_t IndiceReg = 1; IndiceReg<LongueurTableauListesPointsFiltreesC2D; IndiceReg++)//pour toutes les régions
   {
      CurvilinearRegion *r = TableauListesPointsFiltreesC2D[IndiceReg];

      if (!r->merged)
      {
         IPoint2D bb2min, bb2max;
         BuildImageFromList(plan2, r, bb2min, bb2max);

         // tester recouvrement p1 et p2, et recopier dans p1
         int MaxOverlapIndex = 0;
         int Sens = MeasureOverlap(Plan1, Plan2, HistoSurfIm1, HistoInter, Larg, Haut, IndiceReg, SeuilRecou, MaxOverlapIndex);

         if (Sens == 12) // on garde la 1 et on vire la 2
            r->merged = true;  // keep r1 and mark region r2 as removed (r2 is smaller than overlapping r1)

         if (Sens == 21) // recopie région du plan2 dans plan1 // if (Surf2>Surf1)
         {
            // mark MaxOverlapIndex region as merged
            TableauListesPointsFiltreesC2D[MaxOverlapIndex-1]->merged = true;
            // we will keep r2, first render it to plan1
            for (int u = 0; u < Taille; u++)
            {
               // remove pixels of overlapped region (MaxOverlapIndex)
               if (Plan1[u] == MaxOverlapIndex) Plan1[u] = 0;

               // copy new region to Plan1 with it's index
               if (Plan2[u] != 0) Plan1[u] = IndiceReg + 1;
            }
         } // si recouv

         if (Sens == 0) 
         {
            // no overlap... add region to the plan1
            for (int u=0; u<Taille; u++)
               if (Plan2[u] != 0 && Plan1[u] == 0) 
                  Plan1[u] = IndiceReg + 1;
         } // si peu de recou
         plan2->clear(bb2min.y, bb2max.y, bb2min.x, bb2max.x);
      }// si les listes ne sont pas vide (fin if)
   }//fin for IndiceReg
   // cleanup
   delete [] HistoSurfIm1;
   delete [] HistoInter;
   delete plan1;
   delete plan2;
}

template <typename PointType>
PointType LineProjection(const PointType &A, const PointType &B, const PointType &C, double *pr)
{
	// P is the projection of C onto AB
	PointType P;
	double r, sqnorme;
	sqnorme = (B.x-A.x)*(B.x-A.x)+(B.y-A.y)*(B.y-A.y);
	r = (A.y-C.y)*(A.y-B.y) - (A.x-C.x)*(B.x-A.x);
	r = r/sqnorme;
	P.x = A.x + r*(B.x-A.x);
	P.y = A.y + r*(B.y-A.y);

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

template <typename PointType>
double LinePosition(const PointType &A, const PointType &B, const PointType &C)
{
	// compute the position from A to the projection of C on AB 
	/*
   	r=0 if P=A
	   r=1 if p=B
	   r<0 if P is on backward extension of AB
	   r>0 if P is on forward ext. of AB
	   0<r<1 if P in interior to AB
	*/
	double r, sqnorme;
	sqnorme = (B.x-A.x)*(B.x-A.x)+(B.y-A.y)*(B.y-A.y);
	r = (B.x-A.x)*(C.x-A.x)+(B.y-A.y)*(C.y-A.y);
	r = r / sqnorme;	
	return r;
}

template <typename PointType>
double Variance(const PointType *pTab, size_t NbPt, double Mean)
{
	double Var=0.0;
	for (size_t n=0;n<NbPt;n++)
	{
		Var += (pTab[n].x - Mean) * (pTab[n].x - Mean);
	}
	Var /= (double)NbPt;
	return Var;
}

template <typename PointType>
double Covariance(const PointType *pTab, size_t NbPt, const PointType &Mean)
{
	double Cov=0.0;
	for (size_t n=0; n<NbPt; n++)
		Cov += (pTab[n].x - Mean.x) * (pTab[n].y - Mean.y);
	Cov /= (double)NbPt;
	return Cov;
}

template <typename PointType>
double Distance(const PointType &p1, const PointType &p2)
{
	return sqrt(double((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)));
}


void InvertTable(STabIndice *Tabi, STabIndice *Tabf,int Taille)
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

int CompareDouble(const void *a,const void *b)
{
	STabIndice *val1=(STabIndice *)a;
	STabIndice *val2=(STabIndice *)b;
	if (val1->r > val2->r) return 1; // 1
	else if (val2->r > val1->r) return -1; //-1
	else return 0;
}

void SectionOptimization(TableOfCurvilinearRegions &TableauListesPointsFiltrees, int Larg, int Haut)
{
   PosDetectCurvi PosNM1,PosN,PosNP1;

int i,k;
int LengthLocalLine=15;
int cp=LengthLocalLine/2; // central point
PosDetectCurvi *pLocalLine;

PosDetectCurvi *pOneCurve;

PosDetectCurvi *pCurveTemp;

DPoint2D GravInit, GravFinal;

 

int kcp;
for (TableOfCurvilinearRegions::iterator it = TableauListesPointsFiltrees.begin(); it != TableauListesPointsFiltrees.end(); it++)// pour toutes les régions
{

CurvilinearRegion::iterator pos = (*it)->begin();

int RegionSize=(*it)->size();
pOneCurve =

new PosDetectCurvi[RegionSize];
pCurveTemp =

new PosDetectCurvi[RegionSize+LengthLocalLine];
for (i=0; i<cp; i++) pCurveTemp[i] = *pos;
for (i=cp; i<RegionSize+cp; i++)
{

pCurveTemp[i] = *pos;

/*if(pos != (*it)->end())*/ pos++;
}

pos--;

for (i=RegionSize+cp; i<RegionSize+LengthLocalLine; i++) pCurveTemp[i] = *pos;
for (k=0; k<RegionSize; k++)
{

GravInit.x=GravInit.y=0.0;

GravFinal.x=GravFinal.y=0.0;

for (i=0; i<LengthLocalLine; i++)
{

GravInit.x+=(

double)pCurveTemp[k+i].PosInit.x;
GravInit.y+=(

double)pCurveTemp[k+i].PosInit.y;
GravFinal.x+=(

double)pCurveTemp[k+i].PosFinal.x;
GravFinal.y+=(

double)pCurveTemp[k+i].PosFinal.y;
}

GravInit.x/=LengthLocalLine;

GravInit.y/=LengthLocalLine;

GravFinal.x/=LengthLocalLine;

GravFinal.y/=LengthLocalLine;

pOneCurve[k].PosInit.x=(

int)GravInit.x;
pOneCurve[k].PosInit.y=(

int)GravInit.y;
pOneCurve[k].PosFinal.x=(

int)GravFinal.x;
pOneCurve[k].PosFinal.y=(

int)GravFinal.y;
pOneCurve[k].PosAxe.x=(GravInit.x+GravFinal.x)/2.0;

pOneCurve[k].PosAxe.y=(GravInit.y+GravFinal.y)/2.0;

pOneCurve[k].Largeur=Distance(GravInit,GravFinal);

pOneCurve[k].PositionMaxi=pCurveTemp[k+cp].PositionMaxi;

pOneCurve[k].ValDist=pCurveTemp[k+cp].ValDist;

}

// copy the modified region in the original list
pos = (*it)->begin();

for (i=0; i<RegionSize; i++)
{

*pos= pOneCurve[i];

if(pos != (*it)->end()) pos++;
}

delete [] pOneCurve;
delete [] pCurveTemp;
}

// for all regions
}
