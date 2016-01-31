#ifndef __CURVLINEAR_H__
#define __CURVLINEAR_H__

#include <ary.h>
#include <dtypes.h>
#include "struct.h"
#include "curvidetect1d.h"
#include "curvidetect2d.h"
#include "petrou.h"
#include <timeutls.h>

struct CCurvLinearParams 
{
	int LargeurEtude;             // length of the half of the section -> it should approximate the maximum width expected in the image
	int LargeurFourier;           // number of Fourier coefficients for comparison
	int LongueurMiniForme;        // the length of the shortest list, processed
	int LargeurMini;              // minimum width of curvilinear segment
   int TailleRegionSortie;
   int SeuilDepass;           
   int SeuilExtrapol;
	double SeuilReponse;          // threshold on the minimal product of left gradient*right gradient...
	double VLocaleLargeur;        // variation of relative change of width between neighbouring sections
   double EdgeImageSmoothing;    // smoothing of the image before edge detection

   bool Do2dConstraints,
        DoFilterCompact,
        DoMergeRegions,
        DoRemoveDuplicates,
        DoUseColor,
        DoFFTOnEdgesOnly;

   CCurvLinearParams ()
      {
         EdgeImageSmoothing = 3.0f;
         LargeurEtude = 40;
         LargeurFourier = 4;
         LongueurMiniForme = 50;
         LargeurMini = 12;
         TailleRegionSortie = 20;
         SeuilReponse = 0.001;
         SeuilDepass = 10;
         SeuilExtrapol = 10;
         VLocaleLargeur = 0.3;

         Do2dConstraints = true;
         DoFilterCompact = true;
         DoMergeRegions = true;
         DoRemoveDuplicates = true;
         DoUseColor = false;
         DoFFTOnEdgesOnly = true;
         // TODO: default values
      }
};

void BuildImageFromList(utls::IAry *plan1, CurvilinearRegion *ListeCurvi, IPoint2D &bbmin, IPoint2D &bbmax);
void SectionOptimization(TableOfCurvilinearRegions &TableauListesPointsFiltrees, int Larg, int Haut);
void MergeConnectedRegions(TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D, int Larg, int Haut);
void FilterCirclesAndSquares(TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D, CCurvLinearParams &p);
void FindOverlappingRegions(TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D, CCurvLinearParams &p);
void TwoDimConstraints(TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D,
                       TableOfCurvilinearRegions &TableauListesPointsFiltrees,
                       CCurvLinearParams &p);

void FilterOverlappingAndShortRegions(
                              TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D,
                              CCurvLinearParams &p, int Larg, int Haut);

template <typename ImageType>
void DetectCurvLinearPoints(ImageType *Image, CCurvLinearParams &p, void (*Callback)(int),
                            utls::BAry **EdgeImage,
                            TableOfCurvilinearRegions &TableauListesPointsFiltrees,
                            TableOfCurvilinearRegions &TableauListesPointsFiltreesC2D)
{   
   //--------------------------------------------------------
   //			Lissage gaussien pour avoir plus de connexité
   //--------------------------------------------------------
   typedef typename AccumulatorDeduction<typename ImageType::value>::value AccumulatorType;
  
   ImageType *ImageGauss = utls::copy_and_gaussian_blur<ImageType, AccumulatorType>(Image, p.EdgeImageSmoothing);
      
   if (Callback)
      Callback(10);
   
   //--------------------------------------------------------
   //			Contour Petrou
   //--------------------------------------------------------
   CPetrou PetrouContour(ImageGauss, 6, 3, 1);

   if (EdgeImage)
   {
      (*EdgeImage) = utls::convert<ImagePetrou, utls::BAry>(PetrouContour.Output);
      // (*EdgeImage)->write_pgm("e:/test.pgm");
   }
   if (Callback)
      Callback(20);
   
   //--------------------------------------------------------
   //			Recopie Image Contour Petrou pour suivi
   //--------------------------------------------------------
   utls::BAry *ImageFreeman = utls::convert<ImagePetrou, utls::BAry>(PetrouContour.Output);

   int w = p.LargeurEtude;
   assert(w>2*p.LargeurFourier); // this should not happen, otherwise TabCoupe will be too small
   int NBPts = 2*(w+w)+1;
   size_t Larg = Image->cols(), Haut = Image->rows();
   typedef typename ImageType::value PixelType;

   PixelType * TabGauche = new PixelType[w+1+p.LargeurFourier];
   PixelType * TabDroite = new PixelType[w+1+p.LargeurFourier];   
   
   IPoint2D * TabCoupe         = new IPoint2D[NBPts];
   IPoint2D * TabCoupeGauche   = new IPoint2D[w+1+p.LargeurFourier];
   for(int i=0;i<NBPts;i++)
   {
      TabCoupe[i].x=0;
      TabCoupe[i].y=0;
   }      
   //--------------------------------------------------------------------
   //			Construction object ccurvidetect1d
   //--------------------------------------------------------------------
   CCurviDetect1D<ImageType> CurviDetect1D(w+1+p.LargeurFourier, p.LargeurFourier, p.LargeurFourier);

   //--------------------------------------------------------
   //			Suivi de contours
   //--------------------------------------------------------
   CCurviDetect2D CurviDetect2D; 
   CurviDetect2D.SuiviFreeman(ImageFreeman);

   if (Callback)
      Callback(40);

   //---------------------------------------------------------
   //				Affichage couleur fenetre résultats
   //---------------------------------------------------------
   size_t NbListe = CurviDetect2D.TabListePoint.size();
   size_t TotalProfiles = 0;
   double gt1=0, gt2=0, gt3=0, gt4=0, gt5=0, gt6=0;
   int CptExtrapolGauche, CptExtrapolDroite;

   for(size_t i=0;i<NbListe;i++)
   {
      CptExtrapolGauche = 0;
      CptExtrapolDroite = 0;
      if (Callback)
         Callback(40+50*i/NbListe);
      if(int(CurviDetect2D.TabListePoint[i]->size())>p.LongueurMiniForme)
      {
         double t1, t2, t3, t4, t5, t6;       
         CurvilinearRegion *ListePointsFiltres = new CurvilinearRegion;// creation pour les listes dépassants la contraintes de largeurs mini.
         CurvilinearRegion *ListePointsFiltresComp = new CurvilinearRegion;

//         ListePointsFiltres->reserve(CurviDetect2D.TabListePoint[i]->size());
//         ListePointsFiltresComp->reserve(CurviDetect2D.TabListePoint[i]->size());
         
         PointList::iterator pos = CurviDetect2D.TabListePoint[i]->begin();

         double LastWidthGauche = -1.0,
                LastWidthDroite = -1.0;

         while (pos != CurviDetect2D.TabListePoint[i]->end())
         {
            t1=get_time();
            int SigneDx = 0;
            
            IPoint2D CurentPoint = *pos;
            float GradMoy_x = PetrouContour.I_x->el[CurentPoint.y][CurentPoint.x];
            float GradMoy_y = PetrouContour.I_y->el[CurentPoint.y][CurentPoint.x];
            
            if(GradMoy_x>0)
               SigneDx = 1; 
            else
               SigneDx = -1;

            if(GradMoy_x==0 && GradMoy_y>=0.0) SigneDx=1;
            if(GradMoy_x==0 && GradMoy_y<0.0) SigneDx=-1;
           
            int Centre;
            CurviDetect2D.RemplitTableauIJCoupeBiDir(TabCoupe,
                                                     CurentPoint.x,CurentPoint.y,
                                                     GradMoy_x,GradMoy_y,
                                                     w,0,SigneDx,&Centre,Haut,Larg);
            t2=get_time();
            TotalProfiles++;
            //----------------------------------------------------------------------
            //			Creation tableau gauche et droite
            //----------------------------------------------------------------------
				for(int j=Centre-w;j<Centre+1+p.LargeurFourier;j++)//Gauche
				{
					TabGauche[Centre-j+p.LargeurFourier] = Image->el[TabCoupe[j].y][TabCoupe[j].x];
					TabCoupeGauche[Centre+p.LargeurFourier-j] = TabCoupe[j];
				}
   			
				for(int j=Centre-p.LargeurFourier, cpti=0; j<Centre+w+1; j++)//droite
				{
					TabDroite[cpti++] = Image->el[TabCoupe[j].y][TabCoupe[j].x];			
				}//fin for j

            t3=get_time();
 
            if (p.DoFFTOnEdgesOnly)
            CurviDetect1D.RespCurviMaxiFast(TabGauche, TabCoupeGauche, LastWidthGauche, ImageFreeman);
            else 
               CurviDetect1D.RespCurviMaxiFast(TabGauche, TabCoupeGauche, LastWidthGauche, 0);
          
            PosDetectCurvi PosGauche;
				PosGauche.PosInit = TabCoupe[Centre+p.LargeurFourier-CurviDetect1D.MaxiCurvi.IptInit]; 
				PosGauche.PosFinal = TabCoupe[Centre+p.LargeurFourier-CurviDetect1D.MaxiCurvi.IptFinal];
				PosGauche.PositionMaxi = TabCoupe[Centre+p.LargeurFourier-CurviDetect1D.MaxiCurvi.Centre];
            PosGauche.ValDist = CurviDetect1D.MaxiCurvi.ValDist;
            PosGauche.Largeur = LastWidthGauche = CurviDetect1D.MaxiCurvi.Largeur;
				PosGauche.PosAxe.x = ((double)PosGauche.PosInit.x+(double)PosGauche.PosFinal.x)/2;
				PosGauche.PosAxe.y = ((double)PosGauche.PosInit.y+(double)PosGauche.PosFinal.y)/2;
            t4=get_time();


            if ((PosGauche.ValDist >= p.SeuilReponse) && (PosGauche.Largeur >= p.LargeurMini))
				{
					ListePointsFiltres->push_back(PosGauche);// ajout à la liste des candidats gardés
					CptExtrapolGauche=0;
				}
				else
				{
					CptExtrapolGauche++;
					if (CptExtrapolGauche > p.SeuilExtrapol)
					{
						if (int(ListePointsFiltres->size()) > p.LongueurMiniForme) 
					      TableauListesPointsFiltrees.push_back(ListePointsFiltres);
						else 
                     delete ListePointsFiltres;
                  ListePointsFiltres = new CurvilinearRegion;// creation pour les listes dépassants la contraintes de largeurs mini.
						CptExtrapolGauche=0;
					}
				}

				PosDetectCurvi PosDroite;
            if (p.DoFFTOnEdgesOnly)
				CurviDetect1D.RespCurviMaxiFast(TabDroite, TabCoupe+Centre+1-p.LargeurFourier, LastWidthDroite, ImageFreeman);
            else
   				CurviDetect1D.RespCurviMaxiFast(TabDroite, TabCoupe+Centre+1-p.LargeurFourier, LastWidthDroite, 0);

				PosDroite.PosInit = TabCoupe[Centre+1-p.LargeurFourier+CurviDetect1D.MaxiCurvi.IptInit];
				PosDroite.PosFinal = TabCoupe[Centre+1-p.LargeurFourier+CurviDetect1D.MaxiCurvi.IptFinal];
				PosDroite.PositionMaxi = TabCoupe[Centre+1-p.LargeurFourier+CurviDetect1D.MaxiCurvi.Centre];
				PosDroite.ValDist = CurviDetect1D.MaxiCurvi.ValDist;
				PosDroite.Largeur = LastWidthDroite = CurviDetect1D.MaxiCurvi.Largeur;
				PosDroite.PosAxe.x=((double)PosDroite.PosInit.x+(double)PosDroite.PosFinal.x)/2;
				PosDroite.PosAxe.y=((double)PosDroite.PosInit.y+(double)PosDroite.PosFinal.y)/2;
            t5=get_time();

            if ((PosDroite.ValDist >= p.SeuilReponse) && (PosDroite.Largeur >= p.LargeurMini))
            {               
               ListePointsFiltresComp->push_back(PosDroite);// ajout à la liste des candidats gardés
					CptExtrapolDroite=0;	
            } else {
					CptExtrapolDroite++;
					if (CptExtrapolDroite > p.SeuilExtrapol)
					{
                  if (int(ListePointsFiltresComp->size()) > p.LongueurMiniForme)
							TableauListesPointsFiltrees.push_back(ListePointsFiltresComp);
                  else 
                     delete ListePointsFiltresComp;
                  ListePointsFiltresComp = new CurvilinearRegion;
						CptExtrapolDroite=0;
					}
				}

            t6=get_time(); pos++; gt1 += t2-t1; gt2 += t3-t2; gt3 += t4-t3; gt4 += t5-t4; gt5 += t6-t5; gt6 += t6-t1;
         }
         if (int(ListePointsFiltres->size())>p.LongueurMiniForme)
            TableauListesPointsFiltrees.push_back(ListePointsFiltres);
         if (int(ListePointsFiltresComp->size())>p.LongueurMiniForme)
            TableauListesPointsFiltrees.push_back(ListePointsFiltresComp);	
      }//fin if nb points liste > 10
      //ajout de la liste en cours au tableau de liste de structure poscurvidetect
   }//for i

   if (Callback)
      Callback(90);

/*
   char buf[100];
   sprintf(buf, "G1: %.5fs, G2: %.5fs, G3: %.5fs, G4: %.5fs, G5: %.5fs FFT: %.5fs\n", gt1,gt2,gt3,gt4,gt5, FFT_time);
   OutputDebugString(buf);
   */
   


   if (p.Do2dConstraints)
      TwoDimConstraints(TableauListesPointsFiltreesC2D, TableauListesPointsFiltrees, p);
   else 
      TableauListesPointsFiltreesC2D = TableauListesPointsFiltrees;

   SectionOptimization(TableauListesPointsFiltreesC2D, Larg, Haut);

   if (p.DoFilterCompact)
      FilterCirclesAndSquares(TableauListesPointsFiltreesC2D, p);
   
   if (p.DoMergeRegions)
   {
      FindOverlappingRegions(TableauListesPointsFiltreesC2D, p);

      // mark all the regions as OK
      for (TableOfCurvilinearRegions::iterator it = TableauListesPointsFiltreesC2D.begin(); it != TableauListesPointsFiltreesC2D.end(); it++)
         (*it)->merged = false;

      // fast merge
      MergeConnectedRegions(TableauListesPointsFiltreesC2D, Larg, Haut);
   }

   // remove (merge) duplicates...
   if (p.DoRemoveDuplicates)
      FilterOverlappingAndShortRegions(TableauListesPointsFiltreesC2D, p, Larg, Haut);

   

   // filter merged regions
   for (TableOfCurvilinearRegions::iterator it = TableauListesPointsFiltreesC2D.begin(); it != TableauListesPointsFiltreesC2D.end(); )
   {
      if ((*it)->merged)
         // delete it...
         it = TableauListesPointsFiltreesC2D.erase(it);
      else 
         it++;
   }

   if (Callback)
      Callback(100);

   delete [] TabCoupe;
   delete [] TabCoupeGauche;
   delete [] TabDroite;
   delete [] TabGauche;
   delete ImageFreeman;
   delete ImageGauss;
/*
   char buf[100];
   sprintf(buf, "Processed lists %d, profiles %d\n", NbListe, TotalProfiles);
   OutputDebugString(buf);*/
}

#endif // __CURVLINEAR_H__
