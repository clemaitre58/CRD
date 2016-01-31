// CurviDetect1D.h: interface for the CCurviDetect1D class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CURVIDETECT1D_H__F29E1257_B2C6_46E2_A051_2ABDD0C76E57__INCLUDED_)
#define AFX_CURVIDETECT1D_H__F29E1257_B2C6_46E2_A051_2ABDD0C76E57__INCLUDED_

#include <ary.h>
#include <dtypes.h>
#include "struct.h"

extern const double gauss_table_sigma_1[];


void InitFenetreHaming(double *Fenetre, int LargeurFiltre);

template<typename SignalType>
void CorrectionHaming(SignalType* Signal, double* FenetreHaming, size_t Largeur)
{
	for(size_t i=0;i<Largeur;i++)
		Signal[i]=Signal[i]*FenetreHaming[i];
}

template <typename SignalType>
bool FFT(int dir, int m, SignalType *x, SignalType *y)
{
//   double time1 = get_time();
	long nn,i,i1,j,k,i2,l,l1,l2;
	SignalType c1,c2,tx,ty,t1,t2,u1,u2,z;

	/* Calculate the number of points */
	nn = 1<<m;

	/* Do the bit reversal */
	i2 = nn >> 1;
	j = 0;
	for (i=0;i<nn-1;i++) {
		if (i < j) {
			tx = x[i];
			ty = y[i];
			x[i] = x[j];
			y[i] = y[j];
			x[j] = tx;
			y[j] = ty;
		}
		k = i2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}

	/* Compute the FFT */
	c1 = -1.0;
	c2 = 0.0;
	l2 = 1;
	for (l=0;l<m;l++) {
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0;
		u2 = 0.0;
		for (j=0;j<l1;j++) {
			for (i=j;i<nn;i+=l2) {
				i1 = i + l1;
				t1 = u1 * x[i1] - u2 * y[i1];
				t2 = u1 * y[i1] + u2 * x[i1];
				x[i1] = x[i] - t1;
				y[i1] = y[i] - t2;
				x[i] += t1;
				y[i] += t2;
			}
			z =  u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		c2 = sqrt((SignalType(1.0) - c1) / 2.0);
		if (dir == 1)
			c2 = -c2;
		c1 = sqrt((SignalType(1.0) + c1) / 2.0);
	}
   return true;
}

double HyperDistanceModuleFourier(utls::RGBValue<double> * aR, utls::RGBValue<double> * aI, 
                                  utls::RGBValue<double> * bR, utls::RGBValue<double> * bI, int LargFiltre);
double HyperDistanceModuleFourier(double * aR, double * aI, double * bR, double * bI, int LargFiltre);

template <typename AryBase>
class CCurviDetect1D
{
   typedef typename AryBase::value                                        ImageSignalType;
   typedef typename AccumulatorDeduction<typename AryBase::value>         AccumulatorType;
   typedef typename AccumulatorType::value                                SignalType;
private:
   int LongueurSignal;
   int LargeurFourier;
   int IptInitFixe;
   int IptFinalMax;

   double ValDist;

   int RacineLargeurFourier;
   SignalType * FGaucheRE;
   SignalType * FDroiteRE;
   SignalType * FCentreGRE;
   SignalType * FCentreDRE;

   SignalType * FGaucheIM;
   SignalType * FDroiteIM;
   SignalType * FCentreGIM;
   SignalType * FCentreDIM;

   int TailleTabDist;

   double * TabDistanceDroite;
   double * TabDistanceGauche;
   double * TabBorSim;

   // Fenetre de correction de Haming
   double * Fenetre;

public:

   //////////////////////////////////////////////////////////////////////
   // Construction/Destruction
   //////////////////////////////////////////////////////////////////////
   CCurviDetect1D::CCurviDetect1D(int ValLongueurSignal, int ValLargeurFourier, int ValIptInitFixe)
   {
      LongueurSignal = ValLongueurSignal;
      LargeurFourier = ValLargeurFourier;
      IptInitFixe = ValIptInitFixe;
      IptFinalMax = LongueurSignal-LargeurFourier-1; //attention 

      RacineLargeurFourier = 0;
      while ((1<<RacineLargeurFourier)<LargeurFourier) RacineLargeurFourier++;

      FGaucheRE = new SignalType[LargeurFourier];
      FDroiteRE = new SignalType[LargeurFourier];
      FCentreGRE = new SignalType[LargeurFourier];
      FCentreDRE = new SignalType[LargeurFourier];

      FGaucheIM = new SignalType[LargeurFourier];
      FDroiteIM = new SignalType[LargeurFourier];
      FCentreGIM = new SignalType[LargeurFourier];
      FCentreDIM = new SignalType[LargeurFourier];

      TailleTabDist = (IptFinalMax-IptInitFixe); // 13_11_2007 pas sur

      TabDistanceDroite = new double [TailleTabDist];
      TabDistanceGauche = new double [TailleTabDist];
      TabBorSim = new double [TailleTabDist];

      // Allocation du tableau qui contiendra la fenetre de hamming
      Fenetre = new double[LargeurFourier];

      // Initialisation de la Fenetre de Haming
      InitFenetreHaming(Fenetre, LargeurFourier);
   }

   virtual CCurviDetect1D::~CCurviDetect1D()
   {
      delete [] FGaucheRE; 
      delete [] FDroiteRE;
      delete [] FCentreGRE;
      delete [] FCentreDRE;

      delete [] FGaucheIM; 
      delete [] FDroiteIM; 
      delete [] FCentreGIM; 
      delete [] FCentreDIM;

      delete [] TabDistanceDroite;
      delete [] TabDistanceGauche;
      delete [] TabBorSim;
      delete [] Fenetre;
   }
   
   RespCurvi MaxiCurvi;
   RespCurvi SecondePasse;

   int * TabLargeur;

   void RespCurviMaxiFast(ImageSignalType * TabCoupe, IPoint2D *TabCoupeXY, double WidthPrior, utls::BAry *EdgeImage)
   {
      //----------------------------------------------------
      // recopie des tab pour fourier + calcul distance
      //----------------------------------------------------
      int i, cptTab=0, cpt=0;
      for (i = IptInitFixe-LargeurFourier;i<IptInitFixe;i++)
      {
         FGaucheRE[cpt] = TabCoupe[i];
         FGaucheIM[cpt]= 0;
         cpt++;
      }
      cpt =0;
      for (i = IptInitFixe+1;i<IptInitFixe+1+LargeurFourier;i++)//modif 09/01/2007
      {
         FCentreGRE[cpt] = TabCoupe[i];
         FCentreGIM[cpt]= 0;
         cpt++;
      }
      CorrectionHaming(FGaucheRE,Fenetre,LargeurFourier);  CorrectionHaming(FGaucheIM,Fenetre,LargeurFourier);
      FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM);

      CorrectionHaming(FCentreGRE,Fenetre,LargeurFourier); CorrectionHaming(FCentreGIM,Fenetre,LargeurFourier);
      FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);

      double MaxDroite = 0.0, 
         MaxBordSim = 0.0,

         MaxGauche = HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,LargeurFourier);
		 //sans fourier
		 //MaxGauche = HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,1);
		int Overlapp = LargeurFourier/4;
      for (int IptFinal=(IptInitFixe+Overlapp+1); IptFinal<IptFinalMax; IptFinal++) //idem modif 21/11/2006
      {
         if (EdgeImage == 0 || EdgeImage->el[TabCoupeXY[IptFinal].y][TabCoupeXY[IptFinal].x] != 0)
         {
            cpt =0;
            for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)
            {
               FDroiteRE[cpt] = TabCoupe[i]; 
               FDroiteIM[cpt]= 0;
               cpt++;
            }

            cpt =0;
            for (i = IptFinal-LargeurFourier;i<IptFinal;i++)//modif 09/01/2007
            {
               FCentreDRE[cpt] = TabCoupe[i];
               FCentreDIM[cpt]= 0;
               cpt++;
            }

            //Correction de Haming
            CorrectionHaming(FDroiteRE,Fenetre,LargeurFourier); CorrectionHaming(FDroiteIM,Fenetre,LargeurFourier);
            //calcul FFT
            FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);

            //Correction de Haming
            CorrectionHaming(FCentreDRE,Fenetre,LargeurFourier); CorrectionHaming(FCentreDIM,Fenetre,LargeurFourier);
            //calcul FFT
            FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);

            TabDistanceDroite[cptTab]  = HyperDistanceModuleFourier(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,LargeurFourier);
            TabBorSim[cptTab]          = HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,LargeurFourier);
			//sans fourier
			/*TabDistanceDroite[cptTab]  = HyperDistanceModuleFourier(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,1);
            TabBorSim[cptTab]          = HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,1);*/
         }
         else
         {
            TabDistanceDroite[cptTab]=0.0;
            TabBorSim[cptTab]=0.0;
         }

         cptTab ++;
      }

      //-----------------------------------------------------------------------------------------------------
      // recherche des max pour l'inversion des fonctions de probalilité "section homogène" et bord similaire
      //-----------------------------------------------------------------------------------------------------

      for(i=0;i<cptTab;i++)
      {
         if (TabBorSim[i]>MaxBordSim)        MaxBordSim = TabBorSim[i];		
         if (TabDistanceDroite[i]>MaxDroite) MaxDroite  = TabDistanceDroite[i];      
      }

      double MaxNorma=utls::max(MaxGauche, utls::max(MaxBordSim, MaxDroite));   
      if (MaxNorma==0) MaxNorma=1.0;

      //--------------------------------------------------------
      //			Fonction export + calcul des couts
      //--------------------------------------------------------   
      double ValDistMax=-1.0;

      MaxiCurvi.ValDist = 0.0f;
      MaxiCurvi.IptInit=0;
      MaxiCurvi.Centre=0;
      MaxiCurvi.IptFinal=0;
      MaxiCurvi.Largeur=0;

      // find the maximum over the valid signal
      for (int IptFinal = IptInitFixe+Overlapp+1,i=0; IptFinal<IptFinalMax; IptFinal++,i++) // modif 21/11/2006 rajout pour limiter les surbalayage
      {
         double width = utls::l2norm2D(
            TabCoupeXY[IptInitFixe].x - TabCoupeXY[IptFinal].x, 
            TabCoupeXY[IptInitFixe].y - TabCoupeXY[IptFinal].y);
         int g;
         if (WidthPrior>0)
         {
            // g = int(((min(width/WidthPrior, WidthPrior/width)-1.0)/0.1))+45; 
            g = int( ( (width - WidthPrior)/91.0)*10.0 )+45;  // 45
            g = utls::max(0, utls::min(g, 90));
         } else 
            g = 45;

         // Sans la section homogene
         ValDist = (MaxGauche/MaxNorma)* 
            (TabDistanceDroite[i]/MaxNorma)*
            ((MaxNorma-TabBorSim[i])/MaxNorma)*
            gauss_table_sigma_1[g];       
         
         if (ValDist>ValDistMax && (EdgeImage == 0 || EdgeImage->el[TabCoupeXY[IptFinal].y][TabCoupeXY[IptFinal].x] != 0))
         {
            ValDistMax = ValDist;
            MaxiCurvi.ValDist = ValDist;
            MaxiCurvi.Largeur = width; // IptFinal-IptInitFixe;
            MaxiCurvi.IptInit = IptInitFixe;
            MaxiCurvi.IptFinal = IptFinal;
            MaxiCurvi.Centre = (IptInitFixe+IptFinal)/2;
         }
      }//fin IptFinal
      assert(i == cptTab);

#ifdef DO_SECOND_PASS
      //-------------------------------------------------------------
      //2° passe
      //-------------------------------------------------------------
      if(MaxiCurvi.IptFinal-IptInitFixe > LargeurFourier+1)
      {
         cptTab=0; cpt = 0;
         for (i = MaxiCurvi.IptFinal+1;i<MaxiCurvi.IptFinal+LargeurFourier+1;i++)
         {
            FDroiteRE[cpt] = TabCoupe[i];
            FDroiteIM[cpt] = 0.0;
            cpt++;
         }

         cpt = 0;
         for (i = MaxiCurvi.IptFinal-LargeurFourier;i<MaxiCurvi.IptFinal;i++)//modif 09/01/2007
         {
            FCentreDRE[cpt] = TabCoupe[i];
            FCentreDIM[cpt] = 0.0;
            cpt++;
         }

         CorrectionHaming(FDroiteRE,Fenetre,LargeurFourier);  CorrectionHaming(FDroiteIM,Fenetre,LargeurFourier);
         FFT(1,RacineLargeurFourier,FDroiteRE,FDroiteIM);

         CorrectionHaming(FCentreDRE,Fenetre,LargeurFourier); CorrectionHaming(FCentreDIM,Fenetre,LargeurFourier);
         FFT(1,RacineLargeurFourier,FCentreDRE,FCentreDIM);

         double MaxDroite = HyperDistanceModuleFourier(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,LargeurFourier),
		 //sans fourier avec ValMoy
	//	    double MaxDroite = HyperDistanceModuleFourier(FDroiteRE,FDroiteIM,FCentreDRE,FCentreDIM,1), 
            MaxBordSim = 0.0,
            MaxGauche = 0.0;

         //----------------------------------------------------
         // recopie des tab pour fourier + calcul distance
         //----------------------------------------------------
         for (int IptFinal=IptInitFixe;IptFinal<MaxiCurvi.IptFinal-LargeurFourier-1;IptFinal++) //idem modif 21/11/2006
         {
            cpt =0;
            for (i = IptFinal-LargeurFourier;i<IptFinal;i++)
            {
               FGaucheRE[cpt] = TabCoupe[i];
               FGaucheIM[cpt] = 0;
               cpt++;
            }

            cpt =0;
            for (i = IptFinal+1;i<IptFinal+LargeurFourier+1;i++)//modif 09/01/2007
            {
               FCentreGRE[cpt] = TabCoupe[i];
               FCentreGIM[cpt] = 0;
               cpt++;
            }

            // Correction de Haming and calcul FFT
            CorrectionHaming(FGaucheRE,Fenetre,LargeurFourier); CorrectionHaming(FGaucheIM,Fenetre,LargeurFourier);
           FFT(1,RacineLargeurFourier,FGaucheRE,FGaucheIM); 

           CorrectionHaming(FCentreGRE,Fenetre,LargeurFourier); CorrectionHaming(FCentreGIM,Fenetre,LargeurFourier);
           FFT(1,RacineLargeurFourier,FCentreGRE,FCentreGIM);

            TabBorSim[cptTab]         = HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,LargeurFourier);
            TabDistanceGauche[cptTab] = HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,LargeurFourier);
			// sans fourier
			/*TabBorSim[cptTab]         = HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FDroiteRE,FDroiteIM,1);
            TabDistanceGauche[cptTab] = HyperDistanceModuleFourier(FGaucheRE,FGaucheIM,FCentreGRE,FCentreGIM,1);*/

            cptTab++;
         }

         //-----------------------------------------------------------------------------------------------------
         // recherche des max pour l'inversion des fonctions de probalilité "section homogène" et bord similaire
         //-----------------------------------------------------------------------------------------------------

         for(i=0;i<cptTab;i++)
         {
            if(TabBorSim[i]>MaxBordSim)        MaxBordSim = TabBorSim[i];
            if(TabDistanceGauche[i]>MaxGauche) MaxGauche  = TabDistanceGauche[i];
         }

         MaxNorma=__max(MaxGauche, __max(MaxBordSim, MaxDroite));   
         if (MaxNorma==0) MaxNorma=1.0;

         //--------------------------------------------------------
         //			Fonction export + calcul des couts
         //--------------------------------------------------------

         i=0; // compteur
         ValDistMax = -1.0;
         SecondePasse.ValDist = 0.0f;
         SecondePasse.IptInit = 0;
         SecondePasse.IptFinal = 0;
         SecondePasse.Largeur = 0;
         SecondePasse.Centre=0;
         for(int IptFinal=IptInitFixe;IptFinal<MaxiCurvi.IptFinal-LargeurFourier-1;IptFinal++) // modif 21/11/2006 rajout pour limiter les surbalayage
         {
            double width = utls::l2norm2D(
               TabCoupeXY[IptFinal].x - TabCoupeXY[MaxiCurvi.IptFinal].x, 
               TabCoupeXY[IptFinal].y - TabCoupeXY[MaxiCurvi.IptFinal].y);         
            int g;

            if (WidthPrior>0)
            {
               // g = int(((min(width/WidthPrior, WidthPrior/width)-1.0)/0.1))+45; 
               g = int( ( (width - WidthPrior)/91.0)*10.0 )+45;  // 45
               g = __max(0, __min(g, 90));
            } else
               g = 45;

            ValDist = (TabDistanceGauche[i]/MaxNorma) * 
               (MaxDroite/MaxNorma)*
               ((MaxNorma-TabBorSim[i])/MaxNorma)*
               gauss_table_sigma_1[g];
            if (ValDistMax<ValDist)
            {
               SecondePasse.ValDist = ValDist;
               SecondePasse.IptInit = IptFinal;
               SecondePasse.IptFinal = MaxiCurvi.IptFinal;
               SecondePasse.Largeur = width;
               SecondePasse.Centre = (SecondePasse.IptFinal+SecondePasse.IptInit)/2;
               ValDistMax=ValDist;
            }
            i++;
         }//fin IptFinal
         MaxiCurvi = SecondePasse;
      }
#endif
   }

   //double MaxGauche2Passe;
   //double MaxDroite2Passe;
   //double MaxHomoSec2Passe;
   //double MaxBordSim2Passe;
   //double TempValDis2Passe;

};

#endif // !defined(AFX_CURVIDETECT1D_H__F29E1257_B2C6_46E2_A051_2ABDD0C76E57__INCLUDED_)
