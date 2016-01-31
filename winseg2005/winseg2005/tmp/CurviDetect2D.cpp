#include "stdafx.h"
#include "conio.h"
#include "STRUCT.H"
#include "Petrou.h"
#include "CurviDetect2D.h"


CCurviDetect2D::CCurviDetect2D(/*void*/)
{
}

CCurviDetect2D::~CCurviDetect2D(void)
{
	
	int i=0;
	if(TabListePoint.GetCount()!=0)
	{
		for(i=0;i<NbListeInArray;i++)
		{
			TabListePoint[i]->RemoveAll();
			delete TabListePoint[i];
		}
	}
}

BOOL CCurviDetect2D::SuiviContour(BYTE *  pImage, BYTE* pImageDest, long Largeur, long Hauteur)
{
long LargIm = Largeur; 
long HautIm = Hauteur;


BOOL PasContour=TRUE;
int essai=0;
Taille  = Largeur*Hauteur;

BordAZero(pImage,Largeur,Hauteur);	// par sécurité

do
{
	BOOL PointDepart = false;
	long i = 0,j = HautIm-1;

	nbelement = 0;

	ListePointContour->RemoveAll();


	do
	{
		if (pImage[(j*LargIm)+i] !=0) 
		{
			FirstPoint.x = i;
			FirstPoint.y = j;
			PointDepart = TRUE;
		}		
		i++;
		if (i == LargIm)
		{
			i=0;
			j--;
		}

	}while(PointDepart==FALSE && j<HautIm);

	if (PointDepart==FALSE) return FALSE; // si image toute noire, on quitte la fonction

	CurentPoint = FirstPoint;
	int r,choix,spp=0;
	r = 1;
	do
	{
		choix=spp+r;
		switch(choix)
		{
			case 0 : CurentPoint.y++; spp=6; break;
			case 1 : CurentPoint.y--; spp=2; break;
			case 2 : CurentPoint.x++; spp=0; break;
			case 3 : CurentPoint.x--; spp=4; break;
			case 4 : CurentPoint.y--; spp=2; break;
			case 5 : CurentPoint.y++; spp=6; break;
			case 6 : CurentPoint.x--; spp=4; break;
			case 7 : CurentPoint.x++; spp=0; break;
		}

		if(pImage[CurentPoint.x + LargIm*CurentPoint.y] !=0)
		{
			ListePointContour->AddTail(CurentPoint);
			r=1;
			nbelement++; // <- Ajout de code ALEX
		}
		else
		{
			r=0;
		}

	}while (FirstPoint != CurentPoint); 

	if (nbelement>1) 
		PasContour=FALSE;
	else
	{
		pImage[(CurentPoint.y*LargIm)+CurentPoint.x] =0;
		essai++;
		if (essai>100) PasContour=FALSE; 
	}

} while (PasContour);


if (!PasContour)
{
	if (m_TabDistance!=NULL) free(m_TabDistance);
	m_TabDistance=(double *) malloc(nbelement*sizeof(double));

	//////////////////////////////////////////////////////////////////////////////////////////////
	//    pour matérialiser le contour 
	//////////////////////////////////////////////////////////////////////////////////////////////

	//POSITION pos;
	pos = ListePointContour->GetHeadPosition();
	long p=0;
	Cmin=Largeur;
	Cmax=0;
	Lmin=Hauteur;
	Lmax=0;
	while(pos)
	{
		CurentPoint = ListePointContour->GetAt(pos);
		p=CurentPoint.x +LargIm*CurentPoint.y;
		if(Cmin>CurentPoint.x) Cmin=CurentPoint.x;
		if(Cmax<CurentPoint.x) Cmax=CurentPoint.x;
		if(Lmin>CurentPoint.y) Lmin=CurentPoint.y;
		if(Lmax<CurentPoint.y) Lmax=CurentPoint.y;
		if (p<0 || p>Taille) p=0;
		pImageDest[p] = 127; 
		ListePointContour->GetNext(pos);
	}
}
	return TRUE;
}

void CCurviDetect2D::BordAZero(BYTE* pImage, long Largeur, long Hauteur)
{
	// --- périmètre de sécurité --------------------------------------------- //

	int k;
	long Taille = Largeur*Hauteur;

	for (k=0;k<Largeur;k++) pImage[k]=0;
	for (k=Taille-Largeur+1;k<Taille;k++) pImage[k]=0;
	for (k=0;k<Hauteur;k++) pImage[k*Largeur]=0;
	/*for (k=0;k<Hauteur;k++) pImage[k*Largeur+1]=0;
	for (k=0;k<Hauteur;k++) pImage[k*Largeur+2]=0;
	for (k=0;k<Hauteur;k++) pImage[k*Largeur+3]=0;*/
	for (k=0;k<Hauteur;k++) pImage[k*Largeur+(Largeur-1)]=0;
}
void CCurviDetect2D::FermetureContourSimple(BYTE * pImage, long Largeur, long Hauteur)
{
	int l,c;
	int ll,cc;
	int p,pp,p1,p2,p3,p4,p5,p6,p7,p8,p9;
	int Somme;
	for (l=1;l<Hauteur-1;l++)
	{
		for (c=1;c<Largeur-1;c++)
		{
			//Somme=0;
			p=Largeur*l+c;
			p5=p;
			p6=p+1;
			p4=p-1;
			p8=p5+Largeur;
			p7=p8-1;
			p9=p8+1;
			p2=p5-Largeur;
			p1=p2-1;
			p3=p2+1;
			/*for (ll=-1;ll<2;ll++)
			{
				for (cc=-1;cc<2;cc++)
				{
					pp=Largeur*(l+ll)+c+cc;
					if (pImage[pp]) Somme++;
				}
			}*/
			if (/*Somme==2 &&*/ pImage[p1]==255 && pImage[p9]==255)
			{
				pImage[p]=255;
			}
			if (/*Somme==2 &&*/  pImage[p3]==255 && pImage[p7]==255)
			pImage[p]=255;

			if (/*Somme==2 &&*/  pImage[p2]==255 && pImage[p8]==255)
				pImage[p]=255;
			if (/*Somme==2 &&*/  pImage[p4]==255 && pImage[p6]==255)
				pImage[p]=255;
		}
	}
}

BOOL CCurviDetect2D::SuiviFreeman(BYTE * pImage, long Largeur, long Hauteur,CMDIFrameWnd *pmd)
{

long LargIm = Largeur; 
long HautIm = Hauteur;
char CodeFreeman=0;
ConnexityDegreeAFreemanCode ConnexityFisrtPoint;
CPoint NewPoint;
ConnexityFisrtPoint.ConnexityDegree=0;
//POSITION AvantDernierPoint;
//POSITION PositionTete;
POSITION PositionQueue;
CPoint Test;
CPoint PointRetour;
BOOL CasPremier = FALSE;

ParcoursFreemanHoraire FreemanHoraire;
ParcoursFreemanAntiHoraire FreemanAntiHoraire;
DirectPrecedent FreemanCasPrecedent;

char CodeFreemanPrevious = 3;
char CodeFreemanPreviousAnti = 7;
long CptPasBlanc = 0;

int * CopieImageInt = (int*) malloc(sizeof(int)*LargIm*HautIm);



BOOL PasContour=TRUE;
BOOL ConditionArret = TRUE;
BOOL Rebouclage = FALSE;
BOOL BoolRetour = FALSE;
BOOL PointDepart;
BOOL CandidatTrouve = FALSE;
bool bFinFreeman=false;
int essai=0;
int NbZero;
int NbMoinsUn =0;
long i = 0,j = 0;


//CImageWnd * pTemp=NULL;	
//pTemp=(CTraitImageWnd * ) new CImageWnd(LargIm,HautIm, COLOR);
//pTemp->Create("Debug color",WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,pmd);

for(long z=0;z<LargIm*HautIm;z++)
{
//pTemp->ImageRGB[z].r = pImage[z];
//pTemp->ImageRGB[z].g = pImage[z];
//pTemp->ImageRGB[z].b = pImage[z];
CopieImageInt[z] = (int)pImage[z];
}

//MessageBox(NULL,"STOP","fin copie",MB_OK);

Taille  = Largeur*Hauteur;

BordAZero(pImage,Largeur,Hauteur);	// par sécurité
//MessageBox(NULL,"STOP","Bord a zero",MB_OK);
//TRACE("Bord a zero ok");
NbListeInArray=0;
//init tableau
//for (int i=0;i<8;i++)
//{
//FreemanHoraire.Code[i] = i;
//}

FreemanHoraire.Code[0]= 0;
FreemanHoraire.Code[1]= 2;
FreemanHoraire.Code[2]= 4;
FreemanHoraire.Code[3]= 6;
FreemanHoraire.Code[4]= 1;
FreemanHoraire.Code[5]= 3;
FreemanHoraire.Code[6]= 5;
FreemanHoraire.Code[7]= 7;

FreemanCasPrecedent.Code[0]= 4;
FreemanCasPrecedent.Code[1]= 5;
FreemanCasPrecedent.Code[2]= 6;
FreemanCasPrecedent.Code[3]= 7;
FreemanCasPrecedent.Code[4]= 0;
FreemanCasPrecedent.Code[5]= 1;
FreemanCasPrecedent.Code[6]= 2;
FreemanCasPrecedent.Code[7]= 3;

/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////recherche du premier point/////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//TRACE("Init Freeman ok");

//MessageBox(NULL,"STOP","Init Freeman",MB_OK);

FirstPoint.x=0;
FirstPoint.y=0;
//do
//{

do
{
	PointDepart = false;
	i=FirstPoint.x;
	j=FirstPoint.y;

	nbelement = 0;

	do
	{
		if (CopieImageInt[(j*LargIm)+i] == 255) 
		{
			FirstPoint.x = i;
			FirstPoint.y = j;
			PointDepart = TRUE;
			//CopieImageInt[(j*LargIm)+i] = -1;
		}		
		i++;
		if (i >= LargIm)
		{
			i=0;
			j++;
		}

	}while(PointDepart==FALSE && j<HautIm);

/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////Etude de la connexité du premier point/////////////////////
/////////////////////////////////////////////////////////////////////////////////////

if(PointDepart)
{
//allocation d'un pointeur pour ajouter une nouvelle liste à l'array


CurentPoint = FirstPoint;

//pTemp->ImageRGB[FirstPoint.x + LargIm*FirstPoint.y].r = 255;
//pTemp->ImageRGB[FirstPoint.x + LargIm*FirstPoint.y].g = 255;
//pTemp->ImageRGB[FirstPoint.x + LargIm*FirstPoint.y].b = 0;
//pTemp->Affichage();
//Sleep(10);

ConnexityFisrtPoint.ConnexityDegree = 0;
CodeFreeman=0;

for(i=0;i<8;i++)
{
	NewPoint = CurentPoint;
	switch(CodeFreeman)
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
	if(CopieImageInt[(NewPoint.y*LargIm)+NewPoint.x]==255)
	{
		ConnexityFisrtPoint.ConnexityDegree++;
		ConnexityFisrtPoint.FreemanCode[CodeFreeman]=TRUE;
	}
	else
	{
		ConnexityFisrtPoint.FreemanCode[CodeFreeman]=FALSE;
	}
	CodeFreeman++;
} // for i - les 8 directions

//CString ChaineTemp;
//ChaineTemp.Format("%d",ConnexityFisrtPoint.ConnexityDegree);
//
//MessageBox(NULL,"STOP : " + ChaineTemp,"Avant Parcours freeman",MB_OK);

/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////Parcours pour le suivi Freeman/////////////////////
/////////////////////////////////////////////////////////////////////////////////////
if(ConnexityFisrtPoint.ConnexityDegree!=0)
{
	ListePointContour = new CList<CPoint,CPoint>;
	ListePointContour->AddTail(FirstPoint);
	
	CurentPoint = FirstPoint;


	do
	{
		
		CodeFreeman=0;
		bFinFreeman=false;
		BoolRetour=false;
		ConditionArret=true;
		NbZero=0;
		NbMoinsUn =0;
		CandidatTrouve = FALSE;

		do	
		{
			
		

			NewPoint = CurentPoint;
			//if ((NewPoint.x == 199)&&(NewPoint.y == 76))
			//{
			//	int stop=0;
			//}
			switch(FreemanHoraire.Code[CodeFreeman])
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
			if(CopieImageInt[NewPoint.x + LargIm*NewPoint.y] ==0) 
			{
				NbZero++;
			}
			if(CopieImageInt[NewPoint.x + LargIm*NewPoint.y] ==-1) 
			{
				NbMoinsUn++;
			}
			
			

			if((CopieImageInt[NewPoint.x + LargIm*NewPoint.y] ==255) /*&& (CopieImageInt[NewPoint.x + LargIm*NewPoint.y] !=-1)*/)
			{
				
					ListePointContour->AddTail(NewPoint);
					CurentPoint = NewPoint;
					nbelement++; // <- Ajout de code ALEX

					CandidatTrouve=TRUE;
					
					/*pTemp->ImageRGB[NewPoint.x + LargIm*NewPoint.y].r = 255;
					pTemp->ImageRGB[NewPoint.x + LargIm*NewPoint.y].g = 0;
					pTemp->ImageRGB[NewPoint.x + LargIm*NewPoint.y].b = 0;*/

					CopieImageInt[NewPoint.x + LargIm*NewPoint.y] =-1;

				/*	pTemp->Affichage();
					Sleep(10);*/

					if(CurentPoint == FirstPoint)
					{
						ConditionArret = FALSE; //on sort si la tete = la queue
						Rebouclage = TRUE; // pas la peine de traiter la connexité
					}

					bFinFreeman=true; // pas la peine de continuer	
				
			}
			else// si on a pas trouvé de point blanc
			{
				if(CodeFreeman == 7)// et qu'on a balayé tous les freemans alors on sort
				{
					if(BoolRetour)
					{
						ConditionArret = FALSE;
					}
					else
					{
						if(ListePointContour->GetCount() > 1)
						{
							BoolRetour=TRUE;
							CopieImageInt[CurentPoint.x + LargIm*CurentPoint.y] =-1;
							PositionQueue = ListePointContour->GetTailPosition();
							CurentPoint = ListePointContour->GetPrev(PositionQueue);
							CurentPoint =ListePointContour->GetPrev(PositionQueue);
							CodeFreeman = -1;
							NbZero=0;
							NbMoinsUn=0;
							i=0;
						}
						else
						{	
							ConditionArret = FALSE;
						}
					} // fin else
				}
			}
			CodeFreeman++;
			if (CodeFreeman >7) bFinFreeman=true;
			
		}while(!bFinFreeman);

		if ((NbZero==7) && (CandidatTrouve==FALSE))//extremité
		{
			CopieImageInt[CurentPoint.x + LargIm*CurentPoint.y] =-1;
			ConditionArret=false;
		}

		if ((NbZero<=7)&&(NbMoinsUn<=7)&&(CandidatTrouve==FALSE))//cas d'un point emprisonné
		{
			CopieImageInt[CurentPoint.x + LargIm*CurentPoint.y] =-1;
			/*pTemp->ImageRGB[CurentPoint.x + LargIm*CurentPoint.y].r = 0;
			pTemp->ImageRGB[CurentPoint.x + LargIm*CurentPoint.y].g = 0;
			pTemp->ImageRGB[CurentPoint.x + LargIm*CurentPoint.y].b = 255;
			pTemp->Affichage();
			Sleep(10);*/

			ConditionArret=false;
		}

		if ((NbZero==8)&&(CandidatTrouve==FALSE))//cas d'un point isolé
		{
			CopieImageInt[CurentPoint.x + LargIm*CurentPoint.y] =-1;
			/*pTemp->ImageRGB[CurentPoint.x + LargIm*CurentPoint.y].r = 0;
			pTemp->ImageRGB[CurentPoint.x + LargIm*CurentPoint.y].g = 255;
			pTemp->ImageRGB[CurentPoint.x + LargIm*CurentPoint.y].b = 0;
			pTemp->Affichage();
			Sleep(10);*/

			ConditionArret=false;
		}
	
	}while (ConditionArret); //condition d'arrêt

	//MessageBox(NULL,"STOP","Après ALLER Parcours freeman",MB_OK);

	NbElementBack=0;

if ((ConnexityFisrtPoint.ConnexityDegree>1)&&(Rebouclage==FALSE))
{
	CurentPoint = FirstPoint;


	do
	{
		
		CodeFreeman=0;
		bFinFreeman=false;
		BoolRetour=false;
		ConditionArret=true;
		NbZero=0;
		NbMoinsUn =0;
		CandidatTrouve = FALSE;

		do	
		{
			
			NewPoint = CurentPoint;
			/*if((CurentPoint.x==195) && (CurentPoint.y==308))
			{
				int Stop = 0;
			}*/
			switch(FreemanHoraire.Code[CodeFreeman])
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
			if(CopieImageInt[NewPoint.x + LargIm*NewPoint.y] ==0) 
			{
				NbZero++;
			}
			if(CopieImageInt[NewPoint.x + LargIm*NewPoint.y] ==-1) 
			{
				NbMoinsUn++;
			}
			if((CopieImageInt[NewPoint.x + LargIm*NewPoint.y] ==255) /*&& (CopieImageInt[NewPoint.x + LargIm*NewPoint.y] !=-1)*/)
			{
				
					ListePointContour->AddHead(NewPoint);
					CurentPoint = NewPoint;
					nbelement++; // <- Ajout de code ALEX
					CandidatTrouve=TRUE;
					
					/*pTemp->ImageRGB[NewPoint.x + LargIm*NewPoint.y].r = 255;
					pTemp->ImageRGB[NewPoint.x + LargIm*NewPoint.y].g = 0;
					pTemp->ImageRGB[NewPoint.x + LargIm*NewPoint.y].b = 0;*/

					CopieImageInt[NewPoint.x + LargIm*NewPoint.y] =-1;

					/*pTemp->Affichage();
					Sleep(10);*/

					if(CurentPoint == FirstPoint)
					{
						CopieImageInt[FirstPoint.x + LargIm*FirstPoint.y] =-1;
						ConditionArret = FALSE; //on sort si la tete = la queue
						Rebouclage = TRUE; // pas la peine de traiter la connexité
					}

					bFinFreeman=true; // pas la peine de continuer	
				
			}
			else// si on a pas trouvé de point blanc
			{
				if(CodeFreeman == 7)// et qu'on a balayé tous les freemans alors on sort
				{
					if(BoolRetour)
					{
						ConditionArret = FALSE;
					}
					else
					{
						if(ListePointContour->GetCount() >1)
						{
							BoolRetour=TRUE;
							CopieImageInt[CurentPoint.x + LargIm*CurentPoint.y] =-1;
							PositionQueue = ListePointContour->GetHeadPosition();
							CurentPoint =ListePointContour->GetNext(PositionQueue);
							CurentPoint =ListePointContour->GetNext(PositionQueue);
							CodeFreeman = -1;
							NbZero=0;
							NbMoinsUn=0;

							i=0;
						}
						else
						{	
							ConditionArret = FALSE;
						}
					}
				}
			}
			CodeFreeman++;

		
		if (CodeFreeman >7) bFinFreeman=true;
		}while(!bFinFreeman);

		if ((NbZero==7) && (CandidatTrouve==FALSE))//extremité
		{
			CopieImageInt[CurentPoint.x + LargIm*CurentPoint.y] =-1;
			ConditionArret=false;
		}

		if ((NbZero<=7)&&(NbMoinsUn<=7)&&(CandidatTrouve==FALSE))//cas d'un point emprisonné
		{
			CopieImageInt[CurentPoint.x + LargIm*CurentPoint.y] =-1;
			/*pTemp->ImageRGB[CurentPoint.x + LargIm*CurentPoint.y].r = 0;
			pTemp->ImageRGB[CurentPoint.x + LargIm*CurentPoint.y].g = 0;
			pTemp->ImageRGB[CurentPoint.x + LargIm*CurentPoint.y].b = 255;
			pTemp->Affichage();
			Sleep(10);*/

			ConditionArret=false;
		}

		if ((NbZero==8)&&(CandidatTrouve==FALSE))//cas d'un point isolé
		{
			CopieImageInt[CurentPoint.x + LargIm*CurentPoint.y] =-1;
			/*pTemp->ImageRGB[CurentPoint.x + LargIm*CurentPoint.y].r = 0;
			pTemp->ImageRGB[CurentPoint.x + LargIm*CurentPoint.y].g = 255;
			pTemp->ImageRGB[CurentPoint.x + LargIm*CurentPoint.y].b = 0;
			pTemp->Affichage();
			Sleep(10);*/

			ConditionArret=false;
		}
		
	
	}while (ConditionArret); //condition d'arrêt; //condition d'arrêt

}//fin traitement cas U
	//MessageBox(NULL,"STOP","Après Retour Parcours freeman",MB_OK);

TabListePoint.Add(ListePointContour);
NbListeInArray++;
}// fin if point isolé
else
{
CopieImageInt[FirstPoint.x + LargIm*FirstPoint.y] =-1;
}

}//fin if(PointDepart)
	//MessageBox(NULL,"STOP","Après if point depart",MB_OK);



}while(PointDepart);

	//MessageBox(NULL,"STOP","Après tant qu'on a des points de depart",MB_OK);



	/*if (nbelement>1) 
		PasContour=FALSE;
	else
	{
		pImage[(CurentPoint.y*LargIm)+CurentPoint.x] =0;
		essai++;
		if (essai>100) PasContour=FALSE; 
	}*/

//} while (PasContour);




	free(CopieImageInt);
	return TRUE;
}

void CCurviDetect2D::RemplitTableauCoupe(ImagePetrou A, CPoint * TabCoupe, int x, int y, float dx, float dy, int w, int * Centre, long Hauteur, long Largeur)
{
	float Pente=0.0F, OrdonneeOrigine=0.0F;
	char equation='0';
	long signe;
	long i,j,l;
	int wx;
	int wy;
	*Centre=w;
	
	int NbPointCoupe =2*w+1;
	
	/*dy = fabs(dy);*/

	if (dx==0) equation = '1';   //la perp est verticale
	else if(fabs(dy)<1) equation = '2';   // la perp est horizontale
	 else equation = '3';   //equation de la forme y=ax+b
    
	 switch (equation)
	{  
	    case '1':
			*Centre=w;
			j=y-w;  // last change
			i=x;
            for (l=0; l<NbPointCoupe; l++)
            {   
				if (i>=A.dim_x) i=A.dim_x-1;
				if (j>=A.dim_y) j=A.dim_y-1;
				if (i<0) i=0;
				if (j<0) j=0;
				TabCoupe[l].x = (unsigned short) i;
                TabCoupe[l].y = (unsigned short) j;
			//	if (TabIJ[l].x==x && TabIJ[l].y==y) *Centre=l;
				j ++;
			}
			break;
		case '2':
			*Centre=w;
			i=x-w;
			j=y;
			for (l=0; l<NbPointCoupe; l++)
			{   
				if (i>=A.dim_x) i=A.dim_x-1;
				if (j>=A.dim_y) j=A.dim_y-1;
				if (i<0) i=0;
				if (j<0) j=0;
 				TabCoupe[l].x = (unsigned short) i;
                TabCoupe[l].y = (unsigned short) j;
			//	if (TabIJ[l].x==x && TabIJ[l].y==y) *Centre=l;
				i ++;
			}
			break;
		case '3':
			Pente = dy/dx; //initialise la pente a
			Pente = fabs(Pente);
			OrdonneeOrigine = y - Pente*x;  //initialise l'ordonnee b de y=a*x+b
 			if (Pente>0.0)
			{

				if (Pente>1.0)	
				{

					j = y-w;
					i =  (long)  (((float)j-OrdonneeOrigine)/Pente);
					for (l=0; l<NbPointCoupe; l++) 
					{   
						if (i>=A.dim_x) i=A.dim_x-1;
						if (j>=A.dim_y) j=A.dim_y-1;
						if (i<0) i=0;
						if (j<0) j=0;
						TabCoupe[l].x = (unsigned short) i;
                		TabCoupe[l].y = (unsigned short) j;
						if (TabCoupe[l].x==x && TabCoupe[l].y==y) *Centre=l;
						j++;   
						i = (long)  (((float)j-OrdonneeOrigine)/Pente);   

					}	
				}
				else
				{	
					i =  x-w;
					j = Pente*(float)i+OrdonneeOrigine;
					for (l=0; l<NbPointCoupe; l++)         
					{
						if (i>=A.dim_x) i=A.dim_x-1;
						if (j>=A.dim_y) j=A.dim_y-1;
						if (i<0) i=0;
						if (j<0) j=0;
						TabCoupe[l].x = (unsigned short) i;
               			TabCoupe[l].y = (unsigned short) j;
						if (TabCoupe[l].x==x && TabCoupe[l].y==y) *Centre=l;
					 	i ++ ;
					 	j = (long)(Pente*(float)i+OrdonneeOrigine);
					}                                   
				}	
			}
			else
			{   
				if (Pente<-1.0)
				{
					j = y+w;
					i =  (long)  (((float)j-OrdonneeOrigine)/Pente);
					for (l=0; l<NbPointCoupe; l++)   
					{  
						if (i>=A.dim_x) i=A.dim_x-1;
						if (j>=A.dim_y) j=A.dim_y-1;
						if (i<0) i=0;
						if (j<0) j=0;
						TabCoupe[l].x = (unsigned short) i;
                		TabCoupe[l].y = (unsigned short) j;
						if (TabCoupe[l].x==x && TabCoupe[l].y==y) *Centre=l;
						j--;   
						i = (long) (((float)j-OrdonneeOrigine)/Pente);   
					}	
				}
				else
				{
					i =  x-w;
					j = Pente*(float)i+OrdonneeOrigine;
					for (l=0; l<NbPointCoupe; l++)     
					{
						if (i>=A.dim_x) i=A.dim_x-1;
						if (j>=A.dim_y) j=A.dim_y-1;
						if (i<0) i=0;
						if (j<0) j=0;
						TabCoupe[l].x = (unsigned short) i;
						TabCoupe[l].y = (unsigned short) j;
						if (TabCoupe[l].x==x && TabCoupe[l].y==y) *Centre=l;
					 	i++ ;
					 	j = (long)(Pente*(float)i+OrdonneeOrigine);
					}	
				}
			
			}	
			break;
			
		default:
			break;
	}
}
void CCurviDetect2D::RemplitTableauIJCoupeBiDir(CPoint *TabIJ,int x,int y,float dx,float dy,int w,int dw,int Dir,int *IptInit, long Hauteur, long Largeur)
{
float Pente=0.0F, OrdonneeOrigine=0.0F;
float Pente2=0.0F, OrdonneeOrigine2=0.0F;
char equation='0';
long signe;
long i,j,l,ip,jp;
int wx;
int wy;

//*IptInit=dw;

*IptInit=w;

//int Largeur;
//int Hauteur;
//Largeur=A.dim_x;
//Hauteur=A.dim_y;
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

} // fin fonction

BOOL CCurviDetect2D::TestCroisement(CPoint i1, CPoint i2, CPoint f1, CPoint f2)
{
	BOOL Cond = TRUE;
	CPoint c1,c2;
	int difx,dify;
	c1.x=(i1.x+f2.x)/2;
	c2.x=(i2.x+f1.x)/2;
	c1.y=(i1.y+f2.y)/2;
	c2.y=(i2.y+f1.y)/2;
	difx=abs(c1.x-c2.x);
	dify=abs(c1.y-c2.y);
	if (difx<3 && dify<3)
	//	if (abs((i1.x+f2.x)-(i2.x+f1.x))<3 && abs((i1.y+f2.y)-(i2.y+f1.y))<3)
	{
		Cond=FALSE;
	}
	return (Cond);
}

void CCurviDetect2D::InvertionPoint(CPoint * xp1, CPoint * xp2)
{
	CPoint Temp;
	Temp.x = xp1->x;
	Temp.y = xp1->y;
	//xp1=xp2;
	// correction jo
	xp1->x=xp2->x;
	xp1->y=xp2->y;
	xp2->x= Temp.x;
	xp2->y= Temp.y;
}

void CCurviDetect2D::RechercheRectEnglobant(CPoint * TabPoint, long * Xmin, long * Xmax, long * Ymin, long * Ymax)
{
	int i=0;
	//init
	*Xmin = TabPoint[0].x;
	*Xmax = TabPoint[0].x;
	*Ymin = TabPoint[0].y;
	*Ymax = TabPoint[0].y;
	//recherche maxi et mini
	for(i=0;i<4;i++)
	{
		if(TabPoint[i].x<*Xmin)//mini en x
		{
			*Xmin=TabPoint[i].x;
		}
		if(TabPoint[i].x>*Xmax)//max en x
		{
			*Xmax=TabPoint[i].x;
		}
		if(TabPoint[i].y<*Ymin)//mini en y
		{
			*Ymin=TabPoint[i].y;
		}
		if(TabPoint[i].y>*Ymax)//max en y
		{
			*Ymax=TabPoint[i].y;
		}
	}
	
}

void CCurviDetect2D::AfficheQuad(CImageWnd * pImage,CRgn * RegionAffichage,long Larg, long * Xmin, long * Xmax, long * Ymin, long * Ymax, int R, int G, int B)
{	
	int i=0,j=0;
	int p;

	for(i=(*Ymin);i<(*Ymax);i++)
	{
		for(j=(*Xmin);j<(*Xmax);j++)
		{
			
			if (RegionAffichage->PtInRegion(j,i))
			{
				p=i*Larg+j;
				pImage->ImageRGB[p].r =R;
				pImage->ImageRGB[p].g = G;
				pImage->ImageRGB[p].b = B;
			}
			
			
		}	
	}
}

void CCurviDetect2D::RempliQuad(int* pPlan,CRgn * RegionAffichage,long Larg, long * Xmin, long * Xmax, long * Ymin, long * Ymax, int NumReg)
{	int i=0,j=0;
	for(i=(*Ymin);i<(*Ymax);i++)
	{
		for(j=(*Xmin);j<(*Xmax);j++)
		{
			if (RegionAffichage->PtInRegion(j,i))
			{
				pPlan[i*Larg+j] =NumReg;
			}
		}	
	}
}

void CCurviDetect2D::MesureRepetabilite(rgb * Img1, rgb * Img2, long Largeur,long Hauteur,long NBRegionCurviImg1,long NBRegionCurviImg2,double * H)
{
//--------------------------------------------------------
//--------------------------------------------------------
//			mesure de repetabilité
//--------------------------------------------------------


// sens 1
long * ScoreUnionAll = (long *)malloc(sizeof(long)*NBRegionCurviImg1*NBRegionCurviImg2);
long *ScoreInterAll = (long *)malloc(sizeof(long)*NBRegionCurviImg1*NBRegionCurviImg2);
long * ScoreUnion= (long *)malloc(sizeof(long)*NBRegionCurviImg1);
long *ScoreInter = (long *)malloc(sizeof(long)*NBRegionCurviImg1);
double *ScoreRepet = (double *)malloc(sizeof(double)*NBRegionCurviImg1);
long UnionAprendre=0, MaxInter=0;
long w,k,i,j,cpt=0, CouleurRegion1,CouleurRegion2;
double ip,jp,zp;
CStdioFile FichierExportRept;
CString ChaineRepet;
for(i=0;i<NBRegionCurviImg1*NBRegionCurviImg2;i++)
{
	ScoreUnionAll[i]=0;
	ScoreInterAll[i]=0;
}

for(i=0;i<NBRegionCurviImg1;i++)
{
	ScoreUnion[i]=0;
	ScoreInter[i]=0;
	ScoreRepet[i]=0.0f;
}
cpt=0;
for(w=1;w<NBRegionCurviImg1+1;w++)// boucle couleur regions Image 1
{
	for(k=1;k<NBRegionCurviImg2+1;k++)// boucle couleur regions Image 2
	{
		for(i=0;i<Hauteur;i++)// ligne colone
		 {
			for(j=0;j<Largeur;j++)
			{
				// calcul des coordonnées dans l'image 2
				/*jp = H[0]*(double)j + H[1]*(double)i +H[2];
				ip = H[3]*(double)j + H[4]*(double)i +H[5];
				zp = H[6]*(double)j + H[7]*(double)i +H[8];

				ip = ip/zp;
				jp = jp/zp;*/
				
				
					CouleurRegion1 = (long)Img1[i*Largeur+j].r;//((long)Img1[j*Largeur+i].r)+((long)Img1[j*Largeur+i].g)<<8+((long)Img1[j*Largeur+i].b)<<16;
					CouleurRegion2 = (long)Img2[i*Largeur+j].r; //((long)Img2[(long)jp*Largeur+(long)ip].r)+((long)Img2[(long)jp*Largeur+(long)ip].g)<<8+((long)Img2[(long)jp*Largeur+(long)ip].b)<<16;
					if (CouleurRegion1==2)
					{
						int stop=0;
					}

					if((CouleurRegion1 == w)&&(CouleurRegion2==k) )
					{
						ScoreInterAll[cpt]++;
					}
					if((CouleurRegion1 == w)||(CouleurRegion2==k))
					{
					ScoreUnionAll[cpt]++;
					}
				
			}//fin for j
		}//fin for i

		cpt++;
	}//fin for k
}//fin for w
cpt=0;

if(FichierExportRept.Open("RRepet.txt",CFile::modeCreate|CFile::modeWrite))
{
	UnionAprendre=ScoreUnionAll[1]; MaxInter=0;
	for(w=1;w<NBRegionCurviImg1+1;w++)// boucle couleur regions Image 1
	{
		for(k=1;k<NBRegionCurviImg2+1;k++)// boucle couleur regions Image 2
		{
			if(ScoreInterAll[cpt]>MaxInter)
			{
				MaxInter=ScoreInterAll[cpt];
				UnionAprendre=ScoreUnionAll[cpt];
			}
				
			cpt++;
		}//fin for k	
		ScoreInter[w-1]=MaxInter;
		ScoreUnion[w-1]=UnionAprendre;
		ScoreRepet[w-1]=(double)MaxInter/(double)UnionAprendre;
		ChaineRepet.Format("%f",ScoreRepet[w-1]);
		FichierExportRept.WriteString(ChaineRepet);
		FichierExportRept.WriteString("\n");
	}//fin for w
	FichierExportRept.Close();
}






free(ScoreUnionAll);
free(ScoreInterAll);
free(ScoreUnion);//bug
free(ScoreInter);
free(ScoreRepet);

}


void CCurviDetect2D::MesureRepetJ(CIma * Ima1 , CIma * Ima2 , CString NomFichier, long NbReg1, long NbReg2)
{
CStdioFile fic;

CString Temp;

fic.Open(NomFichier,CFile::modeWrite|CFile::modeCreate);

int denom,max,l,c,i,j,w,k,memk;
int **HistoInter;
double Mesure;
long * HistoSurfIm1;
long * HistoSurfIm2;


HistoSurfIm1 = (long*) malloc(sizeof(long)*(NbReg1+1));
HistoSurfIm2 = (long*) malloc(sizeof(long)*(NbReg1+1));

HistoInter=(int**) malloc((NbReg1+1)*sizeof(int*));
for (i=0;i<=NbReg1;i++) HistoInter[i]=(int*)malloc((NbReg2+1)*sizeof(int));

for (w=0;w<=NbReg1;w++)
{
	for (k=0;k<=NbReg2;k++)
	{
	HistoInter[w][k]=0;
	}
	HistoSurfIm1[w]=0; 
	HistoSurfIm2[w]=0;
}

for (l=0;l<Ima1->m_Hauteur;l++)
{
	for (c=0;c<Ima1->m_Largeur;c++)
	{
		i=Ima1->m_PImage[Ima1->m_Largeur*l+c].r;
		j=Ima2->m_PImage[Ima2->m_Largeur*l+c].r;
		HistoInter[i][j]++;
		HistoSurfIm1[i]++;
		HistoSurfIm2[j]++;
	}

}

for (w=1;w<=NbReg1;w++)
{
	max=HistoInter[w][1];
	memk=1;
	for (k=1;k<=NbReg2;k++)
	{
		if (max<HistoInter[w][k])
		{
			max=HistoInter[w][k];
			memk=k;
		}
	}

	denom=HistoSurfIm1[w]+HistoSurfIm2[memk]-max;
	if(denom!=0) Mesure=100.0*(double)max/(double)(denom);
	else
	{
		if (max==0) Mesure = 0;
		else Mesure=100.0;
	}

	Temp.Format("%f\n",Mesure);
	fic.WriteString(Temp);
}

fic.Close();

for (i=0;i<NbReg1+1;i++) free(HistoInter[i]);
free(HistoInter);
free(HistoSurfIm1);
free(HistoSurfIm2);

}

int CCurviDetect2D::MesureRecouv(const int * Ima1 ,const int * Ima2 , int Larg,int Haut,long NbReg1, double SeuilRecouv, int *pi,int *TabListeOk)
{

int p,denom,max,l,c,i,j,w,k,memk;
int *HistoInter;
double Mesure;
long * HistoSurfIm1;
long HistoSurfIm2;
double MesureMax=0.0;
int SurfMin=0;
int Sens;
int Taille=Larg*Haut;

*pi=0;

HistoSurfIm1 = (long*) malloc(sizeof(long)*(NbReg1+1));
HistoSurfIm2 =0;

HistoInter=(int*) malloc((NbReg1+1)*sizeof(int));

for (w=0;w<=NbReg1;w++)
{
	HistoInter[w]=0;
	HistoSurfIm1[w]=0; 
}
HistoSurfIm2=0;

for (p=0;p<Taille;p++)
{
		if (Ima2[p]!=0) {HistoInter[Ima1[p]]++;HistoSurfIm2++;}
		HistoSurfIm1[Ima1[p]]++;
}

max=HistoInter[1];
memk=1;
for (w=1;w<=NbReg1;w++)
{
	if (max<HistoInter[w] && TabListeOk[w-1]==1)
	{
		max=HistoInter[w];
		memk=w;
	}
}
if (HistoSurfIm1[memk]==0)
{

	int clam=0;
}
if (HistoSurfIm1[memk]>HistoSurfIm2) 
{
	SurfMin=HistoSurfIm2;
}
else 
{
	SurfMin=HistoSurfIm1[memk];
}
Mesure = (double)max/(double)(SurfMin);

if (Mesure>SeuilRecouv && HistoSurfIm1[memk]>HistoSurfIm2)
{
	Sens=12; // on garde la 1 et on vire la 2
}

if (Mesure>SeuilRecouv && HistoSurfIm1[memk]<=HistoSurfIm2)
{
	Sens=21; //on garde la 2 et on vire la 1
}

if (Mesure<=SeuilRecouv) Sens=0;

free(HistoInter);
free(HistoSurfIm1);

*pi=memk;
return Sens;
}


double CCurviDetect2D::AngleCourbureLocale(Sfpoint P1,Sfpoint P2,Sfpoint P3)
{
	
	// pythagore dans un triangle non rectangle
	double Alpha=0.0;
	double cosb,p1p2,p1p3,p2p3;
	p1p2= sqrt(
	(P1.x-P2.x)*(P1.x-P2.x)+(P1.y-P2.y)*(P1.y-P2.y));
	p1p3= sqrt(
	(P1.x-P3.x)*(P1.x-P3.x)+ (P1.y-P3.y)*(P1.y-P3.y) );
	p2p3= sqrt(
	(P3.x-P2.x)*(P3.x-P2.x)+  (P3.y-P2.y)*(P3.y-P2.y));

	cosb=(p1p2*p1p2 + p2p3*p2p3- p1p3*p1p3);
	cosb=cosb/(2.0*p1p2*p2p3);
	Alpha=acos(cosb);
	return Alpha;
}