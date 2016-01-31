#include "echant.h"


void SousEchantillonnageLigneColonne(SImage Source,SImage Dest)
{
	long cpt,i,j,k;
	cpt=0;
	for (j=0;j<Source.Haut;j+=2)
	{   
		k = 0;     
		for(i=0;i<Dest.Larg;i++)
		{
			Dest.Image[cpt++]=Source.Image[k+j*Source.Larg];
			k +=2;
		}
	}               

}

void SousEchantillonnageLigne(SImage Source,SImage Dest)
{
	long cpt,i,j,k;
	cpt=0;
	for (j=0;j<Source.Haut;j++)
	{   
		k = 0;     
		for(i=0;i<Dest.Larg;i++)
		{
			Dest.Image[cpt++]=Source.Image[k+j*Source.Larg];
			k +=2;
		}
	}               

} 

void SousEchantillonnageColonne(SImage Source,SImage Dest)
{
	long cpt,i,j;
	cpt=0;
	for (j=0;j<Source.Haut-1;j=j+2)
	{        
		for(i=0;i<Source.Larg;i++)
		{
			Dest.Image[cpt++]=Source.Image[i+j*Source.Larg];
		}
	}               

} 
