#include "stdafx.h"
#include "Steger1D.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <conio.h>


CSteger1D::CSteger1D(void)
: TabResult(NULL)
, TabResultD1(NULL)
, TabResultD2(NULL)
{
	Pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
	DeuxPi = 2*Pi;
}

CSteger1D::~CSteger1D(void)
{

	//Desallocation faite dans le faite dans la fonction efface tab
	/*free(GaussianKernel);
	free(GaussianKernelD1);
	free(GaussianKernelD2);
	free(GaussianKernelD3);
	free(TabResult);
	free(TabResultD1);
	free(TabResultD2);*/
}

void CSteger1D::InitKernels(double ValTheta, double ValTaille)
{
	//déclarations locales 
	int i=0;
	//allocation des tableaux pour les filtres gaussiens
	GaussianKernel = (double*)malloc(sizeof(double)*ValTaille);
	GaussianKernelD1 = (double*)malloc(sizeof(double)*ValTaille);
	GaussianKernelD2 = (double*)malloc(sizeof(double)*ValTaille);
	GaussianKernelD3 = (double*)malloc(sizeof(double)*ValTaille);

	for(i=0;i<(ValTaille);i++)
	{
		int x = i-(ValTaille/2);
		GaussianKernel[i] = ((1)/(DeuxPi*ValTheta))
			*(exp((-1*x*x)/(2*ValTheta*ValTheta)));
		GaussianKernelD1[i] = ((-1*x)/(DeuxPi*ValTheta*ValTheta))
			*(exp((-1*x*x)/(2*ValTheta*ValTheta)));
		GaussianKernelD2[i] = (((x*x)-(ValTheta*ValTheta))/(DeuxPi*ValTheta*ValTheta*ValTheta*ValTheta*ValTheta))
			*(exp((-1*x*x)/(2*ValTheta*ValTheta)));
	}
}

// Fonction réalisant la convolution avec les gaussiennes
void CSteger1D::ConvolutionKernels(double * Signal, long TailleSignal, long TailleKernel, double ValTheta)
{
	//déclarations variables locales
	InitKernels(ValTheta,TailleKernel);
	int i=0,j=0,DebutBoucle = TailleKernel;

	InitTabResult(TailleSignal);

	for(i=DebutBoucle;i<TailleSignal;i++)
	{
		for(j=0;j<TailleKernel;j++)
		{
			TabResult[i] += Signal[i-j]*GaussianKernel[j];
			TabResultD1[i] += Signal[i-j]*GaussianKernelD1[j];
			TabResultD2[i] += Signal[i-j]*GaussianKernelD2[j];
		}
	}
	
}

void CSteger1D::EffaceTab(void)
{
	free(TabResult);
	free(TabResultD1);
	free(TabResultD2);
	free(GaussianKernel);
	free(GaussianKernelD1);
	free(GaussianKernelD2);
}

void CSteger1D::InitTabResult(long TailleSignal)
{
	TabResult = (double*)malloc(sizeof(double)*TailleSignal);
	TabResultD1 = (double*)malloc(sizeof(double)*TailleSignal);
	TabResultD2 = (double*)malloc(sizeof(double)*TailleSignal);
	for(int i=0;i<TailleSignal;i++)
	{
		TabResult[i] = 0.0f;
		TabResultD1[i] = 0.0f;
		TabResultD2[i] = 0.0f;
	}
}
