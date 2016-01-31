#pragma once

class CSteger1D
{
public:
	CSteger1D(void);
	~CSteger1D(void);
	// Tableau de la gaussienne 1D
	double *GaussianKernel;
	// Tableau Derivée première 
	double* GaussianKernelD1;
	// Tableau Dérivée seconde
	double *GaussianKernelD2;
	// Tableau pour la dérivée troisième qu'on pourra éventuellement utilisée pour la détection de contours
	double *GaussianKernelD3;
	// Ecart-Type
	double Theta;
	double *ResultatFiltrage;
	void InitKernels(double ValTheta, double ValTaille);
	double Pi;
	double DeuxPi;
	double *TabResult;
	// Fonction réalisant la convolution avec les gaussiennes
	void ConvolutionKernels(double * Signal, long TailleSignal, long TailleKernel,double ValTheta);
	double *TabResultD1;
	double *TabResultD2;
	void EffaceTab(void);
	void InitTabResult(long TailleSignal);
};
