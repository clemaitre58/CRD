// Donnee.cpp: implementation of the CDonneeFloat class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "float.h"
//#include "Decid.h"
#include "DonneeFloat.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CDonneeFloat::CDonneeFloat(CString NomFichier)
{
	int c;
	Erreur=1;
	Tab=NULL;
	Tmin=NULL;
	Tmax=NULL;
	for (c=0;c<256;c++) Tab[c]=NULL;
	NbClass=0;
	LireDonnee(NomFichier);
	NbMesTot=0;

}

CDonneeFloat::CDonneeFloat()
{
//	int c;
	Erreur=1;
	Tab=NULL;
	Tmin=NULL;
	Tmax=NULL;
//	for (c=0;c<256;c++) Tab[c]=NULL;
	NbClass=0;
	NbMesTot=0;

}


CDonneeFloat::~CDonneeFloat()
{
	BYTE c;
	for (c=0;c<NbClass;c++)
	{
		if (Tab[c]!=NULL) free(Tab[c]);
	}
	if (Tab!=NULL) free(Tab);
	if (Tmin!=NULL) free(Tmin);
	if (Tmax!=NULL) free(Tmax);

}


int CDonneeFloat::LireDonnee(CString NomFichier)
{
	DWORD NbMesures;
	BYTE c;
	CFile fi;
	long kmax01=0,i,par;
	float val;
	/* si déjà lu*/
	for (c=0;c<NbClass;c++)
	{
		if (Tab[c]!=NULL) free(Tab[c]);
	}
	if (Tab!=NULL) free(Tab);
	
	/*lecture des donnees et allocation memoire */
	if(fi.Open(NomFichier,CFile::modeRead))
	 {
  		fi.Read(&NbClass,1);
		fi.Read(&NPAR,1);
		Tab=(float **)malloc(NbClass*sizeof(float*));
		Tmin=(float *) malloc(NPAR*sizeof(float));
		Tmax=(float *) malloc(NPAR*sizeof(float));
		for (c=0;c<NbClass;c++)
		{
			/*attention fi.Read(&NbMesures,2);*/
			fi.Read(&NbMesures,4);
			NbMes[c]=NbMesures;
			NbMesTot+=(long)NbMesures;
			kmax[c]=(long)NbMes[c]*(long)NPAR;
			if (kmax[c]>kmax01) kmax01=kmax[c];
			Tab[c]=(float *)malloc(kmax[c]*sizeof(float));
			if (Tab[c]==NULL)   {MessageBox(NULL,"Erreur allocation données","Message",MB_OK);Erreur=-2;}
			fi.Read(Tab[c],kmax[c]*sizeof(float));

			/*essai*/
		//	for(i=0;i<NbMes[c];i++)	for (par=0;par<NPAR;par++) 	if(Tab[c][i*NPAR+par]>0.0) Tab[c][i*NPAR+par]=0.01;
		} // for classes 

		fi.Close();

		for (par=0;par<NPAR;par++)
		{
			Tmax[par]=-FLT_MIN;
			Tmin[par]=FLT_MAX;
			for (c=0;c<NbClass;c++)
			{
				for(i=0;i<NbMes[c];i++)
				{
					val=Tab[c][i*NPAR+par];
					if (_isnan(val) || _isnan(-val) ) 
					{
						val=0; 
						MessageBox(NULL,"attention","infini detecté",MB_OK);
						Tab[c][i*NPAR+par]=0;
					}
					if(val>Tmax[par]) Tmax[par]=val;
					if (val<Tmin[par]) Tmin[par]=val;
				}
			} // pour les classs
			
		}//pour les par

	} // if
	else {MessageBox(NULL,"Erreur fichier","Message",MB_OK);Erreur=-1;}
//	Normalise();
	return Erreur;
}


void CDonneeFloat::NormaliseSimple()
{
	int i,c,par;
	float max=-FLT_MIN;
	float min=FLT_MAX;
	float coef,val;	
	for (par=0;par<NPAR;par++)
	{
		max=-FLT_MIN;
		min=FLT_MAX;
		for (c=0;c<NbClass;c++)
		{

			for(i=0;i<NbMes[c];i++)
			{
				val=Tab[c][i*NPAR+par];
				if(val>max) max=val;
				if (val<min) min=val;
			}
		} // pour les classs

		coef=(max-min)/0.8F;

		for (c=0;c<NbClass;c++)
		{
			for(i=0;i<NbMes[c];i++) 
			{
				if (coef!=0) Tab[c][i*NPAR+par]=0.1F+((Tab[c][i*NPAR+par]-min)/coef);
				else Tab[c][i*NPAR+par]=0.1F;
			}
		}

	} // pour les param

}

void CDonneeFloat::Normalise(float binf, float bsup)
{
	int i,c,par;
	float coef;
	float val;
	CString ch;
	if (bsup<=binf) {MessageBox(NULL,"Attention","erreur normalisation",MB_OK);return;}
	for (par=0;par<NPAR;par++)
	{
		coef=(Tmax[par]-Tmin[par])/(bsup-binf);
		for (c=0;c<NbClass;c++)
		{
			for(i=0;i<NbMes[c];i++) 
			{
				val=Tab[c][i*NPAR+par];
					if (_isnan(val) || _isnan(-val) ) 
					{val=0; MessageBox(NULL,"attention","infini detect",MB_OK);}
			//		if (val==0.0) 
			//		{ ch.Format("%d",i);MessageBox(NULL,"attention",ch,MB_OK);}
				if (val<Tmin[par]) val=Tmin[par];
				else if (val>Tmax[par]) val=Tmax[par];
				else 
				{
					if (coef!=0) val=binf+((val-Tmin[par])/coef);
					else val=binf;
				}
				Tab[c][i*NPAR+par]=val;
			}
		}

	} // pour les param
}

void CDonneeFloat::SauveMinMax(CString NomFichier)
{
	CFile fi;
	if(fi.Open(NomFichier,CFile::modeWrite|CFile::modeCreate))
	{
		fi.Write(&NPAR,sizeof(NPAR));
		fi.Write(Tmin,NPAR*sizeof(float));
		fi.Write(Tmax,NPAR*sizeof(float));
		fi.Close();
	}

}


void CDonneeFloat::LitMinMax(CString NomFichier)
{
	CFile fi;
	if(fi.Open(NomFichier,CFile::modeRead))
	{
		fi.Read(&NPAR,sizeof(NPAR));
		if (Tmin!=NULL) free(Tmin);
		if (Tmax!=NULL) free(Tmax);
		Tmin=(float *) malloc(NPAR*sizeof(float));
		Tmax=(float *) malloc(NPAR*sizeof(float));
		fi.Read(Tmin,NPAR*sizeof(float));
		fi.Read(Tmax,NPAR*sizeof(float));
		fi.Close();
	}

}