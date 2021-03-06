// polymc.cpp : implementation file
//


#include "stdafx.h"
//#include "winseg.h"
#include "polmc.h"




/////////////////////////////////////////////////////////////////////////////
// CPolyMc
//BEGIN_MESSAGE_MAP(CPolyMc, CDocument)
	//{{AFX_MSG_MAP(CPolyMc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
//END_MESSAGE_MAP()


CPolyMc::CPolyMc()
{

polstemp=NULL;
politemp=NULL;
hypi=NULL;
hyps=NULL;
poli=NULL;
pols=NULL;
}


CPolyMc::~CPolyMc()
{

}

void CPolyMc::FinAlloc()
{
int c;
for (c=0;c<NbClass;c++)
{
	free(Va[c]);
	free(poli[c]);
	free(pols[c]);
	free(politemp[c]);
	free(polstemp[c]);
}
free(poli);
free(pols);
free(Va);
free(politemp);
free(polstemp);
if (hypi!=NULL) free(hypi);
if (hyps!=NULL) free(hyps);

}




int CPolyMc::Init(CString NomFich)
{
  CFile fi;
  long kmax01=0,c;
  /*attention short _int*/DWORD NbMesures;

  if(fi.Open(NomFich,CFile::modeRead))
  {
  		fi.Read(&NbClass,1);
		fi.Read(&NPAR,1);
		
		Va=(byte **)malloc(NbClass*sizeof(byte*));
		poli=(byte **)malloc(NbClass*sizeof(byte*));
		pols=(byte **)malloc(NbClass*sizeof(byte*));
		politemp=(byte **)malloc(NbClass*sizeof(byte*));
		polstemp=(byte **)malloc(NbClass*sizeof(byte*));

		for (c=0;c<NbClass;c++)
		{
			/*attention fi.Read(&NbMesures,2);*/
			fi.Read(&NbMesures,4);
			NbMes[c]=NbMesures;
			kmax[c]=(long)NbMes[c]*(long)NPAR;
			if (kmax[c]>kmax01) kmax01=kmax[c];

			Va[c]=(byte *)malloc(kmax[c]);
			if (Va[c]==NULL)   {MessageBox(NULL,"Erreur va","Message",MB_OK);return 0;}
			fi.Read(Va[c],kmax[c]);
			
			poli[c]=(byte *)malloc(kmax[c]);
			if (poli[c]==NULL) {MessageBox(NULL,"Erreur poli","Message",MB_OK);return 0;}
			pols[c]=(byte *)malloc(kmax[c]);
			if (pols[c]==NULL) {MessageBox(NULL,"Erreur pols","Message",MB_OK);return 0;}
			politemp[c]=(byte  *)malloc(kmax[c]);
			if (politemp[c]==NULL) {MessageBox(NULL,"Erreur politemp","Message",MB_OK);return 0;}
			polstemp[c]=(byte  *)malloc(kmax[c]);
    		if (polstemp[c]==NULL) {MessageBox(NULL,"Erreur polstemp","Message",MB_OK);return 0;}

		}//i

	  hypi=(byte    *)malloc(kmax01);
	  if (hypi==NULL) {MessageBox(NULL,"Erreur hypi","Message",MB_OK);return 0;}
	  hyps=(byte    *)malloc(kmax01);
	  if (hyps==NULL) {MessageBox(NULL,"Erreur hyps","Message",MB_OK);return 0;}


	  fi.Close();;

	  
  }else  return 0;
BORD=0;
NPOL=0;
return 1;
}  // INIT
  

void CPolyMc::AffMes()
{
	CString ch1,ch;
	ch1.Format("Nombre de classes :%d\nNombre de parametres  :%d\n\n",NbClass,NPAR);
	for(int i=0;i<NbClass;i++)
	 {
		ch.Format("Nb mesures de classe %d : %ld\n",i,NbMes[i]);
		ch1+=ch;
	}
	MessageBox(NULL,ch1,"Info",MB_OK);
}
  

void CPolyMc::Fin()
{       

 if (poli!=NULL) free(poli);
 if (pols!=NULL) free(pols);

} // fin
  

void CPolyMc::AffPol()
{
  CString ch1,ch;
  ch1.Format("Nombre de classes :%d\nNombre de parametres  :%d\nCoefficient :%u\n",NbClass,NPAR,coef);
  for(int i=0;i<NbClass;i++)
	 {
		ch.Format("Nb pol. classe %d : %ld\n",i,Npol[i]);
		ch1+=ch;
	}
	MessageBox(NULL,ch1,NomPol.Right(NomPol.GetLength()-NomPol.ReverseFind('\\')-1),MB_OK);
} // AffPol


int CPolyMc::LirePol(CString NomIm)
{
	int c;

    CFile fp;
	long i;//,j;
	byte nva;//,nbn;

	CString Cval,CTemp;
  	NomPol=NomIm;
	
	if (fp.Open(NomIm,CFile::modeRead))
	{ 	
		fp.Read(&NbClass,1);
		poli=(byte **)malloc(NbClass*sizeof(byte*));
		pols=(byte **)malloc(NbClass*sizeof(byte*));
		fp.Read(&NPAR,1);
		
		for (c=0;c<NbClass;c++)
		{
	 		fp.Read(&NPOL,2);
			Npol[c]=NPOL;
			kmax[c]=(long)NPOL*(long)NPAR;
			poli[c]=(byte *)malloc(kmax[c]);
			if (poli[c]==NULL) {MessageBox(NULL,"Erreur poli","Message",MB_OK);}
			pols[c]=(byte *)malloc(kmax[c]);
			if (pols[c]==NULL) {MessageBox(NULL,"Erreur pols","Message",MB_OK);return 1;}
			for (i=0;i<(long)NPOL*(long)NPAR;i++)
			{
				fp.Read(&nva,1);
				poli[c][i]=nva;
			} // i 
			for (i=0;i<(long)NPOL*(long)NPAR;i++)
			{
				fp.Read(&nva,1);
				pols[c][i]=nva;
			} // i
		}	  //c
		
		fp.Read(&coef,2);
		fp.Close();
	}
	else return 0;
 	return 1;


}// fin


  
void CPolyMc::sauve_pol(CString NomF)
{
                           
byte nva;                           
long i,c;//,n;                           
FILE *fi;
CString Cval,CTemp;

fi=fopen(LPCSTR(NomF),"wb");
if (fi!=NULL)
{
  fwrite(&NbClass,1,1,fi);
  fwrite(&NPAR,1,1,fi);
  for (c=0;c<NbClass;c++)
  {
	fwrite(&Npol[c],2,1,fi); //Npol[c]
	for (i=0;i<(long)Npol[c]*(long)NPAR;i++)
	{
	nva=poli[c][i];	// poli
	fwrite(&nva,1,1,fi); 
	}
	for (i=0;i<(long)Npol[c]*(long)NPAR;i++)
	{
	nva=pols[c][i];		//pols
	fwrite(&nva,1,1,fi); 
	}
  }// c
  fwrite(&coef,2,1,fi);

 
  fclose(fi);
}
else {MessageBox(NULL,"Erreur fichier","Message",MB_OK);return;}

}

void CPolyMc::sauve_pol_avant(CString NomF)
{
                           
byte nva;                           
long i,c;//,n;                           
FILE *fi;
CString Cval,CTemp;
for (i=0;i<75;i++) texte[i]='p';                           

fi=fopen(LPCSTR(NomF),"wb");
if (fi!=NULL)
{
  fwrite(&NbClass,1,1,fi);
  fwrite(&NPAR,1,1,fi);
  for (c=0;c<NbClass;c++)
  {/*attention fwrite(&NbMes[c],2,1,fi)*/
	fwrite(&NbMes[c],4,1,fi); //Npol[c]
	for (i=0;i<(long)NbMes[c]*(long)NPAR;i++)
	{
	nva=politemp[c][i];	// poli
	fwrite(&nva,1,1,fi); 
	}
	for (i=0;i<(long)NbMes[c]*(long)NPAR;i++)
	{
	nva=polstemp[c][i];		//pols
	fwrite(&nva,1,1,fi); 
	}
  }// c
  fwrite(&coef,2,1,fi);

 
  fclose(fi);
}
else {MessageBox(NULL,"Erreur fichier","Message",MB_OK);return;}

}


/*******************************************************************/

void CPolyMc::calcb(int cl,long kmes)
{ 
  int par;
  byte val,vol;
  byte mes[0x100];
  //mes=malloc(0x100);
  int por,dir,clasop;
  byte dmx,dis;
  long k0;
   
 // if (cl) kmax=kmax0; else kmax=kmax1;
  
  for(par=0;par<NPAR;par++)
    { 
	  hypi[par]=0; hyps[par]=255;
      mes[par]=Va[cl][kmes+par];
	  if (mes[par]==0) mes[par]=1;
	  if (mes[par]==255) mes[par]=254;
	}
  
  k0=0;
  for (clasop=0;clasop<NbClass;clasop++)
  {
  
	  if (clasop!=cl)
	  {
		//printf("\n classe opposee=%d",clasop);

		 for(k0=0;k0<kmax[clasop];k0+=NPAR)   
		{
			 por=0;dir=0;dmx=0;
			 for(par=0;par<NPAR;par++)
			 { val=mes[par];
			   vol=Va[clasop][k0+par];
			   if (vol>val)
				  {dis=vol-val; if(dis>dmx) {dmx=dis; por=par; dir=0;} }
			   else
				  {dis=val-vol; if(dis>dmx) {dmx=dis; por=par; dir=1;} }
			 }//for par
			 //**********************************
			 vol=Va[clasop][k0+por];	 // relect du pt avec dist max
			 if (!dir) { if (hyps[por]>vol) hyps[por]=vol; }
			 else      { if (hypi[por]<vol) hypi[por]=vol; }
		}//for k0
	   } // si autre classe
  }	// pour autres classes
  //************************************
  /////////////////////////////////////////////////////////////////
  //////////ATTENTION MODIF////////////////////////////////////////
  /////////////////////////////////////////////////////////////////
  for(par=0;par<NPAR;par++)
  { val=mes[par];
    vol=hypi[par];
    /*if (vol)*/ hypi[par]=vol+(val-vol)*coef/10;
    vol=hyps[par];
    /*if (vol!=0xFF)*/ hyps[par]=vol-(vol-val)*coef/10;
    politemp[cl][kmes+par]=hypi[par];
    polstemp[cl][kmes+par]=hyps[par];
  }//for par
  //************************************
  //printf("\r%5d",kmes/NPAR);
  //free(mes);
}     // calcb

void CPolyMc::calcb_all()
{ 
  int c;//,par;
  long m0;
  long kmes;
  for (c=0;c<NbClass;c++)
  {
  	//printf("\n classe initiale=%d",c);
	kmes=0;
	for (m0=0;m0<NbMes[c];m0++)
	{
		calcb(c,kmes);
		kmes=kmes+NPAR;
	}
  }	 // pour toutes les classes

 /*
  for (c=0;c<NbClass;c++)
  {
  	printf("\n classe initiale=%d",c);
	for (m0=0;m0<NbMes[c];m0++)
	{
		
		for (par=0;par<NPAR;par++) 
		{
			printf("\n polit[%d][%d][%d]=%d",c,m0,par,politemp[c][m0*NPAR+par]);
			printf("\n polst[%d][%d][%d]=%d",c,m0,par,polstemp[c][m0*NPAR+par]);
		}

	
	}
  }	 // pour toutes les classes

  */
 }   // clalb_all

/*******************************************************************/



int CPolyMc::interne(long k1,long kpol,int c)
{ int par;
  byte bi,bs,b;
  for(par=0;par<NPAR;par++)
  {  bi=poli[c][kpol+par];bs=pols[c][kpol+par];
     b=politemp[c][k1+par];
     if ((b<bi)||(b>bs)) return 0;
     b=polstemp[c][k1+par];
     if ((b<bi)||(b>bs)) return 0;
  }//for par
 // printf("\n un interne c=%d k=%d",c,kpol);
  return 1;
}     // interne

void CPolyMc::etend(long k1,long kpol,int c)
{ int par;
  byte bi,bs,v;
  for(par=0;par<NPAR;par++)
  { 
	bi=poli[c][kpol+par];
	bs=pols[c][kpol+par];
	v=politemp[c][k1+par];
	if (bi<v) hypi[par]=bi; else hypi[par]=v;
	v=polstemp[c][k1+par];
	if (bs>v) hyps[par]=bs; else hyps[par]=v;
  }//for par
}

int CPolyMc::externe(long k0,int cl)
{ int par;
  byte b,bi,bs,i,s;
  for(par=0;par<NPAR;par++)
  { bi=hypi[par];bs=hyps[par];
    b=politemp[cl][k0+par];
    if (b>bi) i=b; else i=bi;
    b=polstemp[cl][k0+par];
    if (b<bs) s=b; else s=bs;
    if (s<=i) return 1;
  }//for par
  return 0;
}

int CPolyMc::externfusion(long k0,int cl)
{ int par;
  byte b,bi,bs,i,s;
  for(par=0;par<NPAR;par++)
  { bi=hypi[par];bs=hyps[par];
    b=poli[cl][k0+par];
    if (b>bi) i=b; else i=bi;
    b=pols[cl][k0+par];
    if (b<bs) s=b; else s=bs;
    if (s<=i) return 1;
  }//for par
  return 0;
}

int CPolyMc::libre(int clas)
{ 
/*	long k0;
	int clop;  
	for (clop=0;clop<NbClass;clop++)
	{
	 if (clop!=clas)
	 {
		for(k0=0;k0<kmax[clop];k0+=NPAR)
		{ 
			if (!externe(k0,clop)) return 0;
		}//for k0
	 } // si cl
	} // for cls
  return 1;*/
	//////////////////////////////////////////////////////////
	////////////////ATTENTION MODIF///////////////////////////
	//////////////////////////////////////////////////////////
	long k0;
	int clop;  


	for (clop=clas+1;clop<NbClass;clop++)
	{
		for(k0=0;k0<kmax[clop];k0+=NPAR)
		{ 
			if (!externe(k0,clop)) return 0;
		}//for k0
	} // for cls
	
	for (clop=0;clop<clas;clop++)
	{
		for(k0=0;k0<Npol[clop]*NPAR;k0+=NPAR)
		{ 
			if (!externfusion(k0,clop)) return 0;
		}//for k0
	} // for cls

  return 1;
}

void CPolyMc::remplace(long kpol,int c)
{ int par;
  for(par=0;par<NPAR;par++)
  {  poli[c][kpol+par]=hypi[par];
     pols[c][kpol+par]=hyps[par];
  }//for par
}

void CPolyMc::ajoute(long k1,int c)
{ int par;
  long knpol=(long)Npol[c]*(long)NPAR;
  for(par=0;par<NPAR;par++)
  {  poli[c][knpol+par]=politemp[c][k1+par];
     pols[c][knpol+par]=polstemp[c][k1+par];
  }//for par
  Npol[c]++;
}

void CPolyMc::calch(long k1,int clas)
{ long kpol,kmaxp;
  //int par;
  //printf("\r%5d",k1/NPAR);
  kmaxp=(long)Npol[clas]*(long)NPAR;
  for(kpol=0;kpol<kmaxp;kpol+=NPAR) {if (interne(k1,kpol,clas)) return;}
  for(kpol=0;kpol<kmaxp;kpol+=NPAR)
  { etend(k1,kpol,clas);// fusion dans Hypi/s
    if (libre(clas)) {remplace(kpol,clas);return;}
  }//for kpol
  ajoute(k1,clas); // change NPOL
}

void CPolyMc::calch_all()
{ 	
	int clas,par;//,m0,c;
	long k1;
  	for (clas=0;clas<NbClass;clas++)
	{
		Npol[clas]=0;
		for (par=0;par<NPAR;par++) {hypi[par]=0;hyps[par]=255;} 
		for(k1=0;k1<kmax[clas];k1+=NPAR) calch(k1,clas);

	}

/*	 
  for (c=0;c<NbClass;c++)
  {
  	//printf("\n classe initiale=%d",c);
//	printf("\n nombre de polytopes=%d",Npol[c]);

	for (m0=0;m0<Npol[c];m0++)
	{
	
		for (par=0;par<NPAR;par++) 
		{
			printf("\n polit[%d][%d][%d]=%d",c,m0,par,poli[c][m0*NPAR+par]);
			printf("\n polst[%d][%d][%d]=%d",c,m0,par,pols[c][m0*NPAR+par]);
		}

	
	}
  }	 // pour toutes les classes
*/
 


}


void CPolyMc::polytope()
{ 
  NPOL=0;
  MessageBox(NULL,"Attention, la fusion a �t� modifi�e !\nLes bornes aussi!!", "WARNING",MB_ICONEXCLAMATION);
  //printf("\n bornes\n");
  calcb_all();
  //if (BORD) printf("\n Attention points sur le bord 0 ou 255 !");
  //printf("\n fusion\n");
//  poli=va0;
//  pols=va1;
  //sauve_pol_avant("testavant.pmc");
  calch_all();
}
  
  

BOOL CPolyMc::Verif_redondance()
{
	long m,mo,el;
	int c,co,par,nb=0;//,par;
	long kmes,kmeso,kmesel;
	CString ch="",ch0;
	BOOL Verif,VerifTot=FALSE;
	for (c=0;c<NbClass;c++)
	{
		kmes=0;
		for(m=0;m<NbMes[c];m++)
		{
			for(co=c+1;co<NbClass;co++)
			{
				kmeso=0;
				for (mo=0;mo<NbMes[co];)
				{
					Verif=TRUE;
					for(par=0;par<NPAR;par++)
					{
						Verif&=(Va[c][kmes+par]==Va[co][kmeso+par]);
					}
					if(Verif) 
					{
						ch0.Format("classe %d et classe %d \n",c,co);
						ch+=ch0;
						nb++;
						///////////// elimination des redondances (on garde la classe la plus proche de 0)
						kmesel=kmeso;
						NbMes[co]--;
						kmax[co]-=NPAR;
						for(el=mo;el<NbMes[co];el++)
						{
							for(par=0;par<NPAR;par++)
							{
								Va[co][kmesel+par]=Va[co][kmesel+NPAR+par];
							}
							kmesel+=NPAR;
						}
						
						///////////////////
					}
					else {kmeso+=NPAR;mo++;}
					VerifTot|=Verif;	
				}// pour mo
			}//pour c0
			kmes+=NPAR;
		}//pour m
	}	 // pour toutes les classes
	
	ch0.Format("nb redondance =%d",nb);
	ch+=ch0;
	MessageBox(NULL,ch," Redondances!",MB_OK);
//	VerifSauvMMC();
	return VerifTot;
}
/*
void CPolyMc::VerifSauvMMC()
{	///sans redondance

	CFile fi;
	long c;
	long NbMesures;

	if(fi.Open("C:\\users\\pierreg\\testauto.mmc",CFile::modeWrite|CFile::modeCreate))
	{
		fi.Write(&NbClass,1);
		fi.Write(&NPAR,1);

		for (c=0;c<NbClass;c++)
		{
			NbMesures=(long)NbMes[c];
			fi.Write(&NbMesures,4);

			kmax[c]=(long)NbMes[c]*(long)NPAR;

			fi.Write(Va[c],kmax[c]);
		}//c

		fi.Close();;


	}
}


int CPolyMc::Verif_recouvr()
{
	int c,clop,k,kop,par,compt=0;
	BOOL verif=FALSE; //verifie si externe
	byte b,bi,bs,i,s;

	CString ch;

	for (c=0;c<NbClass;c++)
	{ 
		for(k=0;k<Npol[c];k+=NPAR)
		{
			for(clop=c+1;clop<NbClass;clop++)
			{
				for(kop=0;kop<Npol[clop];kop+=NPAR)
				{
					verif=FALSE;
					 for(par=0;par<NPAR;par++)
					  { bi=poli[clop][kop+par];bs=pols[clop][kop+par];
						b=poli[c][k+par];
						if (b>bi) i=b; else i=bi;
						b=pols[c][k+par];
						if (b<bs) s=b; else s=bs;
						if (s<=i) {verif = TRUE;break;}
					  }//for par
					 if(!verif) {compt++;ch.Format("classe %d poly n� %d et classe %d poly n� %d",c,k/NPAR,clop,kop/NPAR);MessageBox(NULL,ch,"",MB_OK);}
				}//kop
			}//clop
		}//k
	}//c
	return compt;
}

int CPolyMc::Verifrecouvrtemp()
{
	long k,kop;
	int c,clop,par,compt=0;
	BOOL verif=FALSE; //verifie si externe
	byte b,bi,bs,i,s;

	CString ch;

	for (c=0;c<NbClass;c++)
	{ 
		for(k=0;k<NbMes[c];k+=NPAR)
		{
			for(clop=c+1;clop<NbClass;clop++)
			{
				for(kop=0;kop<NbMes[clop];kop+=NPAR)
				{
					verif=FALSE;
					 for(par=0;par<NPAR;par++)
					  { bi=politemp[clop][kop+par];bs=polstemp[clop][kop+par];
						b=politemp[c][k+par];
						if (b>bi) i=b; else i=bi;
						b=polstemp[c][k+par];
						if (b<bs) s=b; else s=bs;
						if (s<=i) {verif = TRUE;break;}
					  }//for par
					 if(!verif) {compt++;}//ch.Format("classe %d poly n� %d et classe %d poly n� %d",c,k/NPAR,clop,kop/NPAR);MessageBox(NULL,ch,"",MB_OK);}
				}//kop
			}//clop
		}//k
	}//c
	ch.Format("nb recouvrement temp =%d",compt);
	MessageBox(NULL,ch,"recouvrement temp",MB_OK);
	return compt;
}

*/