#include "stdafx.h"
#include "tacq.h"
//#include "carte.h"
#include "mainfrm.h"
// ******************************************************************************************
// ******************************************************************************************
// LA CLASSE CTraitAcqWnd EST FAITE POUR ETRE MODIFIEE
//**********************  CLASSE DERIVEE de CAcqWnd********************
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!
// déclaration des messages!
// !!!!!!!!!!!!!!!!!!!!!!!!!
BEGIN_MESSAGE_MAP(CTraitAcqWnd, CAcqImageWnd)
	// fonctions déclarées dans CImage  
	//{{AFX_MSG_MAP(CTraitAcqWnd)
	ON_COMMAND(ID_ACQURIR_SQUENCE_3AQUISITIONSQUENCECALCUL, OnAcqurirSquence3aquisitionsquencecalcul)
	ON_UPDATE_COMMAND_UI(ID_ACQURIR_SQUENCE_3AQUISITIONSQUENCECALCUL, OnUpdateAcqurirSquence3aquisitionsquencecalcul)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()                            


// contructeur vide                   
CTraitAcqWnd::CTraitAcqWnd(int initLeType,NBCOLOR initColor = NB)
						:CAcqImageWnd(initLeType,initColor)  
{
}
						
CTraitAcqWnd::~CTraitAcqWnd()
{                                
// rajoutes ici
}  

//******************************************************************** 
// exemple de traitement vide !!
void CTraitAcqWnd:: OnAutreTraitement()
{                             
  SetCursor(LoadCursor(0,IDC_WAIT));    // souris sablié  
  // position du traitement 
  SetCursor(LoadCursor(0,IDC_ARROW));  // souris normal
  Invalidate(FALSE); // pour lancer le message paint de la fenêtre
}
//*********************************************************************  
	
// rajoutes tes traitements après et déclares les dans le fichier tacq.h 
// modifies le menu avec APP studio
// n'oublies pas de déclarer les nouveaux messages 


void CTraitAcqWnd::OnAcqurirSquence3aquisitionsquencecalcul() 
{
	DWORD time;
	char ch[10];
	int i,i1,i2,j1,j2;
	DWORD start_adresse1;   
	DWORD start_adresse2;
	MODCNF  *mod;
	BYTE  * Buf,*Buf1;
	Buf1=new BYTE[328*344];
	unsigned __int16 *cufen,**llufen;
	cufen=new unsigned __int16[328*344];
	llufen=new unsigned __int16*[NbMaxImageSeq];
	for(i=0;i<NbMaxImageSeq;i++) llufen[i]=new unsigned __int16[328*344];

	short fram = Return_frame(); 
	short fram2 = Return_frame2();	

    mod =  Return_mod();

	i1=debut_i;
	i2=fin_i;
	j1=debut_j;
	j2=fin_j;


	MessageBox("GO","lancement acquisition");

	icp_acqbits(mod, fram, ICP_SNAP);
	start_adresse1 = icp_get_acq_start_addr(mod);
//	icp_wait_acq (mod, fram);  
	
	icp_acqbits(mod, fram2, ICP_SNAP);
	start_adresse2 = icp_get_acq_start_addr(mod);
//	icp_wait_acq (mod, fram2);  
	
	// initialisation du temps 		
	itx_delta_ms(mod); 	
	
///////////////////////////***********************************///////////////////////////
// fram 2   PONG  
	icp_put_acq_start_addr(mod, start_adresse2); // Mise en place de l'adresse d'écriture
	write_reg(mod, ICP_ACQMD, 2); // declenchement de l'acquisition	      
///////////////////////////***********************************///////////////////////////


	for (i = 1;i<NbMaxImageSeq;i++)
	{
		Buf  = (BYTE  *) Bpile[i-1];
///////////////////////////***********************************///////////////////////////
		//attente fin pong
		while (read_reg(mod, ICP_ACQMD)!=0);	
///////////////////////////***********************************///////////////////////////
		// fram 1   PING  
		
		// lancement PING
		icp_put_acq_start_addr(mod, start_adresse1); // Mise en place de l'adresse d'écriture
		write_reg(mod, ICP_ACQMD, 2); // declenchement de l'acquisition	      
		
		// transfert pong en RAM
		icp_bm_mode (mod, ICP_ENABLE);	    
		icp_read_area (mod, fram2,i1,j1,i2-i1+1,j2-j1+1,(DWORD  *)Buf);  
	
		// traitement sur PONG
	GradientRapide((unsigned __int32 *)Buf,328,344);
	FenetreGlissanteRapide(Buf,cufen,llufen[i-1],328,344,12,12);

		// attente fin ping
		while (read_reg(mod, ICP_ACQMD)!=0);   
		//icp_wait_acq (mod, am_fram1);  
		
		// WAIT 
	//	if (Temps!=0) Wait((long)Temps);
		Buf  = (BYTE  *) Bpile[i];
		
		
		// fram 2   PONG  
		icp_put_acq_start_addr(mod, start_adresse2); // Mise en place de l'adresse d'écriture
		write_reg(mod, ICP_ACQMD, 2); // declenchement de l'acquisition	      
			    
		icp_bm_mode (mod, ICP_ENABLE);
		icp_read_area (mod, fram,i1,j1,i2-i1+1,j2-j1+1,(DWORD  *)Buf); 

		// traitement sur PING
	GradientRapide((unsigned __int32 *)Buf,328,344);
	FenetreGlissanteRapide(Buf,cufen,llufen[i],328,344,12,12);



	//	while (read_reg(mod, ICP_ACQMD)!=0);
		//icp_wait_acq (mod, am_fram2); 
		// WAIT 
	//	if (Temps!=0) Wait((long)Temps);
		i++;
		
	}	

time=(DWORD)itx_delta_ms(mod);	
sprintf(ch,"Image/.1s=%2.3f",(float)100*((float)NbMaxImageSeq/time));
	MessageBox(ch, "Séquence");
//CTraitImageWnd *im;
//im=new CTraitImageWnd(296,296);

	for(i=0;i<NbMaxImageSeq;i++)
		for(int j=0;j<328*328;j++) 
		{
		Bpile[i][j]=(unsigned char)(llufen[i][j]/144);
		//im->Image[j]=(unsigned char) llufen[0][j]/144;	
		}

//im->Create("lum glis 12*12", WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,GetMDIFrame());
//	Invalidate(FALSE);
	delete llufen;
	delete cufen;
	delete Buf1;
}

void CTraitAcqWnd::OnUpdateAcqurirSquence3aquisitionsquencecalcul(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable(SEQUENCE_CREE);	
}

void CTraitAcqWnd::OnDcisionDecisionsurfloat() 
{
	// TODO: Add your command handler code here
	
}
