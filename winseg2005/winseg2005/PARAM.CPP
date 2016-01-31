// poly.cpp : implementation file
//

#include "stdafx.h"
//#include "winseg.h"
#include "param.h"
#include "image.h"
#include "ProgDlg.h"

//////////////////////////////

#define _CRT_SECURE_NO_DEPRECATE
/////////////////////////////////////////////////////////////////////////////
// CParam

unsigned __int8 tabs[0x10000];
unsigned __int8 tabcont[0x10000];

CParam::CParam()
{
	int l,c,vt,contras;
	unsigned long pt=0;
	for(pt=0;pt<=0xFFFF;pt++)
	{	l=pt/0x100;c=pt&0xFF;
		if (l>c) 
		{
			vt=l-c;
			contras=(unsigned char)((float)255*(l-c)/(l+c));
				
		}
		else
		{
			vt=c-l;
			if(l==c) contras=0;
			else contras=(unsigned char)((float)255*(c-l)/(l+c));
		}
		tabs[pt]=unsigned __int8(vt/2);
		tabcont[pt]=contras;
	}

}


CParam::~CParam()
{

}

BOOL CParam::LumPix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
int l,c;      
long Decale[1],p;

char ch[20];

Decode(chaine,Decale);

sprintf(ch,"Luminance pixel: %d",Decale[0]);
CProgressDlg Dlg;
Dlg.Create();
Dlg.SetWindowText(ch);
Dlg.SetRange(0,Hautl);
Dlg.SetStep(1);


 for (l=0;l<Hautl;l++)
 {
	p=Largl*l;
 	for (c=0;c<Largl;c++)
 	{
 		
 	if (p+Decale[0]>=0 && (p+Decale[0])<(Largl*Hautl)) mem[p]=imi[p+Decale[0]];
 	else mem[p]=0;
	p++;
 	}
	Dlg.StepIt();
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}

 } // l

Dlg.DestroyWindow();
return TRUE;
}// lumpix
 	
BOOL CParam::LumMoyGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2],temp,p,q,qot;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	char ch[50];

Decode(chaine,val);
sprintf(ch,"Lum moy fen glis: %d,%d",val[0],val[1]);
CProgressDlg Dlg;
Dlg.Create();
Dlg.SetWindowText(ch);
Dlg.SetRange((int)(val[0]-1)/2,2*Hautl-(int)val[0]/2+1);
Dlg.SetStep(1);

debutl=(int)(val[0]-1)/2;
finl=Hautl-(int)val[0]/2;
debutc=(int)(val[1]-1)/2;
finc=Largl-(int)val[1]/2;	
debutlc=-debutl;
finlc=(int)val[0]/2;
debutcc=-debutc;
fincc=(int)val[1]/2;

qot=val[0]*val[1];

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}

	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}// for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for l

	for (l=debutl;l<finl;l++)
	{
		p=l*Largl+debutc;
		for(c=debutc;c<finc;c++)
		{
			//p=Largl*l+c;
			temp=0;
			for(lc=debutlc;lc<=finlc;lc++)
			{
				q=p+debutcc+Largl*lc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					temp+=imi[q++];
				}//for cc
			}//for lc

			temp/=qot;
			mem[p++]=(byte)temp;
		}//for c
		Dlg.StepIt();
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	}//for l


Dlg.DestroyWindow();
return TRUE;
}//lummoygliss


BOOL CParam::LumMoyFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2],temp,p,q,qot;
	int l,c,lc,cc;
	char ch[40];

	Decode(chaine,val);

sprintf(ch,"Lum moy fen fix: %d,%d",val[0],val[1]);
CProgressDlg Dlg;
Dlg.Create();
Dlg.SetWindowText(ch);
Dlg.SetRange(0,Hautl-val[0]+2);
Dlg.SetStep(val[0]);

qot=val[0]*val[1];
	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}

	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}//for l for c

	}
	for (l=0;l<=Hautl-val[0];l+=val[0])
	{
		p=l*Largl;
		for(c=0;c<=Largl-val[1];c+=val[1])
		{
			temp=0;

			for(lc=0;lc<val[0];lc++)
			{
				q=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					temp+=imi[q++];
				}//for cc
			}//for lc

			temp/=qot;
			
			for(lc=0;lc<val[0];lc++)
			{
				q=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					mem[q++]=(byte)temp;
				}//for cc
			}//for lc
			
			p+=val[1];
		}//for c
		Dlg.StepIt();
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	}//for l

	Dlg.DestroyWindow();
	return TRUE;
}//lummoyfix

BOOL CParam::GradPix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{

int gra,l,c;      
long Decale[1],p;
char ch[50];

Decode(chaine,Decale);

sprintf(ch,"Gradient pixel: %d",Decale[0]);
CProgressDlg Dlg;
Dlg.Create();
Dlg.SetWindowText(ch);
Dlg.SetRange(0,Hautl);
Dlg.SetStep(1);



 for (l=0;l<Hautl;l++)
 {
	 p=Largl*l;
 	for (c=0;c<Largl;c++)
 	{
 	                                                    
 	
 	if (p+Decale[0]>=0 && (p+Decale[0]+Largl+1)<(Largl*Hautl)) 
 	{
 		gra=(int)abs(imi[p+Decale[0]]-imi[p+Decale[0]+Largl+1]) + (int)abs(imi[p+Decale[0]+1]-imi[p+Decale[0]+Largl]);
 	    if(gra<=255)   mem[p++]=(byte) gra;
		else mem[p++]=255;
 	}
 	else mem[p++]=0;
	
 	}//c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
 } // l
 
Dlg.DestroyWindow();
return TRUE;
}// gradpix

BOOL CParam::GradMoyGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{ 
	long val[2],temp,p,q,qot;
	int gra,l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;    
	byte * tampon;
	tampon= new byte[Largl*Hautl];

	Decode(chaine,val);
char ch[40];
sprintf(ch,"Grad moy fen glis: %d,%d",val[0],val[1]);
CProgressDlg Dlg;
Dlg.Create();
Dlg.SetWindowText(ch);
Dlg.SetRange((int)(val[0]-1)/2,2*Hautl-(int)val[0]/2+1);
Dlg.SetStep(1);

debutl=(int)(val[0]-1)/2;
finl=Hautl-(int)val[0]/2;
debutc=(int)(val[1]-1)/2;
finc=Largl-(int)val[1]/2;	
debutlc=-debutl;
finlc=(int)val[0]/2;
debutcc=-debutc;
fincc=(int)val[1]/2;

qot=val[0]*val[1];

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}


	for (l=0;l<Hautl-1;l++)
	{
		p=l*Largl;
		for (c=0;c<Largl-1;c++)
		{                                                  

			gra=(int)abs(imi[p]-imi[p+Largl+1]) + (int)abs(imi[p+1]-imi[p+Largl]);
 			if(gra<=255)   tampon[p++]=(byte) gra;
			else tampon[p++]=255;
		}//c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	} // l

	for (l=debutl;l<finl;l++)
	{
		p=debutc+l*Largl;
		for(c=debutc;c<finc;c++)
		{
			temp=0;

			for(lc=debutlc;lc<=finlc;lc++)
			{ 
				q=p+debutcc+Largl*lc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					temp+=tampon[q++];
				}//for cc
			}//for lc

			temp/=qot;
			mem[p++]=(byte)temp;
		
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();		
	}//for l
	
	Dlg.DestroyWindow();
	delete tampon;
	return TRUE;
}//gradmoyglis

BOOL CParam::GradMoyGlis2(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{ 
	long val[2],temp,p,q,qot;
	int gra,l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;    
	byte * tampon;
	tampon= new byte[Largl*Hautl];

	Decode(chaine,val);
char ch[40];
sprintf(ch,"Grad moy fen glis2: %d,%d",val[0],val[1]);
CProgressDlg Dlg;
Dlg.Create();
Dlg.SetWindowText(ch);
Dlg.SetRange((int)(val[0]-1)/2,2*Hautl-(int)val[0]/2+1);
Dlg.SetStep(1);

debutl=(int)(val[0]-1)/2;
finl=Hautl-(int)val[0]/2;
debutc=(int)(val[1]-1)/2;
finc=Largl-(int)val[1]/2;	
debutlc=-debutl;
finlc=(int)val[0]/2;
debutcc=-debutc;
fincc=(int)val[1]/2;

qot=val[0]*val[1];

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}


	for (l=0;l<Hautl;l++)
	{
		p=l*Largl;
		for (c=0;c<Largl-1;c++)
		{                                                  

			gra=(int)abs(imi[p]-imi[p+Largl+1])/2 + (int)abs(imi[p+1]-imi[p+Largl])/2;
 			tampon[p++]=(byte) gra;
		}//c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	} // l

	for (l=debutl;l<finl;l++)
	{
		p=debutc+l*Largl;
		for(c=debutc;c<finc;c++)
		{
			temp=0;

			for(lc=debutlc;lc<=finlc;lc++)
			{ 
				q=p+debutcc+Largl*lc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					temp+=tampon[q++];
				}//for cc
			}//for lc

			temp/=qot;
			mem[p++]=(byte)temp;
		
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();		
	}//for l
	
	Dlg.DestroyWindow();
	delete tampon;
	return TRUE;
}//gradmoyglis2

BOOL CParam::Gradientrapide(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{ 
	long val[2],qot;
	int stop,l;
	unsigned __int16 temp;
	unsigned __int16* 	cufen;
	unsigned __int16* 	lufen;
	cufen=new unsigned __int16[Largl*Hautl];
	lufen=new unsigned __int16[Largl*Hautl];
	Decode(chaine,val);

	if((val[0]<1)||(val[1]<1)||(val[0]>16)||(val[1]>16)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	if((val[0]%4)||(val[1]%4)) {MessageBox(NULL,"les dimensions des fenêtres ne sont pas multiples de 4!","Impossible",MB_OK|MB_ICONSTOP);return FALSE;}
	char ch[40];
	sprintf(ch,"Grad rapide glis: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,100);
	Dlg.SetStep(50);

	qot=val[0]*val[1]/2;
	int i;
	for(i=0;i<Largl*Hautl;i++)mem[i]=imi[i];

	GradientRap((unsigned __int32 *)mem,Largl,Hautl);
	Dlg.StepIt();
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	FenetreGlissanteRap(mem, cufen,lufen,Largl,Hautl,val[1],val[0]);
	
	stop=int(val[1]/2)+1+Largl*int(val[0]/2);
	l=stop;
	for(i=0;i<Largl*Hautl-stop;i++) 
	{
		if((temp=unsigned __int16(lufen[l++]/qot))>255) mem[i]=255;
		else mem[i]=(unsigned char)temp;
	}
	for(i=Largl*Hautl-stop;i<Largl*Hautl;i++) mem[i]=0;

	Dlg.StepIt();		

	
	Dlg.DestroyWindow();
	delete cufen;
	delete lufen;
	return TRUE;
}//gradientrapide

BOOL CParam::GradientrapideNonNorm(CString chaine, byte  *mem, byte   *imi, long Largl, long Hautl)
{ 
	long val[2],qot;
	int stop,l;
	unsigned __int16 temp;
	unsigned __int16* 	cufen;
	unsigned __int16* 	lufen;
	cufen=new unsigned __int16[Largl*Hautl];
	lufen=new unsigned __int16[Largl*Hautl];
	Decode(chaine,val);

	if((val[0]<1)||(val[1]<1)||(val[0]>16)||(val[1]>16)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	if((val[0]%4)||(val[1]%4)) {MessageBox(NULL,"les dimensions des fenêtres ne sont pas multiples de 4!","Impossible",MB_OK|MB_ICONSTOP);return FALSE;}
	char ch[40];
	sprintf(ch,"Grad rapide NonNorm: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,100);
	Dlg.SetStep(50);

	qot=val[0]*val[1]/4;///6.4
	int i;
	for(i=0;i<Largl*Hautl;i++)mem[i]=imi[i];

	GradientRap((unsigned __int32 *)mem,Largl,Hautl);
	Dlg.StepIt();
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	FenetreGlissanteRap(mem, cufen,lufen,Largl,Hautl,val[1],val[0]);
	
	stop=int(val[1]/2)+1+Largl*int(val[0]/2);
	l=stop;
	for(i=0;i<Largl*Hautl-stop;i++) 
	{
		if((temp=unsigned __int16(lufen[l++]/qot))>255) mem[i]=255;
		else mem[i]=(unsigned char)temp;
	}
	for(i=Largl*Hautl-stop;i<Largl*Hautl;i++) mem[i]=0;


	Dlg.StepIt();		

	
	Dlg.DestroyWindow();
	delete cufen;
	delete lufen;
	return TRUE;
}//gradientrapide

BOOL CParam::Lumirapide(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{ 
	long val[2],qot;
	int i,stop,l;
	unsigned __int16* 	cufen;
	unsigned __int16* 	lufen;
	cufen=new unsigned __int16[Largl*Hautl];
	lufen=new unsigned __int16[Largl*Hautl];
	Decode(chaine,val);

	if((val[0]<1)||(val[1]<1)||(val[0]>16)||(val[1]>16)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	if((val[0]%4)||(val[1]%4)) {MessageBox(NULL,"les dimensions des fenêtres ne sont pas multiples de 4!","Impossible",MB_OK|MB_ICONSTOP);return FALSE;}
	char ch[40];
	sprintf(ch,"Lum rapide glis: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,100);
	Dlg.SetStep(50);

	qot=val[0]*val[1];

	Dlg.StepIt();
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	FenetreGlissanteRap(imi, cufen,lufen,Largl,Hautl,val[1],val[0]);
	
	stop=int(val[1]/2)+Largl*int(val[0]/2);
	l=stop;
	for(i=0;i<Largl*Hautl-stop;i++) mem[i]=unsigned char(lufen[l++]/qot);
	for(i=Largl*Hautl-stop;i<Largl*Hautl;i++) mem[i]=0;

	Dlg.StepIt();		

	
	Dlg.DestroyWindow();
	delete cufen;
	delete lufen;
	return TRUE;
}//luminance rapide

void CParam::GradientRap(unsigned __int32 *lram,int large,int haute)
{
	int l,c,c4;
	unsigned long pth,ptl;
	unsigned __int32 ph,pl,gradph;
	unsigned __int8 newph,newpl,oldph=0,oldpl=0;
	c4=large/4;
	pth=0;ptl=pth+c4;
	for(l=0;l<(haute-1);l++)
	{ for(c=0;c<c4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			_asm
			{
				mov ebx,ph
				mov ecx,pl
				bswap ecx
				bswap ebx
				mov newph,bl
				mov newpl,cl
				xor eax,eax
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,tabs[eax]
				mov ax,cx
				mov cl,tabs[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,tabs[eax]
				mov ax,cx
				mov cl,tabs[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov bl,tabs[eax]
				mov ax,cx
				mov cl,tabs[eax]
				add bl,cl
				ror ebx,8
				ror ecx,8
				;***************
				xchg bl,cl
				mov ax,bx
				mov ah,oldph
				mov bl,tabs[eax]
				mov ax,cx
				mov ah,oldpl
				mov cl,tabs[eax]
				add bl,cl
				ror ebx,8
				;	ror ecx,8
				;***************
				mov al,newph
				mov oldph,al
				mov ah,newpl
				mov oldpl,ah
				bswap ebx
				mov gradph,ebx
				;xchg ah,al	;12->21
				;bswap eax	;1234->4321
			}
			lram[pth]=gradph;
			pth++;ptl++;
		}//for c
	}//for l
	
}

void CParam::FenetreGlissanteRap(unsigned char* bram,unsigned __int16* cufen,
											unsigned __int16* lufen,int large,int haute,int LF,int HF)
{
	unsigned __int32 p,pth,ptl,pb;
	unsigned __int32 ph,pl;
	unsigned __int16 som;
	unsigned __int32 *lram;
	unsigned __int32 som2,*dcufen,*dlufen;
	int l,c,c4,lf4,c2;//44 msc

	lram=(unsigned __int32 *)bram;
	////////colonnes
	c4=large/4;lf4=LF/4;
	lram=(unsigned __int32 *)bram;
	pb=0;ptl=0;pth=ptl+lf4;
	for(l=0;l<haute;l++)
	{	som=0;p=l*large;
		for(c=0;c<LF;c++) 
		{som=som+bram[p++];cufen[pb++]=som;}
		for(c=lf4;c<c4;c++)
		{	ph=lram[pth];pl=lram[ptl];
			som=som+unsigned __int16((ph&0xff)-(pl&0xff));
			cufen[pb++]=som;
			som=som+unsigned __int16(((ph&0xff00)>>8)-((pl&0xff00)>>8));
			cufen[pb++]=som;
			som=som+unsigned __int16(((ph&0xff0000)>>16)-((pl&0xff0000)>>16));
			cufen[pb++]=som;
			som=som+unsigned __int16((ph>>24)-(pl>>24));
			cufen[pb++]=som;
			pth++;ptl++;
		}//for c
		pth+=lf4;ptl+=lf4;
	}//for l

//////////////////////ligne
	dcufen=(unsigned __int32 *)cufen;
	dlufen=(unsigned __int32 *)lufen;
	c2=large/2;
	for(c=0;c<c2;c++)
	{	pth=c;som2=0;ptl=c;
		for(l=0;l<HF;l++) 
		{	som2=som2+dcufen[pth];//attention pas de carry
			dlufen[pth]=som2;
			pth+=c2;
		}//for l
		
		for(l=HF;l<haute;l++)
		{	som2=som2+(dcufen[pth]-dcufen[ptl]);
			dlufen[pth]=som2;
			pth+=c2;ptl+=c2;
		}//for l
	}//for c
}


BOOL CParam::GradMoyFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
long val[2],temp,p,q,qot;
int gra,l,c,lc,cc;    
byte * tampon;
tampon= new byte[Largl*Hautl];
	
Decode(chaine,val);

char ch[40];
sprintf(ch,"Grad moy fen fix: %d,%d",val[0],val[1]);
CProgressDlg Dlg;
Dlg.Create();
Dlg.SetWindowText(ch);
Dlg.SetRange(0,2*Hautl-val[0]+2);
Dlg.SetStep(1);

qot=val[0]*val[1];

if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}

 for (l=0;l<Hautl;l++)
 {
	 p=Largl*l;
 	for (c=0;c<Largl-1;c++)
 	{   
 		
 			gra=(int)abs(imi[p]-imi[p+Largl+1]) + (int)abs(imi[p+1]-imi[p+Largl]);
 			if(gra<=255)   tampon[p++]=(byte) gra;
			else tampon[p++]=255;		
 	}//c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
 } // l
 Dlg.SetStep(val[0]);

for (l=0;l<=Hautl-val[0];l+=val[0])
{
	p=l*Largl;
	for(c=0;c<=Largl-val[1];c+=val[1])
	{
		temp=0;
		
		for(lc=0;lc<val[0];lc++)
		{
			q=p+Largl*lc;
			for(cc=0;cc<val[1];cc++)
			{
				temp+=tampon[q++];
			}//for cc
		}//for lc

		temp/=qot;
		
		for(lc=0;lc<val[0];lc++)
		{
			q=p+Largl*lc;
			for(cc=0;cc<val[1];cc++)
			{
				mem[q++]=(byte) temp;
			}//for cc
		}//for lc
		p+=val[1];
	}//for c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
}//for l

Dlg.DestroyWindow();
delete tampon;
return TRUE;
}//GradMoyFix

BOOL CParam::LumMoyDecLigFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
long val[2],p,q,qot;
long *temp;
int l,c,lo,co;


	Decode(chaine,val);

char ch[50];
sprintf(ch,"Lum moy décal lign fen fix: %d,%d",val[0],val[1]);
CProgressDlg Dlg;
Dlg.Create();
Dlg.SetWindowText(ch);
Dlg.SetRange(0,2*Hautl-val[0]+2);
Dlg.SetStep(1);

qot=val[0]*val[1];
	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Lum moy décal lign fen fix",MB_OK|MB_ICONSTOP);return FALSE;}
	
	
	temp=new long [(int) Largl];
	
	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}// for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l

Dlg.SetStep(val[0]);

	for (l=0;l<=Hautl-val[0];l+=val[0])
	{
		p=l*Largl;
		for(c=0;c<=Largl-val[1];c+=val[1])
		{
			for(lo=0;lo<val[0];lo++)
			{
				q=p+Largl*lo;
				for(co=0;co<val[1];co++)
				{
					if(l==0) {mem[q++]=0;}
					 else {mem[q++]=(byte)temp[c];}
				}//for co
			}//for lo
			
			temp[c]=0; 
			for(lo=0;lo<val[0];lo++)
			{
				q=p+Largl*lo;
				for(co=0;co<val[1];co++)
				{
					temp[c]+=imi[q++];
				}//for co
			}//for lo

			temp[c]/=qot;
			p+=val[1];
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l
	delete temp;
	
	Dlg.DestroyWindow();
	return TRUE;
}//LumMoyDecLigFix



BOOL CParam::Entropie(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[3],ligne,colo,li,co,p=0,k=0,qot,q;
	double d;
	unsigned char *pilnpi;
	word  *cpt;
	word *som;
	int i;
	Decode(chaine,val);
	
char ch[50];
sprintf(ch,"Entropie fix: %d,%d,%d",val[0],val[1],val[2]);
CProgressDlg Dlg;
Dlg.Create();
Dlg.SetWindowText(ch);
Dlg.SetRange(0,2*Hautl-val[0]+2);
Dlg.SetStep(1);


qot=val[0]*val[1];

	if((val[0]<1)||(val[1]<0)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	if(val[2]<=0) {MessageBox(NULL,"erreur dans les parametres","Entropie glis",MB_OK|MB_ICONSTOP);return FALSE;} 

	pilnpi=new unsigned char[val[0]*val[1]+1];
	cpt=new word[val[2]];
	som=new word[((int)Hautl/val[0]+1)*((int)Largl/val[1]+1)];
  
	for(ligne=0;ligne<Hautl;ligne++)
	{	
		p=Largl*ligne;
		for(colo=0;colo<Largl;colo++)
		{
			mem[p]=imi[p++];
		}//for  colo
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for ligne

    pilnpi[0]=0;
    for ( i=1;i<val[0]*val[1]+1;i++)  
    {
    	d=( (253*exp(1.0)*i/(qot))*log(double(i)/(qot)));// variation entre 0 et 255
    	pilnpi[i]=(byte)(fabs(d));
    }
Dlg.SetStep(val[0]);
 for (ligne=0;ligne<=Hautl-val[0];ligne+=val[0])
 {  
   p=Largl*ligne;
   for (colo=0;colo<=Largl-val[1];colo+=val[1])
   {
  	
  	for (i=0;i<val[2];i++) cpt[i]=0; 

  	som[k]=0;
  	
  	for(li=0;li<val[0];li++)
	{
		q=p+Largl*li;
		for(co=0;co<val[1];co++)
		{
			  
			cpt[imi[q++]*val[2]/256]++;
        } // co
	} // li
	for (i=0;i<val[2];i++)
	{
	som[k]=som[k]+pilnpi[cpt[i]];
//	printf("\n pilnpi=%d",pilnpi[cpt[i]]);
	}
//	printf("\n somme=%u",somme[k]);
	som[k]=(word) (som[k]/val[2]);

	for(li=0;li<val[0];li++)
	{
		q=p+Largl*li;
		for(co=0;co<val[1];co++)
		{
			mem[q++]=(unsigned char) som[k];
		}
	}
    k++;
	p+=val[1];
   } // colo
   if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
   Dlg.StepIt();
  } // ligne	
 
 Dlg.DestroyWindow();
 delete pilnpi;
 delete cpt;
 delete som;
 return TRUE;
};// Entrpoie


BOOL CParam::EntropieGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{	
	long val[3],ligne,colo,li,co,p=0,k=0,q,qot;
	int debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	double d;
	unsigned char *pilnpi;
	word  *cpt;
	word som;
	int i;
	Decode(chaine,val);
	char ch[50];
	sprintf(ch,"Entropie glis: %d,%d,%d",val[0],val[1],val[2]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,2*Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	qot=val[0]*val[1];

	if((val[0]<1)||(val[1]<0)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	if(val[2]<=0) {MessageBox(NULL,"erreur dans les parametres","Entropie glis",MB_OK|MB_ICONSTOP);return FALSE;} 

	pilnpi=new unsigned char[val[0]*val[1]+1];
	cpt=new word[val[2]];
	 
	for(ligne=0;ligne<Hautl;ligne++)
	{	
		p=Largl*ligne;
		for(colo=0;colo<Largl;colo++)
		{
			mem[p]=imi[p++];
		}// for colo
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for ligne

    pilnpi[0]=0;

    for (int i=1;i<val[0]*val[1]+1;i++)  
    {
    	d=( (253*exp(1.0)*i/(qot))*log(double(i)/(qot)));// variation entre 0 et 255
    	pilnpi[i]=(byte)(fabs(d));
    }

	for (ligne=debutl;ligne<finl;ligne++)
	{
		p=Largl*ligne+debutc;
		for(colo=debutc;colo<finc;colo++)
		{
			for (i=0;i<val[2];i++) cpt[i]=0;

			
			for(li=debutlc;li<=finlc;li++)
			{
				q=p+Largl*li+debutcc;
				for(co=debutcc;co<=fincc;co++)
				{
					cpt[imi[q++]*val[2]/256]++;
				}//for co
			}//for li
			som=0;
			for (i=0;i<val[2];i++)
			{
				som+=pilnpi[cpt[i]];
			}
			som=(word)(som/val[2]);
			mem[p++]=(unsigned char)som;
			
		}//for colo
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for ligne

Dlg.DestroyWindow();
delete pilnpi;
delete cpt;
return TRUE;
};//EntropieGlis

BOOL CParam::ContrastLoc(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2];
	long p,q;
	int l,c,lc,cc;
	unsigned char temp,max,min;
	double cont;
	Decode(chaine,val);
	char ch[50];
	sprintf(ch,"Contrast fix: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,Hautl-val[0]+2);
	Dlg.SetStep(val[0]);

	if((val[0]<1)||(val[1]<0)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}


	for (l=0;l<=Hautl-val[0];l+=val[0])
	{
		p=l*Largl;
		for(c=0;c<=Largl-val[1];c+=val[1])
		{
			max=0;			
			min=255;
			for(lc=0;lc<val[0];lc++)
			{
				q=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					temp=imi[q++];
					if(temp>=max) max=temp;
					if (temp<=min) min=temp;
				}//for cc
			}//for lc
			if ((max+min)==0) cont=0;
			else cont=255*((double)(max-min)/(max+min));
			

			for(lc=0;lc<val[0];lc++)
			{
				q=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					mem[q++]=(byte) cont;
				}//for cc
			}//for lc
			p+=val[1];
		}//for c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	return TRUE;
}//contrastLoc
BOOL CParam::Max_Min(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2];
	long p,q;
	int l,c,lc,cc;
	unsigned char temp,max,min;
	
	Decode(chaine,val);
char ch[50];
sprintf(ch,"Max_min fix: %d,%d",val[0],val[1]);
CProgressDlg Dlg;
Dlg.Create();
Dlg.SetWindowText(ch);
Dlg.SetRange(0,Hautl-val[0]+2);
Dlg.SetStep(val[0]);
	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}


	for (l=0;l<=Hautl-val[0];l+=val[0])
	{
		p=l*Largl;
		for(c=0;c<=Largl-val[1];c+=val[1])
		{
			max=0;			
			min=255;
			for(lc=0;lc<val[0];lc++)
			{
				q=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					temp=imi[q++];
					if(temp>=max) max=temp;
					if (temp<=min) min=temp;
				}//for cc
			}//for lc
			 
			for(lc=0;lc<val[0];lc++)
			{
				q=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					mem[q++]=(byte) (max-min);
				}//for cc
			}//for lc
			p+=val[1];
		}//for c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	return TRUE;
}//Max-min

BOOL CParam::ContrastGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2],temp,p,q,max,min;
	double cont;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;

	Decode(chaine,val);
	char ch[50];
	sprintf(ch,"Contrast glis: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,2*Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}// for c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for l

	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			temp=0;
			max=0;
			min=255;

			for(lc=debutlc;lc<=finlc;lc++)
			{
				q=p+Largl*lc+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					temp=imi[q++];
					if (temp>=max) max=temp;
					if (temp<=min) min=temp;
				}//for cc
			}//for lc
			if((max==0)&&(min==0)) cont=0;
			else cont=255*((double)(max-min)/(max+min));
			mem[p++]=(byte) cont;
		}//for c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	return TRUE;
}//Contrastgliss

BOOL CParam::ContrastrapideGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2],l,stop;
	unsigned char *tampon;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Contrast rapid glis: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,100);
	Dlg.SetStep(50);

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	tampon=new unsigned char[Largl*Hautl];

	MinimumGlissant(mem,imi,Largl,Hautl,val[1],val[0]);
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	MaximumGlissant(tampon,imi,Largl,Hautl,val[1],val[0]);
	stop=int(val[1]/2)+Largl*int(val[0]/2);
	l=stop;
	Dlg.StepIt();
	int i;
	for( i=0;i<Largl*Hautl-stop;i++) mem[i]=tabcont[0x100*tampon[l]+mem[l++]];
	for(i=Largl*Hautl-stop;i<Largl*Hautl;i++) mem[i]=0;
	Dlg.DestroyWindow();
	delete tampon;
	return TRUE;
}//Contrastrapidgliss



BOOL CParam::Max_MinGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2],temp,p,q,max,min;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;

	Decode(chaine,val);
	char ch[50];
	sprintf(ch,"Max_min glis: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,2*Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);


	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}// for c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for l

	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			temp=0;
			max=0;
			min=255;

			for(lc=debutlc;lc<=finlc;lc++)
			{
				q=p+Largl*lc+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					temp=imi[q++];
					if (temp>=max) max=temp;
					if (temp<=min) min=temp;
				}//for cc
			}//for lc
			 mem[p++]=(byte) (max-min);
		}//for c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	return TRUE;
}//Max-min gliss

BOOL CParam::MedianGlissant(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[3],p,q,r;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	unsigned char *place;
	Decode(chaine,val);

	char ch[50];
	sprintf(ch,"Median glis: %d,%d,%d",val[0],val[1],val[2]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,2*Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	if((val[0]<1)||(val[1]<1)||(val[2]<0)||(val[2]>(val[0]*val[1]))) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	place=new unsigned char[val[0]*val[1]];

	for(l=0;l<Hautl;l++)
	{	
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}// for c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for l

	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			q=0;

			for(lc=debutlc;lc<=finlc;lc++)
			{
				r=p+lc*Largl+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					place[q++]=imi[r++];
				}//for cc
			}//for lc
			qsort(place,(size_t) (val[1]*val[0]),sizeof(unsigned char),comparaison);
			mem[p++]=(byte) place[val[2]];
		}//for c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	delete place;
	return TRUE;
}//Median gliss


BOOL CParam::MedianFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[3],p,q,r;
	int l,c,lc,cc;
	unsigned char *place;
	Decode(chaine,val);
char ch[50];
sprintf(ch,"Median fix: %d,%d",val[0],val[1]);
CProgressDlg Dlg;
Dlg.Create();
Dlg.SetWindowText(ch);
Dlg.SetRange(0,Hautl-val[0]+2);
Dlg.SetStep(val[0]);
	if((val[0]<1)||(val[1]<1)||(val[2]<0)||(val[2]>=(val[0]*val[1]))) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	place=new unsigned char[val[0]*val[1]];

	for (l=0;l<=Hautl-val[0];l+=val[0])
	{
		p=l*Largl;
		for(c=0;c<=Largl-val[1];c+=val[1])
		{
			q=0;			
			for(lc=0;lc<val[0];lc++)
			{
				r=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					place[q++]=imi[r++];
				}//for cc
			}//for lc
			 
			qsort(place,(size_t) (val[1]*val[0]),sizeof(unsigned char),comparaison);

			for(lc=0;lc<val[0];lc++)
			{
				r=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					mem[r++]=place[val[2]];;
				}//for cc
			}//for lc
		p+=val[1];	
		}//for c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	delete place;
	return TRUE;
}//medianfix

BOOL CParam::SobelFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{	
	long val[2],temp,p,q,qot;
	int l,c,lc,cc;
	byte * tampon;
	tampon= new byte[Largl*Hautl];


	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Sobel moy fen fix: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,Hautl-val[0]+2);
	Dlg.SetStep(val[0]);

	qot=val[0]*val[1];

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	
	for(long i=Largl;i<(Largl-1)*Hautl;i++)
	{
		tampon[i]= (unsigned char)(labs(2*(long)imi[i-1]+(long)imi[i+Largl-1]+
				(long)imi[i-Largl-1]-2*(long)imi[i+1]-(long)imi[i+Largl+1]-(long)imi[i-Largl+1])+
  				 labs((long)imi[i+Largl-1]+2*(long)imi[i+Largl]+(long)imi[i+Largl+1]-(long)imi[i-Largl-1]-2*(long)imi[i-Largl]-(long)imi[i-Largl+1]));

	}

	for (l=0;l<=Hautl-val[0];l+=val[0])
	{
		p=l*Largl;
		for(c=0;c<=Largl-val[1];c+=val[1])
		{
			temp=0;
			
			for(lc=0;lc<val[0];lc++)
			{
				q=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					temp+=tampon[q++];
				}//for cc
			}//for lc

			temp=temp/(qot);
			
			for(lc=0;lc<val[0];lc++)
			{
				q=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					mem[q++]=(byte)temp;
				}//for cc
			}//for lc
			p+=val[1];
		}//for c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	delete tampon;
	return TRUE;
}//sobel fix

BOOL CParam::SobelGlissant(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{	
	long val[2],temp,p,q,qot;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	byte * tampon;
	tampon= new byte[Largl*Hautl];	
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Sobel moy fen glis: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);


	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	qot=val[0]*val[1];
	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}

	for(long i=Largl;i<(Largl-1)*Hautl;i++)
	{
		tampon[i]= (unsigned char)(labs(2*(long)imi[i-1]+(long)imi[i+Largl-1]+
				(long)imi[i-Largl-1]-2*(long)imi[i+1]-(long)imi[i+Largl+1]-(long)imi[i-Largl+1])+
  				 labs((long)imi[i+Largl-1]+2*(long)imi[i+Largl]+(long)imi[i+Largl+1]-(long)imi[i-Largl-1]-2*(long)imi[i-Largl]-(long)imi[i-Largl+1]));
	}

	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			temp=0;

			for(lc=debutlc;lc<=finlc;lc++)
			{
				q=p+Largl*lc+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					temp+=tampon[q++];
				}//for cc
			}//for lc

			temp/=qot;
			mem[p++]=(byte)temp;
		}//for c
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	}//for l
Dlg.DestroyWindow();
delete tampon;
return TRUE;
}//sobelglissant

BOOL CParam::ImageEnregistre(byte *mem, long Largl, long Hautl)
{
	CFile Fichier; // variable fichier BITMAP	  
	CFileDialog * FileData;	
	CProgressDlg Dlg;

	// ouverture de la boite de dialogue d'ouverture d'un fichier
	FileData = new CFileDialog (TRUE,"*.bmp","c:\\users\\pierreg\\bmp\\*.bmp",OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,"BITMAP (*.bmp)|*.bmp|BRUTE (*.mem)|*.mem||",NULL);
		
	if (FileData->DoModal()==IDOK)
	{
		if (Fichier.Open(FileData->GetPathName(), CFile::modeRead))
  		 { 
			
  	 		Fichier.Close(); 
			CImageWnd	* pImageWnd;
			pImageWnd = new CImageWnd(FileData->GetPathName(), Largl, Hautl, BMP);    
			//pImageWnd->Create(FileData->GetFileTitle(), WS_CHILD | WS_VISIBLE | WS_OVERLAPPEDWINDOW,NULL);
			if ((pImageWnd->Larg!=Largl)||(pImageWnd->Haut!=Hautl)) 
			{
				MessageBox(NULL,"Erreur dimensions !","Message",MB_OK|MB_ICONEXCLAMATION);
				delete pImageWnd;
				ImageEnregistre(mem,Largl,Hautl);
				return FALSE;
			}
			else 
			{
				Dlg.Create();
				Dlg.SetWindowText("Chargement de "+FileData->GetFileTitle()+".bmp");
				Dlg.SetRange(0,Hautl);
				Dlg.SetStep(1);
				for(long j=0;j<Hautl;j++)
				{
					for(long i=0;i<Largl;i++) 
					{
						mem[i+j*Largl]= pImageWnd->Image[i+j*Largl];	
					}
					Dlg.StepIt();
				}
				delete pImageWnd;
				Dlg.DestroyWindow();
			}
	   }
	   else 
	   {
		   MessageBox(NULL,"Ce fichier n'existe pas","Message",MB_OK|MB_ICONEXCLAMATION);
		   return FALSE;
	   }
  }
  else  return FALSE;
  delete FileData;
  return TRUE;
};

BOOL CParam::RehaussementFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{	
	long val[2],p,q,r,temp,moy,max,min;
	int l,c,lc,cc;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Rehaussement fix: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,Hautl-val[0]+2);
	Dlg.SetStep(val[0]);

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	

	for (l=0;l<=Hautl-val[0];l+=val[0])
	{
		p=l*Largl;
		for(c=0;c<=Largl-val[1];c+=val[1])
		{
			q=0;	
			moy=0;
			max=0;
			min=255;
			for(lc=0;lc<val[0];lc++)
			{
				r=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					moy+=temp=imi[r++];
					if(temp>max) max=temp;
					if(temp<min) min=temp;
					
				}//for cc
			}//for lc
			 
			
			if((moy/(val[0]*val[1]))>((max+min)/2)) temp=max;
			else temp=min;
			for(lc=0;lc<val[0];lc++)
			{
				r=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
						 mem[r++]=(byte) temp;
				}//for cc
			}//for lc
		p+=val[1];
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	return TRUE;
}//rehaussementfix

BOOL CParam::RehaussementGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{	
	long val[2],p,q,r,max,min,temp;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Rehaussement glis: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	

	for(l=0;l<Hautl;l++)
	{	
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}//for l for c
	}
	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			q=0;
			max=0;
			min=255;
			for(lc=debutlc;lc<=finlc;lc++)
			{
				r=p+lc*Largl+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					temp=imi[r++];
					if(temp>max) max=temp;
					if(temp<min) min=temp;
				}//for cc
			}//for lc
			
			if(imi[p]>(max+min)/2)
				mem[p++]=(byte) max;
			else mem[p++]=(byte) min;
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	return TRUE;
}


BOOL CParam::Ouverture(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2],p,r;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	unsigned char *tampon,temp,max,min;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Ouverture: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,2*Hautl);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	tampon=new unsigned char[Largl*Hautl];
	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			tampon[p]=imi[p++];
		}//for l for c
	}

	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			min=255;
			for(lc=debutlc;lc<=finlc;lc++)
			{
				r=p+debutcc+Largl*lc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					temp=imi[r++];
					if(temp<=min) min=temp;
				}//for cc
			}//for lc
			tampon[p++]=min;
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l

	for (l=debutl;l<finl;l++)
		{
			p=Largl*l+debutc;
			for(c=debutc;c<finc;c++)
			{
				max=0;
				for(lc=debutlc;lc<=finlc;lc++)
				{
					r=p+debutcc+Largl*lc;
					for(cc=debutcc;cc<=fincc;cc++)
					{
						temp=tampon[r++];
						if(temp>=max) max=temp;
					}//for cc
				}//for lc
				mem[p++]=max;
			}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
		}//for l

	Dlg.DestroyWindow();
	delete tampon;
	return TRUE;
}//Ouverture

BOOL CParam::OuvertureRapide(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2],l,stop;
	unsigned char *tampon;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Ouverture rapide: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,100);
	Dlg.SetStep(50);

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	tampon=new unsigned char[Largl*Hautl];

	MinimumGlissant(mem,imi,Largl,Hautl,val[1],val[0]);
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	MaximumGlissant(tampon,mem,Largl,Hautl,val[1],val[0]);
	stop=2*(val[1]/2)+Largl*2*(val[0]/2);
	l=stop;
	Dlg.StepIt();
	int i;
	for( i=0;i<Largl*Hautl-stop;i++) mem[i]=tampon[l++];
	for(i=Largl*Hautl-stop;i<Largl*Hautl;i++) mem[i]=0;
	Dlg.DestroyWindow();
	delete tampon;
	return TRUE;
}//Ouverture rapid

BOOL CParam::FermetureRapide(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2],l,stop;
	unsigned char *tampon;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Fermeture rapide: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,100);
	Dlg.SetStep(50);

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	tampon=new unsigned char[Largl*Hautl];

	MaximumGlissant(mem,imi,Largl,Hautl,val[1],val[0]);
	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	MinimumGlissant(tampon,mem,Largl,Hautl,val[1],val[0]);
	stop=2*(val[1]/2)+Largl*2*(val[0]/2);
	l=stop;
	Dlg.StepIt();
	int i;
	for( i=0;i<Largl*Hautl-stop;i++) mem[i]=tampon[l++];
	for(i=Largl*Hautl-stop;i<Largl*Hautl;i++) mem[i]=0;
	Dlg.DestroyWindow();
	delete tampon;
	return TRUE;
}//Ouverture

void CParam::MaximumGlissant(unsigned char *maximumc, unsigned char *imi,int Largl,int Hautl,int LF,int HF)
{
	int q,l;
	unsigned char classe[256];
	unsigned char Max,Maxc,courant,ancien;

	unsigned char *maximuml;

	
	maximuml=new unsigned char[Largl*Hautl];

	for(int i=0;i<64;i++) classe[i]=classe[255-i]=classe[127-i]=classe[128+i]=0;
	//traitement ligne
	q=l=0;
	int i,j;
	//////////////////////////////lignes
	for( j=0;j<Hautl;j++)
	{
		Max=0;
		for( i=0;i<LF;i++)
		{
			courant=imi[q];
			classe[courant]++;
			Max=(((courant) > (Max)) ? (courant) : (Max));
			//Max=__max(courant,Max);
			maximuml[q++]=Max;
			
		}//for i 0->11
		
		for(i=LF;i<Largl;i++)
		{
			courant=imi[q];
			ancien=imi[l++];
			classe[courant]++;
			Max=(((courant) > (Max)) ? (courant) : (Max));
			classe[ancien]--;
			if(ancien==Max) 
			{
				while(!classe[Max])
				{
					Max--;
				}
			}
			maximuml[q++]=Max;
		}//for i LF->Largl
	for(i=0;i<LF;i++) classe[imi[l++]]--;
	}///for j+

	///////////////////////////////colonnes
	for(i=0;i<Largl;i++)
	{
		q=l=i;Maxc=0;
		for(j=0;j<HF;j++)
		{
			courant=maximuml[q];
			classe[courant]++;
			Maxc=(((courant) > (Maxc)) ? (courant) : (Maxc));
			maximumc[q]=Maxc;
			q+=Largl;
		}

		for(j=HF;j<Hautl;j++)
		{
			courant=maximuml[q];
			ancien=maximuml[l];
			classe[courant]++;
			Maxc=(((courant) > (Maxc)) ? (courant) : (Maxc));
			classe[ancien]--;
			if(ancien==Maxc) 
			{
				while(!classe[Maxc])
				{
					Maxc--;
				}
			}
			maximumc[q]=Maxc;
			q+=Largl;
			l+=Largl;
		}
		for(j=0;j<HF;j++) {classe[maximuml[l]]--;l+=Largl;}
	}//for i
	delete maximuml;
}

void CParam::MinimumGlissant(unsigned char *Minimumc,unsigned char *imi,int Largl,int Hautl,int LF,int HF)
{
	SetCursor(LoadCursor(0,IDC_WAIT));
	int q,l;
	unsigned char classe[256];
	unsigned char Min,Minc,courant,ancien;

	unsigned char *Minimuml;

	
	Minimuml=new unsigned char[Largl*Hautl];

	for(int i=0;i<64;i++) classe[i]=classe[255-i]=classe[127-i]=classe[128+i]=0;
	//traitement ligne
	q=l=0;
	int i,j;
	//////////////////////////////lignes
	for( j=0;j<Hautl;j++)
	{
		Min=255;
		for( i=0;i<LF;i++)
		{
			courant=imi[q];
			classe[courant]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			//Min=__Min(courant,Min);
			Minimuml[q++]=Min;
			
		}//for i 0->11
		
		for(i=LF;i<Largl;i++)
		{
			courant=imi[q];
			ancien=imi[l++];
			classe[courant]++;
			Min=(((courant) < (Min)) ? (courant) : (Min));
			classe[ancien]--;
			if(ancien==Min) 
			{
				while(!classe[Min])
				{
					Min++;
				}
			}
			Minimuml[q++]=Min;
		}//for i LF->Largl
	for(i=0;i<LF;i++) classe[imi[l++]]--;
	}///for j+

	///////////////////////////////colonnes
	for(i=0;i<Largl;i++)
	{
		q=l=i;Minc=255;
		for( j=0;j<HF;j++)
		{
			courant=Minimuml[q];
			classe[courant]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			Minimumc[q]=Minc;
			q+=Largl;
		}

		for(j=HF;j<Hautl;j++)
		{
			courant=Minimuml[q];
			ancien=Minimuml[l];
			classe[courant]++;
			Minc=(((courant) < (Minc)) ? (courant) : (Minc));
			classe[ancien]--;
			if(ancien==Minc) 
			{
				while(!classe[Minc])
				{
					Minc++;
				}
			}
			Minimumc[q]=Minc;
			q+=Largl;
			l+=Largl;
		}
		for(j=0;j<HF;j++) {classe[Minimuml[l]]--;l+=Largl;}
	}//for i
	delete Minimuml;
	SetCursor(LoadCursor(0,IDC_ARROW));
}

BOOL CParam::Fermeture(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2],p,q;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	unsigned char *tampon,temp,max,min;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Fermeture: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,2*Hautl);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	tampon=new unsigned char[Largl*Hautl];
	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			tampon[p]=imi[p++];
		}//for l for c
	}

	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			max=0;
			for(lc=debutlc;lc<=finlc;lc++)
			{
				q=p+lc*Largl+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					temp=imi[q++];
					if(temp>=max) max=temp;
				}//for cc
			}//for lc
			
			tampon[p++]=max;
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l

	for (l=debutl;l<finl;l++)
		{
			p=Largl*l+debutc;
			for(c=debutc;c<finc;c++)
			{
				min=255;
				for(lc=debutlc;lc<=finlc;lc++)
				{
					q=p+lc*Largl+debutcc;
					for(cc=debutcc;cc<=fincc;cc++)
					{
						temp=tampon[q++];
						if(temp<=min) min=temp;
					}//for cc
				}//for lc
				
				mem[p++]=min;
			}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
		}//for l

	Dlg.DestroyWindow();
	delete tampon;
	return TRUE;
}//Fermeture

BOOL CParam::EcartmedFen(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[3],p,s,qot,q,r;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	unsigned char *place;
	unsigned char temp;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Ecart/Median fen glis: %d,%d,%d",val[0],val[1],val[2]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	qot=val[0]*val[1];

	if((val[0]<1)||(val[1]<1)||(val[2]<0)||(val[2]>255)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	place=new unsigned char[val[0]*val[1]];

	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}//for l for c
	}
	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			s=0;
			r=0;
			for(lc=debutlc;lc<=finlc;lc++)
			{
				q=p+Largl*lc+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					place[r++]=imi[q++];
				}//for cc
			}//for lc
			qsort(place,(size_t) (val[1]*val[0]),sizeof(unsigned char),comparaison);
			temp=place[(val[0]*val[1]+1)/2];
			for(lc=debutlc;lc<=finlc;lc++)
			{
				q=p+Largl*lc+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					if(labs(imi[q++]-temp)<=val[2]) 
					{
						s++;						
					}
				
				}//for cc
			}//for lc
			
			mem[p++]=(byte) ((255*s)/(qot));
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	delete place;
	return TRUE;
}//EcartmoyFen

BOOL CParam::EcartmedFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[3],p,s,q,qot,t;
	byte r;
	int l,c,lc,cc;
	unsigned char *place;
	unsigned char temp;
	Decode(chaine,val);
	char ch[50];
	sprintf(ch,"Ecart/Median fen fix: %d,%d,%d",val[0],val[1],val[2]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,Hautl-val[0]+2);
	Dlg.SetStep(val[0]);
	
	qot=val[0]*val[1];

	if((val[0]<1)||(val[1]<1)||(val[2]<0)||(val[2]>255)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	place=new unsigned char[val[0]*val[1]];

	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}//for l for c
	}
	for (l=0;l<=Hautl-val[0];l+=val[0])
	{
		p=Largl*l;
		for(c=0;c<=Largl-val[1];c+=val[1])
		{
			s=0;
			t=0;
			for(lc=0;lc<val[0];lc++)
			{
				q=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					place[t++]=imi[q++];
				
				}//for cc
			}//for lc
			qsort(place,(size_t) (val[1]*val[0]),sizeof(unsigned char),comparaison);
			temp=place[(val[0]*val[1]+1)/2];

			for(lc=0;lc<val[0];lc++)
			{
				q=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					if(labs(imi[q++]-temp)<=val[2]) 
					{
						s++;						
					}
				
				}//for cc
			}//for lc
			
			r=(byte) ((255*s)/(qot));
			for(lc=0;lc<val[0];lc++)
			{
				q=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					mem[q++]=r;
				
				}//for cc
			}//for lc
			p+=val[1];
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	delete place;
	return TRUE;
}//EcartmedFix


BOOL CParam::MoyAdap(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[3],p,r,s,q;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Moy. adap.& seuil: %d,%d,%d",val[0],val[1],val[2]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	if((val[0]<1)||(val[1]<1)||(val[2]<0)||(val[2]>255)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	
	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}//for l for c
	}
	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			r=0;
			s=0;

			for(lc=debutlc;lc<=finlc;lc++)
			{
				q=p+Largl*lc+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					if(labs(imi[q]-imi[p])<=val[2]) 
					{
						r+=imi[q];
						s++;
					}
					q++;
				}//for cc
			}//for lc
			
			mem[p++]=(byte) (r/s);
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	return TRUE;
}//MoyAdap

BOOL CParam::MoyAdapAuto(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2],p,r,s,q,t,temp;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	unsigned char *place;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Moy. adap.auto.: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	place=new unsigned char[val[0]*val[1]];

	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}//for l for c
	}

	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			r=0;
			s=0;
			q=0;

			for(lc=debutlc;lc<=finlc;lc++)
			{
				t=p+Largl*lc+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					place[q++]=imi[t++];
				
				}//for cc
			}//for lc
			qsort(place,(size_t) (val[1]*val[0]),sizeof(unsigned char),comparaison);

			temp=place[(val[1]*val[0]+1)/2];
			for(int i=0;i<val[0]*val[1];i++)place[i]= abs(place[i]-temp);//place[i]-place[(val[1]*val[0]+1)/2];
			qsort(place,(size_t) (val[1]*val[0]),sizeof(unsigned char),comparaison);
			
			temp=place[(val[1]*val[0]+1)/2];
			for(lc=debutlc;lc<=finlc;lc++)
			{
				t=p+Largl*lc+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					if(labs(imi[t]-imi[p])<=temp) 
					{
						r+=imi[t];
						s++;
					}
					t++;

				}//for cc
			}//for lc
			mem[p++]=(byte) (r/s);
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	delete place;
	return TRUE;
}//MoyAdapauto

BOOL CParam::MoyTronqueglis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{	
	long val[3],p,q,m;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	long t,r;
	unsigned char *place;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Moy. tronquée glis: %d,%d,%d",val[0],val[1],val[2]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	if((val[0]<1)||(val[1]<1)||(val[2]<0)||(val[2]>(val[0]*val[1]))) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	place=new unsigned char[val[0]*val[1]];

	r=((val[0]*val[1]+1)/2)-val[2];
	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}//for l for c
	}

	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			q=0;
			t=0;
			for(lc=debutlc;lc<=finlc;lc++)
			{
				m=p+Largl*lc+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					place[q++]=imi[m++];
				}//for cc
			}//for lc
			qsort(place,(size_t) (val[1]*val[0]),sizeof(unsigned char),comparaison);

			for(q=0;q<=2*val[2];q++)	t+=place[r+q];

			mem[p++]=(byte)(t/(2*val[2]+1));;
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	delete place;
	return TRUE;
}//Moytronquéeglis

BOOL CParam::MoyTronqueFix(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{	
	long val[3],p,q,m;
	int l,c,lc,cc;
	long t,r;
	unsigned char *place;
	Decode(chaine,val);
	char ch[50];
	sprintf(ch,"Moy. tronquée fix: %d,%d,%d",val[0],val[1],val[2]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,Hautl-val[0]+2);
	Dlg.SetStep(val[0]);
	if((val[0]<1)||(val[1]<1)||(val[2]<0)||(val[2]>(val[0]*val[1]))) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	
	place=new unsigned char[val[0]*val[1]];

	r=((val[0]*val[1])/2)-val[2];
	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}//for l for c
	}
	for (l=0;l<=Hautl-val[0];l+=val[0])
	{
		p=Largl*l;
		for(c=0;c<=Largl-val[1];c+=val[1])
		{
			q=0;
			t=0;
			for(lc=0;lc<val[0];lc++)
			{
				m=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					place[q++]=imi[m++];
				}//for cc
			}//for lc

			qsort(place,(size_t) (val[1]*val[0]),sizeof(unsigned char),comparaison);

			for(q=0;q<=2*val[2];q++)	t+=place[r+q];

			for(lc=0;lc<val[0];lc++)
			{
				m=p+Largl*lc;
				for(cc=0;cc<val[1];cc++)
				{
					mem[m++]=(byte)(t/(2*val[2]+1));
				}//for cc
			}//for lc
			p+=val[1];
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	delete place;
	return TRUE;
}//MoytronquéeFix

BOOL CParam::Maxrepglis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{	
	long val[2],p,q;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	long max,t;
	unsigned char histo[256];
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Max histo glis: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
		
	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}//for l for c
	}

	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			max=0;
			for(int i=0;i<256;i++) histo[i]=0;
			for(lc=debutlc;lc<=finlc;lc++)
			{
				q=p+Largl*lc+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					histo[imi[q]]++;
					t=histo[imi[q]];
					if(t>max) {max=t;mem[p]=imi[q];}
					q++;
				}//for cc
			}//for lc
			p++;
		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	return TRUE;
}//maxrepglis

BOOL CParam::NbMaxrepglis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{	
	long val[2],p,q,qot;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	long max,t;
	unsigned char histo[256];
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"NbMax histo glis: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	qot=val[0]*val[1];
	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
		
	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}//for l for c
	}

	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			max=0;
			for(int i=0;i<256;i++) histo[i]=0;
			for(lc=debutlc;lc<=finlc;lc++)
			{
				q=p+Largl*lc+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					histo[imi[q]]++;
					t=histo[imi[q]];
					if(t>max) {max=t;}
					q++;
				}//for cc
			}//for lc
			mem[p++]=(byte) ((255*max)/(qot));

		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	return TRUE;
}//Nbmaxrepglis

BOOL CParam::EcartTypeFenGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2],p,q,qot;
	int l,c,lc,cc,debutl,finl,debutc,finc,debutlc,finlc,debutcc,fincc;
	long carre,moy;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Ecart typ fen glis: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange((int)(val[0]-1)/2,Hautl-(int)val[0]/2+1);
	Dlg.SetStep(1);

	debutl=(int)(val[0]-1)/2;
	finl=Hautl-(int)val[0]/2;
	debutc=(int)(val[1]-1)/2;
	finc=Largl-(int)val[1]/2;	
	debutlc=-debutl;
	finlc=(int)val[0]/2;
	debutcc=-debutc;
	fincc=(int)val[1]/2;

	qot=val[0]*val[1];
	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
		
	for(l=0;l<Hautl;l++)
	{
		p=Largl*l;
		for(c=0;c<Largl;c++)
		{
			mem[p]=imi[p++];
		}//for l for c
	}

	for (l=debutl;l<finl;l++)
	{
		p=Largl*l+debutc;
		for(c=debutc;c<finc;c++)
		{
			carre=0;
			moy=0;
			for(lc=debutlc;lc<=finlc;lc++)
			{
				q=p+Largl*lc+debutcc;
				for(cc=debutcc;cc<=fincc;cc++)
				{
					carre+=imi[q]*imi[q];
					moy+=imi[q];
					q++;
				}//for cc
			}//for lc
			mem[p++]=(byte)sqrt((double)(carre/qot-moy*moy/(qot*qot)));

		}//for c
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}//for l
	
	Dlg.DestroyWindow();
	return TRUE;
}//ecarttypfenglis


BOOL CParam::MaxFenGlis(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2],l,stop;
	Decode(chaine,val);
	char ch[40];
	sprintf(ch,"Max glis: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,100);
	Dlg.SetStep(50);

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	

	if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
	Dlg.StepIt();
	MaximumGlissant(mem,imi,Largl,Hautl,val[1],val[0]);
	stop=int(val[1]/2)+Largl*int(val[0]/2);
	l=stop;
	Dlg.StepIt();
	int i;
	for( i=0;i<Largl*Hautl-stop;i++) mem[i]=mem[l++];
	for(i=Largl*Hautl-stop;i<Largl*Hautl;i++) mem[i]=0;
	Dlg.DestroyWindow();
	return TRUE;
}
BOOL CParam::DensMax33(CString chaine, byte   *mem, byte   *imi, long Largl, long Hautl)
{
	long val[2];
	Decode(chaine,val);
	int c,l,compt;
	int co=val[1]/2;
	int lo=val[0]/2;
	int nbpix=(2*lo+1)*(2*co+1);
	unsigned char cour;
	BOOL cod ,*tampon=new BOOL[Largl*Hautl];
	char ch[40];
	sprintf(ch,"Densite Max loc3,3: %d,%d",val[0],val[1]);
	CProgressDlg Dlg;
	Dlg.Create();
	Dlg.SetWindowText(ch);
	Dlg.SetRange(0,Hautl);
	Dlg.SetStep(1);

	if((val[0]<1)||(val[1]<1)) {MessageBox(NULL,"erreur de parametre","Message",MB_OK|MB_ICONSTOP);return FALSE;}
	long i;
	for( i=0;i<Largl*Hautl;i++) {tampon[i]=0;mem[i]=0;}

	for( i=Largl;i<Largl*Hautl-Largl;i++)
	{ 
		cod=1;
		cour=imi[i];
		for(c=-1;c<2 && cod;c++)
			for(l=-1;l<2;l++)
			{
				if(imi[i+c+l*Largl]>cour) {cod=0;break;}
			}
		tampon[i]=cod;
	}

	for(long j=lo;j<Hautl-lo;j++)
	{
		for(i=co;i<Largl-co;i++)
		{ 
			compt=0;
			for(c=-co;c<=co;c++)
				for(l=-lo;l<=lo;l++)
				{
					if(tampon[i+c+(l+j)*Largl]) compt++;
				}
			mem[i+j*Largl]=(unsigned char) (255*((float) compt/nbpix));
		}
		if(Dlg.CheckCancelButton()) {MessageBox(NULL,"Les paramètres ne seront pas calculés","Message",MB_OK|MB_ICONSTOP);mem=imi;return FALSE;}
		Dlg.StepIt();
	}
	Dlg.DestroyWindow();
	delete tampon;
	return TRUE;
}
//////////////////////fonction de décodage des charactères et comparaison///////////////////////

void CParam::Decode(CString chaine,long * val)
{
	int Pos,i=0;
	CString Cval;

	Pos=chaine.Find(":");;
	if(Pos==-1) { MessageBox(NULL,"erreur dans les parametres","Message",MB_OK|MB_ICONSTOP);return;}
	Cval=chaine.Right(chaine.GetLength()-Pos-1);
	do
	{
		Pos=Cval.Find(",");
		if(Pos==-1) {val[i]=atol(Cval);}
		else 
		{
			val[i]=atol(Cval.Left(Pos));
			Cval=Cval.Right(Cval.GetLength()-Pos-1);
		}

		i++;

	} while (Pos!=-1);
}

int  comparaison(const void* p,const void* q)
{
  return (int) (*(unsigned char *)p-*(unsigned char *)q);
}


	///////////////////////////FIN//////////////////////////////////




