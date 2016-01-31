#ifndef __STRUCT_H__     
#define __STRUCT_H__  

struct rgb		{unsigned char b;unsigned char g;unsigned char r;};    

struct SMasque
{  
	float *Filtre;
	int Larg;
	int Haut;
	int Centre_col;
	int Centre_lig;
};

struct SFiltre
{
 float *Filtre;
 int Centre;
 int Fin;
}; 

struct SImage
{ 
 	long Larg;
 	long Haut;                   
 	long Lmin,Lmax,Hmin,Hmax;
    unsigned char  * Image; 
}; 

struct SImageRGB
{ 
 	long Larg;
 	long Haut;                   
 	long Lmin,Lmax,Hmin,Hmax;
    rgb  * Image; 
}; 
    
    
struct SImageFloat
{ 
 	long Larg;
 	long Haut;                   
 	long Lmin,Lmax,Hmin,Hmax;
    float  * Image; 
}; 
    
#endif // __STRUCT_H__  