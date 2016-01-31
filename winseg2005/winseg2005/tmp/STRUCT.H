#ifndef __STRUCT_H__     
#define __STRUCT_H__     

enum NBCOLOR {RIEN=0,NB,COLOR};
enum TYPEIMAGE {BMP, MEM, FLOAT32};

enum LETYPE {IMAGE_FIXE,IMAGE_CONTINUE};
       

struct rgb
{
	unsigned char b;
	unsigned char g;
	unsigned char r;
};

struct drgb	
{
	double b;
	double g;
	double r;
};    

struct SMasque
{  
	float *Filtre;
	int Larg;
	int Haut;
	int Centre_col;
	int Centre_lig;
};

struct SImage
{ 
 	long Larg;
 	long Haut;                   
 	long Lmin,Lmax,Hmin,Hmax;
  	unsigned char  * Image; 
};

struct SFiltre
{
 float *Filtre;
 int Centre;
 int Fin;
}; 

struct SImageFloat
{ 
 	long Larg;
 	long Haut;                   
 	long Lmin,Lmax,Hmin,Hmax;
  	float  * Image; 
};

struct SImageRGB
{ 
 	long Larg;
 	long Haut;                   
 	long Lmin,Lmax,Hmin,Hmax;
  	rgb  * Image; 
};

struct SComplexe
{
	float Re;
	float Im;
};

struct ConnexityDegreeAFreemanCode
{
	unsigned char ConnexityDegree;
	BOOL FreemanCode[8];
};

struct ParcoursFreemanHoraire
{
	char Code[8];
};

struct ParcoursFreemanAntiHoraire
{
	char Code[8];
};

struct DirectPrecedent
{
	char Code[8];
};

struct Sfpoint
{
	double x;
	double y;
};

struct STabIndice
{
	double r;
	int i;
};


struct PosDetectCurvi
{
	CPoint PositionMaxi;
	CPoint PosInit;
	CPoint PosFinal;
	double Largeur;
	double ValDist;
	Sfpoint PosAxe;
};

struct SBoundingBox
{
	Sfpoint AxeMin;
	Sfpoint AxeMax;
};
struct RespCurvi
{
	int IptInit;
	int IptFinal;
	double ValDist;
	double Largeur;
	int Centre;
};

typedef float Pixelp;

struct ImagePetrou 
{
  Pixelp **ima;
  int dim_x;
  int dim_y;
} ;

#endif // __STRUCT_H__