//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!			fichiers include personnels													!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include <windows.h>
#include "struct.h"

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!			D�claration des fonctions d'�tiquetage et Centres de Gravit�				!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 long  Etiquett(SImage imS,SImage imEtiqu,unsigned char seuil);
//     imS = structure SImage correspond a l image SOURCE
//     imEtiqu = structure SImage correspond a l image Resultat 
//	   pNbEtiqu = Nb de taches apres etiquettage


 void  CentresGravites(SImage imS,SImage imGravit, long NbEtiqu);
