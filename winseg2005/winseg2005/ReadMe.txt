================================================================================
    BIBLIOTHÈQUE MICROSOFT FOUNDATION CLASS : Vue d'ensemble du projet winseg2005
===============================================================================

L'Assistant Application a créé cette application winseg2005 pour  
vous.  Cette application décrit les principes de base de l'utilisation de  
Microsoft Foundation Classes et vous permet de créer votre application.

Ce fichier contient un résumé du contenu de chacun des fichiers qui 
constituent votre application winseg2005.

winseg2005.vcproj
    Il s'agit du fichier projet principal pour les projets VC++ générés à l'aide 
d'un Assistant Application. 
    Il contient les informations sur la version de Visual C++ qui a généré le 
fichier et  
    des informations sur les plates-formes, configurations et fonctionnalités du 
projet sélectionnées avec 
    l'Assistant Application.

winseg2005.h
    Il s'agit du fichier d'en-tête principal de l'application.  Il contient 
d'autres 
    en-têtes de projet spécifiques (y compris Resource.h) et déclare la 
    classe d'application Cwinseg2005App.

winseg2005.cpp
    Il s'agit du fichier source principal de l'application qui contient la 
    classe d'application Cwinseg2005App.

winseg2005.rc
    Il s'agit de la liste de toutes les ressources Microsoft Windows que le 
    programme utilise.  Elle comprend les icônes, les bitmaps et les curseurs 
qui sont stockés 
    dans le sous-répertoire RES.  Ce fichier peut être modifié directement dans 
Microsoft 
    Visual C++. Vos ressources de projet sont dans 1036.

res\winseg2005.ico
    Il s'agit d'un fichier icône, qui est utilisé comme icône de l'application.  
Cette 
    icône est incluse par le fichier de ressource principal winseg2005.rc.

res\winseg2005.rc2
    Ce fichier contient les ressources qui ne sont pas modifiées par Microsoft  
    Visual C++. Vous devez placer toutes les ressources non modifiables par 
    l'éditeur de ressources dans ce fichier.

/////////////////////////////////////////////////////////////////////////////

Pour la fenêtre frame principale :
    Le projet comprend une interface MFC standard.

MainFrm.h, MainFrm.cpp
    Ces fichiers contiennent la classe de frame CMainFrame 
dérivée de
    CMDIFrameWnd et qui contrôle toutes les fonctionnalités des frames MDI.

res\Toolbar.bmp
    Ce fichier bitmap sert à créer des images en mosaïque pour la barre 
d'outils.
    La barre d'outils et la barre d'état initiales sont construites dans la 
classe  
CMainFrame. Modifiez cette image de barre d'outils avec 
l'éditeur de ressources et 
    actualisez le tableau IDR_MAINFRAME TOOLBAR dans winseg2005.rc 
pour ajouter 
    les boutons de barre d'outils.
/////////////////////////////////////////////////////////////////////////////

Pour la fenêtre frame enfant :

ChildFrm.h, ChildFrm.cpp
    Ces fichiers définissent et implémentent la classe CChildFrame 
    qui prend en charge les fenêtres enfants d'une application MDI.

/////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////

Autres fonctionnalités :

Contrôles ActiveX
    L'application comprend la prise en charge des contrôles ActiveX.

/////////////////////////////////////////////////////////////////////////////

Autres fichiers standard :

StdAfx.h, StdAfx.cpp
    Ces fichiers sont utilisés pour générer un fichier d'en-tête précompilé 
(PCH) 
    nommé winseg2005.pch et un fichier de type précompilé nommé 
Stdafx.obj.

Resource.h
    Il s'agit du ficher d'en-tête standard, qui définit les nouveaux ID de 
ressources.
    Microsoft Visual C++ lit et met à jour ce fichier.

winseg2005.manifest
	Les fichiers manifestes d'application sont utilisés par Windows XP pour 
décrire les dépendances des applications  
	sur des versions spécifiques des assemblys côte à côte. Le chargeur 
utilise ces  
	informations pour charger l'assembly approprié à partir du cache de 
l'assembly ou directement  
	à partir de l'application. Le manifeste de l'application peut être inclus 
pour redistribution  
	comme fichier .manifest externe installé dans le même dossier que 
l'exécutable de l'application  
	ou être inclus dans l'exécutable sous la forme d'une ressource. 
/////////////////////////////////////////////////////////////////////////////

Autres remarques :

L'Assistant Application utilise "TODO:" pour indiquer les parties du code source 
où 
vous devrez ajouter ou modifier du code.

Si votre application utilise les MFC dans une DLL partagée vous devez redistribuer 
les DLL MFC. Si la langue de votre application n'est pas celle du système 
d'exploitation, vous devez également redistribuer le fichier des  ressources 
localisées MFC80XXX.DLL. Pour plus d'informations, consultez la section relative 
à la redistribution des applications Visual C++ dans la documentation MSDN.

/////////////////////////////////////////////////////////////////////////////
