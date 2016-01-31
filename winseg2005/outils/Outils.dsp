# Microsoft Developer Studio Project File - Name="OUTILS" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=OUTILS - Win32 Release
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Outils.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Outils.mak" CFG="OUTILS - Win32 Release"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OUTILS - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "OUTILS - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "OUTILS - Win32 Release"

# PROP BASE Use_MFC 1
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir ".\Release"
# PROP BASE Intermediate_Dir ".\Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 1
# PROP Use_Debug_Libraries 0
# PROP Output_Dir ".\Release"
# PROP Intermediate_Dir ".\Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MT /W3 /GX /Ox /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /FR /YX /c
# ADD CPP /nologo /MT /W3 /GX /Ox /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /FD /c
# SUBTRACT CPP /Fr
# ADD BASE RSC /l 0x40c
# ADD RSC /l 0x40c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "OUTILS - Win32 Debug"

# PROP BASE Use_MFC 1
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir ".\Debug"
# PROP BASE Intermediate_Dir ".\Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 1
# PROP Use_Debug_Libraries 1
# PROP Output_Dir ".\Debug"
# PROP Intermediate_Dir ".\Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /G3 /MTd /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /FR /YX /c
# ADD CPP /nologo /G3 /MTd /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /Fr /YX /FD /c
# ADD BASE RSC /l 0x40c
# ADD RSC /l 0x40c
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "OUTILS - Win32 Release"
# Name "OUTILS - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;hpj;bat;for;f90"
# Begin Source File

SOURCE=.\BSpline.cpp
# End Source File
# Begin Source File

SOURCE=.\COLOR.CPP
# End Source File
# Begin Source File

SOURCE=.\ECHANT.CPP
# End Source File
# Begin Source File

SOURCE=.\ETIQUETT.CPP
# End Source File
# Begin Source File

SOURCE=.\FOURIER.CPP
# End Source File
# Begin Source File

SOURCE=.\INVERSE.CPP
# End Source File
# Begin Source File

SOURCE=.\MIROIR.CPP
# End Source File
# Begin Source File

SOURCE=.\Morpho.cpp
# End Source File
# Begin Source File

SOURCE=.\ROTATION.CPP
# End Source File
# Begin Source File

SOURCE=.\SEGNOP.CPP
# End Source File
# Begin Source File

SOURCE=.\SEUIL.CPP
# End Source File
# Begin Source File

SOURCE=.\VECTEUR.CPP
# End Source File
# Begin Source File

SOURCE=.\WAVELET.CPP
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\BSpline.h
# End Source File
# Begin Source File

SOURCE=.\color.h
# End Source File
# Begin Source File

SOURCE=.\echant.h
# End Source File
# Begin Source File

SOURCE=.\etiquett.h
# End Source File
# Begin Source File

SOURCE=.\fourier.h
# End Source File
# Begin Source File

SOURCE=.\inverse.h
# End Source File
# Begin Source File

SOURCE=.\miroir.h
# End Source File
# Begin Source File

SOURCE=.\Morpho.h
# End Source File
# Begin Source File

SOURCE=.\rotation.h
# End Source File
# Begin Source File

SOURCE=.\segnop.h
# End Source File
# Begin Source File

SOURCE=.\seuil.h
# End Source File
# Begin Source File

SOURCE=.\struct.h
# End Source File
# Begin Source File

SOURCE=.\vecteur.h
# End Source File
# Begin Source File

SOURCE=.\wavelet.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;cnt;rtf;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
