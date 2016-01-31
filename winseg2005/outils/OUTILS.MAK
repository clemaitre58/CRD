# Microsoft Developer Studio Generated NMAKE File, Format Version 4.20
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

!IF "$(CFG)" == ""
CFG=OUTILS - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to OUTILS - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "OUTILS - Win32 Release" && "$(CFG)" != "OUTILS - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "OUTILS.MAK" CFG="OUTILS - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OUTILS - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "OUTILS - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "OUTILS - Win32 Debug"
CPP=cl.exe

!IF  "$(CFG)" == "OUTILS - Win32 Release"

# PROP BASE Use_MFC 1
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 1
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\Release
INTDIR=.\Release

ALL : "$(OUTDIR)\OUTILS.lib"

CLEAN : 
	-@erase "$(INTDIR)\COLOR.OBJ"
	-@erase "$(INTDIR)\ECHANT.OBJ"
	-@erase "$(INTDIR)\ETIQUETT.OBJ"
	-@erase "$(INTDIR)\FOURIER.OBJ"
	-@erase "$(INTDIR)\INVERSE.OBJ"
	-@erase "$(INTDIR)\MIROIR.OBJ"
	-@erase "$(INTDIR)\ROTATION.OBJ"
	-@erase "$(INTDIR)\SEGNOP.OBJ"
	-@erase "$(INTDIR)\SEUIL.OBJ"
	-@erase "$(INTDIR)\VECTEUR.OBJ"
	-@erase "$(INTDIR)\WAVELET.OBJ"
	-@erase "$(OUTDIR)\OUTILS.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /MT /W3 /GX /Ox /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /FR /YX /c
# ADD CPP /nologo /MT /W3 /GX /Ox /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /YX /c
# SUBTRACT CPP /Fr
CPP_PROJ=/nologo /MT /W3 /GX /Ox /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D\
 "_MBCS" /Fp"$(INTDIR)/OUTILS.pch" /YX /Fo"$(INTDIR)/" /c 
CPP_OBJS=.\Release/
CPP_SBRS=.\.
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/OUTILS.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/OUTILS.lib" 
LIB32_OBJS= \
	"$(INTDIR)\COLOR.OBJ" \
	"$(INTDIR)\ECHANT.OBJ" \
	"$(INTDIR)\ETIQUETT.OBJ" \
	"$(INTDIR)\FOURIER.OBJ" \
	"$(INTDIR)\INVERSE.OBJ" \
	"$(INTDIR)\MIROIR.OBJ" \
	"$(INTDIR)\ROTATION.OBJ" \
	"$(INTDIR)\SEGNOP.OBJ" \
	"$(INTDIR)\SEUIL.OBJ" \
	"$(INTDIR)\VECTEUR.OBJ" \
	"$(INTDIR)\WAVELET.OBJ"

"$(OUTDIR)\OUTILS.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "OUTILS - Win32 Debug"

# PROP BASE Use_MFC 1
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 1
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "$(OUTDIR)\OUTILS.lib"

CLEAN : 
	-@erase "$(INTDIR)\COLOR.OBJ"
	-@erase "$(INTDIR)\ECHANT.OBJ"
	-@erase "$(INTDIR)\ETIQUETT.OBJ"
	-@erase "$(INTDIR)\FOURIER.OBJ"
	-@erase "$(INTDIR)\INVERSE.OBJ"
	-@erase "$(INTDIR)\MIROIR.OBJ"
	-@erase "$(INTDIR)\ROTATION.OBJ"
	-@erase "$(INTDIR)\SEGNOP.OBJ"
	-@erase "$(INTDIR)\SEUIL.OBJ"
	-@erase "$(INTDIR)\VECTEUR.OBJ"
	-@erase "$(INTDIR)\WAVELET.OBJ"
	-@erase "$(OUTDIR)\OUTILS.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /G3 /MTd /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /FR /YX /c
# ADD CPP /nologo /G3 /MTd /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /YX /c
# SUBTRACT CPP /Fr
CPP_PROJ=/nologo /G3 /MTd /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS"\
 /D "_MBCS" /Fp"$(INTDIR)/OUTILS.pch" /YX /Fo"$(INTDIR)/" /c 
CPP_OBJS=.\Debug/
CPP_SBRS=.\.
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/OUTILS.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/OUTILS.lib" 
LIB32_OBJS= \
	"$(INTDIR)\COLOR.OBJ" \
	"$(INTDIR)\ECHANT.OBJ" \
	"$(INTDIR)\ETIQUETT.OBJ" \
	"$(INTDIR)\FOURIER.OBJ" \
	"$(INTDIR)\INVERSE.OBJ" \
	"$(INTDIR)\MIROIR.OBJ" \
	"$(INTDIR)\ROTATION.OBJ" \
	"$(INTDIR)\SEGNOP.OBJ" \
	"$(INTDIR)\SEUIL.OBJ" \
	"$(INTDIR)\VECTEUR.OBJ" \
	"$(INTDIR)\WAVELET.OBJ"

"$(OUTDIR)\OUTILS.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

################################################################################
# Begin Target

# Name "OUTILS - Win32 Release"
# Name "OUTILS - Win32 Debug"

!IF  "$(CFG)" == "OUTILS - Win32 Release"

!ELSEIF  "$(CFG)" == "OUTILS - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\FOURIER.CPP
DEP_CPP_FOURI=\
	".\fourier.h"\
	".\struct.h"\
	

"$(INTDIR)\FOURIER.OBJ" : $(SOURCE) $(DEP_CPP_FOURI) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\INVERSE.CPP
DEP_CPP_INVER=\
	".\inverse.h"\
	".\struct.h"\
	

"$(INTDIR)\INVERSE.OBJ" : $(SOURCE) $(DEP_CPP_INVER) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\MIROIR.CPP
DEP_CPP_MIROI=\
	".\miroir.h"\
	".\struct.h"\
	

"$(INTDIR)\MIROIR.OBJ" : $(SOURCE) $(DEP_CPP_MIROI) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\SEGNOP.CPP
DEP_CPP_SEGNO=\
	".\segnop.h"\
	".\struct.h"\
	

"$(INTDIR)\SEGNOP.OBJ" : $(SOURCE) $(DEP_CPP_SEGNO) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\SEUIL.CPP
DEP_CPP_SEUIL=\
	".\seuil.h"\
	".\struct.h"\
	

"$(INTDIR)\SEUIL.OBJ" : $(SOURCE) $(DEP_CPP_SEUIL) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\VECTEUR.CPP
DEP_CPP_VECTE=\
	".\vecteur.h"\
	

"$(INTDIR)\VECTEUR.OBJ" : $(SOURCE) $(DEP_CPP_VECTE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ETIQUETT.CPP

!IF  "$(CFG)" == "OUTILS - Win32 Release"

DEP_CPP_ETIQU=\
	".\etiquett.h"\
	".\seuil.h"\
	".\struct.h"\
	

"$(INTDIR)\ETIQUETT.OBJ" : $(SOURCE) $(DEP_CPP_ETIQU) "$(INTDIR)"


!ELSEIF  "$(CFG)" == "OUTILS - Win32 Debug"

DEP_CPP_ETIQU=\
	".\etiquett.h"\
	".\seuil.h"\
	".\struct.h"\
	

"$(INTDIR)\ETIQUETT.OBJ" : $(SOURCE) $(DEP_CPP_ETIQU) "$(INTDIR)"


!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\ROTATION.CPP
DEP_CPP_ROTAT=\
	".\rotation.h"\
	".\struct.h"\
	

"$(INTDIR)\ROTATION.OBJ" : $(SOURCE) $(DEP_CPP_ROTAT) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\WAVELET.CPP
DEP_CPP_WAVEL=\
	".\struct.h"\
	".\wavelet.h"\
	

"$(INTDIR)\WAVELET.OBJ" : $(SOURCE) $(DEP_CPP_WAVEL) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ECHANT.CPP
DEP_CPP_ECHAN=\
	".\echant.h"\
	".\struct.h"\
	

"$(INTDIR)\ECHANT.OBJ" : $(SOURCE) $(DEP_CPP_ECHAN) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\COLOR.CPP
DEP_CPP_COLOR=\
	".\color.h"\
	".\struct.h"\
	

"$(INTDIR)\COLOR.OBJ" : $(SOURCE) $(DEP_CPP_COLOR) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
