// winseg2005.h : fichier d'en-tête principal pour l'application winseg2005
//
#pragma once

#ifndef __AFXWIN_H__
	#error "incluez 'stdafx.h' avant d'inclure ce fichier pour PCH"
#endif

#include "resource.h"       // symboles principaux


// Cwinseg2005App:
// Consultez winseg2005.cpp pour l'implémentation de cette classe
//

class Cwinseg2005App : public CWinApp
{
public:
	Cwinseg2005App();


// Substitutions
public:
	virtual BOOL InitInstance();
	virtual int ExitInstance();

// Implémentation
protected:
	HMENU m_hMDIMenu;
	HACCEL m_hMDIAccel;

public:
	afx_msg void OnAppAbout();
	afx_msg void OnFileNew();
	DECLARE_MESSAGE_MAP()
};

extern Cwinseg2005App theApp;