// ChildView.h : interface de la classe CChildView
//


#pragma once


// fenêtre CChildView

class CChildView : public CWnd
{
// Construction
public:
	CChildView();

// Attributs
public:

// Opérations
public:

// Substitutions
	protected:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);

// Implémentation
public:
	virtual ~CChildView();

	// Fonctions générées de la table des messages
protected:
	afx_msg void OnPaint();
	DECLARE_MESSAGE_MAP()
};

