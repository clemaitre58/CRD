#pragma once


// Boîte de dialogue CDialogAffine

class CDialogAffine : public CDialog
{
	DECLARE_DYNAMIC(CDialogAffine)

public:
	CDialogAffine(CWnd* pParent = NULL);   // constructeur standard
	virtual ~CDialogAffine();

// Données de boîte de dialogue
	enum { IDD = IDD_DIALOG_PARAMAFFINE };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // Prise en charge de DDX/DDV

	DECLARE_MESSAGE_MAP()
public:
	double m_M1;
	double m_M2;
	double m_M3;
	double m_M4;
	double m_M5;
	double m_M6;
	double m_M7;
};
