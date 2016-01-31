#pragma once
#include "afxcmn.h"


// Boîte de dialogue BoiteParamCurviDetect

class BoiteParamCurviDetect : public CDialog
{
	DECLARE_DYNAMIC(BoiteParamCurviDetect)

public:
	BoiteParamCurviDetect(CWnd* pParent = NULL);   // constructeur standard
	virtual ~BoiteParamCurviDetect();

// Données de boîte de dialogue
	enum { IDD = IDD_BoiteParamCurviDetect };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // Prise en charge de DDX/DDV

	DECLARE_MESSAGE_MAP()
public:
	CSliderCtrl SliderLargeur;
	int ValLargeur;
	CSliderCtrl SliderLargeurFourier;
	int ValLargeurFourier;
	int NbListeMini;
	CSliderCtrl SliderLongueurForme;
	double SeuilReponse;
	CSliderCtrl SliderSeuilReponse;
	afx_msg void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	virtual BOOL OnInitDialog();
protected:
	virtual void OnOK();
public:
	CSliderCtrl SliderLargeurMini;
	int ValLargeurMini;
	CSliderCtrl m_SliderVariationLocaleLargeur;
	double m_EditVariationLocaleLargeur;
	BOOL m_BDebugEdges;
	BOOL m_BDebugTracking;
	BOOL m_BDebugCurviReg;
};
