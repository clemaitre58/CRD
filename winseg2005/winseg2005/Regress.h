// Regress.h: interface for the CRegress class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_REGRESS_H__FCBE424F_3B50_47F0_9465_1555ACEC446D__INCLUDED_)
#define AFX_REGRESS_H__FCBE424F_3B50_47F0_9465_1555ACEC446D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CRegress  
{
public:
	double a,b,r;
	CRegress(POINT *Tab=NULL,int Taille=0);
	virtual ~CRegress();

};

#endif // !defined(AFX_REGRESS_H__FCBE424F_3B50_47F0_9465_1555ACEC446D__INCLUDED_)
