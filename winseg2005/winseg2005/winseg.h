#include "resource.h"

class CMdiApp : public CWinApp
{
	public:
	virtual BOOL PreTranslateMessage(MSG* pMsg);
  
    ~CMdiApp();
	virtual BOOL InitInstance();    

// Implementation

	//{{AFX_MSG(CMdiApp)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};                          

