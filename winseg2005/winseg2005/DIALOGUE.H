#ifndef __DIALOGUE_H__     
#define __DIALOGUE_H__ 


class CDialZone : public CDialog
{ 	
	public:
	long *i1,*j1,*i2,*j2;  
	private:
	CWnd *LeParent;
	    
	protected:
	DECLARE_MESSAGE_MAP();     
		
	public:
	CDialZone(long *,long *,long *,long *,CWnd *pParent);     
	virtual BOOL OnInitDialog( );  	
	
	protected:
	//{{AFX_MSG(CDialZone) 
	afx_msg void OnOK( );
	afx_msg void OnCancel( );

	//}}AFX_MSG  
 	 
};


class CDialSeuil : public CDialog
{
	public:
	unsigned char * val; 
	BOOL *check;
	virtual BOOL OnInitDialog( );
		
	protected:
	DECLARE_MESSAGE_MAP(); 
	public:
	CDialSeuil(unsigned char  *,BOOL *);
	protected:  
	//{{AFX_MSG(CDialSeuil)
	afx_msg void OnOK( );
	afx_msg void OnCheck();  
	//}}AFX_MSG
};    

class CDialSeqlimite : public CDialog
{
	public:
	int * i1;
	int  * i2;
	int   depas;
	virtual BOOL OnInitDialog( );
		
	protected:
	DECLARE_MESSAGE_MAP(); 
	public:
	CDialSeqlimite(int  *,int *);     
	  	
	protected:
	//{{AFX_MSG(CDialSeqlimite)
	afx_msg void OnOK( );     
	//}}AFX_MSG
}; 
   

class CDialSeqCarac : public CDialog
{
	public:
	int * i1;
	unsigned char * i2;
	virtual BOOL OnInitDialog( );
		
	protected:
	DECLARE_MESSAGE_MAP(); 
	public:
	CDialSeqCarac(int *,unsigned char *);     
	//{{AFX_MSG(CDialSeqCarac)
	afx_msg void OnOK( );    
	//}}AFX_MSG
}; 

class CDialDoubleSeuil : public CDialog
{
	public:
	unsigned char * val1, * val2;
	virtual BOOL OnInitDialog( );
		
	protected:
	DECLARE_MESSAGE_MAP(); 
	public:
	CDialDoubleSeuil(unsigned int,unsigned char *,unsigned char *);     
	
	protected:
	//{{AFX_MSG(CDialDoubleSeuil)
	afx_msg void OnOK( );
	//}}AFX_MSG
}; 
#endif // __DIALOGUE_H__   