#include <windows.h>
#include "struct.h"

static HINSTANCE hInstDLL;

int FAR PASCAL LibMain(HANDLE hInstance,WORD wDataSeg,WORD cbHeapSize,LPSTR lpszCmdLine)
{
 wDataSeg;
 lpszCmdLine;
 if (cbHeapSize>0) UnlockData(0);
 hInstDLL = hInstance;
 return TRUE;
}
