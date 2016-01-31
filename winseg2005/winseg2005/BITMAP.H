#include "stdafx.h"
#include "struct.h"

//#include <windows.h>

	
void SetInfoBitmap(BITMAPINFO * info,long L,long H,NBCOLOR NBRGB);

void SetInfoBitmap(BITMAPINFO * infoDest,long L,long H,BITMAPINFO * infoSource);

void SetTeteBitmap(BITMAPFILEHEADER * tete,long L,long H,NBCOLOR NBRGB);
