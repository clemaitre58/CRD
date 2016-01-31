#include "elemfun.h"

namespace utls
{
  // same as fill_triangle except without optimisation for filling big blocks
   template <typename AryBase, typename Point>
   void fill_small_triangle(AryBase *im, const Point &v1, const Point &v2, const Point &v3, typename AryBase::value fill)
   {
      // routine from: Nicolas Capens (http://www.devmaster.net/codespotlight/show.php?id=17)
      // 28.4 fixed-point coordinates 
      const int Y1 = int(16.0f * v1.y);
      const int Y2 = int(16.0f * v2.y);
      const int Y3 = int(16.0f * v3.y);

      const int X1 = int(16.0f * v1.x);
      const int X2 = int(16.0f * v2.x);
      const int X3 = int(16.0f * v3.x);

      // Deltas
      const int DX12 = X1 - X2;
      const int DX23 = X2 - X3;
      const int DX31 = X3 - X1;

      const int DY12 = Y1 - Y2;
      const int DY23 = Y2 - Y3;
      const int DY31 = Y3 - Y1;

      // Fixed-point deltas
      const int FDX12 = DX12 << 4;
      const int FDX23 = DX23 << 4;
      const int FDX31 = DX31 << 4;

      const int FDY12 = DY12 << 4;
      const int FDY23 = DY23 << 4;
      const int FDY31 = DY31 << 4;

      // Bounding rectangle
      int minx = (min(X1, X2, X3) + 0xF) >> 4;
      int maxx = (max(X1, X2, X3) + 0xF) >> 4;
      int miny = (min(Y1, Y2, Y3) + 0xF) >> 4;
      int maxy = (max(Y1, Y2, Y3) + 0xF) >> 4;

      // Check image boundaries
      minx = min(max(im->lb2, minx), im->ub2+1);
      maxx = max(min(maxx, im->ub2+1), im->lb2);

      miny = min(max(im->lb1, miny), im->ub1+1);
      maxy = max(min(maxy, im->ub1+1), im->lb1);   

      const int stride = im->cols();
      typename AryBase::pointer colorBuffer = &im->el[miny][0];//   (char*&)colorBuffer += miny * stride;

      // Half-edge constants
      int C1 = DY12 * X1 - DX12 * Y1;
      int C2 = DY23 * X2 - DX23 * Y2;
      int C3 = DY31 * X3 - DX31 * Y3;

      // Correct for fill convention
      if (DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
      if (DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
      if (DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;

      // Evaluate half-space functions
      int CY1 = C1 + DX12 * (miny<<4) - DY12 * (minx<<4);
      int CY2 = C2 + DX23 * (miny<<4) - DY23 * (minx<<4);
      int CY3 = C3 + DX31 * (miny<<4) - DY31 * (minx<<4);

      for(int y = miny; y < maxy; y++)
      {
         int CX1 = CY1;
         int CX2 = CY2;
         int CX3 = CY3;

         for(int x = minx; x < maxx; x++)
         {
            if(CX1 > 0 && CX2 > 0 && CX3 > 0)
               colorBuffer[x] = fill;

            CX1 -= FDY12;
            CX2 -= FDY23;
            CX3 -= FDY31;
         }
         CY1 += FDX12;
         CY2 += FDX23;
         CY3 += FDX31;
         colorBuffer += stride;
      }
   }

   template <typename AryBase, typename Value, typename Point>
   void accumulate_small_triangle(AryBase *im, const Point &v1, const Point &v2, const Point &v3, Value unit)
   {
      // routine from: Nicolas Capens (http://www.devmaster.net/codespotlight/show.php?id=17)
      // 28.4 fixed-point coordinates 
      const int Y1 = int(16.0f * v1.y);
      const int Y2 = int(16.0f * v2.y);
      const int Y3 = int(16.0f * v3.y);

      const int X1 = int(16.0f * v1.x);
      const int X2 = int(16.0f * v2.x);
      const int X3 = int(16.0f * v3.x);

      // Deltas
      const int DX12 = X1 - X2;
      const int DX23 = X2 - X3;
      const int DX31 = X3 - X1;

      const int DY12 = Y1 - Y2;
      const int DY23 = Y2 - Y3;
      const int DY31 = Y3 - Y1;

      // Fixed-point deltas
      const int FDX12 = DX12 << 4;
      const int FDX23 = DX23 << 4;
      const int FDX31 = DX31 << 4;

      const int FDY12 = DY12 << 4;
      const int FDY23 = DY23 << 4;
      const int FDY31 = DY31 << 4;

      // Bounding rectangle
      int minx = (min(X1, X2, X3) + 0xF) >> 4;
      int maxx = (max(X1, X2, X3) + 0xF) >> 4;
      int miny = (min(Y1, Y2, Y3) + 0xF) >> 4;
      int maxy = (max(Y1, Y2, Y3) + 0xF) >> 4;

      // Check image boundaries
      minx = min(max(im->lb2, minx), im->ub2+1);
      maxx = max(min(maxx, im->ub2+1), im->lb2);

      miny = min(max(im->lb1, miny), im->ub1+1);
      maxy = max(min(maxy, im->ub1+1), im->lb1);   

      const int stride = im->cols();
      typename AryBase::pointer colorBuffer = &im->el[miny][0];//   (char*&)colorBuffer += miny * stride;

      // Half-edge constants
      int C1 = DY12 * X1 - DX12 * Y1;
      int C2 = DY23 * X2 - DX23 * Y2;
      int C3 = DY31 * X3 - DX31 * Y3;

      // Correct for fill convention
      if (DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
      if (DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
      if (DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;

      // Evaluate half-space functions
      int CY1 = C1 + DX12 * (miny<<4) - DY12 * (minx<<4);
      int CY2 = C2 + DX23 * (miny<<4) - DY23 * (minx<<4);
      int CY3 = C3 + DX31 * (miny<<4) - DY31 * (minx<<4);

      for(int y = miny; y < maxy; y++)
      {
         int CX1 = CY1;
         int CX2 = CY2;
         int CX3 = CY3;

         for(int x = minx; x < maxx; x++)
         {
            if(CX1 > 0 && CX2 > 0 && CX3 > 0)            
               colorBuffer[x].push_back(unit);

            CX1 -= FDY12;
            CX2 -= FDY23;
            CX3 -= FDY31;
         }
         CY1 += FDX12;
         CY2 += FDX23;
         CY3 += FDX31;
         colorBuffer += stride;
      }
   }   

   template <typename AryBase, typename Point>
   void fill_triangle(AryBase *im, const Point &v1, const Point &v2, const Point &v3, typename AryBase::value fill)
   {
      // routine from: Nicolas Capens (http://www.devmaster.net/codespotlight/show.php?id=17)
      // 28.4 fixed-point coordinates 
      const int Y1 = int(16.0f * v1.y);
      const int Y2 = int(16.0f * v2.y);
      const int Y3 = int(16.0f * v3.y);

      const int X1 = int(16.0f * v1.x);
      const int X2 = int(16.0f * v2.x);
      const int X3 = int(16.0f * v3.x);

      // Deltas
      const int DX12 = X1 - X2;
      const int DX23 = X2 - X3;
      const int DX31 = X3 - X1;

      const int DY12 = Y1 - Y2;
      const int DY23 = Y2 - Y3;
      const int DY31 = Y3 - Y1;

      // Fixed-point deltas
      const int FDX12 = DX12 << 4;
      const int FDX23 = DX23 << 4;
      const int FDX31 = DX31 << 4;

      const int FDY12 = DY12 << 4;
      const int FDY23 = DY23 << 4;
      const int FDY31 = DY31 << 4;

      // Bounding rectangle
      int minx = (min(X1, X2, X3) + 0xF) >> 4;
      int maxx = (max(X1, X2, X3) + 0xF) >> 4;
      int miny = (min(Y1, Y2, Y3) + 0xF) >> 4;
      int maxy = (max(Y1, Y2, Y3) + 0xF) >> 4;

      // Check image boundaries
      minx = min(max(im->lb2, minx), im->ub2+1);
      maxx = max(min(maxx, im->ub2+1), im->lb2);

      miny = min(max(im->lb1, miny), im->ub1+1);
      maxy = max(min(maxy, im->ub1+1), im->lb1);     

      // Block size, standard 8x8 (must be power of two)
      const int q = 8;

      // Start in corner of 8x8 block
      minx &= ~(q - 1);
      miny &= ~(q - 1);
   
      const int stride = im->cols();
      typename AryBase::pointer colorBuffer = &im->el[miny][0];//   (char*&)colorBuffer += miny * stride;

      // Half-edge constants
      int C1 = DY12 * X1 - DX12 * Y1;
      int C2 = DY23 * X2 - DX23 * Y2;
      int C3 = DY31 * X3 - DX31 * Y3;

      // Correct for fill convention
      if(DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
      if(DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
      if(DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;

      // Loop through blocks
      for(int y = miny; y < maxy; y += q)
      {
         for(int x = minx; x < maxx; x += q)
         {
            // Corners of block
            int x0 = x << 4;
            int x1 = (x + q - 1) << 4;
            int y0 = y << 4;
            int y1 = (y + q - 1) << 4;

            // Evaluate half-space functions
            bool a00 = C1 + DX12 * y0 - DY12 * x0 > 0;
            bool a10 = C1 + DX12 * y0 - DY12 * x1 > 0;
            bool a01 = C1 + DX12 * y1 - DY12 * x0 > 0;
            bool a11 = C1 + DX12 * y1 - DY12 * x1 > 0;
            int a = (a00 << 0) | (a10 << 1) | (a01 << 2) | (a11 << 3);
    
            bool b00 = C2 + DX23 * y0 - DY23 * x0 > 0;
            bool b10 = C2 + DX23 * y0 - DY23 * x1 > 0;
            bool b01 = C2 + DX23 * y1 - DY23 * x0 > 0;
            bool b11 = C2 + DX23 * y1 - DY23 * x1 > 0;
            int b = (b00 << 0) | (b10 << 1) | (b01 << 2) | (b11 << 3);
    
            bool c00 = C3 + DX31 * y0 - DY31 * x0 > 0;
            bool c10 = C3 + DX31 * y0 - DY31 * x1 > 0;
            bool c01 = C3 + DX31 * y1 - DY31 * x0 > 0;
            bool c11 = C3 + DX31 * y1 - DY31 * x1 > 0;
            int c = (c00 << 0) | (c10 << 1) | (c01 << 2) | (c11 << 3);

            // Skip block when outside an edge
            if(a == 0x0 || b == 0x0 || c == 0x0) continue;

            typename AryBase::pointer buffer = colorBuffer;
            int limy = min(maxy, y+q); // don't go bellow the border
            int limx = min(maxx, x+q); // don't go bellow the border

            // Accept whole block when totally covered
            if(a == 0xF && b == 0xF && c == 0xF)
            {
               for(int iy = y; iy < limy; iy++)
               {
                  for(int ix = x; ix < limx; ix++)
                     buffer[ix] = fill;

                  buffer += stride;
               }
            }
            else // Partially covered block
            {
               int CY1 = C1 + DX12 * y0 - DY12 * x0;
               int CY2 = C2 + DX23 * y0 - DY23 * x0;
               int CY3 = C3 + DX31 * y0 - DY31 * x0;

               for(int iy = y; iy < limy; iy++)
               {
                  int CX1 = CY1;
                  int CX2 = CY2;
                  int CX3 = CY3;

                  for(int ix = x; ix < limx; ix++)
                  {
                     if(CX1 > 0 && CX2 > 0 && CX3 > 0)
                        buffer[ix] = fill;

                     CX1 -= FDY12;
                     CX2 -= FDY23;
                     CX3 -= FDY31;
                  }
                  CY1 += FDX12;
                  CY2 += FDX23;
                  CY3 += FDX31;
                  buffer += stride;
               }
            }
         }
         colorBuffer += q * stride;
      }
   }

   template <typename AryBase, typename Point>
   void add_triangle(AryBase *im, const Point &tv1, const Point &tv2, const Point &v3, typename AryBase::value value)
   {
      Point v1, v2;
      // negative => CW, positive => CCW
      if (((tv1.x-tv2.x)*(tv2.y-v3.y) - (tv1.y-tv2.y)*(tv2.x-v3.x))<0)
      {
         v1 = tv1; v2=tv2;
      } else {
         v1 = tv2; v2=tv1;
      }
         
      // routine from: Nicolas Capens (http://www.devmaster.net/codespotlight/show.php?id=17)
      // 28.4 fixed-point coordinates 
      const int Y1 = int(16.0f * v1.y);
      const int Y2 = int(16.0f * v2.y);
      const int Y3 = int(16.0f * v3.y);

      const int X1 = int(16.0f * v1.x);
      const int X2 = int(16.0f * v2.x);
      const int X3 = int(16.0f * v3.x);

      // Deltas
      const int DX12 = X1 - X2;
      const int DX23 = X2 - X3;
      const int DX31 = X3 - X1;

      const int DY12 = Y1 - Y2;
      const int DY23 = Y2 - Y3;
      const int DY31 = Y3 - Y1;

      // Fixed-point deltas
      const int FDX12 = DX12 << 4;
      const int FDX23 = DX23 << 4;
      const int FDX31 = DX31 << 4;

      const int FDY12 = DY12 << 4;
      const int FDY23 = DY23 << 4;
      const int FDY31 = DY31 << 4;

      // Bounding rectangle
      int minx = (min(X1, X2, X3) + 0xF) >> 4;
      int maxx = (max(X1, X2, X3) + 0xF) >> 4;
      int miny = (min(Y1, Y2, Y3) + 0xF) >> 4;
      int maxy = (max(Y1, Y2, Y3) + 0xF) >> 4;
      
      // Check image boundaries
      minx = min(max(im->lb2, minx), im->ub2+1);
      maxx = max(min(maxx, im->ub2+1), im->lb2);

      miny = min(max(im->lb1, miny), im->ub1+1);
      maxy = max(min(maxy, im->ub1+1), im->lb1);     
             

      // Block size, standard 8x8 (must be power of two)
      const int q = 8;

      // Start in corner of 8x8 block
      minx &= ~(q - 1);
      miny &= ~(q - 1);
   
      const int stride = im->cols();
      typename AryBase::pointer colorBuffer = &im->el[miny][0]; //   (char*&)colorBuffer += miny * stride;

      // Half-edge constants
      int C1 = DY12 * X1 - DX12 * Y1;
      int C2 = DY23 * X2 - DX23 * Y2;
      int C3 = DY31 * X3 - DX31 * Y3;

      // Correct for fill convention
      if(DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
      if(DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
      if(DY31 < 0 || (DY31 == 0 && DX31 > 0)) C3++;

      // Loop through blocks
      for(int y = miny; y < maxy; y += q)
      {
         for(int x = minx; x < maxx; x += q)
         {
            // Corners of block
            int x0 = x << 4;
            int x1 = (x + q - 1) << 4;
            int y0 = y << 4;
            int y1 = (y + q - 1) << 4;

            // Evaluate half-space functions
            bool a00 = C1 + DX12 * y0 - DY12 * x0 > 0;
            bool a10 = C1 + DX12 * y0 - DY12 * x1 > 0;
            bool a01 = C1 + DX12 * y1 - DY12 * x0 > 0;
            bool a11 = C1 + DX12 * y1 - DY12 * x1 > 0;
            int a = (a00 << 0) | (a10 << 1) | (a01 << 2) | (a11 << 3);
    
            bool b00 = C2 + DX23 * y0 - DY23 * x0 > 0;
            bool b10 = C2 + DX23 * y0 - DY23 * x1 > 0;
            bool b01 = C2 + DX23 * y1 - DY23 * x0 > 0;
            bool b11 = C2 + DX23 * y1 - DY23 * x1 > 0;
            int b = (b00 << 0) | (b10 << 1) | (b01 << 2) | (b11 << 3);
    
            bool c00 = C3 + DX31 * y0 - DY31 * x0 > 0;
            bool c10 = C3 + DX31 * y0 - DY31 * x1 > 0;
            bool c01 = C3 + DX31 * y1 - DY31 * x0 > 0;
            bool c11 = C3 + DX31 * y1 - DY31 * x1 > 0;
            int c = (c00 << 0) | (c10 << 1) | (c01 << 2) | (c11 << 3);

            // Skip block when outside an edge
            if(a == 0x0 || b == 0x0 || c == 0x0) continue;

            typename AryBase::pointer buffer = colorBuffer;

            int limy = min(maxy, y+q); // don't go bellow the border
            int limx = min(maxx, x+q); // don't go bellow the border

            // Accept whole block when totally covered
            if(a == 0xF && b == 0xF && c == 0xF)
            {
               for(int iy = y; iy < limy; iy++)
               {
                  for(int ix = x; ix < limx; ix++)
                     buffer[ix] += value;

                  buffer += stride;
               }
            }
            else // Partially covered block
            {
               int CY1 = C1 + DX12 * y0 - DY12 * x0;
               int CY2 = C2 + DX23 * y0 - DY23 * x0;
               int CY3 = C3 + DX31 * y0 - DY31 * x0;

               for(int iy = y; iy < limy; iy++)
               {
                  int CX1 = CY1;
                  int CX2 = CY2;
                  int CX3 = CY3;

                  for(int ix = x; ix < limx; ix++)
                  {
                     if(CX1 > 0 && CX2 > 0 && CX3 > 0)
                        buffer[ix] += value;

                     CX1 -= FDY12;
                     CX2 -= FDY23;
                     CX3 -= FDY31;
                  }
                  CY1 += FDX12;
                  CY2 += FDX23;
                  CY3 += FDX31;
                  buffer += stride;
               }
            }
         }
         colorBuffer += q * stride;
      }
   }

   template <typename AryBase, typename Point>
   void fill_quad(AryBase *im,
      const Point &v1, const Point &v2, const Point &v3, const Point &v4, 
      typename AryBase::value fill)
   {
      // routine from: Nicolas Capens (http://www.devmaster.net/codespotlight/show.php?id=17)
      // 28.4 fixed-point coordinates 
      const int Y1 = int(16.0f * v1.y);
      const int Y2 = int(16.0f * v2.y);
      const int Y3 = int(16.0f * v3.y);
      const int Y4 = int(16.0f * v4.y);

      const int X1 = int(16.0f * v1.x);
      const int X2 = int(16.0f * v2.x);
      const int X3 = int(16.0f * v3.x);
      const int X4 = int(16.0f * v4.x);

      // Deltas
      const int DX12 = X1 - X2;
      const int DX23 = X2 - X3;
      const int DX34 = X3 - X4;
      const int DX41 = X4 - X1;

      const int DY12 = Y1 - Y2;
      const int DY23 = Y2 - Y3;
      const int DY34 = Y3 - Y4;
      const int DY41 = Y4 - Y1;

      // Fixed-point deltas
      const int FDX12 = DX12 << 4;
      const int FDX23 = DX23 << 4;
      const int FDX34 = DX34 << 4;
      const int FDX41 = DX41 << 4;

      const int FDY12 = DY12 << 4;
      const int FDY23 = DY23 << 4;
      const int FDY34 = DY34 << 4;
      const int FDY41 = DY41 << 4;

      // Bounding rectangle
      int minx = (utls::min(X1, X2, X3, X4) + 0xF) >> 4;
      int maxx = (utls::max(X1, X2, X3, X4) + 0xF) >> 4;
      int miny = (utls::min(Y1, Y2, Y3, Y4) + 0xF) >> 4;
      int maxy = (utls::max(Y1, Y2, Y3, Y4) + 0xF) >> 4;

      // Check image boundaries
      minx = min(max(im->lb2, minx), im->ub2+1);
      maxx = max(min(maxx, im->ub2+1), im->lb2);

      miny = min(max(im->lb1, miny), im->ub1+1);
      maxy = max(min(maxy, im->ub1+1), im->lb1);     

      // Block size, standard 8x8 (must be power of two)
      const int q = 8;

      // Start in corner of 8x8 block
      minx &= ~(q - 1);
      miny &= ~(q - 1);
   
      const int stride = im->cols();
      typename AryBase::pointer colorBuffer = &im->el[miny][0];//   (char*&)colorBuffer += miny * stride;

      // Half-edge constants
      int C1 = DY12 * X1 - DX12 * Y1;
      int C2 = DY23 * X2 - DX23 * Y2;
      int C3 = DY34 * X3 - DX34 * Y3;
      int C4 = DY41 * X4 - DX41 * Y4;

      // Correct for fill convention
      if(DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
      if(DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
      if(DY34 < 0 || (DY34 == 0 && DX34 > 0)) C3++;
      if(DY41 < 0 || (DY41 == 0 && DX41 > 0)) C4++;

      // Loop through blocks
      for(int y = miny; y < maxy; y += q)
      {
         for(int x = minx; x < maxx; x += q)
         {
            // Corners of block
            int x0 = x << 4;
            int x1 = (x + q - 1) << 4;
            int y0 = y << 4;
            int y1 = (y + q - 1) << 4;

            // Evaluate half-space functions
            bool a00 = C1 + DX12 * y0 - DY12 * x0 > 0;
            bool a10 = C1 + DX12 * y0 - DY12 * x1 > 0;
            bool a01 = C1 + DX12 * y1 - DY12 * x0 > 0;
            bool a11 = C1 + DX12 * y1 - DY12 * x1 > 0;
            int a = (a00 << 0) | (a10 << 1) | (a01 << 2) | (a11 << 3);
    
            bool b00 = C2 + DX23 * y0 - DY23 * x0 > 0;
            bool b10 = C2 + DX23 * y0 - DY23 * x1 > 0;
            bool b01 = C2 + DX23 * y1 - DY23 * x0 > 0;
            bool b11 = C2 + DX23 * y1 - DY23 * x1 > 0;
            int b = (b00 << 0) | (b10 << 1) | (b01 << 2) | (b11 << 3);
    
            bool c00 = C3 + DX34 * y0 - DY34 * x0 > 0;
            bool c10 = C3 + DX34 * y0 - DY34 * x1 > 0;
            bool c01 = C3 + DX34 * y1 - DY34 * x0 > 0;
            bool c11 = C3 + DX34 * y1 - DY34 * x1 > 0;
            int c = (c00 << 0) | (c10 << 1) | (c01 << 2) | (c11 << 3);

            bool d00 = C4 + DX41 * y0 - DY41 * x0 > 0;
            bool d10 = C4 + DX41 * y0 - DY41 * x1 > 0;
            bool d01 = C4 + DX41 * y1 - DY41 * x0 > 0;
            bool d11 = C4 + DX41 * y1 - DY41 * x1 > 0;
            int d = (d00 << 0) | (d10 << 1) | (d01 << 2) | (d11 << 3);

            // Skip block when outside an edge
            if (a == 0x0 || b == 0x0 || c == 0x0 || d == 0x0) continue;

            typename AryBase::pointer buffer = colorBuffer;

            int limy = min(maxy, y+q); // don't go bellow the border
            int limx = min(maxx, x+q); // don't go bellow the border

            // Accept whole block when totally covered
            if(a == 0xF && b == 0xF && c == 0xF && d == 0xF)
            {
               for(int iy = y; iy < limy; iy++)
               {
                  for(int ix = x; ix < limx; ix++)
                     buffer[ix] = fill;
                  buffer += stride;
               }
            }
            else // Partially covered block
            {
               int CY1 = C1 + DX12 * y0 - DY12 * x0;
               int CY2 = C2 + DX23 * y0 - DY23 * x0;
               int CY3 = C3 + DX34 * y0 - DY34 * x0;
               int CY4 = C4 + DX41 * y0 - DY41 * x0;

               for(int iy = y; iy < limy; iy++)
               {
                  int CX1 = CY1;
                  int CX2 = CY2;
                  int CX3 = CY3;
                  int CX4 = CY4;

                  for(int ix = x; ix < limx; ix++)
                  {
                     if(CX1 > 0 && CX2 > 0 && CX3 > 0 && CX4 > 0)
                        buffer[ix] = fill;

                     CX1 -= FDY12;
                     CX2 -= FDY23;
                     CX3 -= FDY34;
                     CX4 -= FDY41;
                  }
                  CY1 += FDX12;
                  CY2 += FDX23;
                  CY3 += FDX34;
                  CY4 += FDX41;
                  buffer += stride;
               }
            }
         }
         colorBuffer += q * stride;
      }
   }

   template <typename AryBase, typename Point>
   void add_quad(AryBase *im,
      const Point &v1, const Point &v2, const Point &v3, const Point &v4, 
      typename AryBase::value value)
   {
      // routine from: Nicolas Capens (http://www.devmaster.net/codespotlight/show.php?id=17)
      // 28.4 fixed-point coordinates 
      const int Y1 = int(16.0f * v1.y);
      const int Y2 = int(16.0f * v2.y);
      const int Y3 = int(16.0f * v3.y);
      const int Y4 = int(16.0f * v4.y);

      const int X1 = int(16.0f * v1.x);
      const int X2 = int(16.0f * v2.x);
      const int X3 = int(16.0f * v3.x);
      const int X4 = int(16.0f * v4.x);

      // Deltas
      const int DX12 = X1 - X2;
      const int DX23 = X2 - X3;
      const int DX34 = X3 - X4;
      const int DX41 = X4 - X1;

      const int DY12 = Y1 - Y2;
      const int DY23 = Y2 - Y3;
      const int DY34 = Y3 - Y4;
      const int DY41 = Y4 - Y1;

      // Fixed-point deltas
      const int FDX12 = DX12 << 4;
      const int FDX23 = DX23 << 4;
      const int FDX34 = DX34 << 4;
      const int FDX41 = DX41 << 4;

      const int FDY12 = DY12 << 4;
      const int FDY23 = DY23 << 4;
      const int FDY34 = DY34 << 4;
      const int FDY41 = DY41 << 4;

      // Bounding rectangle
      int minx = (min(X1, X2, X3, X4) + 0xF) >> 4;
      int maxx = (max(X1, X2, X3, X4) + 0xF) >> 4;
      int miny = (min(Y1, Y2, Y3, Y4) + 0xF) >> 4;
      int maxy = (max(Y1, Y2, Y3, Y4) + 0xF) >> 4;

      // Check image boundaries
      minx = min(max(im->lb2, minx), im->ub2+1);
      maxx = max(min(maxx, im->ub2+1), im->lb2);

      miny = min(max(im->lb1, miny), im->ub1+1);
      maxy = max(min(maxy, im->ub1+1), im->lb1);     

      // Block size, standard 8x8 (must be power of two)
      const int q = 8;

      // Start in corner of 8x8 block
      minx &= ~(q - 1);
      miny &= ~(q - 1);
   
      const int stride = im->cols();
      typename AryBase::pointer colorBuffer = &im->el[miny][0];//   (char*&)colorBuffer += miny * stride;

      // Half-edge constants
      int C1 = DY12 * X1 - DX12 * Y1;
      int C2 = DY23 * X2 - DX23 * Y2;
      int C3 = DY34 * X3 - DX34 * Y3;
      int C4 = DY41 * X4 - DX41 * Y4;

      // Correct for fill convention
      if(DY12 < 0 || (DY12 == 0 && DX12 > 0)) C1++;
      if(DY23 < 0 || (DY23 == 0 && DX23 > 0)) C2++;
      if(DY34 < 0 || (DY34 == 0 && DX34 > 0)) C3++;
      if(DY41 < 0 || (DY41 == 0 && DX41 > 0)) C4++;

      // Loop through blocks
      for(int y = miny; y < maxy; y += q)
      {
         for(int x = minx; x < maxx; x += q)
         {
            // Corners of block
            int x0 = x << 4;
            int x1 = (x + q - 1) << 4;
            int y0 = y << 4;
            int y1 = (y + q - 1) << 4;

            // Evaluate half-space functions
            bool a00 = C1 + DX12 * y0 - DY12 * x0 > 0;
            bool a10 = C1 + DX12 * y0 - DY12 * x1 > 0;
            bool a01 = C1 + DX12 * y1 - DY12 * x0 > 0;
            bool a11 = C1 + DX12 * y1 - DY12 * x1 > 0;
            int a = (a00 << 0) | (a10 << 1) | (a01 << 2) | (a11 << 3);
    
            bool b00 = C2 + DX23 * y0 - DY23 * x0 > 0;
            bool b10 = C2 + DX23 * y0 - DY23 * x1 > 0;
            bool b01 = C2 + DX23 * y1 - DY23 * x0 > 0;
            bool b11 = C2 + DX23 * y1 - DY23 * x1 > 0;
            int b = (b00 << 0) | (b10 << 1) | (b01 << 2) | (b11 << 3);
    
            bool c00 = C3 + DX34 * y0 - DY34 * x0 > 0;
            bool c10 = C3 + DX34 * y0 - DY34 * x1 > 0;
            bool c01 = C3 + DX34 * y1 - DY34 * x0 > 0;
            bool c11 = C3 + DX34 * y1 - DY34 * x1 > 0;
            int c = (c00 << 0) | (c10 << 1) | (c01 << 2) | (c11 << 3);

            bool d00 = C4 + DX41 * y0 - DY41 * x0 > 0;
            bool d10 = C4 + DX41 * y0 - DY41 * x1 > 0;
            bool d01 = C4 + DX41 * y1 - DY41 * x0 > 0;
            bool d11 = C4 + DX41 * y1 - DY41 * x1 > 0;
            int d = (d00 << 0) | (d10 << 1) | (d01 << 2) | (d11 << 3);

            // Skip block when outside an edge
            if (a == 0x0 || b == 0x0 || c == 0x0 || d == 0x0) continue;

            typename AryBase::pointer buffer = colorBuffer;

            int limy = min(maxy, y+q); // don't go bellow the border
            int limx = min(maxx, x+q); // don't go bellow the border

            // Accept whole block when totally covered
            if(a == 0xF && b == 0xF && c == 0xF && d == 0xF)
            {
               for(int iy = y; iy < limy; iy++)
               {
                  for(int ix = x; ix < limx; ix++)
                     buffer[ix] += value;
                  buffer += stride;
               }
            }
            else // Partially covered block
            {
               int CY1 = C1 + DX12 * y0 - DY12 * x0;
               int CY2 = C2 + DX23 * y0 - DY23 * x0;
               int CY3 = C3 + DX34 * y0 - DY34 * x0;
               int CY4 = C4 + DX41 * y0 - DY41 * x0;

               for(int iy = y; iy < limy; iy++)
               {
                  int CX1 = CY1;
                  int CX2 = CY2;
                  int CX3 = CY3;
                  int CX4 = CY4;

                  for(int ix = x; ix < limx; ix++)
                  {
                     if(CX1 > 0 && CX2 > 0 && CX3 > 0 && CX4 > 0)
                        buffer[ix] += value;

                     CX1 -= FDY12;
                     CX2 -= FDY23;
                     CX3 -= FDY34;
                     CX4 -= FDY41;
                  }
                  CY1 += FDX12;
                  CY2 += FDX23;
                  CY3 += FDX34;
                  CY4 += FDX41;
                  buffer += stride;
               }
            }
         }
         colorBuffer += q * stride;
      }
   }

}
