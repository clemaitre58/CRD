#ifndef _UTLS__CONV_H_
#define _UTLS__CONV_H_

#include <assert.h>
#include <elemfun.h>
#include <cmath>
#include <stdio.h>

namespace utls
{
   const int MAX_CONV_BUFFER = 65536;
   /* Same as ConvBuffer, but implemented with loop unrolling for increased
      speed.  This is the most time intensive routine in keypoint detection,
      so deserves careful attention to efficiency.  Loop unrolling simply
      sums 5 multiplications at a time to allow the compiler to schedule
      operations better and avoid loop overhead.  This almost triples
      speed of previous version on a Pentium with gcc. */
   template <typename PixelType, class KernelBase>
   void convolve_buffer_fast(PixelType *buffer, const KernelBase &k, int rsize)
   {
      PixelType *bp;
      typedef typename KernelBase::_ValueType ValueType;
      ValueType sum; const ValueType *kp, *endkp;
         
      if (k.size<9)
      {
         // use kernel specific convolutions for small sizes
         switch (k.size)
         {
         case 3:
            for (int i = 0; i < rsize; i++) 
               buffer[i] = (PixelType)k.Convolve3(buffer+i);
            break;
         case 5:
            for (int i = 0; i < rsize; i++) 
               buffer[i] = (PixelType)k.Convolve5(buffer+i);
            break;
         case 7:
            for (int i = 0; i < rsize; i++) 
               buffer[i] = (PixelType)k.Convolve7(buffer+i);
            break;
         }
      } else {
         for (int i = 0; i < rsize; i++) 
         {
            sum = (ValueType)0.0;
            /* Do 11 multiplications at a time on remaining items. */
            bp = &buffer[i]; kp = &k.kernel[0]; endkp = &k.kernel[k.size];
            while (kp + 10 < endkp)
            {
               sum += ValueType(bp[0]) * kp[0] +  ValueType(bp[1]) * kp[1] + 
                  ValueType(bp[2]) * kp[2] + ValueType(bp[3]) * kp[3] +  ValueType(bp[4]) * kp[4] +
                  ValueType(bp[5]) * kp[5] + ValueType(bp[6]) * kp[6] +  ValueType(bp[7]) * kp[7] +
                  ValueType(bp[8]) * kp[8] + ValueType(bp[9]) * kp[9] +  ValueType(bp[10]) * kp[10];
               bp += 11; kp += 11;
            }
            /* Do 5 multiplications at a time on remaining items. */
            while (kp + 4 < endkp)
            {
               sum += ValueType(bp[0]) * kp[0] +  ValueType(bp[1]) * kp[1] + 
                  ValueType(bp[2]) * kp[2] + ValueType(bp[3]) * kp[3] +  ValueType(bp[4]) * kp[4];
               bp += 5; kp += 5;
            }
            /* Do 2 multiplications at a time on remaining items. */
            while (kp + 1 < endkp) 
            {
               sum += ValueType(bp[0]) * kp[0] +  ValueType(bp[1]) * kp[1];
               bp += 2; kp += 2;
            }
            /* Finish last one if needed. */
            if (kp < endkp)
               sum += ValueType(*bp) * *kp;
            buffer[i] = (PixelType)sum;
         }
      }
   }

   template <class PixelType, class KernelBase>
   void conv_horizontal(const PixelType **src, PixelType **dst, int rows, int cols, const KernelBase &k)
   {
      int r, c, i, halfsize;
      PixelType buffer[MAX_CONV_BUFFER];
      halfsize = k.size / 2;
      assert(cols + k.size < MAX_CONV_BUFFER);

      for (r = 0; r < rows; r++) 
      {
         /* Copy the row into buffer with pixels at ends replicated for
            half the mask size.  This avoids need to check for ends
            within inner loop. */
         for (i = 0; i < halfsize; i++)
            buffer[i] = src[r][0];
         for (i = 0; i < cols; i++)
            buffer[halfsize + i] = src[r][i];
         for (i = 0; i < halfsize; i++)
            buffer[halfsize + cols + i] = src[r][cols - 1];

         convolve_buffer_fast(buffer, k, cols);
         for (c = 0; c < cols; c++)
            dst[r][c] = buffer[c];
      }
   }

   /* Same as ConvHorizontal, but apply to vertical columns of image. */
   template <class PixelType, class KernelBase>
   void conv_vertical(const PixelType **src, PixelType **dst, int rows, int cols, const KernelBase &k)
   {
      int r, c, i, halfsize;
      PixelType buffer[MAX_CONV_BUFFER];
      halfsize = k.size / 2;
      assert(rows + k.size < MAX_CONV_BUFFER);

      for (c = 0; c < cols; c++) 
      {
         for (i = 0; i < halfsize; i++)
            buffer[i] = src[0][c];
         for (i = 0; i < rows; i++)
            buffer[halfsize + i] = src[i][c];
         for (i = 0; i < halfsize; i++)
            buffer[halfsize + rows + i] = src[rows - 1][c];

         convolve_buffer_fast(buffer, k, rows);
         for (r = 0; r < rows; r++)
            dst[r][c] = buffer[r];
      }
   }
}
#endif // _UTLS__CONV_H_
