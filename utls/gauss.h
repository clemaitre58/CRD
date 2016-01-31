#ifndef _UTLS__GAUSS_H_
#define _UTLS__GAUSS_H_

#include <assert.h>
#include <elemfun.h>
#include <cmath>
#include <stdio.h>

namespace utls
{
   const int MAX_KERNEL_SIZE = 1024;

   template<typename PixelType, typename ValueType>
   class GaussianKernel
   {
   public:
      typedef ValueType _ValueType;
      typedef PixelType _PixelType;
      //static const ValueType GaussTruncate = 3.0;
      GaussianKernel(double sigma)
         {
            int i; 
            
            /* calculate proper kernel size - odd, bigger than 3 */
            size = (int)(2.0 * 3.0 * sigma + 1.0);
            if (size<3) size=3;
            if (size % 2 == 0) size++;
            assert(size < MAX_KERNEL_SIZE);
            /* fill in kernel values */
            double x;
            ValueType sum = 0.0;
            /* Fill in kernel values. */
            for (i = 0; i < size; i++) 
            {
               x = i - size / 2;
               kernel[i] = ::exp(double(- x * x / (2.0 * sigma * sigma)));
               sum += kernel[i];
            }
            /* Normalize kernel values to sum to 1.0. */
            for (i = 0; i < size; i++)
               kernel[i] /= ValueType(sum);
         }
   
      /* fast convolution with kernel of size 3 */
      ValueType Convolve3(const PixelType *buffer) const
         {
            assert(size==3);
            // exploit symmetry
            return 
               kernel[0]*ValueType(buffer[0] + buffer[2]) +
               kernel[1]*ValueType(buffer[1]);
         }

      /* fast convolution with kernel of size 5 */
      ValueType Convolve5(const PixelType *buffer) const
         {
            assert(size==5);
            // exploit symmetry
            return 
               kernel[0]*ValueType(buffer[0] + buffer[4]) +
               kernel[1]*ValueType(buffer[1] + buffer[3]) +
               kernel[2]*ValueType(buffer[2]);
         }

      /* fast convolution with kernel of size 7 */
      ValueType Convolve7(const PixelType *buffer) const
         {
            assert(size==7);
            // exploit symmetry
            return (ValueType)
               kernel[0]*ValueType(buffer[0] + buffer[6]) +
               kernel[1]*ValueType(buffer[1] + buffer[5]) +
               kernel[2]*ValueType(buffer[2] + buffer[4]) +
               kernel[3]*ValueType(buffer[3]);
         }
     
   public:
      ValueType kernel[MAX_KERNEL_SIZE];
      int size;
   };
}

#endif // _UTLS__GAUSS_H_
