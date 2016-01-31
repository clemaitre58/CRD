#ifndef __UTLS_MARRAY_H__
#define __UTLS_MARRAY_H__

#include "mex.h"
#include <assert.h>
namespace utls
{
   /* get matlab_class id from type */
  	template<typename T>
   mxClassID mxClassIDFromType(T *array)
   { 
      assert("Unknown type for mxClassFromType."); 
      return mxClassID(-1);
   };

   template<> mxClassID mxClassIDFromType(unsigned char *array) 
   { return mxClassID(mxUINT8_CLASS);}

   template<> mxClassID mxClassIDFromType(unsigned int *array) 
   { return mxClassID(mxUINT32_CLASS);}

   template<> mxClassID mxClassIDFromType(signed char *array) 
   { return mxClassID(mxINT8_CLASS);}

   template<> mxClassID mxClassIDFromType(int *array) 
   { return mxClassID(mxINT32_CLASS);}

   template<> mxClassID mxClassIDFromType(float *array)
   { return mxClassID(mxSINGLE_CLASS);}

   template<> mxClassID mxClassIDFromType(double *array) 
   { return mxClassID(mxDOUBLE_CLASS);}

   template<typename T> 
   struct matlab_array
   {
      T *ptr;
      mxArray *src;

      matlab_array(const mxArray *a, bool convert = false)
         {
            /* check if passed array has excepted type */
            if (!convert)
               assert(mxClassIDFromType(ptr)==mxGetClassID(a));
            src = const_cast<mxArray *>(a);
            reset();
         }
      
      matlab_array(mwSize num_of_dims, mwSize *dims)
         { 
            ptr = 0;
            src = mxCreateNumericArray(num_of_dims, dims, 
                                       mxClassIDFromType(ptr), mxREAL);
            reset();
         }

      inline void reset()
         {
            ptr = (T *)mxGetData(src); 
         }

      inline void store(T value)
         {
            *ptr++ = value;
         }

      inline T read()
         {
            return *ptr++;
         }

      mxArray *array() { return src; }      
   };

   /* convert provided 2d array to matlab matrix of height x width elements
      of type item_class */
   template<typename PixelType, typename MatlabType>
   mxArray *array2d_to_matlab(const PixelType *image, size_t width, 
                              size_t height, size_t phy_width=0)
   {
      mwSize dims[2] = { height, width };
      /* adapt correct result data type */
      matlab_array<MatlabType> result(2, dims);
      if (!phy_width)
         phy_width=width;
      for (size_t x=0; x<width; x++)
         for (size_t y=0; y<height; y++)
            result.store((MatlabType)(image[x + phy_width * y]));
      return result.array();
   }

   template<typename MatlabType, typename PixelType>
   void array2d_from_matlab(const mxArray *src, PixelType **image, 
                            size_t &width, size_t &height, size_t phy_width=0)
   {
      assert(mxGetNumberOfDimensions(src)==2);
      const mwSize *dims = mxGetDimensions(src);
      matlab_array<MatlabType> result(src);      
      height = dims[0]; width = dims[1]; 
      if (!phy_width)
         phy_width=width;
      *image = new PixelType[phy_width*height];
      for (size_t x=0; x<width; x++)
         for (size_t y=0; y<height; y++)
            (*image)[x+phy_width*y] = (PixelType)result.read();
   }

   template<typename PixelType, typename MatlabType>
   mxArray *array3d_to_matlab(const PixelType *image, size_t width, 
                              size_t height, size_t depth)
   {
      // convert to matlab's column based representation
      mwSize dims[3] = { height, width, depth };
      matlab_array<MatlabType> result(3, dims);
      for (size_t z=0; z<depth; z++)
         for (size_t x=0; x<width; x++)
            for (size_t y=0; y<height; y++)
               result.store((MatlabType)image[z + depth * (x + width * y)]);
      return result.array();
   }

   template<typename MatlabType, typename PixelType>
   void array3d_from_matlab(const mxArray *src, PixelType **image, 
                            size_t& width, size_t& height, size_t& depth)
   {
      assert(mxGetNumberOfDimensions(src)==3);
      const mwSize *dims = mxGetDimensions(src);
      matlab_array<MatlabType> result(src);      
      height = dims[0]; width = dims[1]; depth = dims[2];
      *image = new PixelType[width*height*depth];
      for (size_t z=0; z<depth; z++)
         for (size_t x=0; x<width; x++)
            for (size_t y=0; y<height; y++)
            {
               (*image)[z + depth * (x + width * y)] = 
                  (PixelType)result.read();
            }
   }
}

#endif // __UTLS_MARRAY_H__
