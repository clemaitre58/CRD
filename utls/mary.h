#ifndef __UTLS_MARY_H__
#define __UTLS_MARY_H__

#include <mex.h>
#include "marray.h"

namespace utls
{
   /* convert provided Ary (ary.h) to matlab matrix of height x width 
      elements of type item_class */
   template<typename VectorClass, typename MatlabType>
   mxArray *vec_to_matlab(const VectorClass *vec, size_t items)
   {
      mwSize dims[2] = { 1, items };
      matlab_array<MatlabType> result(2, dims);
      for (size_t i=0; i<items; i++)
         result.store((MatlabType)(vec[i]));
      return result.array();
   }

   template<typename MatlabType, typename VectorClass>
   VectorClass *vec_from_matlab(const mxArray *src)
   {
      assert(mxGetNumberOfDimensions(src)==2);
      const mwSize *dims = mxGetDimensions(src); size_t items = dims[0]*dims[1];
      matlab_array<MatlabType> tmp(src, true);
      VectorClass *result = new VectorClass[items];
      for (size_t i=0; i<items; i++)
         result[i] = (VectorClass)tmp.read();
      return result;
   }

   template<typename AryClass, typename MatlabType>
   mxArray *ary_to_matlab(const AryClass *ary, int lb1, int ub1, int lb2, int ub2)
   {
      mwSize dims[2] = { ub1-lb1+1, ub2-lb2+1 };
      /* adapt correct result data type */
      matlab_array<MatlabType> result(2, dims);
      for (int x=lb2; x<=ub2; x++)
         for (int y=lb1; y<=ub1; y++)
            result.store((MatlabType)(ary->el[y][x]));
      return result.array();
   }

   template<typename AryClass, typename MatlabType>
   mxArray *ary_to_matlab(const AryClass *ary)
   {
      return ary_to_matlab<AryClass, MatlabType>(ary, ary->lb1, ary->ub1, ary->lb2, ary->ub2);
   }

   template<typename MatlabType, typename AryClass>
   AryClass *ary_from_matlab(const mxArray *src)
   {
      assert(mxGetNumberOfDimensions(src)==2);
      const mwSize *dims = mxGetDimensions(src);
      matlab_array<MatlabType> tmp(src, true);
      AryClass *result;
      result = new AryClass(0, dims[0]-1, 0, dims[1]-1);
      for (int x=result->lb2; x<=result->ub2; x++)
         for (int y=result->lb1; y<=result->ub1; y++)
            result->el[y][x] = (typename AryClass::value)tmp.read();
      return result;
   }

   template<typename MatlabType, typename AryClass>
   AryClass **array3d_from_matlab(const mxArray *src)
   {
      assert(mxGetNumberOfDimensions(src)==3);
      const mwSize *dims = mxGetDimensions(src);
      matlab_array<MatlabType> tmp(src, true);
		typedef AryClass* AryClassPtr;
      AryClass **result = new AryClassPtr[dims[2]];
      for (mwSize z=0; z<dims[2]; z++)
      {
         result[z] = new AryClass(0, dims[0]-1, 0, dims[1]-1);
         for (int x=result[z]->lb2; x<=result[z]->ub2; x++)
            for (int y=result[z]->lb1; y<=result[z]->ub1; y++)
               result[z]->el[y][x] = (typename AryClass::value)tmp.read();
      }
      return result;
   }

   template<typename MatlabType, typename AryClass>
   AryClass *array3d_from_matlab_interleaved(const mxArray *src)
   {
      assert(mxGetNumberOfDimensions(src)==3);
      const mwSize *dims = mxGetDimensions(src);
      matlab_array<MatlabType> tmp(src, true);
      AryClass *result = new AryClass(0, dims[0]-1, 0, dims[1]-1);
      result = new AryClass(0, dims[0]-1, 0, dims[1]-1);
      for (mwSize z=0; z<dims[2]; z++)
      {
         for (int x=result->lb2; x<=result->ub2; x++)
            for (int y=result->lb1; y<=result->ub1; y++)
            {
               result->el[y][x][z] = (typename AryClass::value::value)tmp.read();
            }
      }
      return result;
   }

   template<typename AryClass, typename MatlabType>
   mxArray *array3d_to_matlab(const AryClass **ary, int lb1, int ub1, int lb2, int ub2, int lb3, int ub3)
   {
      mwSize dims[3]; 
      dims[0] = ub1-lb1+1;
      dims[1] = ub2-lb2+1; 
      dims[2] = ub3-lb3+1;
      /* adapt correct result data type */
      matlab_array<MatlabType> result(3, dims);
      for (int z=lb3; z<=ub3; z++)
         for (int x=lb2; x<=ub2; x++)
            for (int y=lb1; y<=ub1; y++)
               result.store((MatlabType)(ary[z]->el[y][x]));
      return result.array();
   }

   template<typename AryClass, typename MatlabType>
   mxArray *array3d_to_matlab(const AryClass **ary, int rows, int cols, int channels)
   {
      mwSize dims[3]; 
      dims[0] = rows;
      dims[1] = cols; 
      dims[2] = channels;
      /* adapt correct result data type */
      matlab_array<MatlabType> result(3, dims);
      for (int z=0; z<channels; z++)
         for (int x=0; x<cols; x++)
            for (int y=0; y<rows; y++)
               result.store((MatlabType)(ary[z]->el[y][x]));
      return result.array();
   }
}

#endif // __UTLS_MARY_H__
