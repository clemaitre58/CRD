#ifndef _NEWARY_H_
#define _NEWARY_H_

#include <elemfun.h>
#include <cmath>
#include <stdio.h>
#include "dtypes.h"
#include "gauss.h"
#include "conv.h"

namespace utls
{
   template <
      typename PixelType,                        // storage type
      typename ScalarFloatingType = PixelType,   // type for precise math scalar operations (scaling, computing gradient...)
      typename PreciseCoordType   = float>       // precise coordinates type, used in scaling and bilinear interpolation of pixels
   struct Ary
   {
   private:
      typedef Ary<PixelType, ScalarFloatingType, PreciseCoordType> AryBase;

   public:   
      typedef PixelType   value;
      typedef PixelType*  pointer;
      typedef PixelType** row_pointer;
      typedef int         size_type;
      typedef int         coord_type;
   

      typedef ScalarFloatingType scalar_float;
      typedef PreciseCoordType   precise_coord;

      /* construction */
      Ary()
         {
            lb1=lb2=ub1=ub2=0; 
            num_rows=num_cols=0;
            data=0;
         }
      Ary(coord_type fr, coord_type lr,
          coord_type fc, coord_type lc, pointer attach_to = 0)
         {
            cons(fr, lr, fc, lc, attach_to);
         }

      Ary(coord_type nrows, coord_type ncols, pointer attach_to = 0)
         {
            cons(0, nrows-1, 0, ncols-1, attach_to);
         }

      Ary(const Ary &other, bool do_copy=true, bool do_attach = false)
         {
            if (!do_attach) {
               cons(other.lb1, other.ub1, other.lb2, other.ub2, 0);
               if (do_copy)
               copy(&other);
            } else
               cons(other.lb1, other.ub1, other.lb2, other.ub2, other.data);
         }

      ~Ary()
         {
            deallocate();
         }

      /* basics */
      // pointer   ptr() const { return data; }
      pointer   first() const { return data; }
      pointer   last() const { return data+num_rows*num_cols; }
      size_type rows() const { return num_rows; }
      size_type cols() const { return num_cols; }
      size_type size() const { return num_rows*num_cols*sizeof(value); }
      bool isin(coord_type row, coord_type col) const 
         {
            return 
               row>=lb1 && row<=ub1 && col>=lb2 && col<=ub2;
         }

      pointer operator[](int index) { return el[index]; }
      const pointer operator[](int index) const { return el[index]; }

      /* allocation */
      void cons(coord_type firstrow,
                coord_type lastrow, 
                coord_type firstcol, 
                coord_type lastcol, 
                pointer    attach_to)
         {
            num_rows = lastrow - firstrow + 1;         
            num_cols = lastcol - firstcol + 1;
            lb1 = firstrow; ub1 = lastrow;
            lb2 = firstcol; ub2 = lastcol;
            el = new pointer[num_rows]; el = el - firstrow;
            if (!attach_to)
               data = new value[num_rows*num_cols];
            else
               data = attach_to;
            pointer mem = data - firstcol;
            for (coord_type r = firstrow; r <= lastrow; ++r)
            {
               el[r] = mem;            
               mem = mem + num_cols;
            }
         }

      void clear()
         {
            for (int i = 0; i < num_rows*num_cols; i++)
               data[i]=(PixelType)0;
         }

      void clear(int r1, int r2, int c1, int c2)
         {
            for (int r = r1; r<=r2; r++)
            {
               for (int c = c1; c<=c2; c++)
                  el[r][c] = 0;
            }
         }

      void set(const PixelType &value)
         {
            for (int i = 0; i < num_rows*num_cols; i++)
               data[i]=value;
         }
   
      void detach()
         {
            if (el) 
            {
               el = el + lb1;
               delete [] el; 
            }
            el = 0;
            data = 0; 
            num_rows = num_cols = lb1 = lb2 = ub1 = ub2 = 0;
         }

      /* create a copy of an image */
      Ary* copy() const
         {
            AryBase *newary = new AryBase(lb1, ub1, lb1, ub2);
            for (int i = 0; i < num_rows*num_cols; i++)
               newary->data[i] = data[i];
            return newary;
         }

      /* create a copy of an image */
      void copy(const Ary *from)
         {
            for (int i = 0; i < num_rows*num_cols; i++)
            data[i] = from->data[i];
         }

      /* create a copy of an image */
      void copy(const Ary *s, coord_type top, coord_type left)
      {
         int sr, sc, dr, dc;
         PixelType **src = s->el;

         /* it should fit in */
         assert(s->num_rows + top <= num_rows);
         assert(s->num_cols + left <= num_cols);

         for (sr = s->lb1, dr = top; sr <= s->ub1; sr++,dr++)
         {
            pointer src = s->el[sr], dst = el[dr];
            for (sc = s->lb2, dc = left; sc <= s->ub2; sc++, dc++) 
               dst[dc] = src[sc];
         }
         }
   
      /* scale image with a given value */
      void scale(ScalarFloatingType value)
         {
            for (int i = 0; i < num_rows*num_cols; i++)
               data[i] *= value;
         }

      /* scale pixelwise with image im2  and leave result in this image. */
      void scale(Ary *other, ScalarFloatingType factor)
         {
            PixelType *other_data = other->data;
            for (int i = 0; i < num_rows*num_cols; i++)
               data[i] += factor*other_data[i];
         }

      /* offset image by scalar value */
      void add(PixelType value)
         {
            for (int i = 0; i < num_rows*num_cols; i++)
               data[i] += value;
         }

      /* add other image */
      void add(Ary *other)
         {
            PixelType *other_data = other->data;
            for (int i = 0; i < num_rows*num_cols; i++)
               data[i] += other_data[i];
         }
   
      /* offset image by scalar value */
      void subtract(PixelType value)
         {
            for (int i = 0; i < num_rows*num_cols; i++)
               data[i] -= value;
         }

      /* subtract other image */
      void subtract(Ary *other)
         {
            PixelType *other_data = other->data;
            for (int i = 0; i < num_rows*num_cols; i++)
               data[i] -= other_data[i];
         }

      /* add other image */
      void sqrt()
         {
            for (int i = 0; i < num_rows*num_cols; i++)
               data[i] = ::sqrt(data[i]);
         }
   
      void deallocate()
         {
            if (data) 
               delete [] data;
            detach();
         }
  
      void write_pgm(char *filename)
         {
            FILE *f = fopen(filename, "wb+");
            fprintf(f,"P5\n%d %d\n255\n", num_cols, num_rows);
            fclose(f);
            f = fopen(filename, "ab");
            t_byte *buf = new t_byte[num_cols];
            for (int r = lb1; r <= ub1; r++)
            {
               pointer p = el[r]; t_byte *tmp = buf;
               for (int c = lb2; c <= ub2; c++) 
                  (*tmp++) = (t_byte)(*p++);
               fwrite(buf, sizeof(t_byte), num_cols, f);
            }
            fclose(f);
            delete[]buf;
         }

   public:
      coord_type    lb1, lb2, ub1, ub2;
      size_type     num_rows, num_cols;
      pointer       data;
      row_pointer   el;
      /* leave some space for aasociated user's variable */
      int           tag;
      void          *user_data;
   };

   /* one channel */
   typedef Ary<unsigned char>              BAry;
   typedef Ary<int>                        IAry;
   typedef Ary<long unsigned int>          UI64Ary;
   typedef Ary<unsigned int>               LAry;
   typedef Ary<float>                      FAry;
   typedef Ary<double,double,double>       DAry;
   typedef Ary<void *>                     PAry;

   /* three channels */
   typedef Ary<RGBValue<t_byte>, t_byte, RGBValue<float> > BRGBAry;
   typedef Ary<RGBValue<float>, float, RGBValue<float> >   FRGBAry;


   // works only for 0 based Arys
   template<typename AryBase, typename AccumulatorType>
   void gaussian_blur(AryBase *img, double sigma)
   {
      typedef typename AryBase::value PixelType;
      GaussianKernel<PixelType, AccumulatorType> k(sigma);
      conv_horizontal<PixelType>(const_cast<const PixelType **>(img->el), img->el, img->rows(), img->cols(), k);
      conv_vertical<PixelType>(const_cast<const PixelType **>(img->el), img->el, img->rows(), img->cols(), k);
   }

   // works only for 0 based Arys
   template<typename AryBase, typename AccumulatorType>
   AryBase *copy_and_gaussian_blur(const AryBase *img, double sigma)
   {
      typedef typename AryBase::value PixelType;
      AryBase *newimg = new AryBase(img->rows(), img->cols());
      GaussianKernel<PixelType, AccumulatorType> k(sigma);
      conv_horizontal(const_cast<const PixelType **>(img->el), newimg->el, img->rows(), img->cols(), k);
      conv_vertical(const_cast<const PixelType **>(newimg->el), newimg->el, img->rows(), img->cols(), k);
      return newimg;
   }


   template<typename AryBase1, typename AryBase2>
   AryBase2 *convert(AryBase1 *src)
   {
      AryBase2 *dst = new AryBase2(src->lb1, src->ub1, src->lb2, src->ub2, 0);
      for (int r = src->lb1; r <= src->ub1; r++)
      {
         typename AryBase1::pointer s = src->el[r];
         typename AryBase2::pointer d = dst->el[r];
         for (int c = src->lb2; c <= src->ub2; c++) 
            d[c]=(typename AryBase2::value)(s[c]);
      }
      return dst;
   }

   template<typename AryBase1, typename AryBase2>
   void convert(AryBase1 *src, AryBase2 *dst)
   {
      for (int r = src->lb1; r <= src->ub1; r++)
      {
         typename AryBase1::pointer s = src->el[r];
         typename AryBase2::pointer d = dst->el[r];
         for (int c = src->lb2; c <= src->ub2; c++) 
            d[c]=(typename AryBase2::value)(s[c]);
      }     
   }

   template <typename rgb, typename intensity>
   void convert_rgb_to_intensity(rgb &s, intensity &ret)
   {
      ret = intensity((double(s.r)+s.g+s.b)/3.0);
   }

   template<typename AryBase1, typename AryBase2>
   void convert_to_intensity(AryBase1 *src, AryBase2 *dst)
   {
      for (int r = src->lb1; r <= src->ub1; r++)
      {
         typename AryBase1::pointer s = src->el[r];
         typename AryBase2::pointer d = dst->el[r];
         for (int c = src->lb2; c <= src->ub2; c++)
         {
            typename AryBase2::value v; convert_rgb_to_intensity(s[c], v);
            d[c]=v;
         }
      }     
   }

   template<typename AryBase1, typename AryBase2>
   AryBase2 *convert_to_intensity(AryBase1 *src)
   {
      AryBase2 *dst = new AryBase2(src->lb1, src->ub1, src->lb2, src->ub2, 0);
      for (int r = src->lb1; r <= src->ub1; r++)
      {
         typename AryBase1::pointer s = src->el[r];
         typename AryBase2::pointer d = dst->el[r];
         for (int c = src->lb2; c <= src->ub2; c++)
         {
            typename AryBase2::value v; convert_rgb_to_intensity(s[c], v);
            d[c]=v;
         }
      }
      return dst;
   }

   template<typename PixelType, typename ReductionFunctor >
   PixelType reduce_ary(const Ary<PixelType, PixelType> *a, const ReductionFunctor &func, PixelType init_value = 0)
   {
      PixelType value = init_value; int size = a->rows() * a->cols();
      for (int i = 0; i < size; i++)
         value = func(value, a->data[i]);
      return value;
   }
   
   template<typename t> struct maxfunc { t operator()(const t &a, const t &b) const { return max(a,b); } };
   template<typename t> struct minfunc { t operator()(const t &a, const t &b) const { return min(a,b); } };

}

#endif // _NEWARY_H_
