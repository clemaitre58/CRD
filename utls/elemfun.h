#ifndef __UTLS_ELEMFUN_H__
#define __UTLS_ELEMFUN_H__

#include <cmath>

#ifndef _MSC_VER
  #define CONST_TEMPLATE_PARAMETER const
#else
  #define CONST_TEMPLATE_PARAMETER
#endif

#ifndef M_PI                                                                    
  #define M_PI 3.1415926535897932384626433832795028841971693993751              
#endif                                                                          

namespace utls
   {

   template <typename ValueType>
   double l1norm2D(CONST_TEMPLATE_PARAMETER ValueType x, 
      CONST_TEMPLATE_PARAMETER ValueType y) { return fabs(x)+fabs(y); }

   template <typename ValueType>
   double l2norm2D(CONST_TEMPLATE_PARAMETER ValueType x, 
      CONST_TEMPLATE_PARAMETER ValueType y) { return ::sqrt((double)(x*x+y*y)); }


#ifdef _MSC_VER
#undef min
#undef max
#endif
   template <typename ValueType>
   ValueType min(CONST_TEMPLATE_PARAMETER ValueType x, 
      CONST_TEMPLATE_PARAMETER ValueType y) { if (x<y) return x; else return y; }

   template <typename ValueType>
   ValueType max(CONST_TEMPLATE_PARAMETER ValueType x, 
      CONST_TEMPLATE_PARAMETER ValueType y) { if (x>y) return x; else return y; }

   template <typename ValueType>
   ValueType min(
      CONST_TEMPLATE_PARAMETER ValueType x, 
      CONST_TEMPLATE_PARAMETER ValueType y,
      CONST_TEMPLATE_PARAMETER ValueType z) 
      { 
      if (x<y) { if (x<z) return x; } else { if (y<z) return y; }
      return z;
      }

   template <typename ValueType>
   ValueType max(
      CONST_TEMPLATE_PARAMETER ValueType x, 
      CONST_TEMPLATE_PARAMETER ValueType y,
      CONST_TEMPLATE_PARAMETER ValueType z) 
      { 
      if (x>y) { if (x>z) return x; } else { if (y>z) return y; }
      return z;
   }

   template <typename ValueType>
   ValueType min(
      ValueType x, ValueType y,
      ValueType z, ValueType w) 
   {
      if (x>y) x=y;
      if (x>z) x=z;
      if (x>w) x=w;
      return x;
   }

   template <typename ValueType>
   ValueType max(
      ValueType x, ValueType y,
      ValueType z, ValueType w) 
         {
      if (x<y) x=y;
      if (x<z) x=z;
      if (x<w) x=w;
            return x;
         }

   template <typename ValueType>
   ValueType min(
      ValueType x, ValueType y,
      ValueType z, ValueType w,
      ValueType v) 
   {
      if (x>y) x=y;
      if (x>z) x=z;
      if (x>w) x=w;
      if (x>v) x=v;
      return x;
   }

   template <typename ValueType>
   ValueType max(
      ValueType x, ValueType y,
      ValueType z, ValueType w,
      ValueType v)
   {
      if (x<y) x=y;
      if (x<z) x=z;
      if (x<w) x=w;
      if (x<v) x=v;
      return x;
      }

   template <typename ValueType>
   ValueType sqr(CONST_TEMPLATE_PARAMETER ValueType x) { return x*x; }

   template <typename ValueType>
   double dist(CONST_TEMPLATE_PARAMETER ValueType *x, CONST_TEMPLATE_PARAMETER ValueType *y, int n) 
      { 
      double res = 0;
      for (int i=0;i<n;i++)
         res += double(sqr(x[i]-y[i]));
      return sqrt(res);
      }

   template <typename XArray, typename YArray>
   inline void normalPDF (CONST_TEMPLATE_PARAMETER XArray &x, int len, YArray &y, double mi, double sigma)
      {
      for (int f = 0; f < len; f ++)
         {
         double xn = (x [f] - mi) / sigma;
         y [f] = exp (-0.5 * xn * xn)/(sqrt (2*M_PI) * sigma);
         }
      }

   template <typename XArray>                                                      
   inline XArray normalPDF (XArray &x, int len, double mi = 0, double sigma = 1)
      {                                                                               
      XArray v (len);
      normalPDF (x, len, v, mi, sigma);
      return v;
      }

   template <typename ValueType>
   inline ValueType triangleArea2D(CONST_TEMPLATE_PARAMETER ValueType px, CONST_TEMPLATE_PARAMETER ValueType py, 
      CONST_TEMPLATE_PARAMETER ValueType qx, CONST_TEMPLATE_PARAMETER ValueType qy, 
      CONST_TEMPLATE_PARAMETER ValueType rx, CONST_TEMPLATE_PARAMETER ValueType ry)
      {
      // negative => CW, positive => CCW
      return (px-qx)*(qy-ry) - (py-qy)*(qx-rx);
      }

   template <typename Point1, typename Point2, typename Point3>
   inline double triangleArea2D (Point1 &p, Point2 &q, Point3 &r)
      {
      // negative => CW, positive => CCW
      return (p.x - q.x)*(q.y - r.y) - (p.y - q.y)*(q.x - r.x);
      }


   /* Solve the square system of linear equations, Ax=b, where A is given
   in matrix "sq" and b in the vector "solution".  Result is given in
   solution.  Uses Gaussian elimination with pivoting. (from D.Lowe's util)
   */
#ifndef _MSC_VER
   template <typename ValueType, int rowlen>
   void SolveLinearSystem(ValueType *solution, ValueType sq[rowlen][rowlen], int size)
#else
   template <typename ValueType, int rowlen>
   void SolveLinearSystem(ValueType *solution, ValueType sq[3][3], int size)
#endif
      {
      int row, col, c, pivot = 0, i;
      ValueType maxc, coef, temp, mult, val;

      /* Triangularize the matrix. */
      for (col = 0; col < size - 1; col++) {
         /* Pivot row with largest coefficient to top. */
         maxc = -1.0;
         for (row = col; row < size; row++) {
            coef = sq[row][col];
            coef = (coef < 0.0 ? - coef : coef);
            if (coef > maxc) {
               maxc = coef;
               pivot = row;
               }
            }
         if (pivot != col) {
            /* Exchange "pivot" with "col" row (this is no less efficient
            than having to perform all array accesses indirectly). */
            for (i = 0; i < size; i++) {
               temp = sq[pivot][i];
               sq[pivot][i] = sq[col][i];
               sq[col][i] = temp;
               }
            temp = solution[pivot];
            solution[pivot] = solution[col];
            solution[col] = temp;
            }
         /* Do reduction for this column. */
         for (row = col + 1; row < size; row++) {
            mult = sq[row][col] / sq[col][col];
            for (c = col; c < size; c++)	/* Could start with c=col+1. */
               sq[row][c] -= mult * sq[col][c];
            solution[row] -= mult * solution[col];
            }
         }

      /* Do back substitution.  Pivoting does not affect solution order. */
      for (row = size - 1; row >= 0; row--) {
         val = solution[row];
         for (col = size - 1; col > row; col--)
            val -= solution[col] * sq[row][col];
         solution[row] = val / sq[row][row];
         }
      }

   };

#endif // __UTLS_ELEMFUN_H__
