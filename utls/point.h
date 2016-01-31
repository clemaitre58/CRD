#ifndef __UTLS_POINT_H__
#define __UTLS_POINT_H__

#include "elemfun.h"

#define FLOAT_SIMILARITY_TRESHOLD 0.000001
#define DOUBLE_SIMILARITY_TRESHOLD 0.000000001

namespace utls 
{

template <typename ValueType> 
struct Point2D
{
   ValueType x, y;
   Point2D () { x=0; y=0; }
   Point2D (const ValueType X, const ValueType Y) { x=X; y=Y; }
   Point2D (const Point2D &p) { x=p.x; y=p.y; }
   bool operator ==(const Point2D &p) const { return IsSimilar(p); }
   bool operator !=(const Point2D &p) const { return !IsSimilar(p); }
   
   Point2D operator + (const Point2D &p) const { return Point2D (x + p.x, y + p.y); }
   Point2D operator - (const Point2D &p) const { return Point2D (x - p.x, y - p.y); }
   ValueType operator * (const Point2D &p) const { return x * p.x + y * p.y; }

   Point2D operator + (ValueType d) const { return Point2D (x+d, y+d); }
   Point2D operator - (ValueType d) const { return Point2D (x-d, y-d); }
   Point2D operator * (ValueType d) const { return Point2D (x*d, y*d); }
   Point2D operator / (ValueType d) const { return Point2D (x/d, y/d); }
   
   bool operator < (const Point2D &p) const { return x<p.x; }

   double l2norm() const { return l2norm2D(x,y); }

   double distance(const Point2D &p) const { return l2norm2D(x-p.x, y-p.y); }
};

}
#endif // __UTLS_POINT_H__
