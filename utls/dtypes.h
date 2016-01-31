#ifndef _UTLS__DTYPES_H_
#define _UTLS__DTYPES_H_

namespace utls
{
   typedef unsigned char                                          t_byte;

   template <typename Value>
   union RGBValue 
   {
      typedef Value value;
      Value arr [3];
      struct {Value r,g,b;}; 
      Value &operator [](int index) { return arr[index]; }
      const Value &operator [](int index) const { return arr[index]; }

      RGBValue()
      {}

      // define set of basic operations
      RGBValue(Value v)
      {
         r = g = b = v;
      }

      RGBValue(Value ar, Value ag, Value ab)
      {
         r = ar;
         g = ag;
         b = ab;
      }

      RGBValue(const RGBValue &other)
      {
         r = other.r;
         g = other.g;
         b = other.b;
      }

      // scalar assignment
      RGBValue &operator = (const Value v)
      {
         r = g = b = v;
         return *this;
      }

      // RGB assignment
      RGBValue &operator = (const RGBValue &other)
      {
         r = other.r;
         g = other.g;
         b = other.b;
         return *this;
      }

      // basic arithmetic operations
      RGBValue operator + (const RGBValue &other) const
      {
         return RGBValue(r+other.r, g+other.g, b+other.b);
      }
      
      RGBValue &operator += (const RGBValue &other)
      {
         r += other.r; g += other.g; b += other.b; return *this;
      }

      RGBValue operator - (const RGBValue &other) const
      {
         return RGBValue(r-other.r, g-other.g, b-other.b);
      }

      RGBValue operator - () const 
      {
         return RGBValue(-r,-g,-b);
      }

      RGBValue &operator -= (const RGBValue &other)
      {
         r -= other.r; g -= other.g; b -= other.b; return *this;
      }

      RGBValue operator * (const RGBValue &other) const
      {
         return RGBValue(r*other.r, g*other.g, b*other.b);
      }

      RGBValue &operator *= (const RGBValue &other)
      {
         r *= other.r; g *= other.g; b *= other.b; return *this;
      }

      RGBValue operator /(const RGBValue &other) const
      {
         return RGBValue(r/other.r, g/other.g, b/other.b);
      }

      RGBValue &operator /= (const RGBValue &other)
      {
         r /= other.r; g /= other.g; b /= other.b; return *this;
      }

      // scalar arithmetic operations
      RGBValue operator + (const Value &value)
      {
         return RGBValue(r+value, g+value, b+value);
      }

      RGBValue operator - (const Value &value)
      {
         return RGBValue(r-value, g-value, b-value);
      }

      RGBValue operator * (const Value &value)
      {
         return RGBValue(r*value, g*value, b*value);
      }

      RGBValue operator /(const Value &value)
      {
         return RGBValue(r/value, g/value, b/value);         
      }

      template<typename to>
      operator RGBValue<to>() const 
      {
         return RGBValue<to>(to(r), to(g), to(b));
      }


   }; 

   template<typename Value>
   RGBValue<Value> sqrt(const RGBValue<Value> v)
   {
      return RGBValue<Value>(::sqrt(v.r), ::sqrt(v.g), ::sqrt(v.b));
   }

   template<typename Value>
   RGBValue<Value> exp(const RGBValue<Value> v)
   {
      return RGBValue<Value>(::exp(v.r), ::exp(v.g), ::exp(v.b));
   }

   typedef RGBValue<t_byte> t_rgb;
   typedef RGBValue<float>  t_frgb;
}

#endif // _UTLS__DTYPES_H_
