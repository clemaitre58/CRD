#include <vector>
#include "elemfun.h"

template <typename PointIterator, typename PointAdapter>
struct CSpaceHash
{
   typedef PointIterator                     HashPointIterator;
   typedef PointAdapter                      HashPointAdapter;
   typedef std::vector<PointIterator>        HashBin;
   typedef const std::vector<PointIterator>  ConstHashBin;
   typedef std::vector<HashBin *>            HashBinVector;
   typedef std::vector<ConstHashBin *>       ConstHashBinVector;
private:
   int m_width;
   int m_height;
   int m_granularity;
   int m_rows, m_cols;
public:
   HashBinVector m_hash;

   // allocate initial 2D hash
   CSpaceHash (int width, int height, 
               int granularity, double expected_fill): m_width(width), m_height(height), m_granularity(granularity) 
   {
      // compute number of bins
      m_cols = width/granularity +1;
      m_rows = height/granularity +1;
      // compute average fill
      expected_fill *= granularity*granularity;
      size_t bin_size = size_t(ceil(expected_fill));
      // allocate hash
      m_hash.resize(m_rows*m_cols);
      // preallocate bins
      for (typename HashBinVector::iterator it=m_hash.begin();it!=m_hash.end();it++)
      {
         (*it)=new HashBin;
         (*it)->reserve(bin_size);
      }
   };
   
   ~CSpaceHash()
   {
      for (typename HashBinVector::iterator it=m_hash.begin();it!=m_hash.end();it++)
         delete (*it);
   }
   // iterate through the all points and 
   void Hash(const PointIterator &first, const PointIterator &last)
   {
      int w = m_cols, h = m_rows;
      for (PointIterator it = first; it != last; it++)
      {
         int r = PointAdapter::y(it)/m_granularity,
             c = PointAdapter::x(it)/m_granularity;
         // check bin
         assert(r >= 0 && r < m_rows && c >= 0 && c < m_cols);
         m_hash[r*w+c]->push_back(it);
      }
   }

   // return content of bins at least partially covered by the rectangle
   void GetBins(int left, int top, int right, int bottom, ConstHashBinVector &bins)
   {
      int l = left / m_granularity, r = right / m_granularity+1, t = top / m_granularity, b = bottom / m_granularity+1;
      l = utls::max(0,l); t = utls::max(0,t); r = utls::min(r, int(m_cols)); b = utls::min(b, int(m_rows));
      bins.clear();    
      bins.reserve((r-l+1)*(b-t+1));
      t *= m_cols; b *= m_cols;
      for (int j=t; j<b; j+=m_cols)
         for (int i=l; i<r; i++)
            bins.push_back(m_hash[j+i]);
   }
};
