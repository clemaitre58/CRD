#include <mex.h>
#include <elemfun.h>
#include <vector>
#include <matpar.h>
#include <mary.h>
#include <timeutls.h>

#define MAX_DOUBLE 1e100

struct MatchingParameters 
{
   double  threshold;
   int     matching_strategy;
   int     min_consistent;
   MatchingParameters()
   {
      matching_strategy = 0;
      threshold = MAX_DOUBLE;
      min_consistent = 5;
   }
};

void ProcessParameters(const mxArray *mp, MatchingParameters &p)
{
   DPARAM(p, mp, threshold);
   IPARAM(p, mp, matching_strategy);
   IPARAM(p, mp, min_consistent);
}

static inline double distance(double *desc1, double *desc2, double upper_bound, int n)
{
   double dist = 0;
   double *end = desc1 + n;
   // compute distance and stop if it is over upper bound
   while (desc1 + 4 < end)
   {
      dist += (desc1[0] - desc2[0])*(desc1[0] - desc2[0]) +
              (desc1[1] - desc2[1])*(desc1[1] - desc2[1]) +
              (desc1[2] - desc2[2])*(desc1[2] - desc2[2]) +
              (desc1[3] - desc2[3])*(desc1[3] - desc2[3]) +
              (desc1[4] - desc2[4])*(desc1[4] - desc2[4]);
      if (dist > upper_bound)
         return MAX_DOUBLE;
      desc1 += 5;
      desc2 += 5;
   }

   while (desc1 + 1 < end)
   {
      dist += (desc1[0] - desc2[0])*(desc1[0] - desc2[0]) +
              (desc1[1] - desc2[1])*(desc1[1] - desc2[1]);
      if (dist > upper_bound)
         return MAX_DOUBLE;
      desc1 += 2;
      desc2 += 2;
   }
   if (desc1 < end)
      dist += (desc1[0] - desc2[0])*(desc1[0] - desc2[0]);
   return dist;
}

struct ShapeDesc
{
   int    x;
   int    y;
   double s;
   double angle;
   double S[3][3];
};

struct HashEntry
{
   double angle_diff;
   double scale_diff;
   double tx, ty;
};

struct MatchingPair : std::pair<const ShapeDesc*, const ShapeDesc*>
{
   double     T[3][3];
   HashEntry  h;
   int        matching;

   MatchingPair(const MatchingPair &other)
   {
      first = other.first;
      second = other.second;
      matching = 0;
      ComputeHashEntry();
   }

   MatchingPair(const ShapeDesc *f, const ShapeDesc *s)
   {
      first = f;
      second = s;
      matching = 0;
      ComputeHashEntry();
   }

   double NormaliseAngleDifference(double angle, double angle2)
   {
      angle = angle - angle2; // utls::min(angle - angle2 + 2*M_PI, angle2- angle + 2*M_PI);  
      while (angle > M_PI) angle -= M_PI;
      while (angle < 0) angle += M_PI;
      assert(angle>=0 && angle<=double(M_PI));
      return angle;
   }

   void ComputeHashEntry()
   {
      // difference in scaling
      h.scale_diff = first->s/second->s;

      // difference in angle
      h.angle_diff = first->angle-second->angle;

      double c1 = cos(h.angle_diff)/h.scale_diff, s1 = sin(h.angle_diff)/h.scale_diff;

      // h.tx = second->x-( first->x*c1 + first->y*s1);
      // h.ty = second->y-(-first->x*s1 + first->y*c1);

      h.angle_diff = utls::min(
         NormaliseAngleDifference(first->angle, second->angle),
         NormaliseAngleDifference(second->angle, first->angle));
   }
};

template <typename MatrixA, typename MatrixB, typename MatrixResult>
void matrixMultiply2 (CONST_TEMPLATE_PARAMETER MatrixA &a, CONST_TEMPLATE_PARAMETER MatrixB &b, MatrixResult &result)
{
   result[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0];
   result[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1];
   result[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0];
   result[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1];
}

template <typename MatrixA, typename MatrixB, typename MatrixResult>
void matrixMultiply3 (CONST_TEMPLATE_PARAMETER MatrixA &a, CONST_TEMPLATE_PARAMETER MatrixB &b, MatrixResult &result)
{
  result [0][0] = a [0][0] * b [0][0] + a [0][1] * b [1][0] + a [0][2] * b [2][0];
  result [0][1] = a [0][0] * b [0][1] + a [0][1] * b [1][1] + a [0][2] * b [2][1];
  result [0][2] = a [0][0] * b [0][2] + a [0][1] * b [1][2] + a [0][2] * b [2][2];
  result [1][0] = a [1][0] * b [0][0] + a [1][1] * b [1][0] + a [1][2] * b [2][0];
  result [1][1] = a [1][0] * b [0][1] + a [1][1] * b [1][1] + a [1][2] * b [2][1];
  result [1][2] = a [1][0] * b [0][2] + a [1][1] * b [1][2] + a [1][2] * b [2][2];
  result [2][0] = a [2][0] * b [0][0] + a [2][1] * b [1][0] + a [2][2] * b [2][0];
  result [2][1] = a [2][0] * b [0][1] + a [2][1] * b [1][1] + a [2][2] * b [2][1];
  result [2][2] = a [2][0] * b [0][2] + a [2][1] * b [1][2] + a [2][2] * b [2][2];
}

void ComputeTransformationMatrix(ShapeDesc &a)
{
   double c = cos(-a.angle)*a.s, s = sin(-a.angle)*a.s;
   a.S[0][0] = c; a.S[0][1] = s; a.S[0][2] = a.x;
   a.S[0][0] =-s; a.S[0][1] = c; a.S[0][2] = a.y;
   a.S[2][0] = 0; a.S[2][1] = 0; a.S[2][2] = 1;
}

//void ComputeTransformationMatrix(MatchingPair &m)
//{
//   double S[3][3], T[3][3];
//   // inverse transform (p1-> norm. coord)
//   double 
//      c1 = cos(m.first.angle)/m.first.s, 
//      s1 = sin(m.first.angle)/m.first.s;
//   
//   S[0][0] = c1; S[0][1] = s1; S[0][2] = -(m.first.x*c1+m.first.y*s1);
//   S[1][0] =-s1; S[1][1] = c1; S[1][2] = -(-m.first.x*s1+m.first.y*c1);
//   S[2][0] = 0; S[2][1] = 0; S[2][2] = 1;
//
//   // forward transform (norm.coord -> p2)
//   double 
//      c2 = cos(-m.second.angle)*m.second.s, 
//      s2 = sin(-m.second.angle)*m.second.s;
//   T[0][0] = c2; T[0][1] = s2; T[0][2] = m.second.x;
//   T[1][0] =-s2; T[1][1] = c2; T[1][2] = m.second.y;
//   T[2][0] = 0; T[2][1] = 0; T[2][2] = 1;
//
//   // product is transform 
//   matrixMultiply3(S, T, m.T);
//}

void ComputeTransformationMatrix(MatchingPair &m)
{
   // difference in scaling
   double s = m.first->s/m.second->s;
   // difference in angle
   double angle = m.first->angle-m.second->angle;
   double c1 = cos(angle)/s, s1 = sin(angle)/s;
   m.T[0][0] = c1; m.T[0][1] = s1; m.T[0][2] = m.second->x-( m.first->x*c1 + m.first->y*s1);
   m.T[1][0] =-s1; m.T[1][1] = c1; m.T[1][2] = m.second->y-(-m.first->x*s1 + m.first->y*c1);
   m.T[2][0] = 0;  m.T[2][1] = 0;  m.T[2][2] = 1;
}

double RatioOfReprojectionDistance(MatchingPair &m1, MatchingPair &m2)
{
   // compute transformation matrices...
   ComputeTransformationMatrix(m1);
   ComputeTransformationMatrix(m2);

   double scale = m1.first->s+m1.second->s;
   // relative distance of the points in second image
   double 
      rel_dist1 = sqrt(double(utls::sqr(m1.second->x - m2.second->x) + utls::sqr(m1.second->y - m2.second->y)));
   
   // relative distance of the projected points in second image with interchanged transformations
   double 
      x1 = m1.first->x * m2.T[0][0] + m1.first->y * m2.T[0][1] + m2.T[0][2],
      y1 = m1.first->x * m2.T[1][0] + m1.first->y * m2.T[1][1] + m2.T[1][2],
      
      x2 = m2.first->x * m1.T[0][0] + m2.first->y * m1.T[0][1] + m1.T[0][2],
      y2 = m2.first->x * m1.T[1][0] + m2.first->y * m1.T[1][1] + m1.T[1][2];
   double 
      rel_dist2 = sqrt(double(utls::sqr(x1 - x2) + utls::sqr(y1 - y2)));  

   // ideally this ratio will be close to 1
   return utls::max(rel_dist1/rel_dist2, rel_dist2/rel_dist1);
}

template <typename PairIterator>
void GeometryCheck(const PairIterator &begin, const PairIterator &end, MatchingParameters &mp)
{
   for (PairIterator it = begin; it != end; it++)
      it->matching = 0;

   for (PairIterator it = begin; it != end; it++)
   {      
      HashEntry &h1 = it->h;
      for (PairIterator jt = it+1; jt != end; jt++)
      {        
         HashEntry &h2 = jt->h;
         if ((h1.angle_diff - h2.angle_diff) < M_PI/10 && (h2.angle_diff - h1.angle_diff) < M_PI/10 && 
             (h1.scale_diff > h2.scale_diff*0.7)       && (h2.scale_diff > h1.scale_diff*0.7))
         {
            // ok... rotation and scaling is consistent, check reprojection errors with respect to scale of the pair...
            if (RatioOfReprojectionDistance(*it, *jt) < 1.05)
            {
               it->matching ++;
               jt->matching ++;
/*            mexPrintf("angle_diff: %g, scale_diff: %g, tx1: %g, ty1: %g, tx2: %g, ty2: %g\n", 
               fabs(h1.angle_diff - h2.angle_diff), 
               utls::min(h1.scale_diff/h2.scale_diff, h2.scale_diff/h1.scale_diff),
               h1.tx, h1.ty, h2.tx, h2.ty);*/
            }
         }               
      }
   }
}

void ImportDescs(const mxArray *meta_data, const mxArray *desc_data, std::vector<ShapeDesc> &meta, std::vector<double *> &descs)
{
   if (mxIsStruct(meta_data))
   {
      const mwSize *dims = mxGetDimensions(meta_data);
      size_t items = dims[0]*dims[1]; size_t desc_len = 0;
      { 
         const mxArray *tmp = mxGetField(meta_data, 0, "desc"); 
         if (tmp) desc_len = mxGetM(tmp)*mxGetN(tmp); else { mexErrMsgTxt("Missing DESCRIPTION field"); return; }
      }
      double *tmp = mxGetPr(desc_data);
      descs.resize(items);
      meta.resize(items);
      for (size_t i=0; i<items;i++)
      {         
         ShapeDesc &s = meta[i];
         safe_scalar_field_value(meta_data, s.x, "x", i);
         safe_scalar_field_value(meta_data, s.y, "y", i);
         safe_scalar_field_value(meta_data, s.s, "s", i);
         safe_scalar_field_value(meta_data, s.angle, "angle", i);
         descs[i] = tmp; tmp += desc_len;
      }
   }
}

void mexFunction(int nOp, mxArray *Op[], int nIp, const mxArray *Ip[])
{     
   if (nIp>3)
   {
      int d1 = mxGetN(Ip[1]), d2 = mxGetN(Ip[3]);

      int n = (int)mxGetM(Ip[1]);

      if (int(mxGetM(Ip[3])) != n)
      {
         mexErrMsgTxt("Descriptions should have same length = L ( LxM, LxN ) matrices are expected.\n");
         return;
      }
      MatchingParameters p;
      if (nIp>4 && mxIsStruct(Ip[4]))
         ProcessParameters(Ip[4], p);

      std::vector<double *> descs1, descs2;
      std::vector<ShapeDesc> meta1, meta2;

      ImportDescs(Ip[0], Ip[1], meta1, descs1);
      ImportDescs(Ip[2], Ip[3], meta2, descs2);

      Op[0] = mxCreateNumericMatrix(3,utls::min(d1,d2), mxINT32_CLASS, mxREAL);
      Op[1] = mxCreateNumericMatrix(1,d1,mxDOUBLE_CLASS, mxREAL);
      Op[2] = mxCreateNumericMatrix(1,d2,mxDOUBLE_CLASS, mxREAL);
      Op[3] = mxCreateNumericMatrix(1,d1,mxINT32_CLASS, mxREAL);
      Op[4] = mxCreateNumericMatrix(1,d2,mxINT32_CLASS, mxREAL);

      int num_tc=0;
   
      switch (p.matching_strategy)
      {
      case 0:
         {
            int 
               *idx_d1 = (int*)mxGetPr(Op[3]), 
               *idx_d2 = (int*)mxGetPr(Op[4]);
            double
               *minima_d1 = mxGetPr(Op[1]), 
               *minima_d2 = mxGetPr(Op[2]);

            for (int i=0;i<d1;i++)
            {
               idx_d1[i]=-1; minima_d1[i]=MAX_DOUBLE;
            }

            for (int i=0;i<d2;i++)
            {
               idx_d2[i]=-1; minima_d2[i]=MAX_DOUBLE;
            }

            double *dm1, *dm2; 
            int *id1, *id2, last1, last2, first1, first2;

            int chunk_size = 100;
            int cd1, cd2;
            
            cd1 = 1+(d1/chunk_size);
            cd2 = 1+(d2/chunk_size);
            double t1 = get_time();
            for (int chunk1=0; chunk1<cd1; chunk1++)
            {
               for (int chunk2=0; chunk2<cd2; chunk2++)
               {
                  first1 = last1 = chunk1*chunk_size;
                  dm1 = minima_d1 + last1;
                  id1 = idx_d1 + last1;
                  last1 = utls::min(last1 + chunk_size, d1);
                  for (int i = first1; i<last1; i++, dm1++,id1++)
                  {
                     first2 = last2 = chunk2*chunk_size;
                     dm2 = minima_d2 + last2; 
                     id2 = idx_d2 + last2;
                     last2 = utls::min(last2+chunk_size, d2);
                     for (int j = first2; j<last2; j++, dm2++,id2++)
            {
                        double upper_bound = utls::min(utls::max(*dm1, *dm2), p.threshold);
                        double dist = distance(descs1[i], descs2[j], upper_bound, n);
                  if (dist<*dm1)
                  {
                     *dm1=dist; *id1=j;
                  }
                  if (dist<*dm2)
                  {
                     *dm2=dist; *id2=i;
                  }
               }
            }
               }
            }
            std::vector<MatchingPair> pairs;
            int *tc = (int*)mxGetPr(Op[0]);

            id1=idx_d1;
            if (p.min_consistent>0)
            {
            for (int i=0; i<d1; i++,id1++)
                  if (idx_d2[*id1]==i)
                     pairs.push_back(MatchingPair(&meta1[i], &meta2[*id1]));
               mexPrintf("\n  Found %d tentative correspondences in %6.3f sec\n  Checking geometric consistency...\n", pairs.size(), get_time()-t1);
               GeometryCheck(pairs.begin(), pairs.end(), p);
               for (size_t i=0; i<pairs.size(); i++)
               {
                  if (pairs[i].matching >= p.min_consistent)
                  {
                     (*tc++) = pairs[i].first - &meta1[0] + 1;
                     (*tc++) = pairs[i].second - &meta2[0] + 1;
                     (*tc++) = pairs[i].matching;
                     num_tc++;
                  }
               }
            } else {
               for (size_t i=0; i<d1; i++, id1++)
            {
               if (idx_d2[*id1]==i)
               {
                  (*tc++)=i+1;
                  (*tc++)=*id1+1;
                     (*tc++) = 0;
                  num_tc++;
               }
            }
               mexPrintf("\n  Found %d tentative correspondences in %6.3f sec\n", pairs.size(), get_time()-t1);
            }
         }
         break;

      case 1:
         {
            int 
               *idx_d1 = (int*)mxGetPr(Op[3]), 
               *idx_d2 = (int*)mxGetPr(Op[4]);
            double
               *closest_d1 = mxGetPr(Op[1]), 
               *closest_d2 = mxGetPr(Op[2]);

            for (int i=0;i<d1;i++)
            {
               idx_d1[i]=-1; closest_d1[i]=MAX_DOUBLE;
            }

            double *dm1;
            int *id1;

            dm1 = closest_d1;
            id1 = idx_d1;

            int *tc = (int*)mxGetPr(Op[0]);
            int tmp = 0;
            for (int i=0; i<d1; i++,dm1++,id1++)
            {
               double second = MAX_DOUBLE;
               for (int j=0; j<d2; j++)
               {
                  double dist = distance(descs1[i], descs2[j], second, n);
                  tmp++;
                  if (dist < *dm1)
                  {
                     second = *dm1;
                     *dm1 = dist; 
                     *id1 = j;
                  } else
                     // if not closest, check it against second closest...
                     if (dist < second)
                        second = dist;
               }
               if ((*dm1)*1.25 <= second)
               {
                  (*tc++) = i+1;
                  (*tc++) = *id1+1;
                  (*tc++) = 0;
                  num_tc++;
               }
            }
         }
         break;
      case 2:
         {
         }
         break;
      }
      mxSetN(Op[0], num_tc);
   }
}

