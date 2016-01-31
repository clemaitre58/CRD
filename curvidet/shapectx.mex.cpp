#include <mex.h>

#include <cmath>
#include <algorithm>

#include <mary.h>
#include <matpar.h>
#include <curvlinear.h>
#include <timeutls.h>
#include <list>
#include <spacehash.h>
#include <point.h>
#include <elemfun.h>



const int NumberOfOrientations    = 8;
const int NumberOfDistances    = 3;

#define BINARIZE_ORIENTATION
#ifdef BINARIZE_ORIENTATION

const int NumberOfOrientationBins = 4;
const int DescLength              = NumberOfDistances*NumberOfOrientations*NumberOfOrientationBins;

#else
const int DescLength           = NumberOfDistances*NumberOfOrientations;

#endif

struct ShapeDesc
{
   int   x;
   int   y;
   double s;
   double angle;
   double desc[DescLength];
};

typedef std::vector<ShapeDesc> ShapeDescVector;


const char *keypoint_fields[]={ "x", "y", "s", "angle", "desc"};

mxArray *DescToMatlab(double *ivec)
{
   mxArray *res;
   mwSize dims[2] = { DescLength, 1 };   
   res = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
   double *data = (double *)mxGetData(res);
   for (int i=0; i<DescLength; i++)
      *data++ = ivec[i];
   return res;
}

mxArray *ExportDescs(ShapeDescVector &descs)
{
   mxArray *result=mxCreateStructMatrix(1, descs.size(), 5, keypoint_fields);
   int length = 0;
   for (ShapeDescVector::iterator it=descs.begin();it!=descs.end();it++)
   {
      mxSetField(result, length, "x", mxCreateScalarDouble(it->x));
      mxSetField(result, length, "y", mxCreateScalarDouble(it->y));
      mxSetField(result, length, "s", mxCreateScalarDouble(it->s));
      mxSetField(result, length, "angle", mxCreateScalarDouble(it->angle));
      mxSetField(result, length, "desc", DescToMatlab(it->desc));
      length++;
   }
   return result;
}

template <typename Point>
mxArray *ExportPoints(std::vector<Point> &pts)
{
   mxArray *result=mxCreateDoubleMatrix(2,pts.size(), mxREAL);
   double *d=mxGetPr(result);
   for (typename std::vector<Point>::iterator it=pts.begin();it!=pts.end();it++)
   {
      (*d++)=it->x;
      (*d++)=it->y;
   }
   return result;
}

bool ImportCurviLists(const mxArray *input, TableOfCurvilinearRegions &lists)
{
   const mwSize *dims = mxGetDimensions(input);
   int num_dims = mxGetNumberOfDimensions(input);
   if (num_dims>2 || mxGetClassID(input)!=mxCELL_CLASS)
   {
      mexErrMsgTxt("1D or 2D cell array expected!\n");
      return false;
   }
   mwSize cells = mxGetM(input)*mxGetN(input);
   for (mwSize i=0;i<cells;i++)
   {
      mxArray *cell = mxGetCell(input, i);
      const mwSize *dims = mxGetDimensions(cell);
      int num_dims = mxGetNumberOfDimensions(cell);
      if (num_dims!=2 || dims[0]!=7 || mxGetClassID(cell)!=mxDOUBLE_CLASS)
      {
         mexErrMsgTxt("Cells should contain 7xN DOUBLE 2D matrices!\n");
         return false;
      }
      double *d = mxGetPr(cell);
      CurvilinearRegion *list = new CurvilinearRegion;
      list->resize(dims[1]);
      for (CurvilinearRegion::iterator it= list->begin(); it!=list->end(); it++)
      {
         it->PositionMaxi.x = int(*(d++));
         it->PositionMaxi.y = int(*(d++));
         it->PosInit.x      = int(*(d++));
         it->PosInit.y      = int(*(d++));
         it->PosFinal.x     = int(*(d++));
         it->PosFinal.y     = int(*(d++));
         it->Largeur        = int(*(d++));
         // use positionmaxi for the medial axis point and valdist for orientation
         it->PositionMaxi.x = (it->PosInit.x+it->PosFinal.x)/2;
         it->PositionMaxi.y = (it->PosInit.y+it->PosFinal.y)/2;
         it->ValDist        = atan2(double(it->PosFinal.y-it->PosInit.y), double(it->PosFinal.x-it->PosInit.x))+M_PI/2.0f;
      }
      lists.push_back(list);
   }
   return true;
}

struct ShapeContextParams
{
   double context_size_factor; 
   double size_tolerance;
   int   sampling_rate;
   double computed_tolerance_low,
          computed_tolerance_high;
   ShapeContextParams()
   {
      context_size_factor = 6.0f;
      sampling_rate = 10;
      size_tolerance = 0.2;
      computed_tolerance_high = (1.0+size_tolerance);
      computed_tolerance_low = 1.0/computed_tolerance_high;
   }
};

void ProcessParams(const mxArray *params, ShapeContextParams &p)
{
   DPARAM(p, params, context_size_factor);
   IPARAM(p, params, sampling_rate);
   DPARAM(p, params, size_tolerance);
   p.computed_tolerance_high = (1.0+p.size_tolerance);
   p.computed_tolerance_low = 1.0/p.computed_tolerance_high;
}

int last_progress;

void Progress(int p)
{
   if (last_progress!=p)
   {
      for (int i=last_progress/5; i<p/5; i++)
         mexPrintf("#");
      if (p==100)
         mexPrintf("]\n");
      mexEvalString ("drawnow;");
      last_progress=p;
   }
}

struct CurvilinearPointAdapter
{
public:
   static int x(const CurvilinearRegion::iterator &it)
   {
      return it->PositionMaxi.x;
   }
   static int y(const CurvilinearRegion::iterator &it)
   {
      return it->PositionMaxi.y;
   }
   static double s(const CurvilinearRegion::iterator &it)
   {
      return it->Largeur;
   }
   static double angle(const CurvilinearRegion::iterator &it)
   {
      return it->ValDist;
   }
};

typedef CSpaceHash<CurvilinearRegion::iterator, CurvilinearPointAdapter> CurviHash;

class CurvilinearPointFilter
{
public:
   bool operator ()(ShapeDesc &desc, const CurvilinearRegion::iterator &it, ShapeContextParams &p)
   {
      double ratio = (desc.s/it->Largeur);
      return (ratio > p.computed_tolerance_low) && (ratio < p.computed_tolerance_high);
   }
};

//const double DistancesTable[NumberOfDistances+1]={0.2582,0.4472,0.6831,1.0000};
const double DistancesTable[NumberOfDistances+1]={0.15,0.40,0.65,1.00};
const double DistancesTableSquared[NumberOfDistances+1]={0.15*0.15,0.40*0.40,0.65*0.65,1.00};
               
inline static double NormaliseAngleDifference(double angle, double angle2)
{
   angle = angle - angle2; // utls::min(angle - angle2 + 2*M_PI, angle2- angle + 2*M_PI);  
   while (angle > 2.0*M_PI) angle -= 2.0*M_PI;
   while (angle < 0) angle += 2.0*M_PI;
   assert(angle>=0 && angle<=double(2.0*M_PI));
   return angle;
}

inline static double NormaliseAngleDifferenceModulo(double angle, double angle2)
{
   angle = angle - angle2;
   while (angle > M_PI) angle -= M_PI;
   while (angle < 0) angle += M_PI;
   assert(angle>=0 && angle<=double(M_PI));
   return angle;
}

template <typename PointAdapter, typename PointIterator, typename PointFilter>
void AccumulateShapeContext(PointFilter &filter, const PointIterator &begin, const PointIterator &end, ShapeDesc &desc, ShapeContextParams &p)
{
   double context_size = desc.s*p.context_size_factor; context_size *= context_size;
   // binarize all points 
   for (PointIterator it=begin; it!=end; it++)
   {
      // get points and filter them
      if (filter(desc, (*it), p))
      {
         // get coords
         double x = PointAdapter::x(*it) - desc.x;
         double y = PointAdapter::y(*it) - desc.y;
         // convert to logpolar
         double dist  = (x*x+y*y)/context_size;
         // throw away close and distant points
         // a = sqrt([0 1 3 7 15]/15) -> areas are in ratio 1:2:4:8
         // a =
         //     0    0.2582    0.4472    0.6831    1.0000
         if (dist >= DistancesTableSquared[0] && dist<=1.0)
         {
            // binarize distance
            int bin = 0;
            for (int i=1; i<=NumberOfDistances; i++, bin+=NumberOfOrientations)
               if (dist <= DistancesTableSquared[i])
                  break;

            // binarize difference in orientation
            double angle = atan2(y, x); 
            if (0)
            {
               angle = NormaliseAngleDifference(angle, desc.angle);
               angle *= NumberOfOrientations/double(2.0*M_PI);
            } else {
               angle = NormaliseAngleDifferenceModulo(angle, desc.angle);
               angle *= NumberOfOrientations/double(M_PI);
            }
            // the center of the first bin is at 0 -> subtract 0.5
            angle -= 0.5;
            int  iori = int(floor(angle));
            double ori = angle - iori;
            
#ifdef BINARIZE_ORIENTATION
            {
               // compute orientation bin
               double angle = PointAdapter::angle(*it);
               if (0)
               {
                  angle = NormaliseAngleDifference(angle, desc.angle);
                  angle *= NumberOfOrientationBins/double(2.0*M_PI);
               } else {
                  angle = NormaliseAngleDifferenceModulo(angle, desc.angle);
                  angle *= NumberOfOrientationBins/double(M_PI);
               }
               // the center of the first bin is at 0 -> subtract 0.5
               angle -= 0.5;
               int  iorib = int(floor(angle));
               double orib = angle - iorib;

               desc.desc[NumberOfOrientationBins * (bin+(iori+NumberOfOrientations)%NumberOfOrientations) + 
                  (iorib+NumberOfOrientationBins)%NumberOfOrientationBins ] += (1.0-ori)*(1.0-orib);

               desc.desc[NumberOfOrientationBins * (bin+(iori+1)%NumberOfOrientations) + 
                  (iorib+NumberOfOrientationBins)%NumberOfOrientationBins ] += ori*(1.0-orib);

               desc.desc[NumberOfOrientationBins * (bin+(iori+NumberOfOrientations)%NumberOfOrientations) + 
                  (iorib+1)%NumberOfOrientationBins ] += (1.0-ori)*orib;

               desc.desc[NumberOfOrientationBins * (bin+(iori+1)%NumberOfOrientations) + 
                  (iorib+1)%NumberOfOrientationBins ] += ori*orib;
            }           
#else
            desc.desc[bin+(iori+NumberOfOrientations)%NumberOfOrientations ] += 1.0-ori;
            desc.desc[bin+(iori+1)%NumberOfOrientations                    ] += ori;
#endif            
         }      
      }
   }
}

template <typename valuetype>
void NormVector(valuetype *vec, int len)
{
   double norm = 0.0;
   for (int i = 0; i < len; i++) norm += utls::sqr(vec[i]);
   if (norm<1e-5) return;
   norm = 1.0 / sqrt(norm);
   for (int i = 0; i < len; i++) vec[i] *= norm;
}

template <typename valuetype>
void RobustVector(valuetype *vec, int len)
{
   // robustify bins (0-2) ~ 0 (2-5) ~ 1 more ~ 2
   for (int i = 0; i < len; i++) 
      if (vec[i]<4) 
         vec[i]=0;
      else
         if (vec[i]<10)
            vec[i]=0.5;
         else 
            vec[i]=1;
}

template <typename HashType, typename PointIterator, typename PointFilter>
void ShapeContext(HashType &hash, 
                  const PointIterator &begin, 
                  const PointIterator &end, 
                  PointFilter &filter,
                  ShapeContextParams &p, ShapeDescVector &vect)
{
   typedef typename HashType::HashPointAdapter PointAdapter;
   
   typename HashType::ConstHashBinVector bins; bins.reserve(100);
   // sample provided points
   for (typename HashType::HashPointIterator it=begin; it!=end;)
   {
      // get points in the neighbourhood
      vect.push_back(ShapeDesc());
      ShapeDesc &d = vect.back();
      d.x = PointAdapter::x(it);
      d.y = PointAdapter::y(it);
      d.s = PointAdapter::s(it);
      d.angle = PointAdapter::angle(it);
      std::fill(d.desc, d.desc+DescLength, 0.0);
      int dist = int(ceil(d.s*p.context_size_factor*sqrt(2.0)));
      hash.GetBins(d.x-dist,d.y-dist,d.x+dist,d.y+dist, bins);
      // accumulate context
      for (typename HashType::ConstHashBinVector::iterator jt=bins.begin(); jt!=bins.end(); jt++)
         AccumulateShapeContext<PointAdapter>(filter, (*jt)->begin(), (*jt)->end(), d, p);

      // normalize rings...
#ifdef BINARIZE_ORIENTATION      
      NormVector(d.desc,      NumberOfOrientations*NumberOfOrientationBins); 
      NormVector(d.desc +     NumberOfOrientations*NumberOfOrientationBins, NumberOfOrientations*NumberOfOrientationBins); 
      NormVector(d.desc + 2 * NumberOfOrientations*NumberOfOrientationBins, NumberOfOrientations*NumberOfOrientationBins);    
#else
      NormVector(d.desc,      NumberOfOrientations); 
      NormVector(d.desc +     NumberOfOrientations, NumberOfOrientations); 
      NormVector(d.desc + 2 * NumberOfOrientations, NumberOfOrientations);    
#endif
      // scale it to range 0-255
      for (int i=0;i<DescLength;i++)
         d.desc[i]*=255.0;
      for (int i=0;i<p.sampling_rate && it!=end;i++,it++);
   }
}

struct EdgePoint
{
   float x, y, s, angle;
};

typedef std::vector<EdgePoint> EdgePointVector;

struct EdgePointAdapter
{
public:
   static int x(const EdgePointVector::iterator &it)
   {
      return int(it->x);
   }
   static int y(const EdgePointVector::iterator &it)
   {
      return int(it->y);
   }
   static double s(const EdgePointVector::iterator &it)
   {
      return it->s;
   }
   static double angle(const EdgePointVector::iterator &it)
   {
      return it->angle;
   }
};


typedef CSpaceHash<EdgePointVector::iterator, EdgePointAdapter> EdgePointHash;

class EdgePointFilter
{
public:
   bool operator ()(ShapeDesc &desc, const EdgePointVector::iterator &it, ShapeContextParams &p)
   {
      // filter for +-20% of width
      return 
         ((desc.s/it->s) < (1.0+p.size_tolerance)) && 
         ((it->s/desc.s) < (1.0+p.size_tolerance));
   }
};


void ImportEdgePoints(const mxArray *structure, EdgePointVector &v)
{
   int num_pts = mxGetM(structure)*mxGetN(structure);
   v.resize(num_pts);
   for (int id=0; id<num_pts; id++)
   {
      v[id].x = float(mxGetScalar(mxGetField(structure, id, "x")));
      v[id].y =float(mxGetScalar(mxGetField(structure, id, "y")));
      v[id].s =float(mxGetScalar(mxGetField(structure, id, "s")));
      v[id].angle =float(mxGetScalar(mxGetField(structure, id, "angle")));
   }
}


void mexFunction(int nOp, mxArray *Op[], int nIp, const mxArray *Ip[])
{     
   if (nIp>0)
   {
      int num_dims = mxGetNumberOfDimensions(Ip[0]);
      const mwSize *dims = mxGetDimensions(Ip[0]);
      int h, w; h = dims[0]; w = dims[1];
      if (num_dims==2 || (num_dims==3 && dims[2]==3))
      {
         ShapeContextParams p;
         if (nIp>2)
         {
            if (mxIsStruct(Ip[2]))
               ProcessParams(Ip[2],p);
            else
               mexErrMsgTxt("shapectx: Second array should be structure of parameters.\n");
         }

         if (nIp>1)
         {
            if (!mxIsStruct(Ip[1]))
            {
               TableOfCurvilinearRegions clists;
               ImportCurviLists(Ip[1], clists);
            
               double t1, t2, t3;
               mexPrintf("\nImported %d curvilinear segments from image %d x %d:\n  ", clists.size(), w, h);
               mexPrintf("Hashing...  ["); last_progress = 0;
               size_t l = 0;
               t1=get_time();
               CurviHash hash(w, h, 32, 0.2);            
               for (TableOfCurvilinearRegions::iterator it = clists.begin(); it != clists.end(); it++, l++)
               {
                  hash.Hash((*it)->begin(), (*it)->end());
                  Progress(100*l/clists.size());               
               }
               Progress(100);
               t2=get_time();
               mexPrintf("Hashing took %.3fsec\n",t2-t1);

               ShapeDescVector res;
               CurvilinearPointFilter filter;
               mexPrintf("  Contexts... ["); last_progress = 0;
               l = 0;
               for (TableOfCurvilinearRegions::iterator it = clists.begin(); it != clists.end(); it++, l++)
               {
                  ShapeContext(hash, (*it)->begin(), (*it)->end(), filter, p, res);
                  Progress(100*l/clists.size());   
               }
               Progress(100);
               t3=get_time();
               mexPrintf("%d contexts took %.3fsec, overall %.3fsec\n",res.size(), t3-t2, t3-t1);
               Op[0]=ExportDescs(res);
            } else {
               EdgePointVector edges;
               ImportEdgePoints(Ip[1], edges);
            
               double t1, t2, t3;
               mexPrintf("\nImported %d edge points from image %d x %d:\n  ", edges.size(), w, h);
               mexPrintf("Hashing...  ["); last_progress = 0;
               size_t l = 0;
               t1=get_time();
               EdgePointHash hash(w, h, 32, 0.2);
               hash.Hash(edges.begin(), edges.end());
               Progress(100);
               t2=get_time();
               mexPrintf("Hashing took %.3fsec\n",t2-t1);

               ShapeDescVector res;
               EdgePointFilter filter;
               mexPrintf("  Contexts... ["); last_progress = 0;
               l = 0;
               ShapeContext(hash,  edges.begin(), edges.end(), filter, p, res);
               Progress(100);
               t3=get_time();
               mexPrintf("%d contexts took %.3fsec, overall %.3fsec\n", res.size(), t3-t2, t3-t1);
               Op[0]=ExportDescs(res);
            }
         }
      } else 
         mexErrMsgTxt("shapectx: Grayscale image as uint8 2D matrix expected as paramater one.\n");
   }
}
