#include <mex.h>
#include <mary.h>
#include <matpar.h>
#include <curvlinear.h>
#include <timeutls.h>
#include <point.h>
#include <elemfun.h>
#include <raster.h>
using namespace utls;

struct OverlapParams
{
   bool verbose;
   OverlapParams()
   {
      verbose = true;
   }
};

void ProcessParams(const mxArray *params, OverlapParams &p)
{
   BPARAM(p, params, verbose);
}

int last_progress;

void Progress(int p)
{
   if (last_progress!=p)
   {
//      for (int i=last_progress/5; i<p/5; i++) mexPrintf("#");
      if (p==100)
         mexPrintf("]\n");
//      mexEvalString ("drawnow;");
      last_progress=p;
   }
}

bool ImportCurviLists(const mxArray *input, TableOfCurvilinearRegions &lists)
{
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
         it->ValDist        = atan2(double(it->PosFinal.y-it->PosInit.y), double(it->PosFinal.x-it->PosInit.x));
         it->Largeur        = l2norm2D(double(it->PosFinal.x-it->PosInit.x), double(it->PosFinal.y-it->PosInit.y));
      }
      lists.push_back(list);
   }
   return true;
}

struct RegionElement
{
   float             width;
   float             orientation;   
   const CurvilinearRegion *owner;
   RegionElement     *pair;
   RegionElement (float w, float ori, const CurvilinearRegion *o) : width(w), orientation(ori), owner(o), pair(0) {}
};

typedef std::vector<RegionElement> RegionElementVector;
typedef utls::Ary<RegionElementVector, RegionElementVector> RegionElementVectorImage;

void RenderRegion(RegionElementVectorImage *image, CurvilinearRegion *ListeCurvi, IPoint2D &bbmin, IPoint2D &bbmax)
{
   CurvilinearRegion::iterator pos = ListeCurvi->begin(), pos1;
   if (pos != ListeCurvi->end())
   {
      // at least one point?
      pos1 = pos; pos1++;
   } else
      return;
   bbmin.x = min(pos->PosInit.x, pos->PosFinal.x);
   bbmax.x = max(pos->PosInit.x, pos->PosFinal.x);
   bbmin.y = min(pos->PosInit.y, pos->PosFinal.y);
   bbmax.y = max(pos->PosInit.y, pos->PosFinal.y);
   while (pos1 != ListeCurvi->end())
	{
      bbmin.x = min(bbmin.x, pos1->PosInit.x, pos1->PosFinal.x);
      bbmax.x = max(bbmax.x, pos1->PosInit.x, pos1->PosFinal.x);
      bbmin.y = min(bbmin.y, pos1->PosInit.y, pos1->PosFinal.y);
      bbmax.y = max(bbmax.y, pos1->PosInit.y, pos1->PosFinal.y);

      RegionElement r(float(pos->Largeur+pos1->Largeur)/2.0f, float(pos1->ValDist), ListeCurvi);
      if (triangleArea2D(pos->PosInit, pos1->PosInit, pos1->PosFinal) < 0)
         accumulate_small_triangle(image, pos->PosInit, pos1->PosInit, pos1->PosFinal, r);
      else
         accumulate_small_triangle(image, pos->PosInit, pos1->PosFinal, pos1->PosInit, r);

      if (triangleArea2D(pos1->PosFinal, pos->PosFinal, pos->PosInit) < 0)
         accumulate_small_triangle(image, pos1->PosFinal, pos->PosFinal, pos->PosInit, r);
      else
         accumulate_small_triangle(image, pos1->PosFinal, pos->PosInit, pos->PosFinal, r);

      pos++; pos1++;
	}
   // Check image boundaries
   bbmin.x = min(max(image->lb2, bbmin.x), image->ub2+1);
   bbmax.x = max(min(bbmax.x, image->ub2+1), image->lb2);
   bbmin.y = min(max(image->lb1, bbmin.y), image->ub1+1);
   bbmax.y = max(min(bbmax.y, image->ub1+1), image->lb1);
}

void BuildImageFromList(IAry *plan1, CurvilinearRegion *ListeCurvi, IPoint2D &bbmin, IPoint2D &bbmax)
{
   CurvilinearRegion::iterator pos = ListeCurvi->begin(), pos1;
   if (pos != ListeCurvi->end())
   {
      // at least one point?
      pos1 = pos; pos1++;
   } else
      return;
   bbmin.x = min(pos->PosInit.x, pos->PosFinal.x);
   bbmax.x = max(pos->PosInit.x, pos->PosFinal.x);
   bbmin.y = min(pos->PosInit.y, pos->PosFinal.y);
   bbmax.y = max(pos->PosInit.y, pos->PosFinal.y);
   while (pos1 != ListeCurvi->end())
	{
      bbmin.x = min(bbmin.x, pos1->PosInit.x, pos1->PosFinal.x);
      bbmax.x = max(bbmax.x, pos1->PosInit.x, pos1->PosFinal.x);
      bbmin.y = min(bbmin.y, pos1->PosInit.y, pos1->PosFinal.y);
      bbmax.y = max(bbmax.y, pos1->PosInit.y, pos1->PosFinal.y);

      IAry::value v = 1;
      
      if (triangleArea2D(pos->PosInit, pos1->PosInit, pos1->PosFinal) < 0)
         add_triangle(plan1, pos->PosInit, pos1->PosInit, pos1->PosFinal, v);
      else
         add_triangle(plan1, pos->PosInit, pos1->PosFinal, pos1->PosInit, v);

      if (triangleArea2D(pos1->PosFinal, pos->PosFinal, pos->PosInit) < 0)
         add_triangle(plan1, pos1->PosFinal, pos->PosFinal, pos->PosInit, v);
      else
         add_triangle(plan1, pos1->PosFinal, pos->PosInit, pos->PosFinal, v);

      pos++; pos1++;
	}
   // Check image boundaries
   bbmin.x = min(max(plan1->lb2, bbmin.x), plan1->ub2+1);
   bbmax.x = max(min(bbmax.x, plan1->ub2+1), plan1->lb2);
   bbmin.y = min(max(plan1->lb1, bbmin.y), plan1->ub1+1);
   bbmax.y = max(min(bbmax.y, plan1->ub1+1), plan1->lb1);
}

template <typename Matrix>
double det3 (CONST_TEMPLATE_PARAMETER Matrix &m)
{
  return m [0][0] * m [1][1] * m [2][2] +
         m [0][1] * m [1][2] * m [2][0] +
         m [0][2] * m [1][0] * m [2][1] -
         m [0][2] * m [1][1] * m [2][0] -
         m [0][1] * m [1][0] * m [2][2] -
         m [0][0] * m [1][2] * m [2][1];
}

template <typename Matrix, typename MatrixInv>
void inv3 (CONST_TEMPLATE_PARAMETER Matrix &m, MatrixInv &result)
{
  double d = det3 (m);
  result [0][0] = (m [1][1] * m [2][2] - m [1][2] * m [2][1])/d;
  result [0][1] =-(m [0][1] * m [2][2] - m [0][2] * m [2][1])/d;
  result [0][2] = (m [0][1] * m [1][2] - m [0][2] * m [1][1])/d;
  result [1][0] =-(m [1][0] * m [2][2] - m [1][2] * m [2][0])/d;
  result [1][1] = (m [0][0] * m [2][2] - m [0][2] * m [2][0])/d;
  result [1][2] =-(m [0][0] * m [1][2] - m [0][2] * m [1][0])/d;
  result [2][0] = (m [1][0] * m [2][1] - m [1][1] * m [2][0])/d;
  result [2][1] =-(m [0][0] * m [2][1] - m [0][1] * m [2][0])/d;
  result [2][2] = (m [0][0] * m [1][1] - m [0][1] * m [1][0])/d;
}

void TransformPt(IPoint2D &pt, double H[3][3])
{
   double denom  = H[2][0]*pt.x + H[2][1]*pt.y + H[2][2], tmpx = pt.x;
   pt.x = int((H[0][0]*pt.x + H[0][1]*pt.y + H[0][2])/denom);
   pt.y = int((H[1][0]*tmpx + H[1][1]*pt.y + H[1][2])/denom);
}


template <typename Pointer>
void clearPixel(Pointer p)
{
   *p=10;
}

template <>
void clearPixel(RegionElementVectorImage::pointer p)
{
   p->clear();
}

template<typename AryBase>
void RemovePixelsOutsideTestFrame(AryBase *r, AryBase *t, double H[3][3])
{
   double Hinv[3][3]; inv3(H, Hinv);
   typename AryBase::pointer pr = r->first();
   typename AryBase::pointer pt = t->first();
   for (int y=r->lb1; y<=r->ub1; y++)
   {
      double denom  = Hinv[2][1]*y + Hinv[2][2], 
             tmpx   = Hinv[0][1]*y + Hinv[0][2],
             tmpy   = Hinv[1][1]*y + Hinv[1][2];
      for (int x=r->lb2; x<=r->ub2; x++)
      {
         // is point in the other image
         double d  = Hinv[2][0]*x + denom, tx, ty;
         tx = (Hinv[0][0]*x + tmpx)/d;
         ty = (Hinv[1][0]*x + tmpy)/d;
         if (ty < t->lb1 || ty >= t->ub1 || tx < t->lb2 || tx >= t->ub2)
         {
            // it is outside, clear accumulated points
            clearPixel(pr); clearPixel(pt);
         }
         pr++; pt++;
      }
   }
}

void TransformRegion(CurvilinearRegion *ListeCurvi, double H[3][3])
{
   for (CurvilinearRegion::iterator it = ListeCurvi->begin(); it != ListeCurvi->end(); it++)
	{
      TransformPt(it->PosInit, H);
      TransformPt(it->PosFinal, H);
      it->PositionMaxi.x = (it->PosInit.x+it->PosFinal.x)/2;
      it->PositionMaxi.y = (it->PosInit.y+it->PosFinal.y)/2;
      it->Largeur        = l2norm2D(double(it->PosFinal.x-it->PosInit.x), double(it->PosFinal.y-it->PosInit.y));
      it->ValDist        = atan2(double(it->PosFinal.y-it->PosInit.y), double(it->PosFinal.x-it->PosInit.x));
   }
}

double QualityFn(RegionElement a, RegionElement b)
{
   double angle_diff = utls::min(fmod(a.orientation - b.orientation + 2*M_PI, M_PI), fmod(b.orientation - a.orientation + 2*M_PI, M_PI));
   while (angle_diff > M_PI) angle_diff -= M_PI;
   while (angle_diff < 0) angle_diff += M_PI;
   double dist_diff  = (a.width - b.width)/min(a.width,b.width);
  
   return 
      exp(-angle_diff*angle_diff/(M_PI*M_PI/(10*10)))* // reward for small angle difference
      exp(-dist_diff*dist_diff/(0.2*0.2));  // reward for similar width
}

double Match(RegionElementVector *r, RegionElementVector *t)
{
   // deal fast with simple case
   if (r->size()==1 && t->size()==1)
   {
      // average width for match 
      return 4*QualityFn(*r->begin(), *t->begin()) /sqr((*r->begin()).width + (*t->begin()).width);
   }
   // find corresponding region pixels...
   RegionElementVector *shorter, *longer;
   if (t->size() > r->size())
   {
      shorter = t; 
      longer = r;
   } else {
      shorter = r; 
      longer = t;
   }  

   double Quality = 0;
   
   int ssize = shorter->size(),
       lsize = longer->size(),
       total = ssize * lsize, s, l;

   // compute distance table
   double *dist_table = new double[ssize*lsize];
   for (s = 0; s < ssize; s++)
      for (l = 0; l < lsize; l++)      
         dist_table[l*ssize+s] = QualityFn((*shorter)[s], (*longer)[l]);
      
   // iteratively in ssize*ssize*lsize find all pairings
   double max_quality;
   do {
      max_quality = 0;
      size_t i = 0, j;
      for (s = 0; s < total; s++)
      {
         if (max_quality<dist_table[s])
         {
            max_quality = dist_table[s];
            i = s;
         }
      }

      if (max_quality>0)
      {
         j = i/ssize; // get l index...
         i -= j*ssize; // get s index...

         for (s=0; s<ssize; s++) dist_table[j*ssize+s] = 0;
         for (l=0; l<lsize; l++) dist_table[l*ssize+i] = 0;

         (*shorter)[i].pair = &(*longer)[j];
         (*longer)[j].pair = &(*shorter)[i];

         // normalise quality of match by average width of the pair
         Quality += 4*max_quality/(utls::sqr(((*shorter)[i].width+(*longer)[j].width)));
      }
   } while (max_quality>0);
   // repeat until non-zero pairings are found...
   delete dist_table;
   return Quality;
}

double MatchCandidates(RegionElementVectorImage  *ref, RegionElementVectorImage *tst, DAry *ri, DAry *ti)
{
   int l;
   l = last_progress = 0;
   mexPrintf("Matching image...  [");

   RegionElementVectorImage::pointer ref_last = ref->last();
   RegionElementVectorImage::pointer ref_first = ref->first();
   DAry::pointer rptr = ri->first();
   size_t num_pixels = ref_last - ref_first;

   size_t ref_cnt = 0, tst_cnt = 0;
   size_t ref_max = 0, tst_max = 0;

   double cummulative_quality = 0, cummulative_npixels = 0;
   for (RegionElementVectorImage::pointer r = ref_first, t = tst->first(); r < ref_last; r++, t++, rptr++)
   {
      double Q, N = 0;
      
      // find matching quality
      if (r->size()>0)
      {
         // weighted sum of possible matches for this pixel
         for (size_t s = 0; s < r->size(); s++) N += 1.0/sqr(((*r)[s]).width);
      }
      if (t->size()>0)
      {
         // weighted sum of possible matches for this pixel
         for (size_t s = 0; s < t->size(); s++) N += 1.0/sqr(((*t)[s]).width);
      }
   
      if (r->size()>0 && t->size()>0) Q = Match(r, t); else Q = 0; 
      
      *rptr = Q;
      cummulative_quality += Q;
      cummulative_npixels += N;

      Progress(100*(r - ref_first)/num_pixels);
      ref_cnt += r->size(); if (ref_max < r->size()) ref_max = r->size();
      tst_cnt += t->size(); if (tst_max < t->size()) tst_max = t->size();
   }; Progress(100);

   mexPrintf("Total %d (max %d) entries in reference image, and %d (max %d) entries in test image, overlap score %5.2f%%\n", 
      ref_cnt, ref_max, tst_cnt, tst_max, 200*cummulative_quality/cummulative_npixels);

   return 2*cummulative_quality/cummulative_npixels;
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
         OverlapParams p;
         double  H[3][3];
         if (nIp>3 && mxGetClassID(Ip[3])==mxDOUBLE_CLASS)
         {
            const mwSize *dims = mxGetDimensions(Ip[3]);
            int num_dims = mxGetNumberOfDimensions(Ip[3]);
            mexPrintf("dims: %d, [%d %d]\n", num_dims, dims[0], dims[1]);
            if (num_dims == 2 && dims[0] == 3 && dims[1] == 3)
            {
               double *data = mxGetPr(Ip[3]);
               H[0][0] = data[0]; H[0][1] = data[3]; H[0][2] = data[6];
               H[1][0] = data[1]; H[1][1] = data[4]; H[1][2] = data[7];
               H[2][0] = data[2]; H[2][1] = data[5]; H[2][2] = data[8];
            } else 
               mexErrMsgTxt("cd_overlap: Fourth parameter should be double 3x3 homography matrix H: x_r = H * x_t, where x_r is point in reference image.\n");
         } else
            mexErrMsgTxt("cd_overlap: Fourth parameter should be double 3x3 homography matrix H: x_r = H * x_t, where x_r is point in reference image.\n");
         if (nIp>4)
         {
            if (mxIsStruct(Ip[4]))
               ProcessParams(Ip[4],p);
            else
               mexErrMsgTxt("cd_overlap: Third array should be structure of parameters.\n");
         }

         if (nIp>2 && !mxIsStruct(Ip[1]) && !mxIsStruct(Ip[2]))
         {
            TableOfCurvilinearRegions regs1, regs2;
            ImportCurviLists(Ip[1], regs1);
            ImportCurviLists(Ip[2], regs2);
            mexPrintf("\nReference image size %d x %d:\n  ", w, h);
            mexPrintf("Imported %d cregions from reference image and %d cregions from test image\n", regs1.size(), regs2.size());
            // render reference image
            RegionElementVectorImage* ref = new RegionElementVectorImage(0, h-1, 0, w-1);
            mexPrintf("Rendering reference image...  [");
            int l;
            l = last_progress = 0; 
            for (TableOfCurvilinearRegions::iterator it = regs1.begin(); it != regs1.end(); it++, l++)
            {
               IPoint2D bbmin, bbmax;
               RenderRegion(ref, *it, bbmin, bbmax);
               Progress(100*l/regs1.size());               
            }; Progress(100);
            // transform regions to reference frame
            mexPrintf("Transforming...  ["); 
            l = last_progress = 0; 
            for (TableOfCurvilinearRegions::iterator it = regs2.begin(); it != regs2.end(); it++, l++)
            {
               TransformRegion(*it, H);
               Progress(100*l/regs2.size());               
            }; Progress(100);
            // render test image
            RegionElementVectorImage* tst = new RegionElementVectorImage(h,w);
            mexPrintf("Rendering test image...  ["); last_progress = 0;
            l = last_progress = 0; 
            for (TableOfCurvilinearRegions::iterator it = regs2.begin(); it != regs2.end(); it++, l++)
            {
               IPoint2D bbmin, bbmax;
               RenderRegion(tst, *it, bbmin, bbmax);
               Progress(100*l/regs2.size());               
            }; Progress(100);
            RemovePixelsOutsideTestFrame(ref, tst, H);

            DAry *ref_tmp = new DAry(h,w),
                 *tst_tmp = new DAry(h,w);
            double overlap = MatchCandidates(ref, tst, ref_tmp, tst_tmp);
            Op[3] = mxCreateDoubleScalar(overlap);
            IAry 
               *tmp1 = new IAry(h,w), 
               *tmp2 = new IAry(h,w);
            tmp1->clear();
            for (TableOfCurvilinearRegions::iterator it = regs1.begin(); it != regs1.end(); it++, l++)
            {
               IPoint2D bbmin, bbmax;
               BuildImageFromList(tmp1, *it, bbmin, bbmax);
               delete (*it);
            };
            tmp2->clear();
            for (TableOfCurvilinearRegions::iterator it = regs2.begin(); it != regs2.end(); it++, l++)
            {
               IPoint2D bbmin, bbmax;
               BuildImageFromList(tmp2, *it, bbmin, bbmax);
               delete (*it);
            };
            RemovePixelsOutsideTestFrame(tmp1, tmp2, H);
            Op[0] = ary_to_matlab<IAry, int>(tmp1);
            Op[1] = ary_to_matlab<IAry, int>(tmp2);
            Op[2] = ary_to_matlab<DAry, double>(ref_tmp);
            delete tmp1; delete tmp2; delete ref_tmp; delete tst_tmp;
            delete ref;
            delete tst;
         }
      } else 
         mexErrMsgTxt("cd_overlap: Image expected as paramater one.\n");
   }
}
