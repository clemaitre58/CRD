#include <mex.h>
#include <mary.h>
#include <curvlinear.h>
#include <timeutls.h>

mxArray *ExportCurviLists(TableOfCurvilinearRegions &lists)
{
   mxArray *res = mxCreateCellMatrix(1, lists.size());
   for (size_t i=0;i<lists.size();i++)
   {
      mxArray *seg = mxCreateDoubleMatrix(7, lists[i]->size(), mxREAL);
      double *d = mxGetPr(seg);
      for (CurvilinearRegion::iterator it= lists[i]->begin(); it!=lists[i]->end(); it++)
      {
         *(d++) = it->PositionMaxi.x;
         *(d++) = it->PositionMaxi.y;
         *(d++) = it->PosInit.x;
         *(d++) = it->PosInit.y;
         *(d++) = it->PosFinal.x;
         *(d++) = it->PosFinal.y;
         *(d++) = it->Largeur;
      }
      mxSetCell(res, i, seg);
   }
   return res;
}

void ProcessParams(const mxArray *params, CCurvLinearParams &p)
{
   mxArray *tmp;
   // half of the width of the section
   if ((tmp=mxGetField(params, 0, "halfcut_width"))!=0) 
      p.LargeurEtude = (int)mxGetScalar(tmp);

   // width of the Fourier kernel
   if ((tmp=mxGetField(params, 0, "fourier_width"))!=0) 
      p.LargeurFourier = (int)mxGetScalar(tmp);
   
   // shortest edge list that is processed by the detector
   if ((tmp=mxGetField(params, 0, "shortest_edge"))!=0)
      p.LongueurMiniForme = (int)mxGetScalar(tmp);
   
   // minimum width of the curvilinear region 
   if ((tmp=mxGetField(params, 0, "minimum_width"))!=0)
      p.LargeurMini = (int)mxGetScalar(tmp);

   // threshold on response value
   if ((tmp=mxGetField(params, 0, "product_threshold"))!=0)
      p.SeuilReponse = mxGetScalar(tmp);

   // threshold on the variation of the width
   if ((tmp=mxGetField(params, 0, "max_width_variation"))!=0)
      p.VLocaleLargeur = mxGetScalar(tmp);

   // shortest region that is output
   if ((tmp=mxGetField(params, 0, "shortest_segment"))!=0)
      p.TailleRegionSortie = (int)mxGetScalar(tmp);

   // amount of smoothing (sigma) of the image before edge detection
   if ((tmp=mxGetField(params, 0, "smoothing_sigma"))!=0)
      p.EdgeImageSmoothing = mxGetScalar(tmp);

   if ((tmp=mxGetField(params, 0, "do_2d_constraints"))!=0)
      p.Do2dConstraints = mxGetScalar(tmp)>0;
   
   if ((tmp=mxGetField(params, 0, "do_filter_compact"))!=0)
      p.DoFilterCompact = mxGetScalar(tmp)>0;
   
   if ((tmp=mxGetField(params, 0, "do_merge_regions"))!=0)
      p.DoMergeRegions = mxGetScalar(tmp)>0;
   
   if ((tmp=mxGetField(params, 0, "do_remove_duplicates"))!=0)
      p.DoRemoveDuplicates = mxGetScalar(tmp)>0;

   if ((tmp=mxGetField(params, 0, "do_use_color"))!=0)
      p.DoUseColor = mxGetScalar(tmp)>0;
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
      //mexEvalString ("drawnow;");
      last_progress=p;
   }
}

void mexFunction(int nOp, mxArray *Op[], int nIp, const mxArray *Ip[])
{     
   if (nIp>0)
   {
      const mwSize *dims = mxGetDimensions(Ip[0]);
      int num_dims = mxGetNumberOfDimensions(Ip[0]);
      if ((num_dims==2 || (num_dims==3 && dims[2] == 3)) && mxGetClassID(Ip[0])==mxUINT8_CLASS)
      {
         CCurvLinearParams p;
         if (nIp>1)
         {
            if (mxIsStruct(Ip[1]))
               ProcessParams(Ip[1],p);
            else
               mexErrMsgTxt("curvidet: Second array should be structure of parameters.\n");
         }

         last_progress = 0;
         TableOfCurvilinearRegions res1, res2;
         utls::BAry *EdgeImage;

         double t = get_time();
         if (num_dims == 3)
         {
            // color version of detector...
            RGBImage *rgb = utls::array3d_from_matlab_interleaved<byte, RGBImage>(Ip[0]);
            if (p.DoUseColor)
            {
               mexPrintf("\nDetecting curvilinear segments on COLOR image %d x %d:\n  [", rgb->cols(), rgb->rows());
               DetectCurvLinearPoints(rgb, p, 0, &EdgeImage, res1, res2);
            } else {
				mexPrintf("\n il passe par là");
               Image *img = utls::convert_to_intensity<RGBImage, Image>(rgb);
               mexPrintf("\nDetecting curvilinear segments on GRAYSCALE image %d x %d:\n  [", img->cols(), img->rows());
               DetectCurvLinearPoints(img, p, 0, &EdgeImage, res1, res2);
               delete img;
            }
            delete rgb;
         } else {
            Image *img = utls::ary_from_matlab<byte, Image>(Ip[0]);
            mexPrintf("\nDetecting curvilinear segments on image %d x %d:\n  [", img->cols(), img->rows());
            DetectCurvLinearPoints(img, p, 0, &EdgeImage, res1, res2);
            delete img;
         }
         mexPrintf("%d curvilinear segments detected in %.2fsec.\n\n", res2.size(), get_time()-t);
         Op[0] = ExportCurviLists(res1);
         Op[1] = ExportCurviLists(res2);
         Op[2] = utls::ary_to_matlab<utls::BAry, byte>(EdgeImage);
         delete EdgeImage;
      } else 
         mexErrMsgTxt("curvidet: Grayscale image as uint8 2D matrix expected as paramater one.\n");
   }
}
