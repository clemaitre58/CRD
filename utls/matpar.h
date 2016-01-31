#ifndef __UTLS_MATPAR_H__
#define __UTLS_MATPAR_H__

#define DPARAM(params, mp, name)                              \
   if (mxGetField(mp,0, #name)!=0) params. name = mxGetScalar(mxGetField(mp,0, #name));

#define FPARAM(params, mp, name)                              \
   if (mxGetField(mp,0, #name)!=0) params. name = (float)mxGetScalar(mxGetField(mp,0, #name));

#define IPARAM(params, mp, name)                                 \
   if (mxGetField(mp,0, #name)!=0)                          \
      (params. name) = (int)mxGetScalar(mxGetField(mp,0, #name));

#define BPARAM(params, mp, name)                              \
   if (mxGetField(mp,0, #name)!=0)                       \
      params. name = mxGetScalar(mxGetField(mp,0, #name))>0;

#define SPARAM(params, mp, name)                              \
   if (mxGetField(mp,0, #name)!=0)                       \
   {                                                     \
      mxArray *tmp = mxGetField(mp,0, #name);            \
      int size = mxGetN(tmp);                            \
      if (params. name) free(params. name);                      \
      params. name = (char *)malloc((size+1)*sizeof(char));  \
      mxGetString(tmp, params. name, size+1);                \
   }

#define _SPARAM(params, tmp, name)                             \
   {                                                      \
      int size = mxGetN(tmp);                             \
      if (params. name) free(params. name);                       \
      params. name = (char *)malloc((size+1)*sizeof(char));   \
      mxGetString(tmp, params. name, size+1);                 \
   }

template <typename result>
void scalar_field_value(const mxArray *matlab_structure, result &value, const char *field_name, int element_id)
{
   mxArray *tmp = mxGetField(matlab_structure, element_id, field_name);
   if (tmp!=0)                                                                       
      value = (result)mxGetScalar(tmp);
}

template <typename result>
void safe_scalar_field_value(const mxArray *matlab_structure, result &value, const char *field_name, int element_id)
{
   mxArray *tmp = mxGetField(matlab_structure, element_id, field_name);
   if (tmp!=0)                                                                       
      value = (result)mxGetScalar(tmp);
   else {
      char buf[1024];
#ifdef WIN32
#pragma warning ( suppress: 4996 )
      _snprintf(buf, 1023, "Missing field name: %s in element %d", field_name, element_id); buf[1023]=0;
#else
      snprintf(buf, 1023, "Missing field name: %s in element %d", field_name, element_id); buf[1023]=0;
#endif
      mexWarnMsgTxt(buf);
   }
}

template <>
void scalar_field_value(const mxArray *matlab_structure, char *&value, const char *field_name, int element_id)
{
   mxArray *tmp = mxGetField(matlab_structure, element_id, field_name);
   if (tmp!=0)                                                                       
   {
      int size = mxGetN(tmp);        
      if (value!=0) 
         delete(value);
      value = (char *)new char[size+1];
      mxGetString(tmp, value, size+1);
   }
}

template <>
void safe_scalar_field_value(const mxArray *matlab_structure, char *&value, const char *field_name, int element_id)
{
   mxArray *tmp = mxGetField(matlab_structure, element_id, field_name);
   if (tmp!=0)                                                                       
   {
      int size = mxGetN(tmp);        
      if (value!=0) 
         delete(value);
      value = (char *)new char[size+1];
      mxGetString(tmp, value, size+1);
   } else {
      char buf[1024];
#ifdef WIN32
#pragma warning ( suppress: 4996 )
      _snprintf(buf, 1023, "Missing field name: %s in element %d", field_name, element_id); buf[1023]=0;
#else
      snprintf(buf, 1023, "Missing field name: %s in element %d", field_name, element_id); buf[1023]=0;
#endif
      mexWarnMsgTxt(buf);
   }
}

#endif // __UTLS_MATPAR_H__
