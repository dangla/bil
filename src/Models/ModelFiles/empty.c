#if 1
#include "PredefinedMethods.h"

int SetModelProp(Model_t *model)
{
}

#if 1
int ReadMatProp(Material_t *mat,DataFile_t *datafile)
{
  return(0) ;
}

int PrintModelChar(Model_t *model,FILE *ficd)
{
  return(0) ;
}

int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
{
  return(0) ;
}

int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
{
  return(0) ;
}

int ComputeInitialState(Element_t *el)
{
  return(0) ;
}

int  ComputeExplicitTerms(Element_t *el,double t)
{
  return(0) ;
}

int  ComputeImplicitTerms(Element_t *el,double t,double dt)
{
  return(0) ;
}

int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
{
  return(0) ;
}

int  ComputeResidu(Element_t *el,double t,double dt,double *r)
{
}

int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
{
  return(0) ;
}
#endif

#endif
