#ifndef OLDMETHODS_H
#define OLDMETHODS_H

#ifdef MODELINDEX
#define cat(x,y)   x##y
#define xcat(x,y)  cat(x,y)
#define QM         xcat(qm,MODELINDEX)
#define DM         xcat(dm,MODELINDEX)
#define TB         xcat(tb,MODELINDEX)
#define IN         xcat(in,MODELINDEX)
#define EX         xcat(ex,MODELINDEX)
#define MX         xcat(mx,MODELINDEX)
#define RS         xcat(rs,MODELINDEX)
#define CH         xcat(ch,MODELINDEX)
#define CT         xcat(ct,MODELINDEX)
#define SO         xcat(so,MODELINDEX)
#else
#define QM         qm
#define DM         dm
#define TB         tb
#define IN         in
#define EX         ex
#define MX         mx
#define RS         rs
#define CH         ch
#define CT         ct
#define SO         so
#endif


/* Methods */

#include "PredefinedMethods.h"

int SetModelProp(Model_t* model)
{
#ifdef NEQ
  Model_GetNbOfEquations(model) = NEQ ;
#endif
  return(0) ;
}

#include "Models.h"
#include "Mesh.h"
#include "Elements.h"
#include "Nodes.h"
#include "Loads.h"


/* Old Methods */
typedef int  dm_t(int,Material_t*,FILE*) ;
typedef int  qm_t(int,FILE*) ;
typedef void tb_t(Element_t*,IntFct_t*,unsigned int*,int) ;
typedef void in_t(double**,double**,double*,double*,Element_t,int,Symmetry_t) ;
typedef int  ex_t(double**,double**,double*,double*,Element_t,int,Symmetry_t,double) ;
typedef int  mx_t(double**,double**,double**,double*,double*,double*,double*,Element_t,int,Symmetry_t,double,double) ;
typedef void rs_t(double**,double**,double**,double*,double*,double*,double*,Element_t,int,Symmetry_t,double,double) ;
typedef void ch_t(double**,double**,double**,double*,double*,double*,double*,Element_t,int,Symmetry_t,double,double,Load_t) ;
typedef int  ct_t(double**,double**,double**,double*,double*,double*,Element_t,int,Symmetry_t,double,double) ;
typedef int  so_t(double**,double**,double*,double*,double*,Result_t*,Element_t,int,Symmetry_t,double) ;

static dm_t  DM ;
static qm_t  QM ;
static tb_t  TB ;
static in_t  IN ;
static ex_t  EX ;
static mx_t  MX ;
static rs_t  RS ;
static ch_t  CH ;
static ct_t  CT ;
static so_t  SO ;



int ReadMatProp(Material_t *mat,DataFile_t *datafile)
{
  FILE *ficd = DataFile_GetFileStream(datafile) ;
  int dim = Geometry_GetDimension(Model_GetGeometry(Material_GetModel(mat))) ;
  
  return(DM(dim,mat,ficd)) ;
}



int PrintModelChar(Model_t *model,FILE *ficd)
{
  return(QM(0,ficd)) ;
}



int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
{
  Geometry_t *geom = Element_GetGeometry(el) ;
  IntFct_t *fi = IntFcts_GetIntFct(intfcts) ;
  unsigned int n_fi = IntFcts_GetNbOfIntFcts(intfcts) ;
  unsigned int dim = Geometry_GetDimension(geom) ;
  
  TB(el,fi,&n_fi,dim) ;
  
  /* Begin: March 4, 2017 */
  Element_GetNbOfImplicitTerms(el) = el->n_vi ;
  Element_GetNbOfExplicitTerms(el) = el->n_ve ;
  /* End: March 4, 2017 */
  
  return(0) ;
}



int ComputeInitialState(Element_t *el)
{
  Geometry_t *geom = Element_GetGeometry(el) ;
  unsigned int dim = Geometry_GetDimension(geom) ;
  unsigned int neq = Element_GetNbOfEquations(el) ;
  Node_t **no = Element_GetPointerToNode(el) ;
  Symmetry_t sym = Geometry_GetSymmetry(geom) ;
  double *vi = Element_GetImplicitTerm(el) ;
  double *ve = Element_GetExplicitTerm(el) ;
  double *u[Element_MaxNbOfNodes] ;
  double uu[Element_MaxNbOfNodes][Model_MaxNbOfEquations] ;
  double *x[Element_MaxNbOfNodes] ;
  int i ;
  
  for(i = 0 ; i < Element_MaxNbOfNodes ; i++) u[i] = uu[i] ;
  
  for(i = 0 ; i < Element_GetNbOfNodes(el) ; i++) {
    int    j ;
    x[i] = Node_GetCoordinate(no[i]) ;
    for(j = 0 ; j < neq ; j++) {
      int jj = Element_GetNodalUnknownPosition(el,i,j) ;
      if(jj >= 0) {
        u[i][j] = Element_GetCurrentNodalUnknown(el,i)[jj] ;
      }
    }
  }
  
  IN(x,u,vi,ve,*el,dim,sym) ;
  
  for(i = 0 ; i < Element_GetNbOfNodes(el) ; i++) {
    int    j ;
    for(j = 0 ; j < neq ; j++) {
      int jj = Element_GetNodalUnknownPosition(el,i,j) ;
      if(jj >= 0) {
        Element_GetCurrentNodalUnknown(el,i)[jj] = u[i][j] ;
      }
    }
  }
  
  return(0) ;
}



int  ComputeExplicitTerms(Element_t *el,double t)
{
  Geometry_t *geom = Element_GetGeometry(el) ;
  unsigned int dim = Geometry_GetDimension(geom) ;
  unsigned int neq = Element_GetNbOfEquations(el) ;
  Node_t **no = Element_GetPointerToNode(el) ;
  Symmetry_t sym = Geometry_GetSymmetry(geom) ;
  double *vi = Element_GetPreviousImplicitTerm(el) ;
  double *ve = Element_GetExplicitTerm(el) ;
  double *u[Element_MaxNbOfNodes] ;
  double uu[Element_MaxNbOfNodes][Model_MaxNbOfEquations] ;
  double *x[Element_MaxNbOfNodes] ;
  int i ;
  
  for(i = 0 ; i < Element_MaxNbOfNodes ; i++) u[i] = uu[i] ;
  
  for(i = 0 ; i < Element_GetNbOfNodes(el) ; i++) {
    int    j ;
    x[i] = Node_GetCoordinate(no[i]) ;
    for(j = 0 ; j < neq ; j++) {
      int jj = Element_GetNodalUnknownPosition(el,i,j) ;
      if(jj >= 0) {
        u[i][j] = Element_GetPreviousNodalUnknown(el,i)[jj] ;
      }
    }
  }
  
  return(EX(x,u,vi,ve,*el,dim,sym,t)) ;
}



int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
{
  Geometry_t *geom = Element_GetGeometry(el) ;
  unsigned int dim = Geometry_GetDimension(geom) ;
  unsigned int neq = Element_GetNbOfEquations(el) ;
  Node_t **no = Element_GetPointerToNode(el) ;
  Symmetry_t sym = Geometry_GetSymmetry(geom) ;
  double *vi = Element_GetCurrentImplicitTerm(el) ;
  double *ve = Element_GetExplicitTerm(el) ;
  double *vi_n = Element_GetPreviousImplicitTerm(el) ;
  double *u_1[Element_MaxNbOfNodes],*u_n[Element_MaxNbOfNodes] ;
  double uu_1[Element_MaxNbOfNodes][Model_MaxNbOfEquations],uu_n[Element_MaxNbOfNodes][Model_MaxNbOfEquations] ;
  double *x[Element_MaxNbOfNodes] ;
  int i ;
  
  for(i = 0 ; i < Element_MaxNbOfNodes ; i++) { u_1[i] = uu_1[i] ; u_n[i] = uu_n[i] ; }
  
  for(i = 0 ; i < Element_GetNbOfNodes(el) ; i++) {
    int    j ;
    x[i] = Node_GetCoordinate(no[i]) ;
    for(j = 0 ; j < neq ; j++) {
      int jj = Element_GetNodalUnknownPosition(el,i,j) ;
      if(jj >= 0) {
        u_1[i][j] = Element_GetCurrentNodalUnknown(el,i)[jj] ;
        u_n[i][j] = Element_GetPreviousNodalUnknown(el,i)[jj] ;
      }
    }
  }
  
  return(MX(x,u_1,u_n,vi,vi_n,ve,k,*el,dim,sym,dt,t)) ;
}

int  ComputeResidu(Element_t *el,double t,double dt,double *r)
{
  Geometry_t *geom = Element_GetGeometry(el) ;
  unsigned int dim = Geometry_GetDimension(geom) ;
  unsigned int neq = Element_GetNbOfEquations(el) ;
  Node_t **no = Element_GetPointerToNode(el) ;
  Symmetry_t sym = Geometry_GetSymmetry(geom) ;
  double *vi = Element_GetCurrentImplicitTerm(el) ;
  double *ve = Element_GetExplicitTerm(el) ;
  double *vi_n = Element_GetPreviousImplicitTerm(el) ;
  double *u_1[Element_MaxNbOfNodes],*u_n[Element_MaxNbOfNodes] ;
  double uu_1[Element_MaxNbOfNodes][Model_MaxNbOfEquations],uu_n[Element_MaxNbOfNodes][Model_MaxNbOfEquations] ;
  double *x[Element_MaxNbOfNodes] ;
  int i ;
  
  for(i = 0 ; i < Element_MaxNbOfNodes ; i++) { u_1[i] = uu_1[i] ; u_n[i] = uu_n[i] ; }
  
  for(i = 0 ; i < Element_GetNbOfNodes(el) ; i++) {
    int    j ;
    x[i] = Node_GetCoordinate(no[i]) ;
    for(j = 0 ; j < neq ; j++) {
      int jj = Element_GetNodalUnknownPosition(el,i,j) ;
      if(jj >= 0) {
        u_1[i][j] = Element_GetCurrentNodalUnknown(el,i)[jj] ;
        u_n[i][j] = Element_GetPreviousNodalUnknown(el,i)[jj] ;
      }
    }
  }
  
  RS(x,u_1,u_n,vi,vi_n,ve,r,*el,dim,sym,dt,t) ;
  
  return(0) ;
}



int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
{
  Geometry_t *geom = Element_GetGeometry(el) ;
  unsigned int dim = Geometry_GetDimension(geom) ;
  unsigned int neq = Element_GetNbOfEquations(el) ;
  Node_t **no = Element_GetPointerToNode(el) ;
  Symmetry_t sym = Geometry_GetSymmetry(geom) ;
  double *vi = Element_GetCurrentImplicitTerm(el) ;
  double *ve = Element_GetExplicitTerm(el) ;
  double *vi_n = Element_GetPreviousImplicitTerm(el) ;
  double *u_1[Element_MaxNbOfNodes],*u_n[Element_MaxNbOfNodes] ;
  double uu_1[Element_MaxNbOfNodes][Model_MaxNbOfEquations],uu_n[Element_MaxNbOfNodes][Model_MaxNbOfEquations] ;
  double *x[Element_MaxNbOfNodes] ;
  int i ;
  
  for(i = 0 ; i < Element_MaxNbOfNodes ; i++) { u_1[i] = uu_1[i] ; u_n[i] = uu_n[i] ; }
  
  for(i = 0 ; i < Element_GetNbOfNodes(el) ; i++) {
    int    j ;
    x[i] = Node_GetCoordinate(no[i]) ;
    for(j = 0 ; j < neq ; j++) {
      int jj = Element_GetNodalUnknownPosition(el,i,j) ;
      if(jj >= 0) {
        u_1[i][j] = Element_GetCurrentNodalUnknown(el,i)[jj] ;
        u_n[i][j] = Element_GetPreviousNodalUnknown(el,i)[jj] ;
      }
    }
  }
  
  CH(x,u_1,u_n,vi,vi_n,ve,r,*el,dim,sym,dt,t,*cg) ;
  
  return(0) ;
}



int  ComputeImplicitTerms(Element_t *el,double t,double dt)
{
  Geometry_t *geom = Element_GetGeometry(el) ;
  unsigned int dim = Geometry_GetDimension(geom) ;
  unsigned int neq = Element_GetNbOfEquations(el) ;
  Node_t **no = Element_GetPointerToNode(el) ;
  Symmetry_t sym = Geometry_GetSymmetry(geom) ;
  double *vi = Element_GetCurrentImplicitTerm(el) ;
  double *ve = Element_GetExplicitTerm(el) ;
  double *vi_n = Element_GetPreviousImplicitTerm(el) ;
  double *u_1[Element_MaxNbOfNodes],*u_n[Element_MaxNbOfNodes] ;
  double uu_1[Element_MaxNbOfNodes][Model_MaxNbOfEquations],uu_n[Element_MaxNbOfNodes][Model_MaxNbOfEquations] ;
  double *x[Element_MaxNbOfNodes] ;
  int i ;
  
  for(i = 0 ; i < Element_MaxNbOfNodes ; i++) { u_1[i] = uu_1[i] ; u_n[i] = uu_n[i] ; }
  
  for(i = 0 ; i < Element_GetNbOfNodes(el) ; i++) {
    int    j ;
    x[i] = Node_GetCoordinate(no[i]) ;
    for(j = 0 ; j < neq ; j++) {
      int jj = Element_GetNodalUnknownPosition(el,i,j) ;
      if(jj >= 0) {
        u_1[i][j] = Element_GetCurrentNodalUnknown(el,i)[jj] ;
        u_n[i][j] = Element_GetPreviousNodalUnknown(el,i)[jj] ;
      }
    }
  }
  
  return(CT(x,u_1,u_n,vi,vi_n,ve,*el,dim,sym,dt,t)) ;
}



int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
{
  Geometry_t *geom = Element_GetGeometry(el) ;
  unsigned int dim = Geometry_GetDimension(geom) ;
  unsigned int neq = Element_GetNbOfEquations(el) ;
  Node_t **no = Element_GetPointerToNode(el) ;
  Symmetry_t sym = Geometry_GetSymmetry(geom) ;
  double *vi = Element_GetImplicitTerm(el) ;
  double *ve = Element_GetExplicitTerm(el) ;
  double *u[Element_MaxNbOfNodes] ;
  double uu[Element_MaxNbOfNodes][Model_MaxNbOfEquations] ;
  double *x[Element_MaxNbOfNodes] ;
  int i ;
  
  for(i = 0 ; i < Element_MaxNbOfNodes ; i++) u[i] = uu[i] ;
  
  for(i = 0 ; i < Element_GetNbOfNodes(el) ; i++) {
    int    j ;
    x[i] = Node_GetCoordinate(no[i]) ;
    for(j = 0 ; j < neq ; j++) {
      int jj = Element_GetNodalUnknownPosition(el,i,j) ;
      if(jj >= 0) {
        u[i][j] = Element_GetCurrentNodalUnknown(el,i)[jj] ;
      }
    }
  }
  
  return(SO(x,u,vi,ve,s,r,*el,dim,sym,t)) ;
}

#endif
