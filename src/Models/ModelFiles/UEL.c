#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"

/* Choose the numerical method */
#include "FEM.h"

#define TITLE   "Arbitrary UEL model"
#define AUTHORS "S.P.D."

#include "PredefinedMethods.h"

/*
 * The numbers below are arbitrary and serve only as example
 */

/* Nb of equations of the model */
//#define NEQ   (2)     /* Here let's consider an example with 2 equations */

/* Nb of terms per point */
#define NVI   (33)    /*  18 implicit terms per point: 9 stresses, 9 strains */
#define NVE   (0)     /*  2 explicit terms per point */
#define NV0   (0)     /*  2 constant terms per point */


/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)]) 


/* Intern Functions */
static int    pm(const char* s) ;

#if defined(__cplusplus)
  extern "C" {
#endif

extern void uel_(double*,double*,double*,double*,int*,int*,int*,double*,int*,double*,int*,int*,double*,double*,double*,double*,int*,double*,double*,int*,int*,int*,double*,int*,int*,double*,double*,int*,int*,int*,double*,int*,double*,int*,int*,double*) ;

#if defined(__cplusplus)
  }
#endif

/* Intern variables */
static double coef1, coef2,coef3 ;


int pm(const char* s)
{
  if(!strcmp(s,"neq"))  {
    return(99) ;
  } else if(!strcmp(s,"nvi"))  {
    return(98) ;
  } else if(!strncmp(s,"prop",4))  {
    const char* s4 = s + 4 ;
    int i ;
    
    sscanf(s4,"%d",&i) ;
    
    return(i) ;
  } else {
    return(-1) ;
  }
}


int SetModelProp(Model_t* model)
/** Set the model properties 
 *  Return 0 */
{
  /** Number of equations to be solved */
  //Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  {
    int neq = Model_GetNbOfEquations(model) ;
    int i ;
    
    for(i = 0 ; i < neq ; i++) {
      char name[100] ;
      
      sprintf(name,"Equation_%d",i+1);
      Model_CopyNameOfEquation(model,i,name) ;
    }
  }
  
  /** Names of the main (nodal) unknowns */
  {
    int neq = Model_GetNbOfEquations(model) ;
    int i ;
    
    for(i = 0 ; i < neq ; i++) {
      char name[100] ;
      
      sprintf(name,"Unknown_%d",i+1);
      Model_CopyNameOfUnknown(model,i,name) ;
    }
  }
  
  return(0) ;
}


int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 100 ;
  
  Material_ScanProperties(mat,datafile,pm) ;
  
  /* Number of equations */
  {
    int neq = floor(Material_GetProperty(mat)[pm("neq")] + 0.5) ;
    
    if(neq) {
      Material_GetNbOfEquations(mat) = neq ;
    } else {
      int dim = Material_GetDimension(mat) ;
      
      Material_GetNbOfEquations(mat) = dim ;
    }
    
    SetModelProp(Material_GetModel(mat)) ;
  }
  
  /* Number of implicit terms per integration points */
  {
    int nvi = floor(Material_GetProperty(mat)[pm("nvi")] + 0.5) ;
    
    if(!nvi) {
      Material_GetProperty(mat)[pm("nvi")] = NVI ;
    }
  }
  
  return(NbOfProp) ;
}


int PrintModelProp(Model_t* model,FILE *ficd)
/** Print the model properties 
 *  Return the nb of equations */
{
  int neq = Model_GetNbOfEquations(model) ;
  
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(0) ;
  
  printf("\n") ;
  printf("The set of equations is:\n") ;
  printf("\t- First equation      (name of Eq. 1)\n") ;
  printf("\t- Second equation     (name of Eq. 2)\n") ;
  printf("\t- ...                                \n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- First unknown       (name of Unk. 1)\n") ;
  printf("\t- Second unknown      (name of Unk. 2)\n") ;
  printf("\t- ...                                \n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;

  fprintf(ficd,"prop1 = 0.01   # Property 1\n") ;
  fprintf(ficd,"prop2 = 1.e-3  # Property 2\n") ;
  fprintf(ficd,"...                        \n") ;

  return(neq) ;
}


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;

  /** Define the length of tables */
  {
    int nvi = GetProperty("nvi") ;
    
    Element_GetNbOfImplicitTerms(el) = nvi*NbOfIntPoints ;
  }
  Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
  Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ;
  
  return(0) ;
}



int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/** Compute the residu (r) due to loads 
 *  Return 0 if succeeded and -1 if failed */
{
  IntFct_t* fi = Element_GetIntFct(el) ;
  int     nn   = Element_GetNbOfNodes(el) ;
  int     neq  = Element_GetNbOfEquations(el) ;
  int     ndof = nn*neq ;
  FEM_t*  fem  = FEM_GetInstance(el) ;
  int    i ;

  {
    double* r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
    for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
  }
  
  return(0) ;
}


int ComputeInitialState(Element_t* el)
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 *  Return 0 if succeeded and -1 if failed
 */ 
{
  double* vi = Element_GetImplicitTerm(el) ;
  double* ve = Element_GetExplicitTerm(el) ;
  double* v0 = Element_GetConstantTerm(el) ;
  
  /* We can skip if the element is a submanifold, 
   * e.g. a surface in 3D or a line in 2D */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Compute here vi, ve and v0 */
  {
    unsigned int i ;
    
    for(i = 0 ; i < Element_GetNbOfImplicitTerms(el) ; i++) {
      vi[i] = 0. ;
    }
    
  }
  
  return(0) ;
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the explicit terms.
 *  IMPORTANT: if needed use only with the previous values
 *  whatever they are, nodal values or implicit terms.
 *  Return 0 if succeeded and -1 if failed */
{
  double* ve = Element_GetExplicitTerm(el) ;
  /* If you need the implicit terms, use the previous ones */
  double* vi = Element_GetPreviousImplicitTerm(el) ;
  /* If you need the nodal values, use the previous ones */
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  /* If needed ! */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /* Compute here ve */
  {
  }
  
  return(0) ;
}


int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/** Compute the implicit terms 
 *  Return 0 if succeeded and -1 if failed */
{
  double* vi   = Element_GetImplicitTerm(el) ;
  double* vi_n = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  
  /* Compute here vi (with the help of vi_n if needed) */
  {
  }

  return(0) ;
}


int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) 
 *  Return 0 if succeeded and -1 if failed */
{
  double* vi   = Element_GetImplicitTerm(el) ;
  double* vi_n = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  int     nn   = Element_GetNbOfNodes(el) ;
  int     neq  = Element_GetNbOfEquations(el) ;
  int     ndof = nn*neq ;
  double  zero = 0. ;
  int     i ;

  /* Initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  

  /* Compute here the matrix K(i,j) */
  {
    Material_t* mat = Element_GetMaterial(el) ;
    int     dim  = Element_GetDimensionOfSpace(el) ;
    
    for(i = 0 ; i < Element_GetNbOfImplicitTerms(el) ; i++) {
      vi[i] = vi_n[i] ;
    }
    
    double* RHS = NULL ;
    double* AMATRX = k ;
    double* SVARS = vi ;
    double* ENERGY = NULL ;
    int     NDOFEL = ndof ;
    int     NRHS = ndof ;
    int     NSVARS = Element_GetNbOfImplicitTerms(el) ;
    double* PROPS  = Material_GetProperty(mat) ;
    int     NPROPS = Material_GetNbOfProperties(mat) ;
    double* COORDS = Element_ComputeNodalCoordinates(el) ;
    int     MCRD = dim ;
    int     NNODE = nn ;
    double* U = Element_ComputeCurrentNodalUnknowns(el) ;
    double* DU = Element_ComputeIncrementalNodalUnknowns(el) ;
    double* V = DU ;
    double* A = DU ;
    int     JTYPE = 1 ;
    double  TIME[2] = {0,t} ;
    double  DTIME = dt ;
    int     KSTEP = 1 ;
    int     KINC = 1 ;
    int     JELEM = 1 ;
    double  PARAMS[3] = {0,0,0} ;
    int     NDLOAD = 1 ;
    int     JDLTYP[2] = {0,1} ;
    double  ADLMAG[2] = {0,0} ;
    double  PREDEF[1] = {0} ;
    int     NPREDF = 1 ;
    int     LFLAGS[3] = {1,0,2} ;
    int     MLVARX = 1 ;
    double  DDLMAG[1] = {0};
    int     MDLOAD = 0 ;
    double  PNEWDT = 1 ;
    int     JPROPS[1] = {0} ;
    int     NJPROP = 1 ;
    double  PERIOD = 1 ;
    
    uel_(RHS,AMATRX,SVARS,ENERGY,&NDOFEL,&NRHS,&NSVARS,PROPS,&NPROPS,COORDS,&MCRD,&NNODE,U,DU,V,A,&JTYPE,TIME,&DTIME,&KSTEP,&KINC,&JELEM,PARAMS,&NDLOAD,JDLTYP,ADLMAG,PREDEF,&NPREDF,LFLAGS,&MLVARX,DDLMAG,&MDLOAD,&PNEWDT,JPROPS,&NJPROP,&PERIOD) ;
  }
  
  return(0) ;
}


int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/** Compute the residu (r) 
 *  Return 0 if succeeded and -1 if failed */
{
  double* vi   = Element_GetImplicitTerm(el) ;
  double* vi_n = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  int     nn   = Element_GetNbOfNodes(el) ;
  int     neq  = Element_GetNbOfEquations(el) ;
  int     ndof = nn*neq ;
  int     i ;


  /* Initialisation */
  for(i = 0 ; i < ndof ; i++) r[i] = 0. ;
  

  if(Element_IsSubmanifold(el)) return(0) ;

  /* Compute here the residu R(n,i) */
  {
    Material_t* mat = Element_GetMaterial(el) ;
    int     dim  = Element_GetDimensionOfSpace(el) ;
    
    for(i = 0 ; i < Element_GetNbOfImplicitTerms(el) ; i++) {
      vi[i] = vi_n[i] ;
    }
    
    double* RHS = r ;
    double* AMATRX = NULL ;
    double* SVARS = vi ;
    double* ENERGY = NULL ;
    int     NDOFEL = ndof ;
    int     NRHS = ndof ;
    int     NSVARS = Element_GetNbOfImplicitTerms(el) ;
    double* PROPS  = Material_GetProperty(mat) ;
    int     NPROPS = Material_GetNbOfProperties(mat) ;
    double* COORDS = Element_ComputeNodalCoordinates(el) ;
    int     MCRD = dim ;
    int     NNODE = nn ;
    double* U = Element_ComputeCurrentNodalUnknowns(el) ;
    double* DU = Element_ComputeIncrementalNodalUnknowns(el) ;
    double* V = DU ;
    double* A = DU ;
    int     JTYPE = 1 ;
    double  TIME[2] = {0,t} ;
    double  DTIME = dt ;
    int     KSTEP = 1 ;
    int     KINC = 1 ;
    int     JELEM = 1 ;
    double  PARAMS[3] = {0,0,0} ;
    int     NDLOAD = 1 ;
    int     JDLTYP[2] = {0,1} ;
    double  ADLMAG[2] = {0,0} ;
    double  PREDEF[1] = {0} ;
    int     NPREDF = 1 ;
    int     LFLAGS[3] = {1,0,5} ;
    int     MLVARX = 1 ;
    double  DDLMAG[1] = {0};
    int     MDLOAD = 0 ;
    double  PNEWDT = 1 ;
    int     JPROPS[1] = {0} ;
    int     NJPROP = 1 ;
    double  PERIOD = 1 ;
    
    uel_(RHS,AMATRX,SVARS,ENERGY,&NDOFEL,&NRHS,&NSVARS,PROPS,&NPROPS,COORDS,&MCRD,&NNODE,U,DU,V,A,&JTYPE,TIME,&DTIME,&KSTEP,&KINC,&JELEM,PARAMS,&NDLOAD,JDLTYP,ADLMAG,PREDEF,&NPREDF,LFLAGS,&MLVARX,DDLMAG,&MDLOAD,&PNEWDT,JPROPS,&NJPROP,&PERIOD) ;
  }
  
  return(0) ;
}


int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) 
 *  Return the nb of views (scalar, vector or tensor) */
{
  int neq          = Element_GetNbOfEquations(el) ;
  int nvi          = GetProperty("nvi") ;
  int NbOfOutputs  = neq + nvi ;


  if(Element_IsSubmanifold(el)) return(0) ;
  
  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    int    i ;
    
    i = 0  ;
    
    /* Primary unknowns */
    {
      double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
      FEM_t* fem = FEM_GetInstance(el) ;
      int j ;
      
      for(j = 0 ; j < neq ; j++) {
        double unk = FEM_ComputeUnknown(fem,u,intfct,p,j) ;
        char   str[5] ;
        
        sprintf(str,"U%d",j) ;
        
        Result_Store(r + i++,&unk,str,1) ;
      }
    }
    
    /* Secondary variables */
    {
      double* vi = Element_GetCurrentImplicitTerm(el) ;
      int np  = IntFct_GetNbOfPoints(intfct) ;
      int j ;
      
      for(j = 0 ; j < nvi ; j++) {
        double* vj = vi + j ;
        char   str[5] ;
        double sunk = 0 ;
        int k ;
        
        /* Averaging */
        for(k = 0 ; k < np ; k++ , vj += nvi) {
          sunk += vj[0] ;
        }
        
        sunk /= np ;
        
        sprintf(str,"V%d",j) ;
        
        Result_Store(r + i++,&sunk,str,1) ;
      }
    }
    
    if(i != NbOfOutputs) {
      Message_RuntimeError("ComputeOutputs: wrong number of outputs") ;
    }
  }

  return(NbOfOutputs) ;
}


    
/*
void UEL_(double* RHS,double* AMATRX,double* SVARS,double* ENERGY,int* NDOFEL,int* NRHS,int* NSVARS,double* PROPS,int* NPROPS,double* COORDS,int* MCRD,int* NNODE,double* U,double* DU,double* V,double* A,int* JTYPE,double* TIME,double* DTIME,int* KSTEP,int* KINC,int* JELEM,double* PARAMS,int* NDLOAD,int* JDLTYP,double* ADLMAG,double* PREDEF,int* NPREDF,int* LFLAGS,int* MLVARX,double* DDLMAG,int* MDLOAD,double* PNEWDT,int* JPROPS,int* NJPROP,double* PERIOD)
{
}
*/
