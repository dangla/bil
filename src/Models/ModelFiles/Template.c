#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

/* Choose the numerical method */
#include "FEM.h"

#define TITLE "Short title of my model"
#define AUTHORS "Authors"

#include "PredefinedMethods.h"

/*
 * The numbers below are arbitrary and serve only as example
 */

/* Nb of equations of the model */
#define NEQ   (2)     /* Here let's consider an example with 2 equations */

/* Nb of terms per point */
#define NVI   (9)     /*  9 implicit terms per point */
#define NVE   (2)     /*  2 explicit terms per point */
#define NV0   (2)     /*  2 constant terms per point */


/* Indices of equations */
#define IE_Eq1    (0)
#define IE_Eq2    (1)
/* Indices of unknowns */
#define IU_Unk1    (0)
#define IU_Unk2    (1)


/* Value of the nodal unknown (u and el must be used as pointers below) */
#define UNKNOWN(n,i)     Element_GetValueOfNodalUnknown(el,u,n,i)


/* We define some names for nodal unknowns */
#define Unk1(n)          (UNKNOWN(n,IU_Unk1))
#define Unk2(n)          (UNKNOWN(n,IU_Unk2))


/* We define some names for implicit terms (vi must be used as pointer below) */
#define N_1           (vi[0])
#define N_2           (vi[1])
#define W_1           (vi + 2) /* must be a 3D vector */
#define W_2           (vi + 5) /* must be a 3D vector */
#define PHI           (vi[8])


/* We define some names for explicit terms (ve must be used as pointer below) */
#define K_1           (ve[0])
#define K_2           (ve[1])


/* We define some names for constant terms (v0 must be used as pointer below) */
#define PR_1          (v0[0])
#define PR_2          (v0[1])


/* Material Properties 
 * ------------------- */
#define CURVE1(x)    Curve_ComputeValue(Element_GetCurve(el),x)
#define CURVE2(x)    Curve_ComputeValue(Element_GetCurve(el) + 1,x)


/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)]) 


/* Intern Functions */
static int    pm(const char* s) ;


/* Intern variables */
static double coef1, coef2,coef3 ;


int pm(const char* s)
{
  if(strcmp(s,"prop1") == 0)        return (0) ;
  else if(strcmp(s,"prop2") == 0)   return (1) ;
  else if(strcmp(s,"prop3") == 0)   return (2) ;
  else return(-1) ;
}


int SetModelProp(Model_t* model)
/** Set the model properties, return 0.
 *  Warning:
 *  Never call InternationalSystemOfUnits_UseAsLength() or similar
 *  to modify the units because this will also affect other models.
 */
{
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  Model_CopyNameOfEquation(model,IE_Eq1,"first") ;
  Model_CopyNameOfEquation(model,IE_Eq2,"second") ;
  
  /** Names of the main (nodal) unknowns */
  Model_CopyNameOfUnknown(model,IU_Unk1,"x") ;
  Model_CopyNameOfUnknown(model,IU_Unk2,"y") ;
  
  return(0) ;
}


int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 3 ;
  
  Material_ScanProperties(mat,datafile,pm) ;
  
  return(NbOfProp) ;
}


int PrintModelProp(Model_t* model,FILE* ficd)
/** Print the model properties 
 *  Return the nb of equations */
{
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(0) ;
  
  printf("\n") ;
  printf("The set of equations is:\n") ;
  printf("\t- First equation      (name of Eq. 1)\n") ;
  printf("\t- Second equation     (name of Eq. 2)\n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- First unknown       (name of Unk. 1)\n") ;
  printf("\t- Second unknown      (name of Unk. 2)\n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;

  fprintf(ficd,"prop1 = 0.01   # Property 1\n") ;
  fprintf(ficd,"prop2 = 1.e-3  # Property 2\n") ;
  fprintf(ficd,"prop3 = 1.e6   # Property 3\n") ;

  return(NEQ) ;
}


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;

  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
  Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ;
  
  /* Skip the rest of code for basic development.
   * For advanced developments find below 
   * some examples of possible operations */
  
  /* Compute some new interpolation functions */
  {
    int dim = Element_GetDimension(el) ;
    int nn  = Element_GetNbOfNodes(el) ;
    /* Replace "Type1" and "Type2" by existing types */
    int i   = IntFcts_AddIntFct(intfcts,nn,dim,"Type1") ;
    int j   = IntFcts_AddIntFct(intfcts,nn,dim,"Type2") ; /* not used */
    /* here j is equal to i + 1 ! */
    Element_GetIntFct(el) = IntFcts_GetIntFct(intfcts) + i ;
  }
  /* Remove dof on some nodes 
   * (e.g. at the middle of the edges of a triangle) */
  {
    int nn  = Element_GetNbOfNodes(el) ;
    int j = 3 ; /* dof = 3 to be suppressed */
    int i ;
    for(i = nn/2 ; i < nn ; i++) { /* edge nodes */
      Element_GetUnknownPosition(el)[i*NEQ + j]  = -1 ;
      Element_GetEquationPosition(el)[i*NEQ + j] = -1 ;
    }
  }
  return(0) ;
}



int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/** Compute the residu (r) due to loads 
 *  Return 0 if succeeded and -1 if failed */
{
  IntFct_t* fi = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FEM_t* fem = FEM_GetInstance(el) ;
  int    i ;

  {
    double* r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
    for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r1[i] ;
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
    ...
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
    ...
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
    ...
  }

  return(0) ;
}


int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) 
 *  Return 0 if succeeded and -1 if failed */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])

  /* Compute here the matrix K(i,j) */
  {
    ...
  }
  
  return(0) ;
#undef K
}


int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/** Compute the residu (r) 
 *  Return 0 if succeeded and -1 if failed */
{
#define R(n,i)    (r[(n)*NEQ+(i)])

  /* Compute here the residu R(n,i) */
  {
    ...
  }
  
  return(0) ;
#undef R
}


int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) 
 *  Return the nb of views (scalar, vector or tensor) */
{
  int NbOfOutputs  = 3 ;
  double scalar    = 1 ;
  double vector[3] = {1,2,3} ;
  double tensor[9] = {11,12,13,21,22,23,31,32,33} ;
  
  {
    i = 0  ;
    Result_Store(r + i++,&scalar,"NameOfView_x",1) ; /* scalar */
    Result_Store(r + i++,vector,"NameOfView_v",3)  ; /* vector */
    Result_Store(r + i++,tensor,"NameOfView_t",9)  ; /* tensor */
    
    if(i != NbOfOutputs) {
      Message_RuntimeError("ComputeOutputs: wrong number of outputs") ;
    }
  }

  return(NbOfOutputs) ;
}
