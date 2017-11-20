#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

/* Choose the finite volume method */
#include "FVM.h"

#define TITLE "Short title of my model"
#define AUTHORS "Authors"

#include "PredefinedMethods.h"

/*
 * The numbers below are arbitrary and serve only as example
 */

/* Nb of equations of the model */
#define NEQ   (2)     /* Here let's consider an example with 2 equations */


/* Nb of nodes (el must be used below) */
#define NN     Element_GetNbOfNodes(el)


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
#define MassAndFlux(f,i,j)  ((f)[((i)*NN + (j))])

#define NW_1          (vi)
#define NW_1n         (vi_n)
#define N_1(i)        MassAndFlux(NW_1,i,i)
#define N_1n(i)       MassAndFlux(NW_1n,i,i)
#define W_1(i,j)      MassAndFlux(NW_1,i,j)

#define NW_2          (vi   + NN*NN)
#define NW_2n         (vi_n + NN*NN)
#define N_2(i)        MassAndFlux(NW_2,i,i)
#define N_2n(i)       MassAndFlux(NW_2n,i,i)
#define W_2(i,j)      MassAndFlux(NW_2,i,j)

#define PHI(i)        (vi   + 2*NN*NN)[i]
#define PHIn(i)       (vi_n + 2*NN*NN)[i]


/* We define some names for explicit terms (ve must be used as pointer below) */
#define TransferCoefficient(f,n)  ((f) + (n)*NN)

#define K_1        TransferCoefficient(ve,0)
#define K_2        TransferCoefficient(ve,1)


/* We define some names for constant terms (v0 must be used as pointer below) */
#define PR_1          (v0[0])
#define PR_2          (v0[1])


/* Material Properties 
 * ------------------- */
#define PCURVE1      (Element_GetCurve(el) + 0)
#define CURVE1(x)    (Curve_ComputeValue(PCURVE1,x))
#define PCURVE2      (Element_GetCurve(el) + 1)
#define CURVE2(x)    (Curve_ComputeValue(PCURVE2,x))


/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)]) 


/* Intern Functions */
static int    pm(const char* s) ;
static void   GetProperties(Element_t*) ;

static double* ComputeVariables(Element_t*,double**,double**,double*,double,double,int) ;
static void    ComputeSecondaryVariables(Element_t*,double,double,double*) ;


/* Intern parameters */
static double coef1 ;
static double coef2 ;
static double coef3 ;



/* Locally defined intern variables  */
#define NbOfVariables    (NEQ+4)
static double Variables[Element_MaxNbOfNodes][NbOfVariables] ;
static double dVariables[NbOfVariables] ;


/* We define some indices for the local variables */
#define I_Unk1         (IU_Unk1)
#define I_Unk2         (IU_Unk2)
#define I_N_1          (NEQ+0)
#define I_N_2          (NEQ+1)
#define I_PHI          (NEQ+2)
#define I_PHIn         (NEQ+3)
/* etc... */


/* Locally defined intern variable fluxes  */
#define NbOfVariableFluxes    (2)
static double VariableFluxes[Element_MaxNbOfNodes][NbOfVariableFluxes] ;

/* We define some indices for the local variable fluxes */
#define I_W_1           (0)
#define I_W_2           (1)


int pm(const char* s)
{
  if(strcmp(s,"prop1") == 0)        return (0) ;
  else if(strcmp(s,"prop2") == 0)   return (1) ;
  else if(strcmp(s,"prop3") == 0)   return (2) ;
  else return(-1) ;
}


void GetProperties(Element_t* el)
{
  coef1  = GetProperty("prop1") ;
  coef2  = GetProperty("prop2") ;
  coef3  = GetProperty("prop3") ;
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
  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
  
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


int ComputeInitialState(Element_t* el,double t)
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 *  Return 0 if succeeded and -1 if failed
 */ 
{
  double* vi  = Element_GetImplicitTerm(el) ;
  double* ve  = Element_GetExplicitTerm(el) ;
  double* v0  = Element_GetConstantTerm(el) ;
  
  /* Usually we have to skip if the element is a submanifold, 
   * e.g. a surface in 3D or a line in 2D */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    We load some input data
  */
  GetProperties(el) ;
  

  /* Pre-initialization */
  {
    int i ;
    
    for(i = 0 ; i < NN ; i++) {

    }
  }
  
  /* Compute here vi (with the help of vi_n if needed) */
  /* Loop on integration points */
  {
    int i ;
    
    for(i = 0 ; i < NN ; i++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,f_n,t,dt,i) ;

      /* Back up */
      N_1(i)  = x[I_N_1] ;
      N_2(i)  = x[I_N_2] ;
    }
  }
  
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Transfer coefficients */
  ComputeTransferCoefficients(el,u,f) ;

  /* Flux */
  {
    int i ;
    
    for(i = 0 ; i < NN ; i++) {
      int j ;
      
      for(j = i + 1 ; j < NN ; j++) {
        double* w = ComputeVariableFluxes(el,i,j) ;
        

        W_1(i,j)     = w[I_W_1 ] ;
        W_2(i,j)     = w[I_W_1 ] ;

        W_1(j,i)     = - w[I_W_2 ] ;
        W_2(j,i)     = - w[I_W_2 ] ;
      }
    }
  }

  
  return(0) ;
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the (current) explicit terms.
 *  IMPORTANT: if needed use only the previous values
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

  /*
    We load some input data
  */
  GetProperties(el) ;
  
  /*
    Transfer coefficient
  */
  
  ComputeTransferCoefficients(el,u,f) ;
  
  
  return(0) ;
}


int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/** Compute the (current) implicit terms 
 *  Return 0 if succeeded and -1 if failed */
{
  double* vi   = Element_GetCurrentImplicitTerm(el) ;
  double* vi_n = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;

  /*
    We load some input data
  */
  GetProperties(el) ;
  
  /* Compute here vi (with the help of vi_n if needed) */
  /* Loop on integration points */
  {
    int i ;
    
    for(i = 0 ; i < NN ; i++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,f_n,t,dt,i) ;

      /* Back up */
      N_1(i)  = x[I_N_1] ;
      N_2(i)  = x[I_N_2] ;
    }
  }
  
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Transfer coefficients */
  ComputeTransferCoefficients(el,u,f) ;

  /* Flux */
  {
    int i ;
    
    for(i = 0 ; i < NN ; i++) {
      int j ;
      
      for(j = i + 1 ; j < NN ; j++) {
        double* w = ComputeVariableFluxes(el,i,j) ;
        

        W_1(i,j)     = w[I_W_1 ] ;
        W_2(i,j)     = w[I_W_1 ] ;

        W_1(j,i)     = - w[I_W_2 ] ;
        W_2(j,i)     = - w[I_W_2 ] ;
      }
    }
  }

  return(0) ;
}


int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) 
 *  Return 0 if succeeded and -1 if failed */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])

  /*
    We load some input data
  */
  GetProperties(el) ;

  /* Compute here the matrix K(i,j) */
  {
    /* ...*/
  }
  
  return(0) ;
#undef K
}


int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double* f = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;

  
  /*
    Conservation of element 1: (N_1 - N_1n) + dt * div(W_1) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_1,NW_1n,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_1) -= r1[i] ;
    }
  }
  
  /*
    Conservation of element 2: (N_2 - N_2n) + dt * div(W_2) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_2,NW_2n,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_2) -= r1[i] ;
    }
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
    Result_Store(r + i++,vector ,"NameOfView_v",3) ; /* vector */
    Result_Store(r + i++,tensor ,"NameOfView_t",9) ; /* tensor */
    
    if(i != NbOfOutputs) {
      Message_RuntimeError("ComputeOutputs: wrong number of outputs") ;
    }
  }

  return(NbOfOutputs) ;
}



void ComputeTransferCoefficients(Element_t* el,double** u,double* f)
/* Transfer coefficients  */
{
  double* ve = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int i ;

  /* Initialization */
  for(i = 0 ; i < NVE ; i++) ve[i] = 0. ;


  for(i = 0 ; i < nn ; i++) {
    double* x = ComputeVariables(el,u,u,f,0,0,i) ;
    
    {
      double k1 ;
      double k2 ;
        
      K_1[i] = k1 ;
        
      K_2[i] = k2 ;
    }
  }
}



int TangentCoefficients(Element_t* el,double t,double dt,double* c)
/**  Tangent matrix coefficients (c) */
{
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  ObVal_t* obval = Element_GetObjectiveValue(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int    dec = NEQ*NEQ ;
  double dui[NEQ] ;
  FVM_t* fvm   = FVM_GetInstance(el) ;
  double* dist = FVM_ComputeIntercellDistances(fvm) ;
  int    i ;
  
  /* Initialization */
  for(i = 0 ; i < ndof*ndof ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  for(i = 0 ; i < NEQ ; i++) {
    dui[i] =  1.e-2*ObVal_GetValue(obval + i) ;
  }

  
  for(i = 0 ; i < nn ; i++) {
    /* Variables */
    double* xi   = ComputeVariables(el,u,u_n,f_n,t,dt,i) ;
    int k ;

    
    for(k = 0 ; k < NEQ ; k++) {
      double  dui_k  = dui[k] ;
      double* dxi    = ComputeVariableDerivatives(el,t,dt,xi,dui_k,k) ;
    
      /* Content terms at node i */
      {
        double* cii = c + (i*nn + i)*NEQ*NEQ ;
        
        cii[E_1*NEQ   + k] = dxi[I_N_1] ;
        cii[E_2*NEQ   + k] = dxi[I_N_2] ;
      }

      /* Transfer terms from node i to node j: d(wij)/d(ui_k) */
      {
        int j ;
        
        for(j = 0 ; j < nn ; j++) {
          if(j != i) {
            {
              double* cij = c + (i*nn + j)*NEQ*NEQ ;
              double dij  = dist[nn*i + j] ;
              double dtdij = dt/dij ;
              double* dw = ComputeFluxes(el,dxi,i,j) ;
        
              cij[E_1*NEQ    + k] = - dtdij * dw[I_W_1] ;
              cij[E_2*NEQ    + k] = - dtdij * dw[I_W_2] ;
            }
          }
        }
      }
    }
  }

  return(dec) ;
}



double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double t,double dt,int n)
{
  double* x = Variables[n] ;
  
  /* Primary Variables */
  x[IU_Unk1] = Unk1(n) ;
  x[IU_Unk2] = Unk2(n) ;
  
  /* Needed variables to compute secondary components */
  x[I_PHIn] = PHIn(n) ;

  {
    double* v0   = Element_GetConstantTerm(el) ;
    
    x[I_PR_1]  = PR_1 ;
    x[I_PR_2]  = PR_2 ;
  }
  
  ComputeSecondaryVariables(el,t,dt,x) ;
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x)
{
  /* Primary variables */
  double unk1    = x[IU_Unk1] ;
  double unk2    = x[IU_Unk2] ;
  
  /* Compute the secondary variables */
  double phi = ... ;
  double n_1 = ... ;
  double n_2 = ... ;
  
  /* Backup in x */
  x[I_N_1    ] = n_1 ;
  x[I_N_2    ] = n_2 ;
  
  x[I_PHI    ] = phi ;
}



double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double dui,int i)
{
  double* dx = dVariables ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += dui ;
  
  ComputeSecondaryVariables(el,t,dt,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dui ;
  }

  return(dx) ;
}



double* ComputeVariableFluxes(Element_t* el,double** u,int i,int j)
{
  double* grdij = dVariables ;

  /* Gradients */
  {
    int nn = Element_GetNbOfNodes(el) ;
    FVM_t* fvm   = FVM_GetInstance(el) ;
    double* dist = FVM_ComputeIntercellDistances(fvm) ;
    double dij  = dist[nn*i + j] ;
    
    {
      double* xi = Variables[i] ;
      double* xj = Variables[j] ;
      int k ;
    
      for(k = 0 ; k < NbOfVariables ; k++)  {
        grdij[k] = (xj[k] - xi[k]) / dij ;
      }
    }
  }

  /* Fluxes */
  {
    double* w = ComputeFluxes(el,grdij,i,j) ;
    
    return(w) ;
  }
}



double* ComputeFluxes(Element_t* el,double* grdij,int i,int j)
{
  double* va = Element_GetExplicitTerm(el) ;
  double* w  = VariableFluxes[i] ;
  
  {
    /* Gradients */
    double grd_u_1 = grdij[I_Unk1] ;
    double grd_u_2 = grdij[I_Unk2] ;
      
    /* Transfer terms */
    double k1   = 0.5 * (K_1[i] + K_1[j]) ;
    double k2   = 0.5 * (K_2[i] + K_2[j]) ;

    /* Fluxes */
    w[I_W_1]  += - k1 * grd_u_1  ;
    w[I_W_2]  += - k2 * grd_u_2  ;
  }
    
  return(w) ;
}
