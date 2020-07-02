/* General features of the model:
 * Simple diffusion of solute (as sodium or potassium)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

/* The Finite Volume Method */
#include "FVM.h"

#define TITLE   "Bo model (2020)"
#define AUTHORS "Ran-Dangla"

#include "PredefinedMethods.h"




/* Indices of equations/unknowns */
enum {
  E_sulfur   ,
  E_Last
} ;


/* Nb of equations */
#define NbOfEquations      E_Last
#define NEQ                NbOfEquations




/* Value of the nodal unknown (u, u_n and el must be used below) */
#define UNKNOWN(n,i)     Element_GetValueOfNodalUnknown(el,u,n,i)
#define UNKNOWNn(n,i)    Element_GetValueOfNodalUnknown(el,u_n,n,i)


/* Generic names of nodal unknowns */
#define U_Sulfur(n)      (UNKNOWN(n,E_sulfur))
#define Un_Sulfur(n)     (UNKNOWNn(n,E_sulfur))




/* Method chosen at compiling time. 
 * Each equation is associated to a specific unknown.
 * Each unknown can deal with a specific model.
 * Uncomment/comment to let only one unknown per equation */

/* solute: unknown either C or logC */
#define U_C_SO3     U_Sulfur
//#define U_LogC_SO3  U_Sulfur





/* Names of nodal unknowns */
#if defined (U_LogC_SO3)
  #define LogC_SO3(n)    U_Sulfur(n)
  #define LogC_SO3n(n)   Un_Sulfur(n)
  #define C_SO3(n)       (pow(10,LogC_SO3(n)))
  #define C_SO3n(n)      (pow(10,LogC_SO3n(n)))
#elif defined (U_C_SO3)
  #define C_SO3(n)       U_Sulfur(n)
  #define C_SO3n(n)      Un_Sulfur(n)
  #define LogC_SO3(n)    (log10(C_SO3(n)))
  #define LogC_SO3n(n)   (log10(C_SO3n(n)))
#else
  #error "Ambiguous or undefined unknown"
#endif




/* Nb of nodes (el must be used below) */
#define NN     Element_GetNbOfNodes(el)


/* Nb of terms per point */
#define NVE    	  (NN)
#define NVI       (NN*NN + 2*NN)
#define NV0       (0)


/* Names used for implicit terms */
#define MassAndFlux(f,i,j)  ((f)[((i)*NN + (j))])

#define NW_S         (f)
#define NW_Sn        (f_n)
#define N_S(i)       MassAndFlux(NW_S,i,i)
#define N_Sn(i)      MassAndFlux(NW_Sn,i,i)
#define W_S(i,j)     MassAndFlux(NW_S,i,j)

#define N_S_Cry(i)   (f   + NN*NN)[i]
#define N_S_Cryn(i)  (f_n + NN*NN)[i]

#define Phi(i)       (f   + NN*NN + NN)[i]
#define Phin(i)      (f_n + NN*NN + NN)[i]


/* Names used for explicit terms */
#define TransferCoefficient(va,i)  ((va) + (i)*NN)

#define KF_SO3          TransferCoefficient(va,0)


/* Names used for constant terms */
/* nothing */






/* Math constants */
#define Ln10      Math_Ln10




/* Units
 * ----- */
#include "InternationalSystemOfUnits.h"
/* Shorthands of some units */
#define dm    (0.1*InternationalSystemOfUnits_OneMeter)
#define cm    (0.01*InternationalSystemOfUnits_OneMeter)
#define dm2   (dm*dm)
#define dm3   (dm*dm*dm)
#define cm3   (cm*cm*cm)
#define MPa   (1.e6*InternationalSystemOfUnits_OnePascal)
#define GPa   (1.e3*MPa)
#define mol   InternationalSystemOfUnits_OneMole
#define sec   InternationalSystemOfUnits_OneSecond
#define kg    InternationalSystemOfUnits_OneKilogram
#define gr    (0.001*kg)



#define TEMPERATURE  (298)



/* Material properties
 * ------------------- */
#define TortuosityToLiquid            AverageModel
//#define TortuosityToLiquid            TortuosityOhJang
//#define TortuosityToLiquid(phi)            (1)
 




/* To retrieve the material properties */
#define GetProperty(a)                   (Element_GetProperty(el)[pm(a)])



/* Intern Functions */
static Model_ComputePropertyIndex_t  pm ;
static void    GetProperties(Element_t*) ;

static double* ComputeVariables(Element_t*,double**,double*,double,double,int) ;
//static Model_ComputeSecondaryVariables_t    ComputeSecondaryVariables ;
static int     ComputeSecondaryVariables(Element_t*,double,double,double*,double*) ;
static double* ComputeVariableDerivatives(Element_t*,double,double,double*,double,int) ;

static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static double* ComputeVariableFluxes(Element_t*,double**,int,int) ;
static double* ComputeFluxes(Element_t*,double*,int,int) ;
static int     TangentCoefficients(Element_t*,double,double*) ;

static void    ComputePhysicoChemicalProperties(double) ;


static double AverageModel(double) ;


/* Intern variables */
static double phii ;
static double d_so3 ;
static double k_coef ;
static double molarvol ;
static double alpha ;
static double theta ;



#include "DiffusionCoefficientOfMoleculeInWater.h"



void ComputePhysicoChemicalProperties(double TK)
{
  /* Diffusion Coefficient Of Molecules In Water (m2/s) */
  
  //d_gamma  = DiffusionCoefficientOfMoleculeInWater(SO3,TK) ;
  
}






enum {
I_C_SO3  = NEQ   ,

I_N_S     ,
I_N_S_Cry ,

I_Phi          ,
I_Last
} ;


#define NbOfVariables    (I_Last)
static double Variables[Element_MaxNbOfNodes][2*NbOfVariables] ;
static double dVariables[NbOfVariables] ;
#define Variables_n(x)    ((x) + NbOfVariables)




enum {
I_W_S      ,
I_W_Last
} ;


#define NbOfVariableFluxes    (I_W_Last)
static double VariableFluxes[Element_MaxNbOfNodes][NbOfVariableFluxes] ;



int pm(const char *s)
{
       if(strcmp(s,"porosity") == 0)     return (0) ;
  else if(strcmp(s,"diff_coef") == 0)    return (1) ;
  else if(strcmp(s,"kin_coef") == 0)     return (2) ;
  else if(strcmp(s,"molar_vol") == 0)    return (3) ;
  else if(strcmp(s,"alpha") == 0)        return (4) ;
  else if(strcmp(s,"theta") == 0)        return (5) ;
  else return(-1) ;
}




void GetProperties(Element_t* el)
{
  phii     = GetProperty("porosity") ;
  d_so3    = GetProperty("diff_coef") ;
  k_coef   = GetProperty("kin_coef") ;
  molarvol = GetProperty("molar_vol") ;
  alpha    = GetProperty("alpha") ;
  theta    = GetProperty("theta") ;
  
}




int SetModelProp(Model_t* model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_sulfur  ,"sulfur") ;
  
  
#if defined (U_LogC_SO3)
  Model_CopyNameOfUnknown(model,E_sulfur ,"logc_s") ;
#else
  Model_CopyNameOfUnknown(model,E_sulfur ,"c_s") ;
#endif


  Model_GetComputePropertyIndex(model) = pm ;
  
  //Model_GetNbOfVariables(model) = NbOfVariables ;
  //Model_GetNbOfVariableFluxes(model) = NbOfVariableFluxes ;
  //Model_GetComputeSecondaryVariables(model) = ComputeSecondaryVariables ;
  
  ComputePhysicoChemicalProperties(TEMPERATURE) ;
  
  return(0) ;
}




int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_input = 6 ;
  
    
  /* Default initialization */
  /*
  {
    ComputePhysicoChemicalProperties(TEMPERATURE) ;
    
    Material_GetPropertyValue(mat,"d_solute") = d_so3 ;
  }
  */

  Material_ScanProperties(mat,datafile,pm) ;
  
  return(n_input) ;
}



int PrintModelChar(Model_t* model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n") ;
  printf("The set of equations is:\n") ;
  printf("\t- Mass balance of sulfur     (sulfur)\n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- Sulfate concentration       (c_s)\n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;


  fprintf(ficd,"porosity = 0.38   # Porosity\n") ;

  return(NEQ) ;
}



int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
{
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = (Element_IsSubmanifold(el)) ? 0 : NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
  return(0) ;
}



int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/* Residu du aux chargements (r) */
{
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;

  {
    double* r1 = FVM_ComputeSurfaceLoadResidu(fvm,cg,t,dt) ;
    
    for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r1[i] ;
  }
  
  return(0) ;
}



int ComputeInitialState(Element_t* el)
/* Initialise les variables du systeme (f,va) */ 
{
  double* f  = Element_GetImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
  /* Pre-initialization */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
    
      /* Back up */
      N_S_Cry(i) = 0 ;
      Phi(i)     = phii ;
    }
  }
  
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      /* Variables */
      double* x  = ComputeVariables(el,u,f,0,0,i) ;
    
      /* Back up */
      N_S(i)     = x[I_N_S] ;
      N_S_Cry(i) = x[I_N_S_Cry] ;
      Phi(i)     = x[I_Phi] ;
    }
  }
  
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /* Coefficient de transfert */
  ComputeTransferCoefficients(el,u,f) ;


  /* Flux */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = i + 1 ; j < nn ; j++) {
        double* w = ComputeVariableFluxes(el,u,i,j) ;
        
        W_S(i,j)    = w[I_W_S] ;
        
        W_S(j,i)    = - w[I_W_S] ;
      }
    }
  }
  
  return(0) ;
}



int  ComputeExplicitTerms(Element_t* el,double t)
/* Thermes explicites (va)  */
{
  double*  f = Element_GetPreviousImplicitTerm(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  /*
    Coefficients de transfert
  */
  ComputeTransferCoefficients(el,u,f) ;

  return(0) ;
}



int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double* f   = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
  /* Contenus molaires */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      /* Variables */
      double* x      = ComputeVariables(el,u,f_n,t,dt,i) ;
    
      /* Back up */
      N_S(i)     = x[I_N_S] ;
      N_S_Cry(i) = x[I_N_S_Cry] ;
      Phi(i)     = x[I_Phi] ;

      {
        double x_so3  = x[I_C_SO3] ;
      
        if(x_so3 < 0) {
          double x0 = Element_GetNodeCoordinate(el,i)[0] ;
        
          printf("\n") ;
          printf("en x     = %e\n",x0) ;
          printf("x_solute = %e\n",x_so3) ;
          return(1) ;
        }
      }
    }
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;


  /* Flux */
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = i + 1 ; j < nn ; j++) {
        double* w = ComputeVariableFluxes(el,u,i,j) ;
        
        W_S(i,j)    = w[I_W_S] ;
        
        W_S(j,i)    = - w[I_W_S] ;
      }
    }
  }

  return(0) ;
}



int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  double c[ndof*ndof] ;
  int    i ;
  
  /*
    Initialisation 
  */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  TangentCoefficients(el,dt,c) ;
  {
    double* km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }

#if defined (U_LogC_SO3)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,E_sulfur)     *= Ln10*C_SO3(0) ;
      K(i,E_sulfur+NEQ) *= Ln10*C_SO3(1) ;
    }
  }
#endif


  return(0) ;

#undef K
}



int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double* f   = Element_GetCurrentImplicitTerm(el) ;
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
    Conservation of element S: (N_S - N_Sn) + dt * div(W_S) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_S,NW_Sn,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_sulfur) -= r1[i] ;
    }
  }

  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/* Les valeurs exploitees (s) */
{
  double* f = Element_GetCurrentImplicitTerm(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    nso = 5 ;
  int    i ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  

  /* Initialization */
  for(i = 0 ; i < nso ; i++) {
    Result_SetValuesToZero(r + i) ;
  }

  {
    int j = FVM_FindLocalCellIndex(fvm,s) ;
    /* molarites */
    double* x = ComputeVariables(el,u,f,t,0,j) ;


    /* Output quantities */
    i = 0 ;
    
    Result_Store(r + i++,x + I_C_SO3,"sulfate concentration",1) ;
    
    Result_Store(r + i++,x + I_N_S,"sulfur content",1) ;
    
    Result_Store(r + i++,x + I_N_S_Cry,"crystal content",1) ;
    
    Result_Store(r + i++,x + I_Phi,"porosity",1) ;
    
    {
      double w_so3 = W_S(0,1) ;
      
      Result_Store(r + i++,&w_so3,"sulfate flow",1) ;
    }
    
  }
  
  
  if(i != nso) arret("ComputeOutputs") ;
  return(nso) ;
}



void ComputeTransferCoefficients(Element_t* el,double** u,double* f)
/* Termes explicites (va)  */
{
  double* va = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int    i ; 

  /* initialization */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;
  
  
  for(i = 0 ; i < nn ; i++) {
    double* x = ComputeVariables(el,u,f,0,0,i) ;
    
    /* Liquid tortuosity */
    {
      double phi    = x[I_Phi] ;
      double tau    = TortuosityToLiquid(phi) ;
 
      KF_SO3[i]  = d_so3 * tau ;
    }
  }
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


double* ComputeFluxes(Element_t* el,double* grd,int i,int j)
{
  double* va = Element_GetExplicitTerm(el) ;
  double* w  = VariableFluxes[i] ;

  /* Gradients */
  double grd_so3  = grd[I_C_SO3] ;
      
  /* Transfer terms */
  double kf  = 0.5 * (KF_SO3[i]   + KF_SO3[j]) ;

    /* Flux */
  double w_so3  = - kf * grd_so3  ;


  w[I_W_S] = w_so3 ;
    
  return(w) ;
}




int TangentCoefficients(Element_t* el,double dt,double* c)
/**  Tangent matrix coefficients (c) */
{
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  ObVal_t* obval = Element_GetObjectiveValue(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
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
    double* xi         = ComputeVariables(el,u,f_n,0,dt,i) ;
    int k ;
    
    dui[E_sulfur   ] =  1.e-3*ObVal_GetValue(obval + E_sulfur) ;
    
    #if defined (U_LogC_SO3)
    dui[E_sulfur  ] *=  C_SO3n(i) ;
    #endif
    

    for(k = 0 ; k < NEQ ; k++) {
      double dui_k   = dui[k] ;
      double* dxi    = ComputeVariableDerivatives(el,0,dt,xi,dui_k,k) ;
    
      /* Content terms at node i */
      {
        double* cii = c + (i*nn + i)*NEQ*NEQ ;
    
        cii[E_sulfur*NEQ   + k] = dxi[I_N_S] ;
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
      
              cij[E_sulfur*NEQ   + k] = - dtdij * dw[I_W_S] ;
            }
          }
        }
      }
    }
  }

  return(dec) ;
}



double* ComputeVariables(Element_t* el,double** u,double* f_n,double t,double dt,int n)
{
  double* v0 = Element_GetConstantTerm(el) ;
  double* x = Variables[n] ;
  double* x_n = Variables_n(x) ;
  
  /* Primary Variables */
  x[E_sulfur   ] = C_SO3(n) ;
  
  /* Needed variables to compute secondary components */
  x_n[I_N_S_Cry] = N_S_Cryn(n) ;
  x_n[I_Phi]   = Phin(n) ;
  
  ComputeSecondaryVariables(el,t,dt,x_n,x) ;
  return(x) ;
}



double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double dui,int i)
{
  double* x_n = Variables_n(x) ;
  double* dx = dVariables ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }

  /* Needed variables to compute secondary components */
  
  /* We increment the variable as (x + dx) */
  dx[i] += dui ;
  
  ComputeSecondaryVariables(el,t,dt,x_n,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dui ;
  }

  return(dx) ;
}



int  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x_n,double* x)
{
  double x_so3  = x[E_sulfur   ] ;
  
  
  
  /* Backup */
  double x_so3_l = x_so3 ;
  
  /* Porosity */
  double phin = x_n[I_Phi] ;

  double n_s_cryn = x_n[I_N_S_Cry] ;
  double n_s_cry  = n_s_cryn + dt * k_coef * phin * x_so3_l ;
  
  /* Actualize the porosity */
  double phi      = phii - molarvol * n_s_cry ;
  
  
  /* Liquid contents */
  double phi_l    = phi ;
  /* ... as elements */
  double n_so3_l = phi_l*x_so3_l ;
       


  /* Backup */
  
  /* Liquid components */
  x[I_C_SO3] = x_so3 ;
  
  /* Porosity */
  x[I_Phi     ] = phi ;
  
  
  /* Element contents */
  x[I_N_S]     = n_so3_l + n_s_cry ;
  x[I_N_S_Cry] = n_s_cry ;
    
  return(0) ;
}


double AverageModel(double phi)
{
  double tau = 1 - (1 - (1 - theta) * alpha) * (phii - phi) / phii ;
  
  return tau ;
}
