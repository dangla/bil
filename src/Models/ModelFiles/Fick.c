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

#define TITLE   "Non linear Fick's second law (2015)"
#define AUTHORS ""

#include "PredefinedMethods.h"




/* Indices of equations/unknowns */
enum {
  E_solute   ,
  E_Last
} ;


/* Nb of equations */
#define NbOfEquations      E_Last
#define NEQ                NbOfEquations




/* Value of the nodal unknown (u, u_n and el must be used below) */
#define UNKNOWN(n,i)     Element_GetValueOfNodalUnknown(el,u,n,i)
#define UNKNOWNn(n,i)    Element_GetValueOfNodalUnknown(el,u_n,n,i)


/* Generic names of nodal unknowns */
#define U_solute(n)      (UNKNOWN(n,E_solute))
#define Un_solute(n)     (UNKNOWNn(n,E_solute))




/* Method chosen at compiling time. 
 * Each equation is associated to a specific unknown.
 * Each unknown can deal with a specific model.
 * Uncomment/comment to let only one unknown per equation */

/* solute: unknown either C or logC */
#define U_C_solute     U_solute
//#define U_LogC_solute  U_solute





/* Names of nodal unknowns */
#if defined (U_LogC_solute)
  #define LogC_solute(n)    U_solute(n)
  #define LogC_soluten(n)   Un_solute(n)
  #define C_solute(n)       (pow(10,LogC_solute(n)))
  #define C_soluten(n)      (pow(10,LogC_soluten(n)))
#elif defined (U_C_solute)
  #define C_solute(n)       U_solute(n)
  #define C_soluten(n)      Un_solute(n)
  #define LogC_solute(n)    (log10(C_solute(n)))
  #define LogC_soluten(n)   (log10(C_soluten(n)))
#else
  #error "Ambiguous or undefined unknown"
#endif




/* Nb of nodes (el must be used below) */
#define NN     Element_GetNbOfNodes(el)


/* Nb of terms per point */
#define NVE       (NN)
#define NVI       (NN*NN)
#define NV0       (0)


/* Names used for implicit terms */
#define MassAndFlux(f,i,j)  ((f)[((i)*NN + (j))])

#define NW_solute         (f)
#define NW_soluten        (f_n)
#define N_solute(i)       MassAndFlux(NW_solute,i,i)
#define N_soluten(i)      MassAndFlux(NW_soluten,i,i)
#define W_solute(i,j)     MassAndFlux(NW_solute,i,j)


/* Names used for explicit terms */
#define TransferCoefficient(va,i)  ((va) + (i)*NN)

#define KF_solute          TransferCoefficient(va,0)


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




/* To retrieve the material properties */
#define GetProperty(a)                   (Element_GetProperty(el)[pm(a)])



/* Intern Functions */
static Model_ComputePropertyIndex_t  pm ;
static void    GetProperties(Element_t*) ;

static double* ComputeVariables(Element_t*,double**,double*,double,double,int) ;
static int     ComputeSecondaryVariables(Element_t*,double,double,double*) ;
static double* ComputeVariableDerivatives(Element_t*,double,double,double*,double,int) ;

static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static double* ComputeVariableFluxes(Element_t*,double**,int,int) ;
static double* ComputeFluxes(Element_t*,double*,int,int) ;
static int     TangentCoefficients(Element_t*,double,double*) ;

static void    ComputePhysicoChemicalProperties(double) ;


static double TortuosityOhJang(double) ;
static double TortuosityBazantNajjar(double) ;

#define LiquidTortuosity  TortuosityOhJang
//#define LiquidTortuosity  TortuosityBazantNajjar


/* Intern variables */
static double phii ;

static double d_solute ;



#include "DiffusionCoefficientOfMoleculeInWater.h"



void ComputePhysicoChemicalProperties(double TK)
{
  /* Diffusion Coefficient Of Molecules In Water (m2/s) */
  
  d_solute  = DiffusionCoefficientOfMoleculeInWater(Na,TK) ;
  
}






enum {
I_C_solute  = NEQ   ,

I_N_solute     ,

I_Phi          ,
I_Last
} ;


#define NbOfVariables    (I_Last)
static double Variables[Element_MaxNbOfNodes][NbOfVariables] ;
static double dVariables[NbOfVariables] ;




enum {
I_W_solute      ,
I_W_Last
} ;


#define NbOfVariableFluxes    (I_W_Last)
static double VariableFluxes[Element_MaxNbOfNodes][NbOfVariableFluxes] ;



int pm(const char *s)
{
       if(strcmp(s,"porosity") == 0)     return (0) ;
  else if(strcmp(s,"d_solute") == 0)     return (1) ;
  else return(-1) ;
}




void GetProperties(Element_t* el)
{
  phii     = GetProperty("porosity") ;
  d_solute = GetProperty("d_solute") ;
}




int SetModelProp(Model_t* model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_solute  ,"solute") ;
  
  
#if defined (U_LogC_solute)
  Model_CopyNameOfUnknown(model,E_solute ,"logc_na") ;
#else
  Model_CopyNameOfUnknown(model,E_solute ,"c_na") ;
#endif


  Model_GetComputePropertyIndex(model) = pm ;
  
  ComputePhysicoChemicalProperties(TEMPERATURE) ;
  
  return(0) ;
}




int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 2 ;
  
    
  /* Default initialization */
  {
    ComputePhysicoChemicalProperties(TEMPERATURE) ;
    
    Material_GetPropertyValue(mat,"d_solute") = d_solute ;
  }

  Material_ScanProperties(mat,datafile,pm) ;
  
  return(n_donnees) ;
}



int PrintModelChar(Model_t* model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n") ;
  printf("The set of equations is:\n") ;
  printf("\t- Mass balance of solute     (solute)\n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- Solute concentration       (c_na)\n") ;
  
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
  
  {
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      /* Variables */
      double* x  = ComputeVariables(el,u,f,0,0,i) ;
    
      /* Back up */
      N_solute(i) = x[I_N_solute] ;
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
        
        W_solute(i,j)    = w[I_W_solute] ;
        
        W_solute(j,i)    = - w[I_W_solute] ;
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
      N_solute(i) = x[I_N_solute] ;

      {
        double x_solute  = x[I_C_solute] ;
      
        if(x_solute < 0) {
          double x0 = Element_GetNodeCoordinate(el,i)[0] ;
        
          printf("\n") ;
          printf("en x     = %e\n",x0) ;
          printf("x_solute = %e\n",x_solute) ;
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
        
        W_solute(i,j)    = w[I_W_solute] ;
        
        W_solute(j,i)    = - w[I_W_solute] ;
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
  double c[4*NEQ*NEQ] ;
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

#if defined (U_LogC_solute)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,E_solute)     *= Ln10*C_solute(0) ;
      K(i,E_solute+NEQ) *= Ln10*C_solute(1) ;
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
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i = 0 ; i < NEQ*nn ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;

  
  /*
    Conservation of element Na: (N_solute - N_soluten) + dt * div(W_solute) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_solute,NW_soluten,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_solute) -= r1[i] ;
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
  int    nso = 2 ;
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
    
    Result_Store(r + i++,x + I_C_solute,"solute concentration",1) ;
    
    {
      double n_solute = 0.5*(N_solute(0) + N_solute(1)) ;
      
      Result_Store(r + i++,&n_solute,"solute content",1) ;
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
      double tau    = LiquidTortuosity(phi) ;
 
      KF_solute[i]  = d_solute*tau ;
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
  double grd_solute  = grd[I_C_solute] ;
      
  /* Transfer terms */
  double kf  = 0.5 * (KF_solute[i]   + KF_solute[j]) ;

    /* Flux */
  double w_solute  = - kf * grd_solute  ;


  w[I_W_solute] = w_solute ;
    
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
    
    dui[E_solute   ] =  1.e-3*ObVal_GetValue(obval + E_solute) ;
    
    #if defined (U_LogC_solute)
    dui[E_solute  ] *=  C_soluten(i) ;
    #endif
    

    for(k = 0 ; k < NEQ ; k++) {
      double dui_k   = dui[k] ;
      double* dxi    = ComputeVariableDerivatives(el,0,dt,xi,dui_k,k) ;
    
      /* Content terms at node i */
      {
        double* cii = c + (i*nn + i)*NEQ*NEQ ;
    
        cii[E_solute*NEQ   + k] = dxi[I_N_solute] ;
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
      
              cij[E_solute*NEQ   + k] = - dtdij * dw[I_W_solute] ;
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
  
  /* Primary Variables */
  x[E_solute   ] = C_solute(n) ;
  
  /* Needed variables to compute secondary components */
  
  ComputeSecondaryVariables(el,t,dt,x) ;
  return(x) ;
}



double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double dui,int i)
{
  double* dx = dVariables ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }

  /* Needed variables to compute secondary components */
  
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



int  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x)
{
  double x_solute  = x[E_solute   ] ;
  
  
  
  /* Backup */
  double x_solute_l = x_solute ;

  
  
  /* Porosity */
  double phi      = phii ;
  
  
  /* Liquid contents */
  double phi_l    = phi ;
  /* ... as elements */
  double n_solute_l = phi_l*x_solute_l ;
       


  /* Backup */
  
  /* Liquid components */
  x[I_C_solute] = x_solute ;
  
  /* Porosity */
  x[I_Phi     ] = phi ;
  
  
  /* Element contents */
  x[I_N_solute]  = n_solute_l ;
    
  return(0) ;
}






double TortuosityOhJang(double phi)
/* Ref:
 * Byung Hwan Oh, Seung Yup Jang, 
 * Prediction of diffusivity of concrete based on simple analytic equations, 
 * Cement and Concrete Research 34 (2004) 463 - 480.
 * 
 * tau = tau_paste * tau_agg
 * 
 * tau_agg   = ratio of diffusivity in concrete and diffusivity in matrix (cement paste)
 * tau_paste = ratio of diffusivity in cement paste and diffusivity in liquid bulk
 * 
 * tau_paste = (m_p + sqrt(m_p**2 + phi_c/(1 - phi_c) * (Ds/D0)**(1/n)))**n
 * m_p = 0.5 * ((phi_cap - phi_c) + (Ds/D0)**(1/n) * (1 - phi_c - phi_cap)) / (1 - phi_c)
 * phi_cap = capillary porosity of the cement paste 
 * 
 * tau_agg = 1 + V_a/(1/(2*D_i*eps - 1) + (1 - V_a)/3)
 * V_a = Aggregate volume fraction
 * D_i = Diffusivity of the ITZ relative to that of cement paste 
 * eps = thickness ratio of ITZ: t/r_a 
 * t   = thickness of the ITZ
 * r_a = radius of the aggregate 
 * D_i = 7 
 * eps = 0.002 (Concrete)  ; 0.02 (Mortar) 
 * V_a = 0.67  (Concrete)  ; 0.45 (Mortar) 
 * tau_agg = 0.27 (Concrete) ; 0.63 (Mortar)
 */
{
  double v_a     = 0.5 ;
  double phi_cap = (phi > 0) ? (1 - v_a) * phi : 0  ;
  double phi_c = 0.18 ;          /* Critical porosity */
  double n     = 2.7 ;           /* n  = 2.7   (OPC) ; 4.5  (OPC + 10% Silica fume) */
  double ds    = 2.e-4 ;         /* ds = 2.e-4 (OPC) ; 5e-5 (OPC + 10% Silica fume) */
  double dsn   = pow(ds,1/n) ;
  double m_phi = 0.5 * ((phi_cap - phi_c) + dsn * (1 - phi_c - phi_cap)) / (1 - phi_c) ;
  double tau_paste = pow(m_phi + sqrt(m_phi*m_phi + dsn * phi_c/(1 - phi_c)),n) ;
  double tau_agg = 0.27 ;
  double tausat =  tau_paste * tau_agg ;
  
  double tau =  tausat ;
    
  return(tau) ;
}



double TortuosityBazantNajjar(double phi)
/** Liquid totuosity according to Bazant and Najjar */
{
  double iff = (phi < 0.8) ? 2.9e-4*exp(9.95*phi) : phi ;
    
  return(iff) ;
}
