/* General features of the model:
 * Simple diffusion of alkalis (as sodium or potassium):
 * (Na[+],K[+])
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


/* Nb of equations of the model */
#define NEQ    	  (1)


/* Nb of nodes (el must be used below) */
#define NN     Element_GetNbOfNodes(el)


/* Nb of terms per point */
#define NVE    	  (NN)
#define NVI       (NN*NN)
#define NV0       (0)


/* Indices of equations */
#define E_Na      (0)


/* Indices of unknowns */
#define U_C_Na    (0)


#define NOLOG_U   1
#define LOG_U     2
#define Ln10      Math_Ln10
#define U_Na      NOLOG_U

/* Value of the nodal unknown (u and el must be used as pointers below) */
#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWNn(n,i)    (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])



/* We define some names for nodal unknowns */
#if (U_Na == LOG_U)
#define LogC_Na(n)    (UNKNOWN(n,U_C_Na))
#define LogC_Nan(n)   (UNKNOWNn(n,U_C_Na))
#define C_Na(n)       (pow(10,UNKNOWN(n,U_C_Na)))
#define C_Nan(n)	    (pow(10,UNKNOWNn(n,U_C_Na)))
#else
#define C_Na(n)	      (UNKNOWN(n,U_C_Na))
#define C_Nan(n)	    (UNKNOWNn(n,U_C_Na))
#define LogC_Na(n)	  (log10(UNKNOWN(n,U_C_Na)))
#define LogC_Nan(n)	  (log10(UNKNOWNn(n,U_C_Na)))
#endif


/* Names used for implicit terms */
#define MassAndFlux(f,i,j)  ((f)[((i)*NN + (j))])

#define NW_Na         (f)
#define NW_Nan        (f_n)
#define N_Na(i)       MassAndFlux(NW_Na,i,i)
#define N_Nan(i)      MassAndFlux(NW_Nan,i,i)
#define W_Na(i,j)     MassAndFlux(NW_Na,i,j)


/* Names used for explicit terms */
#define TransferCoefficient(va,i)  ((va) + (i)*NN)

#define KF_Na          TransferCoefficient(va,0)


#define TEMPERATURE  (298)


/* To retrieve the material properties */
#define GetProperty(a)                   (Element_GetProperty(el)[pm(a)])



/* Intern Functions */
static int     pm(const char *s) ;
static void    GetProperties(Element_t*) ;

static double* ComputeVariables(Element_t*,double**,double*,double,int) ;
static Model_ComputeSecondaryVariables_t    ComputeSecondaryVariables ;
static double* ComputeVariableDerivatives(Element_t*,double,double*,double,int) ;

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

static double d_na ;


#include "DiffusionCoefficientOfMoleculeInWater.h"


void ComputePhysicoChemicalProperties(double TK)
{
  /* Diffusion Coefficient Of Molecules In Water (m2/s) */
  
  d_na         = DiffusionCoefficientOfMoleculeInWater(Na,TK) ;
  
}

#define NbOfVariables    (4)
static double Variables[Element_MaxNbOfNodes][NbOfVariables] ;
static double dVariables[NbOfVariables] ;

#define I_C_Na         (1)

#define I_N_Na         (2)

#define I_Phi          (3)



#define NbOfVariableFluxes    (1)
static double VariableFluxes[Element_MaxNbOfNodes][NbOfVariableFluxes] ;

#define I_W_Na          (0)



int pm(const char *s)
{
  if(strcmp(s,"porosity") == 0)     return (0) ;
  else return(-1) ;
}


void GetProperties(Element_t* el)
{
  phii     = GetProperty("porosity") ;
}


int SetModelProp(Model_t* model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_Na  ,"sodium") ;
  
  
#if (U_Na == LOG_U)
  Model_CopyNameOfUnknown(model,U_C_Na ,"logc_na") ;
#else
  Model_CopyNameOfUnknown(model,U_C_Na ,"c_na") ;
#endif
  
  Model_GetNbOfVariables(model) = NbOfVariables ;
  Model_GetNbOfVariableFluxes(model) = NbOfVariableFluxes ;
  Model_GetComputeSecondaryVariables(model) = ComputeSecondaryVariables ;
  
  ComputePhysicoChemicalProperties(TEMPERATURE) ;
  
  return(0) ;
}


int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 1 ;

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
  printf("\t- Mass balance of Na     (sodium)\n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- Sodium concentration             (c_na)\n") ;
  
  printf("\n") ;
  printf("PAY ATTENTION to units : \n") ;
  printf("\t length    : dm !\n") ;
  printf("\t time      : s !\n") ;
  printf("\t pressure  : Pa !\n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;


  fprintf(ficd,"porosity = 0.38   # Porosity\n") ;
  fprintf(ficd,"Curves = my_file  # File name: p_c S_l k_rl C/S H/S V_csh\n") ;  

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
      double* x       = ComputeVariables(el,u,f,0,i) ;
    
      /* Back up */
      N_Na(i) = x[I_N_Na] ;
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
        
        W_Na(i,j)    = w[I_W_Na] ;
        
        W_Na(j,i)    = - w[I_W_Na] ;
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
      double* x      = ComputeVariables(el,u,f_n,dt,i) ;
    
      /* Back up */
      N_Na(i) = x[I_N_Na] ;

      {
        double x_na    	  = x[I_C_Na] ;
      
        if(x_na < 0) {
          double x0 = Element_GetNodeCoordinate(el,i)[0] ;
        
          printf("\n") ;
          printf("en x     = %e\n",x0) ;
          printf("x_na     = %e\n",x_na) ;
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
        
        W_Na(i,j)    = w[I_W_Na] ;
        
        W_Na(j,i)    = - w[I_W_Na] ;
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

#if (U_Na == LOG_U)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < 2*NEQ ; i++){
      K(i,U_C_Na)     *= Ln10*C_Na(0) ;
      K(i,U_C_Na+NEQ) *= Ln10*C_Na(1) ;
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
    Conservation of element Na: (N_Na - N_Nan) + dt * div(W_Na) = 0
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_Na,NW_Nan,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Na) -= r1[i] ;
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
    double* x = ComputeVariables(el,u,f,0,j) ;

    double x_na    	  = x[I_C_Na] ;

    double n_Na = 0.5*(N_Na(0) + N_Na(1)) ;
    


    /* quantites exploitees */
    i = 0 ;
    Result_Store(r + i++,&x_na,"x_na",1) ;
    Result_Store(r + i++,&n_Na,"n_Na",1) ;
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
    double* x = ComputeVariables(el,u,f,0,i) ;
    
    /* Liquid tortuosity */
    {
      double phi    = x[I_Phi] ;
      double iff    = LiquidTortuosity(phi) ;
        
      //TORTUOSITY[i] = iff ;
 
      KF_Na[i]     	= d_na*iff ;
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
  double grd_na       = grd[I_C_Na] ;
      
  /* Transfer terms */
  double kf_na   = 0.5 * (KF_Na[i]   + KF_Na[j]) ;

    /* Flux */
  double w_na    = - kf_na * grd_na  ;


  w[I_W_Na] = w_na ;
    
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
    double* xi         = ComputeVariables(el,u,f_n,dt,i) ;
    int k ;
    
    dui[U_C_Na   ] =  1.e-3*ObVal_GetValue(obval + U_C_Na) ;
    
    #if (U_Na == LOG_U)
    dui[U_C_Na  ] *=  C_Nan(i) ;
    #endif
    

    for(k = 0 ; k < NEQ ; k++) {
      double dui_k    = dui[k] ;
      double* dxi    = ComputeVariableDerivatives(el,dt,xi,dui_k,k) ;
    
      /* Content terms at node i */
      {
        double* cii = c + (i*nn + i)*NEQ*NEQ ;
    
        cii[E_Na*NEQ   + k] = dxi[I_N_Na] ;
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
      
              cij[E_Na*NEQ   + k] = - dtdij * dw[I_W_Na] ;
            }
          }
        }
      }
    }
  }

  return(dec) ;
}



double* ComputeVariables(Element_t* el,double** u,double* f_n,double dt,int n)
{
  double* v0 = Element_GetConstantTerm(el) ;
  double* x = Variables[n] ;
  
  /* Primary Variables */
  x[U_C_Na   ] = C_Na(n) ;
  
  /* Needed variables to compute secondary components */
  
  ComputeSecondaryVariables(el,dt,x) ;
  return(x) ;
}


double* ComputeVariableDerivatives(Element_t* el,double dt,double* x,double dui,int i)
{
  double* dx = dVariables ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NEQ ; j++) {
    dx[j] = x[j] ;
  }

  /* Needed variables to compute secondary components */
  
  /* We increment the variable as (x + dx) */
  dx[i] += dui ;
  
  ComputeSecondaryVariables(el,dt,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dui ;
  }

  return(dx) ;
}



void  ComputeSecondaryVariables(Element_t* el,double dt,double* x)
{
  double x_na       = x[U_C_Na   ] ;
  
  
  
  /* Backup */
  double x_na_l = x_na ;

  
  
  /* Porosity */
  double phi      = phii ;
  
  
  /* Liquid contents */
  double phi_l    = phi ;
  /* ... as elements */
  double n_na_l = phi_l*x_na_l ;
       


  /* Backup */
  
  /* Liquid components */
  x[I_C_Na      ] = x_na ;
  
  /* Porosity */
  x[I_Phi     ] = phi ;
  
  
  /* Element contents */
  x[I_N_Na]  = n_na_l ;
    
  return ;
}



double TortuosityOhJang(double phi)
/** Liquid totuosity according to Oh and Jang, CCR2003 */
{
  double phi_cap = 0.5*phi  ;
  double phi_c   = 0.17 ;  /* Percolation capilar porosity */
    
  double n = 2.7 ; 		      /* OPC n = 2.7,        Fly ash n = 4.5 */
  double ds_norm = 1.e-4 ;	/* OPC ds_norm = 1e-4, Fly ash ds_norm = 5e-5 */
  double m_phi = 0.5*(pow(ds_norm,1/n) + phi_cap/(1-phi_c)*(1 - pow(ds_norm,1/n)) - phi_c/(1-phi_c)) ;
  double iff =  pow(m_phi + sqrt(m_phi*m_phi +  pow(ds_norm,1/n)*phi_c/(1-phi_c)),n) ;
    
  return(iff) ;
}



double TortuosityBazantNajjar(double phi)
/** Liquid totuosity according to Bazant and Najjar */
{
  double iff = (phi < 0.8) ? 2.9e-4*exp(9.95*phi) : phi ;
    
  return(iff) ;
}
