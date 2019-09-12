#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "CommonModel.h"

/* Choose the numerical method */
#include "FVM.h"

#define TITLE "Freezing and thawing of concrete with salt"
#define AUTHORS "Zeng"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ      (3)


/* Nb of terms */
#define NVI      (3*NN*NN)
#define NVE      (7*NN)
#define NV0      (0)

/* Indices of equations */
#define E_h2o    (0)
#define E_salt   (1)
#define E_the    (2)


/* Indices of unknowns */
#define U_p_l    (0)
#define U_c_s    (1)
#define U_tem    (2)


/* Value of the nodal unknown (u and el must be used below) */
#define UNKNOWN(n,i)     Element_GetValueOfNodalUnknown(el,u,n,i)
#define UNKNOWNn(n,i)    Element_GetValueOfNodalUnknown(el,u_n,n,i)


/* Names of nodal unknowns */
#define P_l(n)   (UNKNOWN(n,U_p_l))
#define C_s(n)   (UNKNOWN(n,U_c_s))
#define TEM(n)   (UNKNOWN(n,U_tem))

#define P_ln(n)   (UNKNOWNn(n,U_p_l))
#define C_sn(n)   (UNKNOWNn(n,U_c_s))
#define TEM_n(n)  (UNKNOWNn(n,U_tem))


/* Nb of nodes (el must be used below) */
#define NN     Element_GetNbOfNodes(el)


/* Names used for implicit terms */
#define MASSFLUXES(f,i,j)  ((f)[((i)*NN + (j))])

#define M_H2O(i)    MASSFLUXES(f,i,i)
#define W_H2O(i,j)  MASSFLUXES(f,i,j)
#define M_SALT(i)   MASSFLUXES(f + NN*NN,i,i)
#define W_SALT(i,j) MASSFLUXES(f + NN*NN,i,j)
#define S(i)        MASSFLUXES(f + 2*NN*NN,i,i)
#define Q(i,j)      MASSFLUXES(f + 2*NN*NN,i,j)

#define M_H2On(i)    MASSFLUXES(f_n,i,i)
#define M_SALTn(i)   MASSFLUXES(f_n + NN*NN,i,i)
#define S_n(i)       MASSFLUXES(f_n + 2*NN*NN,i,i)



/* Names used for explicit terms */
#define KD_H2O      (va)
#define KD_SALT     (va + NN)
#define KF_SALT     (va + 2*NN)
#define KTH         (va + 3*NN)
#define S_H2O_l     (va + 4*NN)
#define S_H2O_i     (va + 5*NN)
#define S_SALT      (va + 6*NN)



/* Water properties
 * ---------------- */
 #include "WaterViscosity.h"
 #include "Log10ActivityOfWaterInBrine.h"
 /* Molar mass (kg/m3) */
#define M_h2o      (18.e-3)
/* Molar volume of liquid water (m3/mol) */
#define V_h2o      (18.e-6)
/* Molar volume of ice (m3/mol) */
#define V_ice      (19.63e-6)
/* Viscosity (Pa.s) */
#define mu_l       (1.79e-3)
/* Logarithm of water activity in brine */
#define LogActivityOfWater(c_s,T)      activity(c_s,T)
/* Entropy of fusion of ice (J/mol/K) */
#define S_m        (23.54)
/* Melting temperature of ice at atmospheric pressure */
#define T_m        (273.)
/* Atmospheric pressure (Pa) */
#define p_m        (1.e5)
/* Bulk modulus of liquid water (Pa) */
#define K_l        (1.8e9)
/* Bulk modulus of ice (Pa) */
#define K_i        (7.8e9)
/* Specific heat of liquid water (J/kg/K) */
#define C_l        (4180.)
/* Specific heat of ice (J/kg/K) */
#define C_i        (2000.)
/* Thermal conductivity of liquid water (W/m/K) */
#define LAM_l      (0.6)
/* Thermal conductivity of ice (W/m/K) */
#define LAM_i      (2.2)
/* Thermal volumetric expansion coefficient of liquid water (1/K) */
/* #define ALPHA_l    (-298.e-6) */
#define ALPHA_l(theta)    (-68.7e-6 + 13.877e-6*(theta))  /* After Speedy (1987) */
/* Thermal volumetric expansion coefficient of ice (1/K) */
#define ALPHA_i    (155.e-6)



/* Salt properties
 * --------------- */
/* valencies */
#define z_na       (1.)
#define z_cl       (-1.)
#define z_ca       (2.)

/* Partial molar volumes (m3/mole) */
#include "PartialMolarVolumeOfMoleculeInWater.h"
#define V_na       PartialMolarVolumeOfMoleculeInWater(Na)  /* (22.47e-6) */
#define V_cl       PartialMolarVolumeOfMoleculeInWater(Cl)  /* (-0.35e-6) */
#define V_ca       PartialMolarVolumeOfMoleculeInWater(Ca)  /* (-18.7e-6) */

/* Solid molar volumes (m3/mole) */
#define V_nacl     (24.5e-6)
#define V_cacl2    (40.e-6)

/* Molar masses (kg/m3) */
#define M_na       (23.e-3)
#define M_cl       (35.45e-3)
#define M_ca       (40.1e-3)
#define M_cacl2    (110.99e-3)
#define M_nacl     (58.45e-3)

/* Molecular diffusion coefficients (m2/s) */
#define d_ca       (7.92e-10)
#define d_na       (1.33e-9)
#define d_cl       (2.032e-9)
#define d_nacl     (1.e-9)
#define d_cacl2    (1.e-9)


/* Types of salts */
#define NaCl   0
#define CaCl2  1

/* The chosen salt */
#define SALT   CaCl2

#if SALT == CaCl2
#define V_A      V_cl
#define V_C      V_ca
#define V_AC     (V_ca + 2*V_cl)
#define V_salt   V_cacl2
#define M_salt   M_cacl2
#define d_salt   d_cacl2
#elif SALT == NaCl
#define V_A      V_cl
#define V_C      V_na
#define V_AC     (V_na + V_cl)
#define V_salt   V_nacl
#define M_salt   M_nacl
#define d_salt   d_nacl
#else
#error "Salt not available"
#endif



/* Brine properties
 * ---------------- */


/* Material Properties
 * ------------------- */
#define SaturationDegree(x)                Curve_ComputeValue(Element_GetCurve(el),x)
#define dSaturationDegree(x)               Curve_ComputeDerivative(Element_GetCurve(el),x)
#define RelativePermeabilityToLiquid(x)    Curve_ComputeValue(Element_GetCurve(el) + 1,x)



/* Physical constants
 * ------------------ */
/* Perfect gas constant (J/mol/K) */
#define R_g        (8.315)


/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])

/* Fonctions */
static int    pm(const char *s) ;
static void    GetProperties(Element_t*) ;

static double* ComputeVariables(Element_t*,double**,double**,double*,double,double,int) ;
static void    ComputeSecondaryVariables(Element_t*,double,double,double*) ;
static double* ComputeVariableDerivatives(Element_t*,double,double,double*,double,int) ;


static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static double* ComputeVariableFluxes(Element_t*,int,int) ;
static double* ComputeFluxes(Element_t*,double*,int,int) ;

static int    TangentCoefficients(Element_t*,double,double,double*) ;
static double activity(double,double) ;
static double tortuosity(double,double) ;
static double lna_i(double,double,double,double,double,double) ;
static double lng_LinLee(double,double,double,double,double,double) ;
static double lng_TQN(double,double,double,double,double,double,double,double) ;

/* Parameters */
static double lam_l = LAM_l ;
static double lam_i = LAM_i ;
static double rho_h2o_i0 = M_h2o/V_ice ;
static double rho_h2o_l0 = M_h2o/V_h2o ;
static double alpha_i = ALPHA_i ;
static double vr_salt = (V_AC/M_salt - V_h2o/M_h2o) ;
static double phi,k_int,lam_s,C_s,k_s,g_s,alpha_s ;



#define NbOfVariables    (NEQ+12)
static double Variables[Element_MaxNbOfNodes][NbOfVariables] ;
static double dVariables[NbOfVariables] ;


#define I_M_H2O      (NEQ+0)
#define I_M_SALT     (NEQ+1)
#define I_S          (NEQ+2)
#define I_S_H2O_l    (NEQ+3)
#define I_S_H2O_i    (NEQ+4)
#define I_S_SALT     (NEQ+5)
#define I_S_L        (NEQ+6)
#define I_S_I        (NEQ+7)
#define I_P_L        (NEQ+8)
#define I_P_I        (NEQ+9)
#define I_RHO_l      (NEQ+10)
#define I_RHO_i      (NEQ+11)
  
  

#define NbOfVariableFluxes    (3)
static double VariableFluxes[NbOfVariableFluxes] ;

#define I_W_H2O           (0)
#define I_W_SALT          (1)
#define I_Q               (2)





int pm(const char *s)
{
  if(strcmp(s,"phi") == 0) return (0) ;
  else if(strcmp(s,"k_int") == 0) return (1) ;
  else if(strcmp(s,"C_s") == 0) return (2) ;
  else if(strcmp(s,"lam_s") == 0) return (3) ;
  else if(strcmp(s,"k_s") == 0) return (4) ;
  else if(strcmp(s,"g_s") == 0) return (5) ;
  else if(strcmp(s,"alpha_s") == 0) return (6) ;
  else return (-1) ;
}





void GetProperties(Element_t* el)
{
  phi     = GetProperty("phi") ;
  k_int   = GetProperty("k_int") ;
  C_s     = GetProperty("C_s") ;
  lam_s   = GetProperty("lam_s") ;
  k_s     = GetProperty("k_s") ;
  g_s     = GetProperty("g_s") ;
  alpha_s = GetProperty("alpha_s") ;
}



int SetModelProp(Model_t *model)
/** Set the model properties 
 *  Return 0 */
{
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  Model_CopyNameOfEquation(model,E_h2o,"water") ;
  Model_CopyNameOfEquation(model,E_salt,"salt") ;
  Model_CopyNameOfEquation(model,E_the,"the") ;
  
  /** Names of the main (nodal) unknowns */
  Model_CopyNameOfUnknown(model,U_p_l,"p_l") ;
  Model_CopyNameOfUnknown(model,U_c_s,"c_s") ;
  Model_CopyNameOfUnknown(model,U_tem,"tem") ;
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 7 ;
  
  Material_ScanProperties(mat,datafile,pm) ;
  
  return(NbOfProp) ;
}




int PrintModelProp(Model_t *model,FILE *ficd)
/** Print the model properties 
 *  Return the nb of equations */
{ 
  printf(TITLE) ;
  
  if(!ficd) return(0) ;
  
  printf("\n\
This model consists in 3 equations:\n\
\t 1. Mass balance of water  (water)\n\
\t 2. Masse balance of salt  (salt)\n\
\t 2. Entropy balance        (the)\n") ;
  
  
  printf("\n\
The primary unknowns are:\n\
\t 1. The liquid pressure    (p_l)\n\
\t 2. The salt concentration (c_s)\n\
\t 2. The temperature        (tem)\n") ;


  printf("\n\
Example of input data\n\n") ;

  fprintf(ficd,"porosite = 0.2   # Porosite\n") ;
  fprintf(ficd,"k_int = 5.e-21   # Permeabilite intrinseque (m2)\n") ;
  fprintf(ficd,"C_s = 2e+06      # Chaleur volumique du solide\n") ;
  fprintf(ficd,"lam_s = 1.       # Conductivite thermique du solide\n") ;
  fprintf(ficd,"k_s = 3.18e10    # Module de compression du solide\n") ;
  fprintf(ficd,"g_s = 1.91e10    # Module cisaillement du solide\n") ;
  fprintf(ficd,"alpha_s = 54.e-6 # Dilatation thermique vol. du solide\n") ;
  fprintf(ficd,"Curves = myfile # nom du fichier : p_c S_l k_rl tau_l\n") ;

  return(NEQ) ;
}


int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{

  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
  
  return(0) ;
}



int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/** Compute the residu (r) due to loads 
 *  Return 0 if succeeded and -1 if failed */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int    i ;
  

  {
    double *r1 = FVM_ComputeSurfaceLoadResidu(fvm,cg,t,dt) ;
    for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r1[i] ;
    
    /* Thermic */
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_the) {
      R(0,E_the) /= TEM_n(0) ;
    }
  }
  
  return(0) ;
#undef R
}


int ComputeInitialState(Element_t* el)
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 *  Return 0 if succeeded and -1 if failed
 */ 
{
  double* f  = Element_GetImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int    i ;
  
  
  /* We can skip if the element is a submanifold, 
   * e.g. a surface in 3D or a line in 2D */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*Donnees */
  GetProperties(el) ;
  
  /* masses of h2o, salt and entropies */
  for(i = 0 ; i < nn ; i++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,f,0,0,i) ;

    M_H2O(i)   = x[I_M_H2O] ;
    M_SALT(i)  = x[I_M_SALT] ;
    S(i)       = x[I_S] ;
  }

  {
    ComputeTransferCoefficients(el,u,f) ;
  }

  /* Flux */
  {
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = i + 1 ; j < nn ; j++) {
        double* w = ComputeVariableFluxes(el,i,j) ;

        W_H2O(i,j)  =   w[I_W_H2O] ;
        W_SALT(i,j) =   w[I_W_SALT] ;
        Q(i,j)      =   w[I_Q] ;
        
        W_H2O(j,i)  = - w[I_W_H2O] ;
        W_SALT(j,i) = - w[I_W_SALT] ;
        Q(j,i)      = - w[I_Q] ;
      }
    }
  }
  
  return(0) ;
}



int  ComputeExplicitTerms(Element_t* el,double t)
/* Explicit terms (va)  */
{
  double* f = Element_GetPreviousImplicitTerm(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Transfer coefficients
  */
  ComputeTransferCoefficients(el,u,f) ;

  return(0) ;
}




int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double* f   = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int i ;
  

  
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Donnees
  */
  GetProperties(el) ;

  /* masses of h2o, salt and entropies */
  for(i = 0 ; i < nn ; i++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u_n,f_n,t,dt,i) ;

    M_H2O(i)   = x[I_M_H2O] ;
    M_SALT(i)  = x[I_M_SALT] ;
    S(i)       = x[I_S] ;
  }

  /* Flux */
  {
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = i + 1 ; j < nn ; j++) {
        double* w = ComputeVariableFluxes(el,i,j) ;

        W_H2O(i,j)  =   w[I_W_H2O] ;
        W_SALT(i,j) =   w[I_W_SALT] ;
        Q(i,j)      =   w[I_Q] ;
        
        W_H2O(j,i)  = - w[I_W_H2O] ;
        W_SALT(j,i) = - w[I_W_SALT] ;
        Q(j,i)      = - w[I_Q] ;
      }
    }
  }

  return(0) ;
} 





int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  int    i ;
  
  /*
    Initialisation 
  */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  {
    double c[4*NEQ*NEQ] ;
    int dec = TangentCoefficients(el,t,dt,c) ;
    double* km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }


  return(0) ;

#undef K
}



int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ + (i)])
  double* f = Element_GetCurrentImplicitTerm(el) ;
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
    MASS BALANCE FOR H2O : (m_h2o1 - m_h2on) + dt * div(w_h2o1) = 0
  */
  {
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(i == j) {
          g[i*nn + i] = M_H2O(i) - M_H2On(i) ;
        } else {
          g[i*nn + j] = dt * W_H2O(i,j) ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_h2o) -= r1[i] ;
      }
    }
  }
  /*
    MASS BALANCE FOR SALT : (m_salt1 - m_saltn) + dt * div(w_salt1) = 0
  */
  {
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(i == j) {
          g[i*nn + i] = M_SALT(i) - M_SALTn(i) ;
        } else {
          g[i*nn + j] = dt * W_SALT(i,j) ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_salt) -= r1[i] ;
      }
    }
  }
  /*
    ENTROPY BALANCE : (S_1 - S_n) + dt * div(q_1/T_n+S_i*w_i1) = 0
  */
  {
    double g[Element_MaxNbOfNodes*Element_MaxNbOfNodes] ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(i == j) {
          g[i*nn + i] = S(i) - S_n(i) ;
        } else {
          g[i*nn + j] = dt * Q(i,j) ;
        }
      }
    }
    
    {
      double* r1 = FVM_ComputeMassAndFluxResidu(fvm,g) ;
      
      for(i = 0 ; i < nn ; i++) {
        R(i,E_the) -= r1[i] ;
      }
    }
  }
  
  return(0) ;

#undef R
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
{
  double* f = Element_GetCurrentImplicitTerm(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int    nso = 16 ;
  int    i ;
  

  /* Initialization */
  for(i = 0 ; i < nso ; i++) {
    Result_SetValuesToZero(r + i) ;
  }

  /*
    Donnees
  */
  GetProperties(el) ;
  
  /* Quantities */
  {
    int j = FVM_FindLocalCellIndex(fvm,s) ;
    /* Variables */
    double* x = ComputeVariables(el,u,u,f,t,0,j) ;

    /* pressure */
    double p_l  =  x[U_p_l] ;
    /* temperature */
    double tem  =  x[U_tem] ;
    
    /* ice pressure */
    double p_i  = x[I_P_I] ;
    
    /* liquid saturation */
    double s_l  = x[I_S_L] ;
    double s_i  = 1 - s_l ;
    
    /* original Qiang !!!
    double K = k_s*(1 - phi*(1 + 3*k_s/(4*g_s))) ;
    double G = g_s*(1 - 5*phi*(4*g_s + 3*k_s)/(8*g_s + 9*k_s)) ;
    double b_i  = phi*s_i*(1 + 3*k_s/(4*g_s)) ;
    double b_l  = phi*s_l*(1 + 3*k_s/(4*g_s)) ;
    double N_ii = 3*phi*s_i/(4*g_s) ;
    double N_ll = 3*phi*s_l/(4*g_s) ; 
    double alpha_phi_l = alpha_s*(b_l - phi*s_l) ;
    double alpha_phi_i = alpha_s*(b_i - phi*s_i) ;
    double Epsi   = (b_l*p_l + b_i*p_i + alpha_s*(tem - T_m))/K ;
    double phi_l  = b_l*Epsi + p_l*N_ll - alpha_phi_l*(tem - T_m) ;
    double phi_i  = b_i*Epsi + p_i*N_ii - alpha_phi_i*(tem - T_m) ;
    double phi_p  = phi + phi_l + phi_i ;
    
    double Eps    = (b_l*(p_l - p_m) + b_i*(p_i - p_m) + alpha_s*K*(tem - T_m))/K ;
    double Eps1   = (b_l*(p_l - p_m) + b_i*(p_i - p_m))/K ;
    double Eps2   = alpha_s*K*(tem - T_m)/K ;
    double Rho    = (phi_p<1) ? (phi_p*s_l*(p_l - p_m) + phi_p*s_i*(p_i - p_m))/(1-phi_p) : 0 ;
    double Eps_T  = alpha_s*K*(tem - T_m)/K ;
    double Eps_P  = (b_l*(p_l - p_m) + b_i*(p_i - p_m))/K ;
    double v = 0.2 ;
    double Eps_xx = Eps*(1+v)/(1-v);
    double Rho_yy = 2*G*Eps_xx ;
    double Rho_zz = Rho_yy ;
    */
    
    /* coef poroelastiques */
    double K = k_s*(1 - phi*(1 + 0.75*k_s/g_s)) ;
    double G = g_s*(1 - 5*phi*(4*g_s + 3*k_s)/(8*g_s + 9*k_s)) ;
    double k_oedo = K + 4*G/3. ; /* Module oedometrique */
    double b_i  = phi*s_i*(1 + 0.75*k_s/g_s) ;
    double b_l  = phi*s_l*(1 + 0.75*k_s/g_s) ;
    double N_ii = 0.75*phi*s_i/g_s ;
    double N_ll = 0.75*phi*s_l/g_s ; 
    double alpha_phi_l = alpha_s*(b_l - phi*s_l) ;
    double alpha_phi_i = alpha_s*(b_i - phi*s_i) ;
    
    double Eps_T  = alpha_s*K*(tem - T_m)/k_oedo ;
    double Eps_P  = (b_l*(p_l - p_m) + b_i*(p_i - p_m))/k_oedo ;
    
    /* conditions Sig_xx = Eps_yy = Eps_zz = 0 */
    double Eps_xx = Eps_P + Eps_T ;
    double Sig_yy = - 2*G*Eps_xx ;
    /*
    double Sig_zz = Sig_yy ;
    */
      
    double phi_l  = b_l*Eps_xx + N_ll*(p_l - p_m) - alpha_phi_l*(tem - T_m) ;
    double phi_i  = b_i*Eps_xx + N_ii*(p_i - p_m) - alpha_phi_i*(tem - T_m) ;
    double phi_p  = phi + phi_l + phi_i ;
    
    
    
    if(Element_GetCoordinateSystem(el) != CARTESIAN) {
      Message_FatalError("Impossible") ;
    }
    
    i = 0 ;
    /* quantites exploitees */
    Result_Store(r + i++,x + U_p_l,"Liquid pressure",1) ;
    Result_Store(r + i++,x + U_c_s,"Salt concentration",1) ;
    Result_Store(r + i++,x + U_tem,"Temperature",1) ;
    Result_Store(r + i++,x + I_S_L,"Liquid saturation",1) ;
    Result_Store(r + i++,&W_H2O(0,1),"Mass flow of H2O",1) ;
    Result_Store(r + i++,&W_SALT(0,1),"Mass flow of salt",1) ;
    Result_Store(r + i++,x + I_P_I,"Ice pressure",1) ;
    Result_Store(r + i++,x + I_RHO_l,"Solution density",1) ;
    Result_Store(r + i++,x + I_RHO_i,"Ice density",1) ;
    Result_Store(r + i++,&Eps_xx,"strain",1) ;
    Result_Store(r + i++,&Sig_yy,"stress",1) ;
    Result_Store(r + i++,&Eps_T,"strain_temperature",1) ;
    Result_Store(r + i++,&Eps_P,"strain_pressure",1) ;
    Result_Store(r + i++,&phi_p,"porosity",1) ;
    Result_Store(r + i++,&phi_i,"deformation of ice pores",1) ;
    Result_Store(r + i++,&phi_l,"deformation of liquid pores",1) ;
  }

  return (nso) ;
}



void  ComputeTransferCoefficients(Element_t* el,double** u,double* f)
{
  double* va = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int i ;

  /*
    Donnees
  */
  GetProperties(el) ;

  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;

  /* Transfer coefficients */
  for(i = 0 ; i < nn ; i++) {
    double* x = ComputeVariables(el,u,u,f,0,0,i) ;
    
    double p_l = x[U_p_l] ;
    double tem = x[U_tem] ;  
    double c_s = x[U_c_s] ;
    double p_i = x[I_P_I] ;
    double p_c = p_i - p_l ;

    double s_l = x[I_S_L] ;
    double s_i = x[I_S_I] ;
    
    /* Mass densities */
    double rho_salt  = M_salt*c_s ;
    double rho_l     = x[I_RHO_l] ; 
    double rho_h2o_l = rho_l - rho_salt ;

    double k_l   = k_int/mu_l*RelativePermeabilityToLiquid(p_c) ;
    double tau_l = tortuosity(phi,s_l) ;
    double lam_h2o = s_l*lam_l + s_i*lam_i ;

    KD_H2O[i]  = rho_h2o_l*k_l ;
    KD_SALT[i] = rho_salt*k_l ;
    KF_SALT[i] = M_salt*phi*s_l*tau_l*d_salt ;
    KTH[i]     = (lam_s*(1 - 3*phi*(lam_s - lam_h2o)/(3*lam_s - (1 - phi)*(lam_s - lam_h2o))))/tem ;
    S_H2O_l[i] = x[I_S_H2O_l] ;
    S_H2O_i[i] = x[I_S_H2O_i] ;
    S_SALT[i]  = x[I_S_SALT] ;
  }
  
}



double* ComputeVariableFluxes(Element_t* el,int i,int j)
{
  int nn = Element_GetNbOfNodes(el) ;
  double* grdij = dVariables ;
  FVM_t* fvm   = FVM_GetInstance(el) ;
  double* dist = FVM_ComputeIntercellDistances(fvm) ;
  double dij  = dist[nn*i + j] ;


  {
    double* xi  = Variables[i] ;
    double* xj  = Variables[j] ;
    int k ;
      
    for(k = 0 ; k < NbOfVariables ; k++)  {
      grdij[k] = (xj[k] - xi[k])/dij ;
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
  double* w  = VariableFluxes ;

  /* Gradients */
  double grd_p_l      = grd[U_p_l] ;
  double grd_c_s      = grd[U_c_s] ;
  double grd_tem      = grd[U_tem] ;
  
  /* Transfer terms */
  double kd_h2o  = 0.5*(KD_H2O[i]  + KD_H2O[j]) ;
  double kd_salt = 0.5*(KD_SALT[i] + KD_SALT[j]) ;
  double kf_salt = 0.5*(KF_SALT[i] + KF_SALT[j]) ;
  double kth     = 0.5*(KTH[i]     + KTH[j]) ;
  double s_h2o_l = 0.5*(S_H2O_l[i] + S_H2O_l[j]) ;
  double s_salt  = 0.5*(S_SALT[i]  + S_SALT[j]) ;

  /* Fluxes */
  double w_h2o  = - kd_h2o*grd_p_l  + kf_salt*grd_c_s ;
  double w_salt = - kd_salt*grd_p_l - kf_salt*grd_c_s ;
  double q      = - kth*grd_tem + s_h2o_l*w_h2o + s_salt*w_salt ;
  

  w[I_W_H2O ]  = w_h2o ;
  w[I_W_SALT]  = w_salt ;
  w[I_Q]       = q ;
  
  return(w) ;
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
  double dxi[NEQ] ;
  FVM_t* fvm   = FVM_GetInstance(el) ;
  double* dist = FVM_ComputeIntercellDistances(fvm) ;
  int    i ;
  
  /* Initialization */
  for(i = 0 ; i < ndof*ndof ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  for(i = 0 ; i < NEQ ; i++) {
    dxi[i] =  1.e-2*ObVal_GetValue(obval + i) ;
  }

  
  for(i = 0 ; i < nn ; i++) {
    /* Variables */
    double* x   = ComputeVariables(el,u,u_n,f_n,t,dt,i) ;
    int k ;
    
    for(k = 0 ; k < NEQ ; k++) {
      double dxk    = dxi[k] ;
      double* dx    = ComputeVariableDerivatives(el,t,dt,x,dxk,k) ;
    
      /* Content terms at node i */
      {
        double* cii = c + (i*nn + i)*NEQ*NEQ ;
        
        cii[E_h2o*NEQ    + k] = dx[I_M_H2O ] ;
        cii[E_salt*NEQ   + k] = dx[I_M_SALT] ;
        cii[E_the*NEQ    + k] = dx[I_S ] ;
      }
      
      /* Transfer terms from node i to node j: d(w_ij)/d(u_i) */
      {
        int j ;
        
        for(j = 0 ; j < nn ; j++) {
          if(j != i) {
            double* dw = ComputeFluxes(el,dx,i,j) ;
            double* cij = c + (i*nn + j)*NEQ*NEQ ;
            double dij  = dist[nn*i + j] ;
            double dtdij = dt/dij ;
        
            cij[E_h2o*NEQ    + k] = - dtdij * dw[I_W_H2O ] ;
            cij[E_salt*NEQ   + k] = - dtdij * dw[I_W_SALT] ;
            cij[E_the*NEQ    + k] = - dtdij * dw[I_Q ] ;
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
  x[U_p_l] = P_l(n) ;
  x[U_tem] = TEM(n) ;
  x[U_c_s] = C_s(n) ;
  
  /* Needed variables to compute secondary components */
  
  ComputeSecondaryVariables(el,t,dt,x) ;
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x)
{
  /* Primary variables */
  double p_l    = x[U_p_l] ;
  double c_s    = x[U_c_s] ;
  double tem    = x[U_tem] ;
  
  double lna = LogActivityOfWater(c_s,tem) ;
  
  /* Pressures */
  double p_i = p_m + (V_h2o*(p_l - p_m) + S_m*(T_m - tem) + R_g*tem*lna)/V_ice ;
  double p_c = p_i - p_l ;

  /* Saturations */
  double s_l = SaturationDegree(p_c) ;
  double s_i = 1. - s_l ;

  /* Mass densities */
  double alpha_l   = ALPHA_l(tem - T_m) ;
  double rho_h2o_i = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
  double rho_salt  = M_salt*c_s ;
  double rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - vr_salt*rho_salt) ;
  double rho_h2o_l = rho_l - rho_salt ;
  
  /* Mass contents */
  double m_h2o_l = phi*s_l*rho_h2o_l ;
  double m_h2o_i = phi*s_i*rho_h2o_i ;
  double m_salt  = phi*s_l*rho_salt ;
  
  /* Entropies */
  double s_h2o_l = C_l*log(tem/T_m) - alpha_l/rho_h2o_i0*(p_l - p_m) ;
  double s_h2o_i = C_i*log(tem/T_m) - S_m/M_h2o - alpha_i/rho_h2o_i0*(p_i - p_m) ;
  double s_salt  = s_h2o_l ;
  double s_tot   = C_s*log(tem/T_m) + m_h2o_l*s_h2o_l + m_h2o_i*s_h2o_i + m_salt*s_salt ;
  
  
  /* Backup */
  x[I_M_H2O   ] = m_h2o_l + m_h2o_i ;
  x[I_M_SALT  ] = m_salt ;
  x[I_S       ] = s_tot ;
  x[I_S_H2O_l ] = s_h2o_l ;
  x[I_S_H2O_i ] = s_h2o_i ;
  x[I_S_SALT  ] = s_salt ;
  x[I_S_L     ] = s_l ;
  x[I_S_I     ] = s_i ;
  x[I_P_L     ] = p_l ;
  x[I_P_I     ] = p_i ;
  x[I_RHO_l   ] = rho_l ;
  x[I_RHO_i   ] = rho_h2o_i ;
}



double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double dxi,int i)
{
  double* dx = dVariables ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;
  
  ComputeSecondaryVariables(el,t,dt,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}







double activity(double c_s,double tem)
/* activity of water */
{
  double T_0  = T_m ;
  double b0   = sqrt(M_h2o),S0 = pow(M_h2o,1.29) ; /* references */
  /* NaCl (d'apres Lin et Lee) */
  double b_na_nacl = 4.352/b0,b_cl_nacl = 1.827/b0 ;
  double S_na_nacl = 26.448/S0,S_cl_nacl = 19.245/S0 ;
  /* CaCl2 (d'apres Lin et Lee) */
  double b_ca_cacl2 = 3.908/b0,b_cl_cacl2 = 2.085/b0 ;
  double S_ca_cacl2 = 18.321/S0,S_cl_cacl2 = 10.745/S0 ;

  double epsi = 0.0007*(tem - T_0)*(tem - T_0) - 0.3918*(tem - T_0) + 87.663 ;
  double A    = 1398779.816/pow(epsi*tem,1.5)/b0 ;

  /* depend du sel */
  double b_cat ;
  double b_ani ;
  double S_cat ;
  double S_ani ;
  double c_ani ;
  double c_cat ;
  double z_ani ;
  double z_cat ;

  switch(SALT) {
  case(NaCl) : {
    b_cat = b_na_nacl ;
    b_ani = b_cl_nacl ;
    S_cat = S_na_nacl ;
    S_ani = S_cl_nacl ;
    c_ani = c_s ;
    c_cat = c_s ;
    z_ani = z_cl ;
    z_cat = z_na ;
    break ;
  }
  case(CaCl2) : {
    b_cat = b_ca_cacl2 ;
    b_ani = b_cl_cacl2 ;
    S_cat = S_ca_cacl2 ;
    S_ani = S_cl_cacl2 ;
    c_ani = 2*c_s ;
    c_cat = c_s ;
    z_ani = z_cl ;
    z_cat = z_ca ;
    break ;
  }
  default : {
    arret("non prevu") ;
  }
  }
  
  {
  /* concentrations */
  double c_h2o = (1 - c_s*V_salt)/V_h2o ;
  /* molalites*M_h2o */
  double m_ani  = c_ani/c_h2o ;
  double m_cat  = c_cat/c_h2o ;

  /* ion strength */
  double I     =  0.5*(z_ani*z_ani*m_ani + z_cat*z_cat*m_cat);
  
  double II_ani   = lna_i(tem,I,z_ani,b_ani,S_ani,A) ;
  double II_cat   = lna_i(tem,I,z_cat,b_cat,S_cat,A) ;

  /* activity of water */
  double lna_h2o = m_ani*II_ani + m_cat*II_cat ;

  return(lna_h2o) ;

  /* linearised term */
  lna_h2o = - (m_ani + m_cat) ;

  return(lna_h2o) ;
  }
}


double tortuosity(double p,double s)
{
  double tau_l_sat = 0.296e-3*exp(9.95*p)/p ;
  if(s > 0.) return(tau_l_sat/(s*(1 + 625*pow(1 - s,4)))) ;
  else return(0.) ;
}





double lng_LinLee(double T,double I,double z,double b,double S,double A)
/* Le log du coefficient d'activite d'un ion d'apres Lin & Lee */ 
{
  double alpha = 1.29,II = sqrt(I) ;
  double lng ;
  
  lng = - A*(II/(1 + b*II) + 2*log(1 + b*II)/b) + S*pow(I,alpha)/T ;
  
  return(lng*z*z) ;
}

double lng_TQN(double T,double I,double z,double b,double S,double A,double lna_w,double m_t)
/* Le log du coefficient d'activite d'un ion (T.Q Nguyen) :
   lng_i = dGamma/dm_i = (dGamma/dm_i)_I - 0.5*z_i*z_i*(lna_w + m_t)/I 
   lna_w = - m_t - sum_i ( m_i*lng_i ) + Gamma */
{
  double alpha = 1.29,II = sqrt(I) ;
  double lng ;
  
  lng = - A*2*log(1 + b*II)/b + S*pow(I,alpha)/(1+alpha)/T - 0.5*(lna_w + m_t)/I ;
  
  return(lng*z*z) ;
}

double lna_i(double T,double I,double z,double b,double S,double A)
/* Contribution de chaque ion au log de l'activite du solvant 
   lna_w = sum_i ( m_i*lna_i ) (T.Q Nguyen) */ 
{
  double alpha = 1.29,a1 = alpha/(1+alpha),II = sqrt(I) ;
  double lna ;
  
  lna = A*II/(1 + b*II) - a1*S*pow(I,alpha)/T ;
  
  return(-1 + lna*z*z) ;
}
