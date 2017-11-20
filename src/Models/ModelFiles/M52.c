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
#define NVI      (21)
#define NVE      (4)
#define NV0      (0)

/* Indices of equations */
#define E_h2o    (0)
#define E_salt   (1)
#define E_the    (2)


/* Indices of unknowns */
#define I_p_l    (0)
#define I_c_s    (1)
#define I_tem    (2)


/* Value of the nodal unknown (u and el must be used below) */
#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWNn(n,i)    (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])


/* Names of nodal unknowns */
#define P_l(n)   (UNKNOWN(n,I_p_l))
#define C_s(n)   (UNKNOWN(n,I_c_s))
#define TEM(n)   (UNKNOWN(n,I_tem))

#define P_ln(n)   (UNKNOWNn(n,I_p_l))
#define C_sn(n)   (UNKNOWNn(n,I_c_s))
#define TEM_n(n)  (UNKNOWNn(n,I_tem))


/* Names used for implicit terms */
#define M_H2O_l(n)  (f[(n)])
#define M_H2O_i(n)  (f[(n+2)])
#define M_SALT(n)   (f[(n+4)])
#define S(n)        (f[(n+6)])
#define S_H2O_l(n)  (f[(n+8)])
#define S_H2O_i(n)  (f[(n+10)])
#define S_SALT(n)   (f[(n+12)])
#define Rho(n)      (f[(n+14)])
#define Eps(n)      (f[(n+16)])
#define W_H2O       (f[(18)])
#define W_SALT      (f[(19)])
#define Q           (f[(20)])

#define M_H2O_ln(n)  (f_n[(n)])
#define M_H2O_in(n)  (f_n[(n+2)])
#define M_SALTn(n)   (f_n[(n+4)])
#define S_n(n)       (f_n[(n+6)])
#define S_H2O_ln(n)  (f_n[(n+8)])
#define S_H2O_in(n)  (f_n[(n+10)])
#define S_SALTn(n)   (f_n[(n+12)])


/* Names used for explicit terms */
#define KD_H2O      (va[(0)])
#define KD_salt     (va[(1)])
#define KF_salt     (va[(2)])
#define KTH         (va[(3)])


/* valences des ions */
#define z_na       (1.)
#define z_cl       (-1.)
#define z_ca       (2.)

/* volumes molaires liquides (m3/mole) */
#define V_h2o      (18.e-6)
#define V_na       (22.47e-6)
#define V_cl       (-0.35e-6)
#define V_ca       (-18.7e-6)

/* volumes molaires solides (m3/mole) */
#define V_ice      (19.63e-6)
#define V_nacl     (24.5e-6)
#define V_cacl2    (40.e-6)

/* Masses molaires (kg/m3) */
#define M_h2o      (18.e-3)
#define M_na       (23.e-3)
#define M_cl       (35.45e-3)
#define M_ca       (40.1e-3)
#define M_cacl2    (110.99e-3)
#define M_nacl     (58.45e-3)

/* coefficients de diffusion moleculaire (m2/s) */
#define d_ca       (7.92e-10)
#define d_na       (1.33e-9)
#define d_cl       (2.032e-9)
#define d_nacl     (1.e-9)
#define d_cacl2    (1.e-9)

/* viscosite */
#define mu_l       (1.79e-3)  /* Viscosite de l'eau (Pa.s) */

/* grandeurs de reference */
#define p_m        (1.e5)      /* Pression atmospherique (Pa)*/
#define T_m        (273.)      /* Temperature de fusion de la glace (K) */

/* Chaleurs specifiques (J/kg/K) */
#define C_l        (4180.)     /* Chaleur specifique du liquide */
#define C_i        (2000.)     /* Chaleur specifique de la glace */

/* entropie de fusion (J/mol/K) */
#define S_m        (23.54)

/* Conductivites thermiques (W/m/K) */
#define LAM_l      (0.6)
#define LAM_i      (2.2)

/* modules de compression (Pa) */
#define K_l        (1.8e9)
#define K_i        (7.8e9)

/* coefficients de dilatation thermique volumique (1/K) */
/* #define ALPHA_l    (-298.e-6) */
#define ALPHA_l(theta)    (-68.7e-6 + 13.877e-6*(theta))  /* coef secant D'apres Speedy (1987) */
#define ALPHA_i    (155.e-6)

/* constantes physiques */
#define R_g        (8.315)    /* Gaz parfait (J/mol/K) */


/* Type de sel */
#define NaCl   0
#define CaCl2  1

/* Choix du sel */
#define SALT   CaCl2

#if SALT == CaCl2
#define V_salt   V_cacl2
#define M_salt   M_cacl2
#define d_salt   d_cacl2
#elif SALT == NaCl
#define V_salt   V_nacl
#define M_salt   M_nacl
#define d_salt   d_nacl
#else
#define V_salt   V_nacl
#define M_salt   M_nacl
#define d_salt   d_nacl
#error "Type de sel non prevu"
#endif


/* Material Properties */
#define SaturationDegree(x)                Curve_ComputeValue(Element_GetCurve(el),x)
#define dSaturationDegree(x)               Curve_ComputeDerivative(Element_GetCurve(el),x)
#define RelativePermeabilityToLiquid(x)    Curve_ComputeValue(Element_GetCurve(el) + 1,x)


/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])

/* Fonctions */
static int    pm(const char *s) ;
static void   ComputeTransferCoefficients(Element_t*,double**,double*) ;
static int    flux(Element_t*,double**,double*) ;
static int    TangentCoefficients(Element_t*,double,double*) ;
static double activity(double,double) ;
static double tortuosity(double,double) ;
extern double lna_i(double,double,double,double,double,double) ;

/* Parametres */
static double lam_l = LAM_l,lam_i = LAM_i ;
static double rho_h2o_i0 = M_h2o/V_ice,rho_h2o_l0 = M_h2o/V_h2o ;
/* static double alpha_i = ALPHA_i,alpha_l = ALPHA_l ; */
static double alpha_i = ALPHA_i ;
static double v_salt = (V_salt/M_salt - V_h2o/M_h2o) ;
static double phi,k_int,lam_s,C_s,k_s,g_s,alpha_s ;


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
  Model_CopyNameOfUnknown(model,I_p_l,"p_l") ;
  Model_CopyNameOfUnknown(model,I_c_s,"c_s") ;
  Model_CopyNameOfUnknown(model,I_tem,"tem") ;
  
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


int ComputeInitialState(Element_t *el)
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 *  Return 0 if succeeded and -1 if failed
 */ 
{
  double *f  = Element_GetImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int    i ;
  
  
  /* We can skip if the element is a submanifold, 
   * e.g. a surface in 3D or a line in 2D */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*Donnees */
  phi     = GetProperty("phi") ;
  k_int   = GetProperty("k_int") ;
  lam_s   = GetProperty("lam_s") ;
  C_s     = GetProperty("C_s") ;
  k_s     = GetProperty("k_s") ;
  g_s     = GetProperty("g_s") ;
  alpha_s = GetProperty("alpha_s") ;
  
  /* masses of h2o, salt and entropies */
  for(i = 0 ; i < 2 ; i++) {
    double p_l = P_l(i) ;
    double tem = TEM(i) ;
    double c_s = C_s(i) ;
    
    double lna = activity(c_s,tem) ;
    
    /*pressure of ice crystallization */
    double p_i = p_m + (V_h2o*(p_l - p_m) + S_m*(T_m - tem) + R_g*tem*lna)/V_ice ;
    /* capirally pressure */
    double p_c = p_i - p_l ;
    
    /*saturation degree  */
    double s_l = SaturationDegree(p_c) ;
    double s_i = 1. - s_l ;

    /* mass densities */
    double alpha_l   = ALPHA_l(tem - T_m) ;
    double rho_h2o_i = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
    double rho_salt  = M_salt*c_s ;
    double rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - v_salt*rho_salt) ; 
    double rho_h2o_l = rho_l - rho_salt ;

    M_H2O_l(i) = phi*s_l*rho_h2o_l ;
    M_H2O_i(i) = phi*s_i*rho_h2o_i ;
    M_SALT(i)  = phi*s_l*rho_salt ;
    S_H2O_l(i) = C_l*log(tem/T_m) - alpha_l/rho_h2o_i0*(p_l - p_m) ;
    S_H2O_i(i) = C_i*log(tem/T_m) - S_m/M_h2o - alpha_i/rho_h2o_i0*(p_i - p_m)  ;
    S_SALT(i)  = S_H2O_l(i) ;
    S(i)       = C_s*log(tem/T_m)
               + M_H2O_l(i)*S_H2O_l(i)
               + M_H2O_i(i)*S_H2O_i(i)
               + M_SALT(i)*S_SALT(i) ;
  }

  {
    ComputeTransferCoefficients(el,u,f) ;
  }
  
  /* flux */
  {
    flux(el,u,f) ;
  }
  
  return(0) ;
}



int  ComputeExplicitTerms(Element_t *el,double t)
/* Thermes explicites (va)  */
{
  double *f = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int i ;
  
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Coefficients de transfert
  */
  ComputeTransferCoefficients(el,u,f) ;

  return(0) ;
}




int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int i ;
  

  
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Donnees
  */
  phi     = GetProperty("phi") ;
  k_int   = GetProperty("k_int") ;
  lam_s   = GetProperty("lam_s") ;
  C_s     = GetProperty("C_s") ;
  k_s     = GetProperty("k_s") ;
  g_s     = GetProperty("g_s") ;
  alpha_s = GetProperty("alpha_s") ;

  /* masses of h2o, salt and entropies */
  for(i = 0 ; i < 2 ; i++) {
    double p_l = P_l(i) ;
    double tem = TEM(i) ;
    double c_s = C_s(i) ;
    double lna = activity(c_s,tem) ;
    double p_i = p_m + (V_h2o*(p_l - p_m) + S_m*(T_m - tem) + R_g*tem*lna)/V_ice ;
    double p_c = p_i - p_l ;

    double s_l = SaturationDegree(p_c) ;
    double s_i = 1. - s_l ;

    /* mass densities */
    double alpha_l   = ALPHA_l(tem - T_m) ;
    double rho_h2o_i = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
    double rho_salt  = M_salt*c_s ;
    double rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - v_salt*rho_salt);
    double rho_h2o_l = rho_l - rho_salt ;

    M_H2O_l(i) = phi*s_l*rho_h2o_l ;
    M_H2O_i(i) = phi*s_i*rho_h2o_i ;
    M_SALT(i)  = phi*s_l*rho_salt ;
    S_H2O_l(i) = C_l*log(tem/T_m) - alpha_l/rho_h2o_i0*(p_l - p_m) ;
    S_H2O_i(i) = C_i*log(tem/T_m) - S_m/M_h2o - alpha_i/rho_h2o_i0*(p_i - p_m) ;
    S_SALT(i)  = S_H2O_l(i) ;
    S(i)       = C_s*log(tem/T_m)
                + M_H2O_l(i)*S_H2O_l(i)
                + M_H2O_i(i)*S_H2O_i(i)
                + M_SALT(i)*S_SALT(i) ;
  }

  /* flux */
  {
    double *f_n = Element_GetPreviousImplicitTerm(el) ;
    double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
    
    flux(el,u_n,f_n) ;
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
  double c[4*NEQ*NEQ] ;
  int    i ;
  
  /*
    Initialisation 
  */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  TangentCoefficients(el,dt,c) ;
  {
    double *km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }


  return(0) ;

#undef K
}



int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ + (i)])
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double *volume = FVM_ComputeCellVolumes(fvm) ;
  double surf ;
  int    i ;
  double zero = 0. ;
  
  /*
    INITIALISATION DU RESIDU
  */
  for(i = 0 ; i < NEQ*nn ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Boundary Surface Area */
  {
    double *surface = FVM_ComputeCellSurfaceAreas(fvm) ;
    surf = surface[1] ;
  }

  /*
    MASS BALANCE FOR H2O : (m_h2o1 - m_h2on) + dt * div(w_h2o1) = 0
  */
  R(0,E_h2o) -= volume[0]*(M_H2O_l(0) + M_H2O_i(0) - M_H2O_ln(0) - M_H2O_in(0)) + dt*surf*W_H2O ;
  R(1,E_h2o) -= volume[1]*(M_H2O_l(1) + M_H2O_i(1) - M_H2O_ln(1) - M_H2O_in(1)) - dt*surf*W_H2O ;
  /*
    MASS BALANCE FOR SALT : (m_salt1 - m_saltn) + dt * div(w_salt1) = 0
  */
  R(0,E_salt) -= volume[0]*(M_SALT(0) - M_SALTn(0)) + dt*surf*W_SALT ;
  R(1,E_salt) -= volume[1]*(M_SALT(1) - M_SALTn(1)) - dt*surf*W_SALT ;
  /*
    ENTROPY BALANCE : (S_1 - S_n) + dt * div(q_1/T_n+S_i*w_i1) = 0
  */
  R(0,E_the) -= volume[0]*(S(0) - S_n(0)) + dt*surf*Q ;
  R(1,E_the) -= volume[1]*(S(1) - S_n(1)) - dt*surf*Q ;
  
  return(0) ;

#undef R
}



int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
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
  phi     = GetProperty("phi") ;
  k_int   = GetProperty("k_int") ;
  lam_s   = GetProperty("lam_s") ;
  C_s     = GetProperty("C_s") ;
  k_s     = GetProperty("k_s") ;
  g_s     = GetProperty("g_s") ;
  alpha_s = GetProperty("alpha_s") ;
  
  /* Quantities */
  {
    int j = FVM_FindLocalCellIndex(fvm,s) ;

    /* pressure */
    double p_l  =  P_l(j) ;
    /* concentration */
    double c_s  =  C_s(j) ;
    /* temperature */
    double tem  =  TEM(j) ;
    
    /* activity */
    double lna  =  activity(c_s,tem) ;
    
    /* ice pressure */
    double p_i1 = V_h2o*(p_l - p_m)/V_ice ;
    double p_i2 = S_m*(T_m - tem)/V_ice ;
    double p_i3 = R_g*tem*lna/V_ice ;
    double p_i  = p_m + p_i1 + p_i2 + p_i3 ;
    
    /* capillary pressure */
    double p_c  = p_i - p_l ;
    
    /* liquid saturation */
    double s_l  = SaturationDegree(p_c) ;
    double s_i  = 1 - s_l ;
    
    /* mass densities */
    double alpha_l   = ALPHA_l(tem - T_m) ;
    double rho_h2o_i = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
    double rho_salt  = M_salt*c_s ;
    double rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - v_salt*rho_salt) ;
    
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
    Result_Store(r + i++,&p_l,"Liquid pressure",1) ;
    Result_Store(r + i++,&c_s,"Salt concentration",1) ;
    Result_Store(r + i++,&tem,"Temperature",1) ;
    Result_Store(r + i++,&s_l,"Liquid saturation",1) ;
    Result_Store(r + i++,&W_H2O,"Mass flow of H2O",1) ;
    Result_Store(r + i++,&W_SALT,"Mass flow of salt",1) ;
    Result_Store(r + i++,&p_i,"Ice pressure",1) ;
    Result_Store(r + i++,&rho_l,"Solution density",1) ;
    Result_Store(r + i++,&rho_h2o_i,"Ice density",1) ;
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



void  ComputeTransferCoefficients(Element_t *el,double **u,double *f)
{
  double *va = Element_GetExplicitTerm(el) ;
  int i ;

  /*
    Donnees
  */
  phi     = GetProperty("phi") ;
  k_int   = GetProperty("k_int") ;
  lam_s   = GetProperty("lam_s") ;
  C_s     = GetProperty("C_s") ;
  k_s     = GetProperty("k_s") ;
  g_s     = GetProperty("g_s") ;
  alpha_s = GetProperty("alpha_s") ;

  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;

  /* Transfer coefficients */
  for(i = 0 ; i < 2 ; i++) {
    double p_l = P_l(i) ;
    double tem = TEM(i) ;  
    double c_s = C_s(i) ;
    double lna = activity(c_s,tem) ;
    double p_i = p_m + (V_h2o*(p_l - p_m) + S_m*(T_m - tem) + R_g*tem*lna)/V_ice ;
    double p_c = p_i - p_l ;

    double s_l = SaturationDegree(p_c) ;
    double s_i = 1. - s_l ;
    /* mass densities */
    double alpha_l   = ALPHA_l(tem - T_m) ;
    double rho_salt  = M_salt*c_s ;
    double rho_l     = rho_h2o_l0*(1. + (p_l - p_m)/K_l - alpha_l*(tem - T_m) - v_salt*rho_salt) ; 
    double rho_h2o_l = rho_l - rho_salt ;

    double k_l   = k_int/mu_l*RelativePermeabilityToLiquid(p_c) ;
    double tau_l = tortuosity(phi,s_l) ;
    double lam_h2o = s_l*lam_l + s_i*lam_i ;

    KD_H2O  += rho_h2o_l*k_l ;
    KD_salt += rho_salt*k_l ;
    KF_salt += M_salt*phi*s_l*tau_l*d_salt ;
    KTH     += lam_s*(1 - 3*phi*(lam_s - lam_h2o)/(3*lam_s - (1 - phi)*(lam_s - lam_h2o))) ;
  }

  for(i = 0 ; i < NVE ; i++) va[i] *= 0.5 ;

}


int flux(Element_t *el,double **u_n,double *f_n)
{
  double *f = Element_GetImplicitTerm(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  int i ;
  

  if(Element_IsSubmanifold(el)) return(0) ;

  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx = x1 - x0 ;
    /* gradients */
    double grd_p_l = (P_l(1) - P_l(0))/dx ;
    double grd_c_s = (C_s(1) - C_s(0))/dx ;
    double grd_tem = (TEM(1) - TEM(0))/dx ;
    
    double tem     = (TEM_n(0)    + TEM_n(1)   )*0.5 ;
    double s_h2o_l = (S_H2O_ln(0) + S_H2O_ln(1))*0.5 ;
    double s_salt  = (S_SALTn(0)  + S_SALTn(1) )*0.5 ;
    
    /* flux */
    W_H2O  = - KD_H2O*grd_p_l  + KF_salt*grd_c_s ;
    W_SALT = - KD_salt*grd_p_l - KF_salt*grd_c_s ;
    Q      = - KTH/tem*grd_tem + s_h2o_l*W_H2O + s_salt*W_SALT ;
  }

  return(0) ;
} 






int  TangentCoefficients(Element_t *el,double dt,double *c)
/* Tangent matrix coefficients (c) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int    dec = NEQ*NEQ ;
  int    i ;
  
  
  /*
    Initialisation 
  */
  for(i = 0 ; i < ndof*ndof ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;


  /*
    Donnees
  */
  phi     = GetProperty("phi") ;
  k_int   = GetProperty("k_int") ;
  lam_s   = GetProperty("lam_s") ;
  C_s     = GetProperty("C_s") ;
  k_s     = GetProperty("k_s") ;
  g_s     = GetProperty("g_s") ;
  alpha_s = GetProperty("alpha_s") ;

  
  /* termes d'accumulation */
  for(i = 0 ; i < 2 ; i++) {
    /* Content terms at node i */
    double *cii = c + (i*nn + i)*NEQ*NEQ ;
    
    double p_l = P_l(i) ;
    double c_s = C_s(i) ;
    double tem = TEM(i) ;
    double lna = activity(c_s,tem) ;
    double p_i = p_m + (V_h2o*(p_l - p_m) + S_m*(T_m - tem) + R_g*tem*lna)/V_ice ;
    double p_c = p_i - p_l ;
    
    double s_l = SaturationDegree(p_c) ;
    double s_i = 1. - s_l ;

    /* mass densities */
    double alpha_l   = ALPHA_l(tem - T_m) ;
    double rho_h2o_i = rho_h2o_i0*(1. + (p_i - p_m)/K_i - alpha_i*(tem - T_m)) ;
    double rho_salt  = M_salt*c_s ;
    double rho_l     = rho_h2o_l0*(1. - v_salt*rho_salt + (p_l - p_m)/K_l - alpha_l*(tem - T_m)) ;
    double rho_h2o_l = rho_l - rho_salt ;

    /* pre-derivatives */
    double ds_lsdp_c = dSaturationDegree(p_c) ;
    
    /* derivatives with respect to ... */
    /* ... p_l */
    double dp_isdp_l = V_h2o/V_ice ;
    double dp_csdp_l = dp_isdp_l - 1. ;
    double drho_h2o_lsdp_l = rho_h2o_l0/K_l ;
    double drho_h2o_isdp_l = rho_h2o_i0/K_i*dp_isdp_l ;
    double dm_h2o_lsdp_l = phi*(ds_lsdp_c*dp_csdp_l*rho_h2o_l + s_l*drho_h2o_lsdp_l) ;
    double dm_h2o_isdp_l = phi*(-ds_lsdp_c*dp_csdp_l*rho_h2o_i + s_i*drho_h2o_isdp_l) ;
    double dm_saltsdp_l  = phi*ds_lsdp_c*dp_csdp_l*rho_salt ;
    double ds_h2o_lsdp_l = -alpha_l/rho_h2o_i0 ;
    double ds_h2o_isdp_l = -alpha_i/rho_h2o_i0*dp_isdp_l ;
    double ds_saltsdp_l  = ds_h2o_lsdp_l ;
    /* ... c_s */
    double dc_s = 1.e-5,c_s2 = c_s + dc_s ;
    double dlnasdc_s = (activity(c_s2,tem) - lna)/dc_s ;
    double dp_isdc_s = R_g*tem*dlnasdc_s/V_ice ;
    double dp_csdc_s = dp_isdc_s ;
    double drho_h2o_lsdc_s = -rho_h2o_l0*v_salt*M_salt - M_salt ;
    double drho_h2o_isdc_s = rho_h2o_i0/K_i*dp_isdc_s ;
    double dm_h2o_lsdc_s = phi*(ds_lsdp_c*dp_csdc_s*rho_h2o_l + s_l*drho_h2o_lsdc_s) ;
    double dm_h2o_isdc_s = phi*(-ds_lsdp_c*dp_csdc_s*rho_h2o_i + s_i*drho_h2o_isdc_s) ;
    double dm_saltsdc_s  = phi*(ds_lsdp_c*dp_csdc_s*rho_salt + s_l*M_salt) ;
    double ds_h2o_lsdc_s = 0. ;
    double ds_h2o_isdc_s = -alpha_i/rho_h2o_i0*dp_isdc_s ;
    double ds_saltsdc_s = ds_h2o_lsdc_s ;
    /* ... tem */
    double dtem = 1.,tem2 = tem + dtem ;
    double dlnasdtem = (activity(c_s,tem2) - lna)/dtem ;
    double dp_isdtem = (- S_m + R_g*lna + R_g*tem*dlnasdtem)/V_ice ;
    double dp_csdtem = dp_isdtem ;
    double drho_h2o_lsdtem = -rho_h2o_l0*alpha_l ;
    double drho_h2o_isdtem = rho_h2o_i0/K_i*dp_isdtem - rho_h2o_i0*alpha_i ;
    double dm_h2o_lsdtem = phi*(ds_lsdp_c*dp_csdtem*rho_h2o_l + s_l*drho_h2o_lsdtem) ;
    double dm_h2o_isdtem = phi*(-ds_lsdp_c*dp_csdtem*rho_h2o_i + s_i*drho_h2o_isdtem) ;
    double dm_saltsdtem  = phi*ds_lsdp_c*dp_csdtem*rho_salt ;
    double ds_h2o_lsdtem = C_l/tem ;
    double ds_h2o_isdtem = C_i/tem - alpha_i/rho_h2o_i0*dp_isdtem ;
    double ds_saltsdtem  = ds_h2o_lsdtem ;
    /*
      MASS BALANCE FOR H2O : (m_h2o1 - m_h2on) + dt * div(w_h2o1) = 0
    */
    cii[E_h2o*NEQ + I_p_l] = dm_h2o_lsdp_l + dm_h2o_isdp_l ;
    cii[E_h2o*NEQ + I_c_s] = dm_h2o_lsdc_s + dm_h2o_isdc_s ;
    cii[E_h2o*NEQ + I_tem] = dm_h2o_lsdtem + dm_h2o_isdtem ;
    /*
      MASS BALANCE FOR SALT : (m_salt1 - m_saltn) + dt * div(w_salt1) = 0
    */
    cii[E_salt*NEQ + I_p_l] = dm_saltsdp_l ;
    cii[E_salt*NEQ + I_c_s] = dm_saltsdc_s ;
    cii[E_salt*NEQ + I_tem] = dm_saltsdtem ;
    /*
      ENTROPY BALANCE : (S_1 - S_n) + dt * div(q_1/T_n+S_i*w_i1) = 0
    */
    cii[E_the*NEQ + I_p_l] = dm_h2o_lsdp_l*S_H2O_l(i) + dm_h2o_isdp_l*S_H2O_i(i) + dm_saltsdp_l*S_SALT(i) + M_H2O_l(i)*ds_h2o_lsdp_l + M_H2O_i(i)*ds_h2o_isdp_l + M_SALT(i)*ds_saltsdp_l ;
    cii[E_the*NEQ + I_c_s] = dm_h2o_lsdc_s*S_H2O_l(i) + dm_h2o_isdc_s*S_H2O_i(i) + dm_saltsdc_s*S_SALT(i) + M_H2O_l(i)*ds_h2o_lsdc_s + M_H2O_i(i)*ds_h2o_isdc_s + M_SALT(i)*ds_saltsdc_s ;
    cii[E_the*NEQ + I_tem] = C_s/tem + dm_h2o_lsdtem*S_H2O_l(i) + dm_h2o_isdtem*S_H2O_i(i) + dm_saltsdtem*S_SALT(i) + M_H2O_l(i)*ds_h2o_lsdtem + M_H2O_i(i)*ds_h2o_isdtem + M_SALT(i)*ds_saltsdtem ;
  }
    
    
    /*
      termes d'ecoulement
    */
  {
    /* Transfer terms from node i to node j: wij = flux from i to j
     * double *cij = c + (i*nn + j)*NEQ*NEQ 
     * cij is the derivative of wij wrt variables at node i: cij = wij,i.
     * If wij = - Kij*(Pj - Pi)/dij then cij = Kij/dij */
    double *c01 = c + NEQ*NEQ ;
    double *c10 = c + nn*NEQ*NEQ ;
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx = x1 - x0 ;
    double dtdx     = dt/dx ;
    double trd_h2o  = dtdx*KD_H2O ;
    double trf_h2o  = dtdx*(-KF_salt) ;
    double trd_salt = dtdx*KD_salt ;
    double trf_salt = dtdx*KF_salt ;
    double trth     = dtdx*KTH ;

    double tem      = (TEM_n(0)    + TEM_n(1)   )*0.5 ;
    double s_h2o_l  = (S_H2O_ln(0) + S_H2O_ln(1))*0.5 ;
    double s_salt   = (S_SALTn(0)  + S_SALTn(1 ))*0.5 ;
    double cc ;
    /*
      MASS BALANCE FOR H2O : (m_h2o1 - m_h2on) + dt * div(w_h2o1) = 0
    */
    c01[E_h2o*NEQ + I_p_l] = trd_h2o ;
    c10[E_h2o*NEQ + I_p_l] = trd_h2o ;
  
    c01[E_h2o*NEQ + I_c_s] = trf_h2o ;
    c10[E_h2o*NEQ + I_c_s] = trf_h2o ;

    /*
      MASS BALANCE FOR SALT : (m_salt1 - m_saltn) + dt * div(w_salt1) = 0
    */
    c01[E_salt*NEQ + I_p_l]   = trd_salt ;
    c10[E_salt*NEQ + I_p_l]   = trd_salt ;
  
    c01[E_salt*NEQ + I_c_s]   = trf_salt ;
    c10[E_salt*NEQ + I_c_s]   = trf_salt ;

    /*
      ENTROPY BALANCE : (S_1 - S_n) + dt * div(q_1/T_n+S_i*w_i1) = 0
    */
    cc = s_h2o_l*trd_h2o + s_salt*trd_salt ;
    c01[E_the*NEQ + I_p_l]    = cc ;
    c10[E_the*NEQ + I_p_l]    = cc ;

    cc = s_h2o_l*trf_h2o + s_salt*trf_salt ;
    c01[E_the*NEQ + I_c_s]    = cc ;
    c10[E_the*NEQ + I_c_s]    = cc ;

    cc = trth/tem ;
    c01[E_the*NEQ + I_tem]    = cc ;
    c10[E_the*NEQ + I_tem]    = cc ;
  }

  return(dec) ;
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
