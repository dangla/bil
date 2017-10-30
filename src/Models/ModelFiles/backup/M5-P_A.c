#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"
#include "FVM.h"


#define TITLE   "Drying-Wetting (Isothermal 1D)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"

/* Macros */
#define NEQ      (2)
#define NVI      (15)
#define NVE      (6)
#define NV0      (0)


/* Equation index */
#define E_Mass   (0)
#define E_Air    (1)

/* Unknown index*/
#define I_P_L    (0)
#define I_P_A    (1)


#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWNn(n,i)    (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])


/* We define some names for primary unknowns */
#define P_L(n)      (UNKNOWN(n,I_P_L))
#define P_A(n)      (UNKNOWN(n,I_P_A))

#define M_L(n)   (f[(n)])
#define M_V(n)   (f[(n+2)])
#define M_A(n)   (f[(n+4)])
#define M_G(n)   (f[(n+6)])
#define M_T(n)   (f[(n+8)])
#define W_L      (f[(10)])
#define W_V      (f[(11)])
#define W_A      (f[(12)])
#define W_G      (f[(13)])
#define W_T      (f[(14)])

#define M_Ln(n)  (f_n[(n)])
#define M_Vn(n)  (f_n[(n+2)])
#define M_An(n)  (f_n[(n+4)])
#define M_Gn(n)  (f_n[(n+6)])
#define M_Tn(n)  (f_n[(n+8)])

#define K_L      (va[(0)])
#define KD_V     (va[(1)])
#define KF_V     (va[(2)])
#define KD_A     (va[(3)])
#define KF_A     (va[(4)])
#define K_G      (va[(5)])

/* Masses molaires (kg) */
#define M_h2o    (18.e-3)
#define M_air    (28.8e-3)

/* Coefficient de diffusion moleculaire (m2/s) */
#define D_av0    (2.48e-5)

/* Constantes physiques */
#define RT       (2436.)     /* produit de R et T (J/mol) */
                             /* R = 8.3143 J/K/mol cste des gaz parfait */
                             /* T = 293 K */

/* Viscosites (Pa.s) */
#define mu_l     (1.e-3)
#define mu_g     (1.8e-5)

/* masse volumique (kg/m3) */
#define rho_l    (1000.)


/* Liquid gas thermodynamic equilibrium */
#define VaporPressure(p_l)  (p_v0*exp(M_h2o/RT*(p_l - p_l0)/rho_l))


/* Material properties */
#define SATURATION_CURVE           (Element_GetCurve(el))
#define RELATIVEPERMLIQ_CURVE      (Element_GetCurve(el) + 1)
#define RELATIVEPERMGAS_CURVE      (Element_GetCurve(el) + 2)

#define Saturation(x)    (saturation(x,p_c3,SATURATION_CURVE))
#define dSaturation(x)   (dsaturation(x,p_c3,SATURATION_CURVE))
#define RelativePermeabilityToLiquid(x)  (Curve_ComputeValue(RELATIVEPERMLIQ_CURVE,x))
#define RelativePermeabilityToGas(x)  (Curve_ComputeValue(RELATIVEPERMGAS_CURVE,x))
#define TortuosityToGas(f,sg)     (pow(f,aa)*pow(sg,bb))
#define aa                        (1./3)  /* 1/3 Millington, 2.74 */
#define bb                        (7./3)  /* 7/3 Millington,  4.2 */


/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])


/* Functions */
static int     pm(char *s) ;
static void    GetProperties(Element_t*) ;

static double* ComputeComponents(Element_t*,double**,double*,double,int) ;
static void    ComputeSecondaryComponents(Element_t*,double,double*) ;
static double* ComputeComponentDerivatives(Element_t*,double,double*,double,int) ;

static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static void    ComputeFluxes(Element_t*,double**) ;
static double* Fluxes(Element_t*,double*) ;
static int     TangentCoefficients(Element_t*,double,double*) ;

static double  saturation(double,double,Curve_t*) ;
static double  dsaturation(double,double,Curve_t*) ;


/* Parameters */
static double gravite,phi,k_int,p_l0,p_v0,p_a0,p_g0,p_c3 ;


#define NbOfComponents    (NEQ+15)
static double Components[NbOfComponents] ;
static double dComponents[NbOfComponents] ;

#define I_M_L          (NEQ)
#define I_M_V          (NEQ+1)
#define I_M_A          (NEQ+2)
#define I_M_G          (NEQ+3)
#define I_M_T          (NEQ+4)

#define I_RHO_V        (NEQ+5)
#define I_RHO_A        (NEQ+6)
#define I_RHO_G        (NEQ+7)

#define I_S_L          (NEQ+8)

#define I_H_L          (NEQ+9)
#define I_H_G          (NEQ+10)
#define I_COOR_X       (NEQ+11)

#define I_P_V          (NEQ+12)
#define I_P_G          (NEQ+13)
#define I_C_V          (NEQ+14)
  
  

#define NbOfComponentFluxes    (5)
static double ComponentFluxes[NbOfComponentFluxes] ;

#define I_W_L           (0)
#define I_W_V           (1)
#define I_W_A           (2)
#define I_W_G           (3)
#define I_W_T           (4)



int pm(char *s)
{
  if(!strcmp(s,"gravite")) return (0) ;
  else if(!strcmp(s,"phi")) return (1) ;
  else if(!strcmp(s,"k_int")) return (2) ;
  else if(!strcmp(s,"p_l0")) return (3) ;
  else if(!strcmp(s,"p_v0")) return (4) ;
  else if(!strcmp(s,"p_a0")) return (5) ;
  else if(!strcmp(s,"p_c3")) return (6) ;
  else return(-1) ;
}


void GetProperties(Element_t *el)
{
  gravite = GetProperty("gravite") ;
  phi     = GetProperty("phi") ;
  k_int   = GetProperty("k_int") ;
  p_c3    = GetProperty("p_c3") ;
  p_l0    = GetProperty("p_l0") ;
  p_v0    = GetProperty("p_v0") ;
  p_a0    = GetProperty("p_a0") ;
  p_g0    = p_a0 + p_v0 ;
}


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_Mass,"mass") ;
  Model_CopyNameOfEquation(model,E_Air,"air") ;

  Model_CopyNameOfUnknown(model,I_P_L,"p_l") ;
  Model_CopyNameOfUnknown(model,I_P_A,"p_a") ;
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  NbOfProp = 7 ;

  Material_ScanProperties(mat,datafile,pm) ;
  
  return(NbOfProp) ;
}


int PrintModelChar(Model_t *model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n") ;
  printf("The equation to be solved is:\n") ;
  printf("\t- Mass balance of water  (eau)\n") ;
  printf("\t- Mass balance of air    (gas)\n") ;
  
  printf("\n") ;
  printf("The primary unknown is:\n") ;
  printf("\t- Liquid pressure        (p_l)\n") ;
  printf("\t- Air pressure           (p_a)\n") ;

  printf("\n") ;
  printf("Some other informations\n") ;
  printf("Example of input data\n") ;
  printf("\n") ;

  fprintf(ficd,"gravite = -9.81    # Gravity\n") ;
  fprintf(ficd,"phi = 0.3          # Porosity\n") ;
  fprintf(ficd,"k_int = 1.e-20     # Intrinsic Permeability\n") ;
  fprintf(ficd,"p_c3 = 1.e6        # Capillary pressure parameter\n") ;
  fprintf(ficd,"Curves = sol       # File name: p_c S_l k_rl k_rg\n") ;
  fprintf(ficd,"p_l0 = 100000      # Pression de liquide de reference\n") ;
  fprintf(ficd,"p_v0 = 2460        # Pression de vapeur de reference\n") ;
  fprintf(ficd,"p_a0 = 97540       # Pression d\'air de reference\n") ;
  fprintf(ficd,"Curves = my_file   # File name: p_c S_l k_rl k_rg\n") ;
  
  return(NEQ) ;
}


int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
{
  
  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
  
  return(0) ;
}


int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/** Compute the residu (r) due to loads */
{
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  
  {
    double *r1 = FVM_ComputeSurfaceLoadResidu(fvm,cg,t,dt) ;
    int i ;
    
    for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
  }
  
  return(0) ;
}


int ComputeInitialState(Element_t *el)
/* Initialise les variables du systeme (f,va) */ 
{
  double*  f  = Element_GetImplicitTerm(el) ;
  double** u  = Element_ComputePointerToNodalUnknowns(el) ;
  int    i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;

  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
  /* Pre-initialization */
  for(i = 0 ; i < 2 ; i++) {
    double p_l  = P_L(i) ;
    double p_v  = VaporPressure(p_l) ;
    double p_a  = (P_A(i) > 0) ? P_A(i) : (-P_A(i)) - p_v ;
    
    P_A(i)  = p_a ;
  }
  
  
  /* Fluid mass contents */
  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,f,0,i) ;

    /* Back up */
    M_L(i) = x[I_M_L] ;
    M_V(i) = x[I_M_V] ;
    M_A(i) = x[I_M_A] ;
    M_T(i) = x[I_M_T] ;
  }

  /* Transfer coefficients */
  ComputeTransferCoefficients(el,u,f) ;

  /* Fluxes */
  ComputeFluxes(el,u) ;
  
  return(0) ;
}


int  ComputeExplicitTerms(Element_t *el,double t)
/* Explicites terms  */
{
  double*  f = Element_GetPreviousImplicitTerm(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  if(Element_IsSubmanifold(el)) return(0) ;

  
  /*
    Input data
  */
  GetProperties(el) ;
  
  /* Transfer coefficients */
  ComputeTransferCoefficients(el,u,f) ;
  
  return(0) ;
}


int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/* Implicit terms */
{
  double* f  = Element_GetCurrentImplicitTerm(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int    i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;

  
  /*
    Input data
  */
  GetProperties(el) ;
  

  /* masses d'eau et d'air */
  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,f,0,i) ;

    /* Back up */
    M_L(i) = x[I_M_L] ;
    M_V(i) = x[I_M_V] ;
    M_A(i) = x[I_M_A] ;
    M_T(i) = x[I_M_T] ;
    
    {
      double p_a     = x[I_P_A] ;
    
      if(p_a < 0.) {
        printf("pression d\'air = %e\n",p_a) ;
        return(1) ;
      }
    }
  }
  

  /* Fluxes */
  ComputeFluxes(el,u) ;
  
  return(0) ;
}


int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/* Matrix (k) */
{
#define K(i,j)    (k[(i)*2*NEQ + (j)])
  int    nn = Element_GetNbOfNodes(el) ;
  int    ndof = nn*NEQ ;
  double c[Element_MaxNbOfNodes*Element_MaxNbOfNodes*NEQ*NEQ] ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;
  
  /* initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  TangentCoefficients(el,dt,c) ;
  {
    double *km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }
  
  return(0) ;

#undef K
}


#ifdef NOTDEFINED
int  ComputeMatrix1(Element_t *el,double t,double dt,double *k)
/* Matrix (k) */
{
#define K(i,j)    (k[(i)*2*NEQ + (j)])
  double* va = Element_GetExplicitTerm(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int    nn = Element_GetNbOfNodes(el) ;
  int    ndof = nn*NEQ ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double *volume = FVM_ComputeCellVolumes(fvm) ;
  /* Numeros d'ordre des equations et des inconnues locales */
#define NEQ1     (NEQ+1)
#define E_vap    NEQ
#define I_p_v    NEQ
#define E_liq    E_eau

#define KE(i,j)  (ke[(i)*2*NEQ1+(j)])
  double p_c,s_l,s_g,p_v,p_l,p_a,p_g,dslsdpc,pv[2],pa[2],pg[2],dpvsdpl[2] ;
  double rho_v,rho_a,rho_g ;
  double cv[2],ca[2] ;
  double tr_l,trd_v,trf_v,trd_a,trf_a ;
  double dx ;
  double surf ;
  int    i,j,n ;
  double zero = 0.,un = 1. ;
  double ke[4*NEQ1*NEQ1] ;
  
  /* initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /*
    Input data
  */
  GetProperties(el) ;
  

  /*
    INITIALISATION DE LA MATRICE
  */
  for(i=0;i<4*NEQ1*NEQ1;i++) ke[i] = zero ;
  
  /*
    CALCUL DE volume ET DE surf
  */
  
  /* Boundary Surface Area */
  {
    double *area = FVM_ComputeCellSurfaceAreas(fvm) ;
    surf = area[1] ;
  }
  
  /* termes d'accumulation */
  for(i=0;i<2;i++) {
    p_l   = P_L(i) ;
    p_a   = P_A(i) ;
    p_v   = VaporPressure(p_l) ;
    p_g   = p_v + p_a ;
    p_c   = p_g - p_l ;

    pv[i] = p_v ;
    pa[i] = p_a ;
    pg[i] = p_g ;

    s_l   = Saturation(p_c) ;
    s_g   = un - s_l ;

    rho_v  = p_v*M_h2o/RT ;
    rho_a  = p_a*M_air/RT ;
    rho_g  = rho_v + rho_a ;

    cv[i]  = rho_v/rho_g ;
    ca[i]  = un - cv[i] ;

    dslsdpc = dSaturation(p_c) ;
    dpvsdpl[i] = rho_v/rho_l ;
    /*
      CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l) = +A
    */
    KE(E_liq+i*NEQ1,I_P_L+i*NEQ1) += volume[i]*phi*rho_l*(-dslsdpc) ;
    KE(E_liq+i*NEQ1,I_p_v+i*NEQ1) += volume[i]*phi*rho_l*dslsdpc ;
    KE(E_liq+i*NEQ1,I_P_A+i*NEQ1) += volume[i]*phi*rho_l*dslsdpc ;
    /*
      CONSERVATION DE L'EAU VAPEUR  : (m_v1 - m_vn) + dt * div(w_v) = -A
    */
    KE(E_vap+i*NEQ1,I_P_L+i*NEQ1) += volume[i]*phi*rho_v*dslsdpc ;
    KE(E_vap+i*NEQ1,I_p_v+i*NEQ1) += volume[i]*phi*(M_h2o/RT*s_g - rho_v*dslsdpc) ;
    KE(E_vap+i*NEQ1,I_P_A+i*NEQ1) += volume[i]*phi*rho_v*(-dslsdpc) ;
    /*
      CONSERVATION DE L'AIR SEC     : (m_a1 - m_an) + dt * div(w_a) = 0
    */
    KE(E_Air+i*NEQ1,I_P_L+i*NEQ1) += volume[i]*phi*rho_a*dslsdpc ;
    KE(E_Air+i*NEQ1,I_p_v+i*NEQ1) += volume[i]*phi*rho_a*(-dslsdpc) ;
    KE(E_Air+i*NEQ1,I_P_A+i*NEQ1) += volume[i]*phi*(M_air/RT*s_g - rho_a*dslsdpc) ;
  }
  /*
    termes d'ecoulement
  */
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    dx = x1 - x0 ;
  }
  tr_l  = dt*surf/dx*K_L ;
  trd_v = dt*surf/dx*KD_V ;
  trd_a = dt*surf/dx*KD_A ;
  trf_v = dt*surf/dx*KF_V ;
  trf_a = dt*surf/dx*KF_A ;
  /*
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l) = +A
  */
  KE(E_liq,I_P_L)           += + tr_l ;
  KE(E_liq,I_P_L+NEQ1)      += - tr_l ;
  KE(E_liq+NEQ1,I_P_L)      += - tr_l ;
  KE(E_liq+NEQ1,I_P_L+NEQ1) += + tr_l ;
  /*
    CONSERVATION DE L'EAU VAPEUR  : (m_v1 - m_vn) + dt * div(w_v) = -A
  */
  KE(E_vap,I_p_v)      += + trd_v + trf_v*ca[0]*cv[0]/pv[0] ;
  KE(E_vap,I_p_v+NEQ1) += - trd_v - trf_v*ca[1]*cv[1]/pv[1] ;
  KE(E_vap,I_P_A)      += + trd_v + trf_v*(-ca[0]*cv[0]/pa[0]) ;
  KE(E_vap,I_P_A+NEQ1) += - trd_v - trf_v*(-ca[1]*cv[1]/pa[1]) ;
  
  KE(E_vap+NEQ1,I_p_v+NEQ1) += + trd_v + trf_v*ca[1]*cv[1]/pv[1] ;
  KE(E_vap+NEQ1,I_p_v)      += - trd_v - trf_v*ca[0]*cv[0]/pv[0] ;
  KE(E_vap+NEQ1,I_P_A+NEQ1) += + trd_v + trf_v*(-ca[1]*cv[1]/pa[1]) ;
  KE(E_vap+NEQ1,I_P_A)      += - trd_v - trf_v*(-ca[0]*cv[0]/pa[0]) ;
  /*
    CONSERVATION DE L'AIR SEC     : (m_a1 - m_an) + dt * div(w_a) = 0
  */
  KE(E_Air,I_P_A)      += + trd_a + trf_a*ca[0]*cv[0]/pa[0] ;
  KE(E_Air,I_P_A+NEQ1) += - trd_a - trf_a*ca[1]*cv[1]/pa[1] ;
  KE(E_Air,I_p_v)      += + trd_a + trf_a*(-ca[0]*cv[0]/pv[0]) ;
  KE(E_Air,I_p_v+NEQ1) += - trd_a - trf_a*(-ca[1]*cv[1]/pv[1]) ;
  
  KE(E_Air+NEQ1,I_P_A+NEQ1) += + trd_a + trf_a*ca[1]*cv[1]/pa[1] ;
  KE(E_Air+NEQ1,I_P_A)      += - trd_a - trf_a*ca[0]*cv[0]/pa[0] ;
  KE(E_Air+NEQ1,I_p_v+NEQ1) += + trd_a + trf_a*(-ca[1]*cv[1]/pv[1]) ;
  KE(E_Air+NEQ1,I_p_v)      += - trd_a - trf_a*(-ca[0]*cv[0]/pv[0]) ;

  /*
    ELIMINATION DE LA COLONNE I_p_v :
    COLONNE I_P_L += COLONNE I_p_v * dpvsdpl
  */
  for(j=0;j<NEQ1;j++) {
    for(i=0;i<2;i++) {
      KE(j+i*NEQ1,I_P_L+i*NEQ1) += KE(j+i*NEQ1,I_p_v+i*NEQ1) * dpvsdpl[i] ;
    }
    KE(j,I_P_L+NEQ1) += KE(j,I_p_v+NEQ1) * dpvsdpl[1] ;
    KE(j+NEQ1,I_P_L) += KE(j+NEQ1,I_p_v) * dpvsdpl[0] ;
  }
  /*
    ELIMINATION DE LA LIGNE E_vap :
    LIGNE E_liq += LIGNE E_vap
  */
  for(j=0;j<NEQ;j++) {
    for(i=0;i<2;i++) KE(E_liq+i*NEQ1,j+i*NEQ1) += KE(E_vap+i*NEQ1,j+i*NEQ1) ;
    KE(E_liq,j+NEQ1) += KE(E_vap,j+NEQ1) ;
    KE(E_liq+NEQ1,j) += KE(E_vap+NEQ1,j) ;
  }

  /*
    ASSEMBLAGE DANS K
  */
  for(i=0;i<NEQ;i++) {
    for(j=0;j<NEQ;j++) {
      for(n=0;n<2;n++)  K(i+n*NEQ,j+n*NEQ) += KE(i+n*NEQ1,j+n*NEQ1) ;
      K(i,j+NEQ) += KE(i,j+NEQ1) ;
      K(i+NEQ,j) += KE(i+NEQ1,j) ;
    }
  }

  return(0) ;

#undef NEQ1
#undef E_vap
#undef I_p_v
#undef E_liq
#undef K
#undef KE
}
#endif


int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int    ndof = nn*NEQ ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double *volume = FVM_ComputeCellVolumes(fvm) ;
  double surf ;
  int    i ;
  double zero = 0. ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Boundary Surface Area */
  {
    double *area = FVM_ComputeCellSurfaceAreas(fvm) ;
    surf = area[1] ;
  }
  
  /*
    Conservation of total mass
  */
  R(0,E_Mass)  -= volume[0]*(M_T(0)  - M_Tn(0))  + dt*surf*W_T ;
  R(1,E_Mass)  -= volume[1]*(M_T(1)  - M_Tn(1))  - dt*surf*W_T ;
  
  /*
    Conservation of air mass
  */
  R(0,E_Air)  -= volume[0]*(M_A(0)  - M_An(0))  + dt*surf*W_A ;
  R(1,E_Air)  -= volume[1]*(M_A(1)  - M_An(1))  - dt*surf*W_A ;
  
  return(0) ;

#undef R
}


#ifdef NOTDEFINED
int  ComputeResidu1(Element_t *el,double t,double dt,double *r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int    ndof = nn*NEQ ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double *volume = FVM_ComputeCellVolumes(fvm) ;
  /* Numeros d'ordre des equations et des inconnues locales */
#define NEQ1     (NEQ+1)
#define E_vap    NEQ
#define E_liq    E_eau
#define R(n,i)    (r[(n)*NEQ+(i)])
#define RE(n,i)   (re[(n)*NEQ1+(i)])
  double surf ;
  int    i,j ;
  double zero = 0. ;
  double re[NEQ1*2] ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    INITIALISATION DU RESIDU
  */
  for(i=0;i<NEQ1*2;i++) re[i] = zero ;
  
  /* Boundary Surface Area */
  {
    double *area = FVM_ComputeCellSurfaceAreas(fvm) ;
    surf = area[1] ;
  }
  
  /*
    CONSERVATION DE L'EAU LIQUIDE : (m_l1 - m_ln) + dt * div(w_l) = +A
  */
  RE(0,E_liq) -= volume[0]*(M_L(0) - M_Ln(0)) + dt*surf*W_L ;
  RE(1,E_liq) -= volume[1]*(M_L(1) - M_Ln(1)) - dt*surf*W_L ;
  /*
    CONSERVATION DE L'EAU VAPEUR  : (m_v1 - m_vn) + dt * div(w_v) = -A
  */
  RE(0,E_vap) -= volume[0]*(M_V(0) - M_Vn(0)) + dt*surf*W_V ;
  RE(1,E_vap) -= volume[1]*(M_V(1) - M_Vn(1)) - dt*surf*W_V ;
  /*
    CONSERVATION DE L'AIR SEC     : (m_a1 - m_an) + dt * div(w_a) = 0
  */
  RE(0,E_Air) -= volume[0]*(M_A(0) - M_An(0)) + dt*surf*W_A ;
  RE(1,E_Air) -= volume[1]*(M_A(1) - M_An(1)) - dt*surf*W_A ;
  /*
    ELIMINATION DE LA LIGNE E_vap : LIGNE E_liq += LIGNE E_vap
  */
  for(i=0;i<2;i++) RE(i,E_liq) += RE(i,E_vap) ;
  /*
    ASSEMBLAGE DANS R
  */
  for(j=0;j<NEQ;j++) for(i=0;i<2;i++) R(i,j) += RE(i,j) ;
  
  return(0) ;

#undef NEQ1
#undef E_vap
#undef E_liq
#undef R
#undef RE
}
#endif




int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/* Les valeurs exploitees (s) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int    nso = 8 ;
  int    i ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /*
    Input data
  */
  GetProperties(el) ;
  

  /* initialisation */
  for(i = 0 ; i < nso ; i++) {
    Result_SetValuesToZero(r+i) ;
  }

  /* quantites exploitees */
  {
    int    j = FVM_FindLocalCellIndex(fvm,s) ;
    /* Components */
    double *x = ComputeComponents(el,u,f,0,j) ;
    double p_l     = x[I_P_L] ;
    double p_a     = x[I_P_A] ;
    double p_v     = x[I_P_V] ;
    double p_g     = x[I_P_G] ;
    double p_c     = p_g - p_l ;
    double s_l     = Saturation(p_c) ;
    /* flux */
    double w_l  = W_L ;
    double w_v  = W_V ;
    double w_a  = W_A ;
    double w_g  = W_G ;
    
    i = 0 ;
    Result_Store(r + i++,&p_l,"p_l",1) ;
    Result_Store(r + i++,&p_a,"p_a",1) ;
    Result_Store(r + i++,&p_v,"p_v",1) ;
    Result_Store(r + i++,&w_l,"w_l",1) ;
    Result_Store(r + i++,&w_a,"w_a",1) ;
    Result_Store(r + i++,&w_v,"w_v",1) ;
    Result_Store(r + i++,&w_g,"w_g",1) ;
    Result_Store(r + i++,&s_l,"saturation",1) ;
  }

  return(nso) ;
}


void ComputeTransferCoefficients(Element_t *el,double **u,double *f)
/* Termes explicites (va)  */
{
  double *va = Element_GetExplicitTerm(el) ;
  int i ;
  
  /* initialisation */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;
  
  
  /* Transfer coefficients */
  {
    double *x = Components ;
  
    /* Primary Variables */
    x[I_P_L ] = (P_L(0) + P_L(1))*0.5 ;
    x[I_P_A ] = (P_A(0) + P_A(1))*0.5 ;
  
    /* Needed variables to compute secondary components */
    x[I_COOR_X] = 0. ;
    
    ComputeSecondaryComponents(el,0,x) ;
    
    double p_l    = x[I_P_L] ;
    double p_g    = x[I_P_G] ;
    double p_c    = p_g - p_l ;

    double s_l    = x[I_S_L] ;
    double s_g    = 1 - s_l ;

    double rho_v  = x[I_RHO_V] ;
    double rho_a  = x[I_RHO_A] ;
    double rho_g  = x[I_RHO_G] ;

    double c_v    = rho_v/rho_g ;
    double c_a    = 1 - c_v ;
  
    double k_rl  = RelativePermeabilityToLiquid(p_c) ;
    double k_rg  = RelativePermeabilityToGas(p_c) ;
  
    double kh_l  = k_int/mu_l*k_rl ;
    double kh_g  = k_int/mu_g*k_rg ;
  
    double tau   = TortuosityToGas(phi,s_g) ;
    double D_av  = D_av0*p_g0/p_g ;
    double D_eff = phi*s_g*tau*D_av ;
  
    /* Darcy */
    K_L   += rho_l*kh_l ;           /* liquid */
    K_G   += rho_g*kh_g ;           /* gas */
    KD_V  += rho_v*kh_g ;           /* vapor */
    KD_A  += rho_a*kh_g ;           /* air */
    
    /* Fick */
    KF_V  += rho_g*D_eff ;          /* vapor */
    KF_A  += rho_g*D_eff ;          /* air */
    
    /* Barometric diffusion effects */
    KD_V += D_eff*c_v*c_a*(M_air - M_h2o)/RT ; /* vapor */
    KD_A += D_eff*c_v*c_a*(M_h2o - M_air)/RT ; /* air */
  }
  
  /* Averaging */
  /* for(i = 0 ; i < NVE ; i++) va[i] *= 0.5 ; */
}


void ComputeTransferCoefficients1(Element_t *el,double **u,double *f)
/* Termes explicites (va)  */
{
  double *va = Element_GetExplicitTerm(el) ;
  int i ;
  
  /* initialisation */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;
  
  
  /* Transfer coefficients */
  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,f,0,i) ;
    
    double p_l    = x[I_P_L] ;
    double p_g    = x[I_P_G] ;
    double p_c    = p_g - p_l ;

    double s_l    = x[I_S_L] ;
    double s_g    = 1 - s_l ;

    double rho_v  = x[I_RHO_V] ;
    double rho_a  = x[I_RHO_A] ;
    double rho_g  = x[I_RHO_G] ;

    double c_v    = rho_v/rho_g ;
    double c_a    = 1 - c_v ;
  
    double k_rl  = RelativePermeabilityToLiquid(p_c) ;
    double k_rg  = RelativePermeabilityToGas(p_c) ;
  
    double kh_l  = k_int/mu_l*k_rl ;
    double kh_g  = k_int/mu_g*k_rg ;
  
    double tau   = TortuosityToGas(phi,s_g) ;
    double D_av  = D_av0*p_g0/p_g ;
    double D_eff = phi*s_g*tau*D_av ;
  
    /* Darcy */
    K_L   += rho_l*kh_l ;           /* liquid */
    K_G   += rho_g*kh_g ;           /* gas */
    KD_V  += rho_v*kh_g ;           /* vapor */
    KD_A  += rho_a*kh_g ;           /* air */
    
    /* Fick */
    KF_V  += rho_g*D_eff ;          /* vapor */
    KF_A  += rho_g*D_eff ;          /* air */
    
    /* Barometric diffusion effects */
    KD_V += D_eff*c_v*c_a*(M_air - M_h2o)/RT ; /* vapor */
    KD_A += D_eff*c_v*c_a*(M_h2o - M_air)/RT ; /* air */
  }
  
  /* Averaging */
  for(i = 0 ; i < NVE ; i++) va[i] *= 0.5 ;
}




void ComputeFluxes(Element_t *el,double **u)
/* Fluxes */
{
  double *f = Element_GetImplicitTerm(el) ;
  double *grd = dComponents ;


  {
    double *x1 = ComputeComponents(el,u,f,0.,1) ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] = x1[i] ;
  }
  {
    double *x0 = ComputeComponents(el,u,f,0.,0) ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] -= x0[i] ;
  }
  
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx = x1 - x0 ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] /= dx ;
  }
  
  /* Fluxes */
  {
    double *w = Fluxes(el,grd) ;

    W_L     = w[I_W_L] ;
    W_G     = w[I_W_G] ;
    W_V     = w[I_W_V] ;
    W_A     = w[I_W_A] ;
    W_T     = w[I_W_T] ;
  }
}


double* Fluxes(Element_t *el,double *grd)
{
  double* va = Element_GetExplicitTerm(el) ;
  double* w   = ComponentFluxes ;

  /* Gradients */
  double grd_h_l    = grd[I_H_L] ;
  double grd_h_g    = grd[I_H_G] ;
  double grd_c_v    = grd[I_C_V] ;
  double grd_c_a    = - grd[I_C_V] ;
    
    
  /* Flux */
  double w_l = - K_L*grd_h_l ;
  double w_g = - K_G*grd_h_g ;
  double w_v = - KD_V*grd_h_g - KF_V*grd_c_v ;
  double w_a = - KD_A*grd_h_g - KF_A*grd_c_a ;
  double w_t = w_l + w_g ;
   
   
  w[I_W_L ]  = w_l ;
  w[I_W_G ]  = w_g ;
  w[I_W_V ]  = w_v ;
  w[I_W_A ]  = w_a ;
  w[I_W_T ]  = w_t ;
   
  return(w) ;
}


double saturation(double pc,double pc3,Curve_t* curve)
/* Degre de saturation regularise autour de 1 */
{
  double* x = Curve_GetXRange(curve) ;
  double* y = Curve_GetYValue(curve) ;
  double pc1 = x[0] ;
  double sl ;
  
  if(pc >= pc3 || pc1 >= pc3) {
    sl = Curve_ComputeValue(curve,pc) ;
  } else {
    double sl1 = y[0] ;
    double sl3 = Curve_ComputeValue(curve,pc3) ;
    
    sl  = sl1 - (sl1 - sl3) * exp(-(pc - pc3)/(pc1 - pc3)) ;
  }
  
  return(sl) ;
}

double dsaturation(double pc,double pc3,Curve_t* curve)
{
  int   n_i = Curve_GetNbOfPoints(curve) - 1 ;
  double* x = Curve_GetXRange(curve) ;
  double dpc = (x[1] - x[0])/n_i ;
  double sl  = saturation(pc,pc3,curve) ;
  double sl1 = saturation(pc + dpc,pc3,curve) ;
  double dsl = sl1 - sl ;

  return(dsl/dpc) ;
}



int TangentCoefficients(Element_t *el,double dt,double *c)
/**  Tangent matrix coefficients (c) */
{
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  ObVal_t *obval = Element_GetObjectiveValue(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  double* dist = FVM_ComputeIntercellDistances(fvm) ;
  double dij   = dist[1] ;
  double dtdij = dt/dij ;
  int    dec = NEQ*NEQ ;
  int    i ;
  
  /* Initialization */
  for(i = 0 ; i < ndof*ndof ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  for(i = 0 ; i < 2 ; i++) {
    int j = 1 - i ;
    /* Content terms at node i */
    double *cii = c + (i*nn + i)*NEQ*NEQ ;
    /* Transfer terms from node i to node j */
    double *cij = c + (i*nn + j)*NEQ*NEQ ;
    /* Components */
    double *x         = ComputeComponents(el,u,f_n,dt,i) ;
    double dxi[NEQ] ;
    int k ;
    
    for(k = 0 ; k < NEQ ; k++) {
      /* dxi[k] =  1.e-1*ObVal_GetValue(obval + k) ; */
      dxi[k] =  1.e2 ;
    }
    
    for(k = 0 ; k < NEQ ; k++) {
      double dxk    = dxi[k] ;
      double *dx    = ComputeComponentDerivatives(el,dt,x,dxk,k) ;
      double *dw    = Fluxes(el,dx) ;
    
      cii[E_Mass*NEQ   + k] = dx[I_M_T] ;
      cii[E_Air*NEQ    + k] = dx[I_M_A] ;
      
      cij[E_Mass*NEQ   + k] = - dtdij*dw[I_W_T] ;
      cij[E_Air*NEQ    + k] = - dtdij*dw[I_W_A] ;
    }
  }

  return(dec) ;
}





double* ComputeComponents(Element_t *el,double **u,double *f_n,double dt,int n)
{
  double *x = Components ;
  
  /* Primary Variables */
  x[I_P_L ] = P_L(n) ;
  x[I_P_A ] = P_A(n) ;
  
  /* Needed variables to compute secondary components */
  x[I_COOR_X] = Element_GetNodeCoordinate(el,n)[0] ;
    
  ComputeSecondaryComponents(el,dt,x) ;
  return(x) ;
}


double* ComputeComponentDerivatives(Element_t *el,double dt,double *x,double dxi,int i)
{
  double *dx = dComponents ;
  int j ;
  
  /* Primary Variables */
  dx[I_P_L ] = x[I_P_L ] ;
  dx[I_P_A ] = x[I_P_A ] ;
  
  /* Needed variables to compute secondary components */
  dx[I_COOR_X] = x[I_COOR_X] ;
  
  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;
  
  ComputeSecondaryComponents(el,dt,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfComponents ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}



void  ComputeSecondaryComponents(Element_t *el,double dt,double *x)
{
  double p_l     = x[I_P_L] ;
  double p_a     = x[I_P_A] ;
  double z       = x[I_COOR_X];
    
  /* Fluid components */
  double p_v     = VaporPressure(p_l) ;
  double p_g     = p_v + p_a ;
  double p_c     = p_g - p_l ;

  double s_l     = Saturation(p_c) ;
  double s_g     = 1 - s_l ;

  double rho_v   = p_v*M_h2o/RT ;
  double rho_a   = p_a*M_air/RT ;
  double rho_g   = rho_v + rho_a ;
  
  double c_v     = rho_v/rho_g ;

  double m_l     = rho_l*phi*s_l ;
  double m_v     = rho_v*phi*s_g ;
  double m_a     = rho_a*phi*s_g ;
  double m_g     = m_a + m_v ;
  double m_t     = m_l + m_g ;
    
    
  /* Back up */
  
  /* Fluid components */
  x[I_S_L      ] = s_l ;
  x[I_RHO_V    ] = rho_v ;
  x[I_RHO_A    ] = rho_a ;
  x[I_RHO_G    ] = rho_g ;
  x[I_H_L      ] = p_l - rho_l*gravite*z ;
  x[I_H_G      ] = p_g - rho_g*gravite*z ;
  x[I_C_V      ] = c_v ;
  x[I_P_V      ] = p_v ;
  x[I_P_G      ] = p_g ;
  
  /* Fluid mass contents */
  x[I_M_L      ] = m_l ;
  x[I_M_V      ] = m_v ;
  x[I_M_A      ] = m_a ;
  x[I_M_G      ] = m_g ;
  x[I_M_T      ] = m_t ;
}
